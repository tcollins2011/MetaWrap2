"""`metawrap2 benchmark`: its output and its refusals.

The command can download 95 GB and run for hours, so the tests check that it does nothing until
the dataset is named and the size is confirmed.
"""

import os

import pytest

from metawrap2.benchmark import scoring
from metawrap2.commands import benchmark as cmd

BENCHMARKS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks")

#: Captured at import, before the autouse fixture below replaces it, so one test can put the real
#: guard back and assert it.
_REAL_REQUIRE_ENV_MANAGER = cmd._runner.require_env_manager


@pytest.fixture(autouse=True)
def _use_the_shipped_registry(monkeypatch):
    """Point the registry at the repository's own benchmarks/ regardless of the cwd."""
    monkeypatch.setenv("METAWRAP2_BENCHMARKS", os.path.abspath(BENCHMARKS))


@pytest.fixture(autouse=True)
def _dont_require_conda(monkeypatch):
    """Disable the env-manager guard, since these tests fake the executor.

    The guard is tested by :func:`test_run_without_an_env_manager_stops_before_the_first_variant`,
    which restores it.
    """
    monkeypatch.setattr(cmd._runner, "require_env_manager", lambda: None)


def test_score_labels_a_per_sample_module_by_its_sample_directories(tmp_path, capsys):
    """read_qc output sits one directory deeper per sample, and both mates are summed.

    A measurement is keyed on (dataset, comparison, variant, sample, metric), not on the source
    file. Guards against the two mates collapsing into one row, which halved every read_qc figure.
    Two samples x two mates x one read gives four.
    """
    from metawrap2.benchmark import scoring

    root = tmp_path / "work" / "trimmer" / "fastp"
    for sample in ("a", "b"):
        (root / sample).mkdir(parents=True)
        for mate in ("1", "2"):
            (root / sample / ("final_pure_reads_%s.fastq" % mate)).write_text("@r\nACGT\n+\nIIII\n")

    results = tmp_path / "r.tsv"
    scoring.write_tsv(
        [
            scoring.Row(
                dataset="tiny",
                comparison="trimmer",
                variant="fastp",
                metric="cpu_hours",
                value=0.5,
                sample=sample,
            )
            for sample in ("a", "b")
        ],
        str(results),
    )
    rc = cmd.main(
        [
            "score",
            "tiny",
            "trimmer",
            "--work",
            str(tmp_path / "work"),
            "--results",
            str(results),
            "--metric",
            "reads_kept",
        ]
    )
    out = capsys.readouterr().out
    assert rc == 0
    rows = scoring.read_tsv(str(results))
    kept = [r for r in rows if r.metric == "reads_kept"]
    assert sorted(r.sample for r in kept) == ["a", "b"]
    # One row per sample, covering both mates.
    assert len(kept) == 2
    assert all(r.value == 2.0 for r in kept)
    # Two reads and half a CPU-hour per sample: four reads and one CPU-hour in total.
    (line,) = [ln for ln in out.splitlines() if "reads_kept" in ln and "fastp" in ln]
    # fastp is the default, so the line also carries the incumbent marker.
    assert line.split()[:4] == ["fastp", "reads_kept", "4", "1"]


@pytest.mark.conda_envs
def test_run_without_an_env_manager_stops_before_the_first_variant(tmp_path, monkeypatch, capsys):
    """Without an env manager: one message, exit 1, and no variant started."""
    monkeypatch.setattr(cmd._runner, "require_env_manager", _REAL_REQUIRE_ENV_MANAGER)
    monkeypatch.setattr(cmd._runner.shutil, "which", lambda _name: None)

    def explode(*_a, **_k):
        raise AssertionError("started a variant with no env manager to run it in")

    monkeypatch.setattr(cmd._runner._usage, "run_measured", explode)
    rc = cmd.main(["run", "tiny", "binner_set", "-o", str(tmp_path / "work")])
    out = capsys.readouterr().out
    assert rc == 1
    assert "not on PATH" in out
    assert not os.path.exists(tmp_path / "work" / "binner_set.run.json")


# ── chimerism at score time ────────────────────────────────────────────────────────────────────


def _gunc_report(directory, passes, fails):
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "GUNC.progenomes_2.1.maxCSS_level.tsv")
    verdicts = ["p%d\tTrue\n" % i for i in range(passes)]
    verdicts += ["f%d\tFalse\n" % i for i in range(fails)]
    with open(path, "w") as fh:
        fh.write("genome\tpass.GUNC\n")
        fh.writelines(verdicts)


def test_a_score_time_report_is_labelled_with_the_bin_set_it_sits_beside(tmp_path):
    """GUNC reports in ``<bin set>.gunc`` are labelled by bin set, taken from the path.

    The ranked set keeps the plain metric name and the others get a suffix, as AMBER's per-set
    rows do, so GUNC and AMBER numbers line up per bin set.
    """
    _gunc_report(str(tmp_path / "metabat2_bins.gunc"), 4, 0)
    _gunc_report(str(tmp_path / "maxbin2_bins.gunc"), 1, 1)
    rows = cmd._collect_scores(
        "tiny", "binner_set", "v", str(tmp_path), ranked_bin_set="metabat2_bins"
    )
    by_metric = {r.metric: r.value for r in rows}
    assert by_metric["gunc_pass_fraction"] == 1.0
    assert by_metric["gunc_pass_fraction__maxbin2_bins"] == 0.5
    # The ranker sums by metric name; without the suffix both rows would be added, giving 1.5.
    ranked = scoring.value_against_cost(
        rows + [scoring.Row("tiny", "binner_set", "v", "cpu_hours", 2.0)], "gunc_pass_fraction"
    )
    assert ranked[0].value == 1.0


# ── an arm that is not MetaWrap2 ───────────────────────────────────────────────────────────────

#: A registry with one native arm and one external arm. The shipped file has this arm commented
#: out; here it is live, so an external arm ranked beside a native one is covered.
_EXTERNAL_CONFIG = """
schema = 1

[binner_set]
module = "binning"
question = "Which binners should run by default?"
metrics = ["bins_high_quality"]

[[binner_set.variants]]
name = "metabat2_semibin2"
args = ["--metabat2", "--semibin2"]
default = true

[[binner_set.variants]]
name = "mw1"

[binner_set.variants.external]
command = "metawrap binning -o {output}/out -t {threads} -a {assembly} {reads}"
bins = "out/metabat2_bins"
env = "metawrap-1.3.0"
version = "1.3.0"
"""


def _stats_table(path, rows):
    """A MetaWrap2 ``.stats.tsv`` from CheckM2, so scoring reads numbers without running CheckM."""
    with open(path, "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
        fh.writelines("%s\t%d\t%d\t0.5\tcheckm2\t1000\t2000\n" % row for row in rows)


def test_an_external_arm_is_ranked_with_its_cost_like_any_other(tmp_path, monkeypatch, capsys):
    """``score`` ranks an external arm with its cost, and prints the arm's version."""
    registry_dir = tmp_path / "benchmarks"
    registry_dir.mkdir()
    (registry_dir / "configurations.toml").write_text(_EXTERNAL_CONFIG)
    for name in ("datasets.toml", "fixtures"):
        os.symlink(os.path.join(os.path.abspath(BENCHMARKS), name), registry_dir / name)
    monkeypatch.setenv("METAWRAP2_BENCHMARKS", str(registry_dir))

    work = tmp_path / "work"
    for arm, bins in (("metabat2_semibin2", [("bin_1", 95, 1)]), ("mw1", [("bin.1", 80, 4)])):
        root = work / "binner_set" / arm
        root.mkdir(parents=True)
        bin_dir = root / ("%s_bins" % arm)
        bin_dir.mkdir()
        for stem, _completeness, _contamination in bins:
            # A real FASTA: the ranked bin set is the one found on disk, and an empty directory is
            # not counted as a bin set.
            (bin_dir / ("%s.fasta" % stem)).write_text(">contig_1\nACGTACGTAC\n")
        _stats_table(str(root / ("%s_bins.stats.tsv" % arm)), bins)
    (work / "binner_set.run.json").write_text(
        '{"samples": ["sample_1"], "variants": {"metabat2_semibin2": {"returncode": 0}, '
        '"mw1": {"returncode": 0, "external": true, "external_version": "1.3.0"}}}'
    )
    results = tmp_path / "r.tsv"
    scoring.write_tsv(
        [
            scoring.Row("tiny", "binner_set", "metabat2_semibin2", "cpu_hours", 43.8),
            scoring.Row("tiny", "binner_set", "mw1", "cpu_hours", 52.5),
        ],
        str(results),
    )

    rc = cmd.main(["score", "tiny", "binner_set", "--work", str(work), "--results", str(results)])
    out = capsys.readouterr().out

    assert rc == 0
    (external_row,) = [
        line
        for line in out.splitlines()
        if line.strip().startswith("mw1 ") and "bins_high_quality" in line
    ]
    assert "unmeasured" not in external_row
    assert "52.5" in external_row
    # The external arm and its version are printed, not only recorded.
    assert "External arm(s)" in out and "1.3.0" in out


def test_score_caches_and_reuses_scores_per_scorer(tmp_path, monkeypatch):
    """CheckM2 is the default, and a CheckM1 table or cache is not read as CheckM2 scores."""
    from metawrap2 import checkm as checkm_mod
    from metawrap2.benchmark import assess

    bin_set = str(tmp_path / "metabat2_bins")
    os.makedirs(bin_set)
    _stats_table(bin_set + ".stats.tsv", [("bin_1", 95, 1)])
    with open(bin_set + ".stats.tsv") as fh:
        checkm2_table = fh.read()
    with open(bin_set + ".stats.tsv", "w") as fh:
        fh.write(checkm2_table.replace("checkm2", "k__Bacteria"))
    assert cmd._already_assessed(bin_set, ["bins_high_quality"], "checkm1")
    assert not cmd._already_assessed(bin_set, ["bins_high_quality"], "checkm2")

    caches = []
    monkeypatch.setattr(
        checkm_mod,
        "run_checkm_cached",
        lambda path, cache, **kw: caches.append(os.path.basename(cache.path)) or path,
    )
    assess.run_bin_quality([bin_set], threads=1, memory_gb=8)
    assess.run_bin_quality([bin_set], threads=1, memory_gb=8, checkm2=False)
    assert caches == [
        "metabat2_bins.checkm_scores.checkm2.json",
        "metabat2_bins.checkm_scores.checkm1.json",
    ]
