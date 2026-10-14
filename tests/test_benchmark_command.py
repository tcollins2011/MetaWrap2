"""`metawrap2 benchmark` as a user meets it: what it prints, and what it refuses to do.

The refusals are the point. This command can spend 95 GB and an afternoon, so the tests that matter
are the ones asserting it does nothing until told twice - once by naming the dataset, once by
confirming the size.
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
    """These tests fake the executor, so whether conda is installed here is beside the point.

    The guard itself is asserted by
    :func:`test_run_without_an_env_manager_stops_before_the_first_variant`, which turns it back on.
    """
    monkeypatch.setattr(cmd._runner, "require_env_manager", lambda: None)


def test_score_labels_a_per_sample_module_by_its_sample_directories(tmp_path, capsys):
    """read_qc runs once per sample, so its output sits one directory deeper - and sums.

    Both mates are written, which is the case that found a defect: a measurement is keyed on
    (dataset, comparison, variant, sample, metric) and deliberately not on the source file, so two
    rows for one sample collapsed to whichever mate was walked last and every read_qc figure was
    half the reads. Two samples x two mates x one read is four, not two.
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
    # One row per sample, each carrying both mates rather than one row per file.
    assert len(kept) == 2
    assert all(r.value == 2.0 for r in kept)
    # Two reads per sample over half a CPU-hour each: four reads, and one CPU-hour beside them.
    (line,) = [ln for ln in out.splitlines() if "reads_kept" in ln and "fastp" in ln]
    # fastp is the shipped default now, so the line also carries the incumbent marker.
    assert line.split()[:4] == ["fastp", "reads_kept", "4", "1"]


@pytest.mark.conda_envs
def test_run_without_an_env_manager_stops_before_the_first_variant(tmp_path, monkeypatch, capsys):
    """One sentence, exit 1, and nothing spawned - not ten recorded failures and an empty table."""
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
    """Score time runs GUNC per bin set into ``<bin set>.gunc``, next to CheckM's ``<bin set>.checkm``.

    So three binners leave three reports, and the label has to come from the path or they collide.
    The ranked set keeps the plain metric name and the others are suffixed, which is the same rule
    AMBER's per-set rows already follow - and it means the two tools' numbers line up on one set.
    """
    _gunc_report(str(tmp_path / "metabat2_bins.gunc"), 4, 0)
    _gunc_report(str(tmp_path / "maxbin2_bins.gunc"), 1, 1)
    rows = cmd._collect_scores(
        "tiny", "binner_set", "v", str(tmp_path), ranked_bin_set="metabat2_bins"
    )
    by_metric = {r.metric: r.value for r in rows}
    assert by_metric["gunc_pass_fraction"] == 1.0
    assert by_metric["gunc_pass_fraction__maxbin2_bins"] == 0.5
    # Which is the point of the suffix: the ranker is over metric names, so only one is summable.
    # Without it both rows would be added and the variant would rank at a pass fraction of 1.5.
    ranked = scoring.value_against_cost(
        rows + [scoring.Row("tiny", "binner_set", "v", "cpu_hours", 2.0)], "gunc_pass_fraction"
    )
    assert ranked[0].value == 1.0
