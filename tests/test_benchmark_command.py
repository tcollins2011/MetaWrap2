"""`metawrap2 benchmark` as a user meets it: what it prints, and what it refuses to do.

The refusals are the point. This command can spend 95 GB and an afternoon, so the tests that matter
are the ones asserting it does nothing until told twice - once by naming the dataset, once by
confirming the size.
"""

import os

import pytest

from metawrap2 import checkm as _checkm
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


def test_list_names_the_datasets_the_comparisons_and_the_current_defaults(capsys):
    assert cmd.main(["list"]) == 0
    out = capsys.readouterr().out
    assert "cami2_gi" in out
    assert "tiny" in out
    assert "binner_set" in out
    # Someone reading this has to be able to see what a challenger has to beat.
    assert "current default" in out
    # And the download size, before they type fetch.
    assert "GB" in out


def test_list_verbose_shows_the_citation_and_the_gold_standard_files(capsys):
    assert cmd.main(["list", "-v"]) == 0
    out = capsys.readouterr().out
    assert "cite:" in out
    assert "gold:" in out


def test_fetch_prints_the_size_and_downloads_nothing_without_yes(tmp_path, monkeypatch, capsys):
    """Typing a dataset name is not the same as meaning to spend 95 GB."""

    def explode(*_a, **_k):
        raise AssertionError("downloaded without --yes")

    monkeypatch.setattr(cmd._fetch, "_download", explode)
    rc = cmd.main(["fetch", "cami2_gi", "--dest", str(tmp_path / "d")])
    out = capsys.readouterr().out
    assert rc == 1
    assert "--yes" in out
    assert "95." in out  # the size, stated before the confirmation is asked for
    assert not os.path.exists(tmp_path / "d")


def test_fetch_prints_the_licence_and_the_citation_before_asking(tmp_path, capsys):
    """The data is someone else's; the terms travel with it rather than being looked up later."""
    cmd.main(["fetch", "cami2_gi", "--dest", str(tmp_path / "d")])
    out = capsys.readouterr().out
    assert "Cite:" in out
    assert "CC BY 4.0" in out


def test_fetch_refuses_before_confirming_when_it_cannot_fit(tmp_path, monkeypatch, capsys):
    """No point confirming a download that cannot land."""
    from metawrap2 import limits

    monkeypatch.setattr(
        cmd._fetch._limits,
        "disk_state",
        lambda path, budget: limits.DiskState(
            path=str(tmp_path), total_bytes=1000, used_bytes=999, ceiling_bytes=999
        ),
    )
    rc = cmd.main(["fetch", "cami2_gi", "--dest", str(tmp_path / "d"), "--yes"])
    out = capsys.readouterr().out
    assert rc == 1
    assert "disk ceiling" in out


def test_fetch_of_the_committed_fixture_needs_no_confirmation(capsys):
    """It is already in the repository; there is nothing to confirm."""
    assert cmd.main(["fetch", "tiny"]) == 0
    assert "nothing to do" in capsys.readouterr().out


def test_an_unknown_dataset_is_named_alongside_the_known_ones(capsys):
    assert cmd.main(["fetch", "no_such_dataset"]) == 2
    assert "cami2_gi" in capsys.readouterr().out


def test_an_unknown_comparison_is_named_alongside_the_known_ones(capsys):
    assert cmd.main(["score", "tiny", "no_such_comparison"]) == 2
    assert "binner_set" in capsys.readouterr().out


def test_run_on_absent_data_says_how_to_fetch_it_rather_than_fetching_it(tmp_path, capsys):
    rc = cmd.main(
        [
            "run",
            "cami2_gi",
            "binner_set",
            "--data",
            str(tmp_path / "absent"),
            "-o",
            str(tmp_path / "work"),
        ]
    )
    out = capsys.readouterr().out
    assert rc == 1
    assert "benchmark fetch cami2_gi" in out


def test_score_with_no_results_yet_says_so_instead_of_printing_an_empty_table(tmp_path, capsys):
    rc = cmd.main(
        [
            "score",
            "tiny",
            "binner_set",
            "--work",
            str(tmp_path / "work"),
            "--results",
            str(tmp_path / "r.tsv"),
        ]
    )
    out = capsys.readouterr().out
    assert rc == 1
    assert "No results file" in out


def test_score_reads_a_variants_checkm_report_and_ranks_it(tmp_path, capsys):
    """The whole path: tool output on disk in, value per CPU-hour out."""
    from metawrap2.benchmark import scoring

    root = tmp_path / "work" / "binner_set" / "metabat2"
    root.mkdir(parents=True)
    (root / "quality_report.tsv").write_text(
        "Name\tCompleteness\tContamination\nb1\t96.0\t1.0\nb2\t93.0\t2.0\n"
    )
    results = tmp_path / "r.tsv"
    # The cost half, as a run would have recorded it.
    scoring.write_tsv(
        [
            scoring.Row(
                dataset="tiny",
                comparison="binner_set",
                variant="metabat2",
                metric="cpu_hours",
                value=2.0,
            )
        ],
        str(results),
    )

    rc = cmd.main(
        [
            "score",
            "tiny",
            "binner_set",
            "--work",
            str(tmp_path / "work"),
            "--results",
            str(results),
            "--metric",
            "bins_near_complete",
        ]
    )
    out = capsys.readouterr().out
    assert rc == 0
    assert "with cpu_hours beside it" in out
    # Two near-complete bins, and the two CPU-hours they cost reported next to them.
    (line,) = [ln for ln in out.splitlines() if "bins_near_complete" in ln and "metabat2" in ln]
    assert line.split() == ["metabat2", "bins_near_complete", "2", "2"]


@pytest.mark.parametrize("flag,expected", [([], False), (["--checkm2"], True)])
def test_checkm2_reaches_the_assessor_only_when_asked(tmp_path, monkeypatch, flag, expected):
    """``--checkm2`` has to travel score -> _assess_variant -> assess_bins to be worth anything.

    Scoring cami2_gi's binner_set on 2026-10-01 lost two of its four declared metrics because
    CheckM1's pplacer stage died with a multiprocessing ``EOFError`` on every variant, and there was
    no way to ask for CheckM2 from the command line even though the assessor had supported it all
    along. A flag that is parsed and then dropped looks exactly like a flag that works.
    """
    root = tmp_path / "work" / "binner_set" / "metabat2"
    (root / "metabat2_bins").mkdir(parents=True)
    (root / "metabat2_bins" / ("bin.1" + _checkm.BIN_EXTENSION)).write_text(">c1\nACGT\n")

    seen = {}

    def spy(_root, **kw):
        seen.update(kw)
        return [str(root / "metabat2_bins")], {}

    monkeypatch.setattr(cmd._assess, "assess_bins", spy)
    rc = cmd.main(
        ["score", "tiny", "binner_set", "--work", str(tmp_path / "work"), "--metric", "amber_f1_bp"]
        + flag
    )
    assert rc == 0
    assert seen.get("checkm2") is expected


def test_score_ignores_rows_from_another_dataset(tmp_path, capsys):
    """One results file holds every dataset; mixing them would compare across machines."""
    from metawrap2.benchmark import scoring

    results = tmp_path / "r.tsv"
    scoring.write_tsv(
        [
            scoring.Row(
                dataset=name,
                comparison="binner_set",
                variant="metabat2",
                metric=metric,
                value=value,
            )
            for name, metric, value in (
                ("tiny", "bins_near_complete", 2.0),
                ("tiny", "cpu_hours", 1.0),
                ("cami2_gi", "bins_near_complete", 900.0),
                ("cami2_gi", "cpu_hours", 1.0),
            )
        ],
        str(results),
    )
    cmd.main(
        [
            "score",
            "tiny",
            "binner_set",
            "--work",
            str(tmp_path / "work"),
            "--results",
            str(results),
            "--metric",
            "bins_near_complete",
        ]
    )
    out = capsys.readouterr().out
    assert "900" not in out


def test_a_dry_run_prints_the_commands_and_spawns_nothing(tmp_path, monkeypatch, capsys):
    """--dry-run is a global flag the dispatcher strips, so this is also the test that it arrives."""
    from metawrap2 import command as shared

    monkeypatch.setattr(shared.runner, "dry_run", True)

    def explode(*_a, **_k):
        raise AssertionError("ran a module during a dry run")

    monkeypatch.setattr(cmd._runner._usage, "run_measured", explode)
    rc = cmd.main(["run", "tiny", "binner_set", "-o", str(tmp_path / "work")])
    out = capsys.readouterr().out
    assert rc == 0
    assert "metawrap2.cli binning" in out
    assert "--metabat2" in out
    assert not os.path.exists(tmp_path / "work" / "binner_set.run.json")


def test_a_real_run_records_its_cost_in_the_results_tsv(tmp_path, monkeypatch, capsys):
    from metawrap2.benchmark import scoring

    monkeypatch.setattr(
        cmd._runner._usage,
        "run_measured",
        lambda argv, **_kw: (
            0,
            {
                "wall_seconds": 60.0,
                "user_seconds": 3600.0,
                "system_seconds": 0.0,
                "peak_rss_bytes": 1024,
            },
        ),
    )
    results = tmp_path / "r.tsv"
    rc = cmd.main(
        [
            "run",
            "tiny",
            "binner_set",
            "-o",
            str(tmp_path / "work"),
            "--variant",
            "metabat2",
            "--results",
            str(results),
        ]
    )
    assert rc == 0
    assert "Recorded" in capsys.readouterr().out
    rows = scoring.read_tsv(str(results))
    assert {r.metric for r in rows} == {"wall_seconds", "cpu_hours", "peak_rss_bytes"}
    assert next(r.value for r in rows if r.metric == "cpu_hours") == pytest.approx(1.0)


def test_running_one_named_variant_runs_only_that_one(tmp_path, monkeypatch):
    ran = []
    monkeypatch.setattr(
        cmd._runner._usage,
        "run_measured",
        lambda argv, **_kw: (ran.append(argv), (0, {"wall_seconds": 1.0}))[1],
    )
    cmd.main(
        [
            "run",
            "tiny",
            "binner_set",
            "-o",
            str(tmp_path / "work"),
            "--variant",
            "concoct",
            "--results",
            str(tmp_path / "r.tsv"),
        ]
    )
    assert len(ran) == 1
    assert "--concoct" in ran[0]


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


# ── abundance scoring ──────────────────────────────────────────────────────────────────────────


def _abundance_table(path, rows):
    path.write_text(
        "Genomic bins\tsample_0\tsample_1\n"
        + "".join("%s\t%s\t%s\n" % (name, a, b) for name, a, b in rows)
    )


def test_collect_scores_scores_an_abundance_table_against_the_truth_it_is_handed(tmp_path):
    """The one metric not read off the run's own output: the table says what the quantifier thought,
    and only the dataset's read mappings say what was true."""
    _abundance_table(tmp_path / "bin_abundance_table.tsv", [("bin.1", 1, 9), ("bin.2", 9, 1)])
    truth = {"sample_0": {"bin.1": 0.1, "bin.2": 0.9}, "sample_1": {"bin.1": 0.9, "bin.2": 0.1}}
    rows = cmd._collect_scores("tiny", "quantifier", "sylph", str(tmp_path), abundance_truth=truth)
    by_metric = {r.metric: r.value for r in rows}
    assert by_metric["abundance_spearman"] == 1.0
    assert by_metric["abundance_l1_similarity"] == 1.0
    assert all(r.source == "bin_abundance_table.tsv" for r in rows)


def _gunc_report(directory, passes, fails):
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "GUNC.progenomes_2.1.maxCSS_level.tsv")
    verdicts = ["p%d\tTrue\n" % i for i in range(passes)]
    verdicts += ["f%d\tFalse\n" % i for i in range(fails)]
    with open(path, "w") as fh:
        fh.write("genome\tpass.GUNC\n")
        fh.writelines(verdicts)


def test_the_final_sets_gunc_report_is_the_variants_own_figure(tmp_path):
    """GUNC names the database in its filename and never the bins, so the directory is the only
    evidence of which set a report describes - and bin_refinement writes the final set's to one
    defined place."""
    _gunc_report(str(tmp_path / "work_files" / "gunc"), 3, 1)
    rows = cmd._collect_scores(
        "tiny",
        "binner_set",
        "metabat2_maxbin2",
        str(tmp_path),
        ranked_bin_set="metawrap_70_10_bins",
    )
    assert [(r.metric, r.value) for r in rows] == [("gunc_pass_fraction", 0.75)]


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


def test_a_report_in_no_recognised_place_is_kept_but_not_ranked(tmp_path):
    """The conservative end of the rule: a fraction from an unknown set of bins is still a number
    worth having in the table, and must not be added to the variant's own."""
    _gunc_report(str(tmp_path / "work_files" / "gunc"), 4, 0)
    _gunc_report(str(tmp_path / "work_files" / "scratch"), 1, 1)
    rows = cmd._collect_scores("tiny", "binner_set", "v", str(tmp_path), ranked_bin_set="final")
    by_metric = {r.metric: r.value for r in rows}
    assert by_metric["gunc_pass_fraction"] == 1.0
    assert by_metric[os.path.join("gunc_pass_fraction__work_files", "scratch")] == 0.5


def test_a_per_sample_gunc_report_is_still_the_samples_own_figure(tmp_path):
    """The sample name sits between the variant root and work_files, and must not be read as a second
    bin set - otherwise every sample's report would be suffixed and none of them rankable."""
    _gunc_report(str(tmp_path / "sample_0" / "work_files" / "gunc"), 1, 1)
    rows = cmd._collect_scores(
        "tiny", "binner_set", "v", str(tmp_path), per_sample=True, ranked_bin_set="final"
    )
    assert [(r.metric, r.sample, r.value) for r in rows] == [
        ("gunc_pass_fraction", "sample_0", 0.5)
    ]


def test_collect_scores_leaves_an_abundance_table_alone_without_a_truth(tmp_path):
    """Scoring a quantifier against nothing would be a number with no meaning behind it."""
    _abundance_table(tmp_path / "bin_abundance_table.tsv", [("bin.1", 1, 9)])
    assert cmd._collect_scores("tiny", "quantifier", "sylph", str(tmp_path)) == []


def test_an_unscorable_abundance_table_is_reported_and_does_not_stop_the_rest(tmp_path, capsys):
    """A truth naming none of the table's samples is a mismatch worth saying out loud, but the
    variant's CheckM numbers are still worth having."""
    _abundance_table(tmp_path / "bin_abundance_table.tsv", [("bin.1", 1, 9), ("bin.2", 9, 1)])
    rows = cmd._collect_scores(
        "tiny", "quantifier", "sylph", str(tmp_path), abundance_truth={"other": {"bin.1": 1.0}}
    )
    assert rows == []
    assert "could not be scored" in capsys.readouterr().out


def test_the_results_file_is_anchored_to_the_checkout_not_the_current_directory(
    tmp_path, monkeypatch
):
    """A relative default meant scoring from the wrong directory wrote a second results file.

    It reported success, the numbers in it were real, and the committed file was untouched - so the
    regression a committed results TSV exists to make visible was invisible, and the rows had to be
    merged back by hand.
    """
    registry_dir = tmp_path / "elsewhere" / "benchmarks"
    (registry_dir / "results").mkdir(parents=True)
    monkeypatch.setenv("METAWRAP2_BENCHMARKS", str(registry_dir))
    monkeypatch.chdir(tmp_path)

    path = cmd.default_results("cami2_gi")
    assert path == str(registry_dir / "results" / "cami2_gi.tsv")
    assert os.path.isabs(path), "a relative path is what caused this"
