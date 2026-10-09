"""Scoring and running: the tidy table, the parsers, and the ratio the defaults are decided on.

The arithmetic here is the whole basis for "MetaWrap2 picks defaults on accuracy per CPU-hour", so
these tests are mostly about that claim being true of the code: that cost is summed the same way
accuracy is, that a variant with no measured cost is not silently ranked first, and that the
incumbent is visible in the table someone reads to make the decision.
"""

import json
import os

import pytest

from metawrap2.benchmark import registry, runner, scoring


def row(variant, metric, value, sample="", comparison="binner_set", dataset="tiny"):
    return scoring.Row(
        dataset=dataset,
        comparison=comparison,
        variant=variant,
        metric=metric,
        value=value,
        sample=sample,
    )


# ── the tidy table ────────────────────────────────────────────────────────────────────────────


def test_rows_round_trip_through_the_tsv(tmp_path):
    path = str(tmp_path / "r.tsv")
    written = scoring.write_tsv([row("metabat2", "bins_near_complete", 12)], path)
    assert written == 1
    back = scoring.read_tsv(path)
    assert len(back) == 1
    assert back[0].variant == "metabat2"
    assert back[0].value == 12.0


def test_appending_does_not_write_a_second_header(tmp_path):
    """A results file accumulates over runs; a header in the middle makes it unreadable."""
    path = str(tmp_path / "r.tsv")
    scoring.write_tsv([row("a", "bins_total", 1)], path)
    scoring.write_tsv([row("b", "bins_total", 2)], path, append=True)
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    assert len(lines) == 3
    assert lines[0].startswith("dataset\t")
    assert sum(1 for line in lines if line.startswith("dataset\t")) == 1
    assert len(scoring.read_tsv(path)) == 2


def test_a_row_carries_the_unit_for_its_metric_without_being_told(tmp_path):
    """So nobody reading the TSV has to guess whether a duration is seconds or hours."""
    path = str(tmp_path / "r.tsv")
    scoring.write_tsv([row("a", "cpu_hours", 3.5)], path)
    assert scoring.read_tsv(path)[0].unit == "cpu_h"


def test_an_unparseable_value_is_skipped_rather_than_guessed_at(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text(
        "\t".join(scoring.COLUMNS)
        + "\ntiny\tc\tv\t\tbins_total\tn/a\t\t\ntiny\tc\tv\t\tbins_total\t4\t\t\n"
    )
    rows = scoring.read_tsv(str(path))
    assert [r.value for r in rows] == [4.0]


# ── the parsers ───────────────────────────────────────────────────────────────────────────────


CHECKM2 = "Name\tCompleteness\tContamination\nbin_001\t95.5\t1.2\nbin_002\t72.0\t8.0\nbin_003\t40.0\t2.0\n"
CHECKM1 = "Bin Id\tMarker lineage\tCompleteness\tContamination\nbin_001\tk__Bacteria\t95.5\t1.2\n"


def test_bin_quality_reads_a_checkm2_report(tmp_path):
    path = tmp_path / "quality_report.tsv"
    path.write_text(CHECKM2)
    assert scoring.bin_quality(str(path)) == [
        ("bin_001", 95.5, 1.2),
        ("bin_002", 72.0, 8.0),
        ("bin_003", 40.0, 2.0),
    ]


def test_bin_quality_also_reads_a_checkm1_report(tmp_path):
    """The two tools disagree on spelling; a benchmark that knew one would score nothing."""
    path = tmp_path / "checkm.tsv"
    path.write_text(CHECKM1)
    assert scoring.bin_quality(str(path)) == [("bin_001", 95.5, 1.2)]


def test_a_bin_checkm_could_not_judge_is_skipped_not_counted_as_zero(tmp_path):
    """Counting it as zero would make a failed bin look like a bad one."""
    path = tmp_path / "quality_report.tsv"
    path.write_text(CHECKM2 + "bin_004\tNA\tNA\n")
    assert len(scoring.bin_quality(str(path))) == 3


def test_quality_counts_use_the_community_thresholds(tmp_path):
    path = tmp_path / "quality_report.tsv"
    path.write_text(CHECKM2)
    counts = scoring.quality_counts(scoring.bin_quality(str(path)))
    assert counts["bins_near_complete"] == 1  # 95.5 / 1.2
    assert counts["bins_high_quality"] == 2  # and 72.0 / 8.0
    assert counts["bins_medium_quality"] == 2  # 40.0 is below 50
    assert counts["bins_total"] == 3


def test_a_contaminated_but_complete_bin_does_not_count_as_near_complete():
    assert scoring.quality_counts([("b", 99.0, 40.0)])["bins_near_complete"] == 0


def test_amber_metrics_from_a_wide_results_table(tmp_path):
    path = tmp_path / "results.tsv"
    path.write_text(
        "Tool\tF1 score (bp)\tARI (bp)\tAverage purity (bp)\n" "metabat2\t0.81\t0.77\t0.93\n"
    )
    metrics = scoring.amber_metrics(str(path))
    assert metrics["amber_f1_bp"] == 0.81
    assert metrics["amber_ari_bp"] == 0.77
    assert metrics["amber_purity_bp"] == 0.93


def test_amber_metrics_from_a_long_results_table(tmp_path):
    """AMBER has shipped both shapes; reading only one of them would score half the versions."""
    path = tmp_path / "results.tsv"
    path.write_text("tool\tmetric\tvalue\nmetabat2\tF1 score (bp)\t0.66\n")
    assert scoring.amber_metrics(str(path)) == {"amber_f1_bp": 0.66}


def test_an_unrecognised_amber_column_is_ignored_not_guessed_at(tmp_path):
    path = tmp_path / "results.tsv"
    path.write_text("Tool\tSomething New\tF1 score (bp)\nx\t9\t0.5\n")
    assert scoring.amber_metrics(str(path)) == {"amber_f1_bp": 0.5}


def test_gunc_pass_fraction_is_a_fraction_so_variants_are_comparable(tmp_path):
    """Variants produce different numbers of bins; a count of passes would not compare."""
    path = tmp_path / "GUNC.progenomes_2.1.maxCSS_level.tsv"
    path.write_text("genome\tpass.GUNC\nb1\tTrue\nb2\tFalse\nb3\tTrue\nb4\tTrue\n")
    assert scoring.gunc_pass_fraction(str(path)) == 0.75


def test_gunc_returns_none_when_the_file_says_nothing_about_passing(tmp_path):
    path = tmp_path / "GUNC.maxCSS_level.tsv"
    path.write_text("genome\tclade_separation_score\nb1\t0.1\n")
    assert scoring.gunc_pass_fraction(str(path)) is None


def test_assembly_stats_computes_n50_and_total(tmp_path):
    path = tmp_path / "final_assembly.fasta"
    # 100, 60, 40: total 200, so N50 is the contig at which the running sum passes 100 - the 100.
    path.write_text(">a\n%s\n>b\n%s\n>c\n%s\n" % ("A" * 100, "C" * 60, "G" * 40))
    stats = scoring.assembly_stats(str(path))
    assert stats["assembly_total_bp"] == 200
    assert stats["assembly_n50"] == 100
    assert stats["assembly_contigs"] == 3


def test_assembly_stats_reads_a_gzipped_assembly(tmp_path):
    import gzip

    path = tmp_path / "final_assembly.fasta.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">a\n%s\n" % ("A" * 50))
    assert scoring.assembly_stats(str(path))["assembly_total_bp"] == 50


def test_assembly_stats_counts_a_multi_line_contig_once(tmp_path):
    path = tmp_path / "a.fa"
    path.write_text(">a\nAAAA\nAAAA\n")
    stats = scoring.assembly_stats(str(path))
    assert stats["assembly_contigs"] == 1
    assert stats["assembly_total_bp"] == 8


# ── the ranking ───────────────────────────────────────────────────────────────────────────────


def test_the_ranking_is_on_quality_with_cost_reported_beside_it():
    """Both halves are reported; only the quality half decides the order."""
    rows = [
        row("cheap", "bins_near_complete", 10),
        row("cheap", "cpu_hours", 2.0),
        row("dear", "bins_near_complete", 11),
        row("dear", "cpu_hours", 8.0),
    ]
    ranked = scoring.value_against_cost(rows, "bins_near_complete")
    assert [r.variant for r in ranked] == ["dear", "cheap"], "11 bins beats 10 bins"
    assert ranked[0].value == 11 and ranked[0].cost == 8.0


def test_a_variant_that_scores_higher_for_far_more_compute_wins_the_ranking():
    """This assertion used to run the other way, and that is the correction.

    Ranking on quality-per-cost put the cheap arm first and MetaWrap2 shipped it as a default - which
    then recovered fewer MAGs end to end than the tool it replaces. Cost is for the reader to weigh
    against their own budget, and for `--preset` to express; it does not get to reorder the table.
    """
    rows = [
        row("incumbent", "amber_f1_bp", 0.80),
        row("incumbent", "cpu_hours", 1.0),
        row("challenger", "amber_f1_bp", 0.83),
        row("challenger", "cpu_hours", 4.0),
    ]
    ranked = scoring.value_against_cost(rows, "amber_f1_bp", defaults=["incumbent"])
    assert ranked[0].variant == "challenger"
    assert ranked[1].is_default, "the incumbent is still marked, just no longer first"


def test_arms_that_score_identically_are_separated_by_cost():
    """The only place cost decides anything: a tie on the metric."""
    rows = [
        row("dear", "amber_f1_bp", 0.80),
        row("dear", "cpu_hours", 9.0),
        row("cheap", "amber_f1_bp", 0.80),
        row("cheap", "cpu_hours", 1.0),
    ]
    assert [r.variant for r in scoring.value_against_cost(rows, "amber_f1_bp")] == ["cheap", "dear"]


def test_values_and_costs_are_both_summed_across_samples():
    """Ten samples give one row per variant, not ten - and both halves must sum the same way."""
    rows = [
        row("v", "bins_near_complete", 4, sample="sample_0"),
        row("v", "bins_near_complete", 6, sample="sample_1"),
        row("v", "cpu_hours", 1.0, sample="sample_0"),
        row("v", "cpu_hours", 3.0, sample="sample_1"),
    ]
    ranked = scoring.value_against_cost(rows, "bins_near_complete")
    assert ranked[0].value == 10
    assert ranked[0].cost == 4.0


def test_a_fraction_is_averaged_across_samples_rather_than_summed():
    """Summing a fraction does not just inflate it, it ranks by how many samples were measured.

    Two variants that both scored 0.9 purity everywhere would be separated by one having ten samples
    and the other two, and the ten-sample one would sort first at 9.0 against 1.8. The cost side is
    still a total, because ten samples really do cost ten samples' CPU-hours.
    """
    rows = [
        row("v", "amber_purity_bp", 0.9, sample="sample_0"),
        row("v", "amber_purity_bp", 0.7, sample="sample_1"),
        row("v", "cpu_hours", 1.0, sample="sample_0"),
        row("v", "cpu_hours", 3.0, sample="sample_1"),
    ]
    ranked = scoring.value_against_cost(rows, "amber_purity_bp")
    assert ranked[0].value == pytest.approx(0.8)
    assert ranked[0].samples == 2
    assert ranked[0].cost == 4.0


def test_averaging_a_fraction_is_said_out_loud_in_the_table():
    """An unweighted mean is a figure a reader may want to recompute, so the table admits taking one."""
    rows = [
        row("v", "gunc_pass_fraction", 1.0, sample="sample_0"),
        row("v", "gunc_pass_fraction", 0.5, sample="sample_1"),
        row("v", "cpu_hours", 2.0),
    ]
    text = scoring.format_scores(scoring.value_against_cost(rows, "gunc_pass_fraction"))
    assert "averaged over each variant's samples" in text and "v over 2" in text

    # Not said when there is nothing to say: one pooled row is its own mean.
    pooled = [row("v", "gunc_pass_fraction", 0.75), row("v", "cpu_hours", 2.0)]
    assert "averaged over" not in scoring.format_scores(
        scoring.value_against_cost(pooled, "gunc_pass_fraction")
    )


def test_a_variant_with_no_recorded_cost_still_ranks_on_its_quality():
    """A missing cost is a gap in the measurement, and it is no longer a reason to rank an arm last.

    Under the ratio an unmeasured cost could not be ranked at all, so a real result got pushed to the
    bottom of the table by a bookkeeping gap. Ranking on quality, it sits where its bins put it and the
    gap is simply reported.
    """
    rows = [
        row("measured", "bins_total", 5),
        row("measured", "cpu_hours", 5.0),
        row("costless", "bins_total", 500),
    ]
    ranked = scoring.value_against_cost(rows, "bins_total")
    assert [r.variant for r in ranked] == ["costless", "measured"]
    assert ranked[0].cost == 0.0
    assert "unmeasured" in scoring.format_scores(ranked)


def test_the_table_shows_both_sides_and_marks_the_incumbent():
    rows = [row("a", "bins_total", 8), row("a", "cpu_hours", 2.0)]
    text = scoring.format_scores(scoring.value_against_cost(rows, "bins_total", defaults=["a"]))
    assert "current default" in text
    assert "cpu_hours" in text
    assert "value" in text


def test_the_table_does_not_offer_a_ratio_to_rank_by():
    """The footer has to say what the ranking is, because a reader may remember the old one."""
    rows = [row("a", "bins_total", 8), row("a", "cpu_hours", 2.0)]
    text = scoring.format_scores(scoring.value_against_cost(rows, "bins_total"))
    assert "never as a denominator" in text
    assert "per unit" not in text


def test_the_table_says_so_when_there_is_nothing_to_rank():
    assert "nothing to compare" in scoring.format_scores([])


def test_cost_can_be_reported_as_wall_time_instead():
    """wall_seconds rewards a tool for being given more cores, which is why cpu_hours is the default."""
    rows = [row("a", "bins_total", 10), row("a", "wall_seconds", 100.0)]
    # Under the default cost metric this variant recorded no cost at all, and says so.
    assert scoring.value_against_cost(rows, "bins_total")[0].cost == 0.0
    ranked = scoring.value_against_cost(rows, "bins_total", cost_metric="wall_seconds")
    assert ranked[0].cost == 100.0


# ── the runner ────────────────────────────────────────────────────────────────────────────────


FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks")


def _tiny():
    return registry.get_dataset("tiny", os.path.join(FIXTURES, "datasets.toml"))


def _comparison(name):
    return registry.get_comparison(name, os.path.join(FIXTURES, "configurations.toml"))


def test_a_binning_variant_is_fed_the_gold_standard_assembly(tmp_path):
    """CAMI ships one so a binner can be measured without an assembler's variance in front of it."""
    dataset = _tiny()
    comparison = _comparison("binner_set")
    inputs = runner.gather_inputs(
        dataset, comparison, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
    )
    assert inputs.assembly is not None
    assert inputs.assembly.endswith("assembly.fa")
    argv = runner.build_argv(comparison, comparison.variant("metabat2"), inputs, "out", 8, 64)
    assert argv[1:4] == ["-m", "metawrap2.cli", "binning"]
    assert argv[argv.index("-a") + 1] == inputs.assembly
    assert argv[argv.index("-t") + 1] == "8"
    assert len([a for a in argv if a.endswith(".fq.gz")]) == 2


def test_the_global_force_flag_reaches_every_variant(tmp_path, monkeypatch):
    """A comparison that failed once could not be re-run: each variant is its own process, and the
    flags given to ``metawrap2`` applied to the harness and to nothing it started. Every variant
    stopped at "already contains a MetaWrap2 run" and the way forward was deleting directories."""
    from metawrap2 import command as shared

    dataset = _tiny()
    comparison = _comparison("binner_set")
    inputs = runner.gather_inputs(
        dataset, comparison, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
    )
    monkeypatch.setattr(shared.runner, "force", True)
    monkeypatch.setattr(shared.runner, "resume", False)
    assert runner._inherited_global_flags() == ["--force"]

    argv = runner.build_argv(
        comparison,
        comparison.variant("metabat2"),
        inputs,
        "out",
        8,
        64,
        global_flags=runner._inherited_global_flags(),
    )
    # Before the module name, which is where a global flag has to be for the CLI to strip it.
    assert argv[3] == "--force"
    assert argv[4] == "binning"


def test_an_assembly_variant_gets_comma_separated_mate_lists_and_a_memory_grant(tmp_path):
    dataset = _tiny()
    comparison = _comparison("assembler")
    inputs = runner.gather_inputs(
        dataset, comparison, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
    )
    argv = runner.build_argv(comparison, comparison.variant("megahit"), inputs, "out", 8, 64)
    # Sliced past `python -m metawrap2.cli`, whose own -m would otherwise be found first.
    flags = argv[3:]
    assert flags[flags.index("-m") + 1] == "64"
    assert flags[flags.index("-1") + 1].endswith("_1.fq.gz")
    assert flags[flags.index("-2") + 1].endswith("_2.fq.gz")


def test_a_refinement_comparison_without_bins_says_what_to_run_first(tmp_path):
    with pytest.raises(runner.RunnerError, match="--bins-from"):
        runner.gather_inputs(
            _tiny(),
            _comparison("refinement_completeness"),
            str(tmp_path / "data"),
            str(tmp_path / "work"),
            say=lambda _m: None,
        )


def test_a_module_the_runner_cannot_feed_is_refused_rather_than_guessed_at(tmp_path):
    """Numbers from the wrong command are worse than an error."""
    comparison = registry.Comparison(
        name="x",
        module="kraken2",
        variants=(registry.Variant(name="a"), registry.Variant(name="b")),
    )
    with pytest.raises(runner.RunnerError, match="cannot feed"):
        runner.gather_inputs(
            _tiny(), comparison, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
        )


def test_missing_data_points_at_the_fetch_command(tmp_path):
    dataset = registry.Dataset(
        name="absent",
        title="Absent",
        base_url="https://example.invalid/",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )
    with pytest.raises(runner.RunnerError, match="benchmark fetch absent"):
        runner.prepare_reads(dataset, str(tmp_path / "nope"), str(tmp_path), say=lambda _m: None)


def test_interleaved_reads_are_split_once_and_then_reused(tmp_path):
    """Splitting twice would put the second variant's reads on a colder cache than the first's."""
    data = tmp_path / "data" / "extracted" / "2017.12.04_18.45.54_sample_7" / "reads"
    data.mkdir(parents=True)
    record = "@r%d/%d\nACGT\n+\nIIII\n"
    (data / "anonymous_reads.fq").write_text(
        "".join(record % (n, mate) for n in range(3) for mate in (1, 2))
    )
    dataset = registry.Dataset(
        name="il",
        title="Interleaved",
        base_url="https://example.invalid/",
        reads_layout="interleaved",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )
    said = []
    work = str(tmp_path / "work")
    pairs = runner.prepare_reads(dataset, str(tmp_path / "data"), work, say=said.append)
    assert [p[0] for p in pairs] == ["sample_7"]
    assert os.path.isfile(pairs[0][1]) and os.path.isfile(pairs[0][2])
    assert any("splitting" in line for line in said)

    said.clear()
    runner.prepare_reads(dataset, str(tmp_path / "data"), work, say=said.append)
    assert any("already split" in line for line in said)
    assert not any("splitting" in line for line in said)


def test_the_sample_label_drops_the_publication_timestamp():
    """Otherwise every results row carries a 2017 timestamp and no two body sites sort together."""
    assert (
        runner._sample_label("/d/2017.12.04_18.45.54_sample_12/reads/anonymous_reads.fq.gz")
        == "sample_12"
    )


def test_cost_rows_record_cpu_hours_wall_time_and_peak_memory():
    dataset = _tiny()
    comparison = _comparison("binner_set")
    rows = runner.cost_rows(
        dataset,
        comparison,
        comparison.variant("metabat2"),
        {
            "wall_seconds": 1800.0,
            "user_seconds": 3400.0,
            "system_seconds": 200.0,
            "peak_rss_bytes": 8 * 1024**3,
        },
    )
    by_metric = {r.metric: r.value for r in rows}
    assert by_metric["wall_seconds"] == 1800.0
    assert by_metric["cpu_hours"] == pytest.approx(1.0)
    assert by_metric["peak_rss_bytes"] == 8 * 1024**3


def test_an_unmeasurable_platform_still_records_wall_time():
    """usage.run_measured degrades to wall time alone; that must not lose the whole row."""
    dataset = _tiny()
    comparison = _comparison("binner_set")
    rows = runner.cost_rows(
        dataset, comparison, comparison.variant("metabat2"), {"wall_seconds": 5.0}
    )
    assert [r.metric for r in rows] == ["wall_seconds"]


def test_a_dry_run_records_no_run_file_and_spawns_nothing(tmp_path, monkeypatch):
    def explode(*_a, **_k):
        raise AssertionError("ran a variant during a dry run")

    monkeypatch.setattr(runner._usage, "run_measured", explode)
    rows = runner.run_comparison(
        _tiny(),
        _comparison("binner_set"),
        str(tmp_path / "data"),
        str(tmp_path / "work"),
        dry_run=True,
        say=lambda _m: None,
    )
    assert rows == []
    assert not os.path.exists(tmp_path / "work" / "binner_set.run.json")


def test_a_failing_variant_is_recorded_and_the_others_still_run(tmp_path, monkeypatch):
    """One variant failing is a result about that variant; the rest are still worth having."""
    calls = []

    def fake(argv, **_kw):
        calls.append(argv)
        failed = "--maxbin2" in argv
        return (1 if failed else 0), {
            "wall_seconds": 10.0,
            "user_seconds": 20.0,
            "system_seconds": 0.0,
            "peak_rss_bytes": 1024,
        }

    monkeypatch.setattr(runner._usage, "run_measured", fake)
    work = tmp_path / "work"
    rows = runner.run_comparison(
        _tiny(),
        _comparison("binner_set"),
        str(tmp_path / "data"),
        str(work),
        variants=["metabat2", "maxbin2"],
        say=lambda _m: None,
    )
    assert len(calls) == 2
    variants = {r.variant for r in rows}
    assert "metabat2" in variants
    # The failed variant contributes no cost rows, so it cannot be ranked as cheap-and-good.
    assert "maxbin2" not in variants

    with open(work / "binner_set.run.json") as fh:
        record = json.load(fh)
    assert record["variants"]["maxbin2"]["returncode"] == 1
    assert record["variants"]["metabat2"]["returncode"] == 0
    # The grant is recorded next to the results, so a table can be read against what it was given.
    assert "threads" in record["budget"]


def test_every_variant_of_a_comparison_gets_the_same_grant(tmp_path, monkeypatch):
    """Two variants given different budgets are not a comparison."""
    grants = []

    def fake(argv, **_kw):
        grants.append(argv[argv.index("-t") + 1])
        return 0, {"wall_seconds": 1.0}

    monkeypatch.setattr(runner._usage, "run_measured", fake)
    runner.run_comparison(
        _tiny(),
        _comparison("binner_set"),
        str(tmp_path / "data"),
        str(tmp_path / "work"),
        variants=["metabat2", "concoct"],
        say=lambda _m: None,
    )
    assert len(set(grants)) == 1


# ── per-sample modules and shared flags ───────────────────────────────────────────────────────


def _two_sample_reads(tmp_path):
    """A paired two-sample dataset on disk, which the tiny committed fixture is too small to be."""
    data = tmp_path / "data" / "extracted"
    data.mkdir(parents=True)
    for sample in ("a", "b"):
        for mate in (1, 2):
            (data / ("%s_%d.fastq" % (sample, mate))).write_text("@r\nACGT\n+\nIIII\n")
    return registry.Dataset(
        name="pair",
        title="Two samples",
        base_url="https://example.invalid/",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )


def test_a_read_qc_variant_runs_once_per_sample_and_labels_its_cost(tmp_path, monkeypatch):
    """A per-sample module yields one measured run per sample, summing into one ratio."""
    dataset = _two_sample_reads(tmp_path)
    comparison = _comparison("trimmer")
    seen = []

    def fake(argv, **_kw):
        seen.append(argv)
        return 0, {"wall_seconds": 1.0, "user_seconds": 1800.0, "system_seconds": 0.0}

    monkeypatch.setattr(runner._usage, "run_measured", fake)
    work = tmp_path / "work"
    rows = runner.run_comparison(
        dataset,
        comparison,
        str(tmp_path / "data"),
        str(work),
        variants=["fastp"],
        say=lambda _m: None,
    )
    assert len(seen) == 2
    # One sample per invocation, each with its own -1/-2 and its own output directory.
    for argv in seen:
        assert argv.count("-1") == 1 and "," not in argv[argv.index("-1") + 1]
    assert {r.sample for r in rows} == {"a", "b"}
    assert (work / "trimmer" / "fastp" / "a").is_dir()
    assert (work / "trimmer" / "fastp" / "b").is_dir()
    # Half an hour of CPU per sample, summed across both, is one CPU-hour for the variant.
    (scored,) = scoring.value_against_cost(
        rows + [scoring.Row("pair", "trimmer", "fastp", "reads_kept", 500.0, sample="a")],
        "reads_kept",
    )
    assert scored.cost == pytest.approx(1.0)
    assert scored.value == pytest.approx(500.0)


def test_one_sample_failing_discards_the_whole_variant(tmp_path, monkeypatch):
    """A total summed over one sample is not comparable with one summed over two."""
    dataset = _two_sample_reads(tmp_path)
    comparison = _comparison("trimmer")
    calls = []

    def fake(argv, **_kw):
        calls.append(argv)
        return (1 if len(calls) == 1 else 0), {"wall_seconds": 1.0}

    monkeypatch.setattr(runner._usage, "run_measured", fake)
    rows = runner.run_comparison(
        dataset,
        comparison,
        str(tmp_path / "data"),
        str(tmp_path / "work"),
        variants=["fastp"],
        say=lambda _m: None,
    )
    assert rows == []
    # And the remaining samples are not run: they could not contribute a comparable number.
    assert len(calls) == 1


def test_flags_that_are_not_under_test_are_declared_once_for_every_variant(tmp_path):
    """A flag copy-pasted into three variants and forgotten in the fourth is a fake difference."""
    dataset = _two_sample_reads(tmp_path)
    comparison = _comparison("trimmer")
    assert "--skip-bmtagger" in comparison.common_args
    inputs = runner.gather_inputs(
        dataset, comparison, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
    )
    for variant in comparison.variants:
        one = runner.Inputs(pairs=[inputs.pairs[0]])
        argv = runner.build_argv(comparison, variant, one, "out", 8, 64)
        assert "--skip-bmtagger" in argv
        # ... and the thing that *is* under test still differs.
        assert variant.args[-1] in argv


def test_read_qc_refuses_a_whole_cohort_rather_than_silently_scoring_one_sample():
    comparison = _comparison("trimmer")
    inputs = runner.Inputs(pairs=[("a", "a_1.fq", "a_2.fq"), ("b", "b_1.fq", "b_2.fq")])
    with pytest.raises(runner.RunnerError, match="one sample at a time"):
        runner.build_argv(comparison, comparison.variant("fastp"), inputs, "out", 8, 64)


# ── read statistics, so two trimmers answer the same question ─────────────────────────────────


def test_read_stats_counts_reads_bases_and_q30_off_the_fastq_itself(tmp_path):
    """Scored from the surviving reads, not from each trimmer's own self-assessment."""
    path = tmp_path / "final_pure_reads_1.fastq"
    # Second read's qualities are '!' (Q0), so only the first read's four bases clear Q30.
    path.write_text("@r1\nACGTA\n+\nIIIII\n@r2\nACGT\n+\n!!!!\n")
    stats = scoring.read_stats(str(path))
    assert stats["reads_kept"] == 2
    assert stats["bases_kept"] == 9
    assert stats["q30_bases_kept"] == 5


def test_read_stats_reads_gzipped_output_too(tmp_path):
    import gzip

    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")
    assert scoring.read_stats(str(path))["reads_kept"] == 1


def test_scoring_the_same_run_twice_does_not_double_its_numbers(tmp_path):
    """The results TSV is appended to, so the same measurement can land in it more than once."""
    rows = [
        scoring.Row("d", "c", "v", "bins_total", 10.0, source="quality_report.tsv"),
        scoring.Row("d", "c", "v", "cpu_hours", 2.0, source="usage.run_measured"),
    ]
    once = scoring.value_against_cost(rows, "bins_total")
    twice = scoring.value_against_cost(rows + rows, "bins_total")
    assert [r.value for r in twice] == [r.value for r in once] == [10.0]
    assert [r.cost for r in twice] == [r.cost for r in once] == [2.0]


# ── which files count as the reads ────────────────────────────────────────────────────────────


def _interleaved(tmp_path, *, glob_pattern=""):
    """A CAMI2-shaped tree: one real interleaved FASTQ per sample, plus a small decoy beside it."""
    data = tmp_path / "data" / "extracted"
    for n in (0, 1):
        reads = data / ("2017.12.04_18.45.54_sample_%d" % n) / "reads"
        reads.mkdir(parents=True)
        (reads / "anonymous_reads.fq").write_text("@r/1\nACGT\n+\nIIII\n@r/2\nTGCA\n+\nIIII\n")
        # The debugging extract CAMI2 ships next to the real file, under the same sample_N
        # directory and therefore under the same sample label.
        (reads / "subsample_0.fq").write_text("@s/1\nAA\n+\nII\n@s/2\nTT\n+\nII\n")
    return registry.Dataset(
        name="cami_shaped",
        title="Interleaved, with a decoy",
        base_url="https://example.invalid/",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
        reads_layout="interleaved",
        reads_glob=glob_pattern,
    )


def test_reads_glob_picks_the_real_file_over_the_decoy_beside_it(tmp_path):
    dataset = _interleaved(tmp_path, glob_pattern="*/reads/anonymous_reads.fq*")
    found = runner.read_candidates(dataset, str(tmp_path / "data" / "extracted"))
    assert [os.path.basename(p) for p in found] == ["anonymous_reads.fq"] * 2


def test_two_files_under_one_sample_label_is_an_error_not_a_silent_overwrite(tmp_path):
    # Without a reads_glob the suffix scan finds both files, and both resolve to "sample_0".
    # Splitting them in turn would leave the smaller one in place, and a trimmer measured on
    # 3 MB of a 10 GB sample is a number nobody can tell is wrong.
    dataset = _interleaved(tmp_path)
    with pytest.raises(runner.RunnerError) as exc:
        runner.prepare_reads(
            dataset, str(tmp_path / "data"), str(tmp_path / "work"), say=lambda _m: None
        )
    assert "resolve to sample" in str(exc.value)
    assert "reads_glob" in str(exc.value)


def test_a_reads_glob_that_matches_nothing_says_so(tmp_path):
    dataset = _interleaved(tmp_path, glob_pattern="*/reads/not_a_real_name.fq")
    with pytest.raises(runner.RunnerError) as exc:
        runner.read_candidates(dataset, str(tmp_path / "data" / "extracted"))
    assert "matches nothing" in str(exc.value)


def test_a_dry_run_names_the_split_reads_without_writing_them(tmp_path):
    """A dry run on CAMI2 must not spend an hour writing 100 GB before printing six commands."""
    dataset = _interleaved(tmp_path, glob_pattern="*/reads/anonymous_reads.fq*")
    work = tmp_path / "work"
    pairs = runner.prepare_reads(
        dataset, str(tmp_path / "data"), str(work), dry_run=True, say=lambda _m: None
    )
    assert [s for s, _1, _2 in pairs] == ["sample_0", "sample_1"]
    assert not (work / "reads").exists()

    # The same call without dry_run does write them, to exactly the paths it just named.
    real = runner.prepare_reads(dataset, str(tmp_path / "data"), str(work), say=lambda _m: None)
    assert real == pairs
    assert all(os.path.getsize(p) for _s, r1, r2 in real for p in (r1, r2))


@pytest.mark.conda_envs
def test_a_comparison_refuses_to_start_without_an_env_manager(monkeypatch, tmp_path):
    """Twenty failed runs and an empty table is a worse answer than one sentence.

    The runner deliberately records a failed variant and carries on, because one tool failing is a
    result about that tool. An absent mamba is not that: it fails every variant of every comparison
    identically, and only after each has been given its turn.
    """
    from metawrap2.benchmark import runner as runner_mod

    monkeypatch.setattr(runner_mod.shutil, "which", lambda name: None)
    with pytest.raises(runner_mod.RunnerError) as caught:
        runner_mod.require_env_manager()
    assert "not on PATH" in str(caught.value)

    monkeypatch.setattr(runner_mod.shutil, "which", lambda name: "/opt/conda/bin/" + name)
    runner_mod.require_env_manager()  # present: no opinion


# ── scoping a comparison to fewer samples, and co-assembling the ones chosen ───────────────────


def _three_sample_reads(tmp_path):
    """Three paired samples on disk, enough to select a subset of."""
    data = tmp_path / "data" / "extracted"
    data.mkdir(parents=True)
    for sample in ("a", "b", "c"):
        for mate in (1, 2):
            (data / ("%s_%d.fastq" % (sample, mate))).write_text(
                "@%s/%d\nACGT\n+\nIIII\n" % (sample, mate)
            )
    return registry.Dataset(
        name="three",
        title="Three samples",
        base_url="https://example.invalid/",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )


def _inputs_for(tmp_path):
    dataset = _three_sample_reads(tmp_path)
    return runner.gather_inputs(
        dataset,
        _comparison("trimmer"),
        str(tmp_path / "data"),
        str(tmp_path / "work"),
        say=lambda _m: None,
    )


def test_max_samples_takes_a_deterministic_prefix_not_an_arbitrary_one(tmp_path):
    """Two variants must get the same subset, so the subset cannot depend on directory order."""
    inputs = _inputs_for(tmp_path)
    assert [p[0] for p in runner.select_samples(inputs, max_samples=2).pairs] == ["a", "b"]
    # And asking for more than there are is not an error - it is every sample.
    assert len(runner.select_samples(inputs, max_samples=99).pairs) == 3


def test_naming_a_sample_that_is_not_there_says_what_is(tmp_path):
    """A typo'd label would otherwise silently change what the comparison measured."""
    inputs = _inputs_for(tmp_path)
    with pytest.raises(runner.RunnerError) as caught:
        runner.select_samples(inputs, samples=["a", "sample_typo"])
    message = str(caught.value)
    assert "sample_typo" in message
    assert "a, b, c" in message  # tells you what you could have meant


def test_selecting_samples_by_name_keeps_only_those(tmp_path):
    inputs = _inputs_for(tmp_path)
    assert [p[0] for p in runner.select_samples(inputs, samples=["c", "a"]).pairs] == ["a", "c"]


def test_a_cohort_assembly_is_fed_one_concatenated_pair_not_a_comma_list(tmp_path):
    """`assembly` takes one file per mate; a comma-separated list made it exit immediately.

    This is the regression test for that: the runner built ``-1 a_1.fastq,b_1.fastq``, which the
    module reported as "does not exist" in under a second, so both variants of the comparison
    failed identically and the results table came out empty.
    """
    inputs = _inputs_for(tmp_path)
    work = tmp_path / "work"
    joined = runner.prepare_coassembly(inputs, str(work), say=lambda _m: None)

    assert len(joined.pairs) == 1
    label, r1, r2 = joined.pairs[0]
    assert label == "coassembly"
    # Every selected sample is in there, and the mates stay separated.
    for path, mate in ((r1, 1), (r2, 2)):
        with open(path) as fh:
            expected = "".join("@%s/%d\nACGT\n+\nIIII\n" % (s, mate) for s in ("a", "b", "c"))
            assert fh.read() == expected

    argv = runner.build_argv(
        _comparison("assembler"), _comparison("assembler").variant("megahit"), joined, "out", 8, 64
    )
    assert argv[argv.index("-1") + 1] == r1
    assert "," not in argv[argv.index("-1") + 1]

    # More than one pair reaching build_argv is the bug itself, so it is refused rather than joined.
    with pytest.raises(runner.RunnerError):
        runner.build_argv(
            _comparison("assembler"),
            _comparison("assembler").variant("megahit"),
            inputs,
            "out",
            8,
            64,
        )


def test_the_concatenation_is_built_once_and_reused_by_the_second_variant(tmp_path):
    """Rebuilding it per variant would hand the second one a colder page cache than the first."""
    inputs = _inputs_for(tmp_path)
    work = tmp_path / "work"
    first = runner.prepare_coassembly(inputs, str(work), say=lambda _m: None)
    stamp = os.path.getmtime(first.pairs[0][1])
    said = []
    again = runner.prepare_coassembly(inputs, str(work), say=said.append)
    assert again.pairs[0][1] == first.pairs[0][1]
    assert os.path.getmtime(first.pairs[0][1]) == stamp  # not rewritten
    assert any("already current" in line for line in said)


def test_one_sample_needs_no_concatenation(tmp_path):
    """Copying a 10 GB file onto itself to make a "cohort" of one would be pure waste."""
    inputs = runner.select_samples(_inputs_for(tmp_path), samples=["b"])
    same = runner.prepare_coassembly(inputs, str(tmp_path / "work"), say=lambda _m: None)
    assert same.pairs == inputs.pairs
    assert not os.path.exists(tmp_path / "work" / "coassembly_1.fastq")


def test_a_dry_run_names_the_concatenation_without_writing_it(tmp_path):
    inputs = _inputs_for(tmp_path)
    said = []
    runner.prepare_coassembly(inputs, str(tmp_path / "work"), dry_run=True, say=said.append)
    assert any("would concatenate" in line for line in said)
    assert not os.path.exists(tmp_path / "work" / "coassembly_1.fastq")


def test_the_run_record_names_the_samples_measured_not_the_concatenation(tmp_path, monkeypatch):
    """A cohort module's rows carry no sample label, so the record is the only place the scope is.

    Without it a three-sample co-assembly and a ten-sample one are indistinguishable in the results
    TSV - and `--max-samples` makes that a real possibility rather than a hypothetical one.
    """
    monkeypatch.setattr(
        runner._usage, "run_measured", lambda argv, **_kw: (0, {"wall_seconds": 1.0})
    )
    work = tmp_path / "work"
    runner.run_comparison(
        _three_sample_reads(tmp_path),
        _comparison("assembler"),
        str(tmp_path / "data"),
        str(work),
        variants=["megahit"],
        max_samples=2,
        say=lambda _m: None,
    )
    with open(work / "assembler.run.json") as fh:
        record = json.load(fh)
    assert record["samples"] == ["a", "b"]


def test_ambers_own_column_names_are_the_ones_read(tmp_path):
    """The prose headings ("F1 score (bp)") are on AMBER's HTML page; its TSV is snake_case.

    Reading only the prose spelling left a results.tsv that parsed to nothing at all, which the
    scorer reported as "nothing to score" - a real run of a working tool, recorded as absent.
    """
    path = tmp_path / "results.tsv"
    path.write_text(
        "Sample\tTool\taccuracy_bp\tadjusted_rand_index_bp\tf1_score_bp\tf1_score_bp_cami1\t"
        "percentage_of_assigned_bps\tprecision_avg_bp\trecall_avg_bp\n"
        "tiny\tconcoct_bins\t0.61\t0.316\t0.861\t0.5\t1.0\t0.876\t0.847\n"
    )
    metrics = scoring.amber_metrics_by_tool(str(path))["concoct_bins"]
    assert metrics == {
        "amber_accuracy_bp": 0.61,
        "amber_ari_bp": 0.316,
        "amber_f1_bp": 0.861,
        "amber_assigned_bp_fraction": 1.0,
        "amber_purity_bp": 0.876,
        "amber_completeness_bp": 0.847,
    }, "CAMI 1's F1 is a different definition and must not land in the same column as this one"


# ── quant_bins in the harness ──────────────────────────────────────────────────────────────────


def test_a_quantifier_comparison_without_bins_says_what_to_run_first(tmp_path):
    """There is nothing to quantify without bins, and inventing some would measure the wrong thing."""
    with pytest.raises(runner.RunnerError, match="--bins-from"):
        runner.gather_inputs(
            _tiny(),
            _comparison("quantifier"),
            str(tmp_path / "data"),
            str(tmp_path / "work"),
            say=lambda _m: None,
        )


def test_a_quantifier_variant_gets_the_bins_and_every_sample_but_no_assembly(tmp_path):
    """All three quantifiers work against the bins. Handing salmon alone a whole-assembly
    pseudo-transcriptome would compare three tools against two different references."""
    bins = tmp_path / "bins" / "metawrap_70_10_bins"
    bins.mkdir(parents=True)
    (bins / "bin.1.fa").write_text(">c1\nACGT\n")
    comparison = _comparison("quantifier")
    inputs = runner.gather_inputs(
        _tiny(),
        comparison,
        str(tmp_path / "data"),
        str(tmp_path / "work"),
        bins_from=str(tmp_path / "bins"),
        say=lambda _m: None,
    )
    argv = runner.build_argv(comparison, comparison.variant("sylph"), inputs, "out", 8, 64)
    assert argv[1:4] == ["-m", "metawrap2.cli", "quant_bins"]
    assert argv[argv.index("-b") + 1] == str(bins)
    assert "-a" not in argv
    assert len([a for a in argv if a.endswith(".fq.gz")]) == 2


def test_one_bin_set_prefers_the_consolidated_set_over_the_raw_ones(tmp_path):
    """Quantifying the raw sets as well would count the same genome several times over."""
    root = tmp_path / "refined"
    for name in ("metawrap_70_10_bins", "metabat2_bins", "concoct_bins"):
        (root / name).mkdir(parents=True)
        (root / name / "bin.1.fa").write_text(">c\nACGT\n")
    assert runner.one_bin_set(str(root)) == str(root / "metawrap_70_10_bins")


def test_one_bin_set_takes_a_directory_of_bins_as_it_stands(tmp_path):
    root = tmp_path / "bins"
    root.mkdir()
    (root / "bin.1.fa").write_text(">c\nACGT\n")
    assert runner.one_bin_set(str(root)) == str(root)


def test_one_bin_set_refuses_to_choose_between_overlapping_raw_sets(tmp_path):
    """Picking one arbitrarily would attribute a binner's numbers to whichever sorted first."""
    root = tmp_path / "raw"
    for name in ("metabat2_bins", "concoct_bins"):
        (root / name).mkdir(parents=True)
        (root / name / "bin.1.fa").write_text(">c\nACGT\n")
    with pytest.raises(runner.RunnerError):
        runner.one_bin_set(str(root))


def test_the_abundance_truth_files_are_keyed_by_the_same_sample_label_as_the_reads(tmp_path):
    """The table's columns carry that label, and the join between estimate and truth is the label."""
    extracted = tmp_path / "data" / "extracted"
    for stamp, number in (("2017.12.04_18.45.54", 0), ("2017.12.04_19.01.02", 10)):
        reads = extracted / ("%s_sample_%d" % (stamp, number)) / "reads"
        reads.mkdir(parents=True)
        (reads / "reads_mapping.tsv.gz").write_bytes(b"")
    dataset = registry.Dataset(
        name="cami",
        title="CAMI-shaped",
        base_url="https://example.invalid/",
        abundance_truth_glob="*/reads/reads_mapping.tsv.gz",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )
    found = runner.abundance_truth_files(dataset, str(tmp_path / "data"))
    assert sorted(found) == ["sample_0", "sample_10"]


def test_a_dataset_without_an_abundance_truth_says_so_rather_than_scoring_nothing(tmp_path):
    dataset = registry.Dataset(
        name="bare",
        title="Bare",
        base_url="https://example.invalid/",
        files=(registry.Remote(path="a.tar.gz", md5="0" * 32, size=1),),
    )
    with pytest.raises(runner.RunnerError, match="abundance_truth_glob"):
        runner.abundance_truth_files(dataset, str(tmp_path / "data"))


def test_cami2_declares_where_its_per_sample_read_truth_lives():
    """Without this key the quantifier comparison cannot be scored on the one dataset that has it."""
    dataset = registry.get_dataset("cami2_gi", os.path.join(FIXTURES, "datasets.toml"))
    assert dataset.abundance_truth_glob == "*/reads/reads_mapping.tsv.gz"


def test_a_rescore_that_reads_a_different_file_replaces_the_old_number(tmp_path):
    """A re-score can change *which* bin set is the variant's answer, not just re-read the same one.

    Fixing bin_dirs so a consolidated set one level down is found did exactly that, and because the
    dedup key included the source file the stale row sat beside the new one and the two were summed:
    COMEBin's 178 high-quality MAGs read as 356. Same dataset/comparison/variant/sample/metric is one
    measurement whatever file it came from.
    """
    rows = [
        scoring.Row("d", "c", "v", "bins_high_quality", 178.0, source="comebin_bins.stats.tsv"),
        scoring.Row("d", "c", "v", "cpu_hours", 39.6, source="usage.run_measured"),
        # The same arm, re-scored later on its consolidated set.
        scoring.Row(
            "d", "c", "v", "bins_high_quality", 183.0, source="metawrap_70_10_bins.stats.tsv"
        ),
    ]
    (scored,) = scoring.value_against_cost(rows, "bins_high_quality")
    assert scored.value == 183.0, "the later row is the answer, not the sum"
    assert scored.cost == pytest.approx(39.6)


def test_per_bin_set_numbers_still_coexist_because_they_are_different_metrics(tmp_path):
    """The suffix is what keeps two bin sets' figures apart; dedup must not need the source for that."""
    rows = [
        scoring.Row(
            "d", "c", "v", "bins_high_quality", 183.0, source="metawrap_70_10_bins.stats.tsv"
        ),
        scoring.Row("d", "c", "v", "bins_high_quality__metabat2_bins", 37.0, source="a.stats.tsv"),
        scoring.Row("d", "c", "v", "cpu_hours", 52.8),
    ]
    (scored,) = scoring.value_against_cost(rows, "bins_high_quality")
    assert scored.value == 183.0


def test_a_failed_variant_is_not_scored_from_its_partial_output(tmp_path, capsys):
    """The runner records a non-zero exit; scoring never read it.

    A variant that died after its first binner wrote a directory still has bins on disk, so CheckM
    and AMBER scored them and the number was written to the committed TSV - and because the runner
    discards a failed variant's cost rows, the arm appeared with a quality value and cost
    "unmeasured", which is also what a finished arm whose rusage was lost looks like.
    """
    import json

    from metawrap2.commands import benchmark as cmd

    work = tmp_path / "work"
    (work / "binner_set" / "metabat2").mkdir(parents=True)
    (work / "binner_set" / "metabat2" / "metabat2_bins").mkdir()
    (work / "binner_set.run.json").write_text(
        json.dumps({"variants": {"metabat2": {"returncode": 1, "runs": []}}})
    )
    rc = cmd.main(
        [
            "score",
            "tiny",
            "binner_set",
            "--work",
            str(work),
            "--results",
            str(tmp_path / "r.tsv"),
        ]
    )
    out = capsys.readouterr().out
    assert "Not scored" in out and "metabat2" in out
    assert "--score-failed" in out
    # Nothing was scored, so there is nothing to rank and no results file to write.
    assert rc == 1
    assert not (tmp_path / "r.tsv").exists()
