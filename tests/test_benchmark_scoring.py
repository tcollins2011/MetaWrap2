"""Tests for benchmark scoring and running: the tidy table, the parsers, and the ranking.

Checks that cost is aggregated the same way as accuracy, that a variant with no measured cost is
not ranked first, and that the current default is marked in the table.
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


# ── the parsers ───────────────────────────────────────────────────────────────────────────────


CHECKM2 = "Name\tCompleteness\tContamination\nbin_001\t95.5\t1.2\nbin_002\t72.0\t8.0\nbin_003\t40.0\t2.0\n"
CHECKM1 = "Bin Id\tMarker lineage\tCompleteness\tContamination\nbin_001\tk__Bacteria\t95.5\t1.2\n"


def test_quality_counts_use_the_community_thresholds(tmp_path):
    """Both CheckM1 and CheckM2 column layouts are read. Bins with NA scores are skipped, so a
    bin CheckM could not assess is not counted as a poor one."""
    path = tmp_path / "quality_report.tsv"
    path.write_text(CHECKM2 + "bin_004\tNA\tNA\n")
    counts = scoring.quality_counts(scoring.bin_quality(str(path)))
    assert counts["bins_near_complete"] == 1  # 95.5 / 1.2
    assert counts["bins_high_quality"] == 2  # and 72.0 / 8.0
    assert counts["bins_medium_quality"] == 2  # 40.0 is below 50
    assert counts["bins_total"] == 3
    # A contaminated but complete bin is not near-complete.
    assert scoring.quality_counts([("b", 99.0, 40.0)])["bins_near_complete"] == 0

    v1 = tmp_path / "checkm.tsv"
    v1.write_text(CHECKM1)
    assert scoring.bin_quality(str(v1)) == [("bin_001", 95.5, 1.2)]


def test_ambers_own_column_names_are_the_ones_read(tmp_path):
    """AMBER's TSV uses snake_case column names; the prose headings are only on its HTML page.

    Reading the prose names parsed nothing, and the scorer reported "nothing to score".
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


# ── the ranking ───────────────────────────────────────────────────────────────────────────────


def test_a_variant_that_scores_higher_for_far_more_compute_wins_the_ranking():
    """Variants are ranked on quality alone; cost is reported but does not reorder the table.

    Ranking on quality per cost once made a cheaper variant the default, and it recovered fewer
    MAGs end to end than the tool it replaced. Cost trade-offs are expressed through `--preset`.
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

    # Counts are summed across samples, and so is cost: ten samples give one row per variant.
    rows = [
        row("v", "bins_near_complete", 4, sample="sample_0"),
        row("v", "bins_near_complete", 6, sample="sample_1"),
        row("v", "cpu_hours", 1.0, sample="sample_0"),
        row("v", "cpu_hours", 3.0, sample="sample_1"),
    ]
    ranked = scoring.value_against_cost(rows, "bins_near_complete")
    assert ranked[0].value == 10
    assert ranked[0].cost == 4.0

    # A fraction is averaged across samples, since a sum would depend on the sample count.
    # Cost is still a total over samples.
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


def test_a_rescore_that_reads_a_different_file_replaces_the_old_number(tmp_path):
    """A re-score from a different source file replaces the earlier row instead of adding to it.

    Rows are deduplicated on dataset/comparison/variant/sample/metric, not on the source file.
    When the source file was part of the key, COMEBin's 178 high-quality MAGs were summed to 356.
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


def test_a_failed_variant_is_not_scored_from_its_partial_output(tmp_path, capsys):
    """A variant with a non-zero exit in the run record is not scored.

    A failed variant can leave partial bins on disk. Scoring them gave a quality value with cost
    "unmeasured", which looked like a finished variant with missing rusage.
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


# ── the runner ────────────────────────────────────────────────────────────────────────────────


FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks")


def _tiny():
    return registry.get_dataset("tiny", os.path.join(FIXTURES, "datasets.toml"))


def _comparison(name):
    return registry.get_comparison(name, os.path.join(FIXTURES, "configurations.toml"))


def test_a_failing_variant_is_recorded_and_the_others_still_run(tmp_path, monkeypatch):
    """A failing variant is recorded and the remaining variants still run."""
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
    # The failed variant contributes no cost rows, so it cannot rank as cheap.
    assert "maxbin2" not in variants

    with open(work / "binner_set.run.json") as fh:
        record = json.load(fh)
    assert record["variants"]["maxbin2"]["returncode"] == 1
    assert record["variants"]["metabat2"]["returncode"] == 0
    # The resource grant is recorded with the results.
    assert "threads" in record["budget"]


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


def test_a_cohort_assembly_is_fed_one_concatenated_pair_not_a_comma_list(tmp_path):
    """A co-assembly gets one concatenated file per mate, since `assembly` takes one file per mate.

    The runner used to pass ``-1 a_1.fastq,b_1.fastq``, which the module rejected as missing.
    """
    inputs = _inputs_for(tmp_path)
    work = tmp_path / "work"
    joined = runner.prepare_coassembly(inputs, str(work), say=lambda _m: None)

    assert len(joined.pairs) == 1
    label, r1, r2 = joined.pairs[0]
    assert label == "coassembly"
    # Every selected sample is included, and the mates stay separate.
    for path, mate in ((r1, 1), (r2, 2)):
        with open(path) as fh:
            expected = "".join("@%s/%d\nACGT\n+\nIIII\n" % (s, mate) for s in ("a", "b", "c"))
            assert fh.read() == expected

    argv = runner.build_argv(
        _comparison("assembler"), _comparison("assembler").variant("megahit"), joined, "out", 8, 64
    )
    assert argv[argv.index("-1") + 1] == r1
    assert "," not in argv[argv.index("-1") + 1]

    # build_argv refuses more than one pair instead of joining them.
    with pytest.raises(runner.RunnerError):
        runner.build_argv(
            _comparison("assembler"),
            _comparison("assembler").variant("megahit"),
            inputs,
            "out",
            8,
            64,
        )
