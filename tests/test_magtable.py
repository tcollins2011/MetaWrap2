"""Tests for the per-genome table that joins quality, taxonomy and coverage.

The two failures checked most: a genome lost because two files spell its name differently, and an
unmeasured cell written as a zero.
"""

from __future__ import annotations

import os
import shutil
from collections import Counter
from pathlib import Path

from metawrap2 import magtable as M

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "report")
#: Every module output directory in the fixture run. See fixtures/report/SOURCE.md.
ALL_PARTS = ("REFINEMENT", "CLASSIFY", "QUANT")


def _run(tmp_path, *parts: str) -> list:
    """A run directory holding only *parts* of the fixture run, so inputs can be taken away."""
    root = tmp_path / "run"
    root.mkdir(exist_ok=True)
    for part in parts:
        shutil.copytree(os.path.join(FIXTURES, part), root / part, dirs_exist_ok=True)
    return [str(root)]


def _table(tmp_path, *parts: str) -> M.Table:
    return M.build(M.discover(_run(tmp_path, *(parts or ALL_PARTS))))


def _by_genome(table: M.Table) -> dict:
    return {row[0]: dict(zip(table.columns, row)) for row in table.rows}


#: What `checkm.write_stats` writes above the rows.
STATS_HEADER = "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"


def _study(tmp_path, bins_per_sample: dict, gunc: dict | None = None) -> list:
    """A per-sample study, laid out as ``metawrap2 run`` writes one.

    One ``BIN_REFINEMENT/<sample>/`` per sample, each with its own consolidated bin set and
    ``work_files/gunc/`` report. Every sample numbers its bins from ``bin_001``.
    """
    root = tmp_path / "study" / "BIN_REFINEMENT"
    for sample, count in bins_per_sample.items():
        directory = root / sample
        directory.mkdir(parents=True)
        with open(directory / "metawrap_70_10_bins.stats.tsv", "w") as fh:
            fh.write(STATS_HEADER)
            fh.writelines(
                "bin_%03d\t95.0\t1.0\t0.41\tClostridiales\t40000\t2500000\n" % number
                for number in range(1, count + 1)
            )
        verdicts = (gunc or {}).get(sample)
        if verdicts is not None:
            report = directory / "work_files" / "gunc"
            report.mkdir(parents=True)
            with open(report / "GUNC.progenomes_2.1.maxCSS_level.tsv", "w") as fh:
                fh.write("genome\ttaxonomic_level\tclade_separation_score\tpass.GUNC\n")
                fh.writelines(
                    "%s\tspecies\t%.2f\t%s\n" % (name, css, "True" if css <= 0.45 else "False")
                    for name, css in verdicts.items()
                )
    return [str(tmp_path / "study")]


# ── the join ─────────────────────────────────────────────────────────────────────────────────


def test_a_genome_is_joined_across_three_spellings_of_its_name(tmp_path):
    """``bin_001`` in the stats table, ``bin_001.fasta`` in taxonomy, ``bin_001.fa`` in GUNC.

    One genome named three ways must join into one row.
    """
    row = _by_genome(_table(tmp_path))["bin_001"]
    assert row["completeness"] == "98.65"  # from the stats table
    assert row["taxonomy"].startswith("Bacillota;")  # from bin_taxonomy.tsv
    assert row["gunc_pass"] == "yes"  # from GUNC's own report
    assert row["abundance.sampleA"] == "143.178035"  # from the abundance matrix


def test_a_bin_named_with_a_dot_is_not_collapsed_to_bin(tmp_path):
    """``bin.2`` is a whole genome's name; ``os.path.splitext`` would make it ``bin``.

    metaWRAP-era bin sets are named this way; splitext would merge the whole set into one row.
    """
    rows = _by_genome(_table(tmp_path))
    assert "bin.2" in rows
    assert "bin" not in rows
    assert rows["bin.2"]["taxonomy"].startswith("Bacteroidota;")


def test_a_genome_only_in_the_abundance_matrix_still_gets_a_row(tmp_path):
    """A genome present in any input gets a row.

    ``bin_004`` was quantified but never scored or classified, as happens when quantification
    covered the pre-refinement set.
    """
    row = _by_genome(_table(tmp_path))["bin_004"]
    assert row["abundance.sampleB"] == "4.200000"
    assert row["completeness"] == M.BLANK


def test_the_consolidated_bin_set_is_preferred_over_a_binners_raw_table(tmp_path):
    """The consolidated stats table is chosen from the four in a bin_refinement directory.

    Alphabetical order would pick CONCOCT's. The decoy in the fixture names genomes no other
    input mentions.
    """
    inputs = M.discover(_run(tmp_path, "REFINEMENT"))
    (bin_set,) = inputs.bin_sets
    assert os.path.basename(bin_set.stats or "") == "metawrap_70_10_bins.stats.tsv"
    assert "concoct_42" not in _by_genome(M.build(inputs))


def test_the_unbinned_catch_all_is_not_a_genome(tmp_path):
    """concoct_bins/unbinned.fasta (leftover contigs, not a MAG) gets no row."""
    stats = os.path.join(FIXTURES, "REFINEMENT", "concoct_bins.stats.tsv")
    assert "unbinned" not in M.read_stats(stats)
    assert "concoct_42" in M.read_stats(stats)


# ── blank is not zero ────────────────────────────────────────────────────────────────────────


def test_a_measured_zero_survives_as_a_zero(tmp_path):
    """A measured contamination of 0 is written as 0."""
    row = _by_genome(_table(tmp_path))["bin.2"]
    assert row["contamination"] == "0.0"
    assert row["abundance.sampleA"] == "0.0"


def test_an_unscored_genome_has_blank_scores_and_no_tier(tmp_path):
    """A row CheckM wrote but could not score gets blank score cells, not zeros."""
    row = _by_genome(_table(tmp_path))["bin_003"]
    assert row["completeness"] == M.BLANK
    assert row["contamination"] == M.BLANK
    assert row["quality_tier"] == M.BLANK
    assert row["mimag_tag"] == M.BLANK


def test_an_uncounted_contig_count_is_blank_rather_than_zero(tmp_path):
    """A genome absent from the membership table has a blank contig count, not zero."""
    rows = _by_genome(_table(tmp_path))
    assert rows["bin_001"]["contigs"] == "2"
    assert rows["bin_003"]["contigs"] == M.BLANK


def test_a_bin_gunc_did_not_score_is_blank_rather_than_failed(tmp_path):
    """A bin GUNC did not score gets a blank gunc_pass, neither yes nor no."""
    rows = _by_genome(_table(tmp_path))
    assert rows["bin.2"]["gunc_pass"] == "no"  # GUNC did score this one, and it failed
    assert rows["bin_003"]["gunc_pass"] == M.BLANK
    assert rows["bin_003"]["gunc_css"] == M.BLANK


def test_a_classified_but_unplaced_genome_is_blank_not_missing(tmp_path):
    """classify_bins writes a row with no lineage for a bin it could not place.

    This differs from a bin never classified; the two are told apart by which inputs the table
    reports as absent.
    """
    table = _table(tmp_path)
    assert _by_genome(table)["bin_003"]["taxonomy"] == M.BLANK
    assert "taxonomy" not in table.absent  # the file was there; the cell was empty


# ── every input is optional ──────────────────────────────────────────────────────────────────


def test_a_run_that_skipped_classification_and_quantification_still_gets_a_table(tmp_path):
    """Every input is optional: a missing one leaves its columns blank and is named in ``absent``.

    This separates "not classified" from "classified as nothing".
    """
    table = _table(tmp_path, "REFINEMENT")
    assert [row[0] for row in table.rows] == ["bin_001", "bin.2", "bin_003"]
    assert table.abundance_samples == []
    assert sorted(table.absent) == ["abundance", "taxonomy"]
    assert all(row[table.columns.index("taxonomy")] == M.BLANK for row in table.rows)


def test_a_run_with_only_an_abundance_matrix_is_still_a_table(tmp_path):
    """With only quantification, there is still one row per quantified genome."""
    table = _table(tmp_path, "QUANT")
    assert [row[0] for row in table.rows] == ["bin.2", "bin_001", "bin_004"]
    assert sorted(table.absent) == ["contigs", "gunc", "stats", "taxonomy"]


def test_no_inputs_at_all_is_an_empty_table_rather_than_an_error(tmp_path):
    """Building from a run with none of the inputs is not an error."""
    table = M.build(M.discover([str(tmp_path)]))
    assert table.rows == []
    assert sorted(table.absent) == sorted(M.INPUT_KINDS)


def test_an_older_runs_suffixes_are_still_found(tmp_path):
    """The older ``.stats`` and ``.tab`` suffixes are still read.

    This repo's reference run uses them.
    """
    root = tmp_path / "old"
    (root / "REFINE").mkdir(parents=True)
    (root / "REFINE" / "metawrap_50_10_bins.stats").write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin.1\t91.0\t2.0\t0.40\tClostridiales\t50000\t2000000\n"
    )
    (root / "CLASSIFY").mkdir()
    (root / "CLASSIFY" / "bin_taxonomy.tab").write_text("bin.1.fa\tBacillota\n")
    table = M.build(M.discover([str(root)]))
    assert _by_genome(table)["bin.1"]["taxonomy"] == "Bacillota"


# ── the tier column and the file that is written ─────────────────────────────────────────────


def test_each_genome_gets_one_tier_rather_than_every_tier_it_clears(tmp_path):
    """Each genome gets one tier, unlike the benchmark's non-exclusive counts.

    Otherwise one genome would appear in several bars of the chart.
    """
    assert M.quality_tier(98.65, 1.21) == "near_complete"
    assert M.quality_tier(75.4, 0.0) == "high_quality"
    assert M.quality_tier(20.0, 1.0) == M.BELOW_TIERS
    assert M.quality_tier(99.0, None) == M.BLANK


def test_the_tier_chart_counts_every_genome_including_the_unmeasured_ones(tmp_path):
    """Unscored genomes are counted under NOT_MEASURED, not dropped."""
    table = _table(tmp_path)
    counts = M.tier_counts(table)
    assert sum(count for _label, count in counts) == len(table.rows)
    assert dict(counts)[M.NOT_MEASURED] == 2  # bin_003 was not scored, bin_004 not even seen


def test_the_written_tsv_is_the_table_with_its_header(tmp_path):
    """The written TSV has one header, one row per genome, and keeps blank cells."""
    table = _table(tmp_path)
    path = M.write(table, str(tmp_path / "out" / M.TABLE_NAME))
    lines = Path(path).read_text().splitlines()
    assert lines[0].split("\t") == table.columns
    assert len(lines) == len(table.rows) + 1
    # Blank cells must survive writing; a writer that skipped empty trailing cells would shorten
    # the row.
    assert lines[3].split("\t") == table.rows[2]
    assert lines[3].endswith("\t\t")


# ── a study is more than one bin set ─────────────────────────────────────────────────────────


def test_every_samples_bin_set_is_in_the_table_not_just_the_first(tmp_path):
    """Every sample's bin set in a per-sample study is read, not only the first.

    `find_stats` once returned one path, and a 62-genome study was reported as 24.
    """
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 24, "sample_7": 15, "sample_13": 23})))
    assert table.bin_sets == 3
    assert len(table.rows) == 62
    counts = Counter(row[0] for row in table.rows)
    assert counts == {"sample_6": 24, "sample_13": 23, "sample_7": 15}


def test_more_than_one_bin_set_carries_a_sample_column_in_front_of_the_genome(tmp_path):
    """The sample column, from the directory holding the stats table, comes before the genome."""
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 2, "sample_7": 2})))
    assert table.columns[:2] == [M.SAMPLE_COLUMN, "genome"]
    assert sorted({row[0] for row in table.rows}) == ["sample_6", "sample_7"]


def test_a_single_bin_set_gets_no_sample_column(tmp_path):
    """A single-bin-set run has no sample column."""
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 3})))
    assert table.bin_sets == 1
    assert M.SAMPLE_COLUMN not in table.columns
    assert table.columns[0] == "genome"


def test_two_samples_each_with_a_bin_001_are_two_rows_not_one(tmp_path):
    """Rows are keyed on sample and name, since every sample numbers its bins from bin_001."""
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 1, "sample_7": 1})))
    assert [(row[0], row[1]) for row in table.rows] == [
        ("sample_6", "bin_001"),
        ("sample_7", "bin_001"),
    ]


def test_each_bin_set_is_paired_with_its_own_gunc_report(tmp_path):
    """Each sample's genomes are joined to that sample's GUNC report.

    A single report joined onto every sample gave 23 wrong chimerism calls and an extra genome on
    the reference study. Here sample_6's bin_001 passes and sample_7's fails, so a crossed join
    gives the wrong value.
    """
    study = _study(
        tmp_path,
        {"sample_6": 1, "sample_7": 1},
        gunc={"sample_6": {"bin_001": 0.02}, "sample_7": {"bin_001": 0.91}},
    )
    table = M.build(M.discover(study))
    rows = {(row[0], row[1]): dict(zip(table.columns, row)) for row in table.rows}
    assert rows[("sample_6", "bin_001")]["gunc_pass"] == "yes"
    assert rows[("sample_7", "bin_001")]["gunc_pass"] == "no"


def test_a_sample_whose_gunc_did_not_run_is_blank_rather_than_the_neighbours_verdict(tmp_path):
    """A sample without a GUNC report gets blank GUNC cells, not another sample's values."""
    study = _study(tmp_path, {"sample_6": 1, "sample_7": 1}, gunc={"sample_6": {"bin_001": 0.02}})
    table = M.build(M.discover(study))
    at = table.columns.index("gunc_pass")
    by_sample = {row[0]: row[at] for row in table.rows}
    assert by_sample["sample_6"] == "yes"
    assert by_sample["sample_7"] == M.BLANK


def test_the_intermediate_tables_in_work_files_are_not_bin_sets(tmp_path):
    """Only the consolidated set is read, not the intermediate sets bin_refinement leaves.

    bin_refinement leaves seven per sample (binsA, binsAB, binsO.to_score, ...); reading them
    would count every genome several times.
    """
    study = _study(tmp_path, {"sample_6": 2})
    work = tmp_path / "study" / "BIN_REFINEMENT" / "sample_6" / "work_files"
    work.mkdir(exist_ok=True)
    for stem in ("binsA", "binsAB", "binsO.to_score"):
        with open(work / ("%s.stats.tsv" % stem), "w") as fh:
            fh.write(STATS_HEADER + "bin_001\t50.0\t9.0\t0.40\troot\t1000\t900000\n")
    table = M.build(M.discover(study))
    assert table.bin_sets == 1
    assert len(table.rows) == 2


def test_a_study_wide_table_matching_no_bin_is_recorded_rather_than_guessed(tmp_path):
    """A cohort's taxonomy is keyed on dereplicated names, which match no sample's bin_001.

    The cells stay blank and the bin set is named in `unjoined`, so it can be told apart from an
    input that was never written.
    """
    study = _study(tmp_path, {"sample_6": 2, "sample_7": 2})
    classify = tmp_path / "study" / "CLASSIFY"
    classify.mkdir()
    (classify / "bin_taxonomy.tsv").write_text("derep_cluster_1.fasta\tBacillota;Clostridia\n")
    table = M.build(M.discover(study))
    assert ("taxonomy", "sample_6") in table.unjoined
    assert ("taxonomy", "sample_7") in table.unjoined
    assert "taxonomy" not in table.absent  # the file is there; it is the join that is not
    at = table.columns.index("taxonomy")
    assert all(row[at] == M.BLANK for row in table.rows if row[0])


def test_a_genome_named_only_by_a_study_wide_table_has_a_blank_sample(tmp_path):
    """A genome from an unknown bin set still gets a row, with the sample left blank."""
    study = _study(tmp_path, {"sample_6": 1, "sample_7": 1})
    classify = tmp_path / "study" / "CLASSIFY"
    classify.mkdir()
    (classify / "bin_taxonomy.tsv").write_text("derep_1.fasta\tBacillota\n")
    table = M.build(M.discover(study))
    extra = [row for row in table.rows if row[1] == "derep_1"]
    assert len(extra) == 1
    assert extra[0][0] == M.BLANK
