"""The per-genome table, which is the row that did not exist anywhere before.

Quality was in one file, taxonomy in another and coverage in a third, and the join between them is
where this can go wrong in ways a reader cannot see. These tests pin the two failures that would
make the table worse than no table at all: a genome lost because two files spell its name
differently, and an unmeasured cell written as a zero.
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

    One ``BIN_REFINEMENT/<sample>/`` per sample, each with its own consolidated bin set and its own
    ``work_files/gunc/`` report - which is the real layout of the study this was written against,
    and the one a single-bin-set reader got wrong. Every sample numbers its bins from ``bin_001``.
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

    One genome, three files, three conventions. Without a shared stripping rule this is three
    genomes with two thirds of the columns blank in each - and nothing on the row says so.
    """
    row = _by_genome(_table(tmp_path))["bin_001"]
    assert row["completeness"] == "98.65"  # from the stats table
    assert row["taxonomy"].startswith("Bacillota;")  # from bin_taxonomy.tsv
    assert row["gunc_pass"] == "yes"  # from GUNC's own report
    assert row["abundance.sampleA"] == "143.178035"  # from the abundance matrix


def test_a_bin_named_with_a_dot_is_not_collapsed_to_bin(tmp_path):
    """``bin.2`` is a whole genome's name; ``os.path.splitext`` would make it ``bin``.

    Every metaWRAP-era bin set is named this way, so a join that used splitext would merge an
    entire bin set into one row called ``bin`` and keep whichever file was read last.
    """
    rows = _by_genome(_table(tmp_path))
    assert "bin.2" in rows
    assert "bin" not in rows
    assert rows["bin.2"]["taxonomy"].startswith("Bacteroidota;")


def test_a_genome_only_in_the_abundance_matrix_still_gets_a_row(tmp_path):
    """No input may drop a genome, or the table silently answers a narrower question.

    ``bin_004`` was quantified but never scored or classified - the state of every bin in a run
    where quantification covered the pre-refinement set.
    """
    row = _by_genome(_table(tmp_path))["bin_004"]
    assert row["abundance.sampleB"] == "4.200000"
    assert row["completeness"] == M.BLANK


def test_the_consolidated_bin_set_is_preferred_over_a_binners_raw_table(tmp_path):
    """A bin_refinement directory holds four stats tables and alphabetical order picks CONCOCT's.

    The decoy in the fixture names genomes nothing else mentions, so choosing it would produce a
    report about a different set of genomes with no page saying which.
    """
    inputs = M.discover(_run(tmp_path, "REFINEMENT"))
    (bin_set,) = inputs.bin_sets
    assert os.path.basename(bin_set.stats or "") == "metawrap_70_10_bins.stats.tsv"
    assert "concoct_42" not in _by_genome(M.build(inputs))


def test_the_unbinned_catch_all_is_not_a_genome(tmp_path):
    """`binning` writes concoct_bins/unbinned.fasta on purpose; it is leftovers, not a MAG.

    A row for it would be counted in every per-genome total on the report's pages.
    """
    stats = os.path.join(FIXTURES, "REFINEMENT", "concoct_bins.stats.tsv")
    assert "unbinned" not in M.read_stats(stats)
    assert "concoct_42" in M.read_stats(stats)


# ── blank is not zero ────────────────────────────────────────────────────────────────────────


def test_a_measured_zero_survives_as_a_zero(tmp_path):
    """A contamination of 0 is a finding and has to read as one.

    The guard below only matters because this one holds: if unmeasured cells became zeros, a
    genuine 0.0 would be indistinguishable from a bin CheckM never scored.
    """
    row = _by_genome(_table(tmp_path))["bin.2"]
    assert row["contamination"] == "0.0"
    assert row["abundance.sampleA"] == "0.0"


def test_an_unscored_genome_has_blank_scores_and_no_tier(tmp_path):
    """CheckM wrote the row but could not score it, so nothing is written where a score goes.

    This project has been burned by exactly this in three other tools' tables: a 0 completeness
    reads as a total failure of the genome rather than of the measurement, and it gets believed.
    """
    row = _by_genome(_table(tmp_path))["bin_003"]
    assert row["completeness"] == M.BLANK
    assert row["contamination"] == M.BLANK
    assert row["quality_tier"] == M.BLANK
    assert row["mimag_tag"] == M.BLANK


def test_an_uncounted_contig_count_is_blank_rather_than_zero(tmp_path):
    """A genome absent from the membership table has an unknown contig count, not zero contigs.

    Zero contigs is an empty FASTA, which is a different and much more alarming result.
    """
    rows = _by_genome(_table(tmp_path))
    assert rows["bin_001"]["contigs"] == "2"
    assert rows["bin_003"]["contigs"] == M.BLANK


def test_a_bin_gunc_did_not_score_is_blank_rather_than_failed(tmp_path):
    """Writing it as passing defeats the check; writing it as failed condemns it unmeasured."""
    rows = _by_genome(_table(tmp_path))
    assert rows["bin.2"]["gunc_pass"] == "no"  # GUNC did score this one, and it failed
    assert rows["bin_003"]["gunc_pass"] == M.BLANK
    assert rows["bin_003"]["gunc_css"] == M.BLANK


def test_a_classified_but_unplaced_genome_is_blank_not_missing(tmp_path):
    """classify_bins writes a row with no lineage for a bin it could not place.

    That is a different statement from the bin never having been classified, and the only way to
    tell them apart on a report page is which inputs the table reports as absent.
    """
    table = _table(tmp_path)
    assert _by_genome(table)["bin_003"]["taxonomy"] == M.BLANK
    assert "taxonomy" not in table.absent  # the file was there; the cell was empty


# ── every input is optional ──────────────────────────────────────────────────────────────────


def test_a_run_that_skipped_classification_and_quantification_still_gets_a_table(tmp_path):
    """The table is the point; the optional inputs are columns, not preconditions.

    A missing module leaves its column blank for every genome and is named in ``absent``, so a
    reader is told the difference between "not classified" and "classified as nothing".
    """
    table = _table(tmp_path, "REFINEMENT")
    assert [row[0] for row in table.rows] == ["bin_001", "bin.2", "bin_003"]
    assert table.abundance_samples == []
    assert sorted(table.absent) == ["abundance", "taxonomy"]
    assert all(row[table.columns.index("taxonomy")] == M.BLANK for row in table.rows)


def test_a_run_with_only_an_abundance_matrix_is_still_a_table(tmp_path):
    """Nothing scored, nothing classified - and still one row per genome quantified."""
    table = _table(tmp_path, "QUANT")
    assert [row[0] for row in table.rows] == ["bin.2", "bin_001", "bin_004"]
    assert sorted(table.absent) == ["contigs", "gunc", "stats", "taxonomy"]


def test_no_inputs_at_all_is_an_empty_table_rather_than_an_error(tmp_path):
    """The report asks for a table before it knows whether there is one to build."""
    table = M.build(M.discover([str(tmp_path)]))
    assert table.rows == []
    assert sorted(table.absent) == sorted(M.INPUT_KINDS)


def test_an_older_runs_suffixes_are_still_found(tmp_path):
    """``.stats`` and ``.tab`` predate the normalised suffixes and the output still exists.

    A reader that only knew ``.stats.tsv`` would call this repo's own reference run absent in its
    entirety, which is the one answer a report must not give about a run that finished.
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
    """The benchmark counts tiers non-exclusively; a table needs one label per row.

    Counted the benchmark's way, the same genome would appear in three bars of the report's chart
    and the bars would sum to more than the number of genomes.
    """
    assert M.quality_tier(98.65, 1.21) == "near_complete"
    assert M.quality_tier(75.4, 0.0) == "high_quality"
    assert M.quality_tier(20.0, 1.0) == M.BELOW_TIERS
    assert M.quality_tier(99.0, None) == M.BLANK


def test_the_tier_chart_counts_every_genome_including_the_unmeasured_ones(tmp_path):
    """A chart that dropped unscored genomes would make a failed CheckM look like fewer MAGs."""
    table = _table(tmp_path)
    counts = M.tier_counts(table)
    assert sum(count for _label, count in counts) == len(table.rows)
    assert dict(counts)[M.NOT_MEASURED] == 2  # bin_003 was not scored, bin_004 not even seen


def test_the_written_tsv_is_the_table_with_its_header(tmp_path):
    """One header, one row per genome, tab separated, and no cell invented on the way out."""
    table = _table(tmp_path)
    path = M.write(table, str(tmp_path / "out" / M.TABLE_NAME))
    lines = Path(path).read_text().splitlines()
    assert lines[0].split("\t") == table.columns
    assert len(lines) == len(table.rows) + 1
    # The blank cells have to survive being written: a row of a TSV is where they are easiest to
    # lose, because a writer that skipped empty trailing cells would shorten the row.
    assert lines[3].split("\t") == table.rows[2]
    assert lines[3].endswith("\t\t")


# ── a study is more than one bin set ─────────────────────────────────────────────────────────


def test_every_samples_bin_set_is_in_the_table_not_just_the_first(tmp_path):
    """A per-sample study recovered 62 genomes and the table said 24 - one sample's bin set.

    `find_stats` returned one path, so two of three samples were missing and nothing on the page
    said so. A reader would conclude the study recovered a third of what it did.
    """
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 24, "sample_7": 15, "sample_13": 23})))
    assert table.bin_sets == 3
    assert len(table.rows) == 62
    counts = Counter(row[0] for row in table.rows)
    assert counts == {"sample_6": 24, "sample_13": 23, "sample_7": 15}


def test_more_than_one_bin_set_carries_a_sample_column_in_front_of_the_genome(tmp_path):
    """The sample is the directory the stats table sat in, and it has to come before the genome."""
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 2, "sample_7": 2})))
    assert table.columns[:2] == [M.SAMPLE_COLUMN, "genome"]
    assert sorted({row[0] for row in table.rows}) == ["sample_6", "sample_7"]


def test_a_single_bin_set_gets_no_sample_column(tmp_path):
    """A per-run report should not grow a field that says the same thing on every row."""
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 3})))
    assert table.bin_sets == 1
    assert M.SAMPLE_COLUMN not in table.columns
    assert table.columns[0] == "genome"


def test_two_samples_each_with_a_bin_001_are_two_rows_not_one(tmp_path):
    """Every sample numbers its bins from bin_001, so the name alone cannot key a row.

    Keyed on the name, each sample's bin_001 overwrites the last and a study collapses into the
    distinct names of its largest bin set - with no duplicate left to notice.
    """
    table = M.build(M.discover(_study(tmp_path, {"sample_6": 1, "sample_7": 1})))
    assert [(row[0], row[1]) for row in table.rows] == [
        ("sample_6", "bin_001"),
        ("sample_7", "bin_001"),
    ]


def test_each_bin_set_is_paired_with_its_own_gunc_report(tmp_path):
    """One GUNC report was found for the whole study and joined onto every sample's genomes.

    They join because both sets are numbered from bin_001, so the verdicts are wrong rather than
    missing: on the reference study that was 23 wrong chimerism calls plus a phantom genome from
    the larger bin set. Here sample_6's bin_001 passes and sample_7's fails, so a crossed join
    shows up as the wrong verdict rather than as no verdict.
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
    """Half a study scored is not a reason to lend one sample's verdicts to another."""
    study = _study(tmp_path, {"sample_6": 1, "sample_7": 1}, gunc={"sample_6": {"bin_001": 0.02}})
    table = M.build(M.discover(study))
    at = table.columns.index("gunc_pass")
    by_sample = {row[0]: row[at] for row in table.rows}
    assert by_sample["sample_6"] == "yes"
    assert by_sample["sample_7"] == M.BLANK


def test_the_intermediate_tables_in_work_files_are_not_bin_sets(tmp_path):
    """bin_refinement leaves seven per sample - binsA, binsAB, binsO.to_score - beside one answer.

    They are the arithmetic of consolidation, and counting them would report every genome several
    times over and inflate the study's yield.
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

    There is no join to make, so the cells stay blank and the bin set is named in `unjoined` - the
    page then says which, because a present-but-inapplicable input leaves the same blank cells as
    an input that was never written.
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
    """It is still a row - nothing drops a genome - but which bin set it came from is unknown."""
    study = _study(tmp_path, {"sample_6": 1, "sample_7": 1})
    classify = tmp_path / "study" / "CLASSIFY"
    classify.mkdir()
    (classify / "bin_taxonomy.tsv").write_text("derep_1.fasta\tBacillota\n")
    table = M.build(M.discover(study))
    extra = [row for row in table.rows if row[1] == "derep_1"]
    assert len(extra) == 1
    assert extra[0][0] == M.BLANK
