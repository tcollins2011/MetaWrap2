"""MIMAG quality tags: the tier boundaries, the unverified-rRNA state, and the stats columns."""

import os

from metawrap2 import checkm, mimag
from metawrap2.mimag import MimagTag, RnaEvidence

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "mimag")

#: A genome meeting every criterion, so a test can change one of them.
COMPLETE = RnaEvidence(rrna_5s=True, rrna_16s=True, rrna_23s=True, trna_types=18)


def _lines(path):
    """The lines of a written stats table, so no test leaves a file handle open."""
    with open(path) as fh:
        return fh.read().splitlines()


def _stats(tmp_path, rows, name="metawrap_70_10_bins.stats.tsv"):
    """A .stats.tsv with the columns metawrap2.checkm.write_stats produces."""
    path = tmp_path / name
    path.write_text(
        "\t".join(checkm.STATS_HEADER)
        + "\n"
        + "".join("%s\t%s\t%s\t0.500\tBacteria\t20000\t3000000\n" % row for row in rows)
    )
    return str(path)


# --- the tier boundaries ---------------------------------------------------------------------


def test_eighteen_trna_types_pass_and_seventeen_do_not():
    """The 18/20 tRNA boundary."""
    at_18 = mimag.assess("bin_001", 95.0, 1.0, RnaEvidence(True, True, True, 18))
    at_17 = mimag.assess("bin_002", 95.0, 1.0, RnaEvidence(True, True, True, 17))
    assert at_18.tag is MimagTag.HIGH
    assert at_17.tag is MimagTag.MEDIUM


def test_high_quality_needs_completeness_strictly_above_ninety():
    """Bowers et al. require ">90%", so 90.00 is not high quality even with every rRNA."""
    assert mimag.assess("b", 90.01, 1.0, COMPLETE).tag is MimagTag.HIGH
    assert mimag.assess("b", 90.0, 1.0, COMPLETE).tag is MimagTag.MEDIUM


def test_high_quality_needs_contamination_strictly_below_five():
    """Bowers et al. require "<5%", so 5.00 is not high quality. The benchmark's 90/5 tier uses
    <=5, and the two must stay separate."""
    assert mimag.assess("b", 99.0, 4.99, COMPLETE).tag is MimagTag.HIGH
    assert mimag.assess("b", 99.0, 5.0, COMPLETE).tag is MimagTag.MEDIUM


def test_medium_quality_starts_at_fifty_percent_complete_inclusive():
    """MIMAG medium is ">=50%" while high is ">90%", so a 50.0% bin is medium quality."""
    assert mimag.assess("b", 50.0, 9.0).tag is MimagTag.MEDIUM
    assert mimag.assess("b", 49.99, 9.0).tag is MimagTag.LOW


def test_contamination_of_ten_percent_or_more_is_not_any_mimag_tier():
    """Every MIMAG tier, including low quality, requires <10% contamination, so a 30%
    contaminated bin gets no tier."""
    assert mimag.assess("b", 95.0, 10.0, COMPLETE).tag is MimagTag.NOT_MIMAG
    assert mimag.assess("b", 20.0, 30.0).tag is MimagTag.NOT_MIMAG
    # 9.99% is inside the ceiling: medium, since it is over high quality's 5%.
    assert mimag.assess("b", 95.0, 9.99, COMPLETE).tag is MimagTag.MEDIUM
    assert mimag.assess("b", 20.0, 9.99).tag is MimagTag.LOW


def test_missing_completeness_is_unscored_rather_than_zero():
    """A bin CheckM could not score is tagged unscored, not low quality."""
    assert mimag.assess("b", None, 1.0, COMPLETE).tag is MimagTag.UNSCORED
    assert mimag.assess("b", 95.0, None, COMPLETE).tag is MimagTag.UNSCORED


# --- unmeasured versus measured-absent rRNA --------------------------------------------------


def test_unmeasured_rrna_is_distinct_from_rrna_measured_absent():
    """Two 95/1 genomes: one searched with no 16S found, one never searched. Neither is high
    quality, and only the first is downgraded to medium; the second is unverified."""
    searched = mimag.assess("found_none", 95.0, 1.0, RnaEvidence(True, False, True, 20))
    never_looked = mimag.assess("no_evidence", 95.0, 1.0, None)

    assert searched.tag is MimagTag.MEDIUM
    assert never_looked.tag is MimagTag.HIGH_UNVERIFIED
    assert searched.tag is not never_looked.tag
    assert searched.evidence.describe_rrna() == "5S+23S"
    assert never_looked.evidence.describe_rrna() == mimag.UNKNOWN


def test_rrna_searched_and_entirely_absent_reads_as_none_not_unknown():
    """A genome searched by Infernal with no rRNA found reads "none", not "unknown"."""
    searched = RnaEvidence(rrna_5s=False, rrna_16s=False, rrna_23s=False, trna_types=0)
    assert searched.describe_rrna() == "none"
    assert RnaEvidence().describe_rrna() == mimag.UNKNOWN


def test_measured_rrna_with_unmeasured_trnas_stays_unverified():
    """With Infernal only, rRNAs are present and tRNAs uncounted, so the genome stays
    unverified."""
    infernal_only = RnaEvidence(rrna_5s=True, rrna_16s=True, rrna_23s=True, trna_types=None)
    assert mimag.assess("b", 97.0, 2.0, infernal_only).tag is MimagTag.HIGH_UNVERIFIED


def test_a_measured_absence_decides_even_when_the_rest_is_unmeasured():
    """A measured missing 16S rules out high quality even when tRNAs were not counted."""
    partial = RnaEvidence(rrna_5s=None, rrna_16s=False, rrna_23s=None, trna_types=None)
    assert mimag.assess("b", 97.0, 2.0, partial).tag is MimagTag.MEDIUM


def test_a_bin_below_the_high_tier_is_not_marked_unverified():
    """Only 90/5 bins can be unverified; a 60% complete bin is medium quality."""
    assert mimag.assess("b", 60.0, 4.0, None).tag is MimagTag.MEDIUM
    assert mimag.assess("b", 20.0, 4.0, None).tag is MimagTag.LOW


# --- parsing the evidence the pipeline already writes ----------------------------------------


def test_bakta_tsv_gives_both_rrna_and_trna_evidence():
    """Bakta's TSV gives both rRNA and tRNA evidence in one file, so it is preferred over the GFF
    and the standalone tools."""
    found = mimag.parse_bakta_tsv(os.path.join(FIXTURES, "bakta", "bin_001", "bin_001.tsv"))
    assert found.rrna == (True, True, True)
    assert found.trna_types == 18
    assert mimag.assess("bin_001", 98.0, 1.0, found).tag is MimagTag.HIGH


def test_trnas_are_counted_as_amino_acids_not_as_features():
    """MIMAG counts amino acids (18 of 20). bin_001 has 19 tRNA rows, one of them tRNA-Xxx, which
    names no amino acid and is not counted."""
    found = mimag.parse_bakta_tsv(os.path.join(FIXTURES, "bakta", "bin_001", "bin_001.tsv"))
    assert found.trna_types == 18


def test_bakta_gff_agrees_with_the_tsv_on_the_same_genome():
    """The GFF fallback and the TSV parser give the same verdict for one annotation."""
    tsv = mimag.parse_bakta_tsv(os.path.join(FIXTURES, "bakta", "bin_001", "bin_001.tsv"))
    gff = mimag.parse_bakta_gff(os.path.join(FIXTURES, "bakta", "bin_001", "bin_001.gff3"))
    assert (gff.rrna, gff.trna_types) == (tsv.rrna, tsv.trna_types)


def test_bakta_tsv_records_a_missing_rrna_as_measured_absent():
    """bin_002 has no 16S row, which is returned as False (measured absent), not None."""
    found = mimag.parse_bakta_tsv(os.path.join(FIXTURES, "bakta", "bin_002", "bin_002.tsv"))
    assert found.rrna == (True, False, True)
    assert found.trna_types == 17


def test_infernal_does_not_mistake_the_eukaryotic_5_8s_family_for_5s(tmp_path):
    """RF00002 (5_8S_rRNA, eukaryotic) is not counted as a 5S hit by a substring match."""
    source = os.path.join(FIXTURES, "rna", "bin_003.infernal.tblout")
    with open(source) as fh:
        kept = [line for line in fh if line.startswith("#") or "5_8S_rRNA" in line]
    only_5_8s = tmp_path / "eukaryote.infernal.tblout"
    only_5_8s.write_text("".join(kept))

    assert len(kept) == 3, "the fixture should hold exactly one 5_8S_rRNA row"
    assert mimag.parse_infernal_tblout(str(only_5_8s)).rrna_5s is False


def test_infernal_alone_leaves_the_trna_criterion_unmeasured():
    """Infernal evidence reports no tRNA count (not 0), since Infernal does not search tRNAs."""
    found = mimag.parse_infernal_tblout(os.path.join(FIXTURES, "rna", "bin_003.infernal.tblout"))
    assert found.rrna == (True, True, True)
    assert found.trna_types is None


def test_trnascan_ignores_undetermined_and_pseudogene_rows():
    """bin_003 has 18 amino-acid rows plus an Undet and a Pseudo row, which are not counted."""
    assert mimag.parse_trnascan(os.path.join(FIXTURES, "rna", "bin_003.trnascan.tsv")) == 18


def test_the_fallback_pair_is_merged_into_one_evidence_record():
    """bin_003's Infernal (rRNA) and tRNAscan-SE (tRNA) outputs merge into one record that meets
    the criterion."""
    found = mimag.evidence_for_genome(FIXTURES, "bin_003")
    assert found is not None
    assert found.rrna == (True, True, True)
    assert found.trna_types == 18
    assert mimag.assess("bin_003", 96.0, 2.0, found).tag is MimagTag.HIGH


def test_a_genome_with_no_annotation_output_yields_no_evidence_record():
    """Returns None, not an empty RnaEvidence, so "nothing ran" differs from "found nothing";
    collect_evidence omits the key."""
    assert mimag.evidence_for_genome(FIXTURES, "bin_never_annotated") is None
    assert mimag.collect_evidence(FIXTURES, ["bin_001", "bin_never_annotated"]).keys() == {
        "bin_001"
    }


# --- the stats table ------------------------------------------------------------------------


def test_stats_columns_are_appended_after_the_frozen_base_columns(tmp_path):
    """MIMAG columns are appended after the base columns, which are read by position
    (choose_best_bin uses 0,1,2,5) and pinned by a golden test."""
    path = _stats(tmp_path, [("bin_001", "95.00", "1.000"), ("bin_002", "60.00", "4.000")])
    before = [line.split("\t") for line in _lines(path)]

    mimag.annotate_stats(path)

    after = [line.split("\t") for line in _lines(path)]
    assert after[0][: len(checkm.STATS_HEADER)] == list(checkm.STATS_HEADER)
    assert after[0][len(checkm.STATS_HEADER) :] == list(mimag.STATS_COLUMNS)
    for old, new in zip(before, after):
        assert new[: len(old)] == old


def test_stats_rows_are_tagged_from_their_own_completeness_and_contamination(tmp_path):
    """Each row is tagged from its own completeness and contamination."""
    path = _stats(
        tmp_path,
        [
            ("bin_001", "95.00", "1.000"),
            ("bin_002", "60.00", "4.000"),
            ("bin_003", "30.00", "2.000"),
            ("bin_004", "95.00", "30.00"),
        ],
    )
    tags = {a.bin: a.tag for a in mimag.annotate_stats(path)}
    assert tags == {
        "bin_001": MimagTag.HIGH_UNVERIFIED,
        "bin_002": MimagTag.MEDIUM,
        "bin_003": MimagTag.LOW,
        "bin_004": MimagTag.NOT_MIMAG,
    }


def test_an_absent_rrna_source_leaves_the_column_unknown_rather_than_failing_the_table(tmp_path):
    """Without annotate output the rRNA column is "unknown" and the table is still written."""
    path = _stats(tmp_path, [("bin_001", "95.00", "1.000")])
    mimag.annotate_stats(path, evidence=None)
    header, row = (line.split("\t") for line in _lines(path))
    cells = dict(zip(header, row))
    assert cells["mimag_tag"] == "high_quality_unverified_rrna"
    assert cells["mimag_rrna"] == "unknown"
    assert cells["mimag_trnas"] == ""


def test_annotating_a_stats_table_twice_does_not_duplicate_the_columns(tmp_path):
    """Re-running the tagging step (e.g. on --resume) does not add the columns again, as in
    metawrap2.gunc.annotate_stats."""
    path = _stats(tmp_path, [("bin_001", "95.00", "1.000")])
    mimag.annotate_stats(path)
    mimag.annotate_stats(path)
    header = _lines(path)[0].split("\t")
    assert header.count("mimag_tag") == 1
    assert len(header) == len(checkm.STATS_HEADER) + len(mimag.STATS_COLUMNS)


def test_tag_mimag_resolves_evidence_from_an_annotate_output_directory(tmp_path):
    """checkm.tag_mimag finds Bakta output by genome name in an annotate directory."""
    path = _stats(tmp_path, [("bin_001", "98.00", "1.000"), ("bin_002", "95.00", "1.000")])
    tags = {a.bin: a.tag for a in checkm.tag_mimag(path, evidence_dir=FIXTURES)}
    assert tags == {"bin_001": MimagTag.HIGH, "bin_002": MimagTag.MEDIUM}


def test_the_summary_line_shows_high_and_unverified_high_side_by_side(tmp_path):
    """The summary printed by bin_refinement shows both the high and unverified-high counts."""
    path = _stats(tmp_path, [("bin_001", "95.00", "1.000"), ("bin_002", "60.00", "4.000")])
    line = mimag.summarise(mimag.annotate_stats(path))
    assert "0 high quality," in line
    assert "1 high quality with rRNA/tRNA unmeasured" in line
    assert "1 medium quality" in line
    assert "low quality" not in line
