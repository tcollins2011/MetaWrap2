"""Tests for community-level parsers: SingleM's two tables, Nonpareil's .npo, PPanGGOLiN and
IQ-TREE's Newick.

SingleM's ``read_fraction`` is a percentage (``read_fraction.py`` computes
``account / metagenome_size * 100``, capped at 100). Guessing the unit from the magnitude reads a
host-dominated sample's 0.8% as 80%.

PPanGGOLiN's statistics table starts with two ``#`` settings lines, which a plain DictReader takes
as the header. IQ-TREE writes both support values (``)97.5/100:``) under ``--alrt``, which a
single-number pattern does not match.

All fixtures are vendored tool output; see the SOURCE.md beside each.
"""

from __future__ import annotations

import csv
import os

import pytest

from metawrap2 import funcabund, nonpareil, pangenome, phylogeny, recovery

FIXTURES = os.path.dirname(__file__) + "/fixtures"
APPRAISE_SUMMARY = os.path.join(FIXTURES, "appraise", "appraise_summary.tsv")
PANGENOME = os.path.join(FIXTURES, "pangenome")
TREEFILE = os.path.join(FIXTURES, "phylogeny", "core.treefile")
NPO = os.path.join(FIXTURES, "nonpareil", "HumanGut.npo")

#: `singlem microbial_fraction`'s table, in read_fraction.py's layout: one host-dominated sample at
#: 0.80% and one typical at 85.30%, which a magnitude heuristic would treat differently.
MICROBIAL_FRACTION = (
    "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\n"
    "host_dominated\t8000000\t1000000000\t0.80\n"
    "typical\t853000000\t1000000000\t85.30\n"
)


# ── SingleM's read fraction, which is a percentage ────────────────────────────────────────────


def test_singlems_read_fraction_is_already_a_percentage(tmp_path):
    """0.80 in the file is 0.8% of the reads, i.e. 0.008.

    Values are always divided by 100.
    """
    path = tmp_path / "microbial_fraction.tsv"
    path.write_text(MICROBIAL_FRACTION)
    fractions = funcabund.parse_singlem_fraction(str(path))
    assert fractions["host_dominated"] == pytest.approx(0.008)
    assert fractions["typical"] == pytest.approx(0.853)
    # The same rule is applied to both rows.
    assert fractions["typical"] / fractions["host_dominated"] == pytest.approx(85.30 / 0.80)


def test_the_read_fraction_column_has_had_three_names_and_all_are_accepted(tmp_path):
    """``read_fraction``, ``prokaryotic_fraction``, ``bacterial_archaeal_fraction`` across releases.

    An unknown column name raises KeyError naming the accepted ones, instead of dropping the
    samples.
    """
    assert funcabund.SINGLEM_FRACTION_COLUMNS == (
        "read_fraction",
        "prokaryotic_fraction",
        "bacterial_archaeal_fraction",
    )
    for name in funcabund.SINGLEM_FRACTION_COLUMNS:
        path = tmp_path / ("%s.tsv" % name)
        path.write_text(MICROBIAL_FRACTION.replace("read_fraction", name))
        assert funcabund.parse_singlem_fraction(str(path))["typical"] == pytest.approx(0.853)

    renamed = tmp_path / "renamed.tsv"
    renamed.write_text(MICROBIAL_FRACTION.replace("read_fraction", "microbial_fraction"))
    with pytest.raises(KeyError) as caught:
        funcabund.parse_singlem_fraction(str(renamed))
    assert "read_fraction" in str(caught.value), "the message names the columns it looked for"


def test_a_fraction_above_one_hundred_is_refused_rather_than_rescaled(tmp_path):
    """A value above 100 raises ValueError, since SingleM caps its output at 100."""
    path = tmp_path / "microbial_fraction.tsv"
    path.write_text(MICROBIAL_FRACTION.replace("85.30", "8530.0"))
    with pytest.raises(ValueError) as caught:
        funcabund.parse_singlem_fraction(str(path))
    assert "above 100%" in str(caught.value)


# ── SingleM's appraise summary, which carries its own summary rows ────────────────────────────


def test_singlems_own_total_and_average_rows_are_not_samples():
    """``total`` and ``average`` rows in the ``sample`` column are skipped.

    The fixture has two real rows and four summary rows.
    """
    assert recovery.APPRAISAL_SUMMARY_ROWS == ("total", "average")
    with open(APPRAISE_SUMMARY, newline="") as fh:
        raw = [row["sample"] for row in csv.DictReader(fh, delimiter="\t")]
    assert raw.count("total") == 2 and raw.count("average") == 2

    rows = recovery.parse_appraisal(APPRAISE_SUMMARY)
    assert {row.sample for row in rows} == {"ERR011347_1"}
    assert len(rows) == 2, "one row per domain for the one real sample"
    bacteria = next(row for row in rows if row.domain == "d__Bacteria")
    # 23 of 199 coverage units.
    assert bacteria.percent_binned == pytest.approx(11.56, abs=0.01)


def test_an_absent_assembled_coverage_column_is_none_and_not_zero():
    """``coverage_assembled`` appears only when an assembly was appraised alongside the genomes.

    None means no assembly was given; 0.0 would mean none of the given assembly was binned. The
    column is absent in the common case.
    """
    header = set()
    with open(APPRAISE_SUMMARY, newline="") as fh:
        header.update(next(csv.reader(fh, delimiter="\t")))
    assert "coverage_assembled" not in header

    rows = recovery.parse_appraisal(APPRAISE_SUMMARY)
    assert all(row.assembled is None for row in rows)
    assert all(row.percent_assembled_not_binned is None for row in rows)
    # The absent value is written as blank, not zero.
    assert recovery.appraisal_row(rows[0])["coverage_assembled_not_binned"] == ""
    assert recovery.appraisal_row(rows[0])["percent_community_assembled_not_binned"] == ""


# ── Nonpareil's .npo ──────────────────────────────────────────────────────────────────────────


def test_npo_metadata_keys_are_nonpareils_own_single_letter_names():
    """Keys such as ``# @R``, ``# @L``, ``# @overlap``, ``# @divide`` are read as Nonpareil's.

    ``R`` is the read count and ``L`` the mean read length; the file does not say so. Swapping
    them would rescale the effort axis by about 50,000 on this sample.
    """
    curve = nonpareil.read_npo(NPO, label="HumanGut")
    assert curve.reads == pytest.approx(4158944.0)
    assert curve.read_length == pytest.approx(77.680)
    assert curve.overlap == pytest.approx(50.0)
    assert curve.log_sample == pytest.approx(0.70)
    assert curve.version == "2.40"
    # No k-mer size in the header, so this is the alignment kernel.
    assert curve.kernel == "alignment"
    assert 0.0 < curve.C <= 1.0, "coverage is a fraction, not a percentage"


def test_ppanggolins_statistics_table_opens_with_two_comment_lines():
    """``#soft_core=0.95`` and ``#duplication_margin=0.05`` come *before* the header row.

    A plain csv.DictReader takes the first settings line as the header, so ``Genome_name`` and
    ``Genes`` lookups return None.
    """
    path = os.path.join(PANGENOME, "genomes_statistics.tsv")
    with open(path, newline="") as fh:
        first_two = [next(fh), next(fh)]
    assert all(line.startswith("#") for line in first_two)
    with open(path, newline="") as fh:
        naive = csv.DictReader(fh, delimiter="\t")
        assert "Genome_name" not in (
            naive.fieldnames or []
        ), "the naive read sees the settings line"

    genes = pangenome.read_genome_statistics(path)
    assert pangenome.STATS_GENOME == "Genome_name" and pangenome.STATS_GENES == "Genes"
    assert len(genes) == 3
    assert genes["metawrap_50_10_bins__bin_017"] == pytest.approx(2394.0)


def test_partition_families_are_counted_per_group_and_a_missing_run_is_none():
    """``partitions/persistent.txt`` is one family name per line - a group-level count.

    The statistics table's ``Persistent_families`` is per genome, so summing it would count
    shared families once per genome. From the partition files: 1,509 persistent, 3,406 shell,
    2,446 cloud.
    """
    summary = pangenome.parse_partitions(PANGENOME, "group_1", reliable=True)
    assert summary.families == {"persistent": 1509, "shell": 3406, "cloud": 2446}
    assert summary.total_families == 7361
    assert summary.mean_genes == pytest.approx(3214.0)
    # Openness is a percentage of families, not a fraction.
    assert summary.openness == pytest.approx(100.0 - 1509 / 7361 * 100.0)
    assert summary.shares["persistent"] + summary.openness == pytest.approx(100.0)

    # A directory PPanGGOLiN never wrote to gives None, not zero families.
    assert pangenome.parse_partitions(str(PANGENOME) + "-never-ran", "group_2", True) is None


# ── IQ-TREE's support values ──────────────────────────────────────────────────────────────────


def test_iqtree_writes_two_support_values_and_the_bootstrap_is_the_second():
    """In ``)100/100:`` the values are SH-aLRT then UFBoot; the threshold applies to UFBoot.

    A single-number pattern found no branches in the combined form, and the first value is
    SH-aLRT, a different test.
    """
    summary = phylogeny.support_summary(TREEFILE)
    assert summary.branches == 1
    assert summary.supported == 1
    assert summary.fraction_supported == pytest.approx(1.0)
    assert summary.threshold == phylogeny.SUPPORTED_BRANCH == 95.0

    # Both IQ-TREE forms; the second value is the one compared.
    values = phylogeny._SUPPORT.findall("(a:0.1,(b:0.1,c:0.1)80.5/92:0.3,d:0.1);")
    assert values == [("80.5", "92")]
    assert phylogeny._SUPPORT.findall("(a:0.1,(b:0.1,c:0.1)100:0.3,d:0.1);") == [("100", "")]


def test_a_tree_with_no_branches_reports_none_support_rather_than_zero(tmp_path):
    """No tree, or one with no internal branches, gives None, not 0.0.

    0.0 would rank a variant last on resolution. A three-taxon unrooted tree has no internal branch.
    """
    absent = phylogeny.support_summary(str(tmp_path / "iqtree-never-ran.treefile"))
    assert absent.branches == 0
    assert absent.supported == 0
    assert absent.fraction_supported is None

    star = tmp_path / "star.treefile"
    star.write_text("(a:0.1,b:0.1,c:0.1);\n")
    assert phylogeny.support_summary(str(star)).fraction_supported is None
