"""Tests for parsing inStrain and TRACS tables: column names, units, and genome names.

Complements ``tests/test_strains.py``. A renamed column parses to None, which reads as "not in both
samples", so a schema mistake turns a comparison table into an apparent absence (all 42
comparisons of the first run).

The fixtures are vendored inStrain 1.10.0 and TRACS output over three human gut metagenomes and
this repo's 26 refined bins; see ``tests/fixtures/strains/SOURCE.md``.
"""

from __future__ import annotations

import csv
import os

import pytest

from metawrap2 import strains

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "strains")
PROFILE = os.path.join(FIXTURES, "genome_info.tsv")
COMPARE = os.path.join(FIXTURES, "genomeWide_compare.tsv")
TRACS = os.path.join(FIXTURES, "tracs_transmission_distances.csv")
SAMPLES = ("ERR011347", "ERR011348")


def _header(path, delimiter="\t"):
    with open(path, newline="") as fh:
        return next(csv.reader(fh, delimiter=delimiter))


# ── inStrain's column names, which are not guessable ─────────────────────────────────────────


def test_instrain_profile_columns_are_spelled_the_way_instrain_spells_them():
    """``breadth_minCov`` and ``breadth_expected``: camelCase in one, suffix order reversed in the other.

    Names are taken from the real `genome_info.tsv` header, since inStrain's documentation names
    some differently. A wrong name gives None and the genome looks undetected.
    """
    header = set(_header(PROFILE))
    for key, column in strains.PROFILE_COLUMNS.items():
        assert column in header, "%s -> %r is not a column inStrain writes" % (key, column)
    assert strains.PROFILE_COLUMNS["breadth_min_cov"] == "breadth_minCov"
    assert strains.PROFILE_COLUMNS["expected_breadth"] == "breadth_expected"
    assert strains.PROFILE_COLUMNS["snvs"] == "SNV_count"
    assert strains.PROFILE_COLUMNS["consensus_snvs"] == "SNS_count"


def test_a_missing_identity_column_is_an_error_rather_than_an_empty_table(tmp_path):
    """A missing `popANI` column raises KeyError instead of parsing to None.

    Without the check, 42 rows of None read as "different strain over too little of the genome".
    """
    header = set(_header(COMPARE))
    for key, column in strains.COMPARE_COLUMNS.items():
        assert column in header, "%s -> %r is not a column inStrain writes" % (key, column)
    # `population_SNPs` here, `SNV_count` in the profile table: inStrain's two outputs name the
    # same quantity differently.
    assert strains.COMPARE_COLUMNS["population_snvs"] == "population_SNPs"
    assert strains.COMPARE_COLUMNS["a"] == "name1"

    renamed = tmp_path / "genomeWide_compare.tsv"
    with open(COMPARE) as src:
        renamed.write_text(src.read().replace("popANI", "popANI_reference", 1))
    with pytest.raises(KeyError) as caught:
        strains.parse_comparison(str(renamed))
    assert "genomeWide_compare.tsv" in str(caught.value)

    profile = tmp_path / "genome_info.tsv"
    profile.write_text("coverage\tbreadth\n1.0\t0.5\n")
    with pytest.raises(KeyError):
        strains.parse_profile(str(profile), "S1")


# ── the unit, applied unconditionally ─────────────────────────────────────────────────────────


def test_the_three_fraction_columns_are_scaled_whatever_their_magnitude():
    """Fractions are scaled to percentages per column, not by magnitude.

    popANI is ~0.99999 and percent_compared ranges from 0.001 to 1.0. A rule such as "above 1 is
    already a percentage" would leave 0.001 as 0.001% and scale 1.0 to 100%.
    """
    assert strains.FRACTION_COLUMNS == ("pop_ani", "con_ani", "percent_compared")
    assert strains._as_percent("1.0") == pytest.approx(100.0)
    assert strains._as_percent("0.001") == pytest.approx(0.1)
    assert strains._as_percent("0.9999994540295621") == pytest.approx(99.99994540295621)

    rows = {c.genome: c for c in strains.parse_comparison(COMPARE)}
    with open(COMPARE, newline="") as fh:
        raw = {row["genome"]: row for row in csv.DictReader(fh, delimiter="\t")}
    for genome, comparison in rows.items():
        assert comparison.pop_ani == pytest.approx(float(raw[genome]["popANI"]) * 100.0)


def test_a_blank_cell_is_none_and_not_a_zero(tmp_path):
    """inStrain leaves cells empty, and writes ``nan``, for quantities it did not compute.

    These parse to None ("not measured"), not zero; a zero breadth would put a genome below the
    floor instead of outside the comparison.
    """
    assert strains._as_float("") is None
    assert strains._as_float("nan") is None
    assert strains._as_float("0") == 0.0

    path = tmp_path / "genome_info.tsv"
    path.write_text(
        "genome\tcoverage\tbreadth\tnucl_diversity\tSNV_count\n" "bin.13\t\t0.0\tnan\t\n"
    )
    profile = strains.parse_profile(str(path), "S1")[0]
    assert profile.coverage is None
    assert profile.nucleotide_diversity is None
    assert profile.snvs is None
    assert profile.breadth == 0.0


# ── identifiers: the dot in a metaWRAP bin name ───────────────────────────────────────────────


def test_a_metawrap_bin_name_keeps_its_numeric_suffix():
    """`bin.13` is a bin name, not `bin` with a `.13` extension.

    ``os.path.splitext("bin.13")`` returns ``("bin", ".13")``, which would collapse a 26-bin set
    onto the single key `bin`.
    """
    assert os.path.splitext("bin.13") == ("bin", ".13")
    assert strains._strip_extension("bin.13") == "bin.13"
    assert strains._strip_extension("bin.13.fa") == "bin.13"
    assert strains._strip_extension("bin.13.fasta") == "bin.13"
    # The real tables keep the dotted names.
    assert "bin.11" in {p.genome for p in strains.parse_profile(PROFILE, "S1")}
    assert "bin.1" in {c.genome for c in strains.parse_comparison(COMPARE)}


def test_tracs_column_names_carry_their_spaces_and_the_file_is_comma_separated():
    """``SNP distance``, ``sites considered``, ``expected K`` - spaces, and a `.csv` not a `.tsv`.

    Read as a TSV, the header is one field and every row parses to None, reported as ABSENT for
    every pair.
    """
    header = _header(TRACS, delimiter=",")
    assert header[0] == "sampleA"
    for key, column in strains.TRACS_COLUMNS.items():
        assert column in header, "%s -> %r is not a column TRACS writes" % (key, column)
    assert strains.TRACS_COLUMNS["snp_distance"] == "SNP distance"
    assert strains.TRACS_COLUMNS["sites"] == "sites considered"


def test_the_genome_is_split_off_by_a_known_sample_name_not_at_an_underscore():
    """TRACS labels a row ``<sample>_<genome>``, and genomes are called things like `bin_001`.

    Splitting at an underscore cannot work (the last gives sample `ERR011347_bin` and genome
    `001`). Known sample names are tried longest first so one that prefixes another cannot claim
    its rows.
    """
    assert strains._split_label("ERR011347_bin.1", SAMPLES) == ("ERR011347", "bin.1")
    assert strains._split_label("ERR011347_bin_001", SAMPLES) == ("ERR011347", "bin_001")
    assert strains._split_label("S1_g", ("S1", "S1b")) == ("S1", "g")
    assert strains._split_label("S1b_g", ("S1", "S1b")) == ("S1b", "g")
    # With no known sample matching, the genome is reported as unknown.
    assert strains._split_label("unexpected_label", SAMPLES) == ("unexpected_label", "")


def test_an_empty_tracs_table_is_no_rows_rather_than_a_zero_distance(tmp_path):
    """A missing or header-only TRACS table gives no rows; a blank row is ABSENT, not distance 0.

    TRACS skips any reference under 25% covered, so missing rows are common. A distance of 0 would
    be a transmission call.
    """
    assert strains.parse_tracs(str(tmp_path / "never-ran.csv")) == []

    header_only = tmp_path / "transmission_distances.csv"
    with open(TRACS) as src:
        header_only.write_text(src.readline())
    assert strains.parse_tracs(str(header_only), samples=SAMPLES) == []

    blank = tmp_path / "blank_row.csv"
    with open(TRACS) as src:
        lines = src.read().splitlines()
    blank.write_text(lines[0] + "\n" + "ERR011347_bin.1,ERR011348_bin.1,NA,,NA,NA,,,combinedbin\n")
    row = strains.parse_tracs(str(blank), samples=SAMPLES)[0]
    assert row.snp_distance is None
    assert row.sites_considered is None
    assert row.verdict is strains.StrainVerdict.ABSENT
