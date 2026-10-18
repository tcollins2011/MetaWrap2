"""inStrain's and TRACS' tables: the column names, the units, and where the genome name comes from.

Complements ``tests/test_strains.py``, which pins the two numbers the first real run got wrong.
These pin the *schema* those numbers depend on, because a renamed column here is not an error: it
parses to None, None reads as "not in both samples", and a whole comparison table comes back as a
biological absence. All 42 comparisons of the first real run did exactly that.

The fixtures are the vendored real ones - inStrain 1.10.0 and TRACS output over three human gut
metagenomes and this repo's 26 refined bins; see ``tests/fixtures/strains/SOURCE.md``.
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

    Every one of these is read off the real `genome_info.tsv` header rather than from inStrain's
    documentation, because the documentation names some of them differently. A miss is silent - the
    field becomes None and the genome looks undetected.
    """
    header = set(_header(PROFILE))
    for key, column in strains.PROFILE_COLUMNS.items():
        assert column in header, "%s -> %r is not a column inStrain writes" % (key, column)
    assert strains.PROFILE_COLUMNS["breadth_min_cov"] == "breadth_minCov"
    assert strains.PROFILE_COLUMNS["expected_breadth"] == "breadth_expected"
    assert strains.PROFILE_COLUMNS["snvs"] == "SNV_count"
    assert strains.PROFILE_COLUMNS["consensus_snvs"] == "SNS_count"


def test_a_missing_identity_column_is_an_error_rather_than_an_empty_table(tmp_path):
    """The one column that is checked explicitly, because losing it loses the whole result.

    Without the guard, `popANI` going missing yields 42 rows of None, which the verdict logic reads
    as "different strain over too little of the genome" - a complete, plausible, wrong answer.
    """
    header = set(_header(COMPARE))
    for key, column in strains.COMPARE_COLUMNS.items():
        assert column in header, "%s -> %r is not a column inStrain writes" % (key, column)
    # `population_SNPs` here against `SNV_count` in the profile table: inStrain is inconsistent with
    # itself across its two outputs, so neither table's names can be inferred from the other's.
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
    """popANI ~0.99999 and percent_compared from 0.001 to 1.0 - so a magnitude rule cannot work.

    "Above 1 must already be a percentage" is right for popANI and catastrophic for
    percent_compared, where it would leave 0.001 as 0.001% for one genome and scale another's 1.0
    to 100%. The scaling is a property of the column, so it is listed by column.
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

    A zero here is a measurement: zero coverage, zero diversity, zero SNVs. None is "not measured".
    The verdict logic distinguishes them, and a zero breadth would put a genome below the floor
    rather than outside the comparison.
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

    ``os.path.splitext("bin.13")`` returns ``("bin", ".13")``, which collapses every bin in a set
    onto the single key `bin` - so a 26-bin set joins as one genome and 25 profiles are lost into
    whichever row happened to be last.
    """
    assert os.path.splitext("bin.13") == ("bin", ".13")
    assert strains._strip_extension("bin.13") == "bin.13"
    assert strains._strip_extension("bin.13.fa") == "bin.13"
    assert strains._strip_extension("bin.13.fasta") == "bin.13"
    # And the real tables come back with the dotted names intact.
    assert "bin.11" in {p.genome for p in strains.parse_profile(PROFILE, "S1")}
    assert "bin.1" in {c.genome for c in strains.parse_comparison(COMPARE)}


def test_tracs_column_names_carry_their_spaces_and_the_file_is_comma_separated():
    """``SNP distance``, ``sites considered``, ``expected K`` - spaces, and a `.csv` not a `.tsv`.

    Read as a TSV the whole header is one field and every row parses to None, which the verdict
    logic reports as ABSENT: "these samples do not share this genome", for every pair.
    """
    header = _header(TRACS, delimiter=",")
    assert header[0] == "sampleA"
    for key, column in strains.TRACS_COLUMNS.items():
        assert column in header, "%s -> %r is not a column TRACS writes" % (key, column)
    assert strains.TRACS_COLUMNS["snp_distance"] == "SNP distance"
    assert strains.TRACS_COLUMNS["sites"] == "sites considered"


def test_the_genome_is_split_off_by_a_known_sample_name_not_at_an_underscore():
    """TRACS labels a row ``<sample>_<genome>``, and genomes are called things like `bin_001`.

    Splitting at the last underscore gives the sample `ERR011347_bin` and the genome `001`;
    splitting at the first gives the sample `ERR011347` only by luck. The known sample names are
    tried longest-first so one that prefixes another cannot claim its rows.
    """
    assert strains._split_label("ERR011347_bin.1", SAMPLES) == ("ERR011347", "bin.1")
    assert strains._split_label("ERR011347_bin_001", SAMPLES) == ("ERR011347", "bin_001")
    assert strains._split_label("S1_g", ("S1", "S1b")) == ("S1", "g")
    assert strains._split_label("S1b_g", ("S1", "S1b")) == ("S1b", "g")
    # With no known sample matching, the genome is reported unknown rather than guessed at.
    assert strains._split_label("unexpected_label", SAMPLES) == ("unexpected_label", "")


def test_an_empty_tracs_table_is_no_rows_rather_than_a_zero_distance(tmp_path):
    """TRACS skips any reference under 25% covered, so "no row" is its commonest honest answer.

    A row invented with distance 0 would be a transmission call; a row with distance None and no
    sites is ABSENT. Both are "no evidence", and only one of them is reported as a shared strain.
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
