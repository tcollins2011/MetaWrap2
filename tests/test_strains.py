"""Tests for strain-level comparison: the breadth floor and the units of inStrain's identities.

Both fixtures are inStrain 1.10.0 output over three human gut metagenomes and this repo's 26
refined bins (see tests/fixtures/strains/SOURCE.md). The first run of this module reported all 42
comparisons as "not in both" because of two schema mistakes.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import strains

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "strains")
COMPARE = os.path.join(FIXTURES, "genomeWide_compare.tsv")


# ── units ────────────────────────────────────────────────────────────────────────────────────


def test_identities_are_stored_as_fractions_and_parsed_as_percentages():
    """inStrain's table says 0.99999945; its documentation and paper say 99.999%.

    Parsed as written and compared with the published threshold, every comparison in this fixture
    would be a different strain over too little of the genome, as the first run reported.
    """
    comparisons = {c.genome: c for c in strains.parse_comparison(COMPARE)}
    bin1 = comparisons["bin.1"]
    assert bin1.pop_ani == pytest.approx(99.9999454, abs=1e-5)
    assert bin1.percent_compared == pytest.approx(87.33, abs=0.01)
    # Above the threshold and above the breadth floor, so a real same-strain call.
    assert bin1.verdict is strains.StrainVerdict.SAME
    # `percent_compared`, not `percent_genome_compared`. A wrong name parses to None, which reads
    # as "not in both".
    assert strains.COMPARE_COLUMNS["percent_compared"] == "percent_compared"


# ── the breadth guard ────────────────────────────────────────────────────────────────────────


def test_a_perfect_popani_over_a_sliver_of_a_genome_is_undetermined():
    """High popANI over a small fraction of the genome is UNDETERMINED, not a shared strain.

    In this fixture 13 comparisons have popANI >= 99.999% over 0.1-46.8% of the genome. Without
    the breadth floor there would be 17 shared strains instead of 4.
    """
    comparisons = strains.parse_comparison(COMPARE)
    high_identity = [c for c in comparisons if (c.pop_ani or 0) >= strains.SAME_STRAIN_POPANI]
    assert len(high_identity) == 17
    supported = [c for c in high_identity if c.verdict is strains.StrainVerdict.SAME]
    refused = [c for c in high_identity if c.verdict is strains.StrainVerdict.UNDETERMINED]
    assert len(supported) == 4
    assert len(refused) == 13
    assert max(c.percent_compared for c in refused) < strains.MIN_PERCENT_COMPARED


# ── TRACS, the second backend ────────────────────────────────────────────────────────────────

TRACS = os.path.join(FIXTURES, "tracs_transmission_distances.csv")
SAMPLES = ("ERR011347", "ERR011348")


def test_the_genome_cannot_come_from_the_msa_column():
    """TRACS truncates the reference name at the first dot: `combinedbin` for both bin.1 and bin.4.

    Reading that column gave the genome `bin` for every row, merging two genomes.
    """
    with open(TRACS) as fh:
        rows = list(__import__("csv").DictReader(fh))
    assert {row["MSA file"] for row in rows} == {"combinedbin"}
    # Recovered from the sample labels instead.
    assert {row.genome for row in strains.parse_tracs(TRACS, samples=SAMPLES)} == {"bin.1", "bin.4"}
    by_genome = {row.genome: row for row in strains.parse_tracs(TRACS, samples=SAMPLES)}
    assert by_genome["bin.1"].snp_distance == pytest.approx(5.0)
    assert by_genome["bin.1"].sites_considered == 2006904
