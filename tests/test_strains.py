"""Strain-level comparison: the breadth guard, and the units inStrain stores its identities in.

Both fixtures are real inStrain 1.10.0 output over three human gut metagenomes and this repo's 26
refined bins (see tests/fixtures/strains/SOURCE.md). They exist because the first real run of this
module reported every one of 42 comparisons as "not in both" - two schema mistakes that each look
like a plausible biological answer.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import strains

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "strains")
COMPARE = os.path.join(FIXTURES, "genomeWide_compare.tsv")


# ── the units, which are the whole trap ──────────────────────────────────────────────────────


def test_identities_are_stored_as_fractions_and_parsed_as_percentages():
    """inStrain's table says 0.99999945; its documentation and paper say 99.999%.

    Parsed literally against the published threshold, every comparison in this fixture is a
    different strain over too little of the genome to tell - which is exactly what the first real
    run reported, with no error anywhere.
    """
    comparisons = {c.genome: c for c in strains.parse_comparison(COMPARE)}
    bin1 = comparisons["bin.1"]
    assert bin1.pop_ani == pytest.approx(99.9999454, abs=1e-5)
    assert bin1.percent_compared == pytest.approx(87.33, abs=0.01)
    # Above the threshold and above the breadth floor, so a real same-strain call.
    assert bin1.verdict is strains.StrainVerdict.SAME
    # `percent_compared`, not `percent_genome_compared` - and the wrong name is not an error: a
    # missing column parses to None, which reads as "not in both", so the whole table came back as
    # an absence rather than as a parse failure.
    assert strains.COMPARE_COLUMNS["percent_compared"] == "percent_compared"


# ── the breadth guard ────────────────────────────────────────────────────────────────────────


def test_a_perfect_popani_over_a_sliver_of_a_genome_is_undetermined():
    """Two samples sharing 2% of a genome at 100% identity is not evidence of a shared strain.

    Measured on this fixture: 13 comparisons have popANI >= 99.999% over 0.1-46.8% of the genome.
    Without the floor the module would claim 17 shared strains where 4 are supportable.
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

    The first version of this parser read that column and returned the genome `bin` for every row -
    two different genomes collapsed onto one name, with nothing to say so.
    """
    with open(TRACS) as fh:
        rows = list(__import__("csv").DictReader(fh))
    assert {row["MSA file"] for row in rows} == {"combinedbin"}
    # Recovered from the sample labels instead.
    assert {row.genome for row in strains.parse_tracs(TRACS, samples=SAMPLES)} == {"bin.1", "bin.4"}
    by_genome = {row.genome: row for row in strains.parse_tracs(TRACS, samples=SAMPLES)}
    assert by_genome["bin.1"].snp_distance == pytest.approx(5.0)
    assert by_genome["bin.1"].sites_considered == 2006904
