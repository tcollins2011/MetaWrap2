"""Viral contigs: quality tier and interpretation of CheckV completeness.

The fixtures are five contigs of output from geNomad 1.12.0 and CheckV 1.1.1 on a CAMI2 assembly.
Real output is used because both schemas were initially misread: CheckV's confidence is in
parentheses in ``completeness_method``, and CheckV's ``provirus`` column is a literal Yes/No
unrelated to geNomad's ``topology``.

Covered:

* ``HMM-based (lower-bound)`` means "at least this complete", not an estimate. 35 of 161 contigs
  on a real assembly have it, mostly viruses with no close relative.
* The tier is taken from CheckV, which applies MIUViG including terminal-repeat evidence for
  Complete genomes, and is not re-derived from the number.
* A contig CheckV did not score is undetermined, not low quality (16 of the 161).
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import viruses

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "viruses")
GENOMAD = os.path.join(FIXTURES, "virus_summary.tsv")
CHECKV = os.path.join(FIXTURES, "quality_summary.tsv")


@pytest.fixture(scope="module")
def merged():
    return viruses.merge("assembly", viruses.parse_genomad(GENOMAD), viruses.parse_checkv(CHECKV))


# ── parsing real output ──────────────────────────────────────────────────────────────────────


def test_a_lower_bound_is_not_treated_as_an_estimate(merged):
    """An HMM-based completeness figure is treated as a lower bound."""
    assert len(merged) == 5
    bounded = [c for c in merged if c.completeness_confidence == viruses.LOWER_BOUND]
    assert bounded, "the fixture has to contain one for this to test anything"
    assert all(not c.completeness_is_reliable for c in bounded)
    assert all(c.completeness is not None for c in bounded), "the number is still reported"

    # The tier is CheckV's own determination, not re-derived from the completeness number.
    by_contig = {c.contig: viruses.tier_of(c) for c in merged}
    assert by_contig["k141_62224_length_6665_cov_7.3708"] is viruses.Tier.UNDETERMINED
    assert by_contig["k141_21474_length_8868_cov_5.0000"] is viruses.Tier.LOW_QUALITY
    # A contig CheckV never scored is undetermined, not low quality (16 of 161 on a real assembly).
    unscored = [c for c in merged if c.completeness is None]
    assert unscored
    for contig in unscored:
        assert viruses.tier_of(contig) is viruses.Tier.UNDETERMINED
        assert not contig.completeness_is_reliable


def test_the_provirus_flag_comes_from_checkv_not_from_topology():
    """The provirus flag is CheckV's Yes/No column. geNomad's ``topology`` describes terminal
    repeats; reading it as a provirus marker labelled every linear contig a provirus."""
    genomad = {"c1": {"seq_name": "c1", "length": "5000", "topology": "No terminal repeats"}}
    not_provirus = viruses.merge("g", genomad, {})[0]
    assert not not_provirus.provirus

    checkv = {"c1": {"contig_id": "c1", "provirus": "Yes", "proviral_length": "3000"}}
    is_provirus = viruses.merge("g", genomad, checkv)[0]
    assert is_provirus.provirus
    assert is_provirus.trimmed_length == 3000


# ── vOTU clustering ──────────────────────────────────────────────────────────────────────────


def _edge(ref, query, ani, af):
    from metawrap2.modules.dereplicate import Edge

    return Edge(ref, query, ani, af, af)


def test_the_viral_coverage_threshold_keeps_a_shared_region_from_merging():
    """95% ANI over 50% of the sequence is a shared region, not one population. Genome
    thresholds (95/10) would merge these; viral thresholds (95/85) do not."""
    lengths = {"A": 40000, "B": 38000}
    assert len(viruses.build_votus(["A", "B"], [_edge("A", "B", 97.0, 50.0)], lengths)) == 2
    assert len(viruses.build_votus(["A", "B"], [_edge("A", "B", 97.0, 90.0)], lengths)) == 1
