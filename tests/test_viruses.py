"""Viral contigs: the tier, and the three ways a completeness number can mean less than it looks.

The fixtures are real output, five contigs' worth, from running geNomad 1.12.0 and CheckV 1.1.1 over
a CAMI2 assembly. Real rather than invented because both schemas were guessed wrong the first time:
the confidence CheckV reports is inside a parenthesis in ``completeness_method``, and CheckV's
``provirus`` column is a literal Yes/No that has nothing to do with geNomad's ``topology``. Parsing
against invented fixtures would have agreed with the guesses.

What the tests are about:

* **A lower bound is not an estimate.** ``HMM-based (lower-bound)`` means "at least this complete".
  35 of 161 contigs on a real assembly land there, and treating those as estimates under-reports
  completeness for precisely the viruses with no close relative - the novel ones.
* **CheckV's tier is CheckV's.** It applies MIUViG including the terminal-repeat evidence that makes
  a genome Complete rather than merely 100% estimated, so re-deriving a tier from the number could
  only disagree with the tool that produced it.
* **A contig CheckV never scored is undetermined, not low quality.** 16 of the 161 are like that.
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
    """The one that matters: an HMM-based figure is a floor, not a measurement."""
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
    """CheckV's ``provirus`` is a literal Yes/No. geNomad's ``topology`` describes terminal repeats,
    and reading it as a provirus marker labelled every linear contig one."""
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
    """95% ANI over 50% of the sequence is a shared gene, not the same population. The genome
    thresholds (95/10) would merge these; the viral ones (95/85) must not."""
    lengths = {"A": 40000, "B": 38000}
    assert len(viruses.build_votus(["A", "B"], [_edge("A", "B", 97.0, 50.0)], lengths)) == 2
    assert len(viruses.build_votus(["A", "B"], [_edge("A", "B", 97.0, 90.0)], lengths)) == 1
