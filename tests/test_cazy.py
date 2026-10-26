"""CAZyme calls from dbCAN: the agreement between its three methods, and the units it reports in.

Checked:

* dbCAN runs three methods; the number agreeing on each call is kept;
* dbCAN-sub reports subfamilies where the others report families, and a subfamily counts as
  agreement with its family;
* a CBM binds a carbohydrate without necessarily acting on one, so it is not counted as catalytic.

The fixtures are real dbCAN 5.2.9 output - see tests/fixtures/cazy/SOURCE.md.
"""

from __future__ import annotations

import os

from metawrap2 import cazy

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "cazy")
OVERVIEW = os.path.join(FIXTURES, "overview.tsv")
MULTIDOMAIN = os.path.join(FIXTURES, "overview_multidomain.tsv")


# ── parsing real output ──────────────────────────────────────────────────────────────────────


def test_a_real_bins_calls_are_parsed():
    calls = cazy.parse_overview(OVERVIEW, "bin_006")
    assert len(calls) == 65
    assert all(call.genome == "bin_006" for call in calls)
    # 29 of 65 calls are supported by two or more methods.
    assert sum(1 for call in calls if call.is_confident) == 29
    # A CBM-only protein is never catalytic.
    cbm_only = [call for call in calls if call.classes and set(call.classes) == {"CBM"}]
    assert cbm_only, "the fixture has CBM-only proteins"
    assert all(not call.is_catalytic for call in cbm_only)
    # This bin: 24 glycoside hydrolases in 18 families, 6 of them confident.
    by_class = {s.cazy_class: s for s in cazy.summarise_classes(calls, "bin_006")}
    assert (by_class["GH"].genes, by_class["GH"].families) == (24, 18)
    assert by_class["GH"].confident_genes == 6


def test_coordinates_are_stripped_but_the_family_is_kept():
    """A subfamily call counts as agreement with its family, and coordinates never block a match."""
    calls = {call.gene: call for call in cazy.parse_overview(OVERVIEW, "bin_006")}
    call = calls["k141_11967_length_9582_cov_7.0000_3"]
    # dbCAN wrote GT119(45-389) and GT119_e82(39-389); coordinates differ per method, so they are
    # stripped before comparing.
    assert "GT119" in call.families
    assert call.families_agreed == ("GT119",)
    assert all("(" not in family for family in call.families)
    # A dbCAN-sub suffix such as `_e82` still agrees with the family.
    assert cazy._agreed(("GH13",), ("GH13",), ("GH13_e82",)) == ("GH13",)


def test_a_gene_with_no_ec_on_either_domain_has_no_ec():
    """A multi-domain gene with no EC on any domain (dbCAN writes `-|-`) gets an empty `ec`.

    dbCAN joins per-domain fields with `|` and writes `-` for an unannotated domain; 2 of the 11
    rows in the fixture are like this. Whole `-` members are dropped, so the trailing dash of a
    partial EC number (`3.2.1.-`) is kept.
    """
    calls = {c.gene: c for c in cazy.parse_overview(MULTIDOMAIN, "bin_001")}
    assert calls["k141_10594_length_49463_cov_14.0000_27"].ec == ""
    assert calls["k141_10594_length_49463_cov_14.0000_17"].ec == "3.2.1.-:2"
    assert sum(1 for c in calls.values() if c.ec) == 4
