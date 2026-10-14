"""CAZyme calls from dbCAN: the agreement between its three methods, and the units it reports in.

Every test here is about a way the output misleads if it is flattened into one "CAZy family" column:

* dbCAN runs three independent methods and a third of the calls on a real bin come from one of them;
* dbCAN-sub reports subfamilies where the others report families, so naive comparison makes every
  subfamily call look like a disagreement;
* a CBM binds a carbohydrate and need not act on one, so summing it with the glycoside hydrolases
  inflates a count of degradative enzymes.

The fixtures are real dbCAN 5.2.9 output - see tests/fixtures/cazy/SOURCE.md.
"""

from __future__ import annotations

import os

from metawrap2 import cazy

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "cazy")
OVERVIEW = os.path.join(FIXTURES, "overview.tsv")


# ── parsing real output ──────────────────────────────────────────────────────────────────────


def test_a_real_bins_calls_are_parsed():
    calls = cazy.parse_overview(OVERVIEW, "bin_006")
    assert len(calls) == 65
    assert all(call.genome == "bin_006" for call in calls)
    # 29 of 65 have two or more methods behind them. That ratio is the reason the column exists.
    assert sum(1 for call in calls if call.is_confident) == 29
    # A binding module binds a carbohydrate and need not act on one, so it is never catalytic.
    cbm_only = [call for call in calls if call.classes and set(call.classes) == {"CBM"}]
    assert cbm_only, "the fixture has CBM-only proteins"
    assert all(not call.is_catalytic for call in cbm_only)
    # Real numbers from this bin: 24 glycoside hydrolases in 18 families, 6 of them confident.
    by_class = {s.cazy_class: s for s in cazy.summarise_classes(calls, "bin_006")}
    assert (by_class["GH"].genes, by_class["GH"].families) == (24, 18)
    assert by_class["GH"].confident_genes == 6


def test_coordinates_are_stripped_but_the_family_is_kept():
    """A subfamily call counts as agreement with its family, and coordinates never block a match."""
    calls = {call.gene: call for call in cazy.parse_overview(OVERVIEW, "bin_006")}
    call = calls["k141_11967_length_9582_cov_7.0000_3"]
    # dbCAN wrote GT119(45-389) and GT119_e82(39-389); the coordinates are per method, the family
    # is not, so comparing the strings would never find agreement.
    assert "GT119" in call.families
    assert call.families_agreed == ("GT119",)
    assert all("(" not in family for family in call.families)
    # dbCAN-sub exists to add the `_e82` suffix; treating it as a different call defeats that.
    assert cazy._agreed(("GH13",), ("GH13",), ("GH13_e82",)) == ("GH13",)
