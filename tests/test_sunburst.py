"""The static wheel, which exists because Krona writes no image.

`ktImportText` produces an interactive HTML file and nothing else, so until this module there was no
picture of the functional hierarchy anywhere on disk - not in a report, not in the release notes,
not in a paper figure. These tests pin the parts that would make such a picture wrong rather than
ugly: that ring totals are conserved, that the branches KEGG cannot place are still drawn, and that
a branch keeps its colour between two wheels so they can be compared side by side.
"""

from __future__ import annotations

import os

from metawrap2 import sunburst as S

#: What `hierarchy.write_krona_input` writes: a weight, then one column per level. Two biological
#: branches, one scaffolding branch, and one path that stops at level B.
KRONA_TEXT = (
    "40\tMetabolism\tCarbohydrate metabolism\tGlycolysis / Gluconeogenesis\n"
    "25\tMetabolism\tCarbohydrate metabolism\tCitrate cycle (TCA cycle)\n"
    "15\tMetabolism\tAmino acid metabolism\tLysine biosynthesis\n"
    "12\tNot Included in Pathway or Brite\n"
    "8\tGenetic Information Processing\tTranslation\tRibosome\n"
)


def _write(tmp_path, text=KRONA_TEXT):
    path = tmp_path / "genome.bin_001.all.txt"
    path.write_text(text)
    return str(path)


# ── reading the file Krona is given ──────────────────────────────────────────────────────────


def test_the_wheel_is_drawn_from_the_same_file_krona_is_given(tmp_path):
    """Read back rather than recomputed, so the picture cannot drift from the wheel beside it."""
    rows = S.read_krona_text(_write(tmp_path))
    assert len(rows) == 5
    assert rows[0] == (
        40.0,
        ("Metabolism", "Carbohydrate metabolism", "Glycolysis / Gluconeogenesis"),
    )
    assert rows[3] == (12.0, ("Not Included in Pathway or Brite",))


def test_a_malformed_line_is_skipped_and_a_zero_is_not_a_wedge(tmp_path):
    """A zero weight is a statement that a branch is absent; a malformed line is not a zero."""
    rows = S.read_krona_text(
        _write(tmp_path, "nonsense\tMetabolism\n0\tMetabolism\n5\tMetabolism\n")
    )
    assert rows == [(5.0, ("Metabolism",))]


# ── conserving the total, which is what a percentage on this plot depends on ─────────────────


def test_every_ring_sums_to_the_same_total(tmp_path):
    """A path that stops short of a ring is counted at the depth it reaches.

    Without that, each ring outward would quietly total less than the one inside it and every
    percentage on the outer rings would be an overestimate - the same failure the hierarchy's own
    distribution guards against.
    """
    rows = S.read_krona_text(_write(tmp_path))
    total = sum(weight for weight, _ in rows)
    for depth in (1, 2, 3):
        ring = S.collapse_ring(rows, depth)
        assert sum(weight for _, weight in ring) == total


def test_the_branches_kegg_cannot_place_are_still_drawn(tmp_path):
    """Dropping scaffolding to make a prettier wheel inflates every remaining percentage."""
    ring = dict(S.collapse_ring(S.read_krona_text(_write(tmp_path)), 1))
    assert ("Not Included in Pathway or Brite",) in ring
    assert ring[("Not Included in Pathway or Brite",)] == 12.0


def test_children_are_grouped_under_their_parent_not_sorted_by_weight(tmp_path):
    """A child wedge has to sit inside the parent it belongs to, or the plot is decorative."""
    ring = S.collapse_ring(S.read_krona_text(_write(tmp_path)), 2)
    parents = [key[0] for key, _ in ring]
    # Each parent appears in one contiguous run.
    assert parents == sorted(parents, key=parents.index)
    assert len(set(parents)) == 3


# ── comparability between wheels ─────────────────────────────────────────────────────────────


def test_a_branch_keeps_its_colour_between_two_wheels():
    """Two samples' wheels are meant to be read side by side, so colour cannot follow rank."""
    assert S.colour_for("Metabolism") == S.colour_for("Metabolism")
    assert S.colour_for("Metabolism") != S.colour_for("Human Diseases")
    # An unknown branch still gets a stable colour within one plot.
    table: dict = {}
    first = S.colour_for("Some Future Branch", table)
    assert S.colour_for("Some Future Branch", table) == first


def test_level_a_order_is_fixed_rather_than_taken_from_the_data(tmp_path):
    """KEGG's own order, so the same branch is in the same place in every wheel."""
    rows = S.read_krona_text(_write(tmp_path))
    names = [key[0] for key, _ in S.collapse_ring(rows, 1)]
    assert names.index("Metabolism") < names.index("Genetic Information Processing")
    assert names[-1] == "Not Included in Pathway or Brite"


# ── the file it writes ───────────────────────────────────────────────────────────────────────


def test_a_png_is_written_beside_the_krona_input(tmp_path):
    path = _write(tmp_path)
    out = S.draw_from_krona_text(path)
    assert out == os.path.splitext(path)[0] + ".png"
    assert os.path.getsize(out) > 5000, "a wheel this size should not be a near-empty PNG"


def test_an_empty_distribution_is_no_wheel_rather_than_an_error(tmp_path):
    """A genome with no KO at all is a real result, not a failure of the plot."""
    assert S.draw([], str(tmp_path / "x.png")) is None
    assert S.draw_from_krona_text(_write(tmp_path, "0\tMetabolism\n")) is None


def test_a_long_branch_name_is_wrapped_rather_than_clipped():
    wrapped = S._wrap("Environmental Information Processing")
    assert "\n" in wrapped
    assert wrapped.replace("\n", " ") == "Environmental Information Processing"
