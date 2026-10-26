"""Tests for the KEGG module definition parser and completeness scoring.

Grammar tests: a space is AND, a comma is OR, a plus is a complex, a minus is non-essential, and a
leading minus marks a whole step non-essential.

Snapshot tests pin facts about the vendored module table (for example, how many modules have one
top-level step, or that M00851 is a flat OR of 21 KOs). `pathways` reports per-module flags based
on these, so a new snapshot that changes them must be reviewed. See
src/metawrap2/data/kegg/SNAPSHOT.md.
"""

from __future__ import annotations

import pytest

from metawrap2 import kegg

# ── grammar ──────────────────────────────────────────────────────────────────────────────────


def test_a_comma_binds_tighter_than_a_space():
    """``K1 K2,K3`` parses as "K1, then either K2 or K3".

    The other reading, "either K1 and K2, or K3", would make a genome with only K3 look complete.
    """
    tree = kegg.parse_definition("K00001 K00002,K00003")
    assert len(kegg.top_level_steps(tree)) == 2
    routes = sorted(kegg._required_sets(tree), key=sorted)
    assert routes == [frozenset({"K00001", "K00002"}), frozenset({"K00001", "K00003"})]

    # ``K11023 -K11024 K11025``: a leading minus marks a whole step non-essential. Seven modules
    # use this form; parsing the minus as a binary operator rejects them.
    tree = kegg.parse_definition("K11023 -K11024 K11025")
    assert len(kegg.top_level_steps(tree)) == 3
    assert kegg._required_sets(tree) == [frozenset({"K11023", "K11025"})]


# ── scoring ──────────────────────────────────────────────────────────────────────────────────


def _module(definition: str, module_class: str = "Pathway modules; test") -> kegg.Module:
    return kegg.Module(
        module="M99999",
        name="test",
        module_class=module_class,
        definition=definition,
        tree=kegg.parse_definition(definition),
    )


def test_a_half_assembled_complex_is_a_missing_step_not_a_partial_one():
    """A complex with two of three subunits scores 2/3 pathwise and 0 stepwise.

    The complex is absent, which only the stepwise score shows; both scores are reported.
    """
    module = _module("K00001+K00002+K00003")
    score = kegg.score_module(module, {"K00001", "K00002"})
    assert score.completeness_pathwise == pytest.approx(2 / 3)
    assert score.completeness_stepwise == 0.0


# ── the vendored snapshot ────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def modules():
    return kegg.load_modules()


def test_every_definition_in_the_snapshot_parses(modules):
    """All 573 snapshot definitions parse; an unparsed module would score 0 for every genome."""
    assert len(modules) == 573
