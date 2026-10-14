"""The KEGG module definition parser, and what completeness means.

Two kinds of test here, and the second kind is the reason this file is long.

The first kind pins the grammar: that a space is AND, a comma is OR, a plus is a complex, a minus is
non-essential, and that a leading minus marks a whole step non-essential. Those are small and
obvious.

The second kind pins **facts about the vendored module table** - how many modules have one top-level
step, how many contain a reaction with no KO, that M00851 is a flat OR of 21 KOs. Those look like
tests of someone else's data, and they are not: `pathways` reports a flag per module for each of
them, so a reader discounts a completeness figure *because* of these numbers. If a new snapshot
changes them, the flags mean something different and the change has to be noticed. See
src/metawrap2/data/kegg/SNAPSHOT.md.
"""

from __future__ import annotations

import pytest

from metawrap2 import kegg

# ── grammar ──────────────────────────────────────────────────────────────────────────────────


def test_a_comma_binds_tighter_than_a_space():
    """``K1 K2,K3`` is "K1, then either K2 or K3" - the reading the definitions are written for.

    Getting this backwards would read it as "either K1 and K2, or K3", which turns a two-step
    pathway into a one-step one and makes a genome carrying only K3 look complete.
    """
    tree = kegg.parse_definition("K00001 K00002,K00003")
    assert len(kegg.top_level_steps(tree)) == 2
    routes = sorted(kegg._required_sets(tree), key=sorted)
    assert routes == [frozenset({"K00001", "K00002"}), frozenset({"K00001", "K00003"})]

    # ``K11023 -K11024 K11025``: a leading minus marks a whole step non-essential. Seven real
    # modules are written this way, and reading it as a binary operator rejected all seven.
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
    """Stepwise says 0, pathwise says 2/3. Reporting only pathwise hides that nothing works.

    This is the clearest case for emitting both: a two-of-three-subunit complex scores 0.67 pathwise
    and 0.0 stepwise, and the biology is that the complex is absent.
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
    """No silent skips: a definition this cannot read is a module that would score 0 for everyone."""
    assert len(modules) == 573
