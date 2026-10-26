"""Parsing IQ-TREE support values, model and alignment shape.

Fixtures are an IQ-TREE 3.1.4 run over a PPanGGOLiN core-gene alignment of four MAGs from a CAMI2
cohort, including a tree with `100/100` support that an earlier parser counted as zero branches.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import phylogeny as phylo

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "phylogeny")
TREEFILE = os.path.join(FIXTURES, "core.treefile")
REPORT = os.path.join(FIXTURES, "core.iqtree")
ALIGNMENT = os.path.join(FIXTURES, "core_genome_alignment.aln")


# ── support values ───────────────────────────────────────────────────────────────────────────


def test_combined_support_values_are_parsed_and_the_bootstrap_is_the_one_read(tmp_path):
    """Combined `SH-aLRT/UFBoot` labels (e.g. `)97.5/100:`, written with --alrt) are parsed.

    A single-number pattern matched none of them, and a four-taxon tree with one 100/100 branch was
    counted as having zero branches.
    """
    support = phylo.support_summary(TREEFILE)
    assert support.branches == 1
    assert support.supported == 1
    assert support.fraction_supported == pytest.approx(1.0)

    # The UFBoot value decides, not the SH-aLRT: SUPPORTED_BRANCH is calibrated for UFBoot.
    path = tmp_path / "t.treefile"
    # SH-aLRT 100 but UFBoot 80: not supported at this threshold.
    path.write_text("(a:0.1,(b:0.2,c:0.1)100/80:0.3,d:0.2);")
    support = phylo.support_summary(str(path))
    assert support.branches == 1
    assert support.supported == 0


# ── the model ────────────────────────────────────────────────────────────────────────────────


def test_the_chosen_model_is_read_out_of_iqtrees_report():
    """The model IQ-TREE chose is read from its report and recorded in the summary row."""
    assert phylo.model_from_report(REPORT) == "GTR+F+I"
    alignment = phylo.read_alignment(ALIGNMENT)
    assert (alignment.size, alignment.columns) == (4, 400)
    row = phylo.summary_row("t", alignment, "GTR+F+I", phylo.support_summary(TREEFILE))
    assert row["model"] == "GTR+F+I"
    assert row["support_threshold"] == "95"
    assert row["branches_supported"] == "1"
