"""Pangenomes: what the unit is, which file holds the answer, and what a low persistent share means.

The fixtures are real PPanGGOLiN 2.3.2 output over cluster_0027 of a three-sample CAMI2 cohort -
the cluster whose genomes turned out not to be one species, which is the interesting case - plus
three real rows of `dereplicate`'s cluster table.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import pangenome as pan

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "pangenome")


# ── the partitioning, from the right file ────────────────────────────────────────────────────


def test_families_are_counted_from_the_partition_files():
    """The per-genome table's `Persistent_families` is per genome, so summing it over a group
    multiplies every shared family by the number of genomes that carry it."""
    summary = pan.parse_partitions(FIXTURES, "cluster_0027", reliable=False)
    assert summary is not None
    assert summary.families[pan.PERSISTENT] == 1509
    assert summary.families[pan.SHELL] == 3406
    assert summary.families[pan.CLOUD] == 2446
    assert summary.total_families == 7361
    # PPanGGOLiN's statistics file opens with two `#` settings lines; handed to a plain DictReader
    # it keys every row by `#soft_core=0.95` and the pangenome reads as having no genes.
    assert summary.genomes == 3
    assert summary.mean_genes == pytest.approx((4347 + 2394 + 2901) / 3)
    genes = pan.read_genome_statistics(os.path.join(FIXTURES, "genomes_statistics.tsv"))
    assert genes["metawrap_50_10_bins__bin_017"] == pytest.approx(2394)


# ── what a low persistent share says about the clustering ────────────────────────────────────


def test_a_low_persistent_share_is_evidence_about_the_cluster():
    """Genomes of one species share most of their gene content. These do not.

    Measured across seven real clusters: six at 68-96% persistent and this one at 21%.
    """
    summary = pan.parse_partitions(FIXTURES, "cluster_0027", reliable=False)
    assert summary is not None
    shares = summary.shares
    assert shares[pan.PERSISTENT] == pytest.approx(20.5, abs=0.5)
    assert pan.suspect_groups([summary]) == [summary]
    assert "grouping is the thing to check" in pan.SUSPECT_CLUSTER_CAVEAT
