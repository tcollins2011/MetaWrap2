"""Tests for pangenome partition counts and the low-persistent-share check.

The fixtures are PPanGGOLiN 2.3.2 output for cluster_0027 of a three-sample CAMI2 cohort (a
cluster whose genomes are not one species), plus three rows of `dereplicate`'s cluster table.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import pangenome as pan

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "pangenome")


# ── the partitioning, from the right file ────────────────────────────────────────────────────


def test_families_are_counted_from_the_partition_files():
    """Families are counted from the partition files. `Persistent_families` in the per-genome
    table is per genome, so summing it would count shared families once per genome."""
    summary = pan.parse_partitions(FIXTURES, "cluster_0027", reliable=False)
    assert summary is not None
    assert summary.families[pan.PERSISTENT] == 1509
    assert summary.families[pan.SHELL] == 3406
    assert summary.families[pan.CLOUD] == 2446
    assert summary.total_families == 7361
    # PPanGGOLiN's statistics file starts with two `#` settings lines, which a plain DictReader
    # would take as the header.
    assert summary.genomes == 3
    assert summary.mean_genes == pytest.approx((4347 + 2394 + 2901) / 3)
    genes = pan.read_genome_statistics(os.path.join(FIXTURES, "genomes_statistics.tsv"))
    assert genes["metawrap_50_10_bins__bin_017"] == pytest.approx(2394)


# ── what a low persistent share says about the clustering ────────────────────────────────────


def test_a_low_persistent_share_is_evidence_about_the_cluster():
    """A cluster with a low persistent share is flagged as suspect.

    Across seven clusters, six were 68-96% persistent and this one 21%.
    """
    summary = pan.parse_partitions(FIXTURES, "cluster_0027", reliable=False)
    assert summary is not None
    shares = summary.shares
    assert shares[pan.PERSISTENT] == pytest.approx(20.5, abs=0.5)
    assert pan.suspect_groups([summary]) == [summary]
    assert "Check the grouping, starting with" in pan.SUSPECT_CLUSTER_CAVEAT
