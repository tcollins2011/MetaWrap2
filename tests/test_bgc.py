"""Parsing GECCO's cluster table, including truncation and the empty-genome placeholder row.

Both fixtures are GECCO 0.11.0 output (see tests/fixtures/bgc/SOURCE.md). In the empty one,
``--force-tsv`` has written a placeholder row, which must not be read as a cluster of zero genes.
"""

from __future__ import annotations

import os

from metawrap2 import bgc

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "bgc")
CLUSTERS = os.path.join(FIXTURES, "clusters.tsv")
EMPTY = os.path.join(FIXTURES, "clusters_empty.tsv")

#: The fixture's contigs and the lengths in their names. GECCO does not report contig lengths, and
#: truncation cannot be judged without them.
LENGTHS = {
    "k141_17778_length_9234_cov_12.0000": 9234,
    "k141_23349_length_59489_cov_10.2261": 59489,
    "k141_45765_length_125938_cov_13.0000": 125938,
    "k141_55485_length_23761_cov_10.0000": 23761,
}


def test_the_placeholder_row_of_a_genome_with_no_clusters_is_not_a_cluster():
    """The blank placeholder row from `--force-tsv` parses to no clusters."""
    assert bgc.parse_clusters(EMPTY, "bin_006", LENGTHS) == []


# ── truncation ───────────────────────────────────────────────────────────────────────────────


def test_clusters_at_a_contig_edge_are_flagged_as_truncated():
    """Clusters within EDGE_WINDOW of a contig end are flagged as truncated."""
    clusters = bgc.parse_clusters(CLUSTERS, "bin_003", LENGTHS)
    assert len(clusters) == 4
    assert (clusters[0].start, clusters[0].end) == (392, 4427)
    # start 392 on a 9,234 bp contig: 391 from the near end, 4,807 from the far one.
    assert clusters[0].distance_to_contig_end == 391
    truncated = [c for c in clusters if c.is_truncated]
    assert len(truncated) == 2
    # Both are within EDGE_WINDOW of their contig's end.
    assert all(c.distance_to_contig_end is not None for c in truncated)
    assert all(c.distance_to_contig_end < bgc.EDGE_WINDOW for c in truncated)


# ── Unknown type ─────────────────────────────────────────────────────────────────────────────


# ── genes and domains ────────────────────────────────────────────────────────────────────────


# ── the per-genome summary ───────────────────────────────────────────────────────────────────
