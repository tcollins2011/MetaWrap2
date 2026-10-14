"""Biosynthetic gene clusters: GECCO's table, and the two ways a cluster count misleads.

Both fixtures are real GECCO 0.11.0 output (see tests/fixtures/bgc/SOURCE.md). The empty one is the
one that matters: with ``--force-tsv`` GECCO writes a placeholder row rather than a header alone, and
read literally that is one unclassified cluster of zero genes - a BGC that does not exist.
"""

from __future__ import annotations

import os

from metawrap2 import bgc

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "bgc")
CLUSTERS = os.path.join(FIXTURES, "clusters.tsv")
EMPTY = os.path.join(FIXTURES, "clusters_empty.tsv")

#: The contigs the fixture's clusters sit on, with the lengths their own names carry. Passed to the
#: parser because GECCO does not report them, and without them no cluster can be called truncated.
LENGTHS = {
    "k141_17778_length_9234_cov_12.0000": 9234,
    "k141_23349_length_59489_cov_10.2261": 59489,
    "k141_45765_length_125938_cov_13.0000": 125938,
    "k141_55485_length_23761_cov_10.0000": 23761,
}


def test_the_placeholder_row_of_a_genome_with_no_clusters_is_not_a_cluster():
    """`--force-tsv` writes a row of blanks, which is how "no clusters" used to read as one.

    The flag is still passed, because without it a genome with no clusters writes no file at all and
    "found nothing" becomes indistinguishable from "the run failed".
    """
    assert bgc.parse_clusters(EMPTY, "bin_006", LENGTHS) == []


# ── truncation: the thing a cluster count hides ──────────────────────────────────────────────


def test_clusters_at_a_contig_edge_are_flagged_as_truncated():
    """A cluster is a neighbourhood, so one at a contig edge continues into unassembled sequence."""
    clusters = bgc.parse_clusters(CLUSTERS, "bin_003", LENGTHS)
    assert len(clusters) == 4
    assert (clusters[0].start, clusters[0].end) == (392, 4427)
    # start 392 on a 9,234 bp contig: 391 from the near end, 4,807 from the far one.
    assert clusters[0].distance_to_contig_end == 391
    truncated = [c for c in clusters if c.is_truncated]
    assert len(truncated) == 2
    # Both sit within one gene's worth of sequence of their contig's end.
    assert all(c.distance_to_contig_end is not None for c in truncated)
    assert all(c.distance_to_contig_end < bgc.EDGE_WINDOW for c in truncated)


# ── Unknown is a weak classification, not an unreal cluster ──────────────────────────────────


# ── genes and domains, by this pipeline's own ids ────────────────────────────────────────────


# ── the per-genome summary ───────────────────────────────────────────────────────────────────
