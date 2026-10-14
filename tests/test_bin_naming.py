"""Tests for the one bin-naming convention: bin_001.fasta, for every binner.

Each binner names its output its own way - metaBAT2 counts from 1, MaxBin2 writes zero-padded
``.001.fasta``, CONCOCT emits its own cluster ids - so without normalising, the same number means a
different thing in each folder. The properties worth pinning are that the convention is uniform,
that it counts from 1, and that padding makes filename order and bin order the same thing (which
matters because refined bins are numbered in sorted filename order).
"""

from __future__ import annotations

from metawrap2.modules import binning


def test_renumber_orders_numerically_not_lexicographically(tmp_path):
    """metaBAT2's bin.2 must stay ahead of its bin.10, so bin contents keep their order."""
    bins = tmp_path / "bins"
    bins.mkdir()
    for number in (1, 2, 10):
        (bins / ("bin.%d.fa" % number)).write_text(">contig_from_%d\nACGT\n" % number)
    binning.renumber_bins(str(bins))
    assert (bins / "bin_001.fasta").read_text().startswith(">contig_from_1")
    assert (bins / "bin_002.fasta").read_text().startswith(">contig_from_2")
    assert (bins / "bin_003.fasta").read_text().startswith(">contig_from_10")

    # CONCOCT numbers clusters from 0; bins are numbered from 1 like every other binner's. Cluster
    # ids are not necessarily contiguous and are not renumbered, so the number still identifies
    # which CONCOCT cluster the bin came from.
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>c2\nTTTT\n>c3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\nc2,1\nc3,5\n")
    out = tmp_path / "concoct_bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert (out / "bin_001.fasta").read_text() == ">c1\nACGT\n"
    assert (out / "bin_002.fasta").read_text() == ">c2\nTTTT\n"
    assert (out / "bin_006.fasta").read_text() == ">c3\nGGGG\n"
