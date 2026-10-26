"""Tests for the one bin-naming convention: bin_001.fasta, for every binner.

Binners name output differently (metaBAT2 counts from 1, MaxBin2 writes zero-padded
``.001.fasta``, CONCOCT uses cluster ids). Checked: one convention, counting from 1, with padding so
filename order matches bin order (refined bins are numbered in sorted filename order).
"""

from __future__ import annotations

from metawrap2.modules import binning


def test_renumber_orders_numerically_not_lexicographically(tmp_path):
    """metaBAT2's bin.2 is numbered before bin.10."""
    bins = tmp_path / "bins"
    bins.mkdir()
    for number in (1, 2, 10):
        (bins / ("bin.%d.fa" % number)).write_text(">contig_from_%d\nACGT\n" % number)
    binning.renumber_bins(str(bins))
    assert (bins / "bin_001.fasta").read_text().startswith(">contig_from_1")
    assert (bins / "bin_002.fasta").read_text().startswith(">contig_from_2")
    assert (bins / "bin_003.fasta").read_text().startswith(">contig_from_10")

    # CONCOCT clusters count from 0 and bins from 1. Non-contiguous cluster ids are kept (plus
    # one), so the number identifies the CONCOCT cluster.
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>c2\nTTTT\n>c3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\nc2,1\nc3,5\n")
    out = tmp_path / "concoct_bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert (out / "bin_001.fasta").read_text() == ">c1\nACGT\n"
    assert (out / "bin_002.fasta").read_text() == ">c2\nTTTT\n"
    assert (out / "bin_006.fasta").read_text() == ">c3\nGGGG\n"
