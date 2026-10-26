"""Multi-split binning: splitting a catalogue's bins back per sample.

Per-sample binning bins each assembly separately with a cross-sample depth table. Multi-split bins
the concatenated catalogue and splits the clusters afterwards, so the model sees every sample's
contigs; Söding's analysis reports 19% more bins and 79% more high-quality bins this way.

SemiBin2 (``multi_easy_bin -s``) and VAMB (``-o``) split for themselves when given the separator.
VAMB's default separator is ``C``, which splits any contig id containing a capital C, and
``-o ''`` is only correct for a single assembly.

The generic split is checked for unprefixed contig names, for a bin from several samples becoming
several bins, and for reporting input that is not a catalogue.
"""

from __future__ import annotations

from metawrap2 import binsplit
from metawrap2.modules import binning


def _write_bin(path, names_and_seqs):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, seq) for name, seq in names_and_seqs)
    return str(path)


# ── splitting a whole bin set ────────────────────────────────────────────────────────────────


def test_split_bins_are_also_offered_as_one_flat_bin_set(tmp_path):
    """Split bins are also linked into one flat folder, for tools that take a bin set as a folder.

    The benchmark scorer collects top-level `*_bins` directories; without the flat view it scored
    the unsplit catalogue bins.
    """
    bins = tmp_path / "metabat2_bins"
    bins.mkdir()
    (bins / "bin_001.fasta").write_text(
        ">sample_0~c1\nACGTACGTAC\n>sample_1~c2\nTTTTTTTTTT\n>sample_0~c3\nGGGGGGGGGG\n"
    )
    (bins / "bin_002.fasta").write_text(">sample_1~c9\nAAAACCCCGG\n")
    nested = tmp_path / "metabat2_bins_split"
    report = binsplit.split_bin_set(str(bins), str(nested))
    assert report.bins_in == 2
    assert sorted(report.labels) == ["sample_0", "sample_1"]

    flat = tmp_path / "metabat2_split_bins"
    linked = binsplit.flatten_split_bins(str(nested), str(flat))
    names = sorted(p.name for p in flat.iterdir())
    assert linked == report.pieces_out == 3
    # Each piece keeps its name, which includes the sample.
    assert names == ["bin_001.sample_0.fasta", "bin_001.sample_1.fasta", "bin_002.sample_1.fasta"]
    # The contents are the per-sample pieces with the catalogue prefix removed.
    assert ">c2" in (flat / "bin_001.sample_1.fasta").read_text()


def test_binsplit_takes_a_required_value_so_it_cannot_eat_the_read_files():
    """`--binsplit` requires a value, so it does not consume a following read file (as it did
    with nargs="?")."""
    args = binning._parse_args(
        ["-a", "cat.fa", "-o", "out", "--metabat2", "--binsplit", "~", "r_1.fastq"]
    )
    assert args.binsplit == "~"
    assert args.reads == ["r_1.fastq"]

    off = binning._parse_args(["-a", "a.fa", "-o", "out", "--metabat2", "r_1.fastq"])
    assert off.binsplit == ""
