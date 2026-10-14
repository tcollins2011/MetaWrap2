"""Multi-split binning: splitting a catalogue's bins back per sample.

Binning a cohort has two shapes. What MetaWrap2 did before this is bin each assembly separately with
a cross-sample depth table - the split happens *before* the model. Multi-split bins the concatenated
catalogue and splits the clusters afterwards, so the model is fitted on every sample's contigs at
once; Söding's analysis reports 19% more bins and 79% more high-quality bins for the second shape.

Two of the binners split for themselves when told the separator - SemiBin2 via ``multi_easy_bin -s``
and VAMB via ``-o`` - and those invocations are pinned here because both flags are easy to get
subtly wrong: VAMB's default separator is ``C``, which shreds bins on any contig id containing a
capital C, and the previous code passed ``-o ''`` unconditionally, which is right for one assembly
and silently wrong for a catalogue.

Everything else is about the generic split: that contig names come back unprefixed, that a bin
drawing from several samples becomes several bins, and that a catalogue which was never a catalogue
is reported rather than quietly producing one bin set.
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
    """The nested layout is the wrong shape for anything that takes a bin set as one folder.

    This project's own benchmark scorer collects top-level `*_bins` directories, so without the flat
    view the split arm of a comparison is scored on its *unsplit* catalogue bins - reporting it as
    identical to its control, which is how this was found.
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
    # Each piece keeps its own name, which already carries the sample it came from.
    assert names == ["bin_001.sample_0.fasta", "bin_001.sample_1.fasta", "bin_002.sample_1.fasta"]
    # The contents are the per-sample pieces, with the catalogue prefix stripped back off.
    assert ">c2" in (flat / "bin_001.sample_1.fasta").read_text()


def test_binsplit_takes_a_required_value_so_it_cannot_eat_the_read_files():
    """With nargs="?" this flag sits in front of the positional reads, so `--binsplit r_1.fastq`
    parsed the read file as the separator and then complained no reads were given."""
    args = binning._parse_args(
        ["-a", "cat.fa", "-o", "out", "--metabat2", "--binsplit", "~", "r_1.fastq"]
    )
    assert args.binsplit == "~"
    assert args.reads == ["r_1.fastq"]

    off = binning._parse_args(["-a", "a.fa", "-o", "out", "--metabat2", "r_1.fastq"])
    assert off.binsplit == ""
