"""The two read filters that feed `reassemble_bins`, and the three ways they read a SAM wrongly.

Neither script had a test, which is how all three survived: both decided whether a record had aligned
from the presence of an ``NM`` tag rather than from the FLAG, neither excluded secondary or
supplementary alignments, and a missing ``NM`` scored as 0 - a perfect match - rather than as unknown.

Every SAM record below is written out in full rather than built by a helper, because what is being
tested is the reading of specific columns and flag bits, and a helper that set them would be the same
assumption twice.
"""

from __future__ import annotations

import os

from metawrap2.scripts import filter_nanopore_reads_for_bin_reassembly as nano
from metawrap2.scripts import filter_reads_for_bin_reassembly as short

# QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL tags...
_COLS = "%s\t%d\tc1\t100\t60\t4M\t=\t200\t104\t%s\t%s"


def _bins(tmp_path):
    """Two bins, one contig each, named the way every binner's output is named."""
    (tmp_path / "bins").mkdir()
    (tmp_path / "bins" / "bin_001.fasta").write_text(">c1 len=5000\nACGTACGTAC\n")
    (tmp_path / "bins" / "bin_002.fasta").write_text(">c2 len=4000\nTTTTTTTTTT\n")
    return str(tmp_path / "bins")


def _read(path):
    with open(path) as fh:
        return fh.read()


def test_a_supplementary_alignment_does_not_displace_the_mate_it_is_filed_next_to(tmp_path):
    """The failure is not that a junk record is emitted - it is that a real pair disappears.

    A supplementary record carries the first-mate bit, so the old pairing took it for the forward
    mate, dropped the real one, and then yielded it against the next reverse mate: two different
    reads, recruited on where neither of them aligned. Here that puts the forward mate on c2 and the
    reverse on c1, which disagree, so the pair is silently discarded - and the hard-clipped SEQ would
    have gone into the FASTQ had they agreed.
    """
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:0",
        "r9\t%d\tc2\t50\t60\t2H2M\t=\t60\t12\tGT\tJJ\tNM:i:0" % (0x800 | 0x40),
        _COLS % ("r1", 0x80, "ACGT", "JJJJ") + "\tNM:i:0",
    ]
    bins = _bins(tmp_path)
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), bins, out, strict=2, permissive=5)
    assert counts["recruited"] == 1
    assert _read(os.path.join(out, "bin_001.strict_1.fastq")) == "@r1/1\nACGT\n+\nJJJJ\n"
    assert not os.path.exists(os.path.join(out, "bin_002.strict_1.fastq"))

    # Alignment status comes off the FLAG, not the presence of NM: an aligner that writes no NM
    # used to recruit nothing, so every bin reassembled from no reads.
    sam = [_COLS % ("r1", 0x40, "ACGT", "JJJJ"), _COLS % ("r1", 0x80, "ACGT", "JJJJ")]
    counts = short.filter_reads(iter(sam), bins, str(tmp_path / "no_nm"), strict=2, permissive=5)
    assert counts["recruited"] == 1 and counts["permissive"] == 1

    # A missing NM on an aligned mate is unknown, not 0 (a perfect match): permissive, never strict.
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:5",
        _COLS % ("r1", 0x80, "ACGT", "JJJJ"),
    ]
    counts = short.filter_reads(iter(sam), bins, str(tmp_path / "unk"), strict=6, permissive=99)
    assert counts["recruited"] == 1 and counts["permissive"] == 1 and counts["strict"] == 0


def test_a_long_read_secondary_alignment_is_not_turned_back_into_a_fastq_record(tmp_path):
    """The one that was reaching SPAdes. ``minimap2 -ax map-ont`` emits secondary alignments by
    default, writes their SEQ as ``*``, and puts NM on them - so the old NM-presence test let each
    one through as a one-base read. A long read crossing a repeat produces several.
    """
    sam = [
        "r1\t0\tc1\t100\t60\t4M\t*\t0\t0\tACGT\tJJJJ\tNM:i:30",
        "r1\t%d\tc1\t900\t0\t4M\t*\t0\t0\t*\t*\tNM:i:40" % 0x100,
        "r1\t%d\tc1\t950\t60\t2H2M\t*\t0\t0\tGT\tJJ\tNM:i:1" % 0x800,
    ]
    bins = _bins(tmp_path)
    out = str(tmp_path / "out")
    counts = nano.filter_reads(iter(sam), bins, out)
    assert counts["recruited"] == 1
    assert _read(os.path.join(out, "bin_001.nanopore.fastq")) == "@r1/1\nACGT\n+\nJJJJ\n"

    # An unmapped long read placed at a contig: RNAME passes the reference check, only FLAG says no.
    sam = ["r1\t%d\tc1\t100\t0\t*\t*\t0\t0\tACGT\tJJJJ" % 0x4]
    assert nano.filter_reads(iter(sam), bins, str(tmp_path / "unmapped"))["recruited"] == 0
