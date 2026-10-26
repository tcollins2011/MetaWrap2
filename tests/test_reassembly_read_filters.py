"""Tests for the short- and long-read filters that feed `reassemble_bins`.

Covers three past bugs: alignment status was taken from the presence of an ``NM`` tag instead of
the FLAG, secondary and supplementary alignments were not excluded, and a missing ``NM`` scored as
0 (a perfect match) instead of unknown.

SAM records are written out in full, not built by a helper, so the columns and flag bits under test
are explicit.
"""

from __future__ import annotations

import os

from metawrap2.scripts import filter_nanopore_reads_for_bin_reassembly as nano
from metawrap2.scripts import filter_reads_for_bin_reassembly as short

# QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL tags...
_COLS = "%s\t%d\tc1\t100\t60\t4M\t=\t200\t104\t%s\t%s"


def _bins(tmp_path):
    """Two bins, one contig each, named as binner output is named."""
    (tmp_path / "bins").mkdir()
    (tmp_path / "bins" / "bin_001.fasta").write_text(">c1 len=5000\nACGTACGTAC\n")
    (tmp_path / "bins" / "bin_002.fasta").write_text(">c2 len=4000\nTTTTTTTTTT\n")
    return str(tmp_path / "bins")


def _read(path):
    with open(path) as fh:
        return fh.read()


def test_a_supplementary_alignment_does_not_displace_the_mate_it_is_filed_next_to(tmp_path):
    """A supplementary record between two mates does not replace the forward mate.

    The supplementary record carries the first-mate bit. Previously it replaced the real forward
    mate and was paired with the next reverse mate; here the two disagree on bin (c2 and c1), so
    the pair was dropped.
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

    # Alignment status comes from the FLAG; with the NM check, an aligner that writes no NM
    # recruited no reads.
    sam = [_COLS % ("r1", 0x40, "ACGT", "JJJJ"), _COLS % ("r1", 0x80, "ACGT", "JJJJ")]
    counts = short.filter_reads(iter(sam), bins, str(tmp_path / "no_nm"), strict=2, permissive=5)
    assert counts["recruited"] == 1 and counts["permissive"] == 1

    # A missing NM on an aligned mate is unknown, so the pair is permissive, never strict.
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:5",
        _COLS % ("r1", 0x80, "ACGT", "JJJJ"),
    ]
    counts = short.filter_reads(iter(sam), bins, str(tmp_path / "unk"), strict=6, permissive=99)
    assert counts["recruited"] == 1 and counts["permissive"] == 1 and counts["strict"] == 0


def test_a_long_read_secondary_alignment_is_not_turned_back_into_a_fastq_record(tmp_path):
    """Long-read secondary and supplementary alignments are not written to the FASTQ.

    ``minimap2 -ax map-ont`` emits secondary alignments by default, with SEQ ``*`` and an NM tag,
    so the NM-presence check passed each one to SPAdes as a one-base read.
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

    # An unmapped read placed at a contig has an RNAME; only the FLAG marks it unmapped.
    sam = ["r1\t%d\tc1\t100\t0\t*\t*\t0\t0\tACGT\tJJJJ" % 0x4]
    assert nano.filter_reads(iter(sam), bins, str(tmp_path / "unmapped"))["recruited"] == 0
