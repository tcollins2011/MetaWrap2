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


def test_a_clean_pair_is_recruited_to_its_bin_at_both_stringencies(tmp_path):
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:0",
        _COLS % ("r1", 0x80 | 0x10, "ACGT", "JJJJ") + "\tNM:i:1",
    ]
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=2, permissive=5)
    assert counts == {"pairs": 1, "recruited": 1, "strict": 1, "permissive": 1}
    # Reverse-strand mate comes back in the original orientation, which is what an assembler wants.
    assert _read(os.path.join(out, "bin_001.strict_2.fastq")) == "@r1/2\nACGT\n+\nJJJJ\n"


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
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=2, permissive=5)
    assert counts["recruited"] == 1
    assert _read(os.path.join(out, "bin_001.strict_1.fastq")) == "@r1/1\nACGT\n+\nJJJJ\n"
    assert not os.path.exists(os.path.join(out, "bin_002.strict_1.fastq"))


def test_an_unknown_edit_distance_is_permissive_and_never_strict(tmp_path):
    """0 is a *perfect* match, so defaulting to it promoted an unknown to the top of the ranking.

    The pair below carries 5 mismatches on the mate that reports any, and nothing from the mate that
    does not. Summed as 5 it clears a strict cutoff of 6; the point is that nobody wrote the second
    number down, so strict - which is a claim about how well the read matches - cannot be supported.
    """
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:5",
        _COLS % ("r1", 0x80, "ACGT", "JJJJ"),  # aligned, no NM tag
    ]
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=6, permissive=99)
    assert counts["recruited"] == 1
    assert counts["permissive"] == 1
    assert counts["strict"] == 0
    assert not os.path.exists(os.path.join(out, "bin_001.strict_1.fastq"))


def test_an_unmapped_mate_contributes_no_mismatches_rather_than_zero(tmp_path):
    """An unmapped mate placed at its mate's reference - which the spec allows and bwa does.

    Its RNAME says c1, so it cannot be told from an aligned record by reference name; only the FLAG
    says it did not align. Counting a mismatch count for it at all is the bug, and counting it as 0
    made a half-aligned pair look flawless.
    """
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ") + "\tNM:i:1",
        _COLS % ("r1", 0x80 | 0x4, "ACGT", "JJJJ"),
    ]
    assert (
        short.pair_mismatches(sam[0].split("\t"), sam[1].split("\t")) == 1
    ), "the unmapped mate must not contribute, and must not make the sum unknown either"
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=2, permissive=5)
    assert counts["recruited"] == 1
    assert counts["strict"] == 1


def test_alignment_status_comes_off_the_flag_not_the_nm_tag(tmp_path):
    """An aligner that does not write NM has to work, because this module no longer names one.

    The old test for "did this align" was the presence of NM, which is a claim about bwa. Under it
    this pair recruits nothing, so every bin would reassemble from no reads - a failure that produces
    an empty result rather than an error.
    """
    sam = [
        _COLS % ("r1", 0x40, "ACGT", "JJJJ"),
        _COLS % ("r1", 0x80, "ACGT", "JJJJ"),
    ]
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=2, permissive=5)
    assert counts["recruited"] == 1
    assert counts["permissive"] == 1


def test_both_mates_unmapped_is_not_recruited(tmp_path):
    sam = [
        "r1\t%d\t*\t0\t0\t*\t*\t0\t0\tACGT\tJJJJ" % (0x4 | 0x40),
        "r1\t%d\t*\t0\t0\t*\t*\t0\t0\tACGT\tJJJJ" % (0x4 | 0x80),
    ]
    out = str(tmp_path / "out")
    counts = short.filter_reads(iter(sam), _bins(tmp_path), out, strict=2, permissive=5)
    assert counts == {"pairs": 1, "recruited": 0, "strict": 0, "permissive": 0}


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
    out = str(tmp_path / "out")
    counts = nano.filter_reads(iter(sam), _bins(tmp_path), out)
    assert counts["recruited"] == 1
    assert _read(os.path.join(out, "bin_001.nanopore.fastq")) == "@r1/1\nACGT\n+\nJJJJ\n"


def test_an_unmapped_long_read_placed_at_a_contig_is_not_recruited(tmp_path):
    """RNAME is set, so the reference-name check passes; only the FLAG says it did not align."""
    sam = ["r1\t%d\tc1\t100\t0\t*\t*\t0\t0\tACGT\tJJJJ" % 0x4]
    out = str(tmp_path / "out")
    counts = nano.filter_reads(iter(sam), _bins(tmp_path), out)
    assert counts["recruited"] == 0
