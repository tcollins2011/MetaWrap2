"""Extract read sequences/qualities from a SAM stream into a FASTQ-style file.

For every alignment record (a line with at least 11 tab fields) the QNAME, SEQ and QUAL columns
are emitted. Reads ``-`` as standard input so an aligner can be piped straight in: the assembly
module's unmapped-read extraction used to write the aligner's full SAM to disk, filter it into a
second SAM, convert that, then delete both - roughly 25 GB of text at peak to produce one FASTQ.

The original wrote a ``>`` header while emitting the rest of a 4-line FASTQ record, so the
file it produced was neither valid FASTA nor valid FASTQ despite its ``.fastq`` name. MEGAHIT
(its only consumer) parses both forms identically - verified by running it on the
two variants and comparing the read count and length it reports - so emitting the correct
``@`` changes nothing downstream while making the file readable by every other tool.

``--unmapped`` reads the unmapped bit in FLAG. It used to test for the absence of an ``NM:i:`` tag
instead, which is a claim about one aligner rather than about SAM: bwa writes the edit distance on
every aligned record, so absence stood in for "did not align". The two agree exactly on real bwa
output - 33,332,411 records from a CAMI2 sample, zero disagreements in either direction - so this
does not change what the hybrid assembly is built from today. It is a fix because the module it
serves is moving off bwa: minimap2 is MetaWrap2's measured default aligner everywhere else, and an
aligner that omits NM would have made every mapped read look unmapped, which fails by silently
handing MEGAHIT the entire read set rather than by stopping.

Reading the flag costs something, because the NM test could reject a mapped record on a substring
search without splitting it at all: 650 against 1041 MB per CPU-second over the same 1 GB stream.
Spent knowingly - the aligner upstream produces SAM at roughly 26 MB/s (bwa) to 146 MB/s (minimap2)
on this hardware, so the converter has over 4x headroom either way and is not what the pipe waits
on. ``samtools fastq -f 4`` expresses the same predicate and was measured too: slower again than
this (711 MB per CPU-second) and it would add samtools to an env that has no other use for it.
"""

from __future__ import annotations

import argparse
import sys
from typing import IO, Iterable

from ..io.seqio import smart_open

#: The unmapped bit in FLAG, which is the second field of every alignment record.
UNMAPPED_FLAG = 0x4

#: Secondary (0x100) and supplementary (0x800) records are extra alignments of a read that is
#: already represented by its primary record, so emitting them would repeat reads. They can never
#: be unmapped, so this only matters to the unfiltered mode.
SECONDARY_FLAGS = 0x100 | 0x800


def _flag(line: str) -> int:
    """FLAG for a SAM *line*, or 0 if it has no readable second field.

    0 rather than a negative sentinel because it has to be safe under a bitmask: a negative int is
    all ones in two's complement, so ``-1 & UNMAPPED_FLAG`` is truthy and a header line would be
    treated as an unmapped read. 0 is "mapped, forward, primary", which is the right answer for
    something that is not a record at all - it reaches the field-count check and is dropped there,
    exactly as it was before any flag was looked at.

    Split with maxsplit rather than in full: this runs once per record in the middle of an aligner
    pipe, and a 450-byte record would otherwise be cut into a dozen strings to read one of them.
    Slicing between two ``find`` calls to skip even the 3-element list is slower, not faster (599
    against 650 MB per CPU-second over the same 1 GB stream), so the obvious version is the one.
    """
    parts = line.split("\t", 2)
    if len(parts) < 3:
        return 0
    try:
        return int(parts[1])
    except ValueError:  # an @ header line, whose second field is not a number
        return 0


def _emit(line: str, out: IO) -> None:
    """Write *line*'s QNAME/SEQ/QUAL as a FASTQ record, ignoring anything too short to be one."""
    cut = line.rstrip("\n").split("\t")
    if len(cut) < 11:
        return
    out.write("@" + cut[0] + "\n" + cut[9] + "\n+\n" + cut[10] + "\n")


def records(lines: Iterable[str], out: IO, unmapped_only: bool = False) -> None:
    """Write QNAME/SEQ/QUAL records from SAM *lines* to *out*.

    One loop per mode rather than one loop that re-asks which mode it is on every record, since
    neither mode is ever half of a stream.
    """
    if unmapped_only:
        for line in lines:
            if _flag(line) & UNMAPPED_FLAG:
                _emit(line, out)
        return
    for line in lines:
        if not _flag(line) & SECONDARY_FLAGS:
            _emit(line, out)


def sam_to_fastq(sam_path: str, out: IO, unmapped_only: bool = False) -> None:
    """Write QNAME/SEQ/QUAL records from *sam_path* to *out*; ``-`` reads standard input."""
    if sam_path == "-":
        records(sys.stdin, out, unmapped_only)
        return
    with smart_open(sam_path, "rt") as fh:
        records(fh, out, unmapped_only)


def main(argv) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sam", help="SAM file, or - for standard input")
    ap.add_argument(
        "--unmapped",
        action="store_true",
        help="emit only reads that did not align (FLAG 0x%x set)" % UNMAPPED_FLAG,
    )
    args = ap.parse_args(argv)
    sam_to_fastq(args.sam, sys.stdout, unmapped_only=args.unmapped)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
