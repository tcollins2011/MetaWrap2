"""Extract read sequences and qualities from a SAM stream into a FASTQ file.

For every alignment record (a line with at least 11 tab fields) the QNAME, SEQ and QUAL columns
are written. ``-`` reads standard input, so an aligner can be piped in directly. The assembly
module's unmapped-read extraction previously wrote the full SAM to disk, filtered it into a second
SAM and converted that, about 25 GB of text at peak for one FASTQ.

The original wrote a ``>`` header on a 4-line FASTQ record, which is neither valid FASTA nor valid
FASTQ. MEGAHIT, the only consumer, reports the same read count and length for both forms, so
writing ``@`` changes nothing downstream and makes the file readable by other tools.

``--unmapped`` tests the unmapped bit in FLAG. The original tested for a missing ``NM:i:`` tag,
which only works for aligners that, like bwa, write NM on every aligned record. On 33,332,411 bwa
records from a CAMI2 sample the two tests agree on every record. An aligner that omits NM (the
module is moving to minimap2) would make every mapped read look unmapped and pass the whole read
set to MEGAHIT without an error.

Reading the flag is slower than the NM substring test (650 against 1041 MB per CPU-second over the
same 1 GB stream). The upstream aligner produces SAM at about 26 MB/s (bwa) to 146 MB/s (minimap2),
so the converter still has over 4x headroom. ``samtools fastq -f 4`` was also measured: 711 MB per
CPU-second, and it would add samtools to an env that does not otherwise need it.
"""

from __future__ import annotations

import argparse
import sys
from typing import IO, Iterable

from ..io.seqio import smart_open

#: The unmapped bit in FLAG, which is the second field of every alignment record.
UNMAPPED_FLAG = 0x4

#: Secondary (0x100) and supplementary (0x800) records repeat a read already given by its primary
#: record. They are never unmapped, so this only applies to the unfiltered mode.
SECONDARY_FLAGS = 0x100 | 0x800


def _flag(line: str) -> int:
    """FLAG for a SAM *line*, or 0 if it has no readable second field.

    The fallback is 0 because a negative sentinel such as -1 has every bit set, so a header line
    would test as unmapped. A line with flag 0 is later dropped by the field-count check.

    ``maxsplit`` avoids splitting the whole record on every line of the stream. Slicing between
    two ``find`` calls was slower (599 against 650 MB per CPU-second over the same 1 GB stream).
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

    Each mode has its own loop so the mode is not tested on every record.
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
