"""Extract read sequences/qualities from a SAM file into a FASTQ-style file.

Ported from metaWRAP's ``sam_to_fastq.py``. For every alignment record (a line with at
least 11 tab fields) the QNAME, SEQ and QUAL columns are emitted. NOTE: the original wrote
a ``>`` header rather than the FASTQ ``@`` and this behaviour is preserved verbatim, since
MEGAHIT consumes the file downstream and the port must stay byte-compatible.
"""

from __future__ import annotations

import sys
from typing import IO

from ..io.seqio import smart_open


def sam_to_fastq(sam_path: str, out: IO) -> None:
    """Write QNAME/SEQ/QUAL records from *sam_path* to *out*."""
    with smart_open(sam_path, "rt") as fh:
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 11:
                continue
            out.write(">" + cut[0] + "\n" + cut[9] + "\n+\n" + cut[10] + "\n")


def main(argv) -> int:
    sam_to_fastq(argv[0], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
