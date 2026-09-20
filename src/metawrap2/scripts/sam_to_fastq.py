"""Extract read sequences/qualities from a SAM file into a FASTQ-style file.

Ported from metaWRAP's ``sam_to_fastq.py``. For every alignment record (a line with at
least 11 tab fields) the QNAME, SEQ and QUAL columns are emitted.

The original wrote a ``>`` header while emitting the rest of a 4-line FASTQ record, so the
file it produced was neither valid FASTA nor valid FASTQ despite its ``.fastq`` name. MEGAHIT
(its only consumer) parses both forms identically - verified by running it on the
two variants and comparing the read count and length it reports - so emitting the correct
``@`` changes nothing downstream while making the file readable by every other tool.
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
            out.write("@" + cut[0] + "\n" + cut[9] + "\n+\n" + cut[10] + "\n")


def main(argv) -> int:
    sam_to_fastq(argv[0], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
