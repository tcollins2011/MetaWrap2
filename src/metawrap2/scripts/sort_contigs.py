"""Sort FASTA contigs by descending sequence length and rewrap at a fixed column width.

Ported from metaWRAP's ``sort_contigs.py``. Reading is gz-transparent via
:mod:`metawrap2.io.seqio`; output is plain FASTA.
"""

from __future__ import annotations

import sys
import textwrap
from typing import IO

from ..io.seqio import iter_fasta

WRAP_WIDTH = 100


def sort_contigs(fasta_path: str, out: IO, width: int = WRAP_WIDTH) -> None:
    """Write contigs of *fasta_path* to *out*, longest first, wrapped at *width* columns."""
    records = list(iter_fasta(fasta_path))
    for header, seq in sorted(records, key=lambda hs: len(hs[1]), reverse=True):
        out.write(">" + header + "\n")
        out.write(textwrap.fill(seq, width, break_on_hyphens=False) + "\n")


def main(argv) -> int:
    sort_contigs(argv[0], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
