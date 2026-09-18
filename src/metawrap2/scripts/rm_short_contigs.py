"""Drop short contigs from a length-sorted SPAdes FASTA.

Ported from metaWRAP's ``rm_short_contigs.py``. The input is assumed to be sorted longest
first (as metaSPAdes ``scaffolds.fasta`` is) and to use SPAdes-style headers
(``NODE_x_length_<len>_cov_...``): the contig length is field index 3 of the ``_``-split
header. Iteration stops at the first contig shorter than the cutoff, exactly like the
original. Reading is gz-transparent via :mod:`metawrap2.io.seqio`.
"""

from __future__ import annotations

import sys
from typing import IO

from ..io.seqio import iter_fasta


def remove_short_contigs(min_len: int, fasta_path: str, out: IO) -> None:
    """Write contigs of *fasta_path* to *out* until one shorter than *min_len* is reached."""
    for header, seq in iter_fasta(fasta_path):
        if int(header.split("_")[3]) < min_len:
            break
        out.write(">" + header + "\n" + seq + "\n")


def main(argv) -> int:
    remove_short_contigs(int(argv[0]), argv[1], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
