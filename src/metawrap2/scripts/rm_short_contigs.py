"""Drop short contigs from a length-sorted SPAdes FASTA.

Ported from metaWRAP's ``rm_short_contigs.py``. The length is taken from the sequence instead
of field 3 of a ``_``-split SPAdes-style header (``NODE_x_length_<len>_cov_...``).

The old parse raised ValueError on non-SPAdes headers and stopped at the first short contig,
which relied on the input being sorted longest-first. Every record is now filtered, so unsorted
input also works. Reading is gz-transparent via :mod:`metawrap2.io.seqio`.
"""

from __future__ import annotations

import sys
from typing import IO

from ..io.seqio import iter_fasta


def remove_short_contigs(min_len: int, fasta_path: str, out: IO) -> None:
    """Write every contig of *fasta_path* at least *min_len* bases long to *out*."""
    for header, seq in iter_fasta(fasta_path):
        if len(seq) < min_len:
            continue
        out.write(">" + header + "\n" + seq + "\n")


def main(argv) -> int:
    remove_short_contigs(int(argv[0]), argv[1], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
