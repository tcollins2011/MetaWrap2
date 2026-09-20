"""Drop short contigs from a length-sorted SPAdes FASTA.

Ported from metaWRAP's ``rm_short_contigs.py``, with one correctness fix: the length now
comes from the actual sequence rather than from parsing field 3 of a ``_``-split
SPAdes-style header (``NODE_x_length_<len>_cov_...``).

The old parse raised ValueError on any assembler that does not use SPAdes naming, and it
also stopped at the *first* short contig - correct only because metaSPAdes output happens to
be sorted longest-first. Measuring the sequence and filtering every record gives identical
results on sorted input and correct results on unsorted input. Reading is gz-transparent via
:mod:`metawrap2.io.seqio`.
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
