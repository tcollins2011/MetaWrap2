"""Rename MEGAHIT contigs into SPAdes-style headers, drop short ones, sort by length.

Ported from metaWRAP's ``fix_megahit_contig_naming.py``. MEGAHIT headers look like
``k141_0 flag=1 multi=2.0000 len=500``; they are rewritten as
``k141_0_length_<len>_cov_<multi>`` so downstream tools that assume SPAdes naming keep
working. Contigs shorter than the cutoff are removed and the rest are sorted longest first,
wrapped at 100 columns. Reading is gz-transparent via :mod:`metawrap2.io.seqio`.
"""

from __future__ import annotations

import sys
import textwrap
from typing import IO

from ..io.seqio import iter_fasta

WRAP_WIDTH = 100


def fix_naming(fasta_path: str, min_len: int, out: IO, width: int = WRAP_WIDTH) -> None:
    """Rename/filter/sort MEGAHIT contigs from *fasta_path* into *out*."""
    renamed = []
    for header, seq in iter_fasta(fasta_path):
        parts = header.split(" ")
        # Read len=/multi= by name rather than by position: MEGAHIT's header layout has
        # changed between releases, and positional indexing turns that into an IndexError
        # (or, worse, silently reads the wrong field) instead of just working.
        fields = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
        length = fields.get("len") or str(len(seq))
        cov = fields.get("multi", "0")
        if int(length) < min_len:
            continue
        renamed.append((parts[0] + "_length_" + length + "_cov_" + cov, seq))
    for name, seq in sorted(renamed, key=lambda ns: len(ns[1]), reverse=True):
        out.write(">" + name + "\n")
        out.write(textwrap.fill(seq, width, break_on_hyphens=False) + "\n")


def main(argv) -> int:
    fix_naming(argv[0], int(argv[1]), sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
