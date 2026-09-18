#!/usr/bin/env python
"""Shorten long, underscore-heavy contig names so PROKKA will accept them.

Streams a FASTA file: sequence lines pass through unchanged, and header lines are truncated
to their first four underscore-delimited fields once a header long enough to need it is seen
(from then on every header is shortened, matching the original one-shot latching behavior).

Ported faithfully from metaWRAP's ``shorten_contig_names.py``. ``.gz`` input is supported via
the shared sequence reader.
"""

from __future__ import annotations

import sys
from typing import Iterator, List

from ..io.seqio import smart_open


def shorten_lines(path: str) -> Iterator[str]:
    """Yield output lines (no trailing newline) with contig headers shortened as needed."""
    shorten = False
    with smart_open(path, "rt") as fh:
        for line in fh:
            if line[:1] != ">":
                yield line.rstrip()
            else:
                if shorten:
                    yield "_".join(line.rstrip().split("_")[:4])
                elif len(line) > 20 and len(line.split("_")) > 5:
                    yield "_".join(line.rstrip().split("_")[:4])
                    shorten = True
                else:
                    yield line.rstrip()


def main(argv: List[str]) -> int:
    for line in shorten_lines(argv[0]):
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
