#!/usr/bin/env python
"""Shorten long, underscore-heavy contig names so PROKKA will accept them.

Streams a FASTA file. Sequence lines pass through unchanged. Once a header longer than 20
characters with more than five underscore-delimited fields is seen, that header and every later one
are truncated to their first four fields, as in the original.

Port of metaWRAP's ``shorten_contig_names.py``. ``.gz`` input is supported.
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
