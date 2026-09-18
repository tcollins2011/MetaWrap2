"""Filter bmtagger-flagged host reads OUT of a FASTQ file (keep the non-host reads).

Ported from metaWRAP's ``skip_human_reads.py``. Given the bmtagger read-id list and the
original FASTQ, every read whose id is *not* in the list is written out. Reading is
gz-transparent via :mod:`metawrap2.io.seqio`.

Run standalone (``python -m metawrap2.scripts.skip_human_reads list.txt reads.fastq``) or
import :func:`filter_reads` in-process.
"""

from __future__ import annotations

import sys
from typing import IO

from ..io.seqio import smart_open


def load_ids(list_path: str) -> set:
    """Load the set of read ids flagged by bmtagger."""
    ids = set()
    with smart_open(list_path, "rt") as fh:
        for line in fh:
            ids.add(line.strip())
    return ids


def filter_reads(list_path: str, fastq_path: str, out: IO, keep_host: bool = False) -> None:
    """Write reads from *fastq_path* to file-object *out*.

    With ``keep_host=False`` (default) host reads are dropped; with ``keep_host=True`` only
    the host reads are kept (this is what ``select_human_reads`` does). Read identity is the
    first whitespace/``/``-delimited token of the header, matching the original.
    """
    host = load_ids(list_path)
    skip = keep_host  # before the first header we keep nothing extra
    with smart_open(fastq_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                token = line[1:].split("/")[0].split()
                in_host = bool(token) and token[0] in host
                skip = (not in_host) if keep_host else in_host
            if not skip:
                out.write(line.rstrip("\n") + "\n")


def main(argv) -> int:
    filter_reads(argv[0], argv[1], sys.stdout, keep_host=False)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
