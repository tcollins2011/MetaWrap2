"""Select the bmtagger-flagged host reads FROM a FASTQ file (inverse of skip_human_reads).

Ported from metaWRAP's ``select_human_reads.py``. Only the reads whose id appears in the
bmtagger list are written out, so the removed host reads can be kept "for science".
"""

from __future__ import annotations

import sys
from typing import IO

from .skip_human_reads import filter_reads


def select_reads(list_path: str, fastq_path: str, out: IO) -> None:
    """Write only the host reads from *fastq_path* to file-object *out*."""
    filter_reads(list_path, fastq_path, out, keep_host=True)


def main(argv) -> int:
    select_reads(argv[0], argv[1], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
