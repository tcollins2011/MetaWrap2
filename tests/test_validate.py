"""Tests for metawrap2.validate: catching empty and truncated files before a tool sees them.

A present but unusable file can pass through the pipeline unnoticed, so the cases tested are
files that look valid: a library cut off by a full disk, a gzip whose download stopped halfway, a
FASTA ending on a header with no sequence.
"""

from __future__ import annotations

import gzip
import os

from metawrap2 import validate

FASTQ_RECORD = "@read%d\nACGTACGTAC\n+\nIIIIIIIIII\n"


def write(path: str, text: str) -> str:
    with open(path, "w") as fh:
        fh.write(text)
    return path


def write_gz(path: str, text: str) -> str:
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def good_fastq(n: int = 5) -> str:
    return "".join(FASTQ_RECORD % i for i in range(n))


# --- check_fastq --------------------------------------------------------------------------


def test_check_fastq_accepts_a_large_valid_gzip(tmp_path):
    """Regression: a complete, valid gzipped FASTQ larger than the read window must pass.

    The old check counted lines in only the last TAIL_BYTES of the file and required that
    tail-window count to be a multiple of four - which it is not for an arbitrary byte offset,
    so large valid libraries were wrongly rejected as "truncated". The count must be over the
    whole file.
    """
    records = validate.TAIL_BYTES // len(good_fastq(1)) + 5000  # comfortably exceeds the window
    path = write_gz(str(tmp_path / "big.fastq.gz"), good_fastq(records))
    assert validate.check_fastq(path) is None


def test_check_fastq_rejects_a_line_truncation(tmp_path):
    """A file cut on a line boundary to a non-multiple-of-four line count is still caught.

    Also the interrupted download: valid gzip header, stream stops mid-deflate. This is the case
    that used to escape as a CRC error from bwa or metaSPAdes, hours later.
    """
    text = (
        good_fastq(500) + "@read500\nACGTACGTAC\n+\n"
    )  # 3 extra lines -> total not a multiple of 4
    path = write_gz(str(tmp_path / "cut.fastq.gz"), text)
    problem = validate.check_fastq(path)
    assert problem and "not a multiple of 4" in problem

    path = write_gz(str(tmp_path / "part.fastq.gz"), good_fastq(500))
    full = os.path.getsize(path)
    with open(path, "r+b") as fh:
        fh.truncate(full // 2)
    problem = validate.check_fastq(path)
    assert problem and ("truncated or corrupt" in problem or "could not be read" in problem)
