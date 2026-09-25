"""Tests for metawrap2.validate: catching empty and truncated files before a tool sees them.

The point of these checks is that a *present but unusable* file is worse than a missing one,
because it propagates silently through the pipeline. So the cases that matter here are the
plausible-looking ones: the library cut off by a full disk, the gzip whose download stopped
halfway, the FASTA ending on a header with no sequence.
"""

from __future__ import annotations

import gzip
import os

import pytest

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


# --- require_nonempty_dir -----------------------------------------------------------------


def test_require_nonempty_dir_checks_extensions(tmp_path):
    # The real failure this catches: a bins directory holding CheckM leftovers but no FASTAs.
    write(str(tmp_path / "checkm.log"), "...")
    problem = validate.require_nonempty_dir(str(tmp_path), "bin folder", extensions=(".fa",))
    assert problem and "contains no .fa files" in problem
    write(str(tmp_path / "bin.1.fa"), ">c1\nACGT\n")
    assert validate.require_nonempty_dir(str(tmp_path), "bin folder", extensions=(".fa",)) is None


# --- check_fastq --------------------------------------------------------------------------


def test_check_fastq_rejects_a_truncated_gzip(tmp_path):
    """An interrupted download: valid gzip header, stream stops mid-deflate.

    This is the case that used to escape as a CRC error from bwa or metaSPAdes, hours later.
    """
    path = write_gz(str(tmp_path / "part.fastq.gz"), good_fastq(500))
    full = os.path.getsize(path)
    with open(path, "r+b") as fh:
        fh.truncate(full // 2)
    problem = validate.check_fastq(path)
    assert problem and ("truncated or corrupt" in problem or "could not be read" in problem)


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
    """A file cut on a line boundary to a non-multiple-of-four line count is still caught."""
    text = (
        good_fastq(500) + "@read500\nACGTACGTAC\n+\n"
    )  # 3 extra lines -> total not a multiple of 4
    path = write_gz(str(tmp_path / "cut.fastq.gz"), text)
    problem = validate.check_fastq(path)
    assert problem and "not a multiple of 4" in problem


# --- check_fasta --------------------------------------------------------------------------


def test_check_fasta_rejects_a_trailing_header_with_no_sequence(tmp_path):
    problem = validate.check_fasta(write(str(tmp_path / "cut.fa"), ">c1\nACGT\n>c2\n"))
    assert problem and "ends with a header line" in problem


# --- require_inputs -----------------------------------------------------------------------


def test_require_inputs_reports_every_problem_at_once(tmp_path):
    """The whole reason this helper exists: five samples with three bad files is one run.

    Also exercises the accept path (the ok file) and the missing/empty rejections via check_fastq.
    """
    ok = write(str(tmp_path / "ok_1.fastq"), good_fastq())
    empty = write(str(tmp_path / "empty_2.fastq"), "")
    missing = str(tmp_path / "gone_1.fastq")
    with pytest.raises(validate.InputProblem) as excinfo:
        validate.require_inputs(
            [
                (validate.check_fastq, ok, "sample 1 forward reads"),
                (validate.check_fastq, empty, "sample 1 reverse reads"),
                (validate.check_fastq, missing, "sample 2 forward reads"),
            ],
            context="binning",
        )
    message = str(excinfo.value)
    assert "Cannot start binning: 2 input problem(s) found." in message
    assert message.count("\n  - ") == 2
    assert "empty_2.fastq" in message and "gone_1.fastq" in message
    assert "ok_1.fastq" not in message
