"""Tests for metawrap2.validate: catching empty and truncated files before a tool sees them.

The point of these checks is that a *present but unusable* file is worse than a missing one,
because it propagates silently through the pipeline. So the cases that matter here are the
plausible-looking ones: the zero-byte FASTQ, the library cut off by a full disk, the gzip whose
download stopped halfway, the FASTA ending on a header with no sequence.
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


# --- require_file -------------------------------------------------------------------------


def test_require_file_accepts_a_normal_file(tmp_path):
    path = write(str(tmp_path / "reads.fastq"), good_fastq())
    assert validate.require_file(path) is None


def test_require_file_rejects_missing(tmp_path):
    problem = validate.require_file(str(tmp_path / "nope.fastq"))
    assert problem and "does not exist" in problem


def test_require_file_rejects_zero_bytes(tmp_path):
    path = write(str(tmp_path / "empty.fastq"), "")
    problem = validate.require_file(path)
    assert problem and "0 bytes" in problem


def test_require_file_rejects_a_directory(tmp_path):
    problem = validate.require_file(str(tmp_path))
    assert problem and "is a directory" in problem


def test_require_file_rejects_the_empty_string(tmp_path):
    # An unset argument arrives as "" rather than None often enough to be worth pinning.
    problem = validate.require_file("", "assembly")
    assert problem == "no assembly was given"


# --- require_nonempty_dir -----------------------------------------------------------------


def test_require_nonempty_dir_accepts_a_dir_with_files(tmp_path):
    write(str(tmp_path / "bin.1.fa"), ">c1\nACGT\n")
    assert validate.require_nonempty_dir(str(tmp_path)) is None


def test_require_nonempty_dir_rejects_empty(tmp_path):
    empty = tmp_path / "bins"
    empty.mkdir()
    problem = validate.require_nonempty_dir(str(empty))
    assert problem and "is empty" in problem


def test_require_nonempty_dir_rejects_a_file(tmp_path):
    path = write(str(tmp_path / "bins"), "x")
    problem = validate.require_nonempty_dir(path)
    assert problem and "is not a directory" in problem


def test_require_nonempty_dir_checks_extensions(tmp_path):
    # The real failure this catches: a bins directory holding CheckM leftovers but no FASTAs.
    write(str(tmp_path / "checkm.log"), "...")
    problem = validate.require_nonempty_dir(str(tmp_path), "bin folder", extensions=(".fa",))
    assert problem and "contains no .fa files" in problem
    write(str(tmp_path / "bin.1.fa"), ">c1\nACGT\n")
    assert validate.require_nonempty_dir(str(tmp_path), "bin folder", extensions=(".fa",)) is None


# --- check_fastq --------------------------------------------------------------------------


def test_check_fastq_accepts_plain_and_gzipped(tmp_path):
    assert validate.check_fastq(write(str(tmp_path / "a.fastq"), good_fastq())) is None
    assert validate.check_fastq(write_gz(str(tmp_path / "a.fastq.gz"), good_fastq())) is None


def test_check_fastq_rejects_fasta_content(tmp_path):
    path = write(str(tmp_path / "wrong.fastq"), ">contig1\nACGT\n")
    problem = validate.check_fastq(path)
    assert problem and "does not look like FASTQ" in problem


def test_check_fastq_rejects_a_record_cut_in_half(tmp_path):
    # What a full disk leaves behind: a whole number of lines, but not a multiple of four.
    path = write(str(tmp_path / "cut.fastq"), good_fastq(3) + "@read3\nACGTACGTAC\n")
    problem = validate.check_fastq(path)
    assert problem and "ends mid-record" in problem


def test_check_fastq_rejects_a_missing_final_newline(tmp_path):
    path = write(str(tmp_path / "cut.fastq"), good_fastq(2).rstrip("\n"))
    problem = validate.check_fastq(path)
    assert problem and "does not end with a newline" in problem


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


def test_check_fastq_rejects_empty(tmp_path):
    problem = validate.check_fastq(write(str(tmp_path / "e.fastq"), ""))
    assert problem and "0 bytes" in problem


# --- check_fasta --------------------------------------------------------------------------


def test_check_fasta_accepts_plain_and_gzipped(tmp_path):
    fasta = ">c1\nACGTACGT\n>c2\nTTTTGGGG\n"
    assert validate.check_fasta(write(str(tmp_path / "a.fa"), fasta)) is None
    assert validate.check_fasta(write_gz(str(tmp_path / "a.fa.gz"), fasta)) is None


def test_check_fasta_rejects_fastq_content(tmp_path):
    problem = validate.check_fasta(write(str(tmp_path / "wrong.fa"), good_fastq(1)))
    assert problem and "does not look like FASTA" in problem


def test_check_fasta_rejects_a_trailing_header_with_no_sequence(tmp_path):
    problem = validate.check_fasta(write(str(tmp_path / "cut.fa"), ">c1\nACGT\n>c2\n"))
    assert problem and "ends with a header line" in problem


def test_check_fasta_rejects_a_header_only_file(tmp_path):
    # Single line, no newline: the degenerate form of the same truncation.
    problem = validate.check_fasta(write(str(tmp_path / "cut.fa"), ">c1"))
    assert problem and "ends with a header line" in problem


# --- require_inputs -----------------------------------------------------------------------


def test_require_inputs_passes_when_everything_is_fine(tmp_path):
    r1 = write(str(tmp_path / "s_1.fastq"), good_fastq())
    r2 = write(str(tmp_path / "s_2.fastq"), good_fastq())
    validate.require_inputs(
        [(validate.check_fastq, r1, "forward reads"), (validate.check_fastq, r2, "reverse reads")]
    )


def test_require_inputs_reports_every_problem_at_once(tmp_path):
    """The whole reason this helper exists: five samples with three bad files is one run."""
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


def test_require_inputs_treats_none_as_a_missing_argument(tmp_path):
    with pytest.raises(validate.InputProblem, match="no assembly was given"):
        validate.require_inputs([(validate.check_fasta, None, "assembly")])
