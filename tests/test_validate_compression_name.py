"""A compressed read file whose name does not say so is reported before any tool runs.

MetaWrap2 reads by content, but the tools go by the suffix: gzip data named ``.fastq`` passed the
validator and FastQC then failed with a traceback.
"""

from __future__ import annotations

import bz2
import gzip

from metawrap2 import validate

RECORD = "@r1\nACGT\n+\nIIII\n"


def test_gzip_data_named_fastq_is_reported_with_the_fix(tmp_path):
    path = tmp_path / "mislabeled_1.fastq"
    path.write_bytes(gzip.compress(RECORD.encode()))
    problem = validate.check_fastq(str(path), "forward reads (-1)")
    assert problem is not None
    assert "mislabeled_1.fastq is gzip-compressed but not named .gz" in problem
    assert "Rename it to mislabeled_1.fastq.gz" in problem
    assert "forward reads (-1)" in problem


def test_bzip2_data_named_fastq_is_reported(tmp_path):
    path = tmp_path / "reads.fq"
    path.write_bytes(bz2.compress(RECORD.encode()))
    problem = validate.check_fastq(str(path))
    assert problem is not None and "bzip2-compressed but not named .bz2" in problem
    assert "reads.fq.bz2" in problem


def test_bzip2_data_named_gz_is_reported(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(bz2.compress(RECORD.encode()))
    problem = validate.check_fastq(str(path))
    assert problem is not None and "Rename it to reads.fastq.bz2" in problem


def test_correctly_named_files_pass(tmp_path):
    plain = tmp_path / "a.fastq"
    plain.write_text(RECORD)
    gz = tmp_path / "a.fastq.gz"
    gz.write_bytes(gzip.compress(RECORD.encode()))
    bz = tmp_path / "a.fq.bz2"
    bz.write_bytes(bz2.compress(RECORD.encode()))
    for path in (plain, gz, bz):
        assert validate.check_fastq(str(path)) is None
