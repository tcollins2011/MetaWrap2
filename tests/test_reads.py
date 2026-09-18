import gzip

import pytest

from metawrap2.io import reads


def _write_fastq(path, name):
    with open(path, "w") as fh:
        fh.write(f"@{name}\nACGT\n+\nIIII\n")


def test_paired_ok(tmp_path):
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    _write_fastq(r1, "read1/1")
    _write_fastq(r2, "read1/2")
    rs = reads.classify([str(r1), str(r2)])
    assert rs.layout == "paired"
    assert rs.r1 == str(r1) and rs.r2 == str(r2)


def test_paired_mismatch_raises(tmp_path):
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    _write_fastq(r1, "readA/1")
    _write_fastq(r2, "readB/2")
    with pytest.raises(ValueError, match="matching pair"):
        reads.classify([str(r1), str(r2)])


def test_single_requires_one(tmp_path):
    r1 = tmp_path / "s.fastq"
    _write_fastq(r1, "x")
    assert reads.classify([str(r1)], single=True).layout == "single"
    with pytest.raises(ValueError, match="single"):
        reads.classify([str(r1), str(r1)], single=True)


def test_interleaved(tmp_path):
    r1 = tmp_path / "s.fastq"
    _write_fastq(r1, "x")
    assert reads.classify([str(r1)], interleaved=True).layout == "interleaved"


def test_one_file_no_mode_is_helpful_error(tmp_path):
    r1 = tmp_path / "s.fastq"
    _write_fastq(r1, "x")
    with pytest.raises(ValueError, match="--single.*--interleaved|--interleaved"):
        reads.classify([str(r1)])


def test_paired_works_with_gzip(tmp_path):
    r1 = tmp_path / "s_1.fastq.gz"
    r2 = tmp_path / "s_2.fastq.gz"
    with gzip.open(r1, "wt") as fh:
        fh.write("@read1/1\nACGT\n+\nIIII\n")
    with gzip.open(r2, "wt") as fh:
        fh.write("@read1/2\nACGT\n+\nIIII\n")
    assert reads.classify([str(r1), str(r2)]).layout == "paired"
