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






def test_interleaved(tmp_path):
    r1 = tmp_path / "s.fastq"
    _write_fastq(r1, "x")
    assert reads.classify([str(r1)], interleaved=True).layout == "interleaved"




