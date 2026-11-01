"""An interleaved file passed with ``--single-end`` is run as asked, with a warning.

Run as single-end, interleaved data doubles the read count and loses the pairing.
"""

from __future__ import annotations

import argparse
import gzip

from metawrap2.io import reads
from metawrap2.modules import _common


def _fastq(names):
    return "".join("@%s\nACGT\n+\nIIII\n" % name for name in names)


def test_mates_with_slash_suffixes_look_interleaved(tmp_path):
    path = tmp_path / "x.fastq"
    path.write_text(_fastq(["a/1", "a/2", "b/1", "b/2", "c/1"]))
    assert reads.looks_interleaved(str(path))


def test_mates_with_identical_names_look_interleaved(tmp_path):
    path = tmp_path / "x.fastq.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(_fastq(["a 1:N:0:1", "a 2:N:0:1", "b 1:N:0:1", "b 2:N:0:1"]))
    assert reads.looks_interleaved(str(path))


def test_single_end_reads_do_not(tmp_path):
    path = tmp_path / "x.fastq"
    path.write_text(_fastq(["a", "b", "c", "d"]))
    assert not reads.looks_interleaved(str(path))
    one = tmp_path / "one.fastq"
    one.write_text(_fastq(["a"]))
    assert not reads.looks_interleaved(str(one))
    assert not reads.looks_interleaved(str(tmp_path / "missing.fastq"))


def test_read_layout_warns_and_keeps_single_end(tmp_path, monkeypatch):
    path = tmp_path / "x.fastq"
    path.write_text(_fastq(["a/1", "a/2", "b/1", "b/2"]))
    warned = []
    monkeypatch.setattr(_common, "warning", warned.append)
    args = argparse.Namespace(reads_1=str(path), reads_2=None, single_end=True, interleaved=False)
    layout = _common.read_layout(args)
    assert layout.layout == "single"
    (message,) = warned
    assert "looks interleaved" in message and "--interleaved" in message


def test_real_single_end_data_gets_no_warning(tmp_path, monkeypatch):
    path = tmp_path / "x.fastq"
    path.write_text(_fastq(["a", "b", "c", "d"]))
    warned = []
    monkeypatch.setattr(_common, "warning", warned.append)
    args = argparse.Namespace(reads_1=str(path), reads_2=None, single_end=True, interleaved=False)
    assert _common.read_layout(args).layout == "single"
    assert warned == []
