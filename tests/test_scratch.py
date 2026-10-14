"""Tests for metawrap2.scratch: one scratch-directory policy, and the preflight space check."""

from __future__ import annotations

import os

import pytest

from metawrap2 import scratch


class FakeSettings:
    def __init__(self, scratch_dir: str = ""):
        self.scratch_dir = scratch_dir


@pytest.fixture(autouse=True)
def _no_ambient_scratch(monkeypatch):
    # The env var must not leak in from the developer's own shell.
    monkeypatch.delenv("METAWRAP2_SCRATCH", raising=False)


# --- free space and the check -------------------------------------------------------------


def test_check_space_fails_when_the_estimate_exceeds_free_space(tmp_path, monkeypatch):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(scratch, "free_bytes", lambda path: 1024)
    ok, message = scratch.check_space("assembly", [str(reads)], [str(tmp_path)])
    assert ok is False
    assert "may need about" in message and "1.0 KB free" in message

    # The scratch directory itself is a named subdirectory, created on demand.
    path = scratch.scratch_dir(FakeSettings(), str(tmp_path / "study"), "metaspades")
    assert os.path.isdir(path)
    assert path.endswith(os.path.join(".metawrap2_tmp", "metaspades"))
