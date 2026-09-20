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


# --- resolution order ---------------------------------------------------------------------


def test_env_var_overrides_the_config(tmp_path, monkeypatch):
    # env var wins over the configured scratch_dir, which in turn wins over the default
    monkeypatch.setenv("METAWRAP2_SCRATCH", str(tmp_path / "node_local"))
    root = scratch.scratch_root(FakeSettings(str(tmp_path / "fast")), str(tmp_path / "study"))
    assert root == str(tmp_path / "node_local")


def test_scratch_dir_creates_a_named_subdirectory(tmp_path):
    path = scratch.scratch_dir(FakeSettings(), str(tmp_path / "study"), "metaspades")
    assert os.path.isdir(path)
    assert path.endswith(os.path.join(".metawrap2_tmp", "metaspades"))


# --- sizing -------------------------------------------------------------------------------


def test_estimate_uses_the_module_multiplier(tmp_path):
    big = tmp_path / "reads.fastq"
    big.write_bytes(b"x" * (1024**3))  # 1 GB, sparse enough to be cheap on any modern fs
    assembly = scratch.estimate_bytes("assembly", [str(big)])
    binning = scratch.estimate_bytes("binning", [str(big)])
    assert assembly == int(1024**3 * scratch.SPACE_MULTIPLIERS["assembly"])
    assert assembly > binning  # metaSPAdes is the worst case, by a lot


# --- free space and the check -------------------------------------------------------------


def test_check_space_fails_when_the_estimate_exceeds_free_space(tmp_path, monkeypatch):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(scratch, "free_bytes", lambda path: 1024)
    ok, message = scratch.check_space("assembly", [str(reads)], [str(tmp_path)])
    assert ok is False
    assert "may need about" in message and "1.0 KB free" in message
