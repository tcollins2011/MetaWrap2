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


def test_default_scratch_is_inside_the_output_directory(tmp_path):
    root = scratch.scratch_root(FakeSettings(), str(tmp_path / "study"))
    assert root == str(tmp_path / "study" / ".metawrap2_tmp")


def test_config_overrides_the_default(tmp_path):
    root = scratch.scratch_root(FakeSettings(str(tmp_path / "fast")), str(tmp_path / "study"))
    assert root == str(tmp_path / "fast")


def test_env_var_overrides_the_config(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_SCRATCH", str(tmp_path / "node_local"))
    root = scratch.scratch_root(FakeSettings(str(tmp_path / "fast")), str(tmp_path / "study"))
    assert root == str(tmp_path / "node_local")


def test_scratch_root_expands_a_tilde(monkeypatch):
    monkeypatch.setenv("METAWRAP2_SCRATCH", "~/mw2tmp")
    root = scratch.scratch_root(FakeSettings(), "/tmp/study")
    assert root == os.path.join(os.path.expanduser("~"), "mw2tmp")
    assert "~" not in root


def test_settings_without_a_scratch_dir_attribute_still_work(tmp_path):
    # scratch_root takes `object`, so an older Settings must not raise AttributeError.
    root = scratch.scratch_root(object(), str(tmp_path / "study"))
    assert root.endswith(".metawrap2_tmp")


def test_scratch_dir_creates_a_named_subdirectory(tmp_path):
    path = scratch.scratch_dir(FakeSettings(), str(tmp_path / "study"), "metaspades")
    assert os.path.isdir(path)
    assert path.endswith(os.path.join(".metawrap2_tmp", "metaspades"))


def test_scratch_dir_is_idempotent(tmp_path):
    first = scratch.scratch_dir(FakeSettings(), str(tmp_path / "s"), "megahit")
    assert scratch.scratch_dir(FakeSettings(), str(tmp_path / "s"), "megahit") == first


# --- sizing -------------------------------------------------------------------------------


def test_total_size_sums_files_and_directory_contents(tmp_path):
    (tmp_path / "a.txt").write_text("x" * 100)
    sub = tmp_path / "bins"
    sub.mkdir()
    (sub / "bin.1.fa").write_text("y" * 50)
    (sub / "bin.2.fa").write_text("z" * 25)
    assert scratch.total_size([str(tmp_path / "a.txt"), str(sub)]) == 175


def test_total_size_ignores_missing_and_blank_paths(tmp_path):
    assert scratch.total_size([str(tmp_path / "nope"), "", None or ""]) == 0


def test_estimate_uses_the_module_multiplier(tmp_path):
    big = tmp_path / "reads.fastq"
    big.write_bytes(b"x" * (1024**3))  # 1 GB, sparse enough to be cheap on any modern fs
    assembly = scratch.estimate_bytes("assembly", [str(big)])
    binning = scratch.estimate_bytes("binning", [str(big)])
    assert assembly == int(1024**3 * scratch.SPACE_MULTIPLIERS["assembly"])
    assert assembly > binning  # metaSPAdes is the worst case, by a lot


def test_an_unknown_module_gets_a_conservative_default(tmp_path):
    big = tmp_path / "x"
    big.write_bytes(b"x" * (1024**3))
    assert scratch.estimate_bytes("something_new", [str(big)]) == int(1024**3 * 2.0)


def test_a_tiny_input_still_reserves_the_minimum(tmp_path):
    small = tmp_path / "reads.fastq"
    small.write_text("@r\nACGT\n+\nIIII\n")
    assert scratch.estimate_bytes("assembly", [str(small)]) == scratch.MINIMUM_ESTIMATE_BYTES


def test_human_formats_each_unit():
    assert scratch.human(0) == "0.0 B"
    assert scratch.human(1536) == "1.5 KB"
    assert scratch.human(2 * 1024**3) == "2.0 GB"
    assert scratch.human(None) == "unknown"


# --- free space and the check -------------------------------------------------------------


def test_free_bytes_walks_up_to_an_existing_parent(tmp_path):
    # The scratch directory usually does not exist yet when the check runs.
    deep = str(tmp_path / "not" / "yet" / "created")
    assert scratch.free_bytes(deep) == scratch.free_bytes(str(tmp_path))


def test_check_space_passes_when_there_is_room(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    ok, message = scratch.check_space("kraken2", [str(reads)], [str(tmp_path)])
    assert ok is True
    assert "estimated need" in message and "free" in message


def test_check_space_fails_when_the_estimate_exceeds_free_space(tmp_path, monkeypatch):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(scratch, "free_bytes", lambda path: 1024)
    ok, message = scratch.check_space("assembly", [str(reads)], [str(tmp_path)])
    assert ok is False
    assert "may need about" in message and "1.0 KB free" in message


def test_check_space_reports_one_filesystem_once(tmp_path, monkeypatch):
    """Output and scratch are usually on the same filesystem; saying so twice is noise."""
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "study"
    out.mkdir()
    ok, message = scratch.check_space(
        "binning", [str(reads)], [str(out), str(out / ".metawrap2_tmp")]
    )
    assert ok is True
    assert message.count(" free") == 1


def test_check_space_tolerates_an_unreadable_target(tmp_path, monkeypatch):
    monkeypatch.setattr(scratch, "free_bytes", lambda path: None)
    ok, message = scratch.check_space("binning", [], [str(tmp_path)])
    assert ok is True
    assert "no target checked" in message


def test_check_space_ignores_blank_targets(tmp_path):
    ok, _ = scratch.check_space("binning", [], ["", None or ""])
    assert ok is True
