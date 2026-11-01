"""``--dry-run`` records the plan and nothing else.

A dry run of read_qc wrote four ``.done`` markers and the decompressed reads. The real run into the
same ``-o`` then refused to start, and ``--resume`` skipped every step on the dry run's markers and
reported success with no output.
"""

from __future__ import annotations

import gzip
import json
import os
import shlex

import pytest

from metawrap2 import command
from metawrap2.checkpoint import Checkpoint, make_checkpoint
from metawrap2.modules import _common, read_qc


@pytest.fixture(autouse=True)
def _reset_runner(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    command.runner.configure(resume=False, force=False, dry_run=False)
    yield
    command.runner.configure(resume=False, force=False, dry_run=False)


def _marker(out, step):
    return os.path.join(str(out), ".metawrap2", "steps", step + ".done")


# ── step markers ─────────────────────────────────────────────────────────────────────────────


def test_a_dry_run_writes_no_step_marker(tmp_path):
    command.runner.configure(dry_run=True)
    Checkpoint(str(tmp_path), resume=False).done("trimming")
    assert not os.path.exists(_marker(tmp_path, "trimming"))


def test_a_real_marker_records_that_it_was_not_a_dry_run(tmp_path):
    Checkpoint(str(tmp_path), resume=False).done("trimming")
    with open(_marker(tmp_path, "trimming")) as fh:
        assert json.load(fh) == {"status": "done", "dry_run": False}
    command.runner.configure(resume=True)
    assert make_checkpoint(str(tmp_path)).todo("trimming") is False


@pytest.mark.parametrize(
    "content, honoured",
    [
        ("done\n", True),  # written before 2.1.5
        ('{"status": "done", "dry_run": true}\n', False),
        ('{"status": "done"}\n', False),
        ("", False),
        ("not json\n", False),
    ],
)
def test_resume_honours_only_markers_of_real_runs(tmp_path, content, honoured):
    os.makedirs(os.path.dirname(_marker(tmp_path, "trimming")))
    with open(_marker(tmp_path, "trimming"), "w") as fh:
        fh.write(content)
    ckpt = Checkpoint(str(tmp_path), resume=True)
    assert ckpt.todo("trimming") is (not honoured)


def test_a_dry_run_does_not_clear_markers(tmp_path):
    Checkpoint(str(tmp_path), resume=False).done("align_reads")
    command.runner.configure(dry_run=True)
    Checkpoint(str(tmp_path), resume=False).clear("align_reads")
    assert os.path.exists(_marker(tmp_path, "align_reads"))


# ── read_qc end to end ───────────────────────────────────────────────────────────────────────


def _gz_reads(tmp_path):
    paths = []
    for mate in ("1", "2"):
        path = tmp_path / ("t_%s.fastq.gz" % mate)
        with gzip.open(path, "wt") as fh:
            fh.write("@r1/%s\nACGT\n+\nIIII\n" % mate)
        paths.append(str(path))
    return paths


def _fake_tools(monkeypatch, calls):
    def fake_run(cmd, env=None, tool=None, **kwargs):
        calls.append(cmd)
        if command.runner.dry_run:
            return 0
        parts = shlex.split(cmd)
        if tool == "fastp":
            for flag in ("-o", "-O"):
                with open(parts[parts.index(flag) + 1], "w") as fh:
                    fh.write("@r\nACGT\n+\nIIII\n")
        return 0

    monkeypatch.setattr(read_qc, "run", fake_run)
    monkeypatch.setattr(read_qc, "env_for", lambda module, settings: None)


def _argv(reads, out):
    return ["-1", reads[0], "-2", reads[1], "-o", str(out), "-t", "4", "--skip-bmtagger"]


def test_a_read_qc_dry_run_decompresses_nothing_and_marks_nothing(tmp_path, monkeypatch):
    calls = []
    _fake_tools(monkeypatch, calls)
    reads = _gz_reads(tmp_path)
    out = tmp_path / "qc"
    command.runner.configure(dry_run=True)
    assert read_qc.main(_argv(reads, out)) == 0
    assert calls, "the plan is still printed"
    left = set(os.listdir(out))
    assert not {"t_1.fastq", "t_2.fastq", "pre-QC_report", "post-QC_report"} & left
    assert not os.path.exists(os.path.join(str(out), ".metawrap2", "steps"))
    # The commands name the files the real run would decompress.
    assert str(out / "t_1.fastq") in calls[0]


def test_the_real_run_after_a_dry_run_starts_and_does_the_work(tmp_path, monkeypatch):
    """Preview, then run into the same -o: no refusal, and every step runs."""
    calls = []
    _fake_tools(monkeypatch, calls)
    reads = _gz_reads(tmp_path)
    out = tmp_path / "qc"
    command.runner.configure(dry_run=True)
    assert read_qc.main(_argv(reads, out)) == 0
    assert os.path.isfile(out / "run_config.json")

    warned = []
    monkeypatch.setattr(read_qc, "warning", warned.append)
    command.runner.configure(dry_run=False)
    calls.clear()
    assert read_qc.main(_argv(reads, out)) == 0
    assert (out / "final_pure_reads_1.fastq").is_file()
    assert any(c.startswith("fastp ") for c in calls)
    assert not [w for w in warned if "already exists" in w]


def test_resume_after_a_dry_run_runs_every_step(tmp_path, monkeypatch):
    calls = []
    _fake_tools(monkeypatch, calls)
    reads = _gz_reads(tmp_path)
    out = tmp_path / "qc"
    command.runner.configure(dry_run=True)
    assert read_qc.main(_argv(reads, out)) == 0
    command.runner.configure(dry_run=False, resume=True)
    calls.clear()
    assert read_qc.main(_argv(reads, out)) == 0
    assert any(c.startswith("fastp ") for c in calls)
    assert (out / "final_pure_reads_2.fastq").is_file()


def test_a_dry_run_keeps_the_record_of_a_real_run(tmp_path, monkeypatch):
    calls = []
    _fake_tools(monkeypatch, calls)
    reads = _gz_reads(tmp_path)
    out = tmp_path / "qc"
    assert read_qc.main(_argv(reads, out)) == 0
    before = (out / "run_config.json").read_text()
    command.runner.configure(dry_run=True)
    assert read_qc.main(_argv(reads, out)) == 0
    assert (out / "run_config.json").read_text() == before


# ── "already exists" ─────────────────────────────────────────────────────────────────────────


def test_a_fresh_output_folder_is_not_reported_as_existing(tmp_path, monkeypatch):
    calls, warned = [], []
    _fake_tools(monkeypatch, calls)
    monkeypatch.setattr(read_qc, "warning", warned.append)
    reads = _gz_reads(tmp_path)
    assert read_qc.main(_argv(reads, tmp_path / "fresh")) == 0
    assert not [w for w in warned if "already exists" in w]


def test_a_folder_with_earlier_results_is_reported(tmp_path, monkeypatch):
    calls, warned = [], []
    _fake_tools(monkeypatch, calls)
    monkeypatch.setattr(read_qc, "warning", warned.append)
    reads = _gz_reads(tmp_path)
    out = tmp_path / "used"
    out.mkdir()
    (out / "notes.txt").write_text("mine\n")
    assert read_qc.main(_argv(reads, out)) == 0
    assert [w for w in warned if "already exists" in w]


def test_own_logs_do_not_count_as_earlier_results(tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    for name in ("metawrap2.log", "metawrap2.log.1", "run.stdout", "run.stderr"):
        (out / name).write_text("")
    (out / ".metawrap2").mkdir()
    assert _common._has_prior_output(str(out)) is False
    (out / "final_pure_reads_1.fastq").write_text("")
    assert _common._has_prior_output(str(out)) is True
