"""Tests for the durable run history (metawrap2.history) and `metawrap2 history`.

Two properties matter more than the formatting: a history write must never be able to fail a
run, and a run that was killed must still leave a trace. Both are pinned below.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import history
from metawrap2.commands import history_cmd


@pytest.fixture
def hist(tmp_path, monkeypatch):
    """An isolated history file. Never touch the developer's real ~/.metawrap2."""
    monkeypatch.delenv("METAWRAP2_HISTORY", raising=False)
    return str(tmp_path / "history.jsonl")


def add(path, module="binning", status=history.STARTED, run_id="r1", **kw):
    history.record(
        module,
        version="2.1.0",
        run_id=run_id,
        command=["metawrap2", module],
        status=status,
        path=path,
        **kw,
    )


# --- recording -----------------------------------------------------------------------------


def test_record_captures_the_context_a_user_needs(hist):
    add(
        hist,
        output=".",
        databases={"CHECKM_DB": "/db/checkm", "BLASTDB": ""},
        conda_env="metawrap2-binning",
    )
    entry = next(history.read(hist))
    assert entry["module"] == "binning"
    assert entry["version"] == "2.1.0"
    assert entry["command"] == ["metawrap2", "binning"]
    assert entry["output"] == os.path.abspath(".")
    assert entry["conda_env"] == "metawrap2-binning"
    assert entry["databases"] == {"CHECKM_DB": "/db/checkm"}  # blank values dropped
    assert entry["time"] and entry["user"] and entry["host"] and entry["python"]


def test_record_never_raises_when_the_history_cannot_be_written(tmp_path):
    """Losing a history line must never be the reason a six-hour assembly fails."""
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory")
    add(str(blocker / "history.jsonl"))  # parent is a file: makedirs/open both fail


# --- runs(): folding start and finish together ---------------------------------------------


def test_a_killed_run_still_leaves_a_started_row(hist):
    """No outcome entry was ever written - the row must survive anyway, visibly unfinished."""
    add(hist, run_id="killed", status=history.STARTED)
    rows = history.runs(hist)
    assert len(rows) == 1 and rows[0]["status"] == history.STARTED


# --- the CLI -------------------------------------------------------------------------------


def test_cli_prints_a_table(hist, capsys):
    add(hist, module="binning", run_id="r1", status=history.STARTED)
    add(hist, module="binning", run_id="r1", status=history.COMPLETED, duration_seconds=3600)
    assert history_cmd.main(["--path", hist]) == 0
    out = capsys.readouterr().out
    assert "binning" in out and "ok" in out and "60m" in out
    assert "1 run(s) shown" in out
