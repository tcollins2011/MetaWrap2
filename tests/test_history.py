"""Tests for the durable run history (metawrap2.history) and `metawrap2 history`.

Two properties matter more than the formatting: a history write must never be able to fail a
run, and a run that was killed must still leave a trace. Both are pinned below.
"""

from __future__ import annotations

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


# --- runs(): folding start and finish together ---------------------------------------------


def test_a_killed_run_still_leaves_a_started_row(hist, tmp_path, capsys):
    """No outcome entry was ever written - the row must survive anyway, visibly unfinished."""
    add(hist, run_id="killed", status=history.STARTED)
    rows = history.runs(hist)
    assert len(rows) == 1 and rows[0]["status"] == history.STARTED

    # Losing a history line must never be the reason a six-hour assembly fails.
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory")
    add(str(blocker / "history.jsonl"))  # parent is a file: makedirs/open both fail

    # `metawrap2 history` folds a start and a finish into one row.
    add(hist, module="binning", run_id="r1", status=history.STARTED)
    add(hist, module="binning", run_id="r1", status=history.COMPLETED, duration_seconds=3600)
    assert history_cmd.main(["--path", hist]) == 0
    out = capsys.readouterr().out
    assert "binning" in out and "ok" in out and "60m" in out
