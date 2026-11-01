"""The hint printed when GTDB-Tk's classify_wf fails.

GTDB-Tk logs "Using GTDB-Tk reference data version ..." once it has found its data. A failure after
that line points at the environment, so the hint must not send the user to the database.
"""

from __future__ import annotations

import os
import types

import pytest

from metawrap2.command import ToolError
from metawrap2.modules import classify_bins as cb

LOADED = (
    "[18:27:43] INFO: GTDB-Tk v2.4.1\n"
    "[18:27:43] INFO: Using GTDB-Tk reference data version r226: /path/to/GTDBTK\n"
    "[18:27:43] ERROR: Uncontrolled exit resulting from an unexpected error.\n"
)
NOT_LOADED = (
    "[18:27:43] INFO: GTDB-Tk v2.4.1\n"
    "[18:27:43] ERROR: The GTDB-Tk reference data does not exist or is corrupted.\n"
)


def test_a_failure_after_the_data_loaded_points_at_the_environment(tmp_path):
    hint = cb._gtdbtk_hint(LOADED, str(tmp_path))
    assert "not a database problem" in hint
    assert "metawrap2 install-env gtdbtk --force" in hint
    assert "GTDBTK_DATA_PATH" not in hint


def test_a_failure_before_the_data_loaded_keeps_the_database_hint(tmp_path):
    assert cb._gtdbtk_hint(NOT_LOADED, str(tmp_path)) == cb.GTDBTK_DB_HINT
    assert "GTDBTK_DATA_PATH" in cb.GTDBTK_DB_HINT


def test_the_line_is_found_in_gtdbtks_own_log_when_the_tail_lost_it(tmp_path):
    """A long traceback pushes the reference-data line out of the 25-line error tail."""
    (tmp_path / "gtdbtk.log").write_text(LOADED)
    assert cb._gtdbtk_hint("Traceback ...\nValueError: x\n", str(tmp_path)) == cb.GTDBTK_ENV_HINT


def test_a_previous_runs_log_lines_are_ignored(tmp_path):
    """GTDB-Tk appends to gtdbtk.log, so only what this run wrote counts."""
    log = tmp_path / "gtdbtk.log"
    log.write_text(LOADED)
    offset = os.path.getsize(log)
    with open(log, "a") as fh:
        fh.write(NOT_LOADED)
    assert cb._gtdbtk_hint("", str(tmp_path), offset) == cb.GTDBTK_DB_HINT


def test_classify_with_gtdbtk_replaces_the_hint_on_the_raised_error(tmp_path, monkeypatch):
    monkeypatch.setattr(cb, "optional_env_for", lambda *a, **k: "metawrap2-gtdbtk")
    monkeypatch.setattr(cb, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(cb, "comm", lambda *a, **k: None)

    def failing_run(cmd, **kw):
        assert kw["hint"] == cb.GTDBTK_DB_HINT
        raise ToolError(
            tool="gtdbtk", returncode=1, cmd=["gtdbtk"], log_tail=LOADED, hint=kw["hint"]
        )

    monkeypatch.setattr(cb, "run", failing_run)
    args = types.SimpleNamespace(
        bins=str(tmp_path / "bins"), output=str(tmp_path / "out"), threads=1
    )
    ckpt = types.SimpleNamespace(todo=lambda step: True, done=lambda step: None)
    with pytest.raises(ToolError) as caught:
        cb._classify_with_gtdbtk(args, {}, ckpt)
    assert caught.value.hint == cb.GTDBTK_ENV_HINT
    assert "not a database problem" in str(caught.value)
