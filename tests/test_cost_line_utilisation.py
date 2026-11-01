"""The cost line: no cores-busy figure on very short steps, and busy cores against the grant."""

from __future__ import annotations

from metawrap2 import command, usage
from metawrap2.modules import _common


def _measured(wall, cpu):
    return {
        "wall_seconds": wall,
        "peak_rss_bytes": 1024**2,
        "user_seconds": cpu,
        "system_seconds": 0,
    }


def test_a_short_step_prints_na_for_cores_busy():
    """A 1 s run reported 12.8 cores busy; the ratio is rounding at that length."""
    line = usage.summarise(_measured(1, 12.8), threads=32)
    assert "cores busy n/a" in line
    assert "12.8" not in line
    assert "n/a" in usage.summarise(_measured(0, 3.0))


def test_busy_cores_are_given_against_the_grant():
    line = usage.summarise(_measured(100, 1420), threads=32)
    assert "14.2 of 32 cores busy (44%)" in line


def test_without_a_grant_the_plain_figure_is_kept():
    assert "14.2 cores busy" in usage.summarise(_measured(100, 1420))


def test_finish_run_passes_the_modules_thread_grant(tmp_path, monkeypatch, capsys):
    rec = _common.RunRecorder(
        module="read_qc",
        version="test",
        parameters={"threads": 8},
        config_file=None,
        conda_env=None,
        databases={},
        inputs=[],
    )
    monkeypatch.setattr(_common, "self_usage", lambda r: _measured(60, 240), raising=True)
    monkeypatch.delenv(_common.DRIVEN_BY_RUN_ENV, raising=False)
    command.set_dry_run(True)
    try:
        _common.finish_run(rec, str(tmp_path / "out"))
    finally:
        command.set_dry_run(False)
    assert "4.0 of 8 cores busy (50%)" in capsys.readouterr().out
