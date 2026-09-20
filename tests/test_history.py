"""Tests for the durable run history (metawrap2.history) and `metawrap2 history`.

Two properties matter more than the formatting: a history write must never be able to fail a
run, and a run that was killed must still leave a trace. Both are pinned below.
"""

from __future__ import annotations

import json
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


# --- where the history lives ---------------------------------------------------------------


def test_default_history_path_is_under_the_user_directory(monkeypatch):
    monkeypatch.delenv("METAWRAP2_HISTORY", raising=False)
    assert history.history_path().endswith(os.path.join(".metawrap2", "history.jsonl"))
    assert "~" not in history.history_path()


def test_env_var_overrides_the_default(monkeypatch, tmp_path):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "h.jsonl"))
    assert history.history_path() == str(tmp_path / "h.jsonl")


def test_an_explicit_path_wins_over_the_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "env.jsonl"))
    assert history.history_path(str(tmp_path / "explicit.jsonl")).endswith("explicit.jsonl")


# --- recording -----------------------------------------------------------------------------


def test_record_writes_one_json_line_per_call(hist):
    add(hist, run_id="r1")
    add(hist, run_id="r2")
    with open(hist) as fh:
        lines = fh.read().splitlines()
    assert len(lines) == 2
    assert all(json.loads(line)["module"] == "binning" for line in lines)


def test_record_creates_the_parent_directory(tmp_path):
    target = str(tmp_path / "deep" / "dir" / "history.jsonl")
    add(target)
    assert os.path.isfile(target)


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


def test_duration_is_rounded_and_optional(hist):
    add(hist, status=history.COMPLETED, duration_seconds=123.4567)
    assert next(history.read(hist))["duration_seconds"] == 123.5


# --- reading -------------------------------------------------------------------------------


def test_read_of_a_missing_file_is_empty(tmp_path):
    assert list(history.read(str(tmp_path / "nope.jsonl"))) == []


def test_read_skips_a_truncated_final_line(hist):
    """An append interrupted mid-line costs that record, not the whole history."""
    add(hist, run_id="r1")
    add(hist, run_id="r2")
    with open(hist, "a") as fh:
        fh.write('{"module": "bin')
    assert [e["run_id"] for e in history.read(hist)] == ["r1", "r2"]


def test_read_skips_blank_lines_and_non_objects(hist):
    add(hist, run_id="r1")
    with open(hist, "a") as fh:
        fh.write('\n\n["not an object"]\n')
    assert len(list(history.read(hist))) == 1


def test_entries_filters_by_module_status_and_output(hist, tmp_path):
    out = str(tmp_path / "study")
    add(hist, module="binning", run_id="r1", status=history.COMPLETED, output=out)
    add(hist, module="assembly", run_id="r2", status=history.FAILED)
    assert [e["module"] for e in history.entries(hist, module="binning")] == ["binning"]
    assert [e["module"] for e in history.entries(hist, status=history.FAILED)] == ["assembly"]
    assert [e["module"] for e in history.entries(hist, output=out)] == ["binning"]


def test_entries_filters_by_date(hist):
    add(hist, run_id="r1")
    assert history.entries(hist, since="1999-01-01")
    assert history.entries(hist, since="2999-01-01") == []


# --- runs(): folding start and finish together ---------------------------------------------


def test_runs_merges_a_start_and_its_outcome(hist):
    add(hist, run_id="abc", status=history.STARTED)
    add(hist, run_id="abc", status=history.COMPLETED, duration_seconds=90)
    rows = history.runs(hist)
    assert len(rows) == 1
    assert rows[0]["status"] == history.COMPLETED
    assert rows[0]["duration_seconds"] == 90
    assert rows[0]["finished"]


def test_a_killed_run_still_leaves_a_started_row(hist):
    """No outcome entry was ever written - the row must survive anyway, visibly unfinished."""
    add(hist, run_id="killed", status=history.STARTED)
    rows = history.runs(hist)
    assert len(rows) == 1 and rows[0]["status"] == history.STARTED


def test_runs_keeps_chronological_order(hist):
    for i, module in enumerate(["read_qc", "assembly", "binning"]):
        add(hist, module=module, run_id="r%d" % i, status=history.COMPLETED)
    assert [r["module"] for r in history.runs(hist)] == ["read_qc", "assembly", "binning"]


# --- the CLI -------------------------------------------------------------------------------


def test_cli_explains_itself_when_there_is_no_history(tmp_path, capsys):
    target = str(tmp_path / "nope.jsonl")
    assert history_cmd.main(["--path", target]) == 0
    out = capsys.readouterr().out
    assert "No history yet" in out and target in out


def test_cli_says_so_when_filters_match_nothing(hist, capsys):
    add(hist, module="binning", status=history.COMPLETED)
    assert history_cmd.main(["--path", hist, "--module", "assembly"]) == 0
    assert "match those filters" in capsys.readouterr().out


def test_cli_prints_a_table(hist, capsys):
    add(hist, module="binning", run_id="r1", status=history.STARTED)
    add(hist, module="binning", run_id="r1", status=history.COMPLETED, duration_seconds=3600)
    assert history_cmd.main(["--path", hist]) == 0
    out = capsys.readouterr().out
    assert "binning" in out and "ok" in out and "60m" in out
    assert "1 run(s) shown" in out


def test_cli_failed_shows_only_unsuccessful_runs(hist, capsys):
    add(hist, module="binning", run_id="r1", status=history.COMPLETED)
    add(hist, module="assembly", run_id="r2", status=history.FAILED)
    add(hist, module="read_qc", run_id="r3", status=history.INTERRUPTED)
    assert history_cmd.main(["--path", hist, "--failed"]) == 0
    out = capsys.readouterr().out
    assert "assembly" in out and "read_qc" in out and "binning" not in out


def test_cli_limit_keeps_the_most_recent(hist, capsys):
    for i in range(5):
        add(hist, module="mod%d" % i, run_id="r%d" % i, status=history.COMPLETED)
    assert history_cmd.main(["--path", hist, "-n", "2"]) == 0
    out = capsys.readouterr().out
    assert "mod4" in out and "mod3" in out and "mod0" not in out


def test_cli_limit_zero_shows_everything(hist, capsys):
    for i in range(3):
        add(hist, module="mod%d" % i, run_id="r%d" % i, status=history.COMPLETED)
    assert history_cmd.main(["--path", hist, "-n", "0"]) == 0
    assert "3 run(s) shown" in capsys.readouterr().out


def test_cli_commands_prints_a_pasteable_command_line(hist, capsys):
    history.record(
        "binning",
        version="2.1.0",
        run_id="r1",
        status=history.COMPLETED,
        command=["metawrap2", "binning", "-o", "out dir"],
        path=hist,
    )
    assert history_cmd.main(["--path", hist, "--commands"]) == 0
    assert "metawrap2 binning -o 'out dir'" in capsys.readouterr().out


def test_cli_json_is_machine_readable(hist, capsys):
    add(hist, module="binning", status=history.COMPLETED)
    assert history_cmd.main(["--path", hist, "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert rows[0]["module"] == "binning"


def test_cli_stats_totals_per_module(hist, capsys):
    add(hist, module="binning", run_id="r1", status=history.COMPLETED, duration_seconds=60)
    add(hist, module="binning", run_id="r2", status=history.FAILED, duration_seconds=30)
    add(hist, module="assembly", run_id="r3", status=history.COMPLETED, duration_seconds=10)
    assert history_cmd.main(["--path", hist, "--stats"]) == 0
    out = capsys.readouterr().out
    lines = [ln.split() for ln in out.splitlines() if ln.startswith("  binning")]
    assert lines and lines[0][1:4] == ["2", "1", "1"]
    assert "all" in out
