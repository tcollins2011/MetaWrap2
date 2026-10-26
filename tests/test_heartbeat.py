"""The progress heartbeat printed while a step is silent.

Timings use a fake clock. The test that runs a real command uses an interval shorter than
interpreter start-up, so it does not sleep.
"""

from __future__ import annotations

import os
from typing import List, Optional

from metawrap2.command import CommandRunner
from metawrap2.logfilter import Heartbeat, LogFilter, duration

#: A COMEBin training line.
EPOCH = "2024-01-01 00:00:00,000 - NN training: epoch=%d/200, batches=27, loss=3.45, top1=83.06"


class FakeClock:
    """A monotonic clock the test advances by hand."""

    def __init__(self) -> None:
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _beats(beat: Heartbeat, clock: FakeClock, seconds: float, step: float) -> List[str]:
    """Every heartbeat produced over *seconds* of silence, polled every *step* seconds."""
    out: List[str] = []
    for _ in range(int(seconds / step)):
        clock.advance(step)
        line: Optional[str] = beat.due()
        if line is not None:
            out.append(line.rstrip("\n"))
    return out


def test_a_step_that_has_printed_nothing_says_so_rather_than_nothing():
    """With no output yet, the heartbeat says the tool has printed nothing."""
    clock = FakeClock()
    beat = Heartbeat(step="comebin", interval=60.0, clock=clock)
    assert beat.due() is None
    clock.advance(60.0)
    line = beat.due()
    assert line == (
        "[metawrap2] comebin is still running, 1m00s elapsed, and has printed nothing at all yet\n"
    )


def test_a_heartbeat_carries_the_elapsed_time_and_what_the_tool_last_said():
    """The heartbeat gives the elapsed time, the silence and the tool's last line."""
    clock = FakeClock()
    beat = Heartbeat(step="comebin", interval=60.0, clock=clock)
    beat.note(EPOCH % 85 + "\n")
    clock.advance(4 * 3600 + 10 * 60)
    line = beat.due()
    assert line is not None
    assert "comebin is still running, 4h10m elapsed" in line
    assert "epoch=85/200" in line
    assert "silent for 4h10m" in line


def test_heartbeats_do_not_accumulate_and_report_the_whole_silence():
    """One heartbeat per interval, each reporting the silence since the tool's last line."""
    clock = FakeClock()
    beat = Heartbeat(step="SemiBin2", interval=60.0, clock=clock)
    beat.note("Starting clustering\n")
    lines = _beats(beat, clock, seconds=300.0, step=5.0)
    assert len(lines) == 5
    assert "silent for 1m00s" in lines[0]
    assert "silent for 5m00s" in lines[-1]
    # The quoted line is still the last thing the tool said, not the previous heartbeat.
    assert all("Starting clustering" in line for line in lines)


def test_a_tool_that_keeps_printing_is_never_called_silent():
    """A tool printing a line every 27 s gets no heartbeat, even when its lines are collapsed."""
    clock = FakeClock()
    beat = Heartbeat(step="comebin", interval=60.0, clock=clock)
    filt = LogFilter(clock=clock, heartbeat=beat)
    heartbeats = []
    for epoch in range(1, 201):
        filt.feed(EPOCH % epoch + "\n")
        clock.advance(27.0)
        line = beat.due()
        if line is not None:
            heartbeats.append(line)
    assert heartbeats == []


def test_an_interval_of_zero_turns_heartbeats_off():
    """An interval of 0 disables heartbeats."""
    clock = FakeClock()
    beat = Heartbeat(step="comebin", interval=0.0, clock=clock)
    assert not beat.enabled
    assert _beats(beat, clock, seconds=3600.0, step=60.0) == []


def test_the_poll_interval_is_shorter_than_the_heartbeat_interval():
    """The poll interval is at most half the heartbeat interval."""
    assert Heartbeat(interval=60.0).tick_seconds <= 60.0 / 2


def test_a_heartbeat_changes_nothing_in_the_captured_log(tmp_path, capsys):
    """Heartbeats appear on the screen and not in run.stderr."""
    log = tmp_path / "run.stderr"
    runner = CommandRunner(
        run_stderr_path=str(log),
        # Shorter than interpreter start-up, so a heartbeat is due before the first line.
        heartbeat_seconds=0.001,
    )
    assert runner.run(["python3", "-c", "import sys; sys.stderr.write('working\\n')"]) == 0

    shown = capsys.readouterr().err
    assert "is still running" in shown, "a silent step produced no heartbeat at all"
    assert "working" in shown
    written = log.read_text()
    assert "working" in written
    assert "is still running" not in written


def test_the_filter_passes_the_same_lines_through_with_a_heartbeat_attached():
    """A filter with a heartbeat attached emits the same lines as one without."""
    lines = [EPOCH % epoch + "\n" for epoch in range(1, 60)]
    plain = LogFilter(clock=FakeClock())
    withbeat = LogFilter(clock=FakeClock(), heartbeat=Heartbeat(interval=1.0, clock=FakeClock()))
    assert [plain.feed(line) for line in lines] == [withbeat.feed(line) for line in lines]


def test_verbose_logs_still_gets_a_heartbeat():
    """With --verbose-logs (filtering off) the heartbeat still sees each line."""
    clock = FakeClock()
    beat = Heartbeat(step="comebin", interval=60.0, clock=clock)
    filt = LogFilter(enabled=False, clock=clock, heartbeat=beat)
    filt.feed(EPOCH % 85 + "\n")
    clock.advance(120.0)
    line = beat.due()
    assert line is not None and "epoch=85/200" in line


def test_elapsed_time_is_written_the_way_it_is_read():
    """Elapsed times are formatted as hours and minutes, minutes and seconds, or seconds."""
    assert duration(4 * 3600 + 10 * 60) == "4h10m"
    assert duration(75 * 60) == "1h15m"
    assert duration(75) == "1m15s"
    assert duration(9) == "9s"
    # A clock that went backwards between two reads must not print a negative age.
    assert duration(-5) == "0s"


def test_the_interval_comes_from_the_config(tmp_path, monkeypatch):
    """``[settings] heartbeat_seconds`` sets the interval; set_heartbeat_seconds overrides it."""
    from metawrap2 import command

    monkeypatch.setattr(command, "_CONFIGURED_HEARTBEAT", None)
    config = tmp_path / "config.toml"
    config.write_text("[settings]\nheartbeat_seconds = 120\n")
    monkeypatch.setenv("METAWRAP2_CONFIG", str(config))
    assert command.configured_heartbeat_seconds() == 120.0

    monkeypatch.setattr(command, "_CONFIGURED_HEARTBEAT", None)
    config.write_text("[settings]\nheartbeat_seconds = 0\n")
    assert command.configured_heartbeat_seconds() == 0.0
    assert not Heartbeat(interval=0.0).enabled

    # A process-level override, like the other display settings.
    runner = CommandRunner()
    assert runner.heartbeat_seconds is None  # unset means "ask the config"
    monkeypatch.setattr(command, "runner", runner)
    command.set_heartbeat_seconds(15)
    assert runner.heartbeat_seconds == 15.0
    command.set_heartbeat_seconds(None)
    assert runner.heartbeat_seconds is None


def test_the_heartbeat_does_not_ask_whether_the_output_is_a_terminal():
    """Heartbeats do not depend on a terminal, as with the collapser's refresh lines."""
    source = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src", "metawrap2")
    with open(os.path.join(source, "logfilter.py")) as handle:
        assert "isatty" not in handle.read()


def test_the_global_flag_sets_the_interval_and_is_removed_from_the_argv():
    """`--heartbeat-seconds N` and `--heartbeat-seconds=N` set the interval before the module name."""
    from metawrap2 import cli, command

    try:
        assert cli._apply_global_flags(["--heartbeat-seconds", "15", "binning"]) == ["binning"]
        assert command.runner.heartbeat_seconds == 15.0
        assert cli._apply_global_flags(["--heartbeat-seconds=0", "binning"]) == ["binning"]
        assert command.runner.heartbeat_seconds == 0.0
    finally:
        command.set_heartbeat_seconds(None)


def test_a_heartbeat_interval_that_is_not_a_non_negative_number_is_refused():
    """A non-number or a negative interval is an error; 0 turns the heartbeat off."""
    import pytest

    from metawrap2 import cli, command
    from metawrap2.limits import LimitsError

    try:
        for bad in ("abc", "-5"):
            with pytest.raises(LimitsError):
                cli._apply_global_flags(["--heartbeat-seconds", bad, "binning"])
    finally:
        command.set_heartbeat_seconds(None)
