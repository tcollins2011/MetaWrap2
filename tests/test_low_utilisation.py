"""The mid-step low-utilisation warning and the /proc CPU-time parsing it uses.

Thresholds are tested with a fake clock. The tests marked ``needs_proc`` read the real /proc.
"""

from __future__ import annotations

import os
import subprocess
import threading
from typing import Dict, List, Optional, Tuple

import pytest

from metawrap2 import usage

#: CPU time is read from /proc, which only Linux has.
needs_proc = pytest.mark.skipif(not os.path.isdir("/proc"), reason="needs /proc")

#: The thread grant used in these tests.
CORES = 32


class FakeClock:
    def __init__(self) -> None:
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now

    def advance(self, seconds: float) -> None:
        self.now += seconds


def _measured(busy: float, elapsed: float, threads: int = CORES, **kw) -> Optional[str]:
    """What a step with *threads* granted is told after *elapsed* seconds at *busy* cores."""
    clock = FakeClock()
    util = usage.Utilisation(threads=threads, step="comebin", clock=clock, **kw)
    clock.advance(elapsed)
    return util.note(busy * elapsed)


def test_the_warning_carries_the_measured_numbers_and_the_conclusion():
    """The warning gives cores busy, the grant and the elapsed time, and advises -j over -t."""
    line = _measured(9.0, 4 * 3600 + 10 * 60)
    assert line is not None
    assert "9.0 of its 32 granted cores busy over 4h10m" in line
    assert "does not parallelise" in line
    assert "will not shorten it" in line
    assert "(-j)" in line and "(-t)" in line
    # The step is named, for runs with several concurrent steps.
    assert line.startswith("[metawrap2] comebin ")


def test_semibin2s_clustering_phase_is_reported_too():
    """A shorter step (3.5 of 32 cores over 1 h 15 m) is also reported."""
    line = _measured(3.5, 75 * 60)
    assert line is not None
    assert "3.5 of its 32 granted cores busy over 1h15m" in line


def test_a_step_too_short_to_judge_is_left_alone():
    """No warning before five minutes, while tools are in their single-threaded start-up."""
    assert _measured(1.0, 240.0) is None
    assert _measured(1.0, 299.0) is None
    # At five minutes it fires.
    assert _measured(1.0, 300.0) is not None


def test_a_step_using_its_grant_is_left_alone():
    """A long step that keeps more than a third of its grant busy gets no warning."""
    assert _measured(30.0, 9 * 3600) is None
    assert _measured(20.0, 9 * 3600) is None
    # The threshold is a third of the grant.
    assert _measured(CORES * 0.34, 3600.0) is None
    assert _measured(CORES * 0.32, 3600.0) is not None


def test_it_warns_once_per_step_and_not_once_per_sample():
    """The warning is printed once per step, not at every 5 s sample."""
    clock = FakeClock()
    util = usage.Utilisation(threads=CORES, step="comebin", clock=clock)
    clock.advance(600.0)
    lines = [util.note(9.0 * (600.0 + tick * 5)) for tick in range(100)]
    assert len([line for line in lines if line is not None]) == 1


def test_zero_seconds_turns_the_warning_off():
    """after_seconds = 0 disables the warning."""
    clock = FakeClock()
    util = usage.Utilisation(threads=CORES, after_seconds=0.0, clock=clock)
    assert not util.enabled
    clock.advance(9 * 3600)
    assert util.note(1.0) is None


def test_a_one_thread_grant_is_never_underused():
    """A one-thread grant never warns."""
    assert not usage.Utilisation(threads=1).enabled
    assert _measured(0.2, 3600.0, threads=1) is None


def test_the_thresholds_come_from_the_config(tmp_path, monkeypatch):
    """Both thresholds are read from ``[settings]``, and 0 turns the warning off."""
    config = tmp_path / "config.toml"
    config.write_text(
        "[settings]\nlow_utilisation_after_seconds = 30\nlow_utilisation_fraction = 0.9\n"
    )
    monkeypatch.setenv("METAWRAP2_CONFIG", str(config))
    util = usage.configured_utilisation(CORES, "comebin")
    assert util is not None
    assert util.after_seconds == 30.0
    assert util.fraction == 0.9

    config.write_text("[settings]\nlow_utilisation_after_seconds = 0\n")
    off = usage.configured_utilisation(CORES, "comebin")
    assert off is not None and not off.enabled


def test_the_grant_is_read_off_the_step_that_was_given_it():
    """The grant is the last ``-t`` in the step's argv, which is the one ``run`` appends."""
    argv = ["python", "-m", "metawrap2.cli", "binning", "-a", "asm.fa", "-t", "48"]
    assert usage.grant_from_argv(argv) == 48
    assert usage.grant_from_argv(argv + ["-t", "12"]) == 12
    assert usage.grant_from_argv(["python", "-m", "metawrap2.cli", "binning"]) is None
    # -t all is resolved to a number before a step is built; a non-numeric value is ignored.
    assert usage.grant_from_argv(["metawrap2", "binning", "-t", "all"]) is None


def _rows(*triples: Tuple[int, Optional[int], Optional[float]]) -> Dict:
    return {pid: (ppid, cpu) for pid, ppid, cpu in triples}


def test_a_trees_cpu_is_summed_without_counting_a_reaped_child_twice():
    """CPU time is summed over the whole tree under the step, excluding other processes."""
    rows = _rows(
        (100, 1, 2.0),  # the step: mamba
        (200, 100, 30.0),  # the tool
        (300, 200, 8.0),  # a worker
        (400, 1, 999.0),  # another tenant's job, not in the tree
    )
    assert usage.tree_cpu_seconds(100, rows) == 40.0
    assert usage.tree_cpu_seconds(200, rows) == 38.0
    # A root missing from the snapshot (the step has exited) gives None.
    assert usage.tree_cpu_seconds(500, rows) is None


def test_an_unreadable_process_is_skipped_rather_than_counted_as_zero():
    """A process whose CPU time cannot be read is skipped."""
    rows = _rows((100, 1, 5.0), (200, 100, None), (300, 200, 5.0))
    assert usage.tree_cpu_seconds(100, rows) == 10.0


def test_the_columns_survive_a_comm_full_of_spaces_and_brackets():
    """ppid and utime..cstime are read from the right columns when comm has spaces and parens."""
    # A comm with spaces and a nested close-paren before the real one.
    fields = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
    line = "4242 (my (weird) proc) S " + fields
    close = line.rfind(")")
    parsed = line[close + 2 :].split()

    # Index 0 is proc(5) field 3, the state; index 1 is the ppid.
    assert parsed[0] == "S"
    assert usage._ppid_from(parsed) == 1
    # utime..cstime are indices 11..14 of this slice, i.e. the values 11+12+13+14 = 50 ticks.
    assert usage._cpu_seconds_from(parsed) == pytest.approx(50 / usage.CLOCK_TICKS)


@needs_proc
def test_the_live_process_parses_and_its_ppid_is_the_real_one():
    """The live process's /proc row gives ``os.getppid()`` and a non-negative CPU time."""
    parsed = usage._stat_fields(os.getpid())
    assert len(parsed) >= 15
    assert usage._ppid_from(parsed) == os.getppid()
    own = usage._cpu_seconds_from(parsed)
    assert own is not None and own >= 0.0


@needs_proc
def test_a_childs_cpu_counts_both_while_it_runs_and_after_it_is_reaped():
    """A child's CPU time is counted while it runs and after it has been reaped."""
    before = usage.tree_cpu_seconds(os.getpid())
    assert before is not None
    burn = "x = 0\nfor i in range(20000000): x += i\n"
    child = subprocess.Popen(["python3", "-c", burn])
    running = usage.tree_cpu_seconds(os.getpid())
    assert child.wait() == 0
    reaped = usage.tree_cpu_seconds(os.getpid())
    assert running is not None and reaped is not None
    # The child burned at least a second of CPU; the total must still hold it once it is gone.
    assert reaped >= before + 0.5
    assert reaped >= running


@needs_proc
def test_the_warning_comes_out_of_the_sampler_that_was_already_running():
    """The warning comes from the StepWatcher sampler (GpuWatcher is the same class)."""
    assert usage.GpuWatcher is usage.StepWatcher

    clock = FakeClock()
    util = usage.Utilisation(threads=CORES, step="comebin", clock=clock)
    clock.advance(600.0)
    seen: List[str] = []
    arrived = threading.Event()

    def collect(line: str) -> None:
        seen.append(line)
        arrived.set()

    watcher = usage.StepWatcher(os.getpid(), interval=0.01, utilisation=util, on_warning=collect)
    watcher.start()
    try:
        assert arrived.wait(10.0), "the sampler produced no warning"
    finally:
        watcher.stop()
    assert "granted cores busy" in seen[0]


def test_a_standalone_module_watches_cores_on_a_machine_with_no_gpu(tmp_path, monkeypatch):
    """Without a GPU, a module run on its own still starts the sampler for its thread grant."""
    import argparse
    from typing import Set

    from metawrap2.config import load_settings
    from metawrap2.modules import _common

    started = {}

    class FakeWatcher:
        def __init__(self, pid, utilisation=None, **kwargs):
            started["utilisation"] = utilisation
            # finish_run reads both when it writes the cost line.
            self.peak_mib = 0
            self.devices: Set[str] = set()

        def start(self):
            return self

        def stop(self):
            return None

    monkeypatch.setattr(_common._usage, "has_gpu", lambda: False)
    monkeypatch.setattr(_common._usage, "StepWatcher", FakeWatcher)
    monkeypatch.delenv(_common.DRIVEN_BY_RUN_ENV, raising=False)
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))

    args = argparse.Namespace(output=str(tmp_path / "out"), threads=CORES, config=None)
    rec = _common.start_run("binning", args, None, load_settings(None), inputs=[])
    try:
        assert started.get("utilisation") is not None, "no core sampling without a GPU"
        assert started["utilisation"].threads == CORES
        assert started["utilisation"].step == "binning"
    finally:
        _common.finish_run(rec, args.output)


def test_nothing_is_sampled_for_a_single_threaded_step():
    """A single-threaded step on a CPU-only machine gets no utilisation sampler."""
    assert usage.configured_utilisation(1, "read_qc") is None
    assert usage.configured_utilisation(None, "read_qc") is None
