"""Tests for metawrap2.limits: one resource budget for a whole run."""

from __future__ import annotations

import pytest

from metawrap2 import limits


@pytest.fixture(autouse=True)
def _no_ambient_overrides():
    # Overrides are process-global by design (they come from flags parsed before any module runs),
    # so a test that sets one must not leak it into the next. Same for the cached GPU count.
    limits.clear_overrides()
    limits._count_gpus.cache_clear()
    yield
    limits.clear_overrides()
    limits._count_gpus.cache_clear()


def machine(cores: int = 64, memory_gb: int = 256, gpus: int = 0, monkeypatch=None):
    """Pin the measured machine so the arithmetic is testable on any host."""
    monkeypatch.setattr(limits, "machine_cores", lambda: cores)
    monkeypatch.setattr(limits, "machine_memory_bytes", lambda: memory_gb * 1024**3)
    monkeypatch.setattr(limits, "machine_gpus", lambda: gpus)


# --- measuring the machine ----------------------------------------------------------------


@pytest.mark.real_machine
def test_cores_never_zero():
    assert limits.machine_cores() >= 1
    assert limits.machine_memory_bytes() >= 1


# --- resolving the budget -----------------------------------------------------------------


def test_fractions_are_applied_to_the_measured_machine(monkeypatch):
    machine(cores=64, memory_gb=100, monkeypatch=monkeypatch)
    budget = limits.resolve()
    assert budget.threads == 32
    assert budget.memory_gb == pytest.approx(80.0)

    # Absolutes win over fractions.
    budget = limits.resolve(limits.LimitsConfig(max_threads=10, max_memory_gb=20))
    assert budget.threads == 10
    assert budget.memory_gb == pytest.approx(20.0)

    # Absolutes are still clamped to the machine: asking for more than exists is a mistake to
    # correct, not a request to honour into an OOM.
    machine(cores=8, memory_gb=16, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=64, max_memory_gb=256))
    assert budget.threads == 8
    assert budget.memory_gb == pytest.approx(16.0)

    # Per-worker shares divide the budget rather than multiplying it: eight workers add up to
    # the budget, not eight times it.
    machine(cores=96, memory_gb=200, monkeypatch=monkeypatch)
    budget = limits.resolve(limits.LimitsConfig(max_threads=96, max_memory_gb=160))
    worker = budget.per_worker(8)
    assert worker.threads == 12
    assert worker.memory_gb == pytest.approx(20.0)
    assert worker.threads * 8 <= budget.threads
    assert worker.memory_bytes * 8 <= budget.memory_bytes


# --- the disk ceiling ---------------------------------------------------------------------


def fake_usage(monkeypatch, total: int, free: int):
    monkeypatch.setattr(
        limits.shutil,
        "disk_usage",
        lambda _path: type("Usage", (), {"total": total, "free": free, "used": total - free})(),
    )


# --- wiring: the thread, memory and disk ceilings ---------------------------------------


class Args:
    def __init__(self, threads=None):
        self.threads = threads


class MemArgs:
    def __init__(self, memory=None):
        self.memory = memory


def test_resolve_threads_caps_to_the_budget(monkeypatch, tmp_path):
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    machine(cores=64, monkeypatch=monkeypatch)
    settings = Settings(limits=limits.LimitsConfig(max_threads=8))
    args = Args(threads=32)
    assert _common.resolve_threads(args, settings) == 8
    # Modules read args.threads directly after this, so the cap has to land there too.
    assert args.threads == 8

    # The memory ceiling is wired the same way.
    machine(memory_gb=100, monkeypatch=monkeypatch)
    mem_args = MemArgs(memory=500)
    assert _common.resolve_memory_gb(mem_args, Settings()) == 80
    assert mem_args.memory == 80

    # And the disk ceiling. Unlike the size estimate, which only warns, this one is measured - so
    # there is nothing to be wrong about and no reason to start a step the monitor would stop
    # minutes later: a full filesystem stops the step before it starts.
    fake_usage(monkeypatch, total=1000, free=10)
    with pytest.raises(SystemExit):
        _common.check_disk_space("assembly", Settings(), [], str(tmp_path))


def test_a_monitor_stopped_mid_sample_does_not_interrupt_the_next_step(monkeypatch):
    """stop() has to win against a watcher already inside sample().

    sample() is two statvfs calls and an os.stat, which on a stalled mount takes seconds. A step
    that finished in that window has already called stop(), so interrupting afterwards delivers a
    KeyboardInterrupt into whatever the process does next - under `run` or `cohort`, the following
    step, which is then aborted and recorded as interrupted while the step that was being watched
    is reported as a disk stop it never hit.
    """
    import threading

    interrupts = []
    monkeypatch.setattr("_thread.interrupt_main", lambda: interrupts.append(1))

    budget = limits.Budget(
        threads=1,
        memory_bytes=1,
        disk_fraction=0.7,
        gpus=0,
        disk_check_interval=0.01,
        machine_cores=1,
        machine_memory_bytes=1,
        machine_gpus=0,
    )
    monitor = limits.DiskMonitor(["/"], budget)
    sampling = threading.Event()
    release = threading.Event()

    def slow_sample():
        """Over the ceiling, but slow enough that the step finishes first."""
        sampling.set()
        release.wait(5)
        return limits.DiskState(path="/", total_bytes=100, used_bytes=99, ceiling_bytes=70)

    monkeypatch.setattr(monitor, "sample", slow_sample)
    monitor.start()
    assert sampling.wait(5), "the watcher never reached sample()"

    # The step finishing, while sample() is still in flight.
    stopped = threading.Thread(target=monitor.stop)
    stopped.start()
    release.set()
    stopped.join(10)

    assert not interrupts, "interrupted the main thread after the step it was watching had finished"
    assert monitor.tripped is None
