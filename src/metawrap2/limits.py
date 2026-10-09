"""One resource budget for a whole run: threads, memory, disk, and GPUs.

Every knob for "how much of this machine may MetaWrap2 use" used to be local and independent.
``-t`` set threads for one module, ``[settings] threads`` set a default, :mod:`metawrap2.scratch`
estimated disk before a step, and memory was whatever each tool decided for itself. Nothing tied
them together, so nothing could answer the question a user on a shared machine actually asks:
*what fraction of this box will the run take?*

That gap has a sharp edge. A cohort driver running eight samples at once, each passing the
configured thread count down to its own tools, asks for eight times the machine. The individual
numbers all look reasonable; the product does not. This module exists so the budget is resolved
**once**, divided among concurrent workers, and reported.

Three design choices worth stating:

**Fractions, not absolutes, are the default.** ``max_thread_fraction = 0.5`` keeps meaning the
same thing when the run moves to a different machine, while ``max_threads = 96`` quietly becomes
"all of it" on a smaller one. Absolutes are still accepted and win when given, because a
scheduler that handed you exactly 16 cores meant 16, not half of what it can see.

**The machine is measured, not assumed.** Under a container, ``os.cpu_count()`` reports the
host's cores and ``SC_PHYS_PAGES`` the host's RAM, both of which can be far more than the cgroup
allows. Overshooting memory there does not degrade, it OOM-kills. So cgroup limits are read and
the smaller number always wins.

**Hitting the disk ceiling stops the run rather than filling the filesystem.** A stop is
recoverable: the step is recorded as interrupted and ``--resume`` redoes it. A full filesystem is
not, and on shared infrastructure it is everyone's problem, not just this run's.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import threading
from dataclasses import dataclass, replace
from functools import lru_cache
from types import TracebackType
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Type

#: Defaults, as a fraction of the machine. Threads are capped at half because a metagenomic run
#: is rarely the only thing on a workstation, and oversubscribing cores costs more in contention
#: than it gains. Memory goes higher because assemblers either fit or fail - there is no graceful
#: degradation - so leaving a fifth for the OS and page cache is enough. Disk sits between:
#: enough room to work, while leaving a filesystem other things can still use.
DEFAULT_THREAD_FRACTION = 0.5
DEFAULT_MEMORY_FRACTION = 0.8
DEFAULT_DISK_FRACTION = 0.7

#: How often to re-check the filesystem during a long step, in seconds. A step that writes for
#: hours can pass a pre-flight estimate and still fill the disk, so the estimate in
#: :mod:`metawrap2.scratch` is necessary but not sufficient.
DEFAULT_DISK_CHECK_INTERVAL = 60

#: The keys ``[limits]`` understands, mirroring ``SETTINGS_KEYS`` in :mod:`metawrap2.config` so a
#: typo is reported with a suggestion rather than silently ignored.
LIMITS_KEYS: Dict[str, str] = {
    "max_threads": "absolute thread ceiling for the whole run (overrides max_thread_fraction)",
    "max_thread_fraction": "fraction of the machine's cores to use (default %.2f)"
    % DEFAULT_THREAD_FRACTION,
    "max_memory_gb": "absolute memory ceiling in GB (overrides max_memory_fraction)",
    "max_memory_fraction": "fraction of the machine's RAM to use (default %.2f)"
    % DEFAULT_MEMORY_FRACTION,
    "max_disk_fraction": "stop when a filesystem is this full (default %.2f)"
    % DEFAULT_DISK_FRACTION,
    "max_gpus": "how many GPUs the run may use (default: all visible)",
    "disk_check_interval": "seconds between disk checks during a long step (default %d)"
    % DEFAULT_DISK_CHECK_INTERVAL,
}


class LimitsError(Exception):
    """A ``[limits]`` value cannot be used as written."""


class DiskCeilingReached(KeyboardInterrupt):
    """A filesystem passed ``max_disk_fraction`` while a step was running.

    Subclasses KeyboardInterrupt for the same reason :class:`metawrap2.interrupt.Stopped` does:
    every ``finally`` block, manifest update and ``--resume`` path already written for an
    interrupted run then applies unchanged, so stopping cleanly costs no new machinery.
    """


# ── measuring the machine ────────────────────────────────────────────────────────────────


def _read_int(path: str) -> Optional[int]:
    """First integer in *path*, or None if unreadable or the cgroup wrote ``max``."""
    try:
        with open(path) as handle:
            text = handle.read().strip()
    except OSError:
        return None
    if not text or text.startswith("max"):
        return None
    match = re.match(r"-?\d+", text)
    return int(match.group(0)) if match else None


def _cgroup_cpu_quota() -> Optional[int]:
    """Cores the cgroup allows, rounded up, or None if unlimited or unreadable.

    cgroup v2 writes ``quota period`` into one file; v1 splits it across two. Both are checked
    because container runtimes still ship either.
    """
    try:
        with open("/sys/fs/cgroup/cpu.max") as handle:
            parts = handle.read().split()
        if len(parts) == 2 and not parts[0].startswith("max"):
            quota, period = int(parts[0]), int(parts[1])
            if quota > 0 and period > 0:
                return max(1, -(-quota // period))
    except (OSError, ValueError):
        pass
    quota_us = _read_int("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
    period_us = _read_int("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
    if quota_us and quota_us > 0 and period_us and period_us > 0:
        return max(1, -(-quota_us // period_us))
    return None


def machine_cores() -> int:
    """Cores this process may actually run on.

    ``len(os.sched_getaffinity(0))`` beats ``os.cpu_count()`` because it respects CPU pinning - a
    scheduler that gave us 8 of 128 cores should not see 128. The cgroup quota is applied on top,
    since a container can be affinity-unrestricted and still throttled.
    """
    affinity = getattr(os, "sched_getaffinity", None)  # Linux only
    try:
        count = len(affinity(0)) if affinity else (os.cpu_count() or 1)
    except OSError:
        count = os.cpu_count() or 1
    quota = _cgroup_cpu_quota()
    if quota:
        count = min(count, quota)
    return max(1, count)


def machine_memory_bytes() -> int:
    """Total memory this process may use, honouring a cgroup limit if one is set."""
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        total = 0
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        limit = _read_int(path)
        if not limit or limit <= 0:
            continue
        # cgroup v1 writes an absurd sentinel (~2^63) rather than "max" for "unlimited", and a
        # limit above physical RAM is meaningless either way, so only a smaller one counts.
        if not total or limit < total:
            total = limit
    return max(1, total)


def machine_gpus() -> int:
    """Number of usable NVIDIA GPUs.

    ``CUDA_VISIBLE_DEVICES`` is honoured first: if something already narrowed the set, that is
    the answer regardless of what the driver reports. Otherwise ``nvidia-smi --list-gpus`` is the
    only reliable count, and its absence means zero rather than an error - a CPU-only machine is
    a normal machine.
    """
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible is not None:
        stripped = visible.strip()
        if not stripped or stripped == "-1":
            return 0
        return len([token for token in stripped.split(",") if token.strip()])
    return _count_gpus()


@lru_cache(maxsize=1)
def _count_gpus() -> int:
    """``nvidia-smi --list-gpus``, asked once per process.

    Cached because the budget is resolved from more than one place (the run banner, each module's
    thread resolution) and the answer cannot change mid-run, so paying for a subprocess each time
    would be waste.
    """
    executable = shutil.which("nvidia-smi")
    if not executable:
        return 0
    try:
        completed = subprocess.run(
            [executable, "--list-gpus"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return 0
    lines = (completed.stdout or "").splitlines()
    return len([line for line in lines if line.strip().startswith("GPU")])


# ── the configured request ───────────────────────────────────────────────────────────────


@dataclass
class LimitsConfig:
    """What the user asked for in ``[limits]``, before it is measured against the machine."""

    max_threads: Optional[int] = None
    max_thread_fraction: float = DEFAULT_THREAD_FRACTION
    max_memory_gb: Optional[float] = None
    max_memory_fraction: float = DEFAULT_MEMORY_FRACTION
    max_disk_fraction: float = DEFAULT_DISK_FRACTION
    max_gpus: Optional[int] = None
    disk_check_interval: int = DEFAULT_DISK_CHECK_INTERVAL


def _fraction(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise LimitsError("limits.%s must be a number between 0 and 1, got %r" % (name, value))
    if not 0 < number <= 1:
        raise LimitsError("limits.%s must be greater than 0 and at most 1, got %r" % (name, value))
    return number


def _positive_int(value: Any, name: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise LimitsError("limits.%s must be a whole number, got %r" % (name, value))
    if number < 1:
        raise LimitsError("limits.%s must be at least 1, got %r" % (name, value))
    return number


def _positive_number(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise LimitsError("limits.%s must be a number, got %r" % (name, value))
    if number <= 0:
        raise LimitsError("limits.%s must be greater than 0, got %r" % (name, value))
    return number


def _non_negative_int(value: Any, name: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise LimitsError("limits.%s must be a whole number, got %r" % (name, value))
    if number < 0:
        # Zero is meaningful - "this run may not touch the GPUs" - so only negatives are bad.
        raise LimitsError("limits.%s must be 0 or more, got %r" % (name, value))
    return number


def config_from_table(table: Dict[str, Any]) -> LimitsConfig:
    """Build a :class:`LimitsConfig` from a parsed ``[limits]`` table, validating as we go.

    Values are checked here rather than at point of use so a bad number stops the run in its
    first second instead of its first hour.
    """
    config = LimitsConfig()
    if "max_threads" in table:
        config.max_threads = _positive_int(table["max_threads"], "max_threads")
    if "max_thread_fraction" in table:
        config.max_thread_fraction = _fraction(table["max_thread_fraction"], "max_thread_fraction")
    if "max_memory_gb" in table:
        config.max_memory_gb = _positive_number(table["max_memory_gb"], "max_memory_gb")
    if "max_memory_fraction" in table:
        config.max_memory_fraction = _fraction(table["max_memory_fraction"], "max_memory_fraction")
    if "max_disk_fraction" in table:
        config.max_disk_fraction = _fraction(table["max_disk_fraction"], "max_disk_fraction")
    if "max_gpus" in table:
        config.max_gpus = _non_negative_int(table["max_gpus"], "max_gpus")
    if "disk_check_interval" in table:
        config.disk_check_interval = _positive_int(
            table["disk_check_interval"], "disk_check_interval"
        )
    return config


# ── command-line overrides ───────────────────────────────────────────────────────────────
# The budget is read from the config file, but a user on a busy machine needs to shrink it for
# one run without editing that file. Those flags are global (they apply to every module and
# utility), and the config is loaded separately inside each of them, so the override is kept here
# and applied by resolve() - the one function every consumer of the budget already goes through.
# The alternative, threading four extra arguments through every module's parser, would put the
# same value in ten places and let them disagree.

_OVERRIDES: Dict[str, Any] = {}

#: How each override flag is validated. Same checks as the config file: ``--max-threads 0`` should
#: fail the same way ``max_threads = 0`` does.
_OVERRIDE_PARSERS: Dict[str, Callable[[Any, str], Any]] = {
    "max_threads": _positive_int,
    "max_thread_fraction": _fraction,
    "max_memory_gb": _positive_number,
    "max_memory_fraction": _fraction,
    "max_disk_fraction": _fraction,
    "max_gpus": _non_negative_int,
    "disk_check_interval": _positive_int,
}


def set_override(key: str, value: Any) -> None:
    """Override one ``[limits]`` field for this process, from a command-line flag."""
    parser = _OVERRIDE_PARSERS.get(key)
    if parser is None:
        raise LimitsError("%r is not an overridable limit" % key)
    _OVERRIDES[key] = parser(value, key)


def overrides() -> Dict[str, Any]:
    """The overrides in effect, for the provenance record."""
    return dict(_OVERRIDES)


def clear_overrides() -> None:
    """Drop all overrides. For tests, which must not leak a budget into each other."""
    _OVERRIDES.clear()


# ── the resolved budget ──────────────────────────────────────────────────────────────────


def human(size: Optional[float]) -> str:
    if size is None:
        return "unknown"
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if value < 1024 or unit == "PB":
            return "%.1f %s" % (value, unit)
        value /= 1024
    return "%.1f PB" % value


@dataclass(frozen=True)
class Budget:
    """The resolved budget for a run: absolute numbers, decided once.

    Both the ceiling and what it was measured against are kept, because "96 threads" means
    nothing in a provenance record without "of 128 cores" beside it.
    """

    threads: int
    memory_bytes: int
    disk_fraction: float
    gpus: int
    disk_check_interval: int
    machine_cores: int
    machine_memory_bytes: int
    machine_gpus: int

    @property
    def memory_gb(self) -> float:
        return self.memory_bytes / 1024**3

    def per_worker(self, concurrency: int) -> Budget:
        """The slice of this budget one of *concurrency* concurrent workers may use.

        Divides; never multiplies. Integer division floors, with a floor of one thread, so asking
        for more workers than there are threads oversubscribes by at most one thread each rather
        than issuing zero-thread commands. GPUs divide the same way but floor at zero, because
        there is no such thing as a fractional device: with 8 GPUs and 20 workers some workers
        get none and must fall back to CPU.
        """
        if concurrency < 1:
            raise LimitsError("concurrency must be at least 1, got %r" % concurrency)
        return Budget(
            threads=max(1, self.threads // concurrency),
            memory_bytes=max(1, self.memory_bytes // concurrency),
            disk_fraction=self.disk_fraction,
            gpus=self.gpus // concurrency,
            disk_check_interval=self.disk_check_interval,
            machine_cores=self.machine_cores,
            machine_memory_bytes=self.machine_memory_bytes,
            machine_gpus=self.machine_gpus,
        )

    def describe(self) -> str:
        """One line for the run header: the budget, and what it is a share of."""
        parts = [
            "%d/%d threads" % (self.threads, self.machine_cores),
            "%s/%s memory" % (human(self.memory_bytes), human(self.machine_memory_bytes)),
            "disk ceiling %d%%" % round(self.disk_fraction * 100),
        ]
        if self.machine_gpus:
            parts.append("%d/%d GPUs" % (self.gpus, self.machine_gpus))
        return ", ".join(parts)

    def as_dict(self) -> Dict[str, Any]:
        """Flat mapping for the provenance record."""
        return {
            "threads": self.threads,
            "memory_bytes": self.memory_bytes,
            "disk_fraction": self.disk_fraction,
            "gpus": self.gpus,
            "disk_check_interval": self.disk_check_interval,
            "machine_cores": self.machine_cores,
            "machine_memory_bytes": self.machine_memory_bytes,
            "machine_gpus": self.machine_gpus,
        }


#: Pairs where setting one on the command line must unset the other from the config file. Without
#: this, ``--max-thread-fraction 0.25`` would be ignored whenever the config file happened to set
#: an absolute ``max_threads``, because the absolute always wins - and the flag would silently do
#: nothing, which is the failure mode this whole module exists to prevent.
_EXCLUSIVE_PAIRS = (
    ("max_threads", "max_thread_fraction"),
    ("max_memory_gb", "max_memory_fraction"),
)


def _with_overrides(config: LimitsConfig) -> LimitsConfig:
    """*config* with any command-line overrides applied on top."""
    if not _OVERRIDES:
        return config
    updates = dict(_OVERRIDES)
    for absolute, fraction in _EXCLUSIVE_PAIRS:
        if fraction in updates and absolute not in updates:
            updates[absolute] = None
    return replace(config, **updates)


def resolve(config: Optional[LimitsConfig] = None) -> Budget:
    """Turn a :class:`LimitsConfig` into an absolute :class:`Budget` for this machine.

    An absolute setting wins over the matching fraction, but is still clamped to what the machine
    has: asking for 256 GB on a 64 GB box is a mistake worth correcting quietly, not a request
    worth honouring into an OOM kill.

    Command-line overrides are applied last, so a flag beats the config file.
    """
    settings = _with_overrides(config or LimitsConfig())
    cores = machine_cores()
    memory = machine_memory_bytes()
    gpus = machine_gpus()

    if settings.max_threads is not None:
        threads = min(settings.max_threads, cores)
    else:
        threads = max(1, int(cores * settings.max_thread_fraction))

    if settings.max_memory_gb is not None:
        memory_bytes = min(int(settings.max_memory_gb * 1024**3), memory)
    else:
        memory_bytes = max(1, int(memory * settings.max_memory_fraction))

    allowed_gpus = gpus if settings.max_gpus is None else min(settings.max_gpus, gpus)

    return Budget(
        threads=threads,
        memory_bytes=memory_bytes,
        disk_fraction=settings.max_disk_fraction,
        gpus=allowed_gpus,
        disk_check_interval=settings.disk_check_interval,
        machine_cores=cores,
        machine_memory_bytes=memory,
        machine_gpus=gpus,
    )


# ── the disk ceiling ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class DiskState:
    """One filesystem's usage at one moment, relative to the configured ceiling."""

    path: str
    total_bytes: int
    used_bytes: int
    ceiling_bytes: int

    @property
    def fraction_used(self) -> float:
        return self.used_bytes / self.total_bytes if self.total_bytes else 0.0

    @property
    def over(self) -> bool:
        return self.used_bytes >= self.ceiling_bytes

    @property
    def headroom_bytes(self) -> int:
        """Bytes that may still be written before the ceiling. Never negative."""
        return max(0, self.ceiling_bytes - self.used_bytes)

    def describe(self) -> str:
        ceiling_percent = self.ceiling_bytes / self.total_bytes * 100 if self.total_bytes else 0.0
        return "%s is %.0f%% full (ceiling %.0f%%, %s of headroom)" % (
            self.path,
            self.fraction_used * 100,
            ceiling_percent,
            human(self.headroom_bytes),
        )


def _nearest_existing(path: str) -> Optional[str]:
    """*path* or its nearest existing ancestor.

    A step's target directory usually does not exist yet, but the filesystem it will live on
    does, and that is what needs measuring.
    """
    probe = os.path.abspath(path)
    while probe and not os.path.exists(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            return None
        probe = parent
    return probe or None


def disk_state(path: str, budget: Budget) -> Optional[DiskState]:
    """Current usage of the filesystem holding *path*, or None if it cannot be measured."""
    probe = _nearest_existing(path)
    if probe is None:
        return None
    try:
        usage = shutil.disk_usage(probe)
    except OSError:
        return None
    # total - free, not usage.used: most filesystems reserve some space for root, so
    # used + free < total. What matters is what a non-root writer can still consume.
    return DiskState(
        path=probe,
        total_bytes=usage.total,
        used_bytes=usage.total - usage.free,
        ceiling_bytes=int(usage.total * budget.disk_fraction),
    )


def distinct_filesystems(paths: Iterable[str]) -> List[str]:
    """*paths* reduced to one representative per filesystem, in the order first seen.

    Output and scratch are often the same device, and reporting one filesystem twice in a budget
    header is noise.
    """
    seen: Dict[int, str] = {}
    ordered: List[str] = []
    for path in paths:
        if not path:
            continue
        probe = _nearest_existing(path)
        if probe is None:
            continue
        try:
            device = os.stat(probe).st_dev
        except OSError:
            continue
        if device not in seen:
            seen[device] = probe
            ordered.append(probe)
    return ordered


def check_ceiling(paths: Iterable[str], budget: Budget) -> Tuple[bool, str]:
    """Are all of *paths* below the disk ceiling? Returns ``(ok, message)``.

    An unmeasurable filesystem counts as fine, matching :func:`metawrap2.scratch.check_space`:
    refusing to run because a stat failed would be worse than running without the guard.
    """
    states = [disk_state(path, budget) for path in distinct_filesystems(paths)]
    measured = [state for state in states if state is not None]
    if not measured:
        return True, "no filesystem could be measured"
    over = [state for state in measured if state.over]
    if over:
        return False, "; ".join(state.describe() for state in over)
    return True, "; ".join(state.describe() for state in measured)


class DiskMonitor:
    """Watches filesystems during a long step and stops the run at the ceiling.

    :mod:`metawrap2.scratch` estimates before a step starts, which catches the case where there
    was never enough room. It cannot catch the case that happens more often: a step that starts
    with room, writes for six hours, and fills the disk at hour five. This samples while the step
    runs.

    The stop is delivered as a KeyboardInterrupt in the main thread, via
    :func:`_thread.interrupt_main` - the one mechanism that reaches code blocked in a
    ``subprocess`` wait. The existing interrupt path then handles it: the step is recorded as
    interrupted and ``--resume`` picks it up. Raising in the monitor thread would achieve nothing,
    since nothing reads that stack.

    Used as a context manager, so the thread cannot outlive the step it is watching.
    """

    def __init__(
        self,
        paths: Iterable[str],
        budget: Budget,
        on_stop: Optional[Callable[[DiskState], None]] = None,
    ):
        self.paths = [path for path in paths if path]
        self.budget = budget
        self.on_stop = on_stop
        #: Set when the ceiling was hit, so a caller can tell a disk stop from a Ctrl-C.
        self.tripped: Optional[DiskState] = None
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        #: Held while deciding to interrupt, and by stop(), so the two cannot interleave.
        self._lock = threading.Lock()

    def sample(self) -> Optional[DiskState]:
        """The first filesystem found over the ceiling, or None if all are under it."""
        for path in distinct_filesystems(self.paths):
            state = disk_state(path, self.budget)
            if state is not None and state.over:
                return state
        return None

    def _watch(self) -> None:
        import _thread

        while not self._stop.wait(self.budget.disk_check_interval):
            state = self.sample()
            if state is None:
                continue
            # Checked again, after sampling and under the lock, because sample() takes time - two
            # statvfs calls and an os.stat, which on a stalled NFS mount is seconds. A step that
            # finished during that window has already called stop(), and interrupting now would
            # raise KeyboardInterrupt into whatever the process does next: under `run` or `cohort`
            # that is the following step, which would be aborted and recorded as interrupted while
            # the step actually being watched is reported as a disk stop it never hit.
            with self._lock:
                if self._stop.is_set():
                    return
                self.tripped = state
                if self.on_stop is not None:
                    try:
                        self.on_stop(state)
                    except Exception:  # noqa: BLE001,S110 - reporting must not swallow the stop
                        pass
                _thread.interrupt_main()
                return

    def start(self) -> None:
        """Begin watching. A no-op when there is nothing to watch."""
        if self._thread is not None or not self.paths or self.budget.disk_check_interval <= 0:
            return
        self._thread = threading.Thread(
            target=self._watch, name="metawrap2-disk-monitor", daemon=True
        )
        self._thread.start()

    def stop(self) -> None:
        """Stop watching and wait for the thread to notice. Safe to call more than once.

        The flag is set under the same lock the watcher takes before interrupting, so this either
        happens before the watcher has decided - and it then returns without interrupting - or
        after it has already interrupted. There is no window where a step that has finished still
        gets a KeyboardInterrupt delivered on its behalf.
        """
        with self._lock:
            self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5)
            self._thread = None

    # typing.Self is 3.11+, and `from __future__ import annotations` makes this string-evaluated.
    def __enter__(self) -> DiskMonitor:  # noqa: PYI034
        self.start()
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc: Optional[BaseException],
        tb: Optional[TracebackType],
    ) -> None:
        self.stop()


# ── GPU arbitration ─────────────────────────────────────────────────────────────────────


class GpuPool:
    """Hands out GPU device ids so two workers never land on the same card.

    Without this, running one GPU binner per sample in parallel puts every process on device 0,
    where they OOM each other out of a single card's memory while the machine's other GPUs sit
    idle. The pool is what lets a per-sample GPU stage actually use all of them.

    :meth:`acquire` returns None rather than blocking when nothing is free. Waiting for a device
    that may be busy for an hour is worse than running that one sample on CPU, and only the
    caller knows whether a CPU fallback exists.
    """

    def __init__(self, devices: Optional[Iterable[int]] = None, count: Optional[int] = None):
        if devices is not None:
            self._free = list(dict.fromkeys(int(device) for device in devices))
        else:
            self._free = list(range(count if count is not None else machine_gpus()))
        self._lock = threading.Lock()

    @classmethod
    def from_budget(cls, budget: Budget) -> GpuPool:
        return cls(count=budget.gpus)

    def __len__(self) -> int:
        """How many devices are free right now."""
        with self._lock:
            return len(self._free)

    def acquire(self) -> Optional[int]:
        """Take a device id, or None if every device is busy."""
        with self._lock:
            if not self._free:
                return None
            return self._free.pop(0)

    def release(self, device: Optional[int]) -> None:
        """Return a device to the pool. Releasing None is a no-op, so callers need no branch."""
        if device is None:
            return
        with self._lock:
            if device not in self._free:
                self._free.append(device)


def gpu_env(device: Optional[int], base: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Environment for a child process pinned to *device*.

    Pinning uses ``CUDA_VISIBLE_DEVICES`` rather than a tool's own device flag, because not every
    tool has one and because a pinned child sees its card as device 0 - so nothing downstream
    needs to know which physical GPU it got. ``device=None`` means "no GPU": an empty value is
    how CUDA is told there are none, which turns a framework's device probe into a clean CPU
    fallback instead of a crash.
    """
    env = dict(base if base is not None else os.environ)
    env["CUDA_VISIBLE_DEVICES"] = "" if device is None else str(device)
    return env
