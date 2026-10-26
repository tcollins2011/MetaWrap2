"""One resource budget for a whole run: threads, memory, disk, and GPUs.

The budget is resolved once from ``[limits]``, divided among concurrent workers, and reported.
Without it, a cohort driver running N samples at once, each passing the configured thread count to
its own tools, would request N times the machine.

- Fractions are the default, so ``max_thread_fraction = 0.5`` means the same on any machine.
  Absolute values (``max_threads``) are accepted and take precedence, for schedulers that allocate
  a fixed number of cores.
- The machine is measured. In a container ``os.cpu_count()`` and ``SC_PHYS_PAGES`` report the
  host's cores and RAM, which can exceed the cgroup limits, so cgroup limits are read and the
  smaller value is used.
- Reaching the disk ceiling stops the run before the filesystem fills. The step is recorded as
  interrupted and ``--resume`` redoes it.
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

#: Defaults, as a fraction of the machine. Threads are capped at half, since a workstation usually
#: runs other work and oversubscribed cores lose time to contention. Memory is higher because
#: assemblers either fit or fail; a fifth is left for the OS and page cache. Disk stops at 70% so
#: the filesystem stays usable for other work.
DEFAULT_THREAD_FRACTION = 0.5
DEFAULT_MEMORY_FRACTION = 0.8
DEFAULT_DISK_FRACTION = 0.7

#: How often to re-check the filesystem during a long step, in seconds. A step can pass the
#: pre-flight estimate in :mod:`metawrap2.scratch` and still fill the disk hours later.
DEFAULT_DISK_CHECK_INTERVAL = 60

#: The keys ``[limits]`` understands, mirroring ``SETTINGS_KEYS`` in :mod:`metawrap2.config` so a
#: typo is reported with a suggestion.
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
    "gpu_strict": "refuse to start a GPU step on a card something else is already using "
    "(default false, which warns)",
    "disk_check_interval": "seconds between disk checks during a long step (default %d)"
    % DEFAULT_DISK_CHECK_INTERVAL,
}


class LimitsError(Exception):
    """A ``[limits]`` value cannot be used as written."""


class DiskCeilingReached(KeyboardInterrupt):
    """A filesystem passed ``max_disk_fraction`` while a step was running.

    Subclasses KeyboardInterrupt, like :class:`metawrap2.interrupt.Stopped`, so the existing
    ``finally`` blocks, manifest updates and ``--resume`` handling for interrupted runs apply.
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

    cgroup v2 writes ``quota period`` into one file; v1 splits it across two. Both are checked.
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
    """Cores this process may run on.

    ``len(os.sched_getaffinity(0))`` is used over ``os.cpu_count()`` because it respects CPU
    pinning by a scheduler. The cgroup quota is applied on top, since a container can be
    affinity-unrestricted and still throttled.
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
        # cgroup v1 writes a sentinel (~2^63) instead of "max" for unlimited, so only a limit
        # below physical RAM is used.
        if not total or limit < total:
            total = limit
    return max(1, total)


def machine_gpus() -> int:
    """Number of usable NVIDIA GPUs.

    ``CUDA_VISIBLE_DEVICES`` takes precedence when set. Otherwise the count comes from
    ``nvidia-smi --list-gpus``; without nvidia-smi the count is zero.
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

    Cached because the budget is resolved in several places (the run banner, each module's thread
    resolution) and the count does not change mid-run.
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
    gpu_strict: bool = False
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
        # Zero is valid ("use no GPUs"); only negatives are rejected.
        raise LimitsError("limits.%s must be 0 or more, got %r" % (name, value))
    return number


#: Accepted spellings of a boolean. A config file gives ``true`` and a command-line flag gives the
#: string "true"; both go through the same check.
_TRUE = ("1", "true", "yes", "on")
_FALSE = ("0", "false", "no", "off")


def _boolean(value: Any, name: str) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise LimitsError("limits.%s must be true or false, got %r" % (name, value))


def config_from_table(table: Dict[str, Any]) -> LimitsConfig:
    """Build a :class:`LimitsConfig` from a parsed ``[limits]`` table, validating as we go.

    Values are checked here, at startup, instead of at point of use.
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
    if "gpu_strict" in table:
        config.gpu_strict = _boolean(table["gpu_strict"], "gpu_strict")
    if "disk_check_interval" in table:
        config.disk_check_interval = _positive_int(
            table["disk_check_interval"], "disk_check_interval"
        )
    return config


# ── command-line overrides ───────────────────────────────────────────────────────────────
# Global flags that change the budget for one run without editing the config file. Each module
# loads the config separately, so overrides are stored here and applied by resolve(), which every
# consumer of the budget calls.

_OVERRIDES: Dict[str, Any] = {}

#: How each override flag is validated: the same checks as the config file, so ``--max-threads 0``
#: fails like ``max_threads = 0``.
_OVERRIDE_PARSERS: Dict[str, Callable[[Any, str], Any]] = {
    "max_threads": _positive_int,
    "max_thread_fraction": _fraction,
    "max_memory_gb": _positive_number,
    "max_memory_fraction": _fraction,
    "max_disk_fraction": _fraction,
    "max_gpus": _non_negative_int,
    "gpu_strict": _boolean,
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
    """Drop all overrides. Used by tests so one test's budget does not leak into another."""
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

    Both the ceiling and the machine totals it was derived from are kept for the provenance
    record.
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

        Threads and memory are divided with integer division and a minimum of one, so more
        workers than threads oversubscribes by at most one thread each instead of issuing
        zero-thread commands. GPUs divide the same way with a minimum of zero; with more workers
        than GPUs, some workers get none and must fall back to CPU.
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


#: Pairs where setting one on the command line unsets the other from the config file. Otherwise
#: ``--max-thread-fraction 0.25`` would have no effect when the config sets ``max_threads``,
#: because the absolute value takes precedence.
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

    An absolute setting takes precedence over the matching fraction but is clamped to what the
    machine has (256 GB requested on a 64 GB machine gives 64 GB).

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


def gpu_strict(config: Optional[LimitsConfig] = None) -> bool:
    """Whether a GPU step must refuse a card something else is using. See :func:`gpu_preflight`.

    Not a :class:`Budget` field, because it is not divided among workers. It is read through the
    override table like the rest of ``[limits]``, so a flag takes precedence over the config file.
    """
    return _with_overrides(config or LimitsConfig()).gpu_strict


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

    A step's target directory usually does not exist yet; the filesystem it will be on does.
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
    # total - free, not usage.used: most filesystems reserve space for root, so used + free <
    # total. This counts what a non-root writer can still use.
    return DiskState(
        path=probe,
        total_bytes=usage.total,
        used_bytes=usage.total - usage.free,
        ceiling_bytes=int(usage.total * budget.disk_fraction),
    )


def distinct_filesystems(paths: Iterable[str]) -> List[str]:
    """*paths* reduced to one representative per filesystem, in the order first seen.

    Output and scratch are often on the same device; it is reported once.
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
    """Check whether all of *paths* are below the disk ceiling. Returns ``(ok, message)``.

    A filesystem that cannot be measured counts as below the ceiling, as in
    :func:`metawrap2.scratch.check_space`.
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

    :mod:`metawrap2.scratch` estimates before a step starts. This samples disk usage while the
    step runs, for steps that start with enough room and fill the disk later.

    The stop is delivered as a KeyboardInterrupt in the main thread via
    :func:`_thread.interrupt_main`, which reaches code blocked in a ``subprocess`` wait. The
    interrupt path then records the step as interrupted for ``--resume``.

    Used as a context manager, so the thread does not outlive the step it watches.
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
            # Re-checked under the lock because sample() can take seconds on a stalled NFS mount.
            # If the step finished meanwhile, stop() has been called, and interrupting now would
            # abort the next step under `run` or `cohort` instead.
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

        The flag is set under the lock the watcher holds while interrupting, so a finished step
        cannot receive a KeyboardInterrupt after this returns.
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

    Without it, parallel per-sample GPU binners all run on device 0 and run out of its memory
    while other GPUs are idle.

    :meth:`acquire` returns None instead of blocking when nothing is free; the caller decides
    whether to fall back to CPU.
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

    Pinning uses ``CUDA_VISIBLE_DEVICES`` instead of a tool's device flag, because not every tool
    has one, and the pinned child sees its card as device 0. ``device=None`` sets an empty value,
    which tells CUDA there are no GPUs so frameworks fall back to CPU.
    """
    env = dict(base if base is not None else os.environ)
    env["CUDA_VISIBLE_DEVICES"] = "" if device is None else str(device)
    return env


# ── GPU preflight ───────────────────────────────────────────────────────────────────────


#: Utilisation at or above which a card counts as in use, even when the driver lists no process on
#: it. Inside a container ``--query-compute-apps`` reports pids from other namespaces as
#: ``[Not Found]`` or omits them. The threshold is above 0 because a card that has just finished a
#: step can still report some load for a moment.
GPU_BUSY_PERCENT = 25

#: Fields requested from ``nvidia-smi --query-gpu``, in order. ``uuid`` joins this query to
#: ``--query-compute-apps``, which identifies cards by uuid.
_GPU_QUERY_FIELDS = (
    "index",
    "uuid",
    "name",
    "utilization.gpu",
    "memory.total",
    "memory.free",
)


@dataclass(frozen=True)
class GpuProcess:
    """One process the driver says is holding a card."""

    pid: int
    memory_mb: int


@dataclass(frozen=True)
class GpuState:
    """One card as the driver reports it, at one moment."""

    index: int
    uuid: str
    name: str
    utilisation_percent: Optional[int]
    memory_total_mb: Optional[int]
    memory_free_mb: Optional[int]
    processes: Tuple[GpuProcess, ...] = ()

    @property
    def busy(self) -> bool:
        """True when another process holds the card or its utilisation is high."""
        if self.processes:
            return True
        return (self.utilisation_percent or 0) >= GPU_BUSY_PERCENT

    def describe(self) -> str:
        """One line: the card's name, utilisation, free memory and other processes."""
        utilisation = (
            "utilisation unknown"
            if self.utilisation_percent is None
            else "utilisation %d%%" % self.utilisation_percent
        )
        if self.memory_free_mb is None or self.memory_total_mb is None:
            memory = "free memory unknown"
        else:
            memory = "%s of %s free" % (
                human(self.memory_free_mb * 1024**2),
                human(self.memory_total_mb * 1024**2),
            )
        if not self.processes:
            held = "no other process on it"
        else:
            held = "%d other process%s on it: %s" % (
                len(self.processes),
                "" if len(self.processes) == 1 else "es",
                ", ".join(
                    "pid %d holding %s" % (p.pid, human(p.memory_mb * 1024**2))
                    for p in self.processes
                ),
            )
        return "%s (%s), %s, %s, %s" % (
            cuda_name(self.index),
            self.name or "unnamed device",
            utilisation,
            memory,
            held,
        )

    def brief(self) -> str:
        """``cuda:3 at 94%``, for listing several cards in one sentence."""
        return "%s at %s" % (
            cuda_name(self.index),
            (
                "unknown load"
                if self.utilisation_percent is None
                else "%d%%" % self.utilisation_percent
            ),
        )


@dataclass(frozen=True)
class GpuProbe:
    """The driver's answer, or why it could not be queried.

    *problem* is a phrase for the message, such as no ``nvidia-smi`` on PATH or a driver that does
    not respond. Neither stops the step; the tool then reports its own device error.
    """

    states: Tuple[GpuState, ...] = ()
    problem: Optional[str] = None


def cuda_name(index: int) -> str:
    """``cuda:N``, the same form as ``--device``."""
    return "cuda:%d" % index


def _as_int(text: str) -> Optional[int]:
    """The integer in a ``--query-gpu`` cell, or None when the driver could not fill it.

    The driver writes ``[N/A]`` for a field it does not have. That is reported as unknown, since 0%
    would say the card is idle.
    """
    digits = "".join(character for character in text if character.isdigit())
    return int(digits) if digits else None


def parse_gpu_query(text: str) -> List[GpuState]:
    """``nvidia-smi --query-gpu`` CSV into one :class:`GpuState` per card, in index order."""
    states: List[GpuState] = []
    for line in text.splitlines():
        cells = [cell.strip() for cell in line.split(",")]
        if len(cells) < len(_GPU_QUERY_FIELDS) or not cells[0].isdigit():
            continue
        states.append(
            GpuState(
                index=int(cells[0]),
                uuid=cells[1],
                name=cells[2],
                utilisation_percent=_as_int(cells[3]),
                memory_total_mb=_as_int(cells[4]),
                memory_free_mb=_as_int(cells[5]),
            )
        )
    return sorted(states, key=lambda state: state.index)


def query_gpus() -> GpuProbe:
    """The current state of every card, or why the driver could not be queried.

    Makes two ``nvidia-smi`` calls, one for the cards and one for the processes on them; together
    they take a fraction of a second. Only GPU steps call this.

    The process list comes from :func:`metawrap2.usage.gpu_compute_apps` and is joined on the card
    uuid, which is how the driver identifies a card for a compute app.
    """
    executable = shutil.which("nvidia-smi")
    if not executable:
        return GpuProbe(problem="no nvidia-smi on PATH")
    try:
        completed = subprocess.run(
            [
                executable,
                "--query-gpu=%s" % ",".join(_GPU_QUERY_FIELDS),
                "--format=csv,noheader,nounits",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return GpuProbe(problem="nvidia-smi did not answer")
    if completed.returncode != 0:
        return GpuProbe(problem="nvidia-smi exited %d" % completed.returncode)
    states = parse_gpu_query(completed.stdout or "")
    if not states:
        return GpuProbe(problem="nvidia-smi listed no card it could describe")

    # Imported here because only this function needs it, and limits is imported by config, which
    # every module and utility imports.
    from . import usage as _usage

    holders: Dict[str, List[GpuProcess]] = {}
    for pid, uuid, memory_mb in _usage.gpu_compute_apps():
        holders.setdefault(uuid, []).append(GpuProcess(pid=pid, memory_mb=memory_mb))
    return GpuProbe(
        states=tuple(
            replace(state, processes=tuple(holders.get(state.uuid, ()))) for state in states
        )
    )


#: The cost of sharing a card, added to every preflight message about a busy card.
_CONTENTION_COST = (
    "On a shared card the step's wall time inflates, and its GPU usage record covers both jobs, so "
    "its own cost cannot be separated out."
)

#: How to turn the warning into a refusal: a flag for one command, or a config key for every run.
_STRICT_HINT = (
    "`--gpu-strict` (or `[limits] gpu_strict = true`) refuses to start instead; use it for "
    "benchmarks."
)


@dataclass(frozen=True)
class GpuPreflight:
    """The text to print before a GPU step, and whether the step may start.

    *message* is empty when nothing was checked (a CPU run); callers print it only when it is set.
    """

    device: str
    message: str = ""
    ok: bool = True
    contended: bool = False


def _card_count(count: int) -> str:
    return "1 card" if count == 1 else "%d cards" % count


def _card_inventory(states: Tuple[GpuState, ...]) -> str:
    """The range of card names on this host, for the message about a missing card."""
    if not states:
        return "none at all"
    if len(states) == 1:
        return "%s only" % cuda_name(states[0].index)
    return "%s to %s" % (cuda_name(states[0].index), cuda_name(states[-1].index))


def gpu_preflight(device: str, strict: bool = False) -> GpuPreflight:
    """Check the card a step is about to train on, before the step starts.

    :func:`gpu_env` pins the step to the card ``--device cuda:N`` names but does not check it. GPU
    steps can run for hours. With an index above the device count the pin hides every card, so the
    binner fails at its first allocation or trains on the CPU. A card already in use makes the
    step slower and its GPU usage unmeasurable.

    A missing index is refused. A card in use gives a warning, or a refusal when *strict* is set.

    Indices are checked against the driver's list, which is also what the pin selects. On a host
    with mixed cards ``nvidia-smi`` orders by PCI bus id while CUDA's default ``CUDA_DEVICE_ORDER``
    is ``FASTEST_FIRST``, so the count is reliable but the name reported for ``cuda:1`` may not be
    the card the tool uses.
    """
    requested = (device or "").strip()
    if not requested.startswith("cuda"):
        # --device cpu, or --device auto resolved to the CPU: nothing is probed or printed.
        return GpuPreflight(device=requested)
    index: Optional[int] = None
    if requested.startswith("cuda:"):
        tail = requested.split(":", 1)[1].strip()
        if not tail.isdigit():
            # gpu_pin passes a malformed index to the tool unchanged, so it is not checked here.
            return GpuPreflight(device=requested)
        index = int(tail)

    probe = query_gpus()
    if probe.problem is not None:
        return GpuPreflight(
            device=requested,
            message="GPU preflight: %s, so %s cannot be checked before the step. Starting anyway - "
            "the tool reports its own device error if there is one." % (probe.problem, requested),
        )

    if index is not None:
        matching = [state for state in probe.states if state.index == index]
        if not matching:
            return GpuPreflight(
                device=requested,
                message="--device %s names a card this host does not have: the driver reports %s, "
                "%s. That is not a slower run, it is a run with no GPU in it - the pin hides every "
                "card, so the binner either fails at its first allocation or trains on the CPU for "
                "hours. Pick one of the cards above, or --device cpu."
                % (requested, _card_count(len(probe.states)), _card_inventory(probe.states)),
                ok=False,
            )
        state = matching[0]
        line = "GPU preflight: %s" % state.describe()
        if not state.busy:
            return GpuPreflight(device=requested, message=line)
        verdict = (
            "Refusing to start: --gpu-strict is set. Name a free card with --device cuda:N, wait "
            "for this one, or drop --gpu-strict to share it."
            if strict
            else "Starting anyway, since sharing a card is sometimes what you want. " + _STRICT_HINT
        )
        return GpuPreflight(
            device=requested,
            message="%s\n%s is in use. %s %s" % (line, requested, _CONTENTION_COST, verdict),
            ok=not strict,
            contended=True,
        )

    # Bare `cuda` means any GPU. gpu_pin leaves it unpinned, so the binner's framework picks the
    # card. Unlike GpuPool.acquire() there is no CPU fallback, so when every card is busy the step
    # shares one, and the message says so.
    free = [state for state in probe.states if not state.busy]
    inventory = ", ".join(state.brief() for state in probe.states)
    if free:
        # List only which cards are free; per-card utilisation makes the line long on large hosts.
        found = (
            "all %s free" % _card_count(len(probe.states))
            if len(free) == len(probe.states)
            else "%d of %s free (%s)"
            % (
                len(free),
                _card_count(len(probe.states)),
                ", ".join(cuda_name(state.index) for state in free),
            )
        )
        return GpuPreflight(
            device=requested,
            message="GPU preflight: %s. `--device cuda` pins nothing, so the tool's own framework "
            "picks the device - usually cuda:0, free or not. Name a card with --device cuda:N to be "
            "sure of a free one." % found,
        )
    verdict = (
        "Refusing to start: --gpu-strict is set. Wait for a card, or --device cpu."
        if strict
        else "Starting anyway on a shared card. " + _STRICT_HINT
    )
    return GpuPreflight(
        device=requested,
        message="GPU preflight: all %s in use (%s)\n`--device cuda` pins nothing and does not fall "
        "back to the CPU: the framework takes a card anyway - in practice cuda:0 - and shares it. "
        "%s %s" % (_card_count(len(probe.states)), inventory, _CONTENTION_COST, verdict),
        ok=not strict,
        contended=True,
    )
