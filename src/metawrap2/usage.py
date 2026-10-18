"""Measuring what a step actually cost, rather than what it was told it could use.

Memory budgets in MetaWrap2 are guesses. ``-m 64`` is passed to metaSPAdes, the space multipliers
in :mod:`metawrap2.scratch` are multipliers someone chose, and nothing has ever measured whether
any of it is close. Meanwhile every run already records its inputs, its duration and its outputs -
so the one missing number is what the step's peak memory really was.

:func:`run_measured` gets it exactly, per step, by reaping the child with :func:`os.wait4` instead
of :func:`subprocess.call`. That returns the rusage of *that* process, which
``getrusage(RUSAGE_CHILDREN)`` cannot do - it reports a running maximum across every child the
process has ever waited for, so with several steps it tells you only which one was worst.

The measurement covers the step's whole process tree, because ``ru_maxrss`` of a waited-for child
includes its own reaped descendants - which is what we want, since the thing being measured is
``mamba run -n env metaspades.py``, not the shim.

Where ``os.wait4`` does not exist (Windows), this degrades to running the command with no
measurement rather than refusing to run it.

GPU memory is the one cost ``rusage`` cannot see, and it is the one that decides whether a step
runs at all: a card either holds the model or the step dies with an allocator error. So when a
driver is present, a sampler thread records the peak VRAM the step's own process tree held on a
single device, which is the number someone choosing a card needs. Three tools here can use a GPU -
SemiBin2, VAMB and COMEBin, all in the ``learned`` environment - and every figure in the README's
hardware table comes from this measurement rather than from an upstream README.
"""

from __future__ import annotations

import os
import platform
import resource
import shutil
import subprocess
import threading
import time
from typing import IO, Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

#: True when per-process resource usage can be measured here.
CAN_MEASURE = hasattr(os, "wait4")

#: How often the GPU sampler asks the driver. SemiBin2 and COMEBin allocate and free within
#: seconds, and a 30 s interval missed a training peak entirely on a short bin set; a call to
#: nvidia-smi costs ~50 ms, so 5 s is cheap enough to leave on for every step.
GPU_SAMPLE_SECONDS = 5.0

MIB = 1024 * 1024

#: Cached answer to "is there a driver here", so a CPU-only machine does not pay a PATH lookup per
#: step and a GPU machine does not pay one per sample.
_HAS_GPU: Optional[bool] = None


def maxrss_bytes(raw: int) -> int:
    """Convert ``ru_maxrss`` to bytes.

    Linux reports kilobytes; macOS reports bytes. Getting this wrong is a factor of 1024, which
    would make every recorded number useless, so it is decided by platform rather than guessed
    from magnitude.
    """
    if platform.system() == "Darwin":
        return int(raw)
    return int(raw) * 1024


def human(size: Optional[float]) -> str:
    if size is None:
        return "unknown"
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return "%.1f %s" % (value, unit)
        value /= 1024
    return "%.1f TB" % value


def _exit_code(status: int) -> int:
    """Turn a wait() status into the returncode subprocess would have reported."""
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    return status


def parse_compute_apps(text: str) -> List[Tuple[int, str, int]]:
    """``nvidia-smi --query-compute-apps=pid,gpu_uuid,used_memory`` into (pid, device, MiB).

    Rows the driver cannot fill come back as ``[Not Found]`` or ``[N/A]`` rather than being left
    out, and a row whose memory is unreadable is dropped rather than counted as zero: a zero would
    pull a peak down, while a missing row only fails to raise it.
    """
    out: List[Tuple[int, str, int]] = []
    for line in text.splitlines():
        parts = [cell.strip() for cell in line.split(",")]
        if len(parts) < 3:
            continue
        pid_text, device, memory_text = parts[0], parts[1], parts[2]
        digits = "".join(c for c in memory_text if c.isdigit())
        if not pid_text.isdigit() or not digits:
            continue
        out.append((int(pid_text), device, int(digits)))
    return out


def gpu_compute_apps() -> List[Tuple[int, str, int]]:
    """What is on the GPUs right now, or an empty list when there is no driver here."""
    executable = shutil.which("nvidia-smi")
    if not executable:
        return []
    try:
        result = subprocess.run(
            [
                executable,
                "--query-compute-apps=pid,gpu_uuid,used_memory",
                "--format=csv,noheader,nounits",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    return parse_compute_apps(result.stdout.decode("utf-8", "replace"))


def has_gpu() -> bool:
    """Whether to sample at all, decided once per process.

    ``nvidia-smi`` on PATH is the test, not a successful query: a driver that is present but busy
    can fail one call, and deciding from that would turn the sampler off for the rest of the run.
    """
    global _HAS_GPU
    if _HAS_GPU is None:
        _HAS_GPU = shutil.which("nvidia-smi") is not None
    return _HAS_GPU


def _parent_of(pid: int) -> Optional[int]:
    """The parent pid from ``/proc/<pid>/stat``, or None if it cannot be read.

    Read from field 4 after the comm field rather than by splitting the whole line: a process name
    can contain spaces and parentheses ("(sd-pam)"), and splitting on whitespace puts the parent in
    a different column for those.
    """
    try:
        with open("/proc/%d/stat" % pid) as handle:
            data = handle.read()
    except OSError:
        return None
    close = data.rfind(")")
    fields = data[close + 2 :].split() if close != -1 else []
    if len(fields) < 2 or not fields[1].lstrip("-").isdigit():
        return None
    return int(fields[1])


def in_process_tree(pid: int, root: int, limit: int = 64) -> bool:
    """True when *pid* is *root* or one of its descendants.

    Walks parents upward, which needs no bookkeeping as children come and go - and GPU processes
    here are grandchildren at least, since every step runs through ``mamba run``. *limit* stops a
    walk that a pid-reuse race could otherwise make long; the chain to pid 1 is a handful of steps.
    """
    current: Optional[int] = pid
    for _ in range(limit):
        if current is None or current <= 0:
            return False
        if current == root:
            return True
        if current == 1:
            return False
        current = _parent_of(current)
    return False


def peak_device_memory(samples: Iterable[List[Tuple[int, str, int]]], belongs: Any = None) -> int:
    """The largest MiB any single device held for the processes of interest.

    Per device rather than summed across devices, because the question a VRAM figure answers is
    "will this fit on one card". Per sample rather than per process, because two of these tools run
    several workers at once and a card has to hold them together.

    *belongs* decides which pids count; the default counts every process on the GPU, which is what
    a probe outside a run wants.
    """
    peak = 0
    for sample in samples:
        by_device: Dict[str, int] = {}
        for pid, device, memory in sample:
            if belongs is not None and not belongs(pid):
                continue
            by_device[device] = by_device.get(device, 0) + memory
        if by_device:
            peak = max(peak, max(by_device.values()))
    return peak


class GpuWatcher:
    """Sample GPU memory for one child's process tree until told to stop.

    A thread rather than a measurement at the end: the peak is transient, and by the time the child
    exits its allocations are gone. Failures are silent by design - a missing driver, a truncated
    nvidia-smi, a pid that exits between the sample and the /proc read - because a step must not
    fail over its own instrumentation.
    """

    def __init__(self, root_pid: int, interval: float = GPU_SAMPLE_SECONDS) -> None:
        self.root_pid = root_pid
        self.interval = interval
        self.peak_mib = 0
        self.devices: Set[str] = set()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> GpuWatcher:
        self._thread.start()
        return self

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                sample = gpu_compute_apps()
                mine = [row for row in sample if in_process_tree(row[0], self.root_pid)]
                self.peak_mib = max(self.peak_mib, peak_device_memory([mine]))
                self.devices.update(device for _pid, device, _memory in mine)
            except Exception:  # noqa: BLE001,S110 - instrumentation must never fail a step
                pass  # nosec B110 - the next sample is the retry; a log line per 5s would be noise
            self._stop.wait(self.interval)

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=self.interval + 5.0)


def run_measured(
    argv: Sequence[str],
    stdout: Optional[IO[Any]] = None,
    stderr: Optional[int] = None,
    env: Optional[Dict[str, str]] = None,
) -> Tuple[int, Dict[str, Any]]:
    """Run *argv* to completion. Returns ``(returncode, usage)``.

    *usage* holds ``wall_seconds`` always, plus ``peak_rss_bytes``, ``user_seconds`` and
    ``system_seconds`` when this platform can measure them. An unmeasurable platform yields the
    wall time alone, so callers never have to branch on :data:`CAN_MEASURE`.

    *env* is passed to the child only. The study driver uses it to tell a step that it is running
    inside a study, which is how the step knows not to print a second cost line beside the one the
    driver measures here. Setting that in this process instead would leak into every later call in
    the same interpreter - which is exactly what a test suite is.
    """
    started = time.time()
    if not CAN_MEASURE:
        rc = subprocess.call(list(argv), stdout=stdout, stderr=stderr, env=env)
        return rc, {"wall_seconds": round(time.time() - started, 1)}

    process = subprocess.Popen(list(argv), stdout=stdout, stderr=stderr, env=env)
    watcher = GpuWatcher(process.pid).start() if has_gpu() else None
    try:
        _pid, status, rusage = os.wait4(process.pid, 0)
    except ChildProcessError:
        # Something else reaped it (a stray SIGCHLD handler). The exit code is still available
        # through Popen, and losing the measurement must not lose the run.
        if watcher is not None:
            watcher.stop()
        return process.wait(), {"wall_seconds": round(time.time() - started, 1)}
    if watcher is not None:
        watcher.stop()
    # Popen still believes the child is running; tell it the truth so it does not try to reap it
    # again (which would raise) and so .returncode is correct for the caller.
    process.returncode = _exit_code(status)
    usage: Dict[str, Any] = {
        "wall_seconds": round(time.time() - started, 1),
        "peak_rss_bytes": maxrss_bytes(rusage.ru_maxrss),
        "user_seconds": round(rusage.ru_utime, 1),
        "system_seconds": round(rusage.ru_stime, 1),
    }
    # Only when the step actually put something on a GPU. A zero would be a claim - "this step was
    # measured and used no VRAM" - about every CPU step on a machine that happens to have a card.
    if watcher is not None and watcher.peak_mib:
        usage["peak_gpu_memory_bytes"] = watcher.peak_mib * MIB
        usage["gpu_devices"] = len(watcher.devices)
    return process.returncode, usage


def self_peak_rss_bytes() -> int:
    """Peak RSS of the MetaWrap2 process itself, for the record."""
    return maxrss_bytes(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def summarise(usage: Dict[str, Any]) -> str:
    """One line describing a step's cost, for a log message."""
    parts: List[str] = ["%.0fs wall" % usage.get("wall_seconds", 0)]
    if "peak_rss_bytes" in usage:
        parts.append("%s peak memory" % human(usage["peak_rss_bytes"]))
    if "user_seconds" in usage:
        cpu = usage.get("user_seconds", 0) + usage.get("system_seconds", 0)
        wall = usage.get("wall_seconds") or 1
        # Useful on its own: a step using 1.0 cores while given 24 threads is not threading.
        parts.append("%.1f cores busy" % (cpu / wall))
    if "peak_gpu_memory_bytes" in usage:
        # Said per device, because the figure is what one card had to hold.
        parts.append(
            "%s peak VRAM on %d GPU(s)"
            % (human(usage["peak_gpu_memory_bytes"]), usage.get("gpu_devices") or 1)
        )
    return ", ".join(parts)
