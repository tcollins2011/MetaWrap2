"""Measure the resources a step used: wall time, CPU time, peak RSS and peak VRAM.

Memory settings in MetaWrap2 (``-m 64`` for metaSPAdes, the space multipliers in
:mod:`metawrap2.scratch`) are estimates. Runs already record inputs, duration and outputs; this
module adds each step's peak memory.

:func:`run_measured` reaps the child with :func:`os.wait4` instead of :func:`subprocess.call`, which
returns the rusage of that process. ``getrusage(RUSAGE_CHILDREN)`` reports a running maximum over
all children ever waited for, so it cannot separate steps.

``ru_maxrss`` of a waited-for child includes its reaped descendants, so the measurement covers the
step's process tree (``mamba run -n env metaspades.py``, not just the shim).

Where ``os.wait4`` does not exist (Windows), the command runs without measurement.

``rusage`` does not cover GPU memory, which determines whether a step fits on a card. When
``nvidia-smi`` is present, a sampler thread records the peak VRAM the step's process tree held on a
single device. SemiBin2, VAMB and COMEBin (all in the ``learned`` environment) can use a GPU, and
the VRAM figures in the README's hardware table come from this measurement.

The same sampler tracks CPU use while a step runs. :class:`Utilisation` warns when a long step
keeps only a small part of its thread grant busy, as COMEBin and SemiBin2 clustering do for hours:
more threads would not help, and the idle cores are better used for other samples (``-j``).
"""

from __future__ import annotations

import os
import platform
import resource
import shutil
import subprocess
import sys
import threading
import time
from typing import IO, Any, Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

# Same elapsed-time format as the progress heartbeat, since both report on the same step.
from metawrap2.logfilter import duration

#: True when per-process resource usage can be measured here.
CAN_MEASURE = hasattr(os, "wait4")

#: Sampling interval. SemiBin2 and COMEBin allocate and free VRAM within seconds, and a 30 s
#: interval missed a training peak on a short bin set. An nvidia-smi call costs ~50 ms.
SAMPLE_SECONDS = 5.0

#: The interval's name from when the sampler watched only GPUs.
GPU_SAMPLE_SECONDS = SAMPLE_SECONDS

#: How long a step must run before its mean core use is checked (0 turns the check off), and the
#: fraction of its thread grant it must keep busy to avoid a warning. Every tool starts with a
#: single-threaded phase (reading an index, loading a model), so the mean is not meaningful in the
#: first minutes. Measured COMEBin runs kept 28% of their grant busy for hours; above a third,
#: extra threads mostly add contention (see metawrap2.command.THREAD_LIMITS).
LOW_UTILISATION_AFTER_SECONDS = 300.0
LOW_UTILISATION_FRACTION = 1.0 / 3.0

MIB = 1024 * 1024

#: The prefix on MetaWrap2's own lines in a stream of tool output, as in metawrap2.logfilter.
_PREFIX = "[metawrap2]"

#: Clock ticks per second, the unit of CPU time in ``/proc/<pid>/stat``. It is set at kernel build
#: time, so it is read from the system.
try:
    CLOCK_TICKS = float(os.sysconf("SC_CLK_TCK")) or 100.0
except (AttributeError, ValueError, OSError):  # pragma: no cover - non-POSIX
    CLOCK_TICKS = 100.0

#: Cached result of :func:`has_gpu`, so the PATH lookup happens once per process.
_HAS_GPU: Optional[bool] = None


def maxrss_bytes(raw: int) -> int:
    """Convert ``ru_maxrss`` to bytes.

    Linux reports kilobytes and macOS bytes, so the unit is chosen by platform.
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

    Unfilled fields come back as ``[Not Found]`` or ``[N/A]``. Rows with unreadable memory are
    dropped instead of counted as zero.
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
    """Current GPU compute processes, or an empty list when nvidia-smi is not available."""
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

    The test is ``nvidia-smi`` on PATH, not a successful query, because a busy driver can fail a
    single call.
    """
    global _HAS_GPU
    if _HAS_GPU is None:
        _HAS_GPU = shutil.which("nvidia-smi") is not None
    return _HAS_GPU


def _stat_fields(pid: int) -> List[str]:
    """``/proc/<pid>/stat`` from the state field onwards, or [] if it cannot be read.

    Parsed after the closing parenthesis of the comm field, because a process name can contain
    spaces and parentheses ("(sd-pam)"). Index 0 is field 3 of proc(5), the state.
    """
    try:
        with open("/proc/%d/stat" % pid) as handle:
            data = handle.read()
    except OSError:
        return []
    close = data.rfind(")")
    return data[close + 2 :].split() if close != -1 else []


def _parent_of(pid: int) -> Optional[int]:
    """The parent pid from ``/proc/<pid>/stat``, or None if it cannot be read."""
    return _ppid_from(_stat_fields(pid))


def _ppid_from(fields: List[str]) -> Optional[int]:
    """The ppid (proc(5) field 4) out of :func:`_stat_fields` output."""
    if len(fields) < 2 or not fields[1].lstrip("-").isdigit():
        return None
    return int(fields[1])


def _cpu_seconds_from(fields: List[str]) -> Optional[float]:
    """CPU seconds a process and its reaped descendants have used, out of _stat_fields output.

    ``utime + stime + cutime + cstime``, proc(5) fields 14 to 17 (indices 11 to 14 here). A
    parent's ``cutime`` holds only children it has waited for, so in a tree total a live process is
    counted through its own row and a reaped one through its parent's, each once.
    """
    if len(fields) < 15:
        return None
    try:
        ticks = sum(int(fields[index]) for index in (11, 12, 13, 14))
    except ValueError:
        return None
    return ticks / CLOCK_TICKS


def process_rows() -> Dict[int, Tuple[Optional[int], Optional[float]]]:
    """``pid -> (ppid, cpu seconds)`` for every process on this machine, in one pass over /proc.

    One read per process per sample. Walking parents upward per pid, as :func:`in_process_tree`
    does, would re-read the same ancestors for every candidate.
    """
    rows: Dict[int, Tuple[Optional[int], Optional[float]]] = {}
    try:
        names = os.listdir("/proc")
    except OSError:  # pragma: no cover - no /proc (macOS, Windows)
        return rows
    for name in names:
        if not name.isdigit():
            continue
        fields = _stat_fields(int(name))
        if not fields:
            # Exited between the listing and the read; its CPU time is now in its parent's row.
            continue
        rows[int(name)] = (_ppid_from(fields), _cpu_seconds_from(fields))
    return rows


def tree_cpu_seconds(
    root: int, rows: Optional[Dict[int, Tuple[Optional[int], Optional[float]]]] = None
) -> Optional[float]:
    """CPU seconds *root* and its descendants have used so far, or None if /proc cannot say.

    Walks down from *root* through a child index built from *rows*, which reads each process once.
    None when *root* is not in the snapshot: the step has exited, and ``os.wait4`` gives the exact
    figure.
    """
    rows = process_rows() if rows is None else rows
    if root not in rows:
        return None
    children: Dict[int, List[int]] = {}
    for pid, (ppid, _cpu) in rows.items():
        if ppid is not None:
            children.setdefault(ppid, []).append(pid)
    total = 0.0
    pending = [root]
    seen = {root}
    while pending:
        pid = pending.pop()
        cpu = rows.get(pid, (None, None))[1]
        if cpu is not None:
            total += cpu
        for child in children.get(pid, ()):
            # Guards against a cycle in a /proc snapshot taken during pid reuse.
            if child not in seen:
                seen.add(child)
                pending.append(child)
    return total


def in_process_tree(pid: int, root: int, limit: int = 64) -> bool:
    """True when *pid* is *root* or one of its descendants.

    Walks parents upward, so no tracking of child processes is needed. GPU processes are at least
    grandchildren, since steps run through ``mamba run``. *limit* bounds the walk in case of a
    pid-reuse race.
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

    Per device, because the figure is used to judge whether a step fits on one card. Summed over
    processes within a sample, because some tools run several workers on one card at once.

    *belongs* selects which pids count; by default every process on the GPU counts.
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


class Utilisation:
    """Mean cores busy against the thread grant, with a low-utilisation warning.

    Uses the CPU totals from the sampler's /proc pass (see :func:`tree_cpu_seconds`).

    Warns once per step. A repeated line every 5 s would add noise without new information.

    *clock* can be replaced in tests.
    """

    def __init__(
        self,
        threads: int,
        step: str = "",
        after_seconds: float = LOW_UTILISATION_AFTER_SECONDS,
        fraction: float = LOW_UTILISATION_FRACTION,
        enabled: bool = True,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self.threads = threads
        self.step = step
        self.after_seconds = after_seconds
        self.fraction = fraction
        #: after_seconds of 0 turns the check off. So does a grant of one thread, which cannot be
        #: underused.
        self.enabled = enabled and after_seconds > 0 and threads > 1
        self._clock = clock
        self._started = clock()
        self._warned = False

    def note(self, cpu_seconds: float) -> Optional[str]:
        """Record the step's CPU total so far. Returns the warning line, once, when it is due."""
        if not self.enabled or self._warned:
            return None
        elapsed = self._clock() - self._started
        if elapsed < self.after_seconds:
            return None
        busy = cpu_seconds / elapsed
        if busy > self.threads * self.fraction:
            return None
        self._warned = True
        return self.message(busy, elapsed)

    def message(self, busy: float, elapsed: float) -> str:
        """The warning text: what was measured, then what to do."""
        subject = self.step or "this step"
        return (
            "%s %s has kept %.1f of its %d granted cores busy over %s. It is working, not stuck: a "
            "step this far below its grant is in a phase that does not parallelise, so more "
            "threads will not shorten it. Those idle cores are better spent on other samples (-j) "
            "than on this step (-t).\n" % (_PREFIX, subject, busy, self.threads, duration(elapsed))
        )


def warn_on_screen(line: str) -> None:
    """Write a mid-step warning to the screen (stderr), not the step's log.

    stderr because a step's stdout can be a data stream (see :mod:`metawrap2.progress`); flushed so
    it appears while the step is running.
    """
    sys.stderr.write(line if line.endswith("\n") else line + "\n")
    sys.stderr.flush()


class StepWatcher:
    """Sample a running step's process tree until told to stop.

    A thread, because GPU allocations are freed before the child exits and the utilisation warning
    is only useful while the step runs.

    Errors (a missing driver, truncated nvidia-smi output, a pid that exits between the sample and
    the /proc read) are ignored, so a step never fails because of its instrumentation.
    """

    def __init__(
        self,
        root_pid: int,
        interval: float = SAMPLE_SECONDS,
        utilisation: Optional[Utilisation] = None,
        on_warning: Optional[Callable[[str], None]] = None,
    ) -> None:
        self.root_pid = root_pid
        self.interval = interval
        self.utilisation = utilisation
        self.peak_mib = 0
        self.devices: Set[str] = set()
        self._on_warning = on_warning if on_warning is not None else warn_on_screen
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def start(self) -> StepWatcher:
        self._thread.start()
        return self

    def _run(self) -> None:
        # Checked once per step, so a CPU-only machine does not do a PATH lookup every 5 s.
        gpu = has_gpu()
        while not self._stop.is_set():
            try:
                if self.utilisation is not None and self.utilisation.enabled:
                    cpu = tree_cpu_seconds(self.root_pid)
                    if cpu is not None:
                        warning = self.utilisation.note(cpu)
                        if warning:
                            self._on_warning(warning)
                if gpu:
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


#: The class's name from when it watched only GPUs, still used by ``modules/_common.py``.
GpuWatcher = StepWatcher


def grant_from_argv(argv: Sequence[str]) -> Optional[int]:
    """The thread grant in a step's own command line, or None if it carries none.

    ``run`` builds every step as ``metawrap2 <module> ... -t <threads>``, with its own ``-t`` after
    the module's arguments, so the last value wins. Callers that know the grant should pass
    ``threads=`` instead.
    """
    found: Optional[int] = None
    for index, token in enumerate(argv):
        value = ""
        if token in ("-t", "--threads") and index + 1 < len(argv):
            value = argv[index + 1]
        elif token.startswith("--threads="):
            value = token.split("=", 1)[1]
        if value.isdigit() and int(value) > 0:
            found = int(value)
    return found


def configured_utilisation(threads: Optional[int], step: str = "") -> Optional[Utilisation]:
    """A :class:`Utilisation` for a step with *threads* granted, or None when there is nothing to
    measure against.

    The thresholds are read from ``[settings]`` for each step. One TOML parse is negligible next to
    a step's run time, and a config edit takes effect without a cache to clear.
    """
    if not threads or threads < 2:
        return None
    from metawrap2.config import load_settings

    # strict=False: config errors are reported when the run loads the config, and a key that only
    # controls a warning must not fail a step.
    settings = load_settings(strict=False)
    after = settings.low_utilisation_after_seconds
    fraction = settings.low_utilisation_fraction
    return Utilisation(
        threads=threads,
        step=step,
        after_seconds=LOW_UTILISATION_AFTER_SECONDS if after is None else max(0.0, after),
        fraction=LOW_UTILISATION_FRACTION if fraction is None else fraction,
    )


def run_measured(
    argv: Sequence[str],
    stdout: Optional[IO[Any]] = None,
    stderr: Optional[int] = None,
    env: Optional[Dict[str, str]] = None,
    threads: Optional[int] = None,
    step: str = "",
    on_warning: Optional[Callable[[str], None]] = None,
) -> Tuple[int, Dict[str, Any]]:
    """Run *argv* to completion. Returns ``(returncode, usage)``.

    *usage* holds ``wall_seconds`` always, plus ``peak_rss_bytes``, ``user_seconds`` and
    ``system_seconds`` when this platform can measure them. An unmeasurable platform yields the
    wall time alone, so callers never have to branch on :data:`CAN_MEASURE`.

    *env* is passed to the child only. The study driver uses it to tell a step it is running inside
    a study, so the step does not print a second cost line. Setting it in this process would leak
    into later calls in the same interpreter, such as in a test suite.

    *threads* is the step's thread grant, for the utilisation warning; it defaults to
    :func:`grant_from_argv`. *step* names the step in the warning, and *on_warning* is where the
    warning goes (the screen by default, never the step's log).
    """
    started = time.time()
    if not CAN_MEASURE:
        rc = subprocess.call(list(argv), stdout=stdout, stderr=stderr, env=env)
        return rc, {"wall_seconds": round(time.time() - started, 1)}

    process = subprocess.Popen(list(argv), stdout=stdout, stderr=stderr, env=env)
    utilisation = configured_utilisation(
        threads if threads is not None else grant_from_argv(argv), step
    )
    # The sampler runs for the utilisation check on any machine, and for VRAM when there is a GPU.
    sampling = (utilisation is not None and utilisation.enabled) or has_gpu()
    watcher = (
        StepWatcher(process.pid, utilisation=utilisation, on_warning=on_warning).start()
        if sampling
        else None
    )
    try:
        _pid, status, rusage = os.wait4(process.pid, 0)
    except ChildProcessError:
        # Something else reaped it (a stray SIGCHLD handler). Return Popen's exit code without
        # the measurement.
        if watcher is not None:
            watcher.stop()
        return process.wait(), {"wall_seconds": round(time.time() - started, 1)}
    if watcher is not None:
        watcher.stop()
    # Set returncode so Popen does not try to reap the child again (which would raise) and the
    # caller sees the right value.
    process.returncode = _exit_code(status)
    usage: Dict[str, Any] = {
        "wall_seconds": round(time.time() - started, 1),
        "peak_rss_bytes": maxrss_bytes(rusage.ru_maxrss),
        "user_seconds": round(rusage.ru_utime, 1),
        "system_seconds": round(rusage.ru_stime, 1),
    }
    # Only when the step used a GPU; otherwise every CPU step on a GPU machine would report 0 VRAM.
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
        # Shows whether a step uses the threads it was given.
        parts.append("%.1f cores busy" % (cpu / wall))
    if "peak_gpu_memory_bytes" in usage:
        # Per device: the figure is what one card had to hold.
        parts.append(
            "%s peak VRAM on %d GPU(s)"
            % (human(usage["peak_gpu_memory_bytes"]), usage.get("gpu_devices") or 1)
        )
    return ", ".join(parts)
