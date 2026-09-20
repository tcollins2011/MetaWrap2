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
"""

from __future__ import annotations

import os
import platform
import resource
import subprocess
import time
from typing import IO, Any, Dict, List, Optional, Sequence, Tuple

#: True when per-process resource usage can be measured here.
CAN_MEASURE = hasattr(os, "wait4")


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


def run_measured(
    argv: Sequence[str],
    stdout: Optional[IO[Any]] = None,
    stderr: Optional[int] = None,
) -> Tuple[int, Dict[str, Any]]:
    """Run *argv* to completion. Returns ``(returncode, usage)``.

    *usage* holds ``wall_seconds`` always, plus ``peak_rss_bytes``, ``user_seconds`` and
    ``system_seconds`` when this platform can measure them. An unmeasurable platform yields the
    wall time alone, so callers never have to branch on :data:`CAN_MEASURE`.
    """
    started = time.time()
    if not CAN_MEASURE:
        rc = subprocess.call(list(argv), stdout=stdout, stderr=stderr)
        return rc, {"wall_seconds": round(time.time() - started, 1)}

    process = subprocess.Popen(list(argv), stdout=stdout, stderr=stderr)
    try:
        _pid, status, rusage = os.wait4(process.pid, 0)
    except ChildProcessError:
        # Something else reaped it (a stray SIGCHLD handler). The exit code is still available
        # through Popen, and losing the measurement must not lose the run.
        return process.wait(), {"wall_seconds": round(time.time() - started, 1)}
    # Popen still believes the child is running; tell it the truth so it does not try to reap it
    # again (which would raise) and so .returncode is correct for the caller.
    process.returncode = _exit_code(status)
    usage: Dict[str, Any] = {
        "wall_seconds": round(time.time() - started, 1),
        "peak_rss_bytes": maxrss_bytes(rusage.ru_maxrss),
        "user_seconds": round(rusage.ru_utime, 1),
        "system_seconds": round(rusage.ru_stime, 1),
    }
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
    return ", ".join(parts)
