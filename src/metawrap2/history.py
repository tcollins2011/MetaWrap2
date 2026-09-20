"""A durable record of every MetaWrap2 run ever started on this machine.

A run's provenance lives in its own output directory, which answers "how was *this* made?".
It cannot answer "what have I run?", "when did I last run binning, and with what?", or "which
analyses used the database I just replaced" - and it disappears with the output directory.

So every module invocation also appends one line to a history file outside any run:

    ~/.metawrap2/history.jsonl        (or $METAWRAP2_HISTORY)

One JSON object per line, append-only. That format is chosen deliberately: appending a line is
atomic enough to be safe from concurrent runs without locking, a truncated final line costs one
record rather than the file, and the whole thing stays greppable. It is never rotated or pruned
by MetaWrap2 - the user's own record outlives any individual analysis - and a failed write is
ignored, because losing a history line must never be the reason a run fails.

Read it with ``metawrap2 history``.
"""

from __future__ import annotations

import getpass
import json
import os
import platform
import socket
import time
from typing import Any, Dict, Iterator, List, Optional, Sequence

#: Where the history lives. Under the user's MetaWrap2 directory, beside the config, so it
#: survives for as long as MetaWrap2 is installed.
DEFAULT_HISTORY = os.path.join("~", ".metawrap2", "history.jsonl")

_HISTORY_ENV = "METAWRAP2_HISTORY"

STARTED = "started"
COMPLETED = "completed"
FAILED = "failed"
INTERRUPTED = "interrupted"


def history_path(explicit: Optional[str] = None) -> str:
    """The history file to use: explicit, then $METAWRAP2_HISTORY, then the default."""
    path = explicit or os.environ.get(_HISTORY_ENV) or DEFAULT_HISTORY
    return os.path.abspath(os.path.expanduser(path))


def _who() -> Dict[str, str]:
    try:
        user = getpass.getuser()
    except (KeyError, OSError):
        user = "?"
    try:
        host = socket.gethostname()
    except OSError:
        host = "?"
    return {"user": user, "host": host}


def record(
    module: str,
    *,
    version: str,
    run_id: str,
    command: Sequence[str],
    output: Optional[str] = None,
    status: str = STARTED,
    duration_seconds: Optional[float] = None,
    databases: Optional[Dict[str, str]] = None,
    conda_env: Optional[str] = None,
    path: Optional[str] = None,
) -> None:
    """Append one entry. Never raises: a history failure must not fail a run."""
    entry: Dict[str, Any] = {
        "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "module": module,
        "status": status,
        "run_id": run_id,
        "version": version,
        "command": list(command),
        "cwd": os.getcwd(),
        "python": platform.python_version(),
    }
    entry.update(_who())
    if output:
        entry["output"] = os.path.abspath(output)
    if duration_seconds is not None:
        entry["duration_seconds"] = round(duration_seconds, 1)
    if conda_env:
        entry["conda_env"] = conda_env
    if databases:
        entry["databases"] = {k: v for k, v in databases.items() if v}
    target = history_path(path)
    try:
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(target, "a") as fh:
            fh.write(json.dumps(entry, sort_keys=True) + "\n")
    except OSError:
        pass  # deliberately silent: see the module docstring


def read(path: Optional[str] = None) -> Iterator[Dict[str, Any]]:
    """Yield every history entry, oldest first, skipping any unreadable line."""
    target = history_path(path)
    if not os.path.isfile(target):
        return
    with open(target, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except ValueError:
                continue  # a partially-written final line costs one record
            if isinstance(entry, dict):
                yield entry


def entries(
    path: Optional[str] = None,
    module: Optional[str] = None,
    status: Optional[str] = None,
    since: Optional[str] = None,
    output: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Filtered history, oldest first."""
    found = []
    for entry in read(path):
        if module and entry.get("module") != module:
            continue
        if status and entry.get("status") != status:
            continue
        if since and entry.get("time", "") < since:
            continue
        if output and os.path.abspath(output) != entry.get("output"):
            continue
        found.append(entry)
    return found


def runs(path: Optional[str] = None, **filters: object) -> List[Dict[str, Any]]:
    """History collapsed to one row per run, merging each run's start and finish entries.

    A module records a ``started`` entry and then a ``completed``/``failed``/``interrupted``
    one, so that an interrupted or killed run still leaves a trace. For display, the two are
    folded together on ``run_id``.
    """
    merged: Dict[str, Dict[str, Any]] = {}
    order: List[str] = []
    for entry in entries(path, **filters):  # type: ignore[arg-type]
        run_id = str(entry.get("run_id") or entry.get("time") or "")
        if run_id not in merged:
            merged[run_id] = dict(entry)
            order.append(run_id)
            continue
        existing = merged[run_id]
        # A later entry is the outcome: keep its status and duration, keep the original start.
        existing["status"] = entry.get("status", existing.get("status"))
        if "duration_seconds" in entry:
            existing["duration_seconds"] = entry["duration_seconds"]
        existing["finished"] = entry.get("time")
    return [merged[r] for r in order]
