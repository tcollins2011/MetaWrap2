"""A persistent record of every MetaWrap2 run started on this machine.

Per-run provenance is stored in each output directory and is lost with it. In addition, every
module invocation appends one line to a history file outside any run:

    ~/.metawrap2/history.jsonl        (or $METAWRAP2_HISTORY)

The file holds one JSON object per line and is append-only. Appending a line is safe enough for
concurrent runs without locking, and a truncated final line loses only that record. MetaWrap2
never rotates or prunes the file. A failed write is ignored so it cannot fail a run.

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

#: Default history file, in the user's MetaWrap2 directory beside the config.
DEFAULT_HISTORY = os.path.join("~", ".metawrap2", "history.jsonl")

_HISTORY_ENV = "METAWRAP2_HISTORY"

STARTED = "started"
COMPLETED = "completed"
FAILED = "failed"
INTERRUPTED = "interrupted"


def history_path(explicit: Optional[str] = None) -> str:
    """Return the history file path: explicit, then $METAWRAP2_HISTORY, then the default."""
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
        pass  # see the module docstring


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
    """Return filtered history entries, oldest first."""
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
    """Return history with one row per run, merging each run's start and finish entries.

    A module records a ``started`` entry and then a ``completed``/``failed``/``interrupted``
    one, so a killed run still leaves a trace. The two are merged on ``run_id``.
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
