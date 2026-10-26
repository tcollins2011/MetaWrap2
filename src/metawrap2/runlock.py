"""One writer per output directory.

Two runs writing the same ``-o`` directory overwrite each other's intermediates, append to the
same logs and race on the manifest, leaving output from two runs with different parameters
without an error. Common causes are a resubmitted cluster job, a forgotten background run, or a
`--resume` started while the original is still running.

:class:`RunLock` writes a small JSON file (PID, host, user, command, start time) created with
``O_EXCL``, so creation is atomic. A second run finds it, reports the holder and stops.

A lock whose process is gone (same host, PID not alive) is stale and is taken over, which is the
usual case after a reboot mid-run. A lock from a different host cannot be checked and is honoured;
``--force`` overrides it. This guards one filesystem against accidental concurrent runs; it is not
a distributed lock.
"""

from __future__ import annotations

import errno
import getpass
import json
import os
import socket
import sys
import time
from types import TracebackType
from typing import Any, Dict, Optional, Type

#: Lock filename inside the output directory.
LOCK_NAME = ".metawrap2.lock"


class LockHeld(Exception):
    """The output directory is locked by another run."""


def lock_path(output_dir: str) -> str:
    return os.path.join(output_dir, LOCK_NAME)


def _identity() -> Dict[str, Any]:
    try:
        user = getpass.getuser()
    except (KeyError, OSError):
        user = "?"
    try:
        host = socket.gethostname()
    except OSError:
        host = "?"
    return {
        "pid": os.getpid(),
        "host": host,
        "user": user,
        "command": list(sys.argv),
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }


def read(output_dir: str) -> Optional[Dict[str, Any]]:
    """The lock record held on *output_dir*, or None if there is none or it is unreadable."""
    try:
        with open(lock_path(output_dir)) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _process_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, owned by someone else
    except OSError:
        return True  # cannot tell; assume it is there
    return True


def is_stale(record: Dict[str, Any]) -> bool:
    """True if the holding process is definitely gone.

    Only decidable for a lock taken on this host. Returns False when unsure, since treating a live
    run as stale is the costlier error.
    """
    try:
        host = socket.gethostname()
    except OSError:
        return False
    if record.get("host") != host:
        return False
    pid = record.get("pid")
    if not isinstance(pid, int):
        return False
    return not _process_alive(pid)


def describe(record: Dict[str, Any]) -> str:
    return "pid %s on %s as %s, started %s" % (
        record.get("pid", "?"),
        record.get("host", "?"),
        record.get("user", "?"),
        record.get("started", "?"),
    )


class RunLock:
    """Hold the lock on an output directory for the duration of a run.

    Used as a context manager. ``acquire()`` raises :class:`LockHeld` when another live run holds
    it; ``force=True`` (the global ``--force``) takes it anyway. Releasing only removes a lock this
    object created, so a forced second run cannot delete a third run's.
    """

    def __init__(self, output_dir: str, force: bool = False, enabled: bool = True):
        self.output_dir = output_dir
        self.force = force
        self.enabled = enabled
        self.path = lock_path(output_dir)
        self.acquired = False

    def acquire(self) -> Optional[str]:
        """Take the lock. Returns a note to log, or None. Raises LockHeld if it is held."""
        if not self.enabled:
            return None
        os.makedirs(self.output_dir, exist_ok=True)
        payload = json.dumps(_identity(), indent=2, sort_keys=True) + "\n"
        note = None
        while True:
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            except FileExistsError:
                record = read(self.output_dir)
                if record is None:
                    # Unreadable or half-written lock file: replace it.
                    note = "replaced an unreadable lock file at %s" % self.path
                elif is_stale(record):
                    note = "took over a stale lock (%s is no longer running)" % describe(record)
                elif self.force:
                    note = "--force: taking the lock from %s" % describe(record)
                else:
                    raise LockHeld(
                        "%s is already in use by another MetaWrap2 run (%s).\n"
                        "Two runs writing one output directory overwrite each other's "
                        "intermediates, so this run is stopping.\n"
                        "If that run is no longer running, delete %s, or re-run with --force."
                        % (self.output_dir, describe(record), self.path)
                    )
                try:
                    os.remove(self.path)
                except OSError as exc:
                    if exc.errno != errno.ENOENT:
                        raise
                continue
            with os.fdopen(fd, "w") as fh:
                fh.write(payload)
            self.acquired = True
            return note

    def release(self) -> None:
        """Remove the lock, but only if this object took it."""
        if not self.acquired:
            return
        self.acquired = False
        try:
            os.remove(self.path)
        except OSError:
            pass  # already gone, or the directory was removed under us

    # typing.Self needs 3.11+ or typing_extensions; `from __future__ import annotations` makes
    # this annotation a string.
    def __enter__(self) -> RunLock:  # noqa: PYI034
        self.acquire()
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc: Optional[BaseException],
        tb: Optional[TracebackType],
    ) -> None:
        self.release()
