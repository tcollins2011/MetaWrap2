"""One writer per output directory.

Two MetaWrap2 runs writing the same ``-o`` directory interleave: they overwrite each other's
intermediates, append to the same logs, and race on the manifest. The result is not a crash but
something worse - a plausible-looking output directory whose contents came from two runs with
different parameters. It happens easily: a resubmitted cluster job, a forgotten background run,
a `--resume` started while the original is still going.

:class:`RunLock` makes that fail immediately and legibly instead. It writes a small JSON file
holding the PID, host, user, command and start time, created with ``O_EXCL`` so the creation is
atomic. A second run finds it and reports who holds it.

**Stale locks are detected, not waited on.** A lock whose process is gone (same host, PID no
longer alive) is stale and is taken over silently - that is the normal state after a machine is
rebooted mid-run, and making the user delete a file by hand at that point is pointless
ceremony. A lock held on a *different* host cannot be checked, so it is believed; ``--force``
overrides. This is not a distributed lock and does not try to be: it is a guard against the
accident, on one filesystem, and says so when it cannot be sure.
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

    Only decidable when the lock was taken on this host: a PID on another machine says nothing
    about a process here. When in doubt this returns False, because wrongly declaring a live run
    stale is the expensive mistake.
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
    it; ``force=True`` takes it regardless, which is what the global ``--force`` means. Releasing
    only removes a lock this object created, so a forced second run cannot delete a third run's.
    """

    def __init__(self, output_dir: str, force: bool = False, enabled: bool = True):
        self.output_dir = output_dir
        self.force = force
        self.enabled = enabled
        self.path = lock_path(output_dir)
        self.acquired = False

    def acquire(self) -> Optional[str]:
        """Take the lock. Returns a note worth logging, or None. Raises LockHeld if it is held."""
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
                    # Unreadable or half-written: nothing to respect, so replace it.
                    note = "replaced an unreadable lock file at %s" % self.path
                elif is_stale(record):
                    note = "took over a stale lock (%s is no longer running)" % describe(record)
                elif self.force:
                    note = "--force: taking the lock from %s" % describe(record)
                else:
                    raise LockHeld(
                        "%s is already in use by another MetaWrap2 run (%s).\n"
                        "Two runs writing one output directory overwrite each other's "
                        "intermediates, so this one is stopping instead.\n"
                        "If that run is really gone, delete %s, or re-run with --force."
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

    # typing.Self is 3.11+, and pulling in typing_extensions for one annotation is not
    # worth a dependency. `from __future__ import annotations` makes this string-evaluated.
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
