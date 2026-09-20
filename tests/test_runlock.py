"""Tests for the one-writer-per-output-directory lock.

The property that matters is asymmetric: refusing to start when another run is live is cheap and
recoverable, while wrongly declaring a live run dead silently corrupts its output. So every
"cannot tell" case must resolve to "the lock is held".
"""

from __future__ import annotations

import json
import os
import socket

import pytest

from metawrap2 import runlock


def this_host() -> str:
    return socket.gethostname()


def write_lock(output_dir: str, **fields) -> str:
    os.makedirs(output_dir, exist_ok=True)
    record = {
        "pid": 999999,
        "host": this_host(),
        "user": "someone",
        "command": ["metawrap2", "binning"],
        "started": "2026-01-01T00:00:00",
    }
    record.update(fields)
    path = runlock.lock_path(output_dir)
    with open(path, "w") as fh:
        json.dump(record, fh)
    return path


# --- acquiring and releasing ---------------------------------------------------------------


def test_context_manager_acquires_and_releases(tmp_path):
    out = str(tmp_path / "study")
    with runlock.RunLock(out) as lock:
        assert lock.acquired
        assert os.path.isfile(runlock.lock_path(out))
    assert not os.path.exists(runlock.lock_path(out))


# --- contention ---------------------------------------------------------------------------


def test_a_live_lock_blocks_a_second_run(tmp_path):
    out = str(tmp_path / "study")
    write_lock(out, pid=os.getpid())  # this process is certainly alive
    with pytest.raises(runlock.LockHeld) as excinfo:
        runlock.RunLock(out).acquire()
    message = str(excinfo.value)
    assert out in message
    assert str(os.getpid()) in message
    assert "--force" in message  # the message must say how to get past it


def test_a_stale_lock_is_taken_over(tmp_path):
    """A machine rebooted mid-run leaves this. Making the user delete a file by hand is noise."""
    out = str(tmp_path / "study")
    write_lock(out, pid=999999)  # no such process
    note = runlock.RunLock(out).acquire()
    assert note and "stale" in note
    assert runlock.read(out)["pid"] == os.getpid()


def test_a_lock_from_another_host_is_believed(tmp_path):
    """A PID on another machine says nothing about a process here, so it must not be reclaimed."""
    out = str(tmp_path / "study")
    write_lock(out, pid=1, host="some-other-node")
    with pytest.raises(runlock.LockHeld):
        runlock.RunLock(out).acquire()
