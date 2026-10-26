"""Tests for the one-writer-per-output-directory lock.

Refusing to start when another run is live is recoverable; treating a live run as dead corrupts
its output. So every undecidable case must resolve to "the lock is held".
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


# --- contention ---------------------------------------------------------------------------


def test_a_live_lock_blocks_a_second_run(tmp_path):
    out = str(tmp_path / "study")
    write_lock(out, pid=os.getpid())  # this process is certainly alive
    with pytest.raises(runlock.LockHeld) as excinfo:
        runlock.RunLock(out).acquire()
    message = str(excinfo.value)
    assert out in message
    assert str(os.getpid()) in message
    assert "--force" in message  # the message names the override

    # A lock from a dead PID on this host (e.g. after a reboot) is taken over, and released on
    # exit.
    write_lock(out, pid=999999)  # no such process
    with runlock.RunLock(out) as lock:
        assert lock.acquired
        assert runlock.read(out)["pid"] == os.getpid()
    assert not os.path.exists(runlock.lock_path(out))
