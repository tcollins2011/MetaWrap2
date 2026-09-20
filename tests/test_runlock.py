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


def test_acquire_creates_a_lock_describing_this_process(tmp_path):
    out = str(tmp_path / "study")
    lock = runlock.RunLock(out)
    assert lock.acquire() is None  # nothing worth reporting
    record = runlock.read(out)
    assert record is not None
    assert record["pid"] == os.getpid()
    assert record["host"] == this_host()
    assert record["command"]


def test_acquire_creates_the_output_directory(tmp_path):
    out = str(tmp_path / "not" / "yet")
    runlock.RunLock(out).acquire()
    assert os.path.isfile(runlock.lock_path(out))


def test_release_removes_the_lock(tmp_path):
    out = str(tmp_path / "study")
    lock = runlock.RunLock(out)
    lock.acquire()
    lock.release()
    assert not os.path.exists(runlock.lock_path(out))


def test_release_only_removes_a_lock_this_object_took(tmp_path):
    """A run that never acquired must not be able to delete someone else's lock."""
    out = str(tmp_path / "study")
    path = write_lock(out, pid=os.getpid())
    runlock.RunLock(out).release()  # never acquired
    assert os.path.isfile(path)


def test_release_is_idempotent(tmp_path):
    lock = runlock.RunLock(str(tmp_path / "study"))
    lock.acquire()
    lock.release()
    lock.release()


def test_context_manager_acquires_and_releases(tmp_path):
    out = str(tmp_path / "study")
    with runlock.RunLock(out) as lock:
        assert lock.acquired
        assert os.path.isfile(runlock.lock_path(out))
    assert not os.path.exists(runlock.lock_path(out))


def test_context_manager_releases_even_when_the_run_raises(tmp_path):
    out = str(tmp_path / "study")
    with pytest.raises(RuntimeError), runlock.RunLock(out):
        raise RuntimeError("the run failed")
    assert not os.path.exists(runlock.lock_path(out))


def test_disabled_lock_does_nothing(tmp_path):
    out = str(tmp_path / "study")
    assert runlock.RunLock(out, enabled=False).acquire() is None
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


def test_two_locks_cannot_both_be_held(tmp_path):
    out = str(tmp_path / "study")
    first = runlock.RunLock(out)
    first.acquire()
    with pytest.raises(runlock.LockHeld):
        runlock.RunLock(out).acquire()
    first.release()
    runlock.RunLock(out).acquire()  # now free


def test_force_takes_a_live_lock_and_says_so(tmp_path):
    out = str(tmp_path / "study")
    write_lock(out, pid=os.getpid())
    note = runlock.RunLock(out, force=True).acquire()
    assert note and "--force" in note
    assert runlock.read(out)["pid"] == os.getpid()


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


def test_an_unreadable_lock_is_replaced(tmp_path):
    out = str(tmp_path / "study")
    os.makedirs(out)
    with open(runlock.lock_path(out), "w") as fh:
        fh.write("{ this is not json")
    note = runlock.RunLock(out).acquire()
    assert note and "unreadable" in note
    assert runlock.read(out)["pid"] == os.getpid()


def test_an_empty_lock_file_is_replaced(tmp_path):
    out = str(tmp_path / "study")
    os.makedirs(out)
    open(runlock.lock_path(out), "w").close()
    assert runlock.RunLock(out).acquire() is not None


# --- staleness detection ------------------------------------------------------------------


def test_is_stale_only_for_a_dead_pid_on_this_host():
    assert runlock.is_stale({"pid": 999999, "host": this_host()}) is True
    assert runlock.is_stale({"pid": os.getpid(), "host": this_host()}) is False
    assert runlock.is_stale({"pid": 1, "host": "elsewhere"}) is False


@pytest.mark.parametrize("record", [{}, {"host": "x"}, {"pid": "not a number"}])
def test_is_stale_is_false_for_a_record_it_cannot_read(record):
    """Every "cannot tell" resolves to "held" - the safe direction."""
    assert runlock.is_stale(record) is False


def test_read_of_a_missing_lock_is_none(tmp_path):
    assert runlock.read(str(tmp_path)) is None


def test_read_of_a_non_object_is_none(tmp_path):
    os.makedirs(str(tmp_path / "s"))
    with open(runlock.lock_path(str(tmp_path / "s")), "w") as fh:
        fh.write("[1, 2, 3]")
    assert runlock.read(str(tmp_path / "s")) is None


def test_describe_names_who_holds_it():
    text = runlock.describe({"pid": 42, "host": "node7", "user": "amy", "started": "2026-01-01"})
    assert "42" in text and "node7" in text and "amy" in text and "2026-01-01" in text
