"""Orphaned downloads: the curl a killed `install-db` leaves behind, and the next run that resumes it.

From one real session. A run was stopped, and `ps` then showed two curls writing the same file:

    44:06  curl ... -C - -o .../dbs/KOFAM/profiles.tar.gz.part https://.../profiles.tar.gz
    16:02  curl ... -C - -o .../dbs/KOFAM/profiles.tar.gz.part https://.../profiles.tar.gz

The 44-minute one was an orphan - `install-db` died, its curl did not - and the new run found a
partial `.part` and resumed it with ``-C -``. Two writers appending at their own offsets make a
file that is not the download, and the next `tar` or `gzip` fails looking like a bad download
rather than like a collision.

Two halves are pinned here: a download runs in its own process group and that group is killed when
this command stops, and a `.part` another process appears to be writing is refused rather than
resumed. Resume itself stays, for a `.part` that is genuinely stale - it is usually hours of work.

Nothing here starts a process or touches the network: ``subprocess.Popen``, ``os.killpg`` and
``_run`` are replaced.
"""

from __future__ import annotations

import os
import signal
import time
from typing import List

import pytest

from metawrap2.commands import install_db


@pytest.fixture(autouse=True)
def no_leaked_groups():
    """Keep the process-group registry out of the other tests, in both directions.

    It is module state, so a fake pid left in it would be handed to ``os.killpg`` by whatever ran
    next - and a stray pid from an earlier test would make a kill assertion here pass for the
    wrong reason.
    """
    install_db._DOWNLOAD_GROUPS.clear()
    yield
    install_db._DOWNLOAD_GROUPS.clear()


class _FakeChild:
    """What ``subprocess.Popen`` hands back, with the two attributes ``_run`` uses."""

    def __init__(self, pid: int = 4321, returncode: int = 0):
        self.pid = pid
        self.stdout = None
        self._returncode = returncode

    def wait(self) -> int:
        return self._returncode


def test_a_fresh_part_file_is_refused_rather_than_resumed(tmp_path, monkeypatch, capsys):
    """A `.part` someone else is still writing must not be resumed: two writers corrupt it.

    If this regresses, the second run appends from wherever the orphan's bytes currently end and
    the corruption surfaces an hour later as a `tar` failure.
    """
    commands: List[List[str]] = []

    def fake_run(cmd, **kw):
        commands.append(cmd)
        return 0

    monkeypatch.setattr(install_db, "_run", fake_run)
    dest = tmp_path / "profiles.tar.gz"
    part = tmp_path / "profiles.tar.gz.part"
    part.write_bytes(b"half a database")

    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is False
    out = capsys.readouterr().out
    assert "another install-db appears to be downloading this database" in out
    assert str(part) in out, "the message has to name the file to delete"
    assert commands == [], "no second curl may be started against a live .part"


def test_a_stale_part_file_is_still_resumed(tmp_path, monkeypatch):
    """Resume is the useful default and must survive the guard: it is usually hours of download.

    A `.part` nobody has touched for longer than the threshold is the ordinary case - a run that
    was stopped and is being retried - and starting that download again from zero is the cost this
    test exists to prevent.
    """
    dest = tmp_path / "profiles.tar.gz"
    part = tmp_path / "profiles.tar.gz.part"
    part.write_bytes(b"half a database")
    idle = install_db.LIVE_PART_SECONDS + 60
    os.utime(part, (time.time() - idle, time.time() - idle))

    commands: List[List[str]] = []

    def fake_run(cmd, **kw):
        commands.append(cmd)
        return 0

    monkeypatch.setattr(install_db, "_run", fake_run)
    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is True
    assert commands, "a stale .part is resumed, not refused"
    assert "-C" in commands[0], "resuming is what makes a retry cheap"
    assert dest.is_file(), "the finished download replaces the .part"


def test_a_download_runs_in_its_own_process_group(monkeypatch):
    """Downloads must not share this process's group, or they cannot be killed as a unit.

    The orphan is a curl that outlived the command that started it; a tracked group is what lets
    an interrupt reach every process a database's installer started, including the ones a conda-env
    downloader starts behind a shell script.
    """
    seen = {}

    def fake_popen(cmd, **kwargs):
        seen.update(kwargs)
        return _FakeChild()

    monkeypatch.setattr(install_db.subprocess, "Popen", fake_popen)
    assert install_db._run(["curl", "-o", "x.part", "https://example.invalid/x"]) == 0
    assert seen.get("start_new_session") is True


def test_stopping_signals_the_process_group_and_not_just_the_child(monkeypatch):
    """The group, through killpg: signalling only the child leaves its curl writing.

    This is the half that acts. Without it the new session above would make things *worse* - a
    Ctrl-C would no longer reach curl at all - so the two always travel together.
    """
    killed: List[tuple] = []
    monkeypatch.setattr(install_db.os, "killpg", lambda pgid, sig: killed.append((pgid, sig)))

    with install_db._tracked(_FakeChild(pid=4321)):
        install_db._terminate_downloads()

    assert killed == [(4321, signal.SIGTERM)]
    assert not install_db._DOWNLOAD_GROUPS, "a group that has been killed is not killed again"


def test_an_interrupted_run_leaves_no_download_running(tmp_path, monkeypatch):
    """Ctrl-C out of `install-db` must kill the downloads it started, before it returns.

    The reporter's orphan came from a stopped run. The kill also has to happen before the thread
    pool is joined: the workers are waiting on children that no longer share this process's
    session, so an executor shutdown that waited first would hang until the download finished.
    """
    killed: List[int] = []
    monkeypatch.setattr(install_db.os, "killpg", lambda pgid, sig: killed.append(pgid))

    def interrupted(name, target, small, jobs):
        install_db._DOWNLOAD_GROUPS.add(4321)
        raise KeyboardInterrupt

    monkeypatch.setattr(install_db, "_install_capturing", interrupted)
    with pytest.raises(KeyboardInterrupt):
        install_db.main(["card", "--dir", str(tmp_path), "--no-config-update"])

    assert killed == [4321]
