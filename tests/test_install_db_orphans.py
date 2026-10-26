"""Orphaned downloads: a curl left by a killed `install-db`, and a later run resuming its file.

After a run was stopped, `ps` showed two curls writing the same file:

    44:06  curl ... -C - -o .../dbs/KOFAM/profiles.tar.gz.part https://.../profiles.tar.gz
    16:02  curl ... -C - -o .../dbs/KOFAM/profiles.tar.gz.part https://.../profiles.tar.gz

The 44-minute one was an orphan (`install-db` exited, its curl did not), and the new run resumed
the partial `.part` with ``-C -``. Two writers corrupt the file, and the later `tar` or `gzip`
failure looks like a bad download.

Tested here: a download runs in its own process group, which is killed when the command stops, and
a `.part` that another process appears to be writing is refused. A stale `.part` is still resumed.

Nothing here starts a process or touches the network: ``subprocess.Popen``, ``os.killpg`` and
``_run`` are replaced.
"""

from __future__ import annotations

import gzip
import os
import signal
import time
from typing import List

import pytest

from metawrap2.benchmark import fetch as _fetcher
from metawrap2.commands import install_db


@pytest.fixture(autouse=True)
def no_leaked_groups():
    """Keep the process-group registry out of the other tests, in both directions.

    It is module state: a fake pid left in it would reach ``os.killpg`` in a later test, and a pid
    from an earlier test could make a kill assertion pass incorrectly.
    """
    install_db._DOWNLOAD_GROUPS.clear()
    yield
    install_db._DOWNLOAD_GROUPS.clear()


class _FakeChild:
    """Stand-in for a ``subprocess.Popen`` result, with the two attributes ``_run`` uses."""

    def __init__(self, pid: int = 4321, returncode: int = 0):
        self.pid = pid
        self.stdout = None
        self._returncode = returncode

    def wait(self) -> int:
        return self._returncode


def test_a_fresh_part_file_is_refused_rather_than_resumed(tmp_path, monkeypatch, capsys):
    """A `.part` that another process is still writing is not resumed."""
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
    """A `.part` untouched for longer than the threshold is resumed.

    This is the usual retry after a stopped run, and restarting could cost hours of download.

    The stand-in must leave a `.part` that decompresses, because a download is accepted only if it
    does.
    """
    # Force the curl path; aria2c resumes from its control file instead.
    monkeypatch.setattr(_fetcher, "aria2c", lambda: None)
    dest = tmp_path / "profiles.tar.gz"
    part = tmp_path / "profiles.tar.gz.part"
    part.write_bytes(b"half a database")
    idle = install_db.LIVE_PART_SECONDS + 60
    os.utime(part, (time.time() - idle, time.time() - idle))

    commands: List[List[str]] = []

    def fake_run(cmd, **kw):
        commands.append(cmd)
        with open(cmd[cmd.index("-o") + 1], "wb") as fh:
            fh.write(gzip.compress(b"a resumed database"))
        return 0

    monkeypatch.setattr(install_db, "_run", fake_run)
    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is True
    assert commands, "a stale .part is resumed, not refused"
    assert "-C" in commands[0], "resuming is what makes a retry cheap"
    assert dest.is_file(), "the finished download replaces the .part"


def test_a_download_runs_in_its_own_process_group(monkeypatch):
    """Each download runs in its own tracked process group.

    The group lets an interrupt reach every process an installer started, including those a
    conda-env downloader starts from a shell script.
    """
    seen = {}

    def fake_popen(cmd, **kwargs):
        seen.update(kwargs)
        return _FakeChild()

    monkeypatch.setattr(install_db.subprocess, "Popen", fake_popen)
    assert install_db._run(["curl", "-o", "x.part", "https://example.invalid/x"]) == 0
    assert seen.get("start_new_session") is True


def test_stopping_signals_the_process_group_and_not_just_the_child(monkeypatch):
    """Stopping signals the whole group with killpg, not only the child.

    Without this, the separate session would stop Ctrl-C from reaching curl at all.
    """
    killed: List[tuple] = []
    monkeypatch.setattr(install_db.os, "killpg", lambda pgid, sig: killed.append((pgid, sig)))

    with install_db._tracked(_FakeChild(pid=4321)):
        install_db._terminate_downloads()

    assert killed == [(4321, signal.SIGTERM)]
    assert not install_db._DOWNLOAD_GROUPS, "a group that has been killed is not killed again"


def test_an_interrupted_run_leaves_no_download_running(tmp_path, monkeypatch):
    """Ctrl-C in `install-db` kills the downloads it started before returning.

    The kill must come before the thread pool is joined: the workers wait on children in another
    session, so waiting first would hang until the download finished.
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
