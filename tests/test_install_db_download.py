"""Tests for `install-db` downloads, which go through ``benchmark.fetch.download``.

The databases fetched by URL are ``.gz``, ``.tar.gz`` or ``.tar.bz2`` files without published
checksums. Each is test-decompressed before it is accepted, a rejected file is deleted, and a
re-run downloads it again.

Nothing here starts a process or opens a socket: ``install_db._run`` is replaced with a function
that writes the bytes the test chooses.
"""

from __future__ import annotations

import gzip
import os
from typing import Callable, List, Optional, Tuple

import pytest

from metawrap2.benchmark import fetch as _fetcher
from metawrap2.commands import install_db

PAYLOAD = b"ACGT" * 4096


@pytest.fixture(autouse=True)
def curl_and_no_waiting(monkeypatch):
    """Force the curl path and skip retry waits, so results do not depend on aria2c being found."""
    monkeypatch.setattr(_fetcher, "aria2c", lambda: None)
    monkeypatch.setattr(_fetcher.time, "sleep", lambda seconds: None)


def _arrivals(*payloads: Optional[bytes]) -> Tuple[List[List[str]], Callable[..., int]]:
    """Replace ``install_db._run`` with one that writes *payloads*, one per attempt, last sticky."""
    calls: List[List[str]] = []
    queue = list(payloads)

    def run(cmd: List[str], **kw) -> int:
        calls.append(list(cmd))
        payload = queue.pop(0) if len(queue) > 1 else (queue[0] if queue else None)
        if payload is None:
            return 1
        with open(cmd[cmd.index("-o") + 1], "wb") as fh:
            fh.write(payload)
        return 0

    return calls, run


def test_a_corrupt_archive_is_not_accepted_as_a_downloaded_database(tmp_path, monkeypatch, capsys):
    """A download that does not decompress is retried, then fails on the decompression check."""
    good = gzip.compress(PAYLOAD)
    calls, run = _arrivals(good[:-1] + bytes([good[-1] ^ 0xFF]))
    monkeypatch.setattr(install_db, "_run", run)
    dest = tmp_path / "profiles.tar.gz"

    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is False
    out = capsys.readouterr().out
    assert "decompression check failed" in out, out
    assert "does not decompress" in out
    assert len(calls) == _fetcher.ATTEMPTS, "a corrupt download is retried"


def test_a_rejected_download_leaves_nothing_for_the_next_run_to_resume(tmp_path, monkeypatch):
    """A rejected download leaves no `.part` and no file under the real name.

    Otherwise `curl -C -` would resume onto the bad bytes, or `_fetch` would report the file as
    already downloaded.
    """
    good = gzip.compress(PAYLOAD)
    _calls, run = _arrivals(good[:-1] + bytes([good[-1] ^ 0xFF]))
    monkeypatch.setattr(install_db, "_run", run)
    dest = tmp_path / "profiles.tar.gz"

    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is False
    assert sorted(os.listdir(tmp_path)) == []


def test_a_good_archive_lands_under_its_real_name_and_says_what_checked_it(
    tmp_path, monkeypatch, capsys
):
    """A good archive is renamed into place after one attempt, and the log names the check used."""
    good = gzip.compress(PAYLOAD)
    calls, run = _arrivals(good)
    monkeypatch.setattr(install_db, "_run", run)
    dest = tmp_path / "profiles.tar.gz"

    assert install_db._fetch("https://example.invalid/profiles.tar.gz", str(dest)) is True
    assert dest.read_bytes() == good
    assert len(calls) == 1
    assert "verified (decompresses cleanly)" in capsys.readouterr().out


def test_the_download_still_runs_through_install_dbs_own_run(tmp_path, monkeypatch):
    """The download is started through `install_db._run`, with curl resuming via ``-C -``.

    `_run` gives the child its own process group and captures its output per database.
    """
    calls, run = _arrivals(gzip.compress(PAYLOAD))
    monkeypatch.setattr(install_db, "_run", run)

    assert install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz")) is True
    assert calls[0][0] == "curl"
    assert "-C" in calls[0], "curl resumes a stale .part"


def test_a_ranged_downloader_is_used_when_one_is_available(tmp_path, monkeypatch):
    """aria2c with ranged connections is used when `benchmark.fetch.aria2c` finds it."""
    monkeypatch.setattr(_fetcher, "aria2c", lambda: "/opt/conda/bin/aria2c")
    calls: List[List[str]] = []

    def run(cmd: List[str], **kw) -> int:
        calls.append(list(cmd))
        target = os.path.join(cmd[cmd.index("--dir") + 1], cmd[cmd.index("--out") + 1])
        with open(target, "wb") as fh:
            fh.write(gzip.compress(PAYLOAD))
        return 0

    monkeypatch.setattr(install_db, "_run", run)
    assert install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz")) is True
    assert calls[0][0] == "/opt/conda/bin/aria2c"
    assert "--split=%d" % _fetcher.CONNECTIONS in calls[0]


def test_an_installer_fails_at_download_time_rather_than_at_unpacking(tmp_path, monkeypatch):
    """The GUNC installer fails at download time on a corrupt .dmnd.gz and leaves nothing behind."""
    good = gzip.compress(PAYLOAD)
    _calls, run = _arrivals(good[:-1] + bytes([good[-1] ^ 0xFF]))
    monkeypatch.setattr(install_db, "_run", run)
    target = tmp_path / "GUNC"

    assert install_db._install_gunc(str(target), False) is False
    assert sorted(os.listdir(target)) == [], "nothing is left to be reported as already downloaded"
