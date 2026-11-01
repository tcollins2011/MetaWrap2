"""Tests for `install-db`'s disk ceiling, space estimate, download retries and failure messages.

No network and no real filesystem sizes: ``shutil.disk_usage`` is replaced with a fixed disk, and
downloads go through a stub ``_run``.
"""

from __future__ import annotations

from collections import namedtuple
from typing import List

import pytest

from metawrap2 import limits
from metawrap2.benchmark import fetch as _fetcher
from metawrap2.commands import install_db

GB = 1024**3
_Usage = namedtuple("_Usage", "total used free")


@pytest.fixture(autouse=True)
def fixed_disk(monkeypatch):
    """A 1000 GB disk with 400 GB used, and no `--max-*` overrides left from another test."""
    monkeypatch.setattr(
        limits.shutil, "disk_usage", lambda path: _Usage(1000 * GB, 400 * GB, 600 * GB)
    )
    limits.clear_overrides()
    yield
    limits.clear_overrides()


@pytest.fixture
def planning(monkeypatch):
    from metawrap2 import command

    monkeypatch.setattr(command.runner, "dry_run", True)
    monkeypatch.setattr(install_db.shutil, "which", lambda program: "/usr/bin/" + program)


# --- the disk ceiling comes from the config ----------------------------------------------


def test_the_configs_max_disk_fraction_sets_the_ceiling(tmp_path, capsys, planning):
    config = tmp_path / "config.toml"
    config.write_text("[limits]\nmax_disk_fraction = 0.85\n")
    install_db.main(["kegg", "--dir", str(tmp_path / "dbs"), "--config", str(config)])
    assert "ceiling 85%" in capsys.readouterr().out


def test_a_lower_ceiling_in_the_config_is_honoured_too(tmp_path, capsys, planning):
    config = tmp_path / "config.toml"
    config.write_text("[limits]\nmax_disk_fraction = 0.5\n")
    install_db.main(["kegg", "--dir", str(tmp_path / "dbs"), "--config", str(config)])
    assert "ceiling 50%" in capsys.readouterr().out


def test_the_flag_still_beats_the_config(tmp_path, capsys, planning):
    config = tmp_path / "config.toml"
    config.write_text("[limits]\nmax_disk_fraction = 0.85\n")
    limits.set_override("max_disk_fraction", "0.9")
    install_db.main(["kegg", "--dir", str(tmp_path / "dbs"), "--config", str(config)])
    assert "ceiling 90%" in capsys.readouterr().out


def test_without_a_config_the_default_ceiling_applies(tmp_path, capsys, planning, monkeypatch):
    monkeypatch.delenv("METAWRAP2_CONFIG", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    install_db.main(["kegg", "--dir", str(tmp_path / "dbs"), "--config", str(tmp_path / "no")])
    assert "ceiling %.0f%%" % (limits.DEFAULT_DISK_FRACTION * 100) in capsys.readouterr().out


# --- space for an archive that is deleted after unpacking --------------------------------


def test_gtdbtk_needs_its_download_and_its_unpacked_size_at_once():
    need, note = install_db._space_needed(["gtdbtk"], small=False)
    assert need == pytest.approx(220.0)
    assert note == "110 GB download + 110 GB unpacked before the archive is removed"
    line = install_db._disk_line("/tmp", need, note=note)
    assert line.startswith(
        "220.0 GB needed (110 GB download + 110 GB unpacked before the archive is removed), "
    )


@pytest.mark.parametrize("name", ["gunc", "eggnog", "singlem"])
def test_the_other_unpacked_databases_count_their_archive(name):
    db = install_db.DATABASES[name]
    need, note = install_db._space_needed([name], small=False)
    assert need == pytest.approx(
        install_db._gigabytes(db.download) + install_db._gigabytes(db.on_disk)
    )
    assert "before the archive is removed" in note


def test_a_database_without_a_held_archive_needs_its_on_disk_size_only():
    need, note = install_db._space_needed(["kegg"], small=False)
    assert need == pytest.approx(install_db._gigabytes(install_db.DATABASES["kegg"].on_disk))
    assert note == ""


def test_several_databases_name_the_archives_they_hold():
    need, note = install_db._space_needed(["gtdbtk", "gunc", "kegg"], small=False)
    assert need == pytest.approx(110 + 110 + 6.7 + 13 + 4 / 1024)
    assert note.endswith("of archives kept until they are unpacked: gtdbtk, gunc")


def test_a_dry_run_of_gtdbtk_asks_for_twice_its_size(tmp_path, capsys, planning):
    install_db.main(["gtdbtk", "--dir", str(tmp_path / "dbs"), "--config", str(tmp_path / "c")])
    out = capsys.readouterr().out
    assert "220.0 GB needed (110 GB download + 110 GB unpacked" in out


def test_a_need_above_the_headroom_is_warned_about():
    # 1000 GB disk, 400 GB used, 70% ceiling: 300 GB of headroom.
    assert "WARNING" not in install_db._disk_line("/tmp", 220.0)
    assert "WARNING" in install_db._disk_line("/tmp", 400.0)


# --- the SingleM note ---------------------------------------------------------------------


def test_the_singlem_note_does_not_blame_or_credit_aria2c(capsys):
    """`singlem data` downloads with its own client, so aria2c does not change its speed."""
    note = install_db.DATABASES["singlem"].describe
    assert "without aria2c" in note and "about 10 minutes" in note and "2-4 hours" in note
    assert "start it early" not in note
    install_db._print_list()
    assert note in capsys.readouterr().out


# --- download retries and the failure message --------------------------------------------


def _failing_downloads(monkeypatch, rc: int) -> List[List[str]]:
    """aria2c "installed", every run exiting *rc*, and waits recorded instead of slept."""
    calls: List[List[str]] = []
    waits: List[float] = []
    monkeypatch.setattr(_fetcher, "aria2c", lambda: "aria2c")
    monkeypatch.setattr(install_db, "_run", lambda cmd, **kw: calls.append(cmd) or rc)
    monkeypatch.setattr(install_db.time, "sleep", waits.append)
    return calls, waits


def test_an_unreachable_server_is_retried_after_long_waits(tmp_path, monkeypatch, capsys):
    calls, waits = _failing_downloads(monkeypatch, rc=2)
    assert install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz")) is False
    assert waits == list(install_db.UNREACHABLE_DELAYS)
    assert len(calls) == len(install_db.UNREACHABLE_DELAYS) + 1


def test_the_waits_can_be_passed_in(tmp_path, monkeypatch):
    calls, waits = _failing_downloads(monkeypatch, rc=6)
    install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz"), delays=[1, 2])
    assert waits == [1, 2] and len(calls) == 3


def test_other_failures_keep_the_short_retries(tmp_path, monkeypatch):
    calls, waits = _failing_downloads(monkeypatch, rc=3)  # not found: waiting does not help
    install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz"))
    assert len(calls) == _fetcher.ATTEMPTS
    assert waits == [_fetcher.RETRY_DELAY, _fetcher.RETRY_DELAY * 2]


def test_the_exit_code_is_put_into_words(tmp_path, monkeypatch, capsys):
    _failing_downloads(monkeypatch, rc=2)
    install_db._fetch("https://example.invalid/x.gz", str(tmp_path / "x.gz"), delays=[])
    out = capsys.readouterr().out
    assert "aria2c exited 2 (timed out)" in out
    assert "download failed after 1 attempt(s): https://example.invalid/x.gz" in out
    assert "attempt 1 of 1 failed" in out and out.count("attempt 1 of 1") == 1


def test_an_unknown_exit_code_is_left_as_a_number():
    assert install_db._explain("transport check failed: aria2c exited 99") == (
        "transport check failed: aria2c exited 99"
    )
    assert install_db._explain("transport check failed: curl exited 7").endswith(
        "(could not connect to the server)"
    )


def test_a_failed_gunc_install_suggests_skip_gunc(tmp_path, monkeypatch, capsys):
    _failing_downloads(monkeypatch, rc=2)
    monkeypatch.setattr(install_db, "UNREACHABLE_DELAYS", ())
    code = install_db.main(
        ["gunc", "--dir", str(tmp_path / "dbs"), "--config", str(tmp_path / "config.toml")]
    )
    out = capsys.readouterr().out
    assert code == 1
    assert "  FAILED: gunc" in out
    failed_block = out.split("  FAILED: gunc", 1)[1]
    assert "aria2c exited 2 (timed out)" in failed_block
    assert "--skip-gunc" in failed_block and "chimera check" in failed_block


def test_other_failed_databases_do_not_mention_skip_gunc(tmp_path, monkeypatch, capsys):
    _failing_downloads(monkeypatch, rc=2)
    monkeypatch.setattr(install_db, "UNREACHABLE_DELAYS", ())
    install_db.main(["kegg", "--dir", str(tmp_path / "dbs"), "--config", str(tmp_path / "c")])
    assert "--skip-gunc" not in capsys.readouterr().out
