"""CheckM2's DIAMOND reference: that it is installable, and that the installed path reaches the tool.

CheckM2's database is not CheckM1's. They share no files, both are called "the CheckM database", and
until this was in the catalog it was the one database MetaWrap2 told users to install by hand - so
the tests here are mostly about the two of them staying distinct.
"""

import pytest

from metawrap2 import checkm, config
from metawrap2.command import runner
from metawrap2.commands import install_db


def test_install_db_knows_checkm2_and_points_the_key_at_the_file():
    db = install_db.DATABASES["checkm2"]
    assert db.key == "CHECKM2DB"
    # CheckM2 takes --database_path to the .dmnd itself, not to the directory holding it.
    assert db.key_file and db.key_file.endswith(".dmnd")
    assert "checkm2" in install_db.INSTALLERS
    assert "CHECKM2DB" in config.DB_KEYS
    # Only the --checkm2 path needs it, so its absence is not a failed `metawrap2 check`.
    assert "CHECKM2DB" in config.OPTIONAL_DB_KEYS


def test_the_two_checkm_databases_are_separate_entries():
    """Installing one must not look like it installed the other."""
    assert install_db.DATABASES["checkm"].key == "CHECKM_DB"
    assert install_db.DATABASES["checkm2"].key == "CHECKM2DB"
    assert install_db.DATABASES["checkm2"].key != install_db.DATABASES["checkm"].key
    assert "not CheckM1's data" in install_db.DATABASES["checkm2"].describe


def test_the_installer_asks_checkm2_to_fetch_its_own_database(tmp_path, monkeypatch):
    """Not a URL of our own: the Zenodo record CheckM2 pins moves when CheckM2 is upgraded."""
    calls = []
    monkeypatch.setattr(install_db, "_mamba_run", lambda env, cmd: calls.append((env, cmd)) or 0)

    assert install_db._install_checkm2(str(tmp_path / "CHECKM2"), small=False) is False
    env, cmd = calls[0]
    assert env == config.OPTIONAL_ENVS["checkm2"]
    assert cmd[:4] == ["checkm2", "database", "--download", "--path"]


def test_an_already_installed_database_is_not_downloaded_again(tmp_path, monkeypatch):
    def explode(*_a, **_k):
        raise AssertionError("re-downloaded a database that was already there")

    monkeypatch.setattr(install_db, "_mamba_run", explode)
    dmnd = tmp_path / (install_db.DATABASES["checkm2"].key_file or "")
    dmnd.parent.mkdir(parents=True)
    dmnd.write_text("not really diamond, but not empty either")
    assert install_db._install_checkm2(str(tmp_path), small=False) is True


@pytest.fixture
def _capture_commands(monkeypatch):
    seen = []
    monkeypatch.setattr(checkm, "run", lambda cmd, **kw: seen.append(cmd))
    monkeypatch.setattr(runner, "dry_run", True)
    return seen


def test_the_configured_database_is_passed_on_the_command_line(_capture_commands, monkeypatch):
    """Explicitly, rather than trusting whatever state the conda env was left in."""
    monkeypatch.setattr(checkm, "checkm2_database", lambda: "/dbs/CHECKM2/x.dmnd")
    checkm.run_checkm2("bins", threads=4, env=None)
    assert "--database_path /dbs/CHECKM2/x.dmnd" in _capture_commands[0]


def test_no_configured_database_leaves_checkm2_to_find_its_own(_capture_commands, monkeypatch):
    """A hand-installed copy is registered inside CheckM2's env; passing "" would break it."""
    monkeypatch.setattr(checkm, "checkm2_database", lambda: "")
    checkm.run_checkm2("bins", threads=4, env=None)
    assert "--database_path" not in _capture_commands[0]


def test_an_unreadable_config_is_not_checkm2s_problem(monkeypatch):
    def explode(**_kw):
        raise RuntimeError("config is a mess")

    monkeypatch.setattr(config, "load_settings", explode)
    assert checkm.checkm2_database() == ""
