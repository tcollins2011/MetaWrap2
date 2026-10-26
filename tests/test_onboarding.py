"""Onboarding: quickstart, config migration, retired names, and the disk-space check's severity.

Databases range from 550 MB to ~300 GB, so `quickstart` must not download one unless asked, and
never on a non-interactive run.

One test per onboarding behaviour.
"""

from __future__ import annotations

import pytest

from metawrap2.commands import config_cmd, quickstart

# --- quickstart: databases are opt-in ------------------------------------------------------


@pytest.fixture
def calls(monkeypatch):
    """Record what quickstart would invoke, without invoking anything."""
    seen = {"env": None, "db": None, "doctor": None}

    monkeypatch.setattr(quickstart, "_package_manager", lambda: "mamba")
    monkeypatch.setattr(
        quickstart,
        "_run_install_env",
        lambda modules, threads, from_lock: seen.__setitem__("env", (list(modules), threads)) or 0,
    )
    monkeypatch.setattr(
        quickstart,
        "_run_install_db",
        lambda names, directory, small, jobs, config: seen.__setitem__(
            "db", (list(names), small, directory)
        )
        or 0,
    )
    monkeypatch.setattr(
        quickstart,
        "_run_doctor",
        lambda modules, config, databases: seen.__setitem__("doctor", databases) or 0,
    )
    return seen


def test_no_database_is_downloaded_by_default(calls, monkeypatch):
    """An unattended quickstart installs environments and no databases."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert quickstart.main(["--binning-only"]) == 0
    assert calls["env"] is not None, "environments should still be created"
    assert calls["db"] is None, "no database download may be started without being asked"

    # Without a terminal there is no prompt and no download.
    monkeypatch.setattr(
        "builtins.input", lambda *a: pytest.fail("must not prompt without a terminal")
    )
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


# --- the starter config is the single source of truth --------------------------------------


def test_the_example_config_matches_the_starter_config():
    """`metawrap2.toml.example` is generated from DEFAULT_CONFIG.

    They used to be maintained separately, and the example went out of date.
    """
    import pathlib

    example = pathlib.Path(__file__).resolve().parent.parent / "metawrap2.toml.example"
    assert example.read_text() == config_cmd.DEFAULT_CONFIG
