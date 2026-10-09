"""Tests for the things a first-time or migrating user hits: quickstart, config migration,
retired names, and the disk-space check's severity.

The property that matters most here is the one about databases. They run from 550 MB to ~300 GB
and take hours to days, so `quickstart` must never start one without being told to - and in
particular must never be able to do so on a non-interactive run, where nobody is there to say no.

This is a curated subset: one strong test per distinct onboarding contract.
"""

from __future__ import annotations

import sys

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
    """The whole point: an unattended quickstart installs environments and nothing else."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert quickstart.main(["--binning-only"]) == 0
    assert calls["env"] is not None, "environments should still be created"
    assert calls["db"] is None, "no database download may be started without being asked"


def test_a_non_interactive_run_never_prompts(monkeypatch, capsys):
    """A batch script or CI job must not be able to trigger a 300 GB download."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    monkeypatch.setattr(
        "builtins.input", lambda *a: pytest.fail("must not prompt without a terminal")
    )
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


# --- config migration ----------------------------------------------------------------------

SHELL_CONFIG = """\
# Paths to custom pipelines and scripts
mw_path=$(which metawrap)
SOFT=${bin_path}/metawrap-scripts

KRAKEN_DB=/scratch/old/MY_KRAKEN_DB
KRAKEN2_DB=/scratch/db/kraken2
BMTAGGER_DB=/scratch/db/bmtagger
BLASTDB=/scratch/db/NCBI_nt
TAXDUMP=/scratch/db/NCBI_tax
SOMETHING_ELSE=/scratch/db/whatever
"""


def test_a_shell_config_is_translated():
    translated, retired, unrecognised = config_cmd.parse_shell_config(SHELL_CONFIG)
    assert translated == {
        "KRAKEN2_DB": "/scratch/db/kraken2",
        "BMTAGGER_DB": "/scratch/db/bmtagger",
        "BLASTDB": "/scratch/db/NCBI_nt",
        "TAXDUMP": "/scratch/db/NCBI_tax",
    }
    assert retired == ["KRAKEN_DB"]
    assert "SOMETHING_ELSE" in unrecognised


# --- a retired key must not block an unrelated run ----------------------------------------


def test_a_retired_database_key_is_a_note_not_a_failure():
    """It is harmless and only a removed module ever read it; refusing to run binning over it
    would be absurd."""
    from metawrap2.config import check_config

    problems, notes = check_config({"databases": {"KRAKEN_DB": "/x", "TAXDUMP": "/y"}}, "cfg.toml")
    assert problems == []
    assert len(notes) == 1 and "KRAKEN_DB" in notes[0] and "delete the line" in notes[0]


# --- retired CLI names explain themselves --------------------------------------------------


@pytest.mark.parametrize("name", ["kraken", "phylosift"])
def test_a_retired_module_name_explains_itself(name, capsys):
    from metawrap2 import cli

    assert cli.main([name, "-o", "out"]) == 2
    err = capsys.readouterr().err
    assert name in err and ("kraken2" in err or "classify_bins" in err)


# --- the disk-space check warns rather than refuses -----------------------------------------


# --- the starter config is the single source of truth --------------------------------------


def test_the_example_config_matches_the_starter_config():
    """`metawrap2.toml.example` is generated from DEFAULT_CONFIG.

    They were once two hand-maintained copies and the example fell behind, which is worse than
    having no example at all: a user pastes a section that no longer exists.
    """
    import pathlib

    example = pathlib.Path(__file__).resolve().parent.parent / "metawrap2.toml.example"
    assert example.read_text() == config_cmd.DEFAULT_CONFIG


def test_the_starter_config_is_valid_and_leaves_the_budget_at_its_defaults():
    """Everything in [limits] is commented out, so the file must parse and change nothing."""
    tomllib = pytest.importorskip("tomllib" if sys.version_info >= (3, 11) else "tomli")

    from metawrap2.config import check_config
    from metawrap2.limits import LimitsConfig, config_from_table

    data = tomllib.loads(config_cmd.DEFAULT_CONFIG)
    problems, _notes = check_config(data, "metawrap2.toml.example")
    assert problems == []
    assert config_from_table(data.get("limits", {})) == LimitsConfig()


def test_the_starter_config_does_not_pin_threads():
    """A hardcoded `threads = 8` would cap every module on a large machine at 8."""
    tomllib = pytest.importorskip("tomllib" if sys.version_info >= (3, 11) else "tomli")

    assert "threads" not in tomllib.loads(config_cmd.DEFAULT_CONFIG)["settings"]


# --- the shipped binning tier and what installing it takes ---------------------------------


def test_the_shipped_tier_and_the_drivers_default_cannot_disagree():
    """One definition of "the default binners", read by both the module and the pipeline driver.

    The driver used to hold its own literal ("--metabat2", "--maxbin2"), which is how a sheet
    naming --preset anything but fast ended up pointed at a bin directory that was never created.
    """
    from metawrap2.commands.run import DEFAULT_BINNERS
    from metawrap2.modules.binning import DEFAULT_PRESET, PRESETS

    assert DEFAULT_PRESET in PRESETS
    assert tuple(b.lstrip("-") for b in DEFAULT_BINNERS) == PRESETS[DEFAULT_PRESET]


def test_quickstart_installs_whatever_env_the_shipped_tier_needs():
    """The shipped tier reaches into an opt-in env, so an install that skipped it was incomplete.

    quickstart finished reporting everything healthy and the first `metawrap2 run` then stopped at
    binning asking for metawrap2-learned, which is the opposite of what quickstart is for.
    """
    from metawrap2.commands.quickstart import _default_tier_envs
    from metawrap2.config import OPTIONAL_ENVS
    from metawrap2.modules.binning import DEFAULT_PRESET, PRESETS

    needed = {
        OPTIONAL_ENVS[b].removeprefix("metawrap2-")
        for b in PRESETS[DEFAULT_PRESET]
        if b in OPTIONAL_ENVS
    }
    assert set(_default_tier_envs()) == needed
