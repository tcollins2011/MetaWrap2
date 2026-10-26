"""Tests for how `doctor` selects an env by prefix and reports a failed probe.

Guards against a bug where env existence and execution were checked against different conda
installations:

* ``config.conda_env_exists`` searched every directory ``config.env_dirs`` returned: all
  installations in ``~/.conda/environments.txt`` plus five hardcoded ``~/miniforge3/envs``-style
  guesses, which were meant only as a fallback when no conda is on PATH.
* ``doctor._env_run`` ran ``mamba run -n <name>``, and ``-n`` resolves only in the active root.

An env from a second conda installation was then reported as present but failed to run, and
``_probe_tool`` reported ``not on PATH`` instead of libmamba's ``The given prefix does not exist``.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from metawrap2 import config
from metawrap2.commands import doctor

#: A prefix under a conda installation that is not the active one.
FOREIGN_PREFIX = "/opt/shared/miniforge3/envs/metawrap2-amr"
ACTIVE_ROOT = "/home/someone/miniforge3"

# ── fixtures ─────────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def no_real_conda(monkeypatch):
    """Isolate the tests from the host's conda installation.

    ``config.env_prefix`` runs ``conda env list``, so it is stubbed. The module-level prefix cache
    is cleared before and after, so a fake prefix does not leak into later tests.
    """
    doctor._PREFIX_CACHE.clear()
    monkeypatch.setattr(doctor, "env_prefix", lambda env: None)
    monkeypatch.setattr(doctor, "_run_flags", lambda env_manager: [])
    for var in ("CONDA_PREFIX", "CONDA_ROOT", "MAMBA_ROOT_PREFIX"):
        monkeypatch.delenv(var, raising=False)
    yield
    doctor._PREFIX_CACHE.clear()


@pytest.fixture
def mamba_argv(monkeypatch, no_real_conda):
    """Capture the argv `doctor` hands to mamba, answering it without starting a process."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(list(cmd))
        return subprocess.CompletedProcess(cmd, 0, stdout="/prefix/bin/amrfinder\n")

    monkeypatch.setattr(doctor.subprocess, "run", fake_run)
    return calls


# ── selecting the environment ────────────────────────────────────────────────────────────────


def test_env_run_enters_the_env_by_prefix(monkeypatch, mamba_argv):
    """When a prefix resolves, the env is entered with ``-p <prefix>``, not ``-n <name>``.

    ``-n`` looks only in the active root, so envs from other installations would fail to run.
    """
    monkeypatch.setattr(doctor, "env_prefix", lambda env: FOREIGN_PREFIX)
    doctor._env_run("metawrap2-amr", "command -v amrfinder")
    cmd = mamba_argv[0]
    assert "-n" not in cmd
    assert cmd[cmd.index("-p") + 1] == FOREIGN_PREFIX


def test_env_run_falls_back_to_the_name_when_no_prefix_resolves(mamba_argv):
    """With no resolved prefix, the env is entered with ``-n <name>``.

    The env may be in a directory ``env_dirs`` cannot see, such as an ``envs_dirs`` entry in the
    manager's own config, which ``-n`` still finds.
    """
    doctor._env_run("metawrap2-amr", "command -v amrfinder")
    cmd = mamba_argv[0]
    assert "-p" not in cmd
    assert cmd[cmd.index("-n") + 1] == "metawrap2-amr"


def test_the_prefix_is_resolved_once_and_reused_for_every_probe(monkeypatch, no_real_conda):
    """The prefix is resolved once per env and passed to each command in the probe.

    Resolving runs ``conda env list``, and a full `doctor` probes about 90 tools with up to three
    commands each. One resolution also keeps all commands in a probe on the same env.
    """
    resolved = []

    def fake_prefix(env):
        resolved.append(env)
        return FOREIGN_PREFIX

    used = []

    def fake_env_run(env, shell_cmd, timeout=60, prefix=None):
        used.append(prefix)
        return 0, "amrfinder 3.12.8"

    monkeypatch.setattr(doctor, "env_prefix", fake_prefix)
    monkeypatch.setattr(doctor, "_env_run", fake_env_run)
    doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert resolved == ["metawrap2-amr"]
    assert used == [FOREIGN_PREFIX] * len(used) and len(used) >= 2


# ── what a failed probe reports ──────────────────────────────────────────────────────────────


def test_a_nonexistent_prefix_is_reported_as_such_not_as_not_on_path(monkeypatch, no_real_conda):
    """libmamba's "prefix does not exist" error is reported as such, with the path."""
    message = (
        "critical libmamba The given prefix does not exist: "
        '"/opt/shared/miniforge3/envs/metawrap2-amr"\n'
    )
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (1, message))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert status == doctor.MISSING
    assert detail == "env prefix does not exist: /opt/shared/miniforge3/envs/metawrap2-amr"


def test_a_failed_probe_keeps_what_the_subprocess_said(monkeypatch, no_real_conda):
    """Other failures include the last line of the subprocess output in the detail.

    Without it, failures such as a missing mamba, a timeout or a broken wrapper script would all
    read as "not on PATH".
    """
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (127, "mamba not found"))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert status == doctor.MISSING
    assert "mamba not found" in detail


def test_a_tool_that_is_genuinely_absent_still_says_not_on_path(monkeypatch, no_real_conda):
    """A failure with no output (`command -v` not finding the tool) reports "not on PATH"."""
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (1, ""))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert (status, detail) == (doctor.MISSING, "not on PATH")


# ── an env that belongs to a different conda installation ────────────────────────────────────


def test_an_env_outside_the_active_root_is_named_in_the_detail(monkeypatch, no_real_conda):
    """An env outside the active conda root gets a note naming its prefix and the active root.

    `install-env` builds into the active root, so without the note `doctor` and `install-env`
    could refer to different envs with the same name.
    """
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", ACTIVE_ROOT)
    monkeypatch.setattr(doctor, "env_prefix", lambda env: FOREIGN_PREFIX)
    note = doctor._foreign_prefix_note("metawrap2-amr")
    assert FOREIGN_PREFIX in note
    assert ACTIVE_ROOT in note


def test_an_env_in_the_active_root_gets_no_such_note(monkeypatch, no_real_conda):
    """An env inside the active root gets no note."""
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", ACTIVE_ROOT)
    monkeypatch.setattr(doctor, "env_prefix", lambda env: ACTIVE_ROOT + "/envs/metawrap2-amr")
    assert doctor._foreign_prefix_note("metawrap2-amr") == ""


def test_the_active_root_is_derived_from_the_active_env(monkeypatch, no_real_conda):
    """The active root is two levels above ``CONDA_PREFIX`` when that is a named env.

    Comparing against CONDA_PREFIX directly would mark every env except the active one as foreign.
    """
    monkeypatch.setenv("CONDA_PREFIX", ACTIVE_ROOT + "/envs/metawrap2-reads")
    assert config.active_conda_root() == ACTIVE_ROOT
    assert config.prefix_in_active_root(ACTIVE_ROOT + "/envs/metawrap2-amr")
    assert not config.prefix_in_active_root(FOREIGN_PREFIX)


# ── which directories are searched at all ────────────────────────────────────────────────────

#: The five install locations that env_dirs guesses.
GUESSES = ("miniforge3", "miniconda3", "mambaforge", "anaconda3", "micromamba")


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    """A HOME with no conda installs, so expanded guesses point under tmp_path."""
    monkeypatch.setenv("HOME", str(tmp_path))
    for var in ("CONDA_PREFIX", "CONDA_ROOT", "MAMBA_ROOT_PREFIX"):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def test_env_dirs_omits_the_hardcoded_guesses_when_conda_can_be_asked(monkeypatch, fake_home):
    """When conda can list its envs, the hardcoded guesses are not added.

    The guesses made ``conda_env_exists`` find envs in other conda installations.
    """
    root = fake_home / "conda"
    (root / "envs" / "metawrap2-bins").mkdir(parents=True)
    listing = "# conda environments:\nbase  *  %s\nmetawrap2-bins  %s\n" % (
        root,
        root / "envs" / "metawrap2-bins",
    )
    monkeypatch.setattr(
        config.shutil, "which", lambda prog: "/usr/bin/conda" if prog == "conda" else None
    )
    monkeypatch.setattr(
        config.subprocess, "run", lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, listing)
    )
    dirs = config.env_dirs()
    assert str(root / "envs") in dirs
    assert not [d for d in dirs if any(g in d for g in GUESSES)]


def test_env_dirs_uses_the_hardcoded_guesses_when_no_manager_is_on_path(monkeypatch, fake_home):
    """With no conda or mamba on PATH, env_dirs returns the hardcoded guesses.

    This covers containers and batch nodes where conda is on disk but not on PATH until activated.
    """
    monkeypatch.setattr(config.shutil, "which", lambda prog: None)
    dirs = config.env_dirs()
    assert dirs == [os.path.join(str(fake_home), guess, "envs") for guess in GUESSES]
