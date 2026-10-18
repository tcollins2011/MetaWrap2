"""One environment, one prefix: how `doctor` selects an env, and what it says when it cannot.

From a bug report against a fresh install. Existence and execution were answered against different
sets of conda installations:

* ``config.conda_env_exists`` tested every directory ``config.env_dirs`` could name, which includes
  every installation listed in ``~/.conda/environments.txt`` - a machine-global registry - plus five
  hardcoded ``~/miniforge3/envs``-style guesses that were appended unconditionally, although
  ``env_dirs``' own docstring said they were a fallback for a machine with no conda on PATH.
* ``doctor._env_run`` then ran ``mamba run -n <name>``, and ``-n`` resolves only against the active
  conda root.

So an env built by a second conda install was reported as existing and refused to run, `quickstart`
built nothing, and `doctor` reported every module broken. What the user actually saw was
``amrfinder MISSING "not on PATH"``, because ``_probe_tool`` turned every nonzero return into that
one message and discarded the output - which said ``The given prefix does not exist`` and named the
path. The wrong diagnosis cost the reporter 50 minutes.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from metawrap2 import config
from metawrap2.commands import doctor

#: A prefix under a conda installation that is not the active one - the reporter's case.
FOREIGN_PREFIX = "/opt/shared/miniforge3/envs/metawrap2-amr"
ACTIVE_ROOT = "/home/someone/miniforge3"

# ── fixtures ─────────────────────────────────────────────────────────────────────────────────


@pytest.fixture
def no_real_conda(monkeypatch):
    """Keep every test here off the conda installation this machine happens to have.

    ``doctor`` resolves a prefix through ``config.env_prefix``, which shells out to
    ``conda env list``; left alone, these tests would assert about the developer's own envs and
    give a different answer on CI. The per-process prefix cache is cleared on the way in *and* on
    the way out: it is module state, so a fake prefix left in it would be read by any later test
    that probes a tool.
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
    """``-p <prefix>``, never ``-n <name>``, whenever a prefix resolves.

    If this regresses, every env belonging to another conda installation goes back to being
    reported present and failing to run: ``-n`` looks only in the active root.
    """
    monkeypatch.setattr(doctor, "env_prefix", lambda env: FOREIGN_PREFIX)
    doctor._env_run("metawrap2-amr", "command -v amrfinder")
    cmd = mamba_argv[0]
    assert "-n" not in cmd
    assert cmd[cmd.index("-p") + 1] == FOREIGN_PREFIX


def test_env_run_falls_back_to_the_name_when_no_prefix_resolves(mamba_argv):
    """With nothing on disk to point at, ``-n`` is still better than refusing to run.

    An env that resolves to no prefix may still exist somewhere ``env_dirs`` cannot see (a
    ``envs_dirs`` entry only the manager's own config knows about), and ``-n`` is what finds it.
    Dropping this fallback would turn that into a hard failure.
    """
    doctor._env_run("metawrap2-amr", "command -v amrfinder")
    cmd = mamba_argv[0]
    assert "-p" not in cmd
    assert cmd[cmd.index("-n") + 1] == "metawrap2-amr"


def test_the_prefix_is_resolved_once_and_reused_for_every_probe(monkeypatch, no_real_conda):
    """One resolution per env, handed to each command that acts on it.

    Resolving walks ``env_dirs``, which runs ``conda env list``; a full `doctor` probes ~90 tools
    with up to three commands each. Resolving per command is both slow and a chance for two
    commands in one probe to disagree about which env they are in.
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
    """The message libmamba prints names the problem; it must reach the user instead of a guess.

    ``MISSING  not on PATH`` for a tool that is installed and on PATH inside its env is a wrong
    answer that sends the user looking at their env's contents rather than at which env was
    entered.
    """
    message = (
        "critical libmamba The given prefix does not exist: "
        '"/opt/shared/miniforge3/envs/metawrap2-amr"\n'
    )
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (1, message))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert status == doctor.MISSING
    assert detail == "env prefix does not exist: /opt/shared/miniforge3/envs/metawrap2-amr"


def test_a_failed_probe_keeps_what_the_subprocess_said(monkeypatch, no_real_conda):
    """Any other failure carries its last line of output, rather than discarding it.

    The output is the only evidence of what went wrong, and it is already captured. Dropping it
    again would mean every unforeseen failure - no mamba, a timeout, a broken wrapper script -
    reads identically and reads as a statement about PATH.
    """
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (127, "mamba not found"))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert status == doctor.MISSING
    assert "mamba not found" in detail


def test_a_tool_that_is_genuinely_absent_still_says_not_on_path(monkeypatch, no_real_conda):
    """`command -v` failing silently is the ordinary case and must stay a one-word answer.

    The extra detail added above has to be conditional: appending an empty parenthesis to every
    missing tool would make the common row harder to read than it was.
    """
    monkeypatch.setattr(doctor, "_env_run", lambda *a, **kw: (1, ""))
    status, detail = doctor._probe_tool("metawrap2-amr", "amrfinder")
    assert (status, detail) == (doctor.MISSING, "not on PATH")


# ── an env that belongs to a different conda installation ────────────────────────────────────


def test_an_env_outside_the_active_root_is_named_in_the_detail(monkeypatch, no_real_conda):
    """Probing it by prefix makes it work; the user still has to be told it is not theirs.

    `install-env` builds into the active root. If `doctor` silently reports on an env from another
    installation, the two commands discuss different environments under one name, and rebuilding
    changes nothing the status table shows.
    """
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", ACTIVE_ROOT)
    monkeypatch.setattr(doctor, "env_prefix", lambda env: FOREIGN_PREFIX)
    note = doctor._foreign_prefix_note("metawrap2-amr")
    assert FOREIGN_PREFIX in note
    assert ACTIVE_ROOT in note


def test_an_env_in_the_active_root_gets_no_such_note(monkeypatch, no_real_conda):
    """The normal install must stay quiet, or the note is noise everybody learns to ignore."""
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", ACTIVE_ROOT)
    monkeypatch.setattr(doctor, "env_prefix", lambda env: ACTIVE_ROOT + "/envs/metawrap2-amr")
    assert doctor._foreign_prefix_note("metawrap2-amr") == ""


def test_the_active_root_is_derived_from_the_active_env(monkeypatch, no_real_conda):
    """``CONDA_PREFIX`` points at the active *env*; the root is two levels up from a named one.

    Comparing a prefix against CONDA_PREFIX directly would call every env except the currently
    activated one foreign, which is the same false alarm in the other direction.
    """
    monkeypatch.setenv("CONDA_PREFIX", ACTIVE_ROOT + "/envs/metawrap2-reads")
    assert config.active_conda_root() == ACTIVE_ROOT
    assert config.prefix_in_active_root(ACTIVE_ROOT + "/envs/metawrap2-amr")
    assert not config.prefix_in_active_root(FOREIGN_PREFIX)


# ── which directories are searched at all ────────────────────────────────────────────────────

#: The five install locations env_dirs guesses at. Named here so the assertions below are about
#: the guesses rather than about whatever happens to be in the list.
GUESSES = ("miniforge3", "miniconda3", "mambaforge", "anaconda3", "micromamba")


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    """A HOME with no conda installs in it, so an expanded guess is a path this test controls."""
    monkeypatch.setenv("HOME", str(tmp_path))
    for var in ("CONDA_PREFIX", "CONDA_ROOT", "MAMBA_ROOT_PREFIX"):
        monkeypatch.delenv(var, raising=False)
    return tmp_path


def test_env_dirs_omits_the_hardcoded_guesses_when_conda_can_be_asked(monkeypatch, fake_home):
    """With an authoritative listing available, the guesses add nothing but other installations.

    They are what made ``conda_env_exists`` answer for a conda installation the user never chose:
    five directories searched that nobody asked about, on top of a listing that already covers
    every env the manager knows.
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
    """With nothing to ask, the guesses are the only answer there is.

    This is the case the docstring was always describing: a container or a batch node where conda
    exists on disk but is not on PATH until an activation script runs.
    """
    monkeypatch.setattr(config.shutil, "which", lambda prog: None)
    dirs = config.env_dirs()
    assert dirs == [os.path.join(str(fake_home), guess, "envs") for guess in GUESSES]
