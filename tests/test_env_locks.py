"""Tests for the committed conda lockfiles and `metawrap2 doctor --fix`.

Neither of these may talk to conda in a test, so the package manager is stubbed. What is worth
pinning is the contract around it: lockfile naming and content, that a lockfile that fell back to
human-readable output is refused, restoring from a lockfile, and that ``--fix`` rebuilds exactly
the broken environments and nothing else.
"""

from __future__ import annotations

import os

import pytest

from metawrap2.commands import install_env as ie

EXPLICIT_OUTPUT = """# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/samtools-1.20-h50ea8bc_0.conda#abc123
https://conda.anaconda.org/bioconda/linux-64/metabat2-2.15-h986a166_1.tar.bz2#def456
"""


@pytest.fixture
def envs_dir(tmp_path, monkeypatch):
    """Point envs_dir() at a scratch directory so no test writes into the repo."""
    target = tmp_path / "envs"
    (target / "locks").mkdir(parents=True)
    monkeypatch.setenv("METAWRAP2_ENVS_DIR", str(target))
    return str(target)


# --- restoring from a lockfile ------------------------------------------------------------


def _queue_bins(envs_dir, monkeypatch, with_lock: bool):
    """Make `install-env bins` reach the build step, with or without a committed lockfile."""
    with open(os.path.join(envs_dir, "bins.yaml"), "w") as fh:
        fh.write("name: metawrap2-bins\ndependencies: [samtools]\n")
    if with_lock:
        with open(ie.lock_path("bins", "linux-64"), "w") as fh:
            fh.write(EXPLICIT_OUTPUT)
    monkeypatch.setattr(ie, "lock_platform", lambda: "linux-64")
    monkeypatch.setattr(ie, "ensure_mamba", lambda: True)
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: False)
    built = []

    def fake_create(module, env, yaml, force, log_dir, lock=None):
        built.append((env, lock))
        return module, True, os.path.join(log_dir, "x.log")

    monkeypatch.setattr(ie, "_create_one", fake_create)
    return built


def test_from_lock_installs_the_committed_lockfile(envs_dir, monkeypatch):
    """--from-lock was parsed and never used: every "reproducible" install solved the yaml.

    The lockfile reader had its own unit test, which passed, and nothing called it. So this goes
    through the command: naming the module (`binning`) must find the environment's lockfile
    (`bins.linux-64.lock`) and hand it to the build.
    """
    built = _queue_bins(envs_dir, monkeypatch, with_lock=True)
    assert ie.main(["binning", "--from-lock", "--no-test"]) == 0
    assert built == [("metawrap2-bins", ie.lock_path("bins", "linux-64"))]


def test_a_lockfile_is_named_after_the_environment_not_the_name_typed(envs_dir, monkeypatch):
    """`--lock binning` wrote binning.linux-64.lock, which --from-lock (keyed on bins) never reads."""
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: True)
    written = []
    monkeypatch.setattr(
        ie, "write_lockfile", lambda name, env, plat=None: written.append(name) or "x"
    )
    assert ie.main(["binning", "--lock"]) == 0
    assert written == ["bins"]
