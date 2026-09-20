"""Tests for the committed conda lockfiles and `metawrap2 doctor --fix`.

Neither of these may talk to conda in a test, so the package manager is stubbed. What is worth
pinning is the contract around it: lockfile naming and content, that a lockfile that fell back to
human-readable output is refused, restoring from a lockfile, and that ``--fix`` rebuilds exactly
the broken environments and nothing else.
"""

from __future__ import annotations

import os
import subprocess

import pytest

from metawrap2.commands import doctor
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


def read(path: str) -> str:
    with open(path) as fh:
        return fh.read()


def stub_conda(monkeypatch, stdout=EXPLICIT_OUTPUT, returncode=0, exc=None):
    monkeypatch.setattr(
        ie.shutil, "which", lambda name: "/usr/bin/conda" if name == "conda" else None
    )

    def fake_run(cmd, **kwargs):
        if exc is not None:
            raise exc
        return subprocess.CompletedProcess(cmd, returncode, stdout=stdout, stderr="")

    monkeypatch.setattr(ie.subprocess, "run", fake_run)


# --- naming -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "machine,system,expected",
    [
        ("x86_64", "Linux", "linux-64"),
        ("AMD64", "Windows", "win-64"),
        ("arm64", "Darwin", "osx-arm64"),
        ("aarch64", "Linux", "linux-aarch64"),
    ],
)
def test_lock_platform_translates_each_machine(monkeypatch, machine, system, expected):
    monkeypatch.setattr(ie.platform, "machine", lambda: machine)
    monkeypatch.setattr(ie.platform, "system", lambda: system)
    assert ie.lock_platform() == expected


# --- writing a lockfile -------------------------------------------------------------------


def test_write_lockfile_records_explicit_urls_with_checksums(envs_dir, monkeypatch):
    stub_conda(monkeypatch)
    path = ie.write_lockfile("binning", "metawrap2-binning", "linux-64")
    assert path == ie.lock_path("binning", "linux-64")
    content = read(path)
    assert "@EXPLICIT" in content
    assert "samtools-1.20-h50ea8bc_0.conda#abc123" in content  # exact build + md5, not a range


def test_write_lockfile_refuses_output_that_is_not_explicit(envs_dir, monkeypatch):
    """A `conda list` that fell back to human-readable output must not be saved as a lockfile."""
    stub_conda(monkeypatch, stdout="samtools 1.20 h50ea8bc_0 conda-forge\n")
    assert ie.write_lockfile("binning", "metawrap2-binning", "linux-64") is None
    assert not os.path.exists(ie.lock_path("binning", "linux-64"))


# --- restoring from a lockfile ------------------------------------------------------------


def test_create_from_lock_invokes_conda_create_with_the_file(envs_dir, monkeypatch):
    lock = ie.lock_path("binning", "linux-64")
    with open(lock, "w") as fh:
        fh.write(EXPLICIT_OUTPUT)
    monkeypatch.setattr(ie.shutil, "which", lambda n: "/usr/bin/conda" if n == "conda" else None)
    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        return 0

    monkeypatch.setattr(ie, "run", fake_run)
    assert ie.create_from_lock("binning", "metawrap2-binning", "linux-64") is True
    assert seen["cmd"] == [
        "/usr/bin/conda",
        "create",
        "-y",
        "-n",
        "metawrap2-binning",
        "--file",
        lock,
    ]


# --- doctor --fix -------------------------------------------------------------------------


def results(**statuses):
    return [
        {"module": m, "env": "metawrap2-%s" % m, "status": s, "tools": []}
        for m, s in statuses.items()
    ]


def test_fix_rebuilds_only_the_broken_environments(monkeypatch, capsys):
    """Rebuilding a healthy env would throw away a working install for no reason."""
    seen = {}
    monkeypatch.setattr(
        "metawrap2.commands.install_env.main", lambda argv: seen.setdefault("argv", argv) and 0
    )
    monkeypatch.setattr(
        doctor,
        "_check_module",
        lambda m: {"module": m, "env": "metawrap2-" + m, "status": doctor.OK, "tools": []},
    )
    rc = doctor.repair(
        ["read_qc", "assembly", "binning"],
        results(read_qc=doctor.OK, assembly=doctor.MISSING, binning=doctor.BROKEN),
        threads=6,
    )
    assert rc == 0
    argv = seen["argv"]
    assert argv[:2] == ["assembly", "binning"]
    assert "read_qc" not in argv
    assert "--force" in argv and "-t" in argv and "6" in argv and "--no-test" in argv
    assert "All rebuilt environments are healthy" in capsys.readouterr().out
