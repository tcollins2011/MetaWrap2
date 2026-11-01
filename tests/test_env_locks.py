"""Tests for the committed conda lockfiles and `metawrap2 doctor --fix`.

The package manager is stubbed. Checked: lockfile naming and content, refusal of a lockfile that
fell back to human-readable output, restoring from a lockfile, and that ``--fix`` rebuilds only the
broken environments.
"""

from __future__ import annotations

import os
import re

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
    """--from-lock passes the environment's lockfile to the build.

    Tested through the command, since --from-lock was once parsed but never used: naming the module
    (`binning`) must find the environment's lockfile (`bins.linux-64.lock`).
    """
    built = _queue_bins(envs_dir, monkeypatch, with_lock=True)
    assert ie.main(["binning", "--from-lock", "--no-test"]) == 0
    assert built == [("metawrap2-bins", ie.lock_path("bins", "linux-64"))]


def test_a_lockfile_is_named_after_the_environment_not_the_name_typed(envs_dir, monkeypatch):
    """`--lock binning` writes the lockfile under the env name (bins), which --from-lock reads."""
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: True)
    written = []
    monkeypatch.setattr(
        ie, "write_lockfile", lambda name, env, plat=None: written.append(name) or "x"
    )
    assert ie.main(["binning", "--lock"]) == 0
    assert written == ["bins"]


# --- committed lockfiles match their yaml -------------------------------------------------

_ENVS = os.path.join(os.path.dirname(__file__), "..", "src", "metawrap2", "envs")


def _yaml_packages(path):
    """Package names listed under ``dependencies:`` in an env yaml."""
    names, inside = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith("dependencies:"):
                inside = True
            elif inside and line.startswith("  - "):
                names.append(re.split(r"[\s<>=!]", line[4:].strip(), maxsplit=1)[0])
            elif inside and line.strip() and not line.startswith(" "):
                inside = False
    return names


@pytest.mark.parametrize(
    "lock",
    sorted(f for f in os.listdir(os.path.join(_ENVS, "locks")) if f.endswith(".lock")),
)
def test_a_committed_lockfile_has_every_package_its_yaml_lists(lock):
    """A tool added to a yaml must also be in its lockfile.

    flye was added to reads.yaml but not to its lockfile, so `install-env --from-lock` built an env
    without it and `doctor` reported assembly as broken.
    """
    env = lock.split(".")[0]
    with open(os.path.join(_ENVS, "locks", lock)) as fh:
        built = {
            re.match(r"(.+?)-\d", url.rsplit("/", 1)[-1]).group(1)
            for url in fh
            if url.startswith("http")
        }
    missing = [p for p in _yaml_packages(os.path.join(_ENVS, env + ".yaml")) if p not in built]
    assert not missing, "%s is missing %s; refresh it with `install-env %s --lock`" % (
        lock,
        ", ".join(missing),
        env,
    )


def _version_key(version):
    """A comparable key from a conda version: its leading numeric fields ("1.3.3e" -> (1, 3, 3))."""
    key = []
    for field in version.split("."):
        digits = re.match(r"\d+", field)
        if not digits:
            break
        key.append(int(digits.group()))
        if digits.group() != field:
            break
    return tuple(key)


def _satisfies(version, spec):
    """Whether *version* meets a comma-separated conda spec such as ``>=2.4.0,<2.5.0``."""
    ops = {
        ">=": lambda a, b: a >= b,
        "<=": lambda a, b: a <= b,
        ">": lambda a, b: a > b,
        "<": lambda a, b: a < b,
        "==": lambda a, b: a == b,
        "!=": lambda a, b: a != b,
    }
    for clause in spec.split(","):
        op, bound = re.match(r"(>=|<=|>|<|==|!=)?\s*(\S+)", clause.strip()).groups()
        have, want = _version_key(version), _version_key(bound.rstrip(".*"))
        if op is None:
            # A bare version or 1.2.*: the leading fields must match.
            if have[: len(want)] != want:
                return False
        else:
            width = max(len(have), len(want))
            have, want = have + (0,) * (width - len(have)), want + (0,) * (width - len(want))
            if not ops[op](have, want):
                return False
    return True


def test_the_version_check_reads_conda_specs():
    assert _satisfies("2.4.1", ">=2.4.0,<2.5.0")
    assert not _satisfies("2.6.1", ">=2.4.0,<2.5.0")
    assert _satisfies("1.26.4", "<2") and not _satisfies("2.0", "<2")
    assert _satisfies("3.12.14", ">=3.9,<3.13") and not _satisfies("3.14.7", ">=3.9,<3.13")
    assert _satisfies("1.3.3e", ">=1.3.3e") and _satisfies("3.0.6", "3.0.*")


def _lock_versions(lock):
    versions = {}
    with open(os.path.join(_ENVS, "locks", lock)) as fh:
        for url in fh:
            if url.startswith("http"):
                name, version = re.match(
                    r"(.+?)-(\d[^-]*)-[^-]+$", url.rsplit("/", 1)[-1].split("#")[0]
                ).groups()
                versions[name] = version
    return versions


def _yaml_specs(path):
    """``{package: spec}`` for the entries in an env yaml that carry a version spec."""
    specs, inside = {}, False
    with open(path) as fh:
        for line in fh:
            if line.startswith("dependencies:"):
                inside = True
            elif inside and line.startswith("  - "):
                parts = line[4:].split("#")[0].split(None, 1)
                if len(parts) == 2:
                    specs[parts[0]] = parts[1].strip()
            elif inside and line.strip() and not line.startswith(" "):
                inside = False
    return specs


# Lockfiles committed before their yaml gained the pins they break. Each must be regenerated on
# linux-64 with `install-env <env> --force` then `install-env <env> --lock`; then remove it here
# (strict xfail, so a fixed lockfile fails until it is removed).
_STALE_LOCKS: dict = {}


@pytest.mark.parametrize(
    "lock",
    [
        (
            pytest.param(f, marks=pytest.mark.xfail(reason=_STALE_LOCKS[f], strict=True))
            if f in _STALE_LOCKS
            else f
        )
        for f in sorted(f for f in os.listdir(os.path.join(_ENVS, "locks")) if f.endswith(".lock"))
    ],
)
def test_a_committed_lockfile_meets_the_version_pins_of_its_yaml(lock):
    """A pin added to a yaml (python <3.13 for GTDB-Tk) must hold in its lockfile too.

    Otherwise `install-env --from-lock` builds the versions the pin exists to keep out.
    """
    env = lock.split(".")[0]
    versions = _lock_versions(lock)
    broken = [
        "%s %s (lockfile has %s)" % (name, spec, versions[name])
        for name, spec in _yaml_specs(os.path.join(_ENVS, env + ".yaml")).items()
        if name in versions and not _satisfies(versions[name], spec)
    ]
    assert not broken, "%s breaks its yaml's pins: %s" % (lock, "; ".join(broken))
