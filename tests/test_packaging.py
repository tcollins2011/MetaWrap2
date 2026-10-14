"""Tests that the packaging metadata does not drift out of step with pyproject.toml.

MetaWrap2 declares its runtime requirements in four places that no tool cross-checks: the wheel
(`pyproject.toml`), the conda recipe (`conda_pkg/meta.yaml`), the Docker image, and the Apptainer
definition. Three of them had fallen behind - still asking for Python 3.8 after the floor moved to
3.10, and the conda recipe still listing only biopython after the plotting stack became a hard
requirement, which would have produced a conda package that installs and then fails on the first
figure.

None of that is caught by the test suite, CI, or a wheel build, because each file is only read by
the tool that consumes it. These tests read them.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List

import pytest

try:  # Python 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read(path: str) -> str:
    with open(os.path.join(REPO, path)) as fh:
        return fh.read()


@pytest.fixture(scope="module")
def pyproject() -> Dict:
    with open(os.path.join(REPO, "pyproject.toml"), "rb") as fh:
        return _toml.load(fh)


def requirement_names(pyproject: Dict) -> List[str]:
    """The distribution names in [project].dependencies, without version specifiers."""
    names = []
    for spec in pyproject["project"]["dependencies"]:
        name = re.split(r"[<>=!~;\[ ]", spec.strip(), maxsplit=1)[0]
        if name:
            names.append(name.lower())
    return names


# --- the version is stated in three places ------------------------------------------------


def test_the_version_is_the_same_in_every_file_that_states_it():
    """Four files carry it, and a release that bumps three of them ships a lie in the fourth.

    `metawrap2 -v`, the wheel's metadata, the conda recipe and CITATION.cff are all read by
    somebody - the last one by anybody citing this - so they are checked together rather than by
    whoever remembers at release time.
    """
    import re

    from metawrap2 import __version__

    for name, pattern in (
        ("CITATION.cff", r'^version:\s*"([^"]+)"'),
        (os.path.join("conda_pkg", "meta.yaml"), r'^\s*version:\s*"([^"]+)"'),
    ):
        with open(os.path.join(REPO, name)) as fh:
            found = re.search(pattern, fh.read(), re.MULTILINE)
        assert found, "%s states no version" % name
        assert found.group(1) == __version__, "%s says %s, not %s" % (
            name,
            found.group(1),
            __version__,
        )


# --- what ships inside the wheel -----------------------------------------------------------


def test_every_env_yaml_sits_inside_the_package():
    """Anything outside src/metawrap2/ is absent from an installed wheel."""
    assert not os.path.isdir(
        os.path.join(REPO, "envs")
    ), "envs/ at the repo root is not shipped in the wheel; it belongs in src/metawrap2/envs/"
    inside = os.path.join(REPO, "src", "metawrap2", "envs")
    assert os.path.isdir(inside)
    assert [f for f in os.listdir(inside) if f.endswith(".yaml")]


def test_the_pre_commit_hooks_run_the_same_formatter_versions_ci_does(pyproject):
    """Two CI failures came from this: the hook and CI disagreeing about what is correctly formatted.

    CI installs the dev extras and runs black/ruff/isort at the versions pinned there; the commit
    hook fetches its own copy at the `rev` in .pre-commit-config.yaml. When those differ - ruff
    v0.6.9 against 0.16.8, then black 24.8.0 against 26.5.1 - the hook rewrites a file one way,
    CI's check wants it the other, and a push that passed every local gate fails remotely.
    """
    import re

    pins = {}
    for spec in pyproject["project"]["optional-dependencies"]["dev"]:
        found = re.match(r"^\s*([A-Za-z0-9_.-]+)==([^;\s]+)", spec)
        if found:
            pins[found.group(1).lower()] = found.group(2)

    with open(os.path.join(REPO, ".pre-commit-config.yaml")) as fh:
        config = fh.read()
    hooks = {
        "black": r"psf/black\s*\n\s*rev:\s*v?(\S+)",
        "isort": r"pycqa/isort\s*\n\s*rev:\s*v?(\S+)",
        "ruff": r"ruff-pre-commit\s*\n\s*rev:\s*v?(\S+)",
    }
    for tool, pattern in hooks.items():
        found = re.search(pattern, config)
        assert found, "no %s hook found in .pre-commit-config.yaml" % tool
        assert tool in pins, "%s is not pinned with == in pyproject's dev extras" % tool
        assert (
            found.group(1) == pins[tool]
        ), "the %s pre-commit hook is %s but CI installs %s; set both to the same version" % (
            tool,
            found.group(1),
            pins[tool],
        )


#: Imports the dependency check skips: the package itself, and the one version-guarded pair -
#: tomllib is stdlib from 3.11, tomli its declared stand-in on 3.10.
_IMPORT_EXEMPT = ("metawrap2", "tomli", "tomllib")

#: Import name -> distribution name, where they differ.
_DIST_FOR_IMPORT = {
    "Bio": "biopython",
    "yaml": "pyyaml",
    "skbio": "scikit-bio",
    "sklearn": "scikit-learn",
}


def _module_level_imports(path: str) -> set:
    """Top-level names imported when *path* is imported - not those deferred into a function."""
    import ast
    import sys

    with open(path) as fh:
        tree = ast.parse(fh.read())
    names = set()
    pending = list(tree.body)
    while pending:
        node = pending.pop()
        if isinstance(node, (ast.Try, ast.If)):
            # A guarded import (try/except ImportError, or a version check) is still run at import.
            # The except branch is the fallback, so only the body counts as required.
            pending.extend(node.body)
            continue
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.add(node.module.split(".")[0])
    return {n for n in names if n not in sys.stdlib_module_names and n not in _IMPORT_EXEMPT}


def test_every_import_that_runs_at_import_time_is_a_declared_dependency(pyproject):
    """`import metawrap2.genecall` raised ModuleNotFoundError on a clean install.

    pyrodigal was imported at the top of genecall.py and declared nowhere, so `annotate` worked
    only on the machine where it had been installed by hand - and the test collecting genecall
    failed everywhere else. A dependency imported inside a function is a deliberate optional one
    and is not checked here; one imported at module level is not optional, whatever was meant.
    """
    declared = set(requirement_names(pyproject))
    missing = {}
    for root, _dirs, files in os.walk(os.path.join(REPO, "src", "metawrap2")):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(root, name)
            for module in _module_level_imports(path):
                dist = _DIST_FOR_IMPORT.get(module, module).lower()
                if dist not in declared:
                    missing.setdefault(dist, []).append(os.path.relpath(path, REPO))
    assert not missing, "imported at module level but not in [project].dependencies: %s" % missing


def test_the_conda_recipe_declares_what_pyproject_declares(pyproject):
    """meta.yaml says to keep the two lists in step; this is what makes that more than a comment."""
    pip_names = set(requirement_names(pyproject)) - {"tomli"}
    with open(os.path.join(REPO, "conda_pkg", "meta.yaml")) as fh:
        run = fh.read().split("  run:", 1)[1].split("\n\n", 1)[0]
    conda_names = {
        re.sub(r"-base$", "", m.group(1).lower())
        for m in re.finditer(r"^\s*-\s*([A-Za-z0-9_.-]+)", run, re.MULTILINE)
    } - {"python", "tomli"}
    assert pip_names == conda_names, "pyproject only: %s; meta.yaml only: %s" % (
        sorted(pip_names - conda_names),
        sorted(conda_names - pip_names),
    )
