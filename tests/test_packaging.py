"""Tests that the packaging metadata stays consistent with pyproject.toml.

Runtime requirements are declared in four places: `pyproject.toml`, the conda recipe
(`conda_pkg/meta.yaml`), the Docker image and the Apptainer definition. Each is read only by its own
tool, so nothing else catches drift. Three had fallen behind: they still required Python 3.8 after
the minimum became 3.10, and the conda recipe listed only biopython after the plotting stack became
a hard requirement.
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
    """CITATION.cff and the conda recipe state the same version as metawrap2.__version__."""
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
    """The pre-commit hook revs for black, isort and ruff match the dev-extra pins used by CI.

    A mismatch (ruff v0.6.9 against 0.16.8, black 24.8.0 against 26.5.1) caused two CI failures on
    code that passed the local hooks.
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


#: Imports the dependency check skips: the package itself, and tomllib (stdlib from 3.11) with its
#: 3.10 stand-in tomli.
_IMPORT_EXEMPT = ("metawrap2", "tomli", "tomllib")

#: Import name -> distribution name, where they differ.
_DIST_FOR_IMPORT = {
    "Bio": "biopython",
    "yaml": "pyyaml",
    "skbio": "scikit-bio",
    "sklearn": "scikit-learn",
}


def _module_level_imports(path: str) -> set:
    """Top-level names imported when *path* is imported, excluding imports inside functions."""
    import ast
    import sys

    with open(path) as fh:
        tree = ast.parse(fh.read())
    names = set()
    pending = list(tree.body)
    while pending:
        node = pending.pop()
        if isinstance(node, (ast.Try, ast.If)):
            # Guarded imports still run at import time; only the body counts, since the except
            # branch is the fallback.
            pending.extend(node.body)
            continue
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.add(node.module.split(".")[0])
    return {n for n in names if n not in sys.stdlib_module_names and n not in _IMPORT_EXEMPT}


def test_every_import_that_runs_at_import_time_is_a_declared_dependency(pyproject):
    """Every module-level third-party import is listed in [project].dependencies.

    Guards against a repeat of pyrodigal, imported at the top of genecall.py but not declared, which
    made `import metawrap2.genecall` fail on a clean install. Imports inside functions are treated
    as optional and not checked.
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
    """The conda recipe's run requirements match pyproject's dependencies."""
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
