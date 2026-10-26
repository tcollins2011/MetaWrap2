"""Tests that what the code loads at runtime is what the wheel ships.

Guards against a missing data file in the wheel: ``[tool.setuptools.package-data]`` once listed
only the env files, so ``data/kegg/modules_table.tsv`` was not shipped and `pathways` and
`quantify_function` raised FileNotFoundError on wheel installs. ``pip install -e .`` reads from the
checkout, so the test suite did not see it.

These tests read the package-data patterns from pyproject.toml, expand them against the source
tree as setuptools does, and compare the result with the files the package loads. The wheel itself
is built and installed into a clean venv by .github/workflows/wheel-install.yml.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path
from typing import Dict, List, Set

import pytest

try:  # Python 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml

REPO = Path(__file__).resolve().parent.parent
PACKAGE = REPO / "src" / "metawrap2"


@pytest.fixture(scope="module")
def packaged_files() -> Set[Path]:
    """Every file inside the package that the wheel would carry, relative to the package root.

    Expanded from the declared patterns with pathlib's recursive glob, as setuptools does.
    """
    with open(REPO / "pyproject.toml", "rb") as fh:
        pyproject = _toml.load(fh)
    patterns: List[str] = pyproject["tool"]["setuptools"]["package-data"]["metawrap2"]
    assert patterns, "pyproject declares no package-data for metawrap2"
    found: Set[Path] = set()
    for pattern in patterns:
        found.update(p.relative_to(PACKAGE) for p in PACKAGE.glob(pattern) if p.is_file())
    return found


def _non_python_files(directory: Path) -> Set[Path]:
    """Data files under *directory*, relative to the package root. Python is packaged already."""
    files = set()
    for root, dirs, names in os.walk(directory):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for name in names:
            if name.endswith((".py", ".pyc")):
                continue
            files.add((Path(root) / name).relative_to(PACKAGE))
    return files


def _literals_in_the_package() -> Set[str]:
    """Every string constant in the package's own source.

    Read with ast, not imported, because importing would run the helper scripts' module-level
    argv handling.
    """
    literals: Set[str] = set()
    for path in PACKAGE.rglob("*.py"):
        if "vendor" in path.parts:
            continue
        tree = ast.parse(path.read_text(errors="replace"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                literals.add(node.value)
    return literals


def _data_files_the_code_loads() -> Dict[Path, str]:
    """Files under ``data/`` whose name the package's source names, as {relative path: basename}.

    Derived from the source so it stays current. Data files not named in the source are also
    shipped (see the test below); a missing file from this set means a module cannot start.
    """
    literals = _literals_in_the_package()
    return {rel: rel.name for rel in _non_python_files(PACKAGE / "data") if rel.name in literals}


def test_every_data_file_in_the_tree_is_shipped(packaged_files):
    """Every file under data/ is matched by a package-data pattern, so new tables ship."""
    data_files = _non_python_files(PACKAGE / "data")
    assert data_files, "src/metawrap2/data/ is empty; this test is checking nothing"
    missing = sorted(str(f) for f in data_files - packaged_files)
    assert not missing, (
        "under src/metawrap2/data/ but not matched by any [tool.setuptools.package-data] "
        "pattern, so absent from the wheel: %s" % missing
    )


def test_every_data_file_the_code_opens_is_shipped(packaged_files):
    """Every data file named in the package source is in the wheel."""
    loaded = _data_files_the_code_loads()
    assert loaded, (
        "no file under src/metawrap2/data/ is named anywhere in the package source; either the "
        "data directory moved or this derivation stopped working - fix it rather than deleting it"
    )
    missing = sorted(str(f) for f in set(loaded) - packaged_files)
    assert not missing, (
        "loaded at runtime but not shipped, so `pip install metawrap2` gives a module that "
        "raises FileNotFoundError: %s" % missing
    )


def test_the_data_files_the_code_names_are_really_there():
    """The data file paths named by the code exist in the tree, inside the package."""
    from metawrap2 import kegg

    assert os.path.isfile(kegg.MODULES_TABLE), (
        "%s is named by metawrap2.kegg and is not in the checkout" % kegg.MODULES_TABLE
    )
    assert (
        Path(kegg.MODULES_TABLE).resolve().is_relative_to(PACKAGE.resolve())
    ), "the module table is loaded from outside the package, so no packaging change can ship it"


def test_the_env_specs_and_lockfiles_are_still_shipped(packaged_files):
    """Env specs and lockfiles are in the wheel.

    `metawrap2 install-env` reads them from the installed package; without the locks,
    `--from-lock` does not work.
    """
    env_files = _non_python_files(PACKAGE / "envs")
    assert env_files, "src/metawrap2/envs/ is empty; this test is checking nothing"
    missing = sorted(str(f) for f in env_files - packaged_files)
    assert not missing, "env specs/lockfiles absent from the wheel: %s" % missing
