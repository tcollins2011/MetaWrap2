"""Tests that what the code loads at runtime is what the wheel ships.

`pathways` and `quantify_function` were dead on every wheel install for as long as they existed:
``src/metawrap2/data/kegg/modules_table.tsv`` is in the repository, the module table is loaded on
the first line of real work, and ``[tool.setuptools.package-data]`` listed only the env files - so
the file was not in the wheel and ``load_modules`` raised FileNotFoundError at
``site-packages/metawrap2/data/kegg/modules_table.tsv``. Nothing caught it: it works under
``pip install -e .``, which reads straight out of the checkout and is the one mode a developer
uses, and the whole suite passed on a tree where two of twenty-four modules could not start.

These tests cross that boundary without building anything: they read the package-data patterns out
of pyproject.toml, expand them against the source tree the way setuptools does, and compare the
result with the files the package actually reaches for. The wheel itself is built and installed
into a clean venv by .github/workflows/wheel-install.yml, which is the only place the real
install layout can be checked.
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

    Expanded from the declared patterns with pathlib's recursive glob, which is what setuptools
    uses for these entries - so this is the pattern list's own answer, not a second opinion about
    which files ought to ship.
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

    Read with ast rather than imported: importing every module to see which paths it builds would
    run the helper scripts' module-level argv handling, and the question here is only which file
    names the source names.
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

    Derived, not listed: a hardcoded list is the same mistake as the hardcoded pattern list - it
    is correct on the day it is written and silent afterwards. A data file whose name appears
    nowhere in the source is still shipped (see the test below); it is this set that must never
    be missing from the wheel, because a missing one is a module that cannot start.
    """
    literals = _literals_in_the_package()
    return {rel: rel.name for rel in _non_python_files(PACKAGE / "data") if rel.name in literals}


def test_every_data_file_in_the_tree_is_shipped(packaged_files):
    """A new reference table under data/ must ship without anyone remembering pyproject.toml.

    This is the test that fails on the next modules_table.tsv: the KEGG table was added after the
    package-data list was written, the list was not extended, and two modules shipped broken.
    """
    data_files = _non_python_files(PACKAGE / "data")
    assert data_files, "src/metawrap2/data/ is empty; this test is checking nothing"
    missing = sorted(str(f) for f in data_files - packaged_files)
    assert not missing, (
        "under src/metawrap2/data/ but not matched by any [tool.setuptools.package-data] "
        "pattern, so absent from the wheel: %s" % missing
    )


def test_every_data_file_the_code_opens_is_shipped(packaged_files):
    """The failure the user hit: an installed module reaching for a file the wheel left behind."""
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
    """The derivation above proves nothing if the paths it checks do not exist in the tree."""
    from metawrap2 import kegg

    assert os.path.isfile(kegg.MODULES_TABLE), (
        "%s is named by metawrap2.kegg and is not in the checkout" % kegg.MODULES_TABLE
    )
    assert (
        Path(kegg.MODULES_TABLE).resolve().is_relative_to(PACKAGE.resolve())
    ), "the module table is loaded from outside the package, so no packaging change can ship it"


def test_the_env_specs_and_lockfiles_are_still_shipped(packaged_files):
    """The guarantee the package-data list was written for, kept while data/ was added to it.

    `metawrap2 install-env` reads these out of the installed package, so an install that has the
    yamls but not the locks silently loses `--from-lock`, which is the reproducible install.
    """
    env_files = _non_python_files(PACKAGE / "envs")
    assert env_files, "src/metawrap2/envs/ is empty; this test is checking nothing"
    missing = sorted(str(f) for f in env_files - packaged_files)
    assert not missing, "env specs/lockfiles absent from the wheel: %s" % missing
