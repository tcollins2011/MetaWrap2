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

#: Files that also state a runtime requirement, and must therefore agree with pyproject.toml.
PACKAGING_FILES = (
    "conda_pkg/meta.yaml",
    "containers/Dockerfile",
    "containers/metawrap2.def",
)

#: Requirements whose name differs between PyPI and conda-forge. matplotlib is the one that
#: matters: the conda package `matplotlib` pulls in a Qt stack nobody needs in a container.
CONDA_NAME = {"matplotlib": "matplotlib-base"}

#: Dependencies that are conditional and so need not appear everywhere.
OPTIONAL_EVERYWHERE = {"tomli"}


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


def test_the_package_version_matches_the_module(pyproject):
    from metawrap2 import __version__

    assert pyproject["project"]["version"] == __version__


def test_the_conda_recipe_version_matches(pyproject):
    """A recipe left on an old version publishes the new code under the old number."""
    recipe = read("conda_pkg/meta.yaml")
    match = re.search(r'^\s*version:\s*"?([^"\s]+)"?', recipe, re.MULTILINE)
    assert match, "conda_pkg/meta.yaml states no version"
    assert match.group(1) == pyproject["project"]["version"]


# --- the Python floor is stated in four places ---------------------------------------------


def test_pyproject_and_mypy_target_the_same_python(pyproject):
    """mypy checking a newer Python than the package claims to support proves nothing."""
    declared = pyproject["project"]["requires-python"]
    target = pyproject["tool"]["mypy"]["python_version"]
    assert target in declared, "mypy targets %s but the package requires %s" % (target, declared)


def test_ruff_and_black_target_the_same_python(pyproject):
    minimum = re.search(r"(\d+)\.(\d+)", pyproject["project"]["requires-python"])
    assert minimum
    expected = "py%s%s" % (minimum.group(1), minimum.group(2))
    assert pyproject["tool"]["ruff"]["target-version"] == expected
    assert pyproject["tool"]["black"]["target-version"] == [expected]


@pytest.mark.parametrize("path", PACKAGING_FILES)
def test_no_packaging_file_asks_for_an_unsupported_python(path, pyproject):
    minimum = re.search(r"(\d+)\.(\d+)", pyproject["project"]["requires-python"])
    assert minimum
    floor = (int(minimum.group(1)), int(minimum.group(2)))
    text = read(path)
    # The lookbehind matters: without it this matches "biopython>=1.79" too.
    stated = re.findall(r"(?<![\w-])python\b\s*>=\s*(\d+)\.(\d+)", text)
    assert stated, "%s states no Python requirement" % path
    for major, minor in stated:
        found = (int(major), int(minor))
        assert found >= floor, "%s asks for python>=%d.%d but the package requires >=%d.%d" % (
            path,
            found[0],
            found[1],
            floor[0],
            floor[1],
        )


def test_ci_tests_the_python_versions_the_package_claims(pyproject):
    workflow = read(".github/workflows/ci.yml")
    minimum = re.search(r"(\d+)\.(\d+)", pyproject["project"]["requires-python"])
    assert minimum
    floor = "%s.%s" % (minimum.group(1), minimum.group(2))
    matrix = re.search(r"python-version:\s*\[([^\]]+)\]", workflow)
    assert matrix, "no python-version matrix in ci.yml"
    versions = re.findall(r"(\d+\.\d+)", matrix.group(1))
    assert floor in versions, "CI does not test the minimum supported Python (%s)" % floor
    # And it must not claim to test one the package refuses to install on.
    for version in versions:
        assert tuple(int(p) for p in version.split(".")) >= tuple(
            int(p) for p in floor.split(".")
        ), "CI tests Python %s but the package requires >=%s" % (version, floor)


# --- the runtime dependencies ---------------------------------------------------------------


def test_the_conda_recipe_lists_every_runtime_dependency(pyproject):
    """A conda package missing seaborn installs cleanly and then fails on the first figure."""
    recipe = read("conda_pkg/meta.yaml").lower()
    missing = []
    for name in requirement_names(pyproject):
        if name in OPTIONAL_EVERYWHERE:
            continue
        if CONDA_NAME.get(name, name) not in recipe:
            missing.append(name)
    assert missing == [], "conda_pkg/meta.yaml does not require: %s" % ", ".join(missing)


@pytest.mark.parametrize("path", ("containers/Dockerfile", "containers/metawrap2.def"))
def test_the_containers_install_every_runtime_dependency(path, pyproject):
    text = read(path).lower()
    missing = []
    for name in requirement_names(pyproject):
        if name in OPTIONAL_EVERYWHERE:
            continue
        if CONDA_NAME.get(name, name) not in text:
            missing.append(name)
    assert missing == [], "%s does not install: %s" % (path, ", ".join(missing))


def test_matplotlib_comes_from_matplotlib_base_in_conda(pyproject):
    """The plain conda `matplotlib` drags in a Qt stack a headless pipeline never uses."""
    for path in PACKAGING_FILES:
        text = read(path)
        if "matplotlib" in text:
            assert "matplotlib-base" in text, "%s should use matplotlib-base" % path


# --- what ships inside the wheel -----------------------------------------------------------


def test_the_env_specs_and_lockfiles_are_declared_as_package_data(pyproject):
    """They live under src/metawrap2/ so that `pip install metawrap2` can find them."""
    patterns = pyproject["tool"]["setuptools"]["package-data"]["metawrap2"]
    assert any(p.endswith("envs/*.yaml") for p in patterns)
    assert any("locks" in p for p in patterns)


def test_every_env_yaml_sits_inside_the_package():
    """Anything outside src/metawrap2/ is absent from an installed wheel."""
    assert not os.path.isdir(
        os.path.join(REPO, "envs")
    ), "envs/ at the repo root is not shipped in the wheel; it belongs in src/metawrap2/envs/"
    inside = os.path.join(REPO, "src", "metawrap2", "envs")
    assert os.path.isdir(inside)
    assert [f for f in os.listdir(inside) if f.endswith(".yaml")]


def test_the_console_script_points_at_a_real_entry_point(pyproject):
    target = pyproject["project"]["scripts"]["metawrap2"]
    module_path, function = target.split(":")
    import importlib

    assert callable(getattr(importlib.import_module(module_path), function))


def test_project_urls_are_declared(pyproject):
    urls = pyproject["project"]["urls"]
    assert "Homepage" in urls and "Repository" in urls
    assert all(v.startswith("https://") for v in urls.values())
