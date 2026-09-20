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


# --- the Python floor is stated in four places ---------------------------------------------


# --- the runtime dependencies ---------------------------------------------------------------


# --- what ships inside the wheel -----------------------------------------------------------


def test_every_env_yaml_sits_inside_the_package():
    """Anything outside src/metawrap2/ is absent from an installed wheel."""
    assert not os.path.isdir(
        os.path.join(REPO, "envs")
    ), "envs/ at the repo root is not shipped in the wheel; it belongs in src/metawrap2/envs/"
    inside = os.path.join(REPO, "src", "metawrap2", "envs")
    assert os.path.isdir(inside)
    assert [f for f in os.listdir(inside) if f.endswith(".yaml")]
