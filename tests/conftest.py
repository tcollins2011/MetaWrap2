import os
import sys

import pytest

# Allow tests to import the package from src/ without an install step.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


@pytest.fixture(autouse=True)
def _conda_env_policy(request, tmp_path_factory, monkeypatch):
    """Pin ``use_conda_envs`` for every test instead of inheriting the machine's config.

    Most tests fake the tool call and want no env wrapping at all; a handful are *about* how an
    env is resolved or required, and mark themselves ``@pytest.mark.conda_envs``. Left to the
    ambient config, both groups pass on one developer's machine and fail on another's: whichever
    way ``~/.metawrap2/config.toml`` or ``$METAWRAP2_CONFIG`` happens to be set decides.
    """
    wanted = request.node.get_closest_marker("conda_envs") is not None
    path = tmp_path_factory.mktemp("config") / "metawrap2.toml"
    path.write_text("[settings]\nuse_conda_envs = %s\n" % ("true" if wanted else "false"))
    monkeypatch.setenv("METAWRAP2_CONFIG", str(path))
