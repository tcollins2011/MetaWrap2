"""A missing opt-in conda env is reported before the work, not at the point of use.

This is the difference between two failure modes rather than between failing and succeeding. Module
envs were always checked up front by env_for, so "metawrap2-bins does not exist" arrives in the
first second. Opt-in envs were looked up straight out of OPTIONAL_ENVS with no check, so the first
sign of trouble was mamba complaining about a missing prefix at whatever point the tool was
reached - after every read had been aligned for a learned binner, and after CheckM1 and
consolidation had both finished for GUNC.
"""

import pytest

from metawrap2.config import Settings
from metawrap2.modules import _common


@pytest.fixture
def env_manager_present(monkeypatch):
    """Satisfy the other precondition for running in an env, so only the env itself is under test.

    Running a tool in a conda env needs two things: the env, and the manager that enters it. Both are
    checked here, and a test about one of them should not depend on the developer's PATH for the
    other - mamba is typically not on PATH in a non-login shell, which is how pytest runs.
    """
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)


def test_a_missing_env_names_the_install_command_rather_than_the_prefix(
    monkeypatch, env_manager_present
):
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: False)
    said = []

    def fake_error(message, exit_code=1):
        said.append(message)
        raise SystemExit(exit_code)

    # The metawrap2 logger does not propagate to the root logger, so caplog cannot see it.
    monkeypatch.setattr(_common, "error", fake_error)
    with pytest.raises(SystemExit):
        _common.optional_env_for("comebin", Settings(), "--comebin")
    # The yaml name, not the env name and not the module name: that is what install-env takes.
    assert "install-env learned" in said[0]
    assert "--comebin needs it" in said[0]
