"""A missing opt-in conda env is reported before the work starts.

Module envs are checked up front by env_for. Opt-in envs used to be looked up in OPTIONAL_ENVS
without a check, so a missing one surfaced as a mamba error about a missing prefix when the tool
was reached: after all reads were aligned for a learned binner, or after CheckM1 and
consolidation for GUNC.
"""

import pytest

from metawrap2.config import Settings
from metawrap2.modules import _common


@pytest.fixture
def env_manager_present(monkeypatch):
    """Stub the env manager check, so only the env itself is under test.

    mamba is often not on PATH in a non-login shell, which is how pytest runs.
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
    # The yaml name (what install-env takes), not the env or module name.
    assert "install-env learned" in said[0]
    assert "--comebin needs it" in said[0]
