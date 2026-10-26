"""Env name resolution for install-env, and whether `check` looks where a run looks.

From user bug reports against 2.0.1. 2.1.0 grouped conda environments by workflow instead of one
per tool, so 2.0.x names such as ``annotate_bins-bakta`` and ``bin_refinement-checkm2`` no longer
existed:

* ``install-env bin_refinement-checkm2`` raised ``KeyError: 'bin_refinement-checkm2'`` in 2.0.1.
  The name must still work, since 2.0.x docs use it.
* A typo or obsolete name resolved to itself, so install printed "no env file" and the status table
  showed ``metawrap2-<typo>  MISSING  conda env not created``, as if a valid env had failed.
* ``check`` reported bakta present because it was on PATH, while the run looked inside a conda env
  that did not exist and failed with "The given prefix does not exist".
"""

from __future__ import annotations

import os

from metawrap2 import config
from metawrap2.commands import install_env

# ── what install-env accepts ─────────────────────────────────────────────────────────────────


def test_a_2_0_x_per_tool_env_name_still_resolves():
    """Env names used in 2.0.x docs and error messages still resolve.

    The first two were named in two bug reports from 2.0.1 users.
    """
    for legacy, expected in [
        ("annotate_bins-bakta", "metawrap2-annotate"),
        ("bin_refinement-checkm2", "metawrap2-checkm2"),
        ("bin_refinement-gunc", "metawrap2-binqc"),
        ("binning-semibin2", "metawrap2-learned"),
        ("classify_bins-gtdbtk", "metawrap2-gtdbtk"),
    ]:
        assert install_env._env_for(legacy) == expected, legacy
    # An alias must not point to another name that no longer exists.
    for legacy, env in install_env.LEGACY_ENV_ALIASES.items():
        assert os.path.isfile(install_env._env_yaml(env)), "%s -> %s has no yaml" % (legacy, env)


# ── every registry key resolves ─────────────────────────────────────────────────────────────


def test_every_key_in_every_registry_resolves_to_a_real_environment():
    """Every key in all four env registries resolves to a real environment.

    The registries are keyed differently: MODULE_ENVS by module, OPTIONAL_ENVS by tool (semibin2,
    vamb, comebin, amber), OPTIONAL_ENVS_BY_MODULE by env base name, and LEGACY_ENV_ALIASES by the
    2.0.x name. Consumers used to read different subsets: 2.0.1's install-env read only
    MODULE_ENVS, so three envs could not be installed, and an earlier fix read two of the four, so
    `install-env semibin2` was refused. A new registry that resolution does not read fails here.
    """

    unreachable = []
    for registry in config.ENV_REGISTRIES:
        for key in registry:
            if config.resolve_env_name(key) is None:
                unreachable.append(key)
    assert not unreachable, "names in a registry that resolution cannot reach: %r" % (unreachable,)
