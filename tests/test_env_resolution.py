"""Which environment a name installs into, and whether `check` asks the question a run will ask.

Both of these came from user bug reports against 2.0.1, and they have one root cause: 2.1.0 grouped
the conda environments by workflow instead of one per tool, so the names 2.0.x documented -
``annotate_bins-bakta``, ``bin_refinement-checkm2`` - stopped existing. What the reports showed:

* ``install-env bin_refinement-checkm2`` raised ``KeyError: 'bin_refinement-checkm2'`` in 2.0.1. The
  traceback is gone now, but the name still has to *work*, because it is what 2.0.x's own docs say
  to type and what an upgrading user will reach for.
* A typo or an obsolete name used to resolve to itself, so the install printed "no env file" and
  then a status table reporting ``metawrap2-<typo>  MISSING  conda env not created`` - which reads
  as an environment that failed to build rather than a name that was never valid.
* ``check`` reported bakta present because it was on PATH, while the run looked for it inside a
  conda env that did not exist and died with "The given prefix does not exist". The check passed on
  a question nobody had asked.
"""

from __future__ import annotations

import os

from metawrap2 import config
from metawrap2.commands import install_env

# ── what install-env accepts ─────────────────────────────────────────────────────────────────


def test_a_2_0_x_per_tool_env_name_still_resolves():
    """These are what 2.0.x's docs and error messages name, so they have to keep working.

    Two separate bug reports from 2.0.1 users named the first two of these exactly.
    """
    for legacy, expected in [
        ("annotate_bins-bakta", "metawrap2-annotate"),
        ("bin_refinement-checkm2", "metawrap2-checkm2"),
        ("bin_refinement-gunc", "metawrap2-binqc"),
        ("binning-semibin2", "metawrap2-learned"),
        ("classify_bins-gtdbtk", "metawrap2-gtdbtk"),
    ]:
        assert install_env._env_for(legacy) == expected, legacy
    # An alias to a name that is also gone would just move the error.
    for legacy, env in install_env.LEGACY_ENV_ALIASES.items():
        assert os.path.isfile(install_env._env_yaml(env)), "%s -> %s has no yaml" % (legacy, env)


# ── the invariant that makes another instance of this bug impossible ─────────────────────────


def test_every_key_in_every_registry_resolves_to_a_real_environment():
    """The root cause, as an assertion.

    There are four registries mapping something onto a conda env, keyed differently: MODULE_ENVS by
    module, OPTIONAL_ENVS by *tool* (semibin2, vamb, comebin, amber), OPTIONAL_ENVS_BY_MODULE by env
    base name, and LEGACY_ENV_ALIASES by the 2.0.x name. Every consumer used to read a different
    subset, which is the whole bug: 2.0.1's install-env read only MODULE_ENVS, so three envs were
    uninstallable; my first fix read only two of the four, so `install-env semibin2` was refused for
    an env the running pipeline itself uses.

    Resolution now reads all of them, and this asserts there is nothing in any of them it cannot
    reach. A fifth registry would fail here the moment it is added, which is the point.
    """

    unreachable = []
    for registry in config.ENV_REGISTRIES:
        for key in registry:
            if config.resolve_env_name(key) is None:
                unreachable.append(key)
    assert not unreachable, "names in a registry that resolution cannot reach: %r" % (unreachable,)
