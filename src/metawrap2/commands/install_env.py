"""`metawrap2 install-env [module ...] [--all]` - create per-module conda envs with mamba.

MetaWrap2 manages everything through conda environments but drives them with **mamba**
(much faster solves). This installer makes sure mamba is available (bootstrapping it once
via conda if needed), creates each module's small env from ``envs/<module>.yaml``, then
verifies what it built with ``metawrap2 test``.
"""

from __future__ import annotations

import argparse
import os
import shutil
from typing import List

from ..command import run
from ..config import MODULE_ENVS

# envs/ lives at the repo/package-data root: src/metawrap2/commands -> ../../../envs
_ENVS_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "envs"))


def _env_yaml(module: str) -> str:
    return os.path.join(_ENVS_DIR, "%s.yaml" % module)


def ensure_mamba() -> bool:
    """Make sure `mamba` is on PATH; bootstrap it into the base env via conda if not.

    Returns True if mamba is available afterwards. This is the only place conda itself is
    used - everything else goes through mamba.
    """
    if shutil.which("mamba"):
        return True
    if not shutil.which("conda"):
        print("Neither mamba nor conda was found on PATH. Install Miniforge "
              "(https://github.com/conda-forge/miniforge) - it ships mamba - then re-run.")
        return False
    print("mamba not found; installing it into the base environment with conda (one-time)...")
    run(["conda", "install", "-n", "base", "-y", "-c", "conda-forge", "mamba"],
        tool="conda", check=False)
    return shutil.which("mamba") is not None


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(prog="metawrap2 install-env")
    ap.add_argument("modules", nargs="*", help="modules whose envs to create")
    ap.add_argument("--all", action="store_true", help="create envs for every module")
    ap.add_argument("--no-test", action="store_true", help="don't verify envs after creating")
    args = ap.parse_args(argv)

    modules = sorted(MODULE_ENVS) if args.all else args.modules
    if not modules:
        ap.error("specify one or more modules, or --all")

    if not ensure_mamba():
        return 1

    rc = 0
    for module in modules:
        yaml = _env_yaml(module)
        if not os.path.isfile(yaml):
            print("No env file for module '%s' (%s); skipping." % (module, yaml))
            rc = 1
            continue
        print("Creating %s from %s ..." % (MODULE_ENVS[module], yaml))
        # mamba env create; -y where supported. Fall back message handled by the runner.
        if run(["mamba", "env", "create", "-y", "-f", yaml], tool="mamba", check=False) != 0:
            rc = 1

    if not args.no_test:
        from .doctor import main as test_main
        print("\nVerifying the environment(s)...")
        if test_main(list(modules) + ["--tools-only"]) != 0:
            rc = 1
    return rc
