"""`metawrap2 check [module ...]` — preflight tools, conda envs, and databases.

Reports what's missing before a run starts. Exits non-zero if anything a selected module
needs is unavailable.
"""

from __future__ import annotations

import argparse
from typing import List

from ..config import (
    MODULE_ENVS,
    MODULE_TOOLS,
    check_databases,
    check_tools,
    conda_env_exists,
    load_settings,
)

OK = "  OK   "
MISSING = "MISSING"


def _status(ok: bool) -> str:
    return OK if ok else MISSING


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(prog="metawrap2 check")
    ap.add_argument("modules", nargs="*", help="modules to check (default: all)")
    ap.add_argument("--config", help="path to metawrap2.toml")
    args = ap.parse_args(argv)

    settings = load_settings(args.config)
    modules = args.modules or sorted(MODULE_TOOLS)

    missing = 0

    print("\nConda environments:\n")
    for module in modules:
        env = MODULE_ENVS.get(module)
        if not env:
            continue
        exists = conda_env_exists(env) if settings.use_conda_envs else True
        if not exists:
            missing += 1
        note = "" if settings.use_conda_envs else " (conda envs disabled; using PATH)"
        print("  [%s]  %-28s %s%s" % (_status(exists), env, module, note))

    print("\nExternal tools (on PATH):\n")
    for tool, module, present in check_tools(modules):
        # A tool absent from PATH is fine if it lives in the module's conda env.
        print("  [%s]  %-18s (%s)" % (_status(present), tool, module))

    print("\nDatabases:\n")
    for key, path, present, used_by in check_databases(settings):
        shown = path if path else "<not set>"
        print("  [%s]  %-12s %-40s used by %s" % (_status(present), key, shown, used_by))

    if missing:
        print(
            "\n%d conda env(s) missing. Create them with:  metawrap2 install-env %s\n"
            % (missing, " ".join(modules))
        )
        return 1
    print("\nPreflight OK for: %s\n" % ", ".join(modules))
    return 0
