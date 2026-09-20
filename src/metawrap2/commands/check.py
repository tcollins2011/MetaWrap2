"""`metawrap2 check [module ...]` — preflight tools, conda envs, and databases.

Reports what's missing before a run starts, and exits non-zero if a selected module is
missing its conda env or a database it cannot run without.

The "External tools" section is informational only: with conda envs enabled (the default)
each module's tools live inside its env and are *expected* to be absent from the host PATH,
so their status cannot be part of the verdict. Use ``metawrap2 test`` to probe the tools
where they actually live, inside each env.
"""

from __future__ import annotations

import argparse
from typing import List

from ..config import (
    DB_KEYS,
    MODULE_ENVS,
    MODULE_TOOLS,
    OPTIONAL_DB_KEYS,
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

    print("\nExternal tools (on PATH -- informational; see `metawrap2 test`):\n")
    for tool, module, present in check_tools(modules):
        # A tool absent from PATH is fine (expected, even) if it lives in the module's env.
        print("  [%s]  %-18s (%s)" % (_status(present), tool, module))

    # Only databases a *selected* module actually reads count towards the verdict, and
    # opt-in-path databases (Bakta, GTDB-Tk, ...) never do.
    needed = {
        key
        for key, used_by in DB_KEYS.items()
        if key not in OPTIONAL_DB_KEYS and any(m in used_by for m in modules)
    }
    missing_dbs = []
    print("\nDatabases:\n")
    for key, path, present, used_by in check_databases(settings):
        shown = path if path else "<not set>"
        required = key in needed
        note = "" if required else "  (optional for the selected modules)"
        print("  [%s]  %-16s %-40s used by %s%s" % (_status(present), key, shown, used_by, note))
        if required and not present:
            missing_dbs.append(key)

    if missing:
        print(
            "\n%d conda env(s) missing. Create them with:  metawrap2 install-env %s"
            % (missing, " ".join(modules))
        )
    if missing_dbs:
        print(
            "\n%d database(s) the selected modules need are not set or not present: %s"
            "\nSet them under [databases] in metawrap2.toml (see "
            "installation/database_installation.md)." % (len(missing_dbs), ", ".join(missing_dbs))
        )
    if missing or missing_dbs:
        print()
        return 1
    print("\nPreflight OK for: %s\n" % ", ".join(modules))
    return 0
