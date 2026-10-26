"""`metawrap2 check [module ...]`: preflight tools, conda envs, Python extras and databases.

Reports what is missing before a run starts. Exits non-zero if a selected module is missing its
conda env, a tool inside that env, a required Python extra, or a database it needs. With conda envs
enabled (the default), tools are looked up inside each module's env, not on the host PATH.
"""

from __future__ import annotations

import shutil
from typing import List

from ..cliformat import make_parser
from ..config import (
    DB_KEYS,
    MODULE_TOOLS,
    OPTIONAL_DB_KEYS,
    check_databases,
    conda_env_exists,
    explain_unknown_name,
    legacy_module_for,
    load_settings,
    resolve_env_name,
    tool_in_env,
    tools_for_name,
)

OK = "  OK   "
MISSING = "MISSING"


def _status(ok: bool) -> str:
    return OK if ok else MISSING


def main(argv: List[str]) -> int:
    ap = make_parser(prog="metawrap2 check")
    ap.add_argument("modules", nargs="*", help="modules to check (default: all)")
    ap.add_argument("--config", help="path to metawrap2.toml")
    args = ap.parse_args(argv)

    settings = load_settings(args.config)
    modules = args.modules or sorted(MODULE_TOOLS)
    # Reject names that match neither a module nor an env, instead of reporting OK after checking
    # nothing. Modules without an env (`diff_abundance`, `quantify_function` run in this package's
    # interpreter) are accepted via MODULES.
    from ..cli import MODULES

    unknown = [name for name in modules if resolve_env_name(name) is None and name not in MODULES]
    if unknown:
        for name in unknown:
            for line in explain_unknown_name(name):
                print(line)
        return 1
    # A 2.0.x env name (<module>-<tool>) is translated to its module, which has the tool list.
    translated = []
    for name in modules:
        module = legacy_module_for(name)
        if module:
            print(
                "'%s' was an environment name in MetaWrap2 2.0.x; checking the '%s' module, "
                "which now provides that tool." % (name, module)
            )
            translated.append(module)
        else:
            translated.append(name)
    modules = translated

    missing = 0
    missing_tools: List[str] = []
    # Some modules need one of this package's Python extras (`diff_abundance` needs scikit-bio).
    # Uses the same probe as `doctor`.
    from .doctor import MODULE_EXTRAS, _check_extra

    extra_rows = [
        _check_extra(module, MODULE_EXTRAS[module]) for module in modules if module in MODULE_EXTRAS
    ]
    missing_extras = [row for row in extra_rows if row["status"] != "OK"]
    if extra_rows:
        print("\nPython extras:\n")
        for row in extra_rows:
            print(
                "  [%s]  %-18s %s" % (_status(row["status"] == "OK"), row["module"], row["detail"])
            )

    print("\nConda environments:\n")
    for module in modules:
        env = resolve_env_name(module)
        if not env:
            continue
        exists = conda_env_exists(env) if settings.use_conda_envs else True
        if not exists:
            missing += 1
        note = "" if settings.use_conda_envs else " (conda envs disabled; using PATH)"
        print("  [%s]  %-28s %s%s" % (_status(exists), env, module, note))

    where = "in the module's conda env" if settings.use_conda_envs else "on PATH"
    print("\nExternal tools (%s):\n" % where)
    checked_any = False
    for module in modules:
        env = resolve_env_name(module)
        for tool in tools_for_name(module):
            checked_any = True
            if settings.use_conda_envs and env:
                present = bool(tool_in_env(tool, env))
            else:
                present = shutil.which(tool) is not None
            print("  [%s]  %-18s (%s)" % (_status(present), tool, module))
            if not present:
                # Counted separately from missing envs: an env with a missing tool is half-built
                # and needs `install-env --force`.
                missing_tools.append("%s (%s)" % (tool, module))
    if not checked_any:
        # Say so when there was nothing to check.
        print("  (nothing to check: %s names no tools)" % ", ".join(modules))

    # Only databases read by a selected module count towards the verdict; optional databases
    # (Bakta, GTDB-Tk, ...) never do.
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
            "\nInstall them with `metawrap2 install-db` (it writes their paths into your "
            "config), or set the paths under [databases] in metawrap2.toml."
            % (len(missing_dbs), ", ".join(missing_dbs))
        )
    if missing_tools:
        print(
            "\n%d tool(s) are missing from an existing environment: %s"
            "\nThe environment is incomplete; rebuild it with "
            "`metawrap2 install-env %s --force`."
            % (len(missing_tools), ", ".join(missing_tools), " ".join(modules))
        )
    if missing_extras:
        print(
            "\n%d module(s) need a Python extra that is not installed: %s"
            % (len(missing_extras), ", ".join(row["module"] for row in missing_extras))
        )
        for row in missing_extras:
            print("  %s" % row["detail"])
    if missing or missing_dbs or missing_tools or missing_extras:
        print()
        return 1
    print("\nPreflight OK for: %s\n" % ", ".join(modules))
    return 0
