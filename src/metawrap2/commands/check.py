"""`metawrap2 check [module ...]` — preflight tools, conda envs, and databases.

Reports what's missing before a run starts, and exits non-zero if a selected module is
missing its conda env or a database it cannot run without.

The "External tools" section is informational only: with conda envs enabled (the default)
each module's tools live inside its env and are *expected* to be absent from the host PATH,
so their status cannot be part of the verdict. Use ``metawrap2 test`` to probe the tools
where they actually live, inside each env.
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
    # Refused rather than checked. `check bin_refinement-checkm2` used to print "Preflight OK for:
    # bin_refinement-checkm2" - it matched no module, so it checked nothing and declared success,
    # which is the worst of the three possible answers.
    # A module with no conda environment is still a module. `diff_abundance` and
    # `quantify_function` run in this package's own interpreter, so `resolve_env_name` returns None
    # for them - and this guard then refused two real modules as "not a module or an environment",
    # which is the same blind spot that kept them out of the readiness table.
    from ..cli import MODULES

    unknown = [name for name in modules if resolve_env_name(name) is None and name not in MODULES]
    if unknown:
        for name in unknown:
            for line in explain_unknown_name(name):
                print(line)
        return 1
    # A legacy or env-level name resolves to an env but is not a module, so it has no tool list and
    # no databases; report what can be reported rather than inventing a module row for it.
    # A 2.0.x name is <module>-<tool>, so check the module it came from rather than reporting OK
    # about an environment with no tool list attached.
    translated = []
    for name in modules:
        module = legacy_module_for(name)
        if module:
            print(
                "'%s' was an environment name in MetaWrap2 2.0.x; checking the '%s' module, which "
                "is where that tool lives now." % (name, module)
            )
            translated.append(module)
        else:
            translated.append(name)
    modules = translated

    missing = 0
    missing_tools: List[str] = []
    # A module can also be gated on one of this package's own Python extras, and until 2.1.3
    # neither `check` nor `doctor` looked at those - so `diff_abundance`, whose only dependency is
    # scikit-bio, passed a preflight it was not in and then failed on its first real invocation.
    # `doctor` owns the probe; this reuses it rather than growing a second, differently-wrong copy.
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
                # Counted separately from the env rows: a present env with a missing tool is a
                # half-built environment, and telling the user "1 conda env missing" sends them to
                # install-env when what they need is --force.
                missing_tools.append("%s (%s)" % (tool, module))
    if not checked_any:
        # Never silently. `check <something with no tool list>` used to print "Preflight OK" having
        # looked at nothing, which is a worse answer than either yes or no.
        print("  (nothing to check: %s names no tools)" % ", ".join(modules))

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
            "\nInstall them with `metawrap2 install-db`, which writes their paths into your "
            "config, or set them yourself under [databases] in metawrap2.toml."
            % (len(missing_dbs), ", ".join(missing_dbs))
        )
    if missing_tools:
        print(
            "\n%d tool(s) are absent from an environment that exists: %s"
            "\nThat is a half-built environment rather than a missing one; rebuild it with "
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
