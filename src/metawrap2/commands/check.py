"""`metawrap2 check [module ...]`: preflight tools, conda envs, Python extras and databases.

Reports what is missing before a run starts. Exits non-zero if a selected module is missing its
conda env, a tool inside that env, a required Python extra, or a database it needs. A database is
optional only when the module has a flag to run without it, and its label names that flag. With conda envs
enabled (the default), tools are looked up inside each module's env, not on the host PATH.

`check --for "<module command>"` reads one full module command line and reports only the
databases and conda envs that exact command needs, then ``ready: yes`` or ``ready: no``.
"""

from __future__ import annotations

import shutil
from typing import List, Optional

from .. import dbcheck as _dbcheck
from ..cliformat import make_parser
from ..config import (
    DB_KEYS,
    MODULE_TOOLS,
    OPTIONAL_ENV_OWNERS,
    OPTIONAL_ENVS,
    conda_env_exists,
    explain_unknown_name,
    legacy_module_for,
    load_settings,
    resolve_env_name,
    tool_in_env,
    tools_for_name,
)
from ..modules import binning as _binning

OK = "  OK   "
MISSING = "MISSING"


def _status(ok: bool) -> str:
    return OK if ok else MISSING


def _preset(argv: List[str]) -> str:
    """The binning preset *argv* selects, as its current name (``default`` when none is given)."""
    name = None
    for i, token in enumerate(argv):
        if token == "--preset" and i + 1 < len(argv):
            name = argv[i + 1]
        elif token.startswith("--preset="):
            name = token.split("=", 1)[1]
    resolved = _binning.resolve_preset(name) or "default"
    return resolved if resolved in _binning.PRESETS else "default"


def envs_for_command(module: str, argv: List[str]) -> List[str]:
    """The conda envs one invocation of *module* runs in: its own, and any opt-in env it selects.

    An opt-in env is selected by the flag of the same name (``classify_bins --gtdbtk``). CheckM2's
    env is needed exactly when the command reads CHECKM2DB, so its flags are taken from there.
    """
    envs: List[str] = []
    own = resolve_env_name(module)
    if own:
        envs.append(own)
    given = set(_dbcheck._flags(argv))
    for name, owner in OPTIONAL_ENV_OWNERS.items():
        if owner == module and "--" + name in given:
            envs.append(OPTIONAL_ENVS[name])
    if module == "binning" and not given & {"--" + b for b in _binning.BINNERS}:
        # No binner named, so the preset decides, and the default preset includes SemiBin2.
        envs += [OPTIONAL_ENVS[b] for b in _binning.PRESETS[_preset(argv)] if b in OPTIONAL_ENVS]
    if any(need.key == "CHECKM2DB" for need in _dbcheck.active_needs(module, argv)):
        envs.append(OPTIONAL_ENVS["checkm2"])
    return list(dict.fromkeys(envs))


def check_command(text: str, config: Optional[str] = None) -> int:
    """`check --for`: the databases and envs one module command needs, and whether it can run."""
    try:
        module, argv = _dbcheck.split_command(text)
    except ValueError as exc:
        print(str(exc))
        return 1
    settings = load_settings(config)
    missing = 0
    print("\n%s %s\n" % (module, " ".join(argv)) if argv else "\n%s\n" % module)
    entries = _dbcheck.resolve(settings, [(module, argv)])
    if not entries:
        print("  (no database needed)")
    for entry in entries:
        if entry.problem is None:
            print("  OK       %-19s %s" % (entry.key, entry.path))
            continue
        name = _dbcheck.install_name(entry.key)
        if entry.flag or not name:
            fix = entry.problem
        else:
            fix = "-> metawrap2 install-db %s" % name
        if entry.warn_only:
            # Not counted: the module runs without it.
            need = next(n for n in _dbcheck.module_needs(module) if n.key == entry.key)
            print("  NOT SET  %-19s %s (optional: %s)" % (entry.key, fix, need.note))
            continue
        missing += 1
        print("  MISSING  %-19s %s" % (entry.key, fix))
    own = resolve_env_name(module)
    for env in envs_for_command(module, argv):
        if not settings.use_conda_envs:
            print("  OK       %-19s (conda envs disabled; using PATH)" % env)
        elif conda_env_exists(env):
            print("  OK       %-19s conda env" % env)
        else:
            missing += 1
            target = module if env == own else env[len("metawrap2-") :]
            print("  MISSING  %-19s -> metawrap2 install-env %s" % (env, target))
    print("\nready: %s\n" % ("no (%d missing)" % missing if missing else "yes"))
    return 1 if missing else 0


def main(argv: List[str]) -> int:
    ap = make_parser(prog="metawrap2 check")
    ap.add_argument("modules", nargs="*", help="modules to check (default: all)")
    ap.add_argument("--config", help="path to metawrap2.toml")
    ap.add_argument(
        "--for",
        dest="command",
        metavar='"MODULE ARGS"',
        help="check only what one module command needs, given as one quoted command line "
        '(e.g. --for "annotate --kofam -b BINS -o OUT")',
    )
    args = ap.parse_args(argv)

    if args.command is not None:
        if args.modules:
            ap.error("give modules or --for, not both")
        return check_command(args.command, args.config)

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

    # Only the selected modules' databases are listed. One counts towards the verdict unless the
    # module has a flag to run without it; the label names that flag.
    shown_keys = _dbcheck.keys_for(modules)
    db_problems: List[str] = []
    print("\nDatabases:\n")
    if not shown_keys:
        print("  (none: %s reads no database)" % ", ".join(modules))
    for key, status, _detail in _dbcheck.check(settings, shown_keys):
        users = [m for m in modules if any(n.key == key for n in _dbcheck.module_needs(m))]
        label = _dbcheck.label_for(key, modules)
        required = not label
        label = "  " + label if label else ""
        shown = settings.db(key) or "<not set>"
        print(
            "  [%s]  %-19s %-40s used by %s%s"
            % (
                status.center(7) if status != _dbcheck.OK else OK,
                key,
                shown,
                ", ".join(users),
                label,
            )
        )
        if required and status != _dbcheck.OK:
            found = _dbcheck.problem(settings, key)
            if found:
                db_problems.append(found)
    others = len([key for key in DB_KEYS if key not in shown_keys])
    if others and args.modules:
        print(
            "\n  %d other database%s not needed by %s; 'metawrap2 doctor --all' lists them."
            % (others, " is" if others == 1 else "s are", ", ".join(modules))
        )

    if missing:
        print(
            "\n%d conda env(s) missing. Create them with:  metawrap2 install-env %s"
            % (missing, " ".join(modules))
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
    if missing or db_problems or missing_tools or missing_extras:
        failures = list(db_problems)
        if missing:
            failures.append("%d conda env(s) missing" % missing)
        if missing_tools:
            failures.append("%d tool(s) missing" % len(missing_tools))
        if missing_extras:
            failures.append("%d Python extra(s) missing" % len(missing_extras))
        print("\nPreflight FAILED for %s: %s\n" % (", ".join(modules), "; ".join(failures)))
        return 1
    print("\nPreflight OK for: %s\n" % ", ".join(modules))
    return 0
