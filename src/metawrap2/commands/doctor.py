"""`metawrap2 test [module ...]` - check every module's conda env and software, run the
unit tests, and print a readable installed / missing / broken map.

For each module it: confirms the module's conda env exists (several modules share one), then
goes *into* that env and probes each key tool (present on PATH? does it run, or is it broken by e.g. a
missing shared library / import error?). It also runs the Python unit-test suite (which
covers the env-independent core). The result is a per-module status table so a user can see
at a glance what is ready, what needs installing, and what is broken.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

from .. import dbcheck as _dbcheck
from ..cliformat import make_parser
from ..command import _run_flags
from ..config import (
    MODULE_ENVS,
    MODULE_TOOLS,
    OPTIONAL_ENVS_BY_MODULE,
    OPTIONAL_TOOLS_BY_MODULE,
    conda_env_exists,
)

# Signatures in a tool's output that mean it's installed but broken (not just "exits nonzero
# on --version", which many healthy tools do).
_BROKEN_SIGNATURES = re.compile(
    r"error while loading shared libraries|cannot open shared object|ImportError|"
    r"ModuleNotFoundError|Traceback \(most recent call last\)|GLIBC|symbol lookup error|"
    r"No such file or directory|command not found|Segmentation fault",
    re.IGNORECASE,
)

# Preferred health probe per tool (a command that a healthy install answers). Tools not
# listed fall back to trying --version, --help, -h in turn.
_PROBES = {
    "samtools": "samtools --version",
    "bwa": "bwa",  # prints usage on stderr, exit 1 when healthy
    "bowtie2": "bowtie2 --version",
    "spades.py": "spades.py --version",
    "metaspades.py": "metaspades.py --version",
    "megahit": "megahit --version",
    "quast": "quast --version",
    "salmon": "salmon --version",
    "kraken2": "kraken2 --version",
    "prokka": "prokka --version",
    "taxator": "taxator --help",
    "blastn": "blastn -version",
    "trim_galore": "trim_galore --version",
    "fastp": "fastp --version",
    "fastqc": "fastqc --version",
    "metabat2": "metabat2 -h",
    "concoct": "concoct --version",
    "minimap2": "minimap2 --version",
    "checkm": "checkm -h",
    "gunc": "gunc --version",
}

# Status glyphs / labels.
OK, MISSING, BROKEN = "OK", "MISSING", "BROKEN"


def _env_run(env: str, shell_cmd: str, timeout: int = 60) -> Tuple[int, str]:
    """Run a shell command inside a conda env via mamba; return (returncode, output).

    Uses ``bash -c``, not ``bash -lc``: a login shell re-reads the user's profile, which
    rebuilds PATH from scratch and throws away the environment's own bin directory - so every
    tool probe reported MISSING even for a perfectly good env. Flags for ``run`` come from
    :func:`metawrap2.command._run_flags` because ``--no-capture-output`` is conda-only.
    """
    cmd = ["mamba", "run"] + _run_flags("mamba") + ["-n", env, "bash", "-c", shell_cmd]
    try:
        p = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            check=False,
        )
        return p.returncode, p.stdout or ""
    except subprocess.TimeoutExpired:
        return 124, "timed out"
    except FileNotFoundError:
        return 127, "mamba not found"


def _probe_tool(env: str, tool: str) -> Tuple[str, str]:
    """Return (status, detail) for one tool inside *env*: OK / MISSING / BROKEN."""
    rc, _ = _env_run(env, "command -v %s" % tool)
    if rc != 0:
        return MISSING, "not on PATH"
    # present - now see if it actually runs
    probe = _PROBES.get(tool)
    probes = [probe] if probe else ["%s --version" % tool, "%s --help" % tool, "%s -h" % tool]
    last = ""
    for pr in probes:
        rc, out = _env_run(env, pr)
        last = out
        if rc == 0:
            version = out.strip().splitlines()[0][:60] if out.strip() else "runs"
            return OK, version
        if _BROKEN_SIGNATURES.search(out):
            return BROKEN, out.strip().splitlines()[-1][:80] if out.strip() else "runtime error"
    # present but no probe exited 0 and no breakage signature -> assume OK (nonzero --version)
    if _BROKEN_SIGNATURES.search(last):
        return BROKEN, "runtime error"
    return OK, "present"


def _tools_for(module: str, env: str) -> List[str]:
    """Every tool *module* should be able to run, for a module name or for a bare env name.

    Three ways in, because all three are names a user can type at `doctor` and `install-env`:

    * a module (``binning``) - that module's tools.
    * an opt-in env's yaml name (``learned``) - the tools that env is supposed to hold.
    * an env shared by several modules, named by the env (``bins``) - the union of their tools, since
      the env is only healthy if every module that will reach for it can run.

    The third case is why this is not a dictionary lookup. Several modules share an env now, so
    ``bins`` is not a key of MODULE_TOOLS at all, and looking it up directly returned nothing - which
    `doctor` then reported as "0/0 tools OK": a pass that had probed nothing. install-env runs exactly
    that check after building an env, so the one command whose job is to verify a build was the one
    reporting success without looking.
    """
    direct = MODULE_TOOLS.get(module) or OPTIONAL_TOOLS_BY_MODULE.get(module)
    if direct:
        return list(direct)
    shared: List[str] = []
    for name, module_env in MODULE_ENVS.items():
        if module_env != env:
            continue
        for tool in MODULE_TOOLS.get(name, []):
            if tool not in shared:
                shared.append(tool)
    return shared


def _check_module(module: str) -> dict:
    # An opt-in env is named after its yaml, so `doctor bin_refinement-gunc` (and the verification
    # install-env runs after building one) lands here with a name MODULE_ENVS does not hold.
    env = MODULE_ENVS.get(module) or OPTIONAL_ENVS_BY_MODULE.get(module, "metawrap2-%s" % module)
    if not conda_env_exists(env):
        return {
            "module": module,
            "env": env,
            "status": MISSING,
            "detail": "conda env not created",
            "tools": [],
        }
    tools = []
    n_missing = n_broken = 0
    for tool in _tools_for(module, env):
        status, detail = _probe_tool(env, tool)
        tools.append((tool, status, detail))
        if status == MISSING:
            n_missing += 1
        elif status == BROKEN:
            n_broken += 1
    if n_missing or n_broken:
        mod_status = BROKEN
        detail = "%d missing, %d broken of %d tools" % (n_missing, n_broken, len(tools))
    else:
        mod_status = OK
        detail = "%d/%d tools OK" % (len(tools), len(tools))
    return {"module": module, "env": env, "status": mod_status, "detail": detail, "tools": tools}


def _check_optional_envs(modules: List[str]) -> List[dict]:
    """The opt-in envs belonging to *modules*, probed if they exist.

    Reported but never counted towards the exit code: an absent opt-in env is a choice, not a
    fault. A *present but broken* one is worth seeing, which is why they are probed at all - GUNC
    runs by default in bin_refinement, so "the env is there and gunc does not run" is a failure a
    user should hear about before a four-hour refinement does.
    """
    from ..config import OPTIONAL_ENV_OWNERS, OPTIONAL_ENV_TOOLS, OPTIONAL_ENVS

    rows: List[dict] = []
    for flag, env in sorted(OPTIONAL_ENVS.items()):
        owner = OPTIONAL_ENV_OWNERS.get(flag)
        if owner is None or owner not in modules:
            continue
        if not conda_env_exists(env):
            rows.append(
                {
                    "flag": flag,
                    "env": env,
                    "status": MISSING,
                    "detail": "not created (opt-in)",
                    "tools": [],
                }
            )
            continue
        tools = [(t, *_probe_tool(env, t)) for t in OPTIONAL_ENV_TOOLS.get(flag, [flag])]
        unhealthy = [t for t in tools if t[1] != OK]
        rows.append(
            {
                "flag": flag,
                "env": env,
                "status": BROKEN if unhealthy else OK,
                "detail": (
                    "%d of %d tools not healthy" % (len(unhealthy), len(tools))
                    if unhealthy
                    else "%d/%d tools OK" % (len(tools), len(tools))
                ),
                "tools": tools,
            }
        )
    return rows


def _find_tests_dir() -> Optional[str]:
    """Locate the repo tests/ dir (present in a source/editable checkout)."""
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        cand = os.path.join(here, "tests")
        if os.path.isdir(cand) and any(f.startswith("test_") for f in os.listdir(cand)):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None


def _run_unit_tests() -> Tuple[Optional[bool], str]:
    tests_dir = _find_tests_dir()
    if not tests_dir:
        return None, "skipped (run from a source checkout to run unit tests)"
    try:
        p = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", tests_dir],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=os.path.dirname(tests_dir),
            check=False,
        )
    except FileNotFoundError:
        return None, "skipped (pytest not installed)"
    summary = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ""
    return p.returncode == 0, summary


def repair(
    modules: List[str], results: List[dict], threads: int = 4, use_lock: bool = False
) -> int:
    """Recreate the environments that are missing or broken. Returns an exit code.

    Recreating the environment is the only real remedy for a broken one - a half-solved env
    cannot be patched into a working one - so rather than telling the user to go and do it,
    ``--fix`` does exactly that, for exactly the environments that are actually broken. Healthy
    environments are never touched.
    """
    from .install_env import main as install_main

    damaged = [r["module"] for r in results if r["status"] in (MISSING, BROKEN)]
    if not damaged:
        print("\n  Nothing to fix: every environment checked is healthy.\n")
        return 0

    print("\n  Recreating %d environment(s): %s" % (len(damaged), ", ".join(damaged)))
    print("  (healthy environments are left alone)\n")
    argv = list(damaged) + ["--force", "-t", str(threads), "--no-test"]
    if use_lock:
        argv.append("--from-lock")
    rc = install_main(argv)

    print("\n  Re-checking what was rebuilt...\n")
    after = [_check_module(m) for m in damaged]
    still_bad = [r["module"] for r in after if r["status"] != OK]
    for r in after:
        print("  %-16s %-30s %s" % (r["module"], r["env"], r["status"]))
        for tool, st, detail in r["tools"]:
            if st != OK:
                print("      - %-18s %-8s %s" % (tool, st, detail))
    if still_bad:
        print(
            "\n  %d environment(s) are still not healthy: %s\n"
            "  Their build logs are the next thing to look at; the solver output usually "
            "names the conflict.\n" % (len(still_bad), ", ".join(still_bad))
        )
        return 1
    print("\n  All rebuilt environments are healthy.\n")
    return rc


#: Which database keys each module needs. An environment can be perfectly healthy and the module
#: still unable to run, because the database it reads is missing or half-downloaded - so checking
#: envs without checking databases answers only half of "can I run this?".
MODULE_DATABASES: Dict[str, List[str]] = {
    "read_qc": ["BMTAGGER_DB"],
    "kraken2": ["KRAKEN2_DB"],
    "blobology": ["BLASTDB", "TAXDUMP"],
    "classify_bins": ["BLASTDB", "TAXDUMP"],
    "bin_refinement": ["CHECKM_DB", "GUNC_DB"],
    "reassemble_bins": ["CHECKM_DB"],
    "binning": ["CHECKM_DB"],
    "annotate_bins": [],
    "assembly": [],
    "quant_bins": [],
}


def _report_databases(modules: List[str], config: Optional[str], verbose: bool = False) -> int:
    """Print the state of every database the given modules need. Returns how many are unusable.

    Optional databases (the opt-in CheckM2/GTDB-Tk/Bakta paths, and CheckM's data when only the
    binners are being used) are reported but not counted as failures, because a user who never
    asked for that code path has nothing wrong with their installation.
    """
    from ..config import OPTIONAL_DB_KEYS, load_settings

    needed: List[str] = []
    for module in modules:
        for key in MODULE_DATABASES.get(module, []):
            if key not in needed:
                needed.append(key)
    if not needed:
        return 0

    try:
        settings = load_settings(config)
    except Exception as exc:  # noqa: BLE001 - a broken config is itself the finding
        print("\nDatabases: could not read the config - %s" % exc)
        return 1

    rows = _dbcheck.check(settings, sorted(needed))
    print("\nDatabase status")
    print("=" * 74)
    print("  %-18s %-11s %s" % ("database", "status", "detail"))
    print("  " + "-" * 70)
    bad = 0
    for key, status, detail in rows:
        optional = key in OPTIONAL_DB_KEYS
        shown = status if status == _dbcheck.OK or not optional else "%s (optional)" % status
        print("  %-18s %-11s %s" % (key, shown, detail))
        if status != _dbcheck.OK and not optional:
            bad += 1
    print("  " + "-" * 70)
    print(
        "  %d usable, %d unusable. Install with: metawrap2 install-db --list"
        % (sum(1 for _k, st, _d in rows if st == _dbcheck.OK), bad)
    )
    if verbose:
        print(
            "\n  The value in brackets is the database's content fingerprint, also recorded in\n"
            "  each run's run_environment.json - so a database replaced in place is visible as a\n"
            "  change rather than an identical-looking path."
        )
    return bad


def main(argv: List[str]) -> int:
    ap = make_parser(prog="metawrap2 test")
    ap.add_argument("modules", nargs="*", help="modules to test (default: all)")
    ap.add_argument("--no-unit-tests", action="store_true", help="skip the Python unit tests")
    ap.add_argument("--tools-only", action="store_true", help="only probe tools/envs")
    ap.add_argument(
        "--no-databases", action="store_true", help="skip checking the configured databases"
    )
    ap.add_argument("--config", help="config file to read database paths from")
    ap.add_argument("--verbose", action="store_true", help="list every tool and its version")
    ap.add_argument(
        "--fix",
        action="store_true",
        help="recreate the environments reported missing or broken, then re-check "
        "them (healthy environments are left alone)",
    )
    ap.add_argument(
        "--from-lock",
        action="store_true",
        help="with --fix, rebuild from the committed lockfiles rather than by " "solving the yaml",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="with --fix, how many environments to rebuild at once (default 4)",
    )
    args = ap.parse_args(argv)

    modules = args.modules or sorted(MODULE_ENVS)
    results = [_check_module(m) for m in modules]

    print("\nMetaWrap2 module status")
    print("=" * 74)
    print("  %-16s %-30s %-8s %s" % ("module", "conda env", "status", "detail"))
    print("  " + "-" * 70)
    counts = {OK: 0, MISSING: 0, BROKEN: 0}
    for r in results:
        counts[r["status"]] += 1
        print("  %-16s %-30s %-8s %s" % (r["module"], r["env"], r["status"], r["detail"]))
        for tool, st, detail in r["tools"]:
            if st != OK or args.verbose:
                print("      - %-18s %-8s %s" % (tool, st, detail))
    print("  " + "-" * 70)
    print(
        "  %d OK, %d missing, %d broken (of %d modules)"
        % (counts[OK], counts[MISSING], counts[BROKEN], len(results))
    )

    optional = _check_optional_envs(modules)
    if optional:
        print("\nOpt-in environments (absent is a choice, not a fault)\n")
        for r in optional:
            print("  %-16s %-30s %-8s %s" % (r["flag"], r["env"], r["status"], r["detail"]))
            for tool, st, detail in r["tools"]:
                if st != OK or args.verbose:
                    print("      - %-18s %-8s %s" % (tool, st, detail))

    db_bad = 0
    if not args.tools_only and not args.no_databases:
        db_bad = _report_databases(modules, args.config, verbose=args.verbose)

    if args.fix:
        return repair(modules, results, threads=args.threads, use_lock=args.from_lock)

    unit_ok = None
    if not args.tools_only and not args.no_unit_tests:
        unit_ok, summary = _run_unit_tests()
        print("\n  Unit tests: %s" % summary)

    print()
    # exit nonzero if anything is missing/broken or unit tests failed
    bad = counts[MISSING] + counts[BROKEN] + db_bad
    if bad or unit_ok is False:
        return 1
    return 0
