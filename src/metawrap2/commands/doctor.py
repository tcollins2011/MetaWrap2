"""`metawrap2 test [module ...]` - check every module's conda env and software, run the
unit tests, and print a readable installed / missing / broken map.

For each module it: confirms the `metawrap2-<module>` conda env exists, then goes *into*
that env and probes each key tool (present on PATH? does it run, or is it broken by e.g. a
missing shared library / import error?). It also runs the Python unit-test suite (which
covers the env-independent core). The result is a per-module status table so a user can see
at a glance what is ready, what needs installing, and what is broken.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from typing import List, Optional, Tuple

from ..config import MODULE_ENVS, MODULE_TOOLS, conda_env_exists

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
    "fastqc": "fastqc --version",
    "metabat2": "metabat2 -h",
    "concoct": "concoct --version",
    "minimap2": "minimap2 --version",
    "checkm": "checkm -h",
}

# Status glyphs / labels.
OK, MISSING, BROKEN = "OK", "MISSING", "BROKEN"


def _env_run(env: str, shell_cmd: str, timeout: int = 60) -> Tuple[int, str]:
    """Run a shell command inside a conda env via mamba; return (returncode, output)."""
    cmd = ["mamba", "run", "--no-capture-output", "-n", env, "bash", "-lc", shell_cmd]
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, timeout=timeout)
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


def _check_module(module: str) -> dict:
    env = MODULE_ENVS.get(module, "metawrap2-%s" % module)
    if not conda_env_exists(env):
        return {"module": module, "env": env, "status": MISSING,
                "detail": "conda env not created", "tools": []}
    tools = []
    n_missing = n_broken = 0
    for tool in MODULE_TOOLS.get(module, []):
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
        p = subprocess.run([sys.executable, "-m", "pytest", "-q", tests_dir],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                           cwd=os.path.dirname(tests_dir))
    except FileNotFoundError:
        return None, "skipped (pytest not installed)"
    summary = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ""
    return p.returncode == 0, summary


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(prog="metawrap2 test")
    ap.add_argument("modules", nargs="*", help="modules to test (default: all)")
    ap.add_argument("--no-unit-tests", action="store_true", help="skip the Python unit tests")
    ap.add_argument("--tools-only", action="store_true", help="only probe tools/envs")
    ap.add_argument("--verbose", action="store_true", help="list every tool and its version")
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
    print("  %d OK, %d missing, %d broken (of %d modules)"
          % (counts[OK], counts[MISSING], counts[BROKEN], len(results)))

    unit_ok = None
    if not args.tools_only and not args.no_unit_tests:
        unit_ok, summary = _run_unit_tests()
        print("\n  Unit tests: %s" % summary)

    print("")
    # exit nonzero if anything is missing/broken or unit tests failed
    bad = counts[MISSING] + counts[BROKEN]
    if bad or unit_ok is False:
        return 1
    return 0
