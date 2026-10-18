"""`metawrap2 test [module ...]` - check every module's conda env and software, run the
unit tests, and print a readable installed / missing / broken map.

For each module it: confirms the module's conda env exists (several modules share one), then
goes *into* that env and probes each key tool (present on PATH? does it run, or is it broken by e.g. a
missing shared library / import error?). It also runs the Python unit-test suite (which
covers the env-independent core). The result is a per-module status table so a user can see
at a glance what is ready, what needs installing, and what is broken.

Every module ``metawrap2 -h`` lists gets a row, including the two that own no conda env because
they run in this package's own interpreter. One of those is gated on a declared Python extra
(see :data:`MODULE_EXTRAS`), which is probed here too: a dependency this command does not look at
is a dependency it reports as satisfied.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml

from .. import dbcheck as _dbcheck
from ..cliformat import make_parser
from ..command import _run_flags
from ..config import (
    MODULE_ENVS,
    MODULE_TOOLS,
    OPTIONAL_TOOLS_BY_MODULE,
    active_conda_root,
    conda_env_exists,
    env_prefix,
    explain_unknown_name,
    legacy_module_for,
    prefix_in_active_root,
    resolve_env_name,
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

# What conda and libmamba say when the environment they were pointed at is not there. This is the
# message `doctor` used to throw away: the probe returned nonzero, the output was dropped, and the
# row read "MISSING  not on PATH" - a diagnosis of the wrong thing entirely, for an env that exists
# and holds the tool. One report cost the user 50 minutes before the real message was found by
# running mamba by hand.
_NO_SUCH_PREFIX = re.compile(
    r"given prefix does not exist|EnvironmentLocationNotFound|Not a conda environment|"
    r"EnvironmentNameNotFound|Could not find conda environment",
    re.IGNORECASE,
)
# The path libmamba quotes back: `critical libmamba The given prefix does not exist: "/x/envs/y"`.
_PREFIX_IN_MESSAGE = re.compile(r'does not exist:\s*"?([^"\n]+?)"?\s*$', re.MULTILINE)

# Status glyphs / labels.
OK, MISSING, BROKEN = "OK", "MISSING", "BROKEN"

#: env name -> its prefix on disk (or None), resolved once per process. Resolving walks
#: :func:`config.env_dirs`, which shells out to ``conda env list``; a full `doctor` probes ~90
#: tools with up to three commands each, so resolving per call meant hundreds of listings for an
#: answer that cannot change mid-run.
_PREFIX_CACHE: Dict[str, Optional[str]] = {}


def _resolve_prefix(env: str) -> Optional[str]:
    """Where *env* lives, or None if it exists nowhere. Resolved once and reused."""
    if env not in _PREFIX_CACHE:
        _PREFIX_CACHE[env] = env_prefix(env)
    return _PREFIX_CACHE[env]


def _env_run(
    env: str, shell_cmd: str, timeout: int = 60, prefix: Optional[str] = None
) -> Tuple[int, str]:
    """Run a shell command inside a conda env via mamba; return (returncode, output).

    Selects the env by ``-p <prefix>`` whenever a prefix resolves, and only falls back to
    ``-n <env>`` when none does. ``-n`` resolves names against the *active* conda root alone,
    while existence is answered across every installation on the machine (see
    :func:`config.conda_env_exists`), so for an env belonging to another conda install the two
    disagreed: the env was reported present and then every command against it died with
    ``critical libmamba The given prefix does not exist``. With ``-p`` there is one answer, and it
    is the path that was checked for existence in the first place.

    Uses ``bash -c``, not ``bash -lc``: a login shell re-reads the user's profile, which
    rebuilds PATH from scratch and throws away the environment's own bin directory - so every
    tool probe reported MISSING even for a perfectly good env. Flags for ``run`` come from
    :func:`metawrap2.command._run_flags` because ``--no-capture-output`` is conda-only.
    """
    if prefix is None:
        prefix = _resolve_prefix(env)
    target = ["-p", prefix] if prefix else ["-n", env]
    cmd = ["mamba", "run"] + _run_flags("mamba") + target + ["bash", "-c", shell_cmd]
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


def _one_line(out: str, width: int = 60) -> str:
    """The last non-empty line of *out*, trimmed - for putting subprocess output in a table cell."""
    lines = [line.strip() for line in out.strip().splitlines() if line.strip()]
    return lines[-1][:width] if lines else ""


def _probe_tool(env: str, tool: str) -> Tuple[str, str]:
    """Return (status, detail) for one tool inside *env*: OK / MISSING / BROKEN.

    A failed probe keeps its output. The first version of this discarded it and reported
    ``MISSING  not on PATH`` for every nonzero return, which turned "the prefix I was pointed at
    does not exist" - a message that names the problem exactly - into a statement about PATH that
    was simply false. The two cases are now distinguished and the captured line is carried into
    the detail, because a user reading the table should not have to re-run mamba by hand to find
    out what it actually said.
    """
    prefix = _resolve_prefix(env)
    rc, out = _env_run(env, "command -v %s" % tool, prefix=prefix)
    if rc != 0:
        if _NO_SUCH_PREFIX.search(out):
            named = _PREFIX_IN_MESSAGE.search(out)
            path = named.group(1) if named else (prefix or env)
            return MISSING, "env prefix does not exist: %s" % path
        detail = _one_line(out)
        return MISSING, ("not on PATH (%s)" % detail if detail else "not on PATH")
    # present - now see if it actually runs
    probe = _PROBES.get(tool)
    probes = [probe] if probe else ["%s --version" % tool, "%s --help" % tool, "%s -h" % tool]
    last = ""
    for pr in probes:
        rc, out = _env_run(env, pr, prefix=prefix)
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


def _foreign_prefix_note(env: str) -> str:
    """A note naming *env*'s prefix when it belongs to a different conda installation, else "".

    Existence is answered machine-wide, so the env `doctor` is reporting on may be one some other
    conda install built - same name, someone else's directory, possibly someone else's package
    versions. Probing it by prefix makes it work, which is right, but silently working on an env
    the user did not build is its own kind of wrong answer: if `install-env` writes to the active
    root and `doctor` reads a different one, the two commands discuss different environments and
    agree about nothing. Naming the path is what lets the user see that.
    """
    prefix = _resolve_prefix(env)
    if not prefix or prefix_in_active_root(prefix):
        return ""
    return "env is outside the active conda root (%s, root is %s)" % (prefix, active_conda_root())


def _check_module(module: str) -> dict:
    # An opt-in env is named after its yaml, so `doctor bin_refinement-gunc` (and the verification
    # install-env runs after building one) lands here with a name MODULE_ENVS does not hold.
    # config.resolve_env_name, not a fallback to "metawrap2-<whatever was typed>". That fallback is
    # why `doctor bin_refinement-checkm2` used to invent metawrap2-bin_refinement-checkm2 and report
    # it MISSING - an environment that never existed, presented as one that failed to install.
    # No "metawrap2-<whatever was typed>" fallback. That fallback is why this command used to
    # invent an environment for an unknown name and report it MISSING; main() refuses such names
    # before reaching here, and anything that still gets through should be visibly wrong rather
    # than plausibly wrong.
    env = resolve_env_name(module) or "<unknown: %s>" % module
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
    foreign = _foreign_prefix_note(env)
    if foreign:
        detail = "%s; %s" % (detail, foreign)
    return {"module": module, "env": env, "status": mod_status, "detail": detail, "tools": tools}


#: The distribution the extras below belong to, spelled the way pip names it.
DISTRIBUTION = "metawrap2"

#: Modules gated on a declared Python extra rather than on a conda env, and the extra each one
#: needs (a key of ``[project.optional-dependencies]``). What the extra *contains* is deliberately
#: not written here: the package names and version floors are read back out of the distribution's
#: own metadata, so this table and ``pip install 'metawrap2[stats]'`` cannot disagree about what
#: that installs.
#:
#: The defect this exists for: `doctor` probed conda envs and databases only, so `diff_abundance` -
#: whose single dependency is this extra - was neither OK nor broken but *absent from the table*,
#: and a fresh install reported "22 OK, 0 missing, 0 broken (of 22 modules)" to a user whose next
#: command died for want of scikit-bio. The module list is now the CLI's own, so a row count that
#: disagrees with `metawrap2 -h` is no longer possible; this table only says which extra gates what.
MODULE_EXTRAS: Dict[str, str] = {"diff_abundance": "stats"}

#: What the env column says for a module that has no conda env. Two labels rather than one blank,
#: because "nothing to install" and "install it with pip" need different things from the reader.
PYTHON_EXTRA = "(python extra)"
PURE_PYTHON = "(pure Python)"

#: How long to wait for the import probe below. Importing scikit-bio pulls in numpy, scipy and
#: pandas, which is seconds on a cold page cache and must not be mistaken for a hang.
_IMPORT_TIMEOUT = 120


def _module_names() -> List[str]:
    """Every module `metawrap2 -h` lists, in the order it lists them.

    Read off the CLI's own table rather than from ``MODULE_ENVS``, which is keyed by the modules
    that own a conda env and therefore holds 22 of the 24. A list derived from `-h` cannot drift
    from `-h`; a second list written out by hand is exactly how this one drifted.
    """
    from ..cli import MODULES

    return list(MODULES)


def _canonical(name: str) -> str:
    """A distribution name in its comparable form (PEP 503): ``scikit_bio`` == ``scikit-bio``."""
    return re.sub(r"[-_.]+", "-", name).lower()


def _requirement_name(requirement: str) -> str:
    """The distribution name out of a requirement string: ``scikit-bio>=0.7.1`` -> ``scikit-bio``."""
    for index, char in enumerate(requirement):
        if char in "<>=!~ [(;":
            return requirement[:index].strip()
    return requirement.strip()


def _version_floor(requirement: str) -> str:
    """The ``>=`` floor in *requirement*, or "" when it has none."""
    _name, separator, rest = requirement.partition(">=")
    return rest.split(",")[0].strip() if separator else ""


def _release(version: str) -> Optional[Tuple[int, ...]]:
    """*version*'s numeric release fields, or None if it is not plainly numeric."""
    fields = version.split(".")
    if not all(field.isdigit() for field in fields):
        return None
    return tuple(int(field) for field in fields)


def _version_at_least(installed: str, floor: str) -> bool:
    """Is *installed* at least *floor*?

    Not a PEP 440 implementation and deliberately not one: ``packaging`` is not a dependency of
    this package, and every floor these extras carry is a plain dotted number. A version on either
    side that does not parse counts as satisfying the floor, because reporting a module as unusable
    over a version string we failed to read would be a worse answer than the missing extra this
    check exists to find.
    """
    here, wanted = _release(installed), _release(floor)
    if here is None or wanted is None:
        return True
    return here >= wanted


def _checkout_root() -> Optional[str]:
    """The repository root when running from a source checkout, else None."""
    tests_dir = _find_tests_dir()
    return os.path.dirname(tests_dir) if tests_dir else None


def _pyproject_extras() -> Dict[str, List[str]]:
    """``[project.optional-dependencies]`` out of this checkout's pyproject.toml, or {}."""
    root = _checkout_root()
    if not root:
        return {}
    path = os.path.join(root, "pyproject.toml")
    try:
        with open(path, "rb") as handle:
            table = _toml.load(handle).get("project", {}).get("optional-dependencies", {})
    except (OSError, _toml.TOMLDecodeError):
        return {}
    return {name: list(requirements) for name, requirements in table.items()}


def _extra_requirements(extra: str) -> List[str]:
    """The requirement strings declared under *extra*, e.g. ``["scikit-bio>=0.7.1"]``.

    Taken from the installed distribution's metadata, which is what pip resolved and therefore the
    only description of the extra that matches this machine. The checkout's pyproject.toml is
    read as a fallback because an editable install's dist-info is only regenerated on reinstall: in
    a development tree it can predate the extra entirely and declare no extras at all, and a
    `doctor` that cannot see an extra reports the module as ready - the defect, again, one layer
    down.
    """
    from importlib import metadata

    marker = 'extra == "%s"' % extra
    try:
        declared = metadata.requires(DISTRIBUTION) or []
    except metadata.PackageNotFoundError:
        declared = []
    found = [
        requirement.split(";", 1)[0].strip()
        for requirement in declared
        if marker in requirement.replace("'", '"')
    ]
    return found or _pyproject_extras().get(extra, [])


def _installed_version(dist: str) -> Optional[str]:
    """The installed version of *dist*, or None if it is not installed."""
    from importlib import metadata

    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return None


def _import_names(dist: str) -> List[str]:
    """The top-level modules *dist* installs, read from the metadata rather than guessed.

    scikit-bio installs ``skbio`` and biom-format installs ``biom``: the import name is not the
    distribution name, and there is no rule that turns one into the other, so it is looked up.
    """
    from importlib import metadata

    wanted = _canonical(dist)
    return sorted(
        name
        for name, dists in metadata.packages_distributions().items()
        if any(_canonical(provider) == wanted for provider in dists)
    )


def _import_probe(dist: str) -> Tuple[int, str]:
    """Import what *dist* installs, in a subprocess; (returncode, output).

    A subprocess because a half-installed extra is one of the things this command exists to find,
    and an extra that fails at import - scikit-bio built against another numpy's ABI is the usual
    way - would otherwise raise out of `doctor` itself and report nothing about anything. (0, "")
    when the metadata names no importable module, which is not evidence either way.
    """
    names = _import_names(dist)
    if not names:
        return 0, ""
    source = "; ".join("import %s" % name for name in names)
    try:
        probe = subprocess.run(
            [sys.executable, "-c", source],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=_IMPORT_TIMEOUT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return 124, "timed out importing %s" % ", ".join(names)
    return probe.returncode, probe.stdout or ""


def _probe_requirement(requirement: str) -> Tuple[str, str]:
    """Return (status, detail) for one requirement out of an extra: OK / MISSING / BROKEN."""
    dist = _requirement_name(requirement)
    installed = _installed_version(dist)
    if installed is None:
        return MISSING, "not installed"
    floor = _version_floor(requirement)
    if floor and not _version_at_least(installed, floor):
        return MISSING, "%s %s is installed, %s needed" % (dist, installed, requirement)
    rc, out = _import_probe(dist)
    if rc != 0:
        return BROKEN, _one_line(out, 80) or "%s does not import" % dist
    return OK, "%s %s" % (dist, installed)


def _check_extra(module: str, extra: str) -> dict:
    """One table row for a module whose dependency is a Python extra rather than a conda env."""
    install = "pip install '%s[%s]'" % (DISTRIBUTION, extra)
    requirements = _extra_requirements(extra)
    if not requirements:
        # Nothing on this machine says what the extra holds, so nothing has been checked - and a
        # row that reports OK having checked nothing is the shape of the defect this command keeps
        # being bitten by (see `_tools_for`). Reported as missing, with the way out.
        return {
            "module": module,
            "env": PYTHON_EXTRA,
            "status": MISSING,
            "detail": "cannot tell what the %s extra holds; %s installs it" % (extra, install),
            "tools": [],
            "conda": False,
        }
    probed = [(requirement, *_probe_requirement(requirement)) for requirement in requirements]
    unmet = [probe for probe in probed if probe[1] != OK]
    if unmet:
        status = BROKEN if any(probe[1] == BROKEN for probe in probed) else MISSING
        detail = "needs %s: %s" % (", ".join(probe[0] for probe in unmet), install)
    else:
        status = OK
        detail = ", ".join(probe[2] for probe in probed)
    return {
        "module": module,
        "env": PYTHON_EXTRA,
        "status": status,
        "detail": detail,
        "tools": probed,
        "conda": False,
    }


def _check_any(module: str) -> dict:
    """The table row for *module*, whichever kind of dependency it has.

    Three kinds, and until now only the first had a row at all:

    * a conda env - 22 of the 24 modules.
    * a declared Python extra - `diff_abundance` needs scikit-bio and nothing else.
    * nothing - `quantify_function` runs in this package's own interpreter.

    The last two were indistinguishable from each other *and* from a module that does not exist,
    because neither appeared in the table. ``conda`` marks the rows ``--fix`` can act on: handing
    `diff_abundance` to `install-env` would ask it to build an environment that has never existed.
    """
    extra = MODULE_EXTRAS.get(module)
    if extra:
        return _check_extra(module, extra)
    if resolve_env_name(module) is None:
        return {
            "module": module,
            "env": PURE_PYTHON,
            "status": OK,
            "detail": "runs in MetaWrap2's own Python; nothing to install",
            "tools": [],
            "conda": False,
        }
    row = _check_module(module)
    row["conda"] = True
    return row


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
        detail = (
            "%d of %d tools not healthy" % (len(unhealthy), len(tools))
            if unhealthy
            else "%d/%d tools OK" % (len(tools), len(tools))
        )
        foreign = _foreign_prefix_note(env)
        if foreign:
            detail = "%s; %s" % (detail, foreign)
        rows.append(
            {
                "flag": flag,
                "env": env,
                "status": BROKEN if unhealthy else OK,
                "detail": detail,
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

    damaged = [r["module"] for r in results if r["status"] in (MISSING, BROKEN) and r["conda"]]
    # Rows with no conda env to rebuild: a module gated on a Python extra is fixed with pip, not by
    # `install-env`, which would be asked to build an environment that has never existed. They are
    # still part of the verdict, or `--fix` would report "nothing to fix" and exit 0 over a table
    # that says MISSING - the same contradiction between a verdict and its own output that this
    # command keeps being reported for.
    unfixable = [r for r in results if r["status"] in (MISSING, BROKEN) and not r["conda"]]
    for r in unfixable:
        print(
            "\n  %s cannot be fixed by rebuilding an environment: %s" % (r["module"], r["detail"])
        )
    if not damaged:
        if not unfixable:
            print("\n  Nothing to fix: every environment checked is healthy.\n")
        return 1 if unfixable else 0

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
        print("  %-18s %-30s %s" % (r["module"], r["env"], r["status"]))
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
    return rc or (1 if unfixable else 0)


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

    # Every module `metawrap2 -h` offers, not every module that owns a conda env: the two differ by
    # the modules that run in this package's own Python, one of which is gated on a Python extra
    # that nothing used to check.
    modules = args.modules or _module_names()
    known = set(_module_names())
    unknown = [name for name in modules if resolve_env_name(name) is None and name not in known]
    if unknown:
        for name in unknown:
            for line in explain_unknown_name(name):
                print(line)
        return 1
    # A 2.0.x name is <module>-<tool>; report on the module it came from, which has a tool list, in
    # place of the bare environment, which does not and therefore reports "0/0 tools OK".
    for name in modules:
        module = legacy_module_for(name)
        if module:
            print(
                "'%s' was an environment name in MetaWrap2 2.0.x; reporting on the '%s' module, "
                "which is where that tool lives now." % (name, module)
            )
    modules = [legacy_module_for(name) or name for name in modules]
    results = [_check_any(m) for m in modules]

    print("\nMetaWrap2 module status")
    print("=" * 74)
    # 18 for the name column: `quantify_function` is 17 characters and ran straight into the env
    # column at 16, which is also the width `metawrap2 -h` had to widen for the same name.
    print("  %-18s %-30s %-8s %s" % ("module", "dependency", "status", "detail"))
    print("  " + "-" * 70)
    counts = {OK: 0, MISSING: 0, BROKEN: 0}
    for r in results:
        counts[r["status"]] += 1
        print("  %-18s %-30s %-8s %s" % (r["module"], r["env"], r["status"], r["detail"]))
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
            print("  %-18s %-30s %-8s %s" % (r["flag"], r["env"], r["status"], r["detail"]))
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
