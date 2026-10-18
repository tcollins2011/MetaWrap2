"""`metawrap2 install-env [module ...] [--all]` - create per-module conda envs with mamba.

MetaWrap2 manages everything through conda environments but drives them with **mamba**
(much faster solves). This installer makes sure mamba is available (bootstrapping it once
via conda if needed), creates each module's small env from ``envs/<module>.yaml``, then
verifies what it built with ``metawrap2 test``.

Environments are built several at a time (``-t``), since each one solves and downloads
independently. Each env's solver output goes to its own log file rather than the screen,
because concurrent solver output interleaved into one stream is unreadable - but each
environment does announce when it starts and how long it took, so the slowest step of an
install is not half an hour of silence.
"""

from __future__ import annotations

import concurrent.futures
import os
import platform
import shutil
import subprocess
import tempfile
import time
from typing import List, Optional, Tuple

from ..cliformat import make_parser
from ..command import run, runner
from ..config import (
    LEGACY_ENV_ALIASES,
    MODULE_ENVS,
    conda_env_exists,
    env_prefix,
    explain_unknown_name,
    known_envs,
    prefix_in_active_root,
    resolve_env_name,
)

# Env creation is dominated by solving and downloading, both of which are per-env
# independent, so building several at once is close to a linear speedup. Default is modest
# so a laptop on a slow link is not swamped; raise it with -t.
DEFAULT_PARALLEL_ENVS = 4


def format_duration(seconds: float) -> str:
    """``45s``, ``3m12s``, ``1h04m`` - short enough to sit at the end of a progress line."""
    seconds = max(0.0, seconds)
    if seconds < 60:
        return "%.0fs" % seconds
    minutes, secs = divmod(int(seconds), 60)
    if minutes < 60:
        return "%dm%02ds" % (minutes, secs)
    hours, minutes = divmod(minutes, 60)
    return "%dh%02dm" % (hours, minutes)


def envs_dir() -> str:
    """Locate the ``envs/`` directory of yaml environment specs.

    ``$METAWRAP2_ENVS_DIR`` wins, so a user can point at their own edited specs. Otherwise
    walk up from this file: that finds ``envs/`` both in a source checkout and in an
    editable install, and also covers ``envs/`` shipped alongside the installed package.
    """
    override = os.environ.get("METAWRAP2_ENVS_DIR")
    if override:
        return override
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        cand = os.path.join(here, "envs")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    # Last resort: the historical repo-relative guess, so the error names a sensible path.
    return os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "envs"))


def _env_yaml(env: str) -> str:
    """Locate the yaml that builds conda environment *env*.

    The yaml is named after the environment, not after the module: every env is
    ``metawrap2-<basename>`` and every spec is ``<basename>.yaml``. Resolving it this way is what
    lets several modules share one environment - ``coverage`` maps to ``metawrap2-bins`` because it
    runs the same aligner and samtools binning does, and there is no ``coverage.yaml`` to write.
    """
    base = env.removeprefix("metawrap2-")
    return os.path.join(envs_dir(), "%s.yaml" % base)


def _env_for(name: str) -> Optional[str]:
    """config.resolve_env_name, kept as a name because tests and callers here use it.

    Deliberately a one-line delegation rather than its own logic. This file had its own copy, and
    `check` and `doctor` each had another; all three disagreed, in three different ways.
    """
    return resolve_env_name(name)


def _resolve_or_explain(name: str) -> Optional[str]:
    """_env_for, but printing what to type instead when the name is not one."""
    env = _env_for(name)
    if env:
        if name in LEGACY_ENV_ALIASES:
            print(
                "'%s' was an environment name in MetaWrap2 2.0.x. Environments are now grouped by "
                "workflow, so it is part of '%s' - installing that." % (name, env)
            )
        return env
    for line in explain_unknown_name(name):
        print(line)
    return None


def ensure_mamba() -> bool:
    """Make sure `mamba` is on PATH; bootstrap it into the base env via conda if not.

    Returns True if mamba is available afterwards. This is the only place conda itself is
    used - everything else goes through mamba.
    """
    if shutil.which("mamba"):
        return True
    if not shutil.which("conda"):
        print(
            "Neither mamba nor conda was found on PATH. Install Miniforge "
            "(https://github.com/conda-forge/miniforge) - it ships mamba - then re-run."
        )
        return False
    print("mamba not found; installing it into the base environment with conda (one-time)...")
    run(
        ["conda", "install", "-n", "base", "-y", "-c", "conda-forge", "mamba"],
        tool="conda",
        check=False,
    )
    return shutil.which("mamba") is not None


#: Lockfiles are per-platform: a build hash that exists for linux-64 does not for osx-arm64.
def lock_platform() -> str:
    """The conda platform tag for this machine, used in lockfile names."""
    machine = platform.machine().lower()
    system = platform.system().lower()
    arch = {"x86_64": "64", "amd64": "64", "aarch64": "aarch64", "arm64": "arm64"}.get(
        machine, machine
    )
    prefix = {"linux": "linux", "darwin": "osx", "windows": "win"}.get(system, system)
    return "%s-%s" % (prefix, arch)


def lock_path(module: str, plat: Optional[str] = None) -> str:
    """Where a module's lockfile lives: ``envs/locks/<module>.<platform>.lock``."""
    return os.path.join(envs_dir(), "locks", "%s.%s.lock" % (module, plat or lock_platform()))


def write_lockfile(module: str, env: str, plat: Optional[str] = None) -> Optional[str]:
    """Write an exact-build lockfile for an existing env. Returns the path, or None.

    Uses ``conda list --explicit --md5``, which emits the URL and checksum of every package
    actually installed. That file recreates the environment bit-for-bit with
    ``conda create --file``, needs no extra tooling (conda-lock is not required), and is
    readable in a diff - so a change in a pinned build shows up in review.
    """
    manager = shutil.which("conda") or shutil.which("mamba") or shutil.which("micromamba")
    if not manager:
        return None
    try:
        out = subprocess.run(
            [manager, "list", "-n", env, "--explicit", "--md5"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=300,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0 or "@EXPLICIT" not in (out.stdout or ""):
        return None
    path = lock_path(module, plat)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    header = (
        "# MetaWrap2 lockfile for the %s module (%s).\n"
        "# Exact builds of everything in the %s environment, so a run from today can be\n"
        "# reproduced later. Generated with:  metawrap2 install-env %s --lock\n"
        "# Recreate with:  conda create -n %s --file %s\n"
        % (module, plat or lock_platform(), env, module, env, os.path.basename(path))
    )
    with open(path, "w") as fh:
        fh.write(header)
        fh.write(out.stdout)
    return path


def _create_one(
    module: str, env: str, yaml: str, force: bool, log_dir: str, lock: Optional[str] = None
) -> Tuple[str, bool, str]:
    """Create one env, capturing its output to a log file. Returns (module, ok, log path).

    Output goes to a file rather than the screen because several of these run at once and
    interleaved solver output is unreadable. The log path is reported on failure.
    """
    log = os.path.join(log_dir, "install-env-%s.log" % module)
    # Recreating means removing first. `mamba env create --force` is gone in mamba 2.x - the command
    # is micromamba's `create`, which rejects the flag outright ("The following argument was not
    # expected: --force"), so --force failed every rebuild instead of doing one. Only reached when
    # the user asked for a rebuild of an env MetaWrap2 itself built.
    if lock:
        # An explicit lockfile names every package build, so this is an install, not a solve.
        steps = [["mamba", "create", "-y", "-n", env, "--file", lock]]
    else:
        steps = [["mamba", "env", "create", "-y", "-f", yaml]]
    if force:
        steps.insert(0, ["mamba", "env", "remove", "-y", "-n", env])
    rc = 0
    with open(log, "w") as fh:
        for cmd in steps:
            fh.write("$ %s\n\n" % " ".join(cmd))
            fh.flush()
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
            if rc != 0 and cmd[2] == "create":
                break
    return module, rc == 0, log


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 install-env",
        description="Create the per-module conda environments, several at a time.",
    )
    ap.add_argument("modules", nargs="*", help="modules whose envs to create")
    ap.add_argument(
        "--all",
        action="store_true",
        help="create every environment, including the opt-in ones (checkm2, gtdbtk, the learned "
        "binners). Several are large",
    )
    ap.add_argument("--no-test", action="store_true", help="don't verify envs after creating")
    ap.add_argument(
        "--force",
        action="store_true",
        help="recreate envs that already exist (default: leave them alone)",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_PARALLEL_ENVS,
        help="how many environments to build concurrently (default %d). Each "
        "env solves and downloads independently, so this is close to a "
        "linear speedup; use 1 for serial, readable output." % DEFAULT_PARALLEL_ENVS,
    )
    ap.add_argument(
        "--lock",
        action="store_true",
        help="write an exact-build lockfile for each named module's existing "
        "environment into envs/locks/, instead of creating anything",
    )
    ap.add_argument(
        "--from-lock",
        action="store_true",
        help="create the environment(s) from their committed lockfiles rather "
        "than solving the yaml, for a bit-for-bit reproducible install",
    )
    ap.add_argument(
        "--log-dir",
        default=None,
        help="where to write per-env build logs (default: a temp directory)",
    )
    args = ap.parse_args(argv)
    # `--force` is also a global flag, and the CLI strips global flags wherever they appear - so
    # `metawrap2 install-env x --force` never reaches this parser and the flag documented above did
    # nothing at all. Read both: whichever way it was spelled, the user asked for a rebuild.
    force = args.force or runner.force

    # By env name, not by module name: several modules share an env now, so iterating MODULE_ENVS
    # would ask mamba to build metawrap2-reads four times. The names are interchangeable here -
    # _module_env falls through to the name itself for anything that is not a module.
    # --all means every environment, including the opt-in ones. It used to mean every *module* env,
    # so `install-env --all` left checkm2, gtdbtk and the learned binners uninstalled - and the
    # learned env holds SemiBin2, which the shipped default binning tier needs. A flag called --all
    # that silently skips four environments is the kind of thing a user only discovers when a run
    # fails halfway.
    if args.all:
        modules = [env.removeprefix("metawrap2-") for env in sorted(known_envs())]
    else:
        modules = args.modules
    if not modules:
        ap.error("specify one or more modules, or --all")

    if args.lock:
        # Locking reads existing environments; it does not need mamba bootstrapped.
        rc = 0
        for module in modules:
            env = _resolve_or_explain(module)
            if env is None:
                rc = 1
                continue
            if not conda_env_exists(env):
                print("%s does not exist yet; create it before locking it." % env)
                rc = 1
                continue
            # Named after the environment, not the name typed: `--lock binning` and `--lock bins`
            # describe one environment and must write one file, the one --from-lock reads.
            path = write_lockfile(env.removeprefix("metawrap2-"), env)
            if path:
                print("wrote %s" % path)
            else:
                print("could not lock %s (is conda on PATH?)" % env)
                rc = 1
        return rc

    if not ensure_mamba():
        return 1

    rc = 0
    todo = []
    # Every env this run has already decided about, whether it was queued for building or skipped
    # because it exists. It used to hold only the queued ones, because the "already exists" branch
    # below `continue`d before reaching the `add` - so a run over modules that share an env printed
    # "metawrap2-reads already exists" five times and "metawrap2-bins already exists" four, while
    # the build path printed each env once. Same decision, two different-looking outputs.
    seen_envs = set()
    for module in modules:
        env = _resolve_or_explain(module)
        if env is None:
            rc = 1
            continue
        yaml = _env_yaml(env)
        if not os.path.isfile(yaml):
            print("No env file for module '%s' (%s); skipping." % (module, yaml))
            rc = 1
            continue
        if env in seen_envs:
            # Two modules sharing one env: `--all` reaches both, and building it twice at once would
            # have two mamba solves writing the same prefix.
            print("%s is already covered by this run; skipping (%s shares it)." % (env, module))
            continue
        seen_envs.add(env)
        if conda_env_exists(env) and not force:
            # Name the prefix, not just the env. `conda_env_exists` answers across every conda
            # installation on the machine, so "already exists" on its own can mean an env that
            # belongs to a different conda and is not the one this user is about to run against.
            # The path is how they tell; without it the only way to find out was to go looking.
            where = env_prefix(env)
            print(
                "%s already exists at %s; skipping (use --force to recreate)."
                % (env, where or "a location conda did not report")
            )
            if where and not prefix_in_active_root(where):
                print(
                    "    note: that prefix is not under the conda installation on your PATH, so "
                    "this environment was built by a different conda."
                )
            continue
        lock = None
        if args.from_lock:
            candidate = lock_path(env.removeprefix("metawrap2-"))
            if os.path.isfile(candidate):
                lock = candidate
            else:
                print(
                    "No %s lockfile for %s; solving %s instead, so this one is not an exact rebuild."
                    % (lock_platform(), env, os.path.basename(yaml))
                )
        todo.append((module, env, yaml, lock))

    if todo and runner.dry_run:
        # The global --dry-run promises "print the commands without running them", and this command
        # used to ignore it: a dry run solved and built every environment it named.
        for _module, env, yaml, lock in todo:
            if force:
                print("[dry-run] mamba env remove -y -n %s" % env)
            if lock:
                print("[dry-run] mamba create -y -n %s --file %s" % (env, lock))
            else:
                print("[dry-run] mamba env create -y -f %s" % yaml)
        return rc

    if not todo:
        print("Nothing to create.")
    else:
        log_dir = args.log_dir or tempfile.mkdtemp(prefix="metawrap2-install-env-")
        os.makedirs(log_dir, exist_ok=True)
        workers = max(1, min(args.threads, len(todo)))
        print(
            "Creating %d environment(s) with %d worker(s). Logs: %s\n"
            % (len(todo), workers, log_dir)
        )
        for module, env, _yaml, _lock in todo:
            print("  queued  %s  (%s)" % (env, module))
        print()

        total = len(todo)

        def build(
            position: int, module: str, env: str, yaml: str, lock: Optional[str]
        ) -> Tuple[str, bool, str]:
            """Create one env, announcing when it starts and how long it took.

            A conda solve prints nothing to the screen on purpose (its output goes to the per-env
            log, which is only shown on failure), so this step - the slowest in the whole install,
            measured at 29 minutes for eleven environments - used to print the banner and then
            nothing at all until it was over. After five minutes of silence the honest question is
            whether it has hung, and the only way to answer it was to go and find the log.

            The lines interleave when several workers are running; each one carries its env name,
            which is what makes that readable. Flushed because stdout is a pipe under `tee`, `nohup`
            and CI, where a block-buffered progress line arrives after the thing it announces.
            """
            print("  [%d/%d] %-22s solving..." % (position, total, env), flush=True)
            started = time.time()
            result = _create_one(module, env, yaml, force, log_dir, lock)
            _m, ok, log = result
            if ok:
                print(
                    "  [%d/%d] %-22s done in %s"
                    % (position, total, env, format_duration(time.time() - started)),
                    flush=True,
                )
            else:
                print(
                    "  [%d/%d] %-22s FAILED after %s   -> %s"
                    % (position, total, env, format_duration(time.time() - started), log),
                    flush=True,
                )
            return result

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(build, i, m, e, y, lk): (m, e)
                for i, (m, e, y, lk) in enumerate(todo, start=1)
            }
            for future in concurrent.futures.as_completed(futures):
                module, env = futures[future]
                try:
                    _m, ok, _log = future.result()
                except Exception as exc:  # noqa: BLE001 - keep building the other envs
                    # The per-env line never printed in this case, so say so here or the env
                    # vanishes from the output entirely.
                    print("  FAILED  %s   -> %s" % (env, exc))
                    ok = False
                if not ok:
                    rc = 1

    if not args.no_test:
        # Only the names that resolved. Passing an unresolved one through made `doctor` invent a
        # module row for it - "nonsense  metawrap2-nonsense  MISSING  conda env not created" -
        # which reads as an environment that failed to install rather than a name that does not
        # exist, and was the confusing half of a user's bug report.
        # The *resolved* env's own base name, not what the user typed. doctor resolves a module or
        # an env name, and a 2.0.x alias is neither - passing "bin_refinement-checkm2" through made
        # it report MISSING for an environment that had just been confirmed present, and the
        # command exited 1 on a successful install.
        verifiable = []
        for name in modules:
            env = _env_for(name)
            if not env:
                continue
            verifiable.append(name if name in MODULE_ENVS else env.removeprefix("metawrap2-"))
        if verifiable:
            from .doctor import main as test_main

            print("\nVerifying the environment(s)...")
            if test_main(verifiable + ["--tools-only"]) != 0:
                rc = 1
    return rc
