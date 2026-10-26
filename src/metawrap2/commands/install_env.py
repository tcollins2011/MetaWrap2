"""`metawrap2 install-env [module ...] [--all]` - create per-module conda envs with mamba.

MetaWrap2 manages everything through conda environments but drives them with **mamba**
(much faster solves). This installer makes sure mamba is available (bootstrapping it once
via conda if needed), creates each module's small env from ``envs/<module>.yaml``, then
verifies what it built with ``metawrap2 test``.

Environments are built several at a time (``-t``), since each one solves and downloads
independently. Each env's solver output goes to its own log file, and a line is printed when each
env starts and finishes, with its duration.
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

# Solving and downloading are independent per env, so parallel builds scale nearly linearly. The
# default is kept low for slow connections; raise it with -t.
DEFAULT_PARALLEL_ENVS = 4


def format_duration(seconds: float) -> str:
    """Format a duration as ``45s``, ``3m12s`` or ``1h04m``, for progress lines."""
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

    ``$METAWRAP2_ENVS_DIR`` takes precedence, for user-edited specs. Otherwise walk up from this
    file, which finds ``envs/`` in a source checkout, an editable install, or the installed
    package.
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
    # Fallback: the old repo-relative path, so an error names a plausible location.
    return os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "envs"))


def _env_yaml(env: str) -> str:
    """Locate the yaml that builds conda environment *env*.

    Specs are named after the environment: env ``metawrap2-<basename>`` is built from
    ``<basename>.yaml``. This lets modules share an env; ``coverage`` uses ``metawrap2-bins``
    and has no ``coverage.yaml``.
    """
    base = env.removeprefix("metawrap2-")
    return os.path.join(envs_dir(), "%s.yaml" % base)


def _env_for(name: str) -> Optional[str]:
    """Alias for config.resolve_env_name, used by tests and callers here.

    Delegates instead of reimplementing; separate copies here and in `check` and `doctor` had
    diverged.
    """
    return resolve_env_name(name)


def _resolve_or_explain(name: str) -> Optional[str]:
    """Like _env_for, but prints suggestions when the name is unknown."""
    env = _env_for(name)
    if env:
        if name in LEGACY_ENV_ALIASES:
            print(
                "'%s' was an environment name in MetaWrap2 2.0.x. Environments are now grouped by "
                "workflow; it is part of '%s', which will be installed." % (name, env)
            )
        return env
    for line in explain_unknown_name(name):
        print(line)
    return None


def ensure_mamba() -> bool:
    """Make sure `mamba` is on PATH; bootstrap it into the base env via conda if not.

    Returns True if mamba is available afterwards. This is the only place conda itself is used.
    """
    if shutil.which("mamba"):
        return True
    if not shutil.which("conda"):
        print(
            "Neither mamba nor conda was found on PATH. Install Miniforge "
            "(https://github.com/conda-forge/miniforge), which includes mamba, then re-run."
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

    Uses ``conda list --explicit --md5``, which lists the URL and checksum of every installed
    package. ``conda create --file`` recreates the env from it without conda-lock, and changes
    are readable in a diff.
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
        "# Exact builds of every package in the %s environment, for reproducing a run.\n"
        "# Generated with:  metawrap2 install-env %s --lock\n"
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

    Output goes to a file because several builds run at once. The log path is reported on failure.
    """
    log = os.path.join(log_dir, "install-env-%s.log" % module)
    # A rebuild removes the env first: mamba 2.x's `env create` rejects --force ("The following
    # argument was not expected: --force").
    if lock:
        # An explicit lockfile names every package build, so no solve is needed.
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
        help="create every environment, including the opt-in ones (gtdbtk, eggnog, rgi, tracs, "
        "the benchmark's amber). Several of these are large",
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
        help="how many environments to build at once (default %d). Each "
        "env solves and downloads independently. Use 1 for serial "
        "output that is easier to read." % DEFAULT_PARALLEL_ENVS,
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
        help="create the environment(s) from their committed lockfiles instead "
        "of solving the yaml, which reproduces the exact builds",
    )
    ap.add_argument(
        "--log-dir",
        default=None,
        help="where to write per-env build logs (default: a temp directory)",
    )
    args = ap.parse_args(argv)
    # `--force` is also a global flag, which the CLI strips before this parser sees it, so check
    # both.
    force = args.force or runner.force

    # --all iterates env names, not modules, because several modules share an env. It includes
    # the separate envs, both those a default run uses (checkm2, learned) and those a flag selects.
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
            # Named after the env, not the name typed, so `--lock binning` and `--lock bins` write
            # the same file that --from-lock reads.
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
    # Envs already handled in this run, whether queued or skipped as existing, so a shared env is
    # reported once.
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
            # Modules sharing an env; building it twice at once would write one prefix twice.
            print("%s is already covered by this run; skipping (%s shares it)." % (env, module))
            continue
        seen_envs.add(env)
        if conda_env_exists(env) and not force:
            # Print the prefix: `conda_env_exists` searches every conda installation on the
            # machine, so the existing env may belong to a different conda.
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
                    "No %s lockfile for %s; solving %s instead, so the builds may differ."
                    % (lock_platform(), env, os.path.basename(yaml))
                )
        todo.append((module, env, yaml, lock))

    if todo and runner.dry_run:
        # The global --dry-run prints the commands without running them.
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
            """Create one env, printing when it starts and how long it took.

            Solver output goes to the per-env log, and this step can be long (29 minutes for eleven
            environments), so start and finish lines show progress. Each line names its env, since
            workers interleave. Output is flushed because stdout is often a pipe (`tee`, `nohup`,
            CI).
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
                    # build() did not print a result line in this case, so print one here.
                    print("  FAILED  %s   -> %s" % (env, exc))
                    ok = False
                if not ok:
                    rc = 1

    if not args.no_test:
        # Verify only names that resolved; otherwise `doctor` reports an unknown name as a
        # missing env. Module names are passed as typed; anything else (such as a 2.0.x alias
        # like "bin_refinement-checkm2", which doctor does not know) is passed as the resolved
        # env's base name.
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
