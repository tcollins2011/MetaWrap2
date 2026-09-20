"""`metawrap2 install-env [module ...] [--all]` - create per-module conda envs with mamba.

MetaWrap2 manages everything through conda environments but drives them with **mamba**
(much faster solves). This installer makes sure mamba is available (bootstrapping it once
via conda if needed), creates each module's small env from ``envs/<module>.yaml``, then
verifies what it built with ``metawrap2 test``.

Environments are built several at a time (``-t``), since each one solves and downloads
independently. Each env's solver output goes to its own log file rather than the screen,
because concurrent solver output interleaved into one stream is unreadable.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import platform
import shutil
import subprocess
import tempfile
from typing import List, Optional, Tuple

from ..command import run
from ..config import MODULE_ENVS, OPTIONAL_ENVS, conda_env_exists

#: The opt-in envs keyed by the yaml/module name they are built from, so --lock can find them.
OPTIONAL_ENVS_BY_MODULE = {name[len("metawrap2-") :]: name for name in OPTIONAL_ENVS.values()}

# Env creation is dominated by solving and downloading, both of which are per-env
# independent, so building several at once is close to a linear speedup. Default is modest
# so a laptop on a slow link is not swamped; raise it with -t.
DEFAULT_PARALLEL_ENVS = 4


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


def _env_yaml(module: str) -> str:
    return os.path.join(envs_dir(), "%s.yaml" % module)


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


def create_from_lock(module: str, env: str, plat: Optional[str] = None) -> bool:
    """Create an env from its committed lockfile instead of solving the yaml."""
    path = lock_path(module, plat)
    if not os.path.isfile(path):
        return False
    manager = shutil.which("conda") or shutil.which("mamba") or shutil.which("micromamba")
    if not manager:
        return False
    return run([manager, "create", "-y", "-n", env, "--file", path], tool="conda", check=False) == 0


def _create_one(
    module: str, env: str, yaml: str, force: bool, log_dir: str
) -> Tuple[str, bool, str]:
    """Create one env, capturing its output to a log file. Returns (module, ok, log path).

    Output goes to a file rather than the screen because several of these run at once and
    interleaved solver output is unreadable. The log path is reported on failure.
    """
    log = os.path.join(log_dir, "install-env-%s.log" % module)
    cmd = ["mamba", "env", "create", "-y", "-f", yaml]
    if force:
        cmd.append("--force")
    with open(log, "w") as fh:
        fh.write("$ %s\n\n" % " ".join(cmd))
        fh.flush()
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return module, rc == 0, log


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(
        prog="metawrap2 install-env",
        description="Create the per-module conda environments, several at a time.",
    )
    ap.add_argument("modules", nargs="*", help="modules whose envs to create")
    ap.add_argument("--all", action="store_true", help="create envs for every module")
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

    modules = sorted(MODULE_ENVS) if args.all else args.modules
    if not modules:
        ap.error("specify one or more modules, or --all")

    if args.lock:
        # Locking reads existing environments; it does not need mamba bootstrapped.
        rc = 0
        for module in modules:
            env = MODULE_ENVS.get(module) or OPTIONAL_ENVS_BY_MODULE.get(module, module)
            if not conda_env_exists(env):
                print("%s does not exist yet; create it before locking it." % env)
                rc = 1
                continue
            path = write_lockfile(module, env)
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
    for module in modules:
        yaml = _env_yaml(module)
        if not os.path.isfile(yaml):
            print("No env file for module '%s' (%s); skipping." % (module, yaml))
            rc = 1
            continue
        env = MODULE_ENVS[module]
        if conda_env_exists(env) and not args.force:
            print("%s already exists; skipping (use --force to recreate)." % env)
            continue
        todo.append((module, env, yaml))

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
        for module, env, _yaml in todo:
            print("  queued  %s  (%s)" % (env, module))
        print()

        done = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_create_one, m, e, y, args.force, log_dir): (m, e) for m, e, y in todo
            }
            for future in concurrent.futures.as_completed(futures):
                module, env = futures[future]
                done += 1
                try:
                    _m, ok, log = future.result()
                except Exception as exc:  # noqa: BLE001 - keep building the other envs
                    ok, log = False, "%s" % exc
                if ok:
                    print("  [%d/%d] OK      %s" % (done, len(todo), env))
                else:
                    print("  [%d/%d] FAILED  %s   -> %s" % (done, len(todo), env, log))
                    rc = 1

    if not args.no_test:
        from .doctor import main as test_main

        print("\nVerifying the environment(s)...")
        if test_main(list(modules) + ["--tools-only"]) != 0:
            rc = 1
    return rc
