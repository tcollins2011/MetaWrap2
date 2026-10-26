"""`metawrap2 quickstart`: set up a fresh checkout in one command.

Runs the setup steps in order: create the environments, optionally download the databases and
write their paths to the config, then check the result.

    metawrap2 quickstart -t 8                     # environments, then verify
    metawrap2 quickstart --binning-only           # just what genome recovery needs
    metawrap2 quickstart --databases small        # also fetch the capped test databases
    metawrap2 quickstart --databases full         # also fetch the production databases

It calls the same ``install-env``, ``install-db`` and ``doctor`` code as the separate commands,
and stops at the first failing step with the command to resume from.

Databases are only downloaded with ``--databases`` (they range from 550 MB to ~300 GB and take
hours to days). Without it, quickstart installs the environments, then prints what each module
needs, its size and the command to fetch it, and offers to run it when interactive. A
non-interactive run (CI, a batch script) never prompts and never downloads.
"""

from __future__ import annotations

import shutil
import sys
import time
from typing import List, Optional, Sequence, Tuple

from ..cliformat import make_parser
from ..config import MODULE_ENVS
from ..logging import announcement, comm, warning

#: The modules genome recovery needs; `--binning-only` installs only these environments.
BINNING_PATH: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "quant_bins",
)

#: Databases the genome-recovery path needs: CheckM2's and GUNC's for bin_refinement's default
#: scoring and chimerism check, and the host genome for read_qc. CheckM1's data is only needed for
#: --checkm1 (`install-db checkm`). The large Kraken2 and BLAST databases are only used by
#: kraken2/blobology/classify_bins and are left out.
BINNING_DATABASES: Tuple[str, ...] = ("checkm2", "gunc", "host")

#: Parallel database downloads. They are network-bound; more than a few saturates the link.
DEFAULT_DB_JOBS = 3


def _default_tier_envs() -> List[str]:
    """Opt-in env names the default binning tier needs, so an install covers a plain run.

    The default tier is metaBAT2 + SemiBin2; SemiBin2 is in metawrap2-learned with the other
    pytorch binners. Derived from binning.DEFAULT_PRESET, so a tier of classical binners returns
    an empty list.
    """
    from ..config import OPTIONAL_ENVS
    from ..modules.binning import DEFAULT_PRESET, PRESETS

    names = []
    for binner in PRESETS[DEFAULT_PRESET]:
        env = OPTIONAL_ENVS.get(binner)
        if env:
            name = env.removeprefix("metawrap2-")
            if name not in names:
                names.append(name)
    return names


def _default_run_envs(modules: Sequence[str]) -> List[str]:
    """The envs outside MODULE_ENVS that a default run of *modules* uses, as install-env names.

    *modules* empty means a full install. bin_refinement scores with CheckM2, which has its own
    env (config.DEFAULT_OPTIONAL_ENVS); the default binning tier needs metawrap2-learned.
    """
    from ..config import DEFAULT_OPTIONAL_ENVS

    names: List[str] = []
    for key, users in DEFAULT_OPTIONAL_ENVS.items():
        if not modules or any(module in modules for module in users):
            names.append(key)
    if not modules or "binning" in modules:
        names += [n for n in _default_tier_envs() if n not in names]
    return names


def _package_manager() -> Optional[str]:
    for program in ("mamba", "micromamba", "conda"):
        if shutil.which(program):
            return program
    return None


def _elapsed(started: float) -> str:
    seconds = time.time() - started
    if seconds < 90:
        return "%.0fs" % seconds
    if seconds < 5400:
        return "%.0fm" % (seconds / 60)
    return "%.1fh" % (seconds / 3600)


#: quickstart always numbers three steps, because whether step 2 runs may only be decided by a
#: prompt after step 1 finishes.
TOTAL_STEPS = 3


def _step(number: int, title: str) -> None:
    announcement("STEP %d/%d: %s" % (number, TOTAL_STEPS, title))


def _skipped_step(number: int, title: str, why: str) -> None:
    """Announce a step that is not going to run, and why.

    Printed so the step numbering has no gaps (e.g. "STEP 2/3 ... (skipped: --databases none)").
    """
    announcement("STEP %d/%d: %s (skipped: %s)" % (number, TOTAL_STEPS, title, why))


def _unique_envs(names: Sequence[str]) -> List[str]:
    """The distinct conda environments *names* resolve to, in first-seen order.

    Several modules share one environment (metawrap2-reads serves five), so the count of names is
    not the count of installs; a full install names 23 and builds 11.
    """
    from ..config import resolve_env_name

    envs: List[str] = []
    for name in names:
        env = resolve_env_name(name) or name
        if env not in envs:
            envs.append(env)
    return envs


def _print_modules_without_envs() -> None:
    """No-op: `doctor` now reports every module, including those with no conda env.

    Kept as the place to report a future dependency kind (a system library, say).
    """
    return


def _explain_failure(step: str, command: str) -> None:
    warning(
        "%s failed.\n"
        "quickstart runs the ordinary metawrap2 commands in order, so you can continue from\n"
        "this step with:\n"
        "    %s\n"
        "The step's output above shows the error. `metawrap2 doctor` reports what already\n"
        "works; a retry skips what is installed." % (step, command)
    )


def _database_command(names: Sequence[str], jobs: int, small: bool) -> str:
    return "metawrap2 install-db %s%s -t %d" % (
        " ".join(names),
        " --small" if small else "",
        jobs,
    )


def _print_database_guidance(names: Sequence[str], jobs: int) -> None:
    """Print the databases still needed, their sizes and the commands to install them."""
    from .install_db import DATABASES

    print("\n  Databases are installed separately because they are large and you choose where")
    print("  they go. What each one is used for, and its size:\n")
    print("    %-10s %-12s %-10s %s" % ("name", "full size", "--small", "needed by"))
    print("    " + "-" * 76)
    listed = list(DATABASES) if names == ["--all"] else [n for n in names if n in DATABASES]
    for name in listed:
        db = DATABASES[name]
        print("    %-10s %-12s %-10s %s" % (name, db.size, db.small_size or "-", db.used_by))
    print(
        "\n  Small capped variants, enough to test the pipeline end to end:\n"
        "      %s\n"
        "  The production databases:\n"
        "      %s\n"
        "  Both commands write the paths into your config. `metawrap2 install-db --list`\n"
        "  shows the full set; --dir installs them somewhere other than ~/metawrap2_dbs.\n"
        % (
            _database_command(listed or names, jobs, True),
            _database_command(listed or names, jobs, False),
        )
    )


def _print_database_status(config: Optional[str]) -> bool:
    """Print each database key's value and whether it exists. True if none is broken.

    Used by `--skip-envs` ("databases only") in place of the environment table, which would show
    every module MISSING. Keys that are configured but missing and keys that are not set are
    listed separately.
    """
    from ..config import check_databases, load_settings

    rows = check_databases(load_settings(config))
    configured = [(key, value, present) for key, value, present, _used in rows if value]
    unset = [key for key, value, _present, _used in rows if not value]

    print("\nMetaWrap2 database status")
    print("  " + "-" * 70)
    for key, value, present in configured:
        print("  %-22s %-38s [%s]" % (key, value, "ok" if present else "MISSING"))
    if not configured:
        print("  nothing configured yet.")
    if unset:
        print("  not set (%d): %s" % (len(unset), ", ".join(unset)))
    print("  " + "-" * 70)

    broken = [key for key, _value, present in configured if not present]
    if broken:
        warning(
            "These database paths are in your config but are not on disk: %s\n"
            "The path is wrong or the download did not finish. `metawrap2 install-db "
            "<name>` downloads one and updates its path in the config." % ", ".join(broken)
        )
    return not broken


def _announce_non_interactive_default(names: Sequence[str], jobs: int) -> None:
    """Report that an unattended run downloads no databases, and what to pass instead.

    Under `sbatch` or in CI no one can answer a prompt, so nothing is downloaded.
    """
    # The command is on its own line so it does not wrap mid-command.
    comm(
        "No --databases given and stdin is not a terminal; downloading nothing.\n"
        "Pass --databases small or --databases full to fetch them here, or run this later:\n"
        "%s" % _database_command(names, jobs, True)
    )


def _ask_about_databases(names: Sequence[str], jobs: int) -> str:
    """Offer to install databases, but only when a person is there to answer.

    Returns "none", "small" or "full". A non-interactive run (batch script, CI) always gets "none".
    """
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
        _announce_non_interactive_default(names, jobs)
        return "none"
    _print_database_guidance(names, jobs)
    print("  Download databases now?")
    print("    [1] no  - I will do it later (default)")
    print("    [2] the capped --small variants")
    print("    [3] the full production databases - hundreds of GB, hours to days")
    try:
        answer = input("  Choose 1/2/3 [1]: ").strip()
    except (EOFError, KeyboardInterrupt):
        print()
        return "none"
    return {"2": "small", "3": "full"}.get(answer, "none")


def _print_next_steps(threads: int, small: bool) -> None:
    print(
        "\n  A first run, on your own reads:\n"
        "\n    metawrap2 run --example > samples.toml     # then edit it to point at your reads"
        "\n    metawrap2 run samples.toml -o study/ -t %d"
        "\n\n  or one module at a time, e.g.:\n"
        "\n    metawrap2 binning -a assembly.fasta -o binning_out -t %d \\"
        "\n        --metabat2 --maxbin2 sampleA_1.fastq sampleA_2.fastq\n" % (threads, threads)
    )
    if small:
        comm(
            "These are the capped databases, for testing the pipeline and not for published "
            "results. Use `metawrap2 install-db --all` to install the full ones."
        )


def _run_install_env(
    modules: Sequence[str], threads: int, from_lock: bool, log_dir: Optional[str] = None
) -> int:
    from .install_env import main as install_env_main

    argv: List[str] = list(modules) + ["-t", str(threads)]
    if from_lock:
        argv.append("--from-lock")
    if log_dir:
        # install-env only prints its log paths on failure; --log-dir keeps the logs of a
        # successful install too.
        argv += ["--log-dir", log_dir]
    # --no-test: the final doctor step verifies every env, so testing each one here is redundant.
    argv.append("--no-test")
    return install_env_main(argv)


def _run_install_db(
    names: Sequence[str], directory: Optional[str], small: bool, jobs: int, config: Optional[str]
) -> int:
    from .install_db import main as install_db_main

    argv: List[str] = list(names) + ["-t", str(jobs)]
    if small:
        argv.append("--small")
    if directory:
        argv += ["--dir", directory]
    if config:
        argv += ["--config", config]
    return install_db_main(argv)


def _run_doctor(modules: Sequence[str], config: Optional[str], databases: bool) -> int:
    from .doctor import main as doctor_main

    argv: List[str] = list(modules) + ["--no-unit-tests"]
    if not databases:
        argv.append("--no-databases")
    if config:
        argv += ["--config", config]
    return doctor_main(argv)


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 quickstart",
        description="Build MetaWrap2's conda environments, optionally download databases, and "
        "check both. It does not install MetaWrap2's Python extras: `diff_abundance` "
        "needs `pip install 'metawrap2[stats]'`, and `doctor` reports if it is missing.",
        epilog="Every step is an ordinary metawrap2 command; run them separately if you prefer.",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="how many conda environments to build at once (default 4). More is faster if "
        "the machine and network allow it",
    )
    ap.add_argument(
        "--binning-only",
        action="store_true",
        help="install only what genome recovery needs (%s) instead of every module"
        % ", ".join(BINNING_PATH),
    )
    ap.add_argument(
        "--modules", nargs="+", metavar="MODULE", help="install exactly these modules' envs"
    )
    ap.add_argument(
        "--databases",
        choices=("none", "small", "full"),
        default=None,
        help="whether to download databases. 'small' downloads capped variants "
        "(for a first run), 'full' the production ones (hundreds of GB, many hours), "
        "'none' skips them. Without this option, quickstart asks before downloading",
    )
    ap.add_argument(
        "--skip-envs",
        action="store_true",
        help="databases only: skip environment creation, and report database status at the end "
        "instead of module status",
    )
    ap.add_argument("--db-dir", help="where to install databases (default ~/metawrap2_dbs)")
    ap.add_argument(
        "--db-jobs",
        type=int,
        default=DEFAULT_DB_JOBS,
        help="how many databases to download at once (default %d)" % DEFAULT_DB_JOBS,
    )
    ap.add_argument(
        "--from-lock",
        action="store_true",
        help="build environments from the committed lockfiles rather than by solving the yamls",
    )
    ap.add_argument("--config", help="config file to write database paths into")
    ap.add_argument(
        "--log-dir",
        help="keep the per-environment build logs here (default: a temp directory, so the "
        "warnings of a successful install are not kept)",
    )
    args = ap.parse_args(argv)

    if args.modules:
        modules = list(args.modules)
        unknown = [m for m in modules if m not in MODULE_ENVS]
        if unknown:
            print(
                "Unknown module(s): %s\nKnown: %s"
                % (", ".join(unknown), ", ".join(sorted(MODULE_ENVS)))
            )
            return 1
    elif args.binning_only:
        modules = list(BINNING_PATH)
    else:
        # An empty list makes `doctor` report every CLI module, including the two without a
        # conda env.
        modules = []

    # Only the genome-recovery path has a short database list; a full install needs all of them.
    db_names = list(BINNING_DATABASES) if (args.binning_only or args.modules) else ["--all"]

    if not _package_manager():
        warning(
            "No conda, mamba or micromamba on PATH.\n"
            "MetaWrap2 runs each module's tools inside a conda environment, so one of these\n"
            "is required. To install Miniforge:\n"
            "    https://github.com/conda-forge/miniforge#install\n"
            "(To manage the tools yourself, set use_conda_envs = false in metawrap2.toml and\n"
            "MetaWrap2 will run them from your PATH.)"
        )
        return 1

    started = time.time()
    announcement("METAWRAP2 QUICKSTART")

    # The modules' envs (every module env for a full install) plus the separate envs a default
    # run uses: CheckM2 and the learned binners. `doctor` is only asked about the modules, since
    # "checkm2" and "learned" are not modules.
    from ..config import ALL_MODULE_ENVS

    envs_to_build = list(modules) or [env.removeprefix("metawrap2-") for env in ALL_MODULE_ENVS]
    envs_to_build += _default_run_envs(modules)

    if args.skip_envs:
        _skipped_step(1, "conda environments", "--skip-envs")
    else:
        unique = _unique_envs(envs_to_build)
        _step(1, "creating %d conda environment(s)" % len(unique))
        comm(
            "This step is slow. Each environment's solver output goes to its own log, so the\n"
            "lines below are the only progress shown. Eleven environments took 29 minutes at\n"
            "-t 4. Existing environments are skipped on a retry."
        )
        # Without --log-dir, install-env uses its own temp directory.
        if args.log_dir:
            rc_env = _run_install_env(envs_to_build, args.threads, args.from_lock, args.log_dir)
        else:
            rc_env = _run_install_env(envs_to_build, args.threads, args.from_lock)
        if rc_env != 0:
            _explain_failure(
                "Environment creation",
                "metawrap2 install-env %s -t %d" % (" ".join(envs_to_build), args.threads),
            )
            return 1

    choice = args.databases
    if choice is None:
        choice = _ask_about_databases(db_names, args.db_jobs)

    if choice == "none":
        _skipped_step(
            2,
            "databases",
            "--databases none" if args.databases == "none" else "no --databases given",
        )
    else:
        _step(2, "installing the %s databases" % choice)
        if (
            _run_install_db(db_names, args.db_dir, choice == "small", args.db_jobs, args.config)
            != 0
        ):
            _explain_failure(
                "Database installation",
                "metawrap2 install-db %s%s -t %d"
                % (" ".join(db_names), " --small" if choice == "small" else "", args.db_jobs),
            )
            return 1

    if args.skip_envs:
        # `--skip-envs` means "databases only", so check databases, not modules (which would all
        # be missing and give exit 1).
        _step(3, "checking the databases")
        healthy = _print_database_status(args.config)
        announcement("QUICKSTART FINISHED IN %s" % _elapsed(started))
        if not healthy:
            return 1
        if choice == "none":
            _print_database_guidance(db_names, args.db_jobs)
        return 0

    _step(3, "checking what works")
    # Databases are only checked if they were installed in this run; otherwise they would all be
    # reported missing.
    healthy = _run_doctor(modules, args.config, databases=choice != "none") == 0
    _print_modules_without_envs()

    announcement("QUICKSTART FINISHED IN %s" % _elapsed(started))
    if not healthy:
        warning(
            "Some modules are not ready; the table above says which.\n"
            "`metawrap2 doctor --fix` recreates only the broken environments. The environments\n"
            "that did build are still installed."
        )
        return 1

    if choice == "none":
        comm("Environments are ready. No databases were downloaded.")
        _print_database_guidance(db_names, args.db_jobs)
    else:
        comm("Everything checked out.")
    _print_next_steps(args.threads, small=choice == "small")
    return 0
