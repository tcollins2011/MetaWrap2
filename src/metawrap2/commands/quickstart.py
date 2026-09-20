"""`metawrap2 quickstart` - one command from a fresh checkout to a pipeline that runs.

Getting started used to mean four commands and a reading of the docs before anything ran: create
the environments, download the databases, point the config at them, check the result. Each step is
simple, and none of them is what the user came to do. This runs them in the right order.

    metawrap2 quickstart -t 8                     # environments, then verify
    metawrap2 quickstart --binning-only           # just what genome recovery needs
    metawrap2 quickstart --databases small        # also fetch the capped test databases
    metawrap2 quickstart --databases full         # also fetch the production databases

It is a *composition*, not a reimplementation: it calls the same ``install-env``, ``install-db``
and ``doctor`` code paths you would have run by hand, and stops at the first failing step with the
specific command to resume from. Nothing it does is unavailable separately.

**It never downloads a database you did not ask for.** The databases run from 550 MB to ~300 GB
and take hours to days, and where they go is a decision about the machine, not something to
default. So ``--databases`` is required to fetch anything. With no answer given, quickstart
installs the environments, then prints what each module needs, how large it is, and the exact
command - and, if it is running interactively, offers to do it. On a non-interactive run (CI, a
batch script) it never asks and never downloads.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from typing import List, Optional, Sequence, Tuple

from ..config import MODULE_ENVS
from ..logging import announcement, comm, warning

#: The modules genome recovery actually needs. `--binning-only` installs these, which is a much
#: smaller commitment than all thirteen environments and covers the path most users came for.
BINNING_PATH: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "quant_bins",
)

#: Databases the genome-recovery path needs. The Kraken2 and BLAST databases are the enormous
#: ones and are only needed by kraken2/blobology/classify_bins, so they are left out here.
BINNING_DATABASES: Tuple[str, ...] = ("checkm", "host")

#: How many database downloads to run at once. Network-bound, so a few in parallel is close to
#: linear, but more than this mostly just saturates a link.
DEFAULT_DB_JOBS = 3


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


def _step(number: int, title: str) -> None:
    announcement("STEP %d/3: %s" % (number, title))


def _explain_failure(step: str, command: str) -> None:
    warning(
        "%s failed.\n"
        "Nothing quickstart does is special - it runs the ordinary commands in order - so you\n"
        "can pick up exactly where it stopped:\n"
        "    %s\n"
        "The step's own output above says what went wrong. `metawrap2 doctor` reports what is\n"
        "already working, so a retry only redoes what is missing." % (step, command)
    )


def _database_command(names: Sequence[str], jobs: int, small: bool) -> str:
    return "metawrap2 install-db %s%s -t %d" % (
        " ".join(names),
        " --small" if small else "",
        jobs,
    )


def _print_database_guidance(names: Sequence[str], jobs: int) -> None:
    """Say what is still needed, how big it is, and exactly what to type."""
    from .install_db import DATABASES

    print("\n  Databases are installed separately, because they are large and where they go is")
    print("  your decision. What each one is for, and how big:\n")
    print("    %-10s %-12s %-10s %s" % ("name", "full size", "--small", "needed by"))
    print("    " + "-" * 76)
    listed = list(DATABASES) if names == ["--all"] else [n for n in names if n in DATABASES]
    for name in listed:
        db = DATABASES[name]
        print("    %-10s %-12s %-10s %s" % (name, db.size, db.small_size or "-", db.used_by))
    print(
        "\n  Capped-but-real variants, enough to see the pipeline work end to end:\n"
        "      %s\n"
        "  The production databases:\n"
        "      %s\n"
        "  Either way the paths are written into your config for you. `metawrap2 install-db\n"
        "  --list` shows the full set, and --dir puts them somewhere other than ~/metawrap2_dbs.\n"
        % (
            _database_command(listed or names, jobs, True),
            _database_command(listed or names, jobs, False),
        )
    )


def _ask_about_databases(names: Sequence[str], jobs: int) -> str:
    """Offer to install databases, but only when a person is there to answer.

    Returns "none", "small" or "full". A non-interactive run always gets "none": a batch script
    or a CI job must never be able to start a multi-hundred-gigabyte download because nobody was
    watching to say no.
    """
    if not (sys.stdin.isatty() and sys.stdout.isatty()):
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
        "\n        --metabat2 --maxbin2 --concoct sampleA_1.fastq sampleA_2.fastq\n"
        % (threads, threads)
    )
    if small:
        comm(
            "Those are the capped databases - right for checking the pipeline runs, not for "
            "publishing. Use `metawrap2 install-db --all` when you need the real ones."
        )


def _run_install_env(modules: Sequence[str], threads: int, from_lock: bool) -> int:
    from .install_env import main as install_env_main

    argv: List[str] = list(modules) + ["-t", str(threads)]
    if from_lock:
        argv.append("--from-lock")
    # --no-test because the final doctor step verifies everything anyway; testing each env as it
    # is built and then again at the end doubles the slowest part for no extra information.
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
    ap = argparse.ArgumentParser(
        prog="metawrap2 quickstart",
        description="Install everything MetaWrap2 needs and verify it, in one command.",
        epilog="Every step is an ordinary metawrap2 command; run them separately if you prefer.",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="how many conda environments to build at once (default 4). More is faster if "
        "your machine and network can take it",
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
        help="whether to download databases too. 'small' fetches capped-but-real variants "
        "(good for a first run), 'full' the production ones (hundreds of GB, many hours), "
        "'none' skips them. Omit it and nothing is downloaded without asking you first",
    )
    ap.add_argument(
        "--skip-envs", action="store_true", help="skip environment creation (databases only)"
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
        modules = sorted(MODULE_ENVS)

    # Only the genome-recovery path has a short database list; a full install needs all of them.
    db_names = list(BINNING_DATABASES) if (args.binning_only or args.modules) else ["--all"]

    if not _package_manager():
        warning(
            "No conda, mamba or micromamba on PATH.\n"
            "MetaWrap2 runs each module's tools inside a conda environment, so one of them has\n"
            "to exist first. Miniforge is the usual choice:\n"
            "    https://github.com/conda-forge/miniforge#install\n"
            "(If you would rather manage the tools yourself, set use_conda_envs = false in\n"
            "metawrap2.toml and MetaWrap2 will call them straight off your PATH.)"
        )
        return 1

    started = time.time()
    announcement("METAWRAP2 QUICKSTART")

    if not args.skip_envs:
        _step(1, "creating %d conda environment(s)" % len(modules))
        comm("This is the slow part. Existing environments are skipped, so a retry is cheap.")
        if _run_install_env(modules, args.threads, args.from_lock) != 0:
            _explain_failure(
                "Environment creation",
                "metawrap2 install-env %s -t %d" % (" ".join(modules), args.threads),
            )
            return 1

    choice = args.databases
    if choice is None:
        choice = _ask_about_databases(db_names, args.db_jobs)

    if choice != "none":
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

    _step(3, "checking what works")
    # Databases nobody asked for would be reported as missing and make the check look failed, so
    # they are only part of the verdict once they have actually been installed.
    healthy = _run_doctor(modules, args.config, databases=choice != "none") == 0

    announcement("QUICKSTART FINISHED IN %s" % _elapsed(started))
    if not healthy:
        warning(
            "Some of it is not ready - the table above says which.\n"
            "`metawrap2 doctor --fix` recreates exactly the environments that are broken and\n"
            "leaves the healthy ones alone. Everything that did install is still installed."
        )
        return 1

    if choice == "none":
        comm("Environments are ready. No databases were downloaded.")
        _print_database_guidance(db_names, args.db_jobs)
    else:
        comm("Everything checked out.")
    _print_next_steps(args.threads, small=choice == "small")
    return 0
