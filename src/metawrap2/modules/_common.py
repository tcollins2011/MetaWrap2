"""Shared helpers for MetaWrap2 modules.

Modules import the banner helpers, the conda-env lookup and some filesystem, read and run
bookkeeping utilities from here. Everything else stays inline in each module.
"""

from __future__ import annotations

import argparse
import fcntl
import os
import re
import resource
import shutil
import sys
from typing import Any, Dict, List, Optional, Sequence

from .. import __version__
from .. import aligners as _aligners
from .. import command as _cmd
from .. import history as _history
from .. import interrupt as _interrupt
from .. import limits as _limits
from .. import logging as _log
from .. import runlock as _runlock
from .. import scratch as _scratch
from .. import usage as _usage
from .. import validate as _validate
from ..checkpoint import make_checkpoint
from ..command import ToolError, run, runner  # re-exported for modules
from ..config import MODULE_ENVS, OPTIONAL_ENVS, Settings, conda_env_exists
from ..constants import FASTA_EXTENSIONS
from ..io import reads as _reads
from ..io.seqio import fasta_files as fasta_files_in
from ..provenance import RunRecorder

#: The lock held on the current run's output directory, if any. Not a field on RunRecorder,
#: because the recorder is serialised into run_config.json and a RunLock is not JSON-serialisable.
_active_lock: Optional[_runlock.RunLock] = None

#: The disk high-water monitor for the current run, held here for the same reason as the lock.
_active_disk_monitor: Optional[_limits.DiskMonitor] = None

#: The GPU memory sampler for the current run, held here for the same reason.
#:
#: Needed as well as ``usage.run_measured`` because a standalone module measures itself from the
#: inside and ``rusage`` has no GPU figure; without it, the three GPU-capable binners reported no
#: peak VRAM when invoked directly.
_active_gpu_watcher: Optional[_usage.GpuWatcher] = None

#: Environment variable the ``run`` driver sets for every step it starts, so a module can tell a
#: step of a study from a run of its own.
#:
#: It decides who prints the cost line. The driver measures each step from outside with
#: ``os.wait4`` and prints ``finished <step> (Ns wall, ...)``, so a module under ``run`` does not
#: print its own (different) figure. Unset means the module runs standalone.
DRIVEN_BY_RUN_ENV = "METAWRAP2_IN_STUDY"

#: The run whose cost has already been reported, so the cost line is printed once.
#: ``finish_run`` can be called twice per run: pathways calls it before an early ``return 0``
#: inside the ``try`` whose ``finally`` calls it again.
_cost_reported_for: Optional[str] = None

# Banner logging (same look as the original comm/announcement/warning/error).
comm = _log.comm
announcement = _log.announcement
warning = _log.warning
error = _log.error

__all__ = [
    "DRIVEN_BY_RUN_ENV",
    "ToolError",
    "absolutize_paths",
    "add_checkm2_only_args",
    "add_scorer_args",
    "announcement",
    "check_fasta",
    "check_fastq",
    "collect_read_pairs",
    "collect_reads",
    "comm",
    "dry_run",
    "ensure_dir",
    "env_for",
    "error",
    "expect_dir",
    "expect_file",
    "expect_produced",
    "file_size",
    "finish_run",
    "force",
    "fs_copy",
    "fs_move",
    "fs_remove",
    "fs_replace_dir",
    "listdir",
    "make_checkpoint",
    "memory_arg",
    "optional_env_for",
    "read_layout",
    "require_databases",
    "require_file",
    "require_nonempty_dir",
    "resolve_memory_gb",
    "resolve_threads",
    "run",
    "scorer_notes",
    "scorer_preflight",
    "self_usage",
    "start_run",
    "threads_arg",
    "validate_inputs",
    "warning",
]


def _require_env_manager(env: str) -> None:
    """Fail now if the env manager that runs tools inside envs is not on PATH.

    Checked at the start of a run, not when the first env command is built, because a module can
    do real work before then (bin_refinement builds every bin-set combination in process before
    CheckM runs).
    """
    if runner.dry_run:
        return
    if shutil.which(runner.env_manager) is None:
        error(str(_cmd.EnvManagerMissing(runner.env_manager, env)))


def env_for(module: str, settings: Settings) -> Optional[str]:
    """Conda env name for *module*, or None if the user disabled conda envs.

    Fails at the start with the command to create the env if it is missing, instead of letting
    the first tool call fail inside ``mamba run`` with a missing-prefix message.
    """
    if not settings.use_conda_envs:
        return None
    env = MODULE_ENVS.get(module)
    if env:
        _require_env_manager(env)
    if env and not conda_env_exists(env):
        error(
            "The conda environment '%s' for the %s module does not exist.\n"
            "Create it with:\n"
            "    metawrap2 install-env %s\n"
            "or run the module's tools from your PATH by setting "
            "use_conda_envs = false under [settings] in metawrap2.toml." % (env, module, module)
        )
    return env


def optional_env_for(key: str, settings: Settings, needed_by: str) -> Optional[str]:
    """Conda env name for the opt-in tool *key* names, or None if conda envs are off.

    The same check as :func:`env_for`, for opt-in envs. Modules should call it up front: an
    opt-in tool can run late in a pipeline (a binner runs after all reads are aligned), so a
    missing env would otherwise surface after hours of work.

    *needed_by* names what asked for the env in the message. It is passed by the caller because
    ``"--" + key`` is not always a real flag: GUNC runs by default (off with ``--skip-gunc``) and
    CoverM is ``--quantifier coverm``.
    """
    if not settings.use_conda_envs:
        return None
    env = OPTIONAL_ENVS[key]
    _require_env_manager(env)
    if not conda_env_exists(env):
        error(
            "The conda environment '%s' does not exist, and %s needs it.\n"
            "Create it with:\n"
            "    metawrap2 install-env %s\n"
            "or run the tool from your PATH by setting use_conda_envs = false under "
            "[settings] in metawrap2.toml." % (env, needed_by, env[len("metawrap2-") :])
        )
    return env


def require_databases(module: str, settings: Settings, keys: List[str]) -> None:
    """Fail fast if a database *module* cannot run without is unset or absent."""
    problems = []
    for key in keys:
        path = settings.db(key)
        if not path:
            problems.append("%s is not set" % key)
        elif not os.path.exists(path):
            problems.append("%s points at %s, which does not exist" % (key, path))
    if problems:
        error(
            "The %s module needs databases that are not ready:\n    %s\n"
            "Install and configure them with:\n    metawrap2 install-db --list\n"
            "and set their paths under [databases] in metawrap2.toml (metawrap2 install-db "
            "writes them there)." % (module, "\n    ".join(problems))
        )


def add_scorer_args(p, what: str = "bins") -> None:
    """``--checkm1``, ``--checkm2`` and ``--full-tree`` for a module that scores *what*.

    CheckM2 is the default. ``--checkm2`` selected CheckM2 before it became the default, and
    ``--quick`` added ``--reduced_tree`` before that became CheckM1's default; both are kept,
    hidden, so scripts still run (see :func:`scorer_notes`).
    """
    scorer = p.add_mutually_exclusive_group()
    scorer.add_argument(
        "--checkm1",
        dest="checkm2",
        action="store_false",
        help="score %s with CheckM1 instead of CheckM2, on CheckM1's reduced genome tree (about "
        "16 GB of RAM per call). Needs `metawrap2 install-db checkm`" % what,
    )
    scorer.add_argument(
        "--checkm2",
        dest="checkm2",
        action="store_true",
        help="score %s with CheckM2 (the default)" % what,
    )
    p.set_defaults(checkm2=True)
    p.add_argument(
        "--full-tree",
        action="store_true",
        help="with --checkm1, place bins in CheckM1's full genome tree instead of the reduced "
        "tree (about 40 GB of RAM per call)",
    )
    p.add_argument("--quick", action="store_true", help=argparse.SUPPRESS)


def scorer_notes(args) -> None:
    """One line each for scorer options that change nothing in this run."""
    if getattr(args, "quick", False):
        comm(
            "--quick is no longer needed: --checkm1 uses CheckM1's reduced tree by default "
            "(--full-tree for the full tree)."
        )
    if getattr(args, "full_tree", False) and args.checkm2:
        comm("--full-tree applies to CheckM1 only (--checkm1); CheckM2 runs without it.")


class _NoCheckM1(argparse.Action):
    """``--checkm1`` on a module that scores with CheckM2 only: stop with one sentence."""

    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(
            "%s scores bins with CheckM2 only, so --checkm1 is not available." % self.const
        )


def add_checkm2_only_args(p, module: str) -> None:
    """Hidden ``--checkm2`` (a no-op) and ``--checkm1`` (an error) for a CheckM2-only *module*.

    Both were options of the module before it dropped CheckM1, so old scripts get a clear message.
    """
    p.add_argument("--checkm2", action="store_true", help=argparse.SUPPRESS)
    p.set_defaults(checkm2=True)
    p.add_argument(
        "--checkm1",
        nargs=0,
        action=_NoCheckM1,
        const=module,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )


def _require_scorer_database(checkm2: bool, settings: Settings) -> None:
    """Fail at the start if the selected CheckM's database is missing.

    An unset CHECKM2DB only warns: CheckM2 also keeps a database path in its env (from ``checkm2
    database --download``), which a manual install uses. CheckM1's stored path is lost when its
    env is rebuilt (see :func:`metawrap2.checkm.checkm1_data_env`), so CheckM1 needs CHECKM_DB.
    """
    key, name, db = (
        ("CHECKM2DB", "CheckM2", "checkm2") if checkm2 else ("CHECKM_DB", "CheckM1", "checkm")
    )
    path = settings.db(key)
    if path and os.path.exists(path):
        return
    if checkm2 and not path:
        warning(
            "CHECKM2DB is not set, so CheckM2 will use the database path stored in its "
            "environment, if there is one. `metawrap2 install-db checkm2` installs the database "
            "and sets CHECKM2DB."
        )
        return
    problem = (
        "%s points at %s, which does not exist" % (key, path) if path else "%s is not set" % key
    )
    # A dry run warns and keeps printing the plan; a real run stops before any work is done.
    (warning if runner.dry_run else error)(
        "%s scores the bins, but its database is not ready: %s.\n"
        "Install it with:\n    metawrap2 install-db %s\n"
        "or set %s under [databases] in metawrap2.toml.%s"
        % (
            name,
            problem,
            db,
            key,
            "" if checkm2 else " Without --checkm1, bins are scored with CheckM2.",
        )
    )


def scorer_preflight(args, settings: Settings) -> Optional[str]:
    """Check the selected CheckM's env and database up front; return the CheckM2 env.

    The CheckM2 env is None under ``--checkm1`` or with conda envs off. CheckM1 runs in the
    bin_refinement env, which callers already resolve. Checked at the start because CheckM can run
    hours into a module (after SPAdes in reassemble_bins).
    """
    _require_scorer_database(args.checkm2, settings)
    if not args.checkm2:
        return None
    return optional_env_for("checkm2", settings, "CheckM2, which scores bins unless --checkm1")


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _describe_overrides(overrides: dict) -> str:
    """Format overrides as command-line flags, e.g. ``--max-threads 8, --max-gpus 0``."""
    return ", ".join(
        "--%s %s" % (key.replace("_", "-"), value) for key, value in sorted(overrides.items())
    )


def start_run(
    module: str, args, env: Optional[str], settings: Settings, inputs: List[str]
) -> RunRecorder:
    """Begin provenance recording for a module run and wire it into the runner.

    ``inputs`` are the paths this run reads (assembly/bins/reads); any that live inside a
    previous MetaWrap2 output dir let provenance propagate. Call :func:`finish_run` in a
    ``finally`` so even a failed run leaves a record.
    """
    output = getattr(args, "output", None)
    # Refuse to overwrite a previous MetaWrap2 run unless --force/--resume (or --dry-run).
    if (
        output
        and os.path.isfile(os.path.join(output, "run_config.json"))
        and not _cmd.runner.force
        and not _cmd.runner.resume
        and not _cmd.runner.dry_run
    ):
        error(
            "Output directory '%s' already contains a MetaWrap2 run. Re-run with "
            "--resume to continue it, --force to overwrite, or choose a new -o." % output
        )

    # A scheduler's SIGTERM (wall clock, memory limit) or a dropped SSH session's SIGHUP unwinds
    # like Ctrl-C, so the finally blocks run and the outcome is recorded.
    _interrupt.install_handlers()

    params = {k: v for k, v in vars(args).items()}
    inputs = [p for p in inputs if p]

    # The run has one resource budget, resolved and printed here so the user can check it before
    # a long run.
    budget = _limits.resolve(settings.limits)
    comm("Resource budget: %s" % budget.describe())
    overrides = _limits.overrides()
    if overrides:
        comm("Budget overridden on the command line: %s" % _describe_overrides(overrides))

    rec = RunRecorder(
        module=module,
        version=__version__,
        parameters=params,
        config_file=getattr(args, "config", None),
        conda_env=env,
        databases=dict(settings.databases),
        inputs=inputs,
        budget=budget.as_dict(),
    )
    _cmd.set_recorder(rec)

    # One writer per output directory. Taken before anything is written, so a second run fails
    # before it can overwrite the first one's intermediates.
    if output and not _cmd.runner.dry_run:
        lock = _runlock.RunLock(output, force=_cmd.runner.force)
        try:
            note = lock.acquire()
        except _runlock.LockHeld as exc:
            error(str(exc))
        if note:
            warning(note)
        global _active_lock
        _active_lock = lock

    # Tee all command output to run.stdout/run.stderr in the output dir (screen unaffected).
    if output:
        os.makedirs(output, exist_ok=True)
        # One combined, timestamped log of MetaWrap2's messages and all tool output, alongside
        # the raw per-stream files.
        _log.configure_run_log(os.path.join(output, _log.RUN_LOG_NAME))
        stdout_log = os.path.join(output, "run.stdout")
        stderr_log = os.path.join(output, "run.stderr")
        # The logs are appended to within a run (one section per command), and across runs only
        # with --resume. A fresh run truncates them so they describe only this run.
        if not _cmd.runner.resume:
            for path in (stdout_log, stderr_log):
                try:
                    open(path, "w").close()
                except OSError:
                    pass
        _cmd.set_run_logs(stdout_log, stderr_log)
        # Print the live log path once at start, so a slow step can be followed before it ends.
        comm(
            "Live log of this run: %s (MetaWrap2's messages and all tool output, written as "
            "the run goes). The raw tool streams are in run.stdout and run.stderr beside it."
            % os.path.join(output, _log.RUN_LOG_NAME)
        )
    _warn_odd_paths(inputs + ([output] if output else []))

    if output and not _cmd.runner.dry_run:
        # Record the environment before any tool runs, so a run that fails in its first command
        # still has it.
        rec.record_environment(output)

    # Monitor disk use for the whole run. check_disk_space only checks once before a step, and a
    # long step can still fill the filesystem partway through.
    if output and not _cmd.runner.dry_run:
        _start_disk_monitor(module, settings, output, budget)

    # Sample GPU memory over the whole run, since the peak is transient. One nvidia-smi call every
    # five seconds, only on machines with a driver. Cleared first so a second run in the same
    # process (a test suite, or `pathways` calling finish_run twice) does not inherit the peak.
    # Skipped under `run`, whose driver samples the step from outside.
    # The same thread checks CPU utilisation, so it also starts on a machine without a GPU when
    # the step has a thread grant worth watching.
    global _active_gpu_watcher
    _active_gpu_watcher = None
    if not _cmd.runner.dry_run and not os.environ.get(DRIVEN_BY_RUN_ENV):
        utilisation = _usage.configured_utilisation(getattr(args, "threads", None), module)
        if _usage.has_gpu() or utilisation is not None:
            _active_gpu_watcher = _usage.StepWatcher(os.getpid(), utilisation=utilisation).start()

    # A machine-wide run history outside any output directory. Written at start so an
    # interrupted or killed run is still recorded; finish_run appends the outcome.
    if not _cmd.runner.dry_run:
        _history.record(
            module,
            version=__version__,
            run_id=rec.run_id,
            command=rec.command_line,
            output=output,
            status=_history.STARTED,
            databases=settings.databases,
            conda_env=env,
        )
    return rec


def _start_disk_monitor(
    module: str, settings: Settings, output_dir: str, budget: _limits.Budget
) -> None:
    """Begin watching this run's filesystems, unless the space checks were turned off.

    Honours ``--skip-space-check``, like the preflight check.
    """
    if _cmd.runner.skip_space_check or getattr(settings, "skip_space_check", False):
        return

    def report(state: _limits.DiskState) -> None:
        warning(
            "Stopping %s: %s. The run is stopped so it does not fill the filesystem. Free some "
            "space (or set scratch_dir under [settings] to a filesystem with more space) and "
            "re-run with --resume. To allow a fuller disk, raise [limits] "
            "max_disk_fraction." % (module, state.describe())
        )

    global _active_disk_monitor
    _active_disk_monitor = _limits.DiskMonitor(
        [output_dir, _scratch.scratch_root(settings, output_dir)], budget, on_stop=report
    )
    _active_disk_monitor.start()


def _outcome_from_context() -> str:
    """The history status implied by how we are leaving the run.

    Called from a module's ``finally``, so ``sys.exc_info()`` still describes whatever is
    unwinding: nothing (success), a stop signal or Ctrl-C (interrupted), or anything else -
    including the ``SystemExit`` that ``error()`` raises - (failed).
    """
    exc = sys.exc_info()[1]
    if exc is None:
        return _history.COMPLETED
    if isinstance(exc, KeyboardInterrupt):  # covers interrupt.Stopped
        return _history.INTERRUPTED
    if isinstance(exc, SystemExit) and not exc.code:
        return _history.COMPLETED
    return _history.FAILED


def self_usage(rec: RunRecorder) -> Dict[str, Any]:
    """This process's resource use so far, in the form :func:`metawrap2.usage.summarise` reads.

    Matches what the ``run`` driver gets from ``os.wait4``: CPU times summed over RUSAGE_SELF and
    RUSAGE_CHILDREN, and peak RSS as the maximum of the two (a high-water mark, so not summed).
    Unit conversion and wording come from :mod:`metawrap2.usage`, shared with the driver.
    """
    me = resource.getrusage(resource.RUSAGE_SELF)
    kids = resource.getrusage(resource.RUSAGE_CHILDREN)
    measured: Dict[str, Any] = {
        "wall_seconds": round(rec.elapsed(), 1),
        "peak_rss_bytes": max(_usage.self_peak_rss_bytes(), _usage.maxrss_bytes(kids.ru_maxrss)),
        "user_seconds": round(me.ru_utime + kids.ru_utime, 1),
        "system_seconds": round(me.ru_stime + kids.ru_stime, 1),
    }
    # Only when this run used a GPU, so a CPU run on a GPU machine reports no VRAM, not zero.
    watcher = _active_gpu_watcher
    if watcher is not None and watcher.peak_mib:
        measured["peak_gpu_memory_bytes"] = watcher.peak_mib * _usage.MIB
        measured["gpu_devices"] = len(watcher.devices)
    return measured


def finish_run(
    rec: RunRecorder,
    output_dir: str,
    inputs: Optional[List[str]] = None,
    status: Optional[str] = None,
) -> None:
    """Write the provenance files (run_config.json, run_commands.txt, provenance.txt).

    Also closes out this run's history entry and releases the output directory's lock.

    *status* is normally left unset. Modules call this from a ``finally``, so the status is
    derived from the exception currently propagating, if any; a killed run or one that called
    ``error()`` is not recorded as completed.
    """
    if status is None:
        status = _outcome_from_context()
    try:
        rec.finalize(output_dir, input_paths=inputs, capture_versions=not _cmd.runner.dry_run)
    finally:
        # Print the run's cost line for a standalone module (under `run` the driver prints it),
        # e.g. to show "24 threads granted, 1.0 core busy". Printed whatever the outcome, and
        # before configure_run_log(None) so it also goes to metawrap2.log.
        global _cost_reported_for
        # Stop the GPU sampler before composing the line; it already holds the peak.
        if _active_gpu_watcher is not None:
            _active_gpu_watcher.stop()
        if not os.environ.get(DRIVEN_BY_RUN_ENV) and _cost_reported_for != rec.run_id:
            _cost_reported_for = rec.run_id
            comm("%s used %s" % (rec.module, _usage.summarise(self_usage(rec))))
        if not _cmd.runner.dry_run:
            _history.record(
                rec.module,
                version=rec.version,
                run_id=rec.run_id,
                command=rec.command_line,
                output=output_dir,
                status=status,
                duration_seconds=rec.elapsed(),
                conda_env=rec.conda_env,
            )
        _log.configure_run_log(None)
        _cmd.set_recorder(None)
        _cmd.set_run_logs(None, None)
        global _active_lock, _active_disk_monitor
        if _active_lock is not None:
            _active_lock.release()
            _active_lock = None
        # Stop the disk monitor explicitly so it cannot outlive the run and interrupt whatever
        # the process does next.
        if _active_disk_monitor is not None:
            _active_disk_monitor.stop()
            _active_disk_monitor = None


# Characters in a path that tend to break downstream shell pipelines / tools.
_ODD_PATH = re.compile(r"[\s'\"$`|;&<>()\\!*?]")


def _warn_odd_paths(paths: List[str]) -> None:
    for p in paths:
        if p and _ODD_PATH.search(p):
            warning(
                "Path contains characters that can break some tools (spaces or shell "
                "metacharacters): %s. Consider renaming to [A-Za-z0-9._/-]." % p
            )


def collect_read_pairs(reads: List[str]) -> List[tuple]:
    """From a positional list of fastq files, pair up mates into (sample, r1, r2).

    Wrapper over :func:`metawrap2.io.reads.collect_pairs`, which handles the filename conventions
    (``_1/_2``, ``_R1/_R2``, ``_R1_001/_R2_001``, ...). Raises ValueError if a mate is missing.
    Compressed files are fine; the aligners read .gz directly.
    """
    return _reads.collect_pairs(reads)


def collect_reads(reads: List[str]) -> List[tuple]:
    """From a positional list of fastq files, return (sample_name, path) for each."""
    return _reads.collect_singles(reads)


def read_layout(args) -> _reads.ReadSet:
    """Validate the read layout a module was given, honouring --single-end/--interleaved.

    Returns a :class:`metawrap2.io.reads.ReadSet`. Problems are reported here before any tool
    runs.
    """
    single = bool(getattr(args, "single_end", False))
    interleaved = bool(getattr(args, "interleaved", False))
    if single and interleaved:
        error("--single-end and --interleaved are mutually exclusive.")
    files = [f for f in (getattr(args, "reads_1", None), getattr(args, "reads_2", None)) if f]
    try:
        return _reads.classify(files, interleaved=interleaved, single=single)
    except ValueError as exc:
        error(str(exc))
        raise  # unreachable: error() exits, but this tells the type checker so


def resolve_threads(args, settings: Settings) -> int:
    """Final thread count for a run: an explicit ``-t`` wins, else the config, else the budget.

    ``-t`` defaults to ``None`` in every module so that an unset value differs from ``-t 1``.
    ``-t all`` means every core the budget allows. With neither flag nor config ``threads``, the
    budget is used, as for ``-m``.

    The ``[limits]`` thread ceiling is enforced here, since every module's thread count passes
    through this function. A request above the ceiling is capped, not refused.
    """
    requested = getattr(args, "threads", None)
    budget = _limits.resolve(settings.limits)
    if requested is None:
        threads = settings.threads or budget.threads
    else:
        threads = requested
    if threads < 1:
        error("--threads must be at least 1 (got %s)." % threads)
    if threads > budget.threads:
        # Warn only for an explicit -t. A capped config `threads` value is already covered by the
        # budget line start_run prints.
        if requested is not None:
            if threads > budget.machine_cores:
                warning(
                    "Asked for %d threads but this machine has %d cores; using %d."
                    % (threads, budget.machine_cores, budget.threads)
                )
            else:
                warning(
                    "Asked for %d threads but the budget allows %d of this machine's %d cores; "
                    "using %d. Raise [limits] max_thread_fraction or max_threads to use more."
                    % (threads, budget.threads, budget.machine_cores, budget.threads)
                )
        threads = budget.threads
    args.threads = threads
    return threads


def resolve_memory_gb(args, settings: Settings) -> int:
    """Final memory allowance in GB for a run: an explicit ``-m`` wins, else the budget.

    ``-m`` defaults to ``None`` so an unset value can be told apart, as for ``-t``. Unset, it
    resolves to the ``[limits]`` budget. (The earlier fixed per-module values, 24 GB for an
    assembly and 4 GB for binning, left large machines unused while metaSPAdes ran out of memory
    and CheckM ran pplacer single-threaded.)

    An explicit request above the ceiling is capped, matching :func:`resolve_threads`.
    """
    requested = getattr(args, "memory", None)
    budget = _limits.resolve(settings.limits)
    allowed = max(1, int(budget.memory_gb))
    if requested is None:
        args.memory = allowed
        return allowed
    if requested < 1:
        error("--memory must be at least 1 GB (got %s)." % requested)
    memory = requested
    if memory > allowed:
        warning(
            "Asked for %d GB but the budget allows %d GB of this machine's %s; using %d. "
            "Raise [limits] max_memory_fraction or max_memory_gb to use more."
            % (memory, allowed, _limits.human(budget.machine_memory_bytes), allowed)
        )
        memory = allowed
    args.memory = memory
    return memory


def memory_arg(parser) -> None:
    """Add the standard ``-m/--memory`` option (default: the ``[limits]`` memory budget)."""
    help_text = (
        "memory in GB the tools may use (default: the [limits] budget, "
        "%d%% of this machine's RAM)" % round(_limits.DEFAULT_MEMORY_FRACTION * 100)
    )
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=None,
        # argparse %-formats help strings, so a literal percent sign must be doubled.
        help=help_text.replace("%", "%%"),
    )


def threads_arg(parser) -> None:
    """Add the standard ``-t/--threads`` option (default: the config's threads setting)."""
    help_text = (
        "number of threads ('all' for every core the budget allows; default: [settings] "
        "threads in metawrap2.toml, else the [limits] budget, %d%% of this machine's cores)"
        % round(_limits.DEFAULT_THREAD_FRACTION * 100)
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=_threads_value,
        default=None,
        # See memory_arg: argparse %-formats help strings, so the percent sign must be doubled.
        help=help_text.replace("%", "%%"),
    )


def _threads_value(text: str) -> int:
    if str(text).strip().lower() == "all":
        return os.cpu_count() or 1
    try:
        return int(text)
    except ValueError:
        raise ValueError("threads must be an integer or 'all', got %r" % text)


# ── dry-run-aware filesystem helpers ─────────────────────────────────────────────────────
# Under --dry-run no command executes, so the files a module expects to read do not exist.
# These helpers skip filesystem effects and output checks under --dry-run, so the whole plan is
# printed instead of failing at the first step that inspects a tool's output.


def dry_run() -> bool:
    """True if this is a --dry-run (no commands are executed)."""
    return bool(_cmd.runner.dry_run)


def force() -> bool:
    """True if this run was told to overwrite an existing output directory (--force).

    A module that clears a previous run's artifacts must check this, because ``--resume`` reaches
    the same code with an existing output folder. Without the check, resuming bin_refinement after
    its refine step deleted the bin sets it was about to skip rebuilding.
    """
    return bool(_cmd.runner.force)


def expect_file(path: str, message: str, allow_empty: bool = False) -> bool:
    """Assert a tool produced *path*, or fail with *message*. Always passes under --dry-run."""
    if dry_run():
        return True
    if not os.path.isfile(path) or (not allow_empty and not os.path.getsize(path)):
        error(message)
    return True


def expect_dir(path: str, message: str) -> bool:
    """Assert a tool produced directory *path*. Always passes under --dry-run."""
    if dry_run():
        return True
    if not os.path.isdir(path):
        error(message)
    return True


def fs_move(src: str, dest: str) -> None:
    """os.replace, skipped under --dry-run (the source was never created)."""
    if dry_run():
        return
    os.replace(src, dest)


def fs_replace_dir(src: str, dest: str) -> None:
    """``os.replace``, but clearing *dest* first if it is a directory. Skipped under --dry-run.

    ``os.replace`` overwrites a file but fails on a non-empty directory, so moving a directory
    into place fails on every rerun into the same output folder (``--force`` or ``--resume``),
    after the work is done. *dest* is always an earlier run's version, so it is removed first.

    *src* may be a file or a directory, as in the loops that move bin sets and their .stats.tsv
    files together.
    """
    if dry_run():
        return
    if os.path.isdir(dest):
        shutil.rmtree(dest, ignore_errors=True)
    os.replace(src, dest)


def fs_copy(src: str, dest: str) -> None:
    """shutil.copy, skipped under --dry-run."""
    if dry_run():
        return
    shutil.copy(src, dest)


def fs_remove(path: str) -> None:
    """Remove *path* if it exists; skipped under --dry-run."""
    if dry_run():
        return
    if os.path.isfile(path):
        os.remove(path)


def file_size(path: str) -> int:
    """Size of *path*, or a nonzero placeholder under --dry-run so checks pass."""
    if dry_run():
        return 1
    return os.path.getsize(path) if os.path.exists(path) else 0


def listdir(path: str) -> List[str]:
    """os.listdir, or [] under --dry-run when the directory was never created."""
    if not os.path.isdir(path):
        if dry_run():
            return []
        return []
    return os.listdir(path)


def list_bins(path: str, extensions: Sequence[str] = FASTA_EXTENSIONS) -> List[str]:
    """Sorted filenames of the bins in *path* - :func:`fasta_files`, as names rather than paths.

    For modules that build an output name per bin (prokka's per-bin directory, the ``.orig.fa``
    copy, the per-bin read files). Using ``os.listdir`` instead made three of them treat the
    unbinned catch-all as a genome.

    Returns an empty list for a missing directory, since callers have already validated their
    inputs or are running ``--dry-run``.
    """
    return [os.path.basename(p) for p in fasta_files_in(path, extensions)]


#: Argument names that hold a filesystem path, and so should be made absolute up front.
_PATH_ARGS = (
    "output",
    "assembly",
    "bins",
    "bins_a",
    "bins_b",
    "bins_c",
    "reads_1",
    "reads_2",
    "nanopore",
    "bakta_db",
)


def absolutize_paths(args, extra: Optional[List[str]] = None) -> None:
    """Rewrite every path argument on *args* to an absolute path, in place.

    Some tools run with ``cwd=<output dir>`` because they write into the working directory
    (taxator-tk's ``binner``, the plotting helpers). With a relative ``-o`` this gave
    ``cd out && ... out/file``, i.e. ``out/out/file``. Absolute paths avoid that and make the
    recorded provenance commands unambiguous.
    """
    names = list(_PATH_ARGS) + list(extra or [])
    for name in names:
        value = getattr(args, name, None)
        if isinstance(value, str) and value:
            setattr(args, name, os.path.abspath(os.path.expanduser(value)))
    reads = getattr(args, "reads", None)
    if isinstance(reads, list):
        args.reads = [
            os.path.abspath(os.path.expanduser(r)) if isinstance(r, str) else r for r in reads
        ]
    seqs = getattr(args, "seqs", None)
    if isinstance(seqs, list):
        args.seqs = [
            os.path.abspath(os.path.expanduser(s)) if isinstance(s, str) else s for s in seqs
        ]


# ── input validation ─────────────────────────────────────────────────────────────────────
# Wrappers over metawrap2.validate that report through error() (a banner + exit) instead of
# raising, which is how modules report user-facing problems everywhere else.


def validate_inputs(context: str, checks: List[tuple]) -> None:
    """Run *checks* and abort with every problem listed if any fail.

    See :func:`metawrap2.validate.require_inputs`. Skipped under --dry-run, and under
    --skip-validation for the case where a check rejects a good file.
    """
    if dry_run():
        return
    if _cmd.runner.skip_validation:
        warning("Skipping input validation for %s (--skip-validation)." % context)
        return
    try:
        _validate.require_inputs(checks, context=context)
    except _validate.InputProblem as problem:
        error(str(problem))


def expect_produced(path: str, what: str, hint: str = "") -> None:
    """Assert a tool produced a usable *path*, not an empty or truncated one.

    ``expect_file`` only checks existence and non-emptiness; this also validates the format
    for FASTA/FASTQ, so an assembler that exits 0 having written a header and nothing else is
    caught here rather than in the next module.
    """
    if dry_run() or _cmd.runner.skip_validation:
        return
    lower = path.lower()
    if lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2")):
        problem = _validate.check_fastq(path, what)
    elif lower.endswith((".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fasta.gz")):
        problem = _validate.check_fasta(path, what)
    else:
        problem = _validate.require_file(path, what)
    if problem:
        message = (
            "%s\nThis file is an output of the step that just ran. The step failed "
            "without reporting an error." % problem
        )
        if hint:
            message += "\n" + hint
        error(message)


#: Re-exported so modules can build check lists without importing metawrap2.validate.
check_fastq = _validate.check_fastq
check_fasta = _validate.check_fasta
require_file = _validate.require_file
require_nonempty_dir = _validate.require_nonempty_dir


# ── scratch space and the disk-space preflight ───────────────────────────────────────────


def scratch_dir(settings: Settings, output_dir: str, name: str) -> str:
    """A named temporary directory for this run, honouring [settings] scratch_dir."""
    return _scratch.scratch_dir(settings, output_dir, name)


def check_disk_space(module: str, settings: Settings, inputs: List[str], output_dir: str) -> None:
    """Warn - or, if asked, refuse - when the filesystem looks too small for this step.

    metaSPAdes needs many times its input size in scratch. The estimate is a rough multiple of
    the input size (see metawrap2.scratch.SPACE_MULTIPLIERS), not a measurement, and can be wrong
    either way, so by default it only warns. Set ``strict_space_check = true`` under
    ``[settings]`` to make it a hard stop (useful on shared filesystems).

    ``--skip-space-check``, or ``skip_space_check = true``, turns the check off entirely.
    """
    if dry_run() or _cmd.runner.skip_space_check or getattr(settings, "skip_space_check", False):
        return
    targets = [output_dir, _scratch.scratch_root(settings, output_dir)]

    # The ceiling check is a measurement (the filesystem is already fuller than the budget
    # allows), unlike the estimate below, so it always stops the step.
    budget = _limits.resolve(settings.limits)
    under_ceiling, ceiling_message = _limits.check_ceiling(targets, budget)
    if not under_ceiling:
        error(
            "Not starting %s: %s\nThis is the [limits] max_disk_fraction ceiling, checked "
            "against measured free space. Free some space, set scratch_dir to a filesystem "
            "with more space, or raise max_disk_fraction." % (module, ceiling_message)
        )

    ok, message = _scratch.check_space(module, inputs, targets)
    if ok:
        comm("disk space check: %s" % message)
        return
    if getattr(settings, "strict_space_check", False):
        error(
            "Not enough free disk space to start %s:\n    %s\n"
            "The run stops because strict_space_check = true is set under [settings].\n"
            "Pass --skip-space-check to start anyway, or set scratch_dir to a larger filesystem."
            % (module, message.replace("\n", "\n    "))
        )
    warning(
        "%s may not have room to finish:\n    %s\n"
        "This is estimated from the input size and may be wrong; the run continues. "
        "To put temporary files on another filesystem, set scratch_dir under [settings]. "
        "To stop the run in this case, set strict_space_check = true."
        % (module, message.replace("\n", "\n    "))
    )


# --- aligning reads to a reference ---------------------------------------------------------

#: The aligner streams into ``samtools sort`` so the SAM (about 25 GB per sample) never reaches
#: disk. ``pipefail`` because bash otherwise reports only the sort's status, and a failed alignment
#: would give an empty BAM with exit 0.
ALIGN_SORT = "set -o pipefail; {align} | samtools sort -T {tmp} -@ {threads} -O BAM -o {bam} -"

#: Threads given to ``samtools sort``. It runs alongside the aligner and saturates well before
#: it; the benchmark that chose minimap2 used 8.
SORT_THREADS = 8


# Read subsampling for coverage was tried and removed. On CAMI2 GI with metaBAT2 (AMBER), half
# the reads kept 33 of 37 high-quality MAGs for 53% of the step's CPU, a quarter kept 25, a tenth
# kept 12. Binning is under 2% of a default run, so halving it saves about 1% of the pipeline and
# loses 11% of the genomes. Assembly is 86% of the run; thin reads there if anywhere.
#
# Read capping is used in reassemble_bins --max-coverage, which thins each bin's recruited reads to
# 50x of the bin's assembled length, since an abundant organism can give SPAdes hundreds of times
# the coverage it needs.


def _index_is_current(index: str, reference: str, aligner: str, command: str = "") -> bool:
    """Whether the index of *reference* can be reused for *command* without rebuilding.

    True only when every index file for the aligner exists, is non-empty and is at least as new as
    the reference, and the build marker records the same index command. The marker is written only
    after a build finishes, so an index that another process is still writing, or one built with a
    different minimap2 preset, is not reused. The file list comes from INDEX_FILES (bwa writes five
    files under one prefix, minimap2 writes one); a partial set counts as absent.
    """
    paths = [reference + suffix for suffix in _aligners.INDEX_FILES.get(aligner, ())]
    if not paths:
        return False
    try:
        reference_mtime = os.path.getmtime(reference)
        for path in paths:
            if os.path.getsize(path) == 0 or os.path.getmtime(path) < reference_mtime:
                return False
        with open(reference + _aligners.INDEX_MARKER) as fh:
            return fh.read().strip() == command.strip()
    except OSError:
        return False


def _build_or_reuse_index(
    index: str, reference: str, aligner: str, command: str, env: Optional[str]
) -> None:
    """Reuse the index of *reference* if it is current for *command*, else build it.

    A lock file serialises processes that map against the same reference: the first builds the
    index and writes the marker, and the others find it current and reuse it. Without the lock,
    parallel mappings would each rebuild the index while others read it.
    """
    marker = reference + _aligners.INDEX_MARKER
    with open(marker + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if _index_is_current(index, reference, aligner, command):
            comm("reusing the existing %s index of %s" % (aligner, os.path.basename(reference)))
            return
        if os.path.exists(marker):
            os.remove(marker)
        run(command, env=env, tool="%s index" % aligner)
        with open(marker, "w") as fh:
            fh.write(command + "\n")


def align_read_sets(
    reference: str,
    reads: Sequence[str],
    out_dir: str,
    *,
    aligner: str,
    threads: int,
    env: Optional[str] = None,
    single_end: bool = False,
    interleaved: bool = False,
    long_read_type: Optional[str] = None,
) -> List[str]:
    """Map each read set in *reads* to *reference*, leaving one sorted BAM per sample in *out_dir*.

    Returns the BAM paths in sample order. Shared by `binning` (every sample against one assembly)
    and the cohort coverage stage (every sample against the concatenated catalogue).

    Existing BAMs are kept, so an interrupted pass resumes at the sample it stopped on.

    """
    spec = _aligners.ALIGNERS[aligner]
    if long_read_type:
        # Long reads need minimap2's platform preset instead of `sr`, which maps a 12 kb nanopore
        # read poorly and gives low, uneven coverage. The preset comes from metawrap2.longreads.
        from .. import longreads as _long

        if aligner != "minimap2":
            raise ValueError(
                "long reads need --aligner minimap2: bwa and bwa-mem2 have no long-read mode, and "
                "coverage from long reads mapped with them is not usable"
            )
        preset = _long.resolve(long_read_type).minimap2_preset
        spec = dict(spec)
        spec["single"] = spec["single"].replace("-ax sr", "-ax %s" % preset)
        spec["paired"] = spec["single"]
        spec["interleaved"] = spec["single"]
    index = _aligners.index_path(reference, aligner)
    index_cmd = spec["index_cmd"]
    if long_read_type and aligner == "minimap2":
        # The index stores the preset's k-mer and window size; an `sr` index used for `map-ont`
        # mapping gives no error, so build it with the matching preset.
        from .. import longreads as _long

        index_cmd = index_cmd.replace(
            "-x sr", "-x %s" % _long.resolve(long_read_type).minimap2_preset
        )
    command = index_cmd.format(assembly=reference, index=index)
    if _cmd.runner.dry_run:
        run(command, env=env, tool="%s index" % aligner)
    else:
        _build_or_reuse_index(index, reference, aligner, command, env)
    tmp = os.path.join(out_dir, "tmp-samtools")
    sort_threads = min(threads, SORT_THREADS)
    sets = (
        collect_reads(list(reads))
        if (single_end or interleaved)
        else collect_read_pairs(list(reads))
    )

    bams = []
    for entry in sets:
        sample = entry[0]
        bam = os.path.join(out_dir, sample + ".bam")
        bams.append(bam)
        if not dry_run() and os.path.isfile(bam):
            comm("skipping aligning %s reads; %s already exists" % (sample, bam))
            continue
        comm("Aligning and sorting %s reads back to assembly" % sample)
        if single_end:
            align = spec["single"].format(threads=threads, index=index, reads=entry[1])
        elif interleaved:
            align = spec["interleaved"].format(threads=threads, index=index, reads=entry[1])
        else:
            _s, r1, r2 = entry
            align = spec["paired"].format(threads=threads, index=index, r1=r1, r2=r2)
        run(
            ["bash", "-c", ALIGN_SORT.format(align=align, tmp=tmp, threads=sort_threads, bam=bam)],
            env=env,
            tool="%s | samtools sort" % aligner,
        )
    return bams
