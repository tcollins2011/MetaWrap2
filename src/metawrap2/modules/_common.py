"""Shared helpers for MetaWrap2 modules.

Deliberately tiny. Modules import the banner helpers, the conda-env lookup, and a couple of
filesystem/read utilities from here; everything else stays inline in the module so a reader
sees exactly what runs. This is not a framework — just the few things every module repeats.
"""

from __future__ import annotations

import os
import re
import sys
from typing import List, Optional

from .. import __version__
from .. import command as _cmd
from .. import history as _history
from .. import interrupt as _interrupt
from .. import logging as _log
from .. import runlock as _runlock
from .. import scratch as _scratch
from .. import validate as _validate
from ..checkpoint import make_checkpoint
from ..command import ToolError, run  # re-exported for modules
from ..config import MODULE_ENVS, Settings, conda_env_exists
from ..io import reads as _reads
from ..provenance import RunRecorder

#: The lock held on the current run's output directory, if any. Deliberately not a field on
#: RunRecorder: everything on the recorder is serialised into run_config.json, and a RunLock is
#: not JSON-serialisable.
_active_lock: Optional[_runlock.RunLock] = None

# Banner logging (same look as the original comm/announcement/warning/error).
comm = _log.comm
announcement = _log.announcement
warning = _log.warning
error = _log.error

__all__ = [
    "ToolError",
    "absolutize_paths",
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
    "fs_copy",
    "fs_move",
    "fs_remove",
    "listdir",
    "make_checkpoint",
    "read_layout",
    "require_databases",
    "require_file",
    "require_nonempty_dir",
    "resolve_threads",
    "run",
    "start_run",
    "threads_arg",
    "validate_inputs",
    "warning",
]


def env_for(module: str, settings: Settings) -> Optional[str]:
    """Conda env name for *module*, or None if the user disabled conda envs.

    Fails fast with the exact command to fix it when the env has not been created. Without
    this the first tool invocation dies inside ``mamba run`` with a message about a missing
    prefix, which does not tell the user what to do about it - and in a long pipeline that
    can happen well into a run.
    """
    if not settings.use_conda_envs:
        return None
    env = MODULE_ENVS.get(module)
    if env and not conda_env_exists(env):
        error(
            "The conda environment '%s' for the %s module does not exist.\n"
            "Create it with:\n"
            "    metawrap2 install-env %s\n"
            "or run the module's tools straight off your PATH by setting "
            "use_conda_envs = false under [settings] in metawrap2.toml." % (env, module, module)
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
            "then set them under [databases] in metawrap2.toml (metawrap2 install-db "
            "writes them for you)." % (module, "\n    ".join(problems))
        )


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


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

    # A scheduler's SIGTERM (wall clock, memory limit) or a dropped SSH session's SIGHUP now
    # unwinds like Ctrl-C, so the finally blocks below run and the outcome is recorded.
    _interrupt.install_handlers()

    params = {k: v for k, v in vars(args).items()}
    inputs = [p for p in inputs if p]
    rec = RunRecorder(
        module=module,
        version=__version__,
        parameters=params,
        config_file=getattr(args, "config", None),
        conda_env=env,
        databases=dict(settings.databases),
        inputs=inputs,
    )
    _cmd.set_recorder(rec)

    # One writer per output directory. Taken before anything is written, so a second run fails
    # before it can overwrite the first one's intermediates rather than halfway through.
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
        # One combined, timestamped log of everything - MetaWrap2's own messages and every
        # line of every tool - alongside the raw per-stream files.
        _log.configure_run_log(os.path.join(output, _log.RUN_LOG_NAME))
        stdout_log = os.path.join(output, "run.stdout")
        stderr_log = os.path.join(output, "run.stderr")
        # The logs are appended to within a run (one section per command). Across runs that
        # is only wanted when continuing one: a fresh run (--force, or a new output dir) used
        # to silently accumulate the previous run's output in the same file, so the log no
        # longer described the run that produced the data next to it.
        if not _cmd.runner.resume:
            for path in (stdout_log, stderr_log):
                try:
                    open(path, "w").close()
                except OSError:
                    pass
        _cmd.set_run_logs(stdout_log, stderr_log)
    _warn_odd_paths(inputs + ([output] if output else []))

    if output and not _cmd.runner.dry_run:
        # What this run is about to use, written before any tool runs - so a run that dies in
        # its first command still says what it was running with.
        rec.record_environment(output)

    # A durable, machine-wide record of the run, outside any output directory. Written at the
    # start so an interrupted or killed run still leaves a trace; finish_run appends the outcome.
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


def finish_run(
    rec: RunRecorder,
    output_dir: str,
    inputs: Optional[List[str]] = None,
    status: Optional[str] = None,
) -> None:
    """Write the provenance files (run_config.json, run_commands.txt, provenance.txt).

    Also closes out this run's history entry and releases the output directory's lock.

    *status* is normally left unset. Modules call this from a ``finally``, where the happy path
    and the two unhappy ones are indistinguishable from the arguments - so the status is derived
    from whether an exception is currently propagating, and which one. Defaulting to "completed"
    (as this used to) meant a run that was killed, or that called ``error()``, still recorded
    success in its history and told you nothing.
    """
    if status is None:
        status = _outcome_from_context()
    try:
        rec.finalize(output_dir, input_paths=inputs, capture_versions=not _cmd.runner.dry_run)
    finally:
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
        global _active_lock
        if _active_lock is not None:
            _active_lock.release()
            _active_lock = None


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

    Thin wrapper over :func:`metawrap2.io.reads.collect_pairs`, which owns the filename
    conventions (``_1/_2``, ``_R1/_R2``, ``_R1_001/_R2_001``, ...). Raises ValueError with a
    clear message if a mate is missing. Compression is fine - the aligners read .gz directly.
    """
    return _reads.collect_pairs(reads)


def collect_reads(reads: List[str]) -> List[tuple]:
    """From a positional list of fastq files, return (sample_name, path) for each."""
    return _reads.collect_singles(reads)


def read_layout(args) -> _reads.ReadSet:
    """Validate the read layout a module was given, honouring --single-end/--interleaved.

    Returns a :class:`metawrap2.io.reads.ReadSet`. Failures are reported here, up front, with
    an actionable message - not hours later from inside an assembler.
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
    """Final thread count for a run: an explicit ``-t`` wins, else the config, else 1.

    ``-t`` defaults to ``None`` in every module so that "user did not say" is
    distinguishable from "user asked for 1". Accepts ``-t all`` for every core on the
    machine. Without this, a user who set ``threads`` in metawrap2.toml still silently got
    single-threaded assemblies and alignments, which on a metagenome is the difference
    between hours and days.
    """
    requested = getattr(args, "threads", None)
    if requested is None:
        threads = settings.threads or 1
    else:
        threads = requested
    if threads < 1:
        error("--threads must be at least 1 (got %s)." % threads)
    available = os.cpu_count() or 1
    if threads > available:
        warning(
            "Asked for %d threads but this machine has %d cores; using %d."
            % (threads, available, available)
        )
        threads = available
    args.threads = threads
    return threads


def threads_arg(parser) -> None:
    """Add the standard ``-t/--threads`` option (default: the config's threads setting)."""
    parser.add_argument(
        "-t",
        "--threads",
        type=_threads_value,
        default=None,
        help="number of threads ('all' for every core; "
        "default: [settings] threads in metawrap2.toml, else 1)",
    )


def _threads_value(text: str) -> int:
    if str(text).strip().lower() == "all":
        return os.cpu_count() or 1
    try:
        return int(text)
    except ValueError:
        raise ValueError("threads must be an integer or 'all', got %r" % text)


# ── dry-run-aware filesystem helpers ─────────────────────────────────────────────────────
# Under --dry-run no command actually executes, so none of the files a module expects to
# consume exist. Using os.replace / os.path.getsize / bare `error(...)` directly meant every
# module crashed with a raw traceback (or a misleading "something went wrong with X") at the
# first step that inspected a tool's output - so --dry-run only ever showed the first command
# or two of a pipeline, which defeats the point of the flag. These helpers make the
# orchestration skip filesystem effects and post-condition checks when nothing really ran.


def dry_run() -> bool:
    """True if this is a --dry-run (no commands are actually executed)."""
    return bool(_cmd.runner.dry_run)


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


def fs_copy(src: str, dest: str) -> None:
    """shutil.copy, skipped under --dry-run."""
    if dry_run():
        return
    import shutil as _shutil

    _shutil.copy(src, dest)


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

    Some tools have to be run with ``cwd=<output dir>`` because they write into the working
    directory (taxator-tk's ``binner``, the plotting helpers). Combining that with paths built
    from a *relative* ``-o`` silently produced ``cd out && ... out/file`` - i.e. ``out/out/file``
    - and the step failed for everyone who passed a relative output directory, which is what
    the tutorial itself does. Absolute paths make the working directory irrelevant, and also
    make the recorded provenance commands unambiguous.
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

    See :func:`metawrap2.validate.require_inputs`. Skipped under --dry-run, where the point is
    to print a plan rather than to touch the data, and under --skip-validation, the escape hatch
    for the rare case where a check wrongly rejects a good file.
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
    """Assert a tool actually produced a usable *path*, not an empty or truncated one.

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
            "%s\nThis file is an output of the step that just ran, so the step failed "
            "even though it did not report an error." % problem
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

    metaSPAdes on a real library needs many times its input size in scratch, and discovering that
    by filling the disk several hours in is expensive. But the estimate is a rough multiple of the
    input size (see metawrap2.scratch.SPACE_MULTIPLIERS), chosen rather than measured, and it
    cannot know how large scratch gets at its peak, because scratch is deleted when a step
    finishes. It is therefore capable of being wrong in both directions.

    Given that, the default is to say so loudly and carry on. Refusing to start on the strength of
    a guess is the worse of the two failures: the user can see the numbers, and knows things about
    their filesystem that this does not. Set ``strict_space_check = true`` under ``[settings]`` to
    make it a hard stop instead - worth doing on shared infrastructure, where one run filling the
    disk is everyone's problem.

    ``--skip-space-check``, or ``skip_space_check = true``, turns the check off entirely.
    """
    if dry_run() or _cmd.runner.skip_space_check or getattr(settings, "skip_space_check", False):
        return
    targets = [output_dir, _scratch.scratch_root(settings, output_dir)]
    ok, message = _scratch.check_space(module, inputs, targets)
    if ok:
        comm("disk space check: %s" % message)
        return
    if getattr(settings, "strict_space_check", False):
        error(
            "Not enough free disk space to start %s:\n    %s\n"
            "This is a hard stop because strict_space_check = true is set under [settings].\n"
            "Pass --skip-space-check to start anyway, or set scratch_dir somewhere roomier."
            % (module, message.replace("\n", "\n    "))
        )
    warning(
        "%s may not have room to finish:\n    %s\n"
        "That is an estimate from the input size, not a measurement, so it can be wrong - "
        "carrying on anyway. To put temporary files somewhere roomier, set scratch_dir under "
        "[settings]; to make this stop the run in future, set strict_space_check = true."
        % (module, message.replace("\n", "\n    "))
    )
