"""Submitting study steps as Slurm jobs, and following them to the end.

``metawrap2 run --executor slurm`` (and ``cohort``) submits each planned step as one batch job.
The job runs the same ``metawrap2 <module>`` command as the local executor, inside the CPUs,
memory and GPUs Slurm grants it. The driver keeps stage order: it waits for every job of a stage
before it submits the next stage. See ``run_slurm_step`` in :mod:`metawrap2.commands.run`.

- :func:`sbatch_command` builds the sbatch arguments and :func:`batch_script` the script, so
  ``--dry-run`` can print both without writing or submitting anything.
- :func:`submit_job` returns the job ID at once; the driver records it in the manifest before it
  waits, so ``--resume`` can reconnect to a job that outlived the driver.
- :func:`wait_for_job` polls ``squeue`` while the job is queued or running, then reads the final
  state, exit code and resource use from ``sacct``.

The batch script is exported the driver's environment (``--export=ALL``), so conda, ``PATH`` and
the ``METAWRAP2_*`` variables are the same on the compute node. The inputs, databases, Python and
output directory must be on a filesystem the compute nodes can read.
"""

from __future__ import annotations

import getpass
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from .cmdtemplate import Raw, fill_command, quoted_words

SBATCH = "sbatch"
SQUEUE = "squeue"
SACCT = "sacct"
SCANCEL = "scancel"

#: Seconds between state checks of one job. Each check is one squeue call.
POLL_SECONDS = 30
#: How long a job may be gone from squeue before sacct must have its final state. The accounting
#: database can lag the controller by a few minutes on a busy cluster.
ACCOUNTING_TIMEOUT_SECONDS = 600
#: How long squeue or sacct may keep failing (a controller restart, a timeout) before the driver
#: gives up on a job. The job itself is left alone.
QUERY_FAILURE_TIMEOUT_SECONDS = 1800
#: Seconds to wait for one Slurm client command.
COMMAND_TIMEOUT_SECONDS = 120

#: Job states after which the job will not run again. COMPLETED with exit code 0 is the only
#: success.
TERMINAL_STATES = frozenset(
    {
        "COMPLETED",
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
        "REVOKED",
        "SPECIAL_EXIT",
    }
)

#: A preempted job is requeued if the partition allows it, and squeue can miss it for a moment in
#: between. PREEMPTED counts as final only when two checks in a row report it.
_CONFIRM_TWICE = frozenset({"PREEMPTED"})

#: ``sbatch --parsable`` prints ``jobid`` or ``jobid;cluster``.
_JOB_ID = re.compile(r"^([0-9]+)(?:;\S+)?$")

#: The batch script. The command is one ``exec`` line, so the script is easy to read and rerun by
#: hand, and the job's exit code is the module's.
BATCH_SCRIPT = """\
#!/usr/bin/env bash
# MetaWrap2 study step {step}, submitted by `metawrap2 run --executor slurm`.
set -eu
cd {cwd}
export OMP_NUM_THREADS={threads}
exec {command}
"""


class SlurmError(RuntimeError):
    """A Slurm command failed, or a job's outcome cannot be determined."""


class SubmissionRejected(SlurmError):
    """sbatch refused the job, so nothing was submitted."""


class JobInterrupted(Exception):
    """The driver was interrupted while waiting for a job."""


@dataclass(frozen=True)
class SlurmJob:
    """The scheduling options of one job (resources aside)."""

    partition: Optional[str] = None
    account: Optional[str] = None
    qos: Optional[str] = None
    constraint: Optional[str] = None
    #: Unset leaves the partition's default time limit.
    time: Optional[str] = None
    gres: Optional[str] = None
    gpus: Optional[str] = None
    extra_args: Tuple[str, ...] = ()
    poll_seconds: int = POLL_SECONDS

    @property
    def requests_gpu(self) -> bool:
        return bool(self.gpus) or bool(self.gres and self.gres.split(":", 1)[0] == "gpu")


@dataclass
class SlurmResult:
    """The final state of one job."""

    job_id: str
    state: str
    returncode: int
    #: ``wall_seconds``, ``user_seconds``, ``system_seconds`` and ``peak_rss_bytes`` from sacct,
    #: in the form :mod:`metawrap2.usage` uses; a key is absent when sacct did not report it.
    usage: Dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.state == "COMPLETED" and self.returncode == 0

    def describe(self) -> str:
        return "Slurm job %s %s (exit code %d)" % (self.job_id, self.state, self.returncode)


# --- running the client commands ----------------------------------------------------------


def _run(
    argv: Sequence[str], timeout: int = COMMAND_TIMEOUT_SECONDS
) -> subprocess.CompletedProcess[str]:
    """Run a Slurm client command and capture its output."""
    try:
        return subprocess.run(
            list(argv), capture_output=True, text=True, check=False, timeout=timeout
        )
    except subprocess.TimeoutExpired as exc:
        raise SlurmError("%s did not answer within %d s" % (argv[0], timeout)) from exc
    except OSError as exc:
        raise SlurmError("could not run %s: %s" % (argv[0], exc)) from exc


def missing_commands(names: Sequence[str] = (SBATCH, SQUEUE, SACCT, SCANCEL)) -> List[str]:
    """The Slurm client commands that are not on PATH."""
    return [name for name in names if shutil.which(name) is None]


# --- the job script and the sbatch command -------------------------------------------------


def batch_script(command: Sequence[str], step: str, cwd: str, threads: int) -> str:
    """The batch script that runs *command* in *cwd*, every argument shell-quoted."""
    if not command:
        raise ValueError("a Slurm job needs a command")
    if threads < 1:
        raise ValueError("a Slurm job needs at least one thread")
    return fill_command(
        BATCH_SCRIPT,
        step=Raw(step.replace("\n", " ")),
        cwd=cwd,
        threads=threads,
        command=quoted_words(command),
    )


def write_script(path: str, contents: str) -> str:
    """Write *contents* to *path* atomically, as an executable script. Returns the path."""
    script = Path(path).resolve()
    script.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".%s." % script.name, dir=str(script.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(contents)
        os.chmod(temporary, 0o755)
        os.replace(temporary, script)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return str(script)


def _filename_pattern(path: str) -> str:
    """*path* for sbatch ``--output``, where ``%`` starts a replacement such as ``%j``."""
    return path.replace("%", "%%")


def sbatch_command(
    script: str,
    *,
    name: str,
    threads: int,
    memory_gb: int,
    log: str,
    job: SlurmJob,
) -> List[str]:
    """The sbatch arguments that submit *script*. The script comes last: sbatch passes anything
    after it to the script, not to the scheduler."""
    if threads < 1 or memory_gb < 1:
        raise ValueError("a Slurm job needs at least one thread and 1 GB of memory")
    if job.gres and job.gpus:
        raise ValueError("set gres or gpus for a Slurm job, not both")
    argv = [
        SBATCH,
        "--parsable",
        "--job-name=%s" % name,
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=%d" % threads,
        "--mem=%dG" % memory_gb,
        "--export=ALL",
        # stdout and stderr both go to the step's log, after the command line the driver wrote.
        "--output=%s" % _filename_pattern(log),
        "--open-mode=append",
    ]
    for option, value in (
        ("--time", job.time),
        ("--partition", job.partition),
        ("--account", job.account),
        ("--qos", job.qos),
        ("--constraint", job.constraint),
        ("--gres", job.gres),
        ("--gpus", job.gpus),
    ):
        if value:
            argv.append("%s=%s" % (option, value))
    argv += list(job.extra_args)
    argv.append(script)
    return argv


def submit_job(sbatch: Sequence[str]) -> str:
    """Run the *sbatch* command and return the job ID.

    Raises :class:`SubmissionRejected` when sbatch refused the job, and :class:`SlurmError` when
    it is not known whether the job was submitted (sbatch timed out).
    """
    try:
        result = subprocess.run(
            list(sbatch),
            capture_output=True,
            text=True,
            check=False,
            timeout=COMMAND_TIMEOUT_SECONDS,
        )
    except OSError as exc:
        raise SubmissionRejected("could not run sbatch: %s" % exc) from exc
    except subprocess.TimeoutExpired as exc:
        raise SlurmError(
            "sbatch did not answer within %d s; the job may or may not have been submitted"
            % COMMAND_TIMEOUT_SECONDS
        ) from exc
    if result.returncode != 0:
        raise SubmissionRejected(
            "sbatch refused the job: %s" % (result.stderr.strip() or result.stdout.strip())
        )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    match = _JOB_ID.match(lines[-1]) if lines else None
    if match is None:
        raise SlurmError("cannot read a job ID from sbatch's output: %r" % result.stdout.strip())
    return match.group(1)


# --- job state ------------------------------------------------------------------------------


def _state(value: str) -> str:
    """A job state without its modifiers: ``CANCELLED by 123`` and ``CANCELLED+`` are CANCELLED."""
    words = value.strip().split()
    return words[0].rstrip("+").upper() if words else ""


def _exit_code(value: str) -> int:
    """sacct's ``code:signal`` as one exit code, 128 + signal for a job killed by a signal."""
    try:
        code, signal = (int(part) for part in value.strip().split(":", 1))
    except ValueError:
        return 1
    if code == 0 and signal:
        return 128 + signal
    return code


def _cpu_seconds(value: str) -> Optional[float]:
    """sacct's ``[DD-[HH:]]MM:SS[.mmm]`` as seconds, or None when empty."""
    value = value.strip()
    if not value:
        return None
    days = 0
    if "-" in value:
        day_part, value = value.split("-", 1)
        days = int(day_part)
    parts = [float(p) for p in value.split(":")]
    seconds = 0.0
    for part in parts:
        seconds = seconds * 60 + part
    return days * 86400 + seconds


_SIZE_UNITS = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}


def _bytes(value: str) -> Optional[int]:
    """sacct's MaxRSS (``1234K``, ``5.5G``) as bytes, or None when empty."""
    match = re.match(r"^([0-9.]+)([KMGT]?)$", value.strip().upper())
    if not match:
        return None
    return int(float(match.group(1)) * _SIZE_UNITS[match.group(2)])


def queue_state(job_id: str) -> Optional[str]:
    """The job's state in squeue, or None when squeue no longer lists it."""
    result = _run([SQUEUE, "--noheader", "--jobs=%s" % job_id, "--format=%i|%T"])
    if result.returncode != 0:
        # Older releases fail on a job that has left the queue; newer ones print nothing.
        if "invalid job id" in (result.stderr or "").lower():
            return None
        raise SlurmError("squeue failed for job %s: %s" % (job_id, result.stderr.strip()))
    for line in result.stdout.splitlines():
        fields = line.strip().split("|")
        if len(fields) >= 2 and fields[0] == job_id:
            return _state(fields[1])
    return None


SACCT_FORMAT = "JobIDRaw,State,ExitCode,ElapsedRaw,UserCPU,SystemCPU,MaxRSS"


def accounting_result(job_id: str) -> Optional[SlurmResult]:
    """The job's final state from sacct, or None when sacct has no final state for it yet.

    The job's own line gives the state and exit code. Its steps (``.batch``, ``.extern``) give
    the peak memory, and a failed step fails the job even when the job line says COMPLETED.
    """
    result = _run(
        [SACCT, "--noheader", "--parsable2", "--jobs=%s" % job_id, "--format=" + SACCT_FORMAT]
    )
    if result.returncode != 0:
        raise SlurmError("sacct failed for job %s: %s" % (job_id, result.stderr.strip()))
    job: Optional[List[str]] = None
    steps: List[List[str]] = []
    for line in result.stdout.splitlines():
        fields = line.strip().split("|")
        if len(fields) < 7:
            continue
        if fields[0] == job_id:
            job = fields
        elif fields[0].startswith(job_id + "."):
            steps.append(fields)
    if job is None:
        return None
    state = _state(job[1])
    if state not in TERMINAL_STATES:
        return None
    code = _exit_code(job[2])
    if state == "COMPLETED":
        for fields in steps:
            step_state = _state(fields[1])
            if step_state in TERMINAL_STATES and step_state != "COMPLETED":
                state = step_state
                break
    if state != "COMPLETED" and code == 0:
        code = 1
    usage: Dict[str, Any] = {}
    if job[3].strip().isdigit():
        usage["wall_seconds"] = float(job[3])
    user, system = _cpu_seconds(job[4]), _cpu_seconds(job[5])
    if user is not None and system is not None:
        usage["user_seconds"] = round(user, 1)
        usage["system_seconds"] = round(system, 1)
    peaks = [b for b in (_bytes(fields[6]) for fields in [job] + steps) if b]
    if peaks:
        usage["peak_rss_bytes"] = max(peaks)
    return SlurmResult(job_id, state, code, usage)


def wait_for_job(
    job_id: str,
    poll_seconds: int = POLL_SECONDS,
    stop: Optional[threading.Event] = None,
    on_state: Optional[Callable[[str], None]] = None,
) -> SlurmResult:
    """Wait until job *job_id* has finished, and return its final state.

    *on_state* is called with each new queue state (PENDING, RUNNING, ...). Raises
    :class:`JobInterrupted` as soon as *stop* is set, and :class:`SlurmError` when the outcome
    cannot be determined; the job is not cancelled in either case.
    """
    if poll_seconds < 1:
        raise ValueError("the Slurm poll interval must be at least 1 s")
    gone_since: Optional[float] = None
    failing_since: Optional[float] = None
    last_state: Optional[str] = None
    unconfirmed: Optional[str] = None
    while True:
        if stop is not None and stop.is_set():
            raise JobInterrupted(job_id)
        try:
            state = queue_state(job_id)
            final = None
            if state is None or state in TERMINAL_STATES:
                final = accounting_result(job_id)
            failing_since = None
        except SlurmError as exc:
            # The controller or the accounting database is busy or restarting. Keep trying for
            # a while: the job is unaffected.
            now = time.monotonic()
            failing_since = failing_since or now
            if now - failing_since >= QUERY_FAILURE_TIMEOUT_SECONDS:
                raise SlurmError(
                    "could not query job %s for %d s: %s"
                    % (job_id, QUERY_FAILURE_TIMEOUT_SECONDS, exc)
                ) from exc
            state, final = last_state, None
        else:
            if state is not None and state != last_state and on_state is not None:
                on_state(state)
            last_state = state or last_state
            if final is not None:
                if final.state not in _CONFIRM_TWICE or unconfirmed == final.state:
                    return final
                unconfirmed = final.state
            else:
                unconfirmed = None
            if state is None and final is None:
                now = time.monotonic()
                gone_since = gone_since or now
                if now - gone_since >= ACCOUNTING_TIMEOUT_SECONDS:
                    raise SlurmError(
                        "job %s has left squeue, but sacct has no final state for it after %d s. "
                        "Check `sacct -j %s` (job accounting must be enabled)."
                        % (job_id, ACCOUNTING_TIMEOUT_SECONDS, job_id)
                    )
            else:
                gone_since = None
        if stop is not None:
            if stop.wait(poll_seconds):
                raise JobInterrupted(job_id)
        else:
            time.sleep(poll_seconds)


def job_finished(job_id: str) -> bool:
    """Whether the job is known to have ended: squeue lists it as ended or not at all."""
    state = queue_state(job_id)
    return state is None or state in TERMINAL_STATES


def find_job(name: str, since: Optional[str] = None) -> Optional[str]:
    """The ID of this user's most recent job called *name*, or None.

    Used on ``--resume`` when the driver stopped after sbatch but before it recorded the ID.
    squeue is checked first, then sacct from *since* (``YYYY-MM-DDTHH:MM:SS``).
    """
    user = getpass.getuser()
    result = _run([SQUEUE, "--noheader", "--user=%s" % user, "--name=%s" % name, "--format=%i"])
    if result.returncode != 0:
        raise SlurmError("squeue failed: %s" % result.stderr.strip())
    ids = [line.strip() for line in result.stdout.splitlines() if line.strip().isdigit()]
    if not ids:
        argv = [SACCT, "--noheader", "--parsable2", "--user=%s" % user, "--name=%s" % name]
        argv += ["--format=JobIDRaw"] + (["--starttime=%s" % since] if since else [])
        result = _run(argv)
        if result.returncode != 0:
            raise SlurmError("sacct failed: %s" % result.stderr.strip())
        ids = [line.strip() for line in result.stdout.splitlines() if line.strip().isdigit()]
    return max(ids, key=int) if ids else None


def cancel_job(job_id: str) -> None:
    """Ask Slurm to cancel job *job_id*."""
    result = _run([SCANCEL, job_id])
    if result.returncode != 0:
        raise SlurmError("scancel failed for job %s: %s" % (job_id, result.stderr.strip()))
