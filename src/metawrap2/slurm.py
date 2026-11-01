"""Slurm job submission and monitoring for MetaWrap2.

`metawrap2 run` and `metawrap2 cohort` execute each planned module as a separate Slurm batch
job. Independent steps within a stage may run concurrently; the driver waits for every job in
the stage before advancing. The existing workflow planner, stage hooks, and resume logic remain
responsible for execution order.

A job runs the same `metawrap2 <module>` command as the local executor. Slurm allocates its
CPUs and memory, redirects output to the step's log, and reports its final state. The tools
within a module still execute normally inside that allocation, including their conda environments.

Submission and waiting are separate operations. `sbatch --parsable` returns the job ID immediately;
`wait_for_job` polls Slurm until the job finishes. This lets the driver record the job ID before
waiting, cancel an interrupted job, and eventually reconnect to existing jobs on `--resume`.

`squeue` reports pending and running jobs. `sacct` supplies final states and exit codes after jobs
leave the queue. Both commands must be available on the machine running the driver, and accounting
must be enabled for completed-job monitoring.

The batch script uses an absolute working directory and exports the submission environment so
the compute node can find MetaWrap2 and its conda environments. Input files, databases, scripts,
and output directories must be accessible from the compute node.

Limitation: this module does not manage the workflow manifest, verify module outputs, or reconnect
to jobs automatically. Those responsibilities belong to the study driver in commands/run.py.
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

#: Slurm commands used for submission, monitoring, and cancellation.
SBATCH = "sbatch"
SQUEUE = "squeue"
SACCT = "sacct"
SCANCEL = "scancel"

#: How often to check job status. Slurm accounting can lag behind the scheduler, so a job
#: disappearing from squeue is not treated as completed until sacct reports a terminal state.
POLL_SECONDS = 10
ACCOUNTING_TIMEOUT_SECONDS = 120

#: Slurm states that indicate the job is no longer running. COMPLETED is the only successful
#: state, and its exit code must also be zero. Everything else is reported as a failure.
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

#: Slurm sometimes appends a modifier such as "CANCELLED by 12345" or "FAILED+" to a state.
JOB_ID_PATTERN = re.compile(r"^([0-9]+)(?:;[A-Za-z0-9_.-]+)?$")

#: The batch script executes one MetaWrap2 module, not each external command within it.
#: Keeping the script minimal makes the generated command easy to inspect and reproduce.
BATCH_SCRIPT = """\
#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS={threads}

exec {command}
"""


@dataclass(frozen=True)
class SlurmConfig:
    """Resource and scheduling options for one Slurm job."""

    partition: Optional[str] = None
    account: Optional[str] = None
    time: str = "24:00:00"

    gres: Optional[str] = None
    gpus: Optional[str] = None
    constraint: Optional[str] = None
    qos: Optional[str] = None

    output: Optional[str] = None
    error: Optional[str] = None

    poll_seconds: int = POLL_SECONDS


@dataclass(frozen=True)
class SlurmResult:
    """Final outcome of one Slurm job."""

    job_id: str
    state: str
    returncode: int
    message: str = ""


class SlurmError(RuntimeError):
    """Raised when a Slurm command fails or the job's outcome cannot be determined."""


def _run(
    argv: Sequence[str],
    timeout: int = 60,
) -> subprocess.CompletedProcess[str]:
    """Run a Slurm control command, capturing its stdout and stderr."""
    try:
        return subprocess.run(
            list(argv),
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise SlurmError(f"Could not execute {argv[0]}: {exc}") from exc


def _require_command(name: str) -> None:
    """Check that a required Slurm command is on PATH."""
    if shutil.which(name) is None:
        raise SlurmError(
            f"{name} is not on PATH. The Slurm executor requires the Slurm client commands."
        )


def _job_id(output: str) -> str:
    """Extract the job ID from `sbatch --parsable` output."""
    value = output.strip()
    match = JOB_ID_PATTERN.fullmatch(value)
    if match is None:
        raise SlurmError(f"Cannot parse Slurm job ID from sbatch output: {value!r}")
    return match.group(1)


def _state(value: str) -> str:
    """Normalise a Slurm job state, dropping modifiers and trailing explanations."""
    return value.strip().split()[0].rstrip("+").upper() if value.strip() else ""


def _exit_code(value: str) -> int:
    """Extract the process exit code from Slurm's `exit:signal` representation."""
    try:
        exit_code, signal_number = value.strip().split(":", 1)
        code = int(exit_code)
        signal = int(signal_number)
        return code if signal == 0 else (128 + signal if code == 0 else code)
    except (ValueError, TypeError) as exc:
        raise SlurmError(f"Cannot parse Slurm exit code: {value!r}") from exc


def write_script(
    command: Sequence[str],
    path: str,
    threads: int,
) -> str:
    """Atomically write a batch script that executes one MetaWrap2 module."""
    if not command:
        raise ValueError("A Slurm job must have a command.")
    if threads < 1:
        raise ValueError("Slurm threads must be positive.")

    script = Path(path).resolve()
    script.parent.mkdir(parents=True, exist_ok=True)

    contents = BATCH_SCRIPT.format(
        threads=threads,
        command=shlex.join(list(command)),
    )

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{script.name}.",
            suffix=".tmp",
            dir=script.parent,
            delete=False,
        ) as fh:
            temporary = fh.name
            fh.write(contents)

        os.replace(temporary, script)
    finally:
        if temporary is not None and os.path.exists(temporary):
            os.unlink(temporary)

    return str(script)


def submit_job(
    command: Sequence[str],
    *,
    name: str,
    threads: int,
    memory_gb: int,
    script_path: str,
    log_path: str,
    config: SlurmConfig,
    working_dir: Optional[str] = None,
) -> str:
    """Submit one module to Slurm and return its job ID without waiting.

    The job's stdout and stderr are written to *log_path*. The batch script and working
    directory use absolute paths so execution does not depend on Slurm's starting directory.

    The caller must persist the returned ID before waiting or starting another job.
    """
    _require_command(SBATCH)

    if threads < 1 or memory_gb < 1:
        raise ValueError("Slurm threads and memory must be positive.")

    script = write_script(command, script_path, threads)
    log = str(Path(log_path).resolve())
    cwd = str(Path(working_dir or os.getcwd()).resolve())

    Path(log).parent.mkdir(parents=True, exist_ok=True)

    argv = [
        SBATCH,
        "--parsable",
        f"--job-name={name}",
        "--ntasks=1",
        f"--cpus-per-task={threads}",
        f"--mem={memory_gb}G",
        f"--time={config.time}",
        f"--chdir={cwd}",
        "--export=ALL",
    ]

    if config.partition:
        argv.append(f"--partition={config.partition}")
    if config.account:
        argv.append(f"--account={config.account}")

    argv.append(script)

    if config.gres and config.gpus:
        raise ValueError("Specify either gres or gpus, not both.")

    if config.gres:
        argv.append(f"--gres={config.gres}")

    if config.gpus:
        argv.append(f"--gpus={config.gpus}")

    if config.constraint:
        argv.append(f"--constraint={config.constraint}")

    if config.qos:
        argv.append(f"--qos={config.qos}")

    stdout_path = config.output or log_path
    stderr_path = config.error or stdout_path

    stdout_path = os.path.abspath(stdout_path)
    stderr_path = os.path.abspath(stderr_path)

    os.makedirs(os.path.dirname(stdout_path), exist_ok=True)
    os.makedirs(os.path.dirname(stderr_path), exist_ok=True)

    argv += [
        f"--output={stdout_path}",
        f"--error={stderr_path}",
    ]

    result = _run(argv)

    if result.returncode != 0:
        raise SlurmError(
            f"Slurm rejected job {name}: {result.stderr.strip() or result.stdout.strip()}"
        )

    return _job_id(result.stdout)


def _queued(job_id: str) -> bool:
    """Return whether the job is still present in Slurm's active queue.

    A job may disappear from squeue after completing, failing, or being
    cancelled. In that case, return False so wait_for_job() can query sacct
    for the final state.
    """
    result = _run([SQUEUE, "--noheader", "--jobs", job_id, "--format=%i"])

    if result.returncode != 0:
        message = (result.stderr or "").strip()

        # Slurm removes terminal jobs from its active queue.
        # Missing jobs must be checked through sacct, not treated as errors.
        if "invalid job id specified" in message.lower():
            return False

        raise SlurmError(f"Could not query squeue for job {job_id}: {message}")

    return any(line.strip() == str(job_id) for line in result.stdout.splitlines())


def _accounting_result(job_id: str) -> Optional[SlurmResult]:
    """Read the top-level job record from sacct, excluding batch and extern substeps."""
    result = _run(
        [
            SACCT,
            "--noheader",
            "--parsable2",
            "--jobs",
            job_id,
            "--format=JobIDRaw,State,ExitCode",
        ]
    )

    if result.returncode != 0:
        raise SlurmError(f"Could not query sacct for job {job_id}: {result.stderr.strip()}")

    for line in result.stdout.splitlines():
        fields = line.split("|")
        if len(fields) < 3 or fields[0].strip() != job_id:
            continue

        state = _state(fields[1])
        if state not in TERMINAL_STATES:
            return None

        code = _exit_code(fields[2])
        if state != "COMPLETED" and code == 0:
            code = 1

        return SlurmResult(
            job_id=job_id,
            state=state,
            returncode=code,
            message=f"Slurm job {job_id}: {state} (exit code {fields[2].strip()})",
        )

    return None


def wait_for_job(
    job_id: str,
    config: SlurmConfig,
) -> SlurmResult:
    """Block until a submitted or previously recorded Slurm job terminates.

    Poll squeue while the job is active. Once it disappears, query sacct
    for the terminal state.

    If accounting remains unavailable after the grace period, raise
    SlurmError without cancelling the job or assuming that it failed.
    """
    _require_command(SQUEUE)
    _require_command(SACCT)

    if config.poll_seconds < 1:
        raise ValueError("Slurm polling interval must be positive.")

    accounting_since = None

    while True:
        if _queued(job_id):
            accounting_since = None

        else:
            result = _accounting_result(job_id)

            if result is not None:
                # Ensure the accounting result is a terminal state.
                # Do not return RUNNING or PENDING as final results.
                state = result.state.split("+", 1)[0].strip().upper()

                if state not in {
                    "PENDING",
                    "CONFIGURING",
                    "RUNNING",
                    "COMPLETING",
                    "SUSPENDED",
                    "RESIZING",
                    "REQUEUED",
                }:
                    return result

            if accounting_since is None:
                accounting_since = time.monotonic()

            elif time.monotonic() - accounting_since >= ACCOUNTING_TIMEOUT_SECONDS:
                raise SlurmError(
                    f"Job {job_id} is no longer in squeue, but sacct "
                    "has not provided a terminal state within "
                    f"{ACCOUNTING_TIMEOUT_SECONDS} seconds. "
                    "Check the scheduler before resubmitting."
                )

        time.sleep(config.poll_seconds)


def cancel_job(job_id: str) -> None:
    """Request cancellation of a Slurm job.

    The caller should confirm its terminal state before restarting the step. A successful
    scancel request alone does not prove that the job has finished terminating.
    """
    _require_command(SCANCEL)
    result = _run([SCANCEL, job_id])

    if result.returncode != 0:
        raise SlurmError(f"Could not cancel Slurm job {job_id}: {result.stderr.strip()}")


def submit_and_wait(
    command: Sequence[str],
    *,
    name: str,
    threads: int,
    memory_gb: int,
    script_path: str,
    log_path: str,
    config: SlurmConfig,
    working_dir: Optional[str] = None,
) -> SlurmResult:
    """Submit one MetaWrap2 module and block until Slurm reports its final state.

    This convenience function suits simple blocking execution. For resumable studies, prefer
    calling submit_job and wait_for_job separately so the driver can persist the job ID.
    """
    job_id = submit_job(
        command,
        name=name,
        threads=threads,
        memory_gb=memory_gb,
        script_path=script_path,
        log_path=log_path,
        config=config,
        working_dir=working_dir,
    )

    return wait_for_job(job_id, config)
