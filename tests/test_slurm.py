"""Slurm submission and monitoring tests; all scheduler commands are mocked."""

import subprocess

import pytest

from metawrap2 import slurm


# --- generated batch scripts ---------------------------------------------------------------


def test_a_slurm_script_preserves_quoted_command_arguments(tmp_path):
    """Module arguments containing spaces or shell metacharacters must remain one argument."""
    script = tmp_path / "job.sh"
    slurm.write_script(["python", "-m", "metawrap2.cli", "read_qc", "a b; $(touch wrong)"], str(script), 8)
    content = script.read_text()
    assert content.startswith("#!/usr/bin/env bash\n")
    assert "export OMP_NUM_THREADS=8" in content
    assert "exec python -m metawrap2.cli read_qc 'a b; $(touch wrong)'" in content


def test_a_slurm_script_rejects_an_empty_command_or_zero_threads(tmp_path):
    """Do not generate a batch job without a command or a valid CPU allocation."""
    with pytest.raises(ValueError):
        slurm.write_script([], str(tmp_path / "empty.sh"), 2)
    with pytest.raises(ValueError):
        slurm.write_script(["true"], str(tmp_path / "zero.sh"), 0)


# --- sbatch arguments and job identifiers -------------------------------------------------


def test_a_slurm_job_passes_partition_account_resources_and_gpu_flags(tmp_path, monkeypatch):
    """A module can request a bigmem or GPU allocation without changing the workflow driver."""
    calls = []
    monkeypatch.setattr(slurm.shutil, "which", lambda name: "/usr/bin/" + name)

    def fake_run(argv, timeout=60):
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 0, "123456;cluster\n", "")

    monkeypatch.setattr(slurm, "_run", fake_run)
    job_id = slurm.submit_job(
        ["python", "-m", "metawrap2.cli", "binning"],
        name="mw2_binning", threads=32, memory_gb=500,
        script_path=str(tmp_path / "job.sh"), log_path=str(tmp_path / "job.log"),
        config=slurm.SlurmConfig(partition="bigmem", account="mschatz1_bigmem",
                                 time="48:00:00", gres="gpu:1"),
    )
    assert job_id == "123456"
    argv = calls[0]
    for flag in ("--job-name=mw2_binning", "--cpus-per-task=32", "--mem=500G",
                 "--partition=bigmem", "--account=mschatz1_bigmem", "--time=48:00:00",
                 "--gres=gpu:1", "--ntasks=1"):
        assert flag in argv


def test_a_rejected_sbatch_submission_reports_the_scheduler_error(tmp_path, monkeypatch):
    """An invalid account or partition fails before a workflow step is considered running."""
    monkeypatch.setattr(slurm.shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: subprocess.CompletedProcess(
        argv, 1, "", "sbatch: error: invalid account"))
    with pytest.raises(slurm.SlurmError, match="invalid account"):
        slurm.submit_job(
            ["true"], name="mw2_test", threads=1, memory_gb=2,
            script_path=str(tmp_path / "job.sh"), log_path=str(tmp_path / "job.log"),
            config=slurm.SlurmConfig(),
        )


# --- job state and cancellation ------------------------------------------------------------


def test_a_completed_job_requires_a_zero_exit_code(monkeypatch):
    """Do not mistake a COMPLETED Slurm state for success when the process returned nonzero."""
    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: subprocess.CompletedProcess(
        argv, 0, "123|COMPLETED|2:0\n123.batch|COMPLETED|2:0\n", ""))
    result = slurm._accounting_result("123")
    assert result.state == "COMPLETED"
    assert result.returncode == 2


def test_an_out_of_memory_job_cannot_be_reported_successful(monkeypatch):
    """Slurm can report OOM with a zero process exit code, but the step must still fail."""
    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: subprocess.CompletedProcess(
        argv, 0, "123|OUT_OF_MEMORY|0:0\n", ""))
    result = slurm._accounting_result("123")
    assert result.state == "OUT_OF_MEMORY"
    assert result.returncode != 0


def test_a_waited_job_does_not_finish_until_accounting_confirms_its_state(monkeypatch):
    """A job disappearing from squeue is not enough to mark it completed."""
    checks = iter([None, slurm.SlurmResult("123", "COMPLETED", 0)])
    monkeypatch.setattr(slurm, "_require_command", lambda name: None)
    monkeypatch.setattr(slurm, "_queued", lambda job_id: False)
    monkeypatch.setattr(slurm, "_accounting_result", lambda job_id: next(checks))
    monkeypatch.setattr(slurm.time, "sleep", lambda seconds: None)
    result = slurm.wait_for_job("123", slurm.SlurmConfig(poll_seconds=1))
    assert result.returncode == 0


def test_cancelling_a_job_uses_its_exact_slurm_id(monkeypatch):
    """Cancellation should target the submitted job, never an entire user's jobs."""
    calls = []
    monkeypatch.setattr(slurm, "_require_command", lambda name: None)

    def fake_run(argv, timeout=60):
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(slurm, "_run", fake_run)
    slurm.cancel_job("123456")
    assert calls == [["scancel", "123456"]]


def test_queued_treats_expired_job_id_as_not_active(monkeypatch):
    """A job missing from squeue must fall through to sacct."""
    from subprocess import CompletedProcess

    from metawrap2 import slurm

    def fake_run(argv):
        return CompletedProcess(
            args=argv,
            returncode=1,
            stdout="",
            stderr="slurm_load_jobs error: Invalid job id specified",
        )

    monkeypatch.setattr(slurm, "_run", fake_run)

    assert slurm._queued("32766245") is False