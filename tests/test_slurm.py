"""Slurm submission and job-state tests. sbatch, squeue, sacct and scancel are all faked."""

import json
import os
import subprocess
import sys
import threading

import pytest

from metawrap2 import slurm


def _done(stdout="", returncode=0, stderr=""):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


class FakeScheduler:
    """Answers squeue and sacct from a script of polls.

    Each poll is ``(squeue state or None, sacct lines or None)``; the last one repeats.
    ``"error"`` in place of the squeue state makes squeue fail.
    """

    def __init__(self, job_id, polls):
        self.job_id = job_id
        self.polls = list(polls)
        self.calls = []

    def _current(self):
        return self.polls[0] if len(self.polls) == 1 else self.polls.pop(0)

    def run(self, argv, timeout=60):
        self.calls.append(list(argv))
        if argv[0] == slurm.SQUEUE:
            self._poll = self._current()
            state = self._poll[0]
            if state == "error":
                return _done(returncode=1, stderr="slurm_load_jobs error: Socket timed out")
            return _done("%s|%s\n" % (self.job_id, state) if state else "")
        if argv[0] == slurm.SACCT:
            return _done("\n".join(self._poll[1] or []) + "\n")
        raise AssertionError(argv)


@pytest.fixture
def no_sleep(monkeypatch):
    monkeypatch.setattr(slurm.time, "sleep", lambda seconds: None)


# --- the batch script -----------------------------------------------------------------------


def test_the_batch_script_keeps_every_argument_one_word(tmp_path):
    """Spaces, quotes and shell syntax in paths reach the module unchanged."""
    cwd = tmp_path / "study dir"
    cwd.mkdir()
    awkward = ["a b", "$(touch wrong)", "it's", "x;y", "*"]
    command = [sys.executable, "-c", "import json, sys; print(json.dumps(sys.argv[1:]))"]
    script = cwd / "job's script.sh"
    slurm.write_script(str(script), slurm.batch_script(command + awkward, "read_qc:A", str(cwd), 4))
    out = subprocess.run(["bash", str(script)], capture_output=True, text=True, check=True).stdout
    assert json.loads(out) == awkward
    assert not (cwd / "wrong").exists()
    assert os.access(str(script), os.X_OK)
    text = script.read_text()
    assert text.startswith("#!/usr/bin/env bash\n")
    assert "export OMP_NUM_THREADS=4" in text


def test_the_batch_script_runs_in_the_drivers_directory(tmp_path):
    cwd = tmp_path / "with space"
    cwd.mkdir()
    script = tmp_path / "job.sh"
    slurm.write_script(str(script), slurm.batch_script(["pwd"], "s", str(cwd), 1))
    out = subprocess.run(["bash", str(script)], capture_output=True, text=True, check=True).stdout
    assert os.path.realpath(out.strip()) == os.path.realpath(str(cwd))


def test_a_batch_script_needs_a_command_and_a_thread(tmp_path):
    with pytest.raises(ValueError):
        slurm.batch_script([], "s", str(tmp_path), 2)
    with pytest.raises(ValueError):
        slurm.batch_script(["true"], "s", str(tmp_path), 0)


# --- the sbatch command ---------------------------------------------------------------------


def _sbatch(job, log="/study/LOGS/binning.log", **kwargs):
    options = {
        "name": "mw2-abc123-binning",
        "threads": 32,
        "memory_gb": 500,
        "log": log,
        "job": job,
    }
    options.update(kwargs)
    return slurm.sbatch_command("/study/job.sh", **options)


def test_every_sbatch_option_comes_before_the_script():
    """sbatch passes anything after the script to the script, so a --gres there is lost."""
    job = slurm.SlurmJob(
        partition="bigmem",
        account="my_lab",
        time="48:00:00",
        gres="gpu:1",
        qos="long",
        constraint="avx2",
        extra_args=("--exclude=node7",),
    )
    argv = _sbatch(job)
    assert argv[0] == "sbatch" and argv[-1] == "/study/job.sh"
    for flag in (
        "--parsable",
        "--job-name=mw2-abc123-binning",
        "--ntasks=1",
        "--cpus-per-task=32",
        "--mem=500G",
        "--export=ALL",
        "--partition=bigmem",
        "--account=my_lab",
        "--time=48:00:00",
        "--gres=gpu:1",
        "--qos=long",
        "--constraint=avx2",
        "--exclude=node7",
        "--output=/study/LOGS/binning.log",
        "--open-mode=append",
    ):
        assert flag in argv[:-1]
    assert not any(a.startswith("--error") for a in argv)


def test_unset_options_are_left_to_the_cluster():
    argv = _sbatch(slurm.SlurmJob())
    assert not any(
        a.split("=")[0] in ("--time", "--partition", "--account", "--gres", "--gpus") for a in argv
    )


def test_a_percent_sign_in_the_log_path_is_not_a_slurm_pattern():
    argv = _sbatch(slurm.SlurmJob(), log="/runs/100%_test/LOGS/a.log")
    assert "--output=/runs/100%%_test/LOGS/a.log" in argv


def test_gres_and_gpus_together_are_refused():
    with pytest.raises(ValueError):
        _sbatch(slurm.SlurmJob(gres="gpu:1", gpus="1"))


# --- submission -----------------------------------------------------------------------------


def test_submission_returns_the_job_id(monkeypatch):
    monkeypatch.setattr(slurm.subprocess, "run", lambda *a, **k: _done("123456;cluster\n"))
    assert slurm.submit_job(["sbatch", "job.sh"]) == "123456"


def test_a_refused_submission_names_the_reason(monkeypatch):
    monkeypatch.setattr(
        slurm.subprocess,
        "run",
        lambda *a, **k: _done(returncode=1, stderr="sbatch: error: Invalid account"),
    )
    with pytest.raises(slurm.SubmissionRejected, match="Invalid account"):
        slurm.submit_job(["sbatch", "job.sh"])


def test_a_timed_out_submission_is_not_reported_as_refused(monkeypatch):
    """sbatch may have submitted the job, so the step must not be treated as never started."""

    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired("sbatch", 120)

    monkeypatch.setattr(slurm.subprocess, "run", timeout)
    with pytest.raises(slurm.SlurmError) as raised:
        slurm.submit_job(["sbatch", "job.sh"])
    assert not isinstance(raised.value, slurm.SubmissionRejected)


# --- final state from sacct -----------------------------------------------------------------


def _sacct(monkeypatch, *lines):
    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: _done("\n".join(lines) + "\n"))
    return slurm.accounting_result("123")


def test_completed_needs_a_zero_exit_code(monkeypatch):
    result = _sacct(monkeypatch, "123|COMPLETED|2:0|60|||", "123.batch|COMPLETED|2:0|60|||")
    assert result.state == "COMPLETED" and result.returncode == 2 and not result.ok


def test_out_of_memory_with_exit_code_zero_still_fails(monkeypatch):
    result = _sacct(monkeypatch, "123|OUT_OF_MEMORY|0:0|60|||")
    assert result.state == "OUT_OF_MEMORY" and result.returncode != 0 and not result.ok


def test_a_failed_batch_step_fails_a_completed_job(monkeypatch):
    result = _sacct(
        monkeypatch,
        "123|COMPLETED|0:0|60|||",
        "123.batch|OUT_OF_MEMORY|0:125|60|||9G",
        "123.extern|COMPLETED|0:0|60|||",
    )
    assert result.state == "OUT_OF_MEMORY" and not result.ok


@pytest.mark.parametrize(
    "state,expected",
    [
        ("CANCELLED by 1000", "CANCELLED"),
        ("CANCELLED+", "CANCELLED"),
        ("TIMEOUT", "TIMEOUT"),
        ("NODE_FAIL", "NODE_FAIL"),
        ("PREEMPTED", "PREEMPTED"),
    ],
)
def test_ended_states_are_failures(monkeypatch, state, expected):
    result = _sacct(monkeypatch, "123|%s|0:0|60|||" % state)
    assert result.state == expected and result.returncode == 1 and not result.ok


def test_a_job_killed_by_a_signal_gets_128_plus_the_signal(monkeypatch):
    assert _sacct(monkeypatch, "123|FAILED|0:9|60|||").returncode == 137


@pytest.mark.parametrize("state", ["PENDING", "CONFIGURING", "RUNNING", "REQUEUED", "SUSPENDED"])
def test_a_job_that_has_not_ended_has_no_result(monkeypatch, state):
    assert _sacct(monkeypatch, "123|%s|0:0|0|||" % state) is None


def test_the_cost_of_a_job_comes_from_sacct(monkeypatch):
    result = _sacct(
        monkeypatch,
        "123|COMPLETED|0:0|3725|01:02:03|00:10.500|",
        "123.batch|COMPLETED|0:0|3725|01:02:03|00:10.500|2G",
        "123.extern|COMPLETED|0:0|3725|00:00:00|00:00:00|1024K",
    )
    assert result.ok
    assert result.usage == {
        "wall_seconds": 3725.0,
        "user_seconds": 3723.0,
        "system_seconds": 10.5,
        "peak_rss_bytes": 2 * 1024**3,
    }


def test_a_cost_over_a_day_is_parsed(monkeypatch):
    result = _sacct(monkeypatch, "123|COMPLETED|0:0|90000|1-01:00:00|00:00:01|")
    assert result.usage["user_seconds"] == 90000.0


# --- queue state ----------------------------------------------------------------------------


def test_a_job_gone_from_squeue_is_not_an_error(monkeypatch):
    """Older releases fail with 'Invalid job id' once a job has left the queue."""
    monkeypatch.setattr(
        slurm,
        "_run",
        lambda argv, timeout=60: _done(
            returncode=1, stderr="slurm_load_jobs error: Invalid job id specified"
        ),
    )
    assert slurm.queue_state("32766245") is None


# --- waiting --------------------------------------------------------------------------------


def test_waiting_follows_the_job_through_its_states(monkeypatch, no_sleep):
    fake = FakeScheduler(
        "123",
        [
            ("PENDING", None),
            ("CONFIGURING", None),
            ("RUNNING", None),
            ("RUNNING", None),
            (None, ["123|RUNNING|0:0|10|||"]),  # sacct lags the queue
            (None, ["123|COMPLETED|0:0|10|||"]),
        ],
    )
    monkeypatch.setattr(slurm, "_run", fake.run)
    seen = []
    result = slurm.wait_for_job("123", 5, on_state=seen.append)
    assert result.ok
    assert seen == ["PENDING", "CONFIGURING", "RUNNING"]


def test_an_ended_job_still_in_squeue_is_read_from_sacct(monkeypatch, no_sleep):
    fake = FakeScheduler("123", [("TIMEOUT", ["123|TIMEOUT|0:15|10|||"])])
    monkeypatch.setattr(slurm, "_run", fake.run)
    result = slurm.wait_for_job("123", 5)
    assert result.state == "TIMEOUT" and not result.ok
    assert len([c for c in fake.calls if c[0] == slurm.SQUEUE]) == 1


def test_a_preempted_job_that_is_requeued_is_followed_to_the_end(monkeypatch, no_sleep):
    fake = FakeScheduler(
        "123",
        [
            (None, ["123|PREEMPTED|0:0|10|||"]),
            ("PENDING", None),
            ("RUNNING", None),
            (None, ["123|COMPLETED|0:0|10|||"]),
        ],
    )
    monkeypatch.setattr(slurm, "_run", fake.run)
    assert slurm.wait_for_job("123", 5).ok


def test_a_preempted_job_that_stays_preempted_fails(monkeypatch, no_sleep):
    fake = FakeScheduler("123", [(None, ["123|PREEMPTED|0:0|10|||"])])
    monkeypatch.setattr(slurm, "_run", fake.run)
    result = slurm.wait_for_job("123", 5)
    assert result.state == "PREEMPTED" and not result.ok


def test_a_busy_controller_does_not_stop_the_wait(monkeypatch, no_sleep):
    fake = FakeScheduler(
        "123",
        [("error", None), ("error", None), ("RUNNING", None), (None, ["123|COMPLETED|0:0|1|||"])],
    )
    monkeypatch.setattr(slurm, "_run", fake.run)
    assert slurm.wait_for_job("123", 5).ok


def test_a_controller_that_keeps_failing_ends_the_wait(monkeypatch, no_sleep):
    monkeypatch.setattr(slurm, "QUERY_FAILURE_TIMEOUT_SECONDS", 0)
    fake = FakeScheduler("123", [("error", None)])
    monkeypatch.setattr(slurm, "_run", fake.run)
    with pytest.raises(slurm.SlurmError, match="could not query job 123"):
        slurm.wait_for_job("123", 5)


def test_a_job_with_no_accounting_record_ends_the_wait(monkeypatch, no_sleep):
    monkeypatch.setattr(slurm, "ACCOUNTING_TIMEOUT_SECONDS", 0)
    fake = FakeScheduler("123", [(None, [])])
    monkeypatch.setattr(slurm, "_run", fake.run)
    with pytest.raises(slurm.SlurmError, match="sacct -j 123"):
        slurm.wait_for_job("123", 5)


def test_a_stop_request_ends_the_wait_without_cancelling(monkeypatch):
    fake = FakeScheduler("123", [("RUNNING", None)])
    monkeypatch.setattr(slurm, "_run", fake.run)
    stop = threading.Event()
    stop.set()
    with pytest.raises(slurm.JobInterrupted):
        slurm.wait_for_job("123", 5, stop=stop)
    assert all(call[0] != slurm.SCANCEL for call in fake.calls)


def test_the_wait_sleeps_for_the_poll_interval(monkeypatch):
    slept = []
    monkeypatch.setattr(slurm.time, "sleep", slept.append)
    fake = FakeScheduler("123", [("RUNNING", None), (None, ["123|COMPLETED|0:0|1|||"])])
    monkeypatch.setattr(slurm, "_run", fake.run)
    slurm.wait_for_job("123", 45)
    assert slept == [45]
    with pytest.raises(ValueError):
        slurm.wait_for_job("123", 0)


# --- cancelling and finding jobs ------------------------------------------------------------


def test_cancelling_a_job_uses_its_exact_id(monkeypatch):
    calls = []
    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: calls.append(list(argv)) or _done())
    slurm.cancel_job("123456")
    assert calls == [["scancel", "123456"]]


def test_a_job_is_found_by_name_in_the_queue_then_in_accounting(monkeypatch):
    calls = []

    def fake_run(argv, timeout=60):
        calls.append(list(argv))
        if argv[0] == slurm.SQUEUE:
            return _done("")
        return _done("41\n41.batch\n57\n")

    monkeypatch.setattr(slurm, "_run", fake_run)
    assert slurm.find_job("mw2-abc123-assembly", since="2026-01-02T03:04:05") == "57"
    assert "--name=mw2-abc123-assembly" in calls[0]
    assert "--starttime=2026-01-02T03:04:05" in calls[1]

    monkeypatch.setattr(slurm, "_run", lambda argv, timeout=60: _done("99\n"))
    assert slurm.find_job("mw2-abc123-assembly") == "99"
