"""The study driver's Slurm executor. sbatch, squeue, sacct and scancel are all faked."""

import json
import os
import threading
import time

import pytest

import metawrap2.config as config_mod
from metawrap2 import command, limits, slurm
from metawrap2.commands import run as run_cmd
from metawrap2.manifest import COMPLETED, FAILED, INTERRUPTED
from metawrap2.slurm import SlurmResult


@pytest.fixture(autouse=True)
def _reset_runner():
    command.runner.configure(resume=False, force=False, dry_run=False)
    yield
    command.runner.configure(resume=False, force=False, dry_run=False)
    command.set_heartbeat_seconds(None)
    limits.clear_overrides()


def _slurm_run(tmp_path, slurm_table=None, overrides=None, settings=None):
    settings = settings or config_mod.Settings(slurm=slurm_table or {})
    return run_cmd.SlurmRun(settings, overrides or {}, tag="abc123", cwd=str(tmp_path))


class FakeSlurm:
    """Stands in for metawrap2.slurm's submit, wait, cancel and lookup calls."""

    def __init__(self, monkeypatch, results=None, job_ids=None):
        self.submitted = []
        self.waited = []
        self.cancelled = []
        self.results = dict(results or {})
        self.job_ids = list(job_ids or ["1001", "1002", "1003", "1004"])
        self.finished = {}
        monkeypatch.setattr(slurm, "submit_job", self.submit)
        monkeypatch.setattr(slurm, "wait_for_job", self.wait)
        monkeypatch.setattr(slurm, "cancel_job", self.cancelled.append)
        monkeypatch.setattr(slurm, "job_finished", lambda job_id: self.finished.get(job_id, True))
        monkeypatch.setattr(slurm, "find_job", lambda name, since=None: None)

    def submit(self, sbatch):
        self.submitted.append(list(sbatch))
        return self.job_ids.pop(0)

    def wait(self, job_id, poll_seconds=30, stop=None, on_state=None):
        self.waited.append(job_id)
        if on_state is not None:
            on_state("RUNNING")
        result = self.results.get(job_id, SlurmResult(job_id, "COMPLETED", 0))
        if isinstance(result, BaseException):
            raise result
        return result


def _step(tmp_path, name="assembly", module="assembly", argv=None):
    out = tmp_path / "study" / "ASSEMBLY"
    return run_cmd.Step(name, module, argv or ["-o", str(out)], outputs=[str(out / "final.fa")])


def _state(tmp_path, resume=False):
    return run_cmd.StepState(run_cmd.Layout(str(tmp_path / "study")), resume=resume)


def _run(tmp_path, step, state, slurm_run, flags=(), dry_run=False):
    layout = run_cmd.Layout(str(tmp_path / "study"))
    return run_cmd.run_slurm_step(step, layout, None, list(flags), dry_run, state, slurm_run)


# --- options: config, command line and per-module resources --------------------------------


def test_module_options_override_the_slurm_table_and_the_command_line_overrides_both(tmp_path):
    path = tmp_path / "metawrap2.toml"
    path.write_text(
        "[slurm]\n"
        'partition = "shared"\n'
        'account = "my_lab"\n'
        "threads = 8\n"
        "memory_gb = 32\n"
        "\n"
        "[slurm.modules.assembly]\n"
        'partition = "bigmem"\n'
        "threads = 32\n"
        "memory_gb = 500\n"
    )
    settings = config_mod.load_settings(str(path))
    assembly = config_mod.slurm_options(settings, "assembly", {"time": "48:00:00", "qos": None})
    qc = config_mod.slurm_options(settings, "read_qc")
    assert (assembly["partition"], assembly["threads"], assembly["memory_gb"]) == (
        "bigmem",
        32,
        500,
    )
    assert assembly["time"] == "48:00:00" and "qos" not in assembly
    assert assembly["account"] == "my_lab"
    assert (qc["partition"], qc["threads"], qc["memory_gb"]) == ("shared", 8, 32)
    assert config_mod.slurm_options(settings, "assembly", {"partition": "priority"})[
        "partition"
    ] == ("priority")


def test_a_gpu_request_in_one_form_replaces_the_other():
    """[slurm] gres and a module's gpus would otherwise reach sbatch together, which it refuses."""
    settings = config_mod.Settings(slurm={"gres": "gpu:1", "modules": {"binning": {"gpus": 2}}})
    assert config_mod.slurm_options(settings, "binning") == {"gpus": 2}
    assert config_mod.slurm_options(settings, "binning", {"gres": "gpu:a100:1"}) == {
        "gres": "gpu:a100:1"
    }


def test_slurm_config_mistakes_are_reported():
    def problems(table):
        return config_mod.check_config({"slurm": table})[0]

    assert any("partiton" in p for p in problems({"partiton": "bigmem"}))
    assert any("memroy_gb" in p for p in problems({"modules": {"assembly": {"memroy_gb": 500}}}))
    assert any("asembly" in p for p in problems({"modules": {"asembly": {"threads": 4}}}))
    assert any("positive whole number" in p for p in problems({"threads": "8"}))
    assert any("positive whole number" in p for p in problems({"memory_gb": 0}))
    assert any("list of strings" in p for p in problems({"extra_args": "--exclude=n1"}))
    assert any("not both" in p or "use one" in p for p in problems({"gres": "gpu:1", "gpus": 1}))
    assert any("poll_seconds" in p for p in problems({"modules": {"binning": {"poll_seconds": 5}}}))
    assert problems({"partition": "a", "threads": 4, "modules": {"binning": {"gpus": 1}}}) == []


def test_a_gpu_module_keeps_its_allocation_separate_from_cpu_modules(tmp_path):
    slurm_run = _slurm_run(
        tmp_path,
        {
            "partition": "shared",
            "threads": 8,
            "memory_gb": 32,
            "modules": {
                "binning": {"partition": "gpu", "threads": 16, "memory_gb": 128, "gres": "gpu:1"}
            },
        },
    )
    b_threads, b_memory, gpu = slurm_run.resources(run_cmd.Step("binning", "binning", []))
    q_threads, q_memory, cpu = slurm_run.resources(run_cmd.Step("read_qc:A", "read_qc", []))
    assert (b_threads, b_memory, gpu.partition, gpu.gres) == (16, 128, "gpu", "gpu:1")
    assert (q_threads, q_memory, cpu.partition, cpu.gres) == (8, 32, "shared", None)


def test_a_deep_binning_step_asks_for_a_gpu(tmp_path):
    """--preset deep trains on cuda; without a GPU in the job its GPU check stops the step."""
    slurm_run = _slurm_run(tmp_path)
    deep = run_cmd.Step("binning", "binning", ["--preset", "deep"])
    default = run_cmd.Step("binning", "binning", [])
    auto = run_cmd.Step("binning", "binning", ["--semibin2", "--device", "auto"])
    assert slurm_run.resources(deep)[2].gres == "gpu:1"
    assert slurm_run.resources(default)[2].gres is None
    assert slurm_run.resources(auto)[2].gres is None
    configured = _slurm_run(tmp_path, {"modules": {"binning": {"gpus": "a100:2"}}})
    job = configured.resources(deep)[2]
    assert (job.gres, job.gpus) == (None, "a100:2")


def test_slurm_jobs_are_not_capped_by_this_machines_budget(tmp_path):
    """The test machine has 32 cores and 128 GB; a job on a compute node may ask for more."""
    slurm_run = _slurm_run(tmp_path, overrides={"threads": 64, "memory_gb": 500})
    threads, memory, _job = slurm_run.resources(_step(tmp_path))
    assert (threads, memory) == (64, 500)


def test_the_command_line_overrides_a_modules_resources(tmp_path):
    slurm_run = _slurm_run(
        tmp_path,
        {"modules": {"assembly": {"threads": 32, "memory_gb": 500}}},
        overrides={"threads": 16, "memory_gb": None},
    )
    assert slurm_run.resources(_step(tmp_path))[:2] == (16, 500)


def test_a_job_name_is_one_word(tmp_path):
    step = run_cmd.Step("read_qc:sample 1;x", "read_qc", [])
    assert _slurm_run(tmp_path).job_name(step) == "mw2-abc123-read_qc_sample_1_x"


def test_unset_resources_use_the_defaults(tmp_path):
    assert _slurm_run(tmp_path).resources(_step(tmp_path))[:2] == (8, 32)
    settings = config_mod.Settings(threads=12)
    assert _slurm_run(tmp_path, settings=settings).resources(_step(tmp_path))[:2] == (12, 32)


# --- one step -------------------------------------------------------------------------------


def test_the_job_id_is_recorded_before_the_wait(tmp_path, monkeypatch):
    """A driver that stops while waiting leaves the ID in the manifest for --resume."""
    state = _state(tmp_path)
    seen = {}
    fake = FakeSlurm(monkeypatch)

    def wait(job_id, *args, **kwargs):
        seen["record"] = dict(state.manifest.step("assembly").data)
        return SlurmResult(job_id, "COMPLETED", 0, {"wall_seconds": 30.0})

    monkeypatch.setattr(slurm, "wait_for_job", wait)
    _, status, _log = _run(tmp_path, _step(tmp_path), state, _slurm_run(tmp_path))
    assert status == COMPLETED
    assert seen["record"]["status"] == "running"
    assert seen["record"]["scheduler_job_id"] == "1001"
    assert seen["record"]["executor"] == "slurm"
    assert seen["record"]["grant"]["job_name"] == "mw2-abc123-assembly"
    record = state.manifest.step("assembly").data
    assert record["status"] == COMPLETED and record["usage"] == {"wall_seconds": 30.0}
    assert len(fake.submitted) == 1


def test_the_job_script_and_log_are_in_the_study(tmp_path, monkeypatch):
    fake = FakeSlurm(monkeypatch)
    step = _step(tmp_path, argv=["-o", str(tmp_path / "study dir" / "ASSEMBLY")])
    _, _status, log = _run(tmp_path, step, _state(tmp_path), _slurm_run(tmp_path))
    sbatch = fake.submitted[0]
    script = sbatch[-1]
    assert script == os.path.join(
        str(tmp_path / "study"), ".metawrap2", "run_steps", "slurm", "assembly.sh"
    )
    with open(script) as fh:
        text = fh.read()
    assert "'%s'" % (tmp_path / "study dir" / "ASSEMBLY") in text
    assert log == str(tmp_path / "study" / "LOGS" / "assembly.log")
    assert "--output=%s" % log in sbatch
    with open(log) as fh:
        lines = fh.read().splitlines()
    assert lines[0].startswith("$ ") and "assembly" in lines[0]
    assert lines[1].startswith("# submitted with: sbatch ")
    assert lines[-1].startswith("[metawrap2] Slurm job 1001 COMPLETED (exit code 0)")


def test_a_failed_job_is_not_marked_completed(tmp_path, monkeypatch):
    FakeSlurm(monkeypatch, results={"1001": SlurmResult("1001", "OUT_OF_MEMORY", 1)})
    state = _state(tmp_path)
    _, status, log = _run(tmp_path, _step(tmp_path), state, _slurm_run(tmp_path))
    assert status == FAILED and state.manifest.step("assembly").status == FAILED
    with open(log) as fh:
        assert "OUT_OF_MEMORY" in fh.read()


def test_a_refused_submission_fails_the_step(tmp_path, monkeypatch):
    """A refusal is certain, so the step is failed, not left running for ever."""
    FakeSlurm(monkeypatch)

    def refuse(sbatch):
        raise slurm.SubmissionRejected("sbatch refused the job: Invalid account")

    monkeypatch.setattr(slurm, "submit_job", refuse)
    state = _state(tmp_path)
    _, status, log = _run(tmp_path, _step(tmp_path), state, _slurm_run(tmp_path))
    assert status == FAILED and state.manifest.step("assembly").status == FAILED
    with open(log) as fh:
        assert "Invalid account" in fh.read()


def test_an_uncertain_submission_is_found_by_name_on_resume(tmp_path, monkeypatch):
    fake = FakeSlurm(monkeypatch)

    def timeout(sbatch):
        raise slurm.SlurmError("sbatch did not answer")

    monkeypatch.setattr(slurm, "submit_job", timeout)
    state = _state(tmp_path)
    with pytest.raises(slurm.SlurmError):
        _run(tmp_path, _step(tmp_path), state, _slurm_run(tmp_path))
    assert state.manifest.step("assembly").status == "running"

    looked_up = []
    monkeypatch.setattr(slurm, "submit_job", fake.submit)
    monkeypatch.setattr(
        slurm, "find_job", lambda name, since=None: looked_up.append((name, since)) or "777"
    )
    resumed = _state(tmp_path, resume=True)
    _, status, _log = _run(tmp_path, _step(tmp_path), resumed, _slurm_run(tmp_path))
    assert status == COMPLETED and fake.submitted == [] and fake.waited == ["777"]
    assert looked_up[0][0] == "mw2-abc123-assembly"
    assert resumed.manifest.step("assembly").data["scheduler_job_id"] == "777"


def _left_running(tmp_path, step, job_id="555", flags=()):
    """A manifest in which an earlier run submitted *step* and then stopped."""
    state = _state(tmp_path)
    slurm_run = _slurm_run(tmp_path)
    threads, memory_gb, _job = slurm_run.resources(step)
    argv = step.command(threads, None, list(flags), memory_gb=memory_gb)
    state.start(step, argv, grant={"executor": "slurm", "job_name": slurm_run.job_name(step)})
    state.manifest.record_scheduler_job(step.name, "slurm", job_id)
    record = state.manifest.data["steps"][step.name]
    record["started"] = "2026-01-02T03:04:05"
    state.manifest.save()


def test_resume_reconnects_to_a_running_job(tmp_path, monkeypatch):
    step = _step(tmp_path)
    _left_running(tmp_path, step)
    fake = FakeSlurm(monkeypatch)
    state = _state(tmp_path, resume=True)
    _, status, _log = _run(tmp_path, step, state, _slurm_run(tmp_path), flags=["--resume"])
    assert status == COMPLETED
    assert fake.submitted == [] and fake.waited == ["555"]
    # The duration runs from the original submission.
    started = time.mktime(time.strptime("2026-01-02T03:04:05", "%Y-%m-%dT%H:%M:%S"))
    duration = state.manifest.step("assembly").data["duration_seconds"]
    assert duration >= time.time() - started - 5


def test_resume_does_not_resubmit_a_job_that_finished_while_the_driver_was_gone(
    tmp_path, monkeypatch
):
    step = _step(tmp_path)
    _left_running(tmp_path, step)
    fake = FakeSlurm(monkeypatch)
    state = _state(tmp_path, resume=True)
    _run(tmp_path, step, state, _slurm_run(tmp_path))
    assert fake.submitted == []
    assert state.manifest.step("assembly").status == COMPLETED


@pytest.mark.parametrize("ended", ["TIMEOUT", "CANCELLED", "NODE_FAIL", "PREEMPTED", "FAILED"])
def test_resume_resubmits_a_job_that_failed_while_the_driver_was_gone(tmp_path, monkeypatch, ended):
    step = _step(tmp_path)
    _left_running(tmp_path, step)
    fake = FakeSlurm(monkeypatch, results={"555": SlurmResult("555", ended, 1)})
    state = _state(tmp_path, resume=True)
    _, status, _log = _run(tmp_path, step, state, _slurm_run(tmp_path))
    assert status == COMPLETED
    assert fake.waited == ["555", "1001"] and len(fake.submitted) == 1
    assert state.manifest.step("assembly").data["scheduler_job_id"] == "1001"


def test_resume_refuses_a_second_job_while_one_with_another_command_runs(tmp_path, monkeypatch):
    _left_running(tmp_path, _step(tmp_path, argv=["-o", "elsewhere"]))
    fake = FakeSlurm(monkeypatch)
    fake.finished["555"] = False
    with pytest.raises(slurm.SlurmError, match="scancel 555"):
        _run(tmp_path, _step(tmp_path), _state(tmp_path, resume=True), _slurm_run(tmp_path))
    assert fake.submitted == []
    fake.finished["555"] = True
    _run(tmp_path, _step(tmp_path), _state(tmp_path, resume=True), _slurm_run(tmp_path))
    assert len(fake.submitted) == 1


def test_without_resume_a_left_job_is_not_reused(tmp_path, monkeypatch):
    step = _step(tmp_path)
    _left_running(tmp_path, step)
    fake = FakeSlurm(monkeypatch)
    _run(tmp_path, step, _state(tmp_path), _slurm_run(tmp_path))
    assert fake.waited == ["1001"]


def test_ctrl_c_cancels_the_job_and_records_the_step_as_interrupted(tmp_path, monkeypatch):
    fake = FakeSlurm(monkeypatch)
    slurm_run = _slurm_run(tmp_path)

    def wait(job_id, poll_seconds=30, stop=None, on_state=None):
        assert stop is slurm_run.stop
        raise slurm.JobInterrupted(job_id)

    monkeypatch.setattr(slurm, "wait_for_job", wait)
    state = _state(tmp_path)
    _, status, _log = _run(tmp_path, _step(tmp_path), state, slurm_run)
    assert status == INTERRUPTED and fake.cancelled == ["1001"]
    assert state.manifest.step("assembly").status == INTERRUPTED


def test_a_step_not_yet_submitted_at_ctrl_c_is_not_submitted(tmp_path, monkeypatch):
    fake = FakeSlurm(monkeypatch)
    slurm_run = _slurm_run(tmp_path)
    slurm_run.stop.set()
    _, status, _log = _run(tmp_path, _step(tmp_path), _state(tmp_path), slurm_run)
    assert status == INTERRUPTED and fake.submitted == []


def test_a_dry_run_prints_the_sbatch_command_and_writes_nothing(tmp_path, monkeypatch, capsys):
    fake = FakeSlurm(monkeypatch)
    slurm_run = _slurm_run(tmp_path, {"partition": "bigmem", "time": "48:00:00"})
    step = _step(tmp_path)
    layout = run_cmd.Layout(str(tmp_path / "study"))
    _, status, _log = run_cmd.run_slurm_step(
        step, layout, None, ["--dry-run"], True, None, slurm_run
    )
    out = capsys.readouterr().out
    assert status == COMPLETED and fake.submitted == []
    assert not (tmp_path / "study").exists()
    sbatch_line, job_line = out.splitlines()[:2]
    assert sbatch_line.startswith("[plan] sbatch --parsable ")
    assert "--partition=bigmem" in sbatch_line and "--time=48:00:00" in sbatch_line
    assert sbatch_line.endswith("assembly.sh")
    assert "job runs:" in job_line and "--dry-run assembly" in job_line


# --- a study's stages -----------------------------------------------------------------------


def _quiet(monkeypatch):
    for name in ("announcement", "comm", "warning"):
        monkeypatch.setattr(run_cmd, name, lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "audit_outputs", lambda *a: None)


def _execute(tmp_path, steps, slurm_run, jobs=2, dry_run=False, state=None):
    return run_cmd.execute(
        steps,
        run_cmd.Layout(str(tmp_path / "study")),
        run_cmd.Study([]),
        0,
        jobs,
        None,
        [],
        state or _state(tmp_path),
        dry_run,
        slurm=slurm_run,
    )


def test_every_job_of_a_stage_is_checked_before_the_next_stage(tmp_path, monkeypatch):
    """A failed first sample must not be hidden by a later one that succeeds."""
    calls = []
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(run_cmd, "_after_group", lambda *a: calls.append("advanced"))
    monkeypatch.setattr(
        run_cmd, "error", lambda message: (_ for _ in ()).throw(RuntimeError(message))
    )

    def fake_step(step, *args):
        calls.append(step.name)
        return step, FAILED if step.name == "read_qc:A" else COMPLETED, "log"

    monkeypatch.setattr(run_cmd, "run_slurm_step", fake_step)
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", [], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", [], "read_qc"),
    ]
    with pytest.raises(RuntimeError, match="1 of 2 step"):
        _execute(tmp_path, steps, _slurm_run(tmp_path))
    assert "read_qc:A" in calls and "read_qc:B" in calls
    assert "advanced" not in calls


def test_stages_keep_their_order_and_each_job_gets_its_modules_resources(tmp_path, monkeypatch):
    calls = []
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(
        run_cmd, "_after_group", lambda step, *a: calls.append(("after", step.module))
    )

    def fake_step(step, layout, config, flags, dry_run, state, slurm_run):
        threads, memory, job = slurm_run.resources(step)
        calls.append(("run", step.name, threads, memory, job.partition))
        return step, COMPLETED, "log"

    monkeypatch.setattr(run_cmd, "run_slurm_step", fake_step)
    slurm_run = _slurm_run(
        tmp_path,
        {
            "threads": 8,
            "memory_gb": 32,
            "partition": "shared",
            "modules": {"assembly": {"threads": 32, "memory_gb": 500, "partition": "bigmem"}},
        },
    )
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", [], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", [], "read_qc"),
        run_cmd.Step("assembly", "assembly", []),
    ]
    assert _execute(tmp_path, steps, slurm_run, dry_run=True) == 0
    qc = [c for c in calls if c[0] == "run" and c[1].startswith("read_qc")]
    assert len(qc) == 2 and all(c[2:] == (8, 32, "shared") for c in qc)
    assembly = ("run", "assembly", 32, 500, "bigmem")
    assert calls.index(("after", "read_qc")) < calls.index(assembly)


def test_an_unknown_job_outcome_stops_the_study_after_the_other_jobs(tmp_path, monkeypatch):
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    after = []
    monkeypatch.setattr(run_cmd, "_after_group", lambda *a: after.append(1))
    fake = FakeSlurm(
        monkeypatch,
        results={"1001": slurm.SlurmError("sacct has no final state")},
    )
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", ["-o", "A"], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", ["-o", "B"], "read_qc"),
    ]
    state = _state(tmp_path)
    assert _execute(tmp_path, steps, _slurm_run(tmp_path), jobs=1, state=state) == 1
    assert sorted(fake.waited) == ["1001", "1002"] and after == []
    # The unknown one stays running, so --resume reconnects to it.
    assert state.manifest.step("read_qc:A").status == "running"
    assert state.manifest.step("read_qc:B").status == COMPLETED


def test_ctrl_c_during_a_stage_stops_every_waiting_step(tmp_path, monkeypatch):
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    seen_stop = []

    def fake_step(step, layout, config, flags, dry_run, state, slurm_run):
        seen_stop.append(slurm_run.stop.wait(5))
        return step, INTERRUPTED, "log"

    def interrupted(futures):
        raise KeyboardInterrupt

    monkeypatch.setattr(run_cmd, "run_slurm_step", fake_step)
    monkeypatch.setattr(run_cmd.concurrent.futures, "as_completed", interrupted)
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", [], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", [], "read_qc"),
    ]
    assert _execute(tmp_path, steps, _slurm_run(tmp_path)) == 130
    assert seen_stop == [True, True]


def test_a_local_resume_warns_about_a_job_left_in_slurm(tmp_path, monkeypatch):
    step = _step(tmp_path)
    _left_running(tmp_path, step)
    warned = []
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "warning", warned.append)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(run_cmd, "run_step", lambda s, *a: (s, COMPLETED, "log"))
    run_cmd.execute(
        [step],
        run_cmd.Layout(str(tmp_path / "study")),
        run_cmd.Study([]),
        4,
        1,
        None,
        [],
        _state(tmp_path, resume=True),
        False,
    )
    assert any("job 555" in w for w in warned)


# --- the command line -----------------------------------------------------------------------

SHEET = """
[settings]
steps = ["read_qc", "assembly", "binning"]

[modules]
binning = "--preset deep"

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq.gz"
r2 = "raw/B_2.fastq.gz"
"""


def _main(tmp_path, monkeypatch, *extra, config=None):
    sheet = tmp_path / "samples.toml"
    sheet.write_text(SHEET)
    toml = tmp_path / "metawrap2.toml"
    toml.write_text(config or '[slurm]\npartition = "shared"\n')
    monkeypatch.setattr(run_cmd, "preflight_databases", lambda *args: None)
    monkeypatch.chdir(tmp_path)
    argv = [str(sheet), "-o", "study", "--config", "metawrap2.toml", "--executor", "slurm"]
    return run_cmd.main(argv + list(extra))


def test_the_slurm_executor_reaches_execute(tmp_path, monkeypatch):
    seen = {}

    def fake_execute(steps, layout, study, threads, jobs, config, flags, state, dry, **kwargs):
        seen.update(kwargs, jobs=jobs, config=config, flags=flags)
        return 0

    monkeypatch.setattr(run_cmd, "execute", fake_execute)
    monkeypatch.setattr(slurm, "missing_commands", list)
    assert _main(tmp_path, monkeypatch, "-t", "64", "-m", "500", "-j", "4") == 0
    slurm_run = seen["slurm"]
    assert isinstance(slurm_run, run_cmd.SlurmRun)
    assert slurm_run.overrides["threads"] == 64 and slurm_run.overrides["memory_gb"] == 500
    assert seen["jobs"] == 4
    # Each job is given the config by absolute path, for nodes without this home directory.
    assert seen["config"] == str(tmp_path / "metawrap2.toml")
    assert slurm_run.cwd == str(tmp_path)


def test_a_slurm_dry_run_plans_every_job_and_creates_nothing(tmp_path, monkeypatch, capsys):
    """Global flags reach each job's command, and nothing is written or submitted."""
    monkeypatch.setattr(slurm, "missing_commands", lambda: ["sbatch"])  # not needed to plan
    monkeypatch.setattr(slurm, "submit_job", lambda sbatch: pytest.fail("submitted"))
    from metawrap2 import cli

    cli._apply_global_flags(["--heartbeat-seconds", "120", "--max-disk-fraction", "0.9", "run"])
    assert _main(tmp_path, monkeypatch, "--dry-run", "--slurm-account", "my_lab", "-j", "2") == 0
    out = capsys.readouterr().out
    assert not (tmp_path / "study").exists()
    plans = [line for line in out.splitlines() if line.startswith("[plan] sbatch")]
    jobs = [line for line in out.splitlines() if "job runs:" in line]
    assert len(plans) == len(jobs) == 4  # read_qc A and B, assembly, binning
    assert all("--account=my_lab" in p and "--partition=shared" in p for p in plans)
    assert all("--heartbeat-seconds 120" in j and "--max-disk-fraction 0.9" in j for j in jobs)
    binning = next(p for p in plans if "binning" in p)
    assert "--gres=gpu:1" in binning
    assert not any("--gres" in p for p in plans if "binning" not in p)


def test_the_slurm_executor_needs_the_slurm_commands(tmp_path, monkeypatch):
    monkeypatch.setattr(slurm, "missing_commands", lambda: ["sbatch", "sacct"])
    with pytest.raises(SystemExit):
        _main(tmp_path, monkeypatch)


def test_a_slurm_study_records_the_cost_from_sacct(tmp_path, monkeypatch):
    """The study's provenance.txt totals the CPU time sacct reports for each job."""
    monkeypatch.setattr(slurm, "missing_commands", list)
    usage = {"wall_seconds": 3600.0, "user_seconds": 7000.0, "system_seconds": 200.0}
    FakeSlurm(
        monkeypatch,
        results={j: SlurmResult(j, "COMPLETED", 0, dict(usage)) for j in ("1001", "1002")},
    )
    monkeypatch.setattr(run_cmd, "_after_group", lambda *a: None)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    sheet = SHEET.replace('steps = ["read_qc", "assembly", "binning"]', 'steps = ["read_qc"]')
    (tmp_path / "samples.toml").write_text(sheet)
    monkeypatch.setattr(run_cmd, "preflight_databases", lambda *args: None)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "metawrap2.toml").write_text("")
    rc = run_cmd.main(["samples.toml", "-o", "study", "--executor", "slurm", "-j", "2"])
    assert rc == 0
    provenance = (tmp_path / "study" / "provenance.txt").read_text()
    assert "CPU-hours:    4.00" in provenance
    manifest = json.loads(
        (tmp_path / "study" / ".metawrap2" / "run_steps" / "manifest.json").read_text()
    )
    assert {s["scheduler_job_id"] for s in manifest["steps"].values()} == {"1001", "1002"}


def test_the_local_executor_is_still_the_default(tmp_path, monkeypatch):
    seen = {}

    def fake_execute(steps, layout, study, threads, jobs, config, flags, state, dry, **kwargs):
        seen.update(kwargs)
        return 0

    monkeypatch.setattr(run_cmd, "execute", fake_execute)
    monkeypatch.setattr(run_cmd, "preflight_databases", lambda *args: None)
    sheet = tmp_path / "samples.toml"
    sheet.write_text(SHEET)
    assert run_cmd.main([str(sheet), "-o", str(tmp_path / "study"), "--dry-run"]) == 0
    assert "slurm" not in seen and "gpus" in seen


def test_the_threads_running_jobs_are_not_left_behind(tmp_path, monkeypatch):
    """Every worker thread has ended when execute returns."""
    _quiet(monkeypatch)
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(run_cmd, "_after_group", lambda *a: None)
    FakeSlurm(monkeypatch)
    before = threading.active_count()
    steps = [run_cmd.Step("read_qc:%s" % s, "read_qc", ["-o", s], "read_qc") for s in "AB"]
    assert _execute(tmp_path, steps, _slurm_run(tmp_path)) == 0
    assert threading.active_count() == before
