"""Regression tests for the MetaWrap2 study driver's blocking Slurm executor.

These tests mock scheduler operations; they never submit jobs to a real Slurm cluster.
They cover configuration, per-step resources, stage barriers, and CLI dispatch.
"""

import pytest

from metawrap2 import config as config_mod
from metawrap2.commands import run as run_cmd
from metawrap2.manifest import COMPLETED, FAILED
from metawrap2.slurm import SlurmResult

# --- Slurm configuration and per-module resource resolution ---------------------------------


def test_a_slurm_config_uses_module_resources_and_cli_overrides(tmp_path):
    """Module settings override global defaults, while explicit CLI settings take precedence.

    Assembly should be able to request a big-memory node without moving read QC there too.
    """
    path = tmp_path / "metawrap2.toml"
    path.write_text(
        "[slurm]\n"
        'partition = "shared"\n'
        'account = "mschatz1"\n'
        "threads = 8\n"
        "memory_gb = 32\n"
        "\n"
        "[slurm.modules.assembly]\n"
        'partition = "bigmem"\n'
        "threads = 32\n"
        "memory_gb = 500\n"
    )

    settings = config_mod.load_settings(str(path))
    assembly = config_mod.resolve_slurm_settings(settings, "assembly", {"time": "48:00:00"})
    qc = config_mod.resolve_slurm_settings(settings, "read_qc")

    assert (assembly["partition"], assembly["threads"], assembly["memory_gb"]) == (
        "bigmem",
        32,
        500,
    )
    assert assembly["time"] == "48:00:00"
    assert (qc["partition"], qc["threads"], qc["memory_gb"]) == ("shared", 8, 32)
    assert (
        config_mod.resolve_slurm_settings(settings, "assembly", {"partition": "priority"})[
            "partition"
        ]
        == "priority"
    )


def test_an_unknown_slurm_setting_is_rejected_before_submission():
    """Misspelled Slurm keys must not silently fall back to the wrong partition or resources."""
    problems, _ = config_mod.check_config({"slurm": {"partiton": "bigmem"}})
    assert any("partiton" in problem for problem in problems)

    problems, _ = config_mod.check_config({"slurm": {"modules": {"assembly": {"memroy_gb": 500}}}})
    assert any("memroy_gb" in problem for problem in problems)


def test_a_gpu_module_keeps_its_allocation_separate_from_cpu_modules():
    """Selecting GPU binning must not send unrelated read-QC jobs to a GPU partition."""
    settings = config_mod.Settings(
        slurm={
            "partition": "shared",
            "threads": 8,
            "memory_gb": 32,
            "modules": {
                "binning": {
                    "partition": "gpu",
                    "threads": 16,
                    "memory_gb": 128,
                    "gres": "gpu:1",
                }
            },
        }
    )

    binning = run_cmd.Step("binning", "binning", [])
    read_qc = run_cmd.Step("read_qc:A", "read_qc", [])
    b_threads, b_memory, gpu = run_cmd.slurm_resources_for(binning, settings, {})
    q_threads, q_memory, cpu = run_cmd.slurm_resources_for(read_qc, settings, {})

    assert (b_threads, b_memory, gpu.partition, gpu.gres) == (16, 128, "gpu", "gpu:1")
    assert (q_threads, q_memory, cpu.partition, cpu.gres) == (8, 32, "shared", None)


def test_slurm_resources_use_requested_fallbacks_without_clamping_to_the_login_node():
    """Unconfigured Slurm jobs use the supplied default grant, not the login node's budget.

    This covers the resource-budget fix: a 32-CPU, 500-GB compute job must remain possible even
    when the controller is running on a much smaller login node.
    """
    settings = config_mod.Settings(slurm={})
    step = run_cmd.Step("assembly", "assembly", [])

    threads, memory, _ = run_cmd.slurm_resources_for(
        step, settings, {}, default_threads=32, default_memory_gb=500
    )

    assert (threads, memory) == (32, 500)


def test_module_resources_override_fallbacks_even_with_multiple_jobs():
    """A per-module Slurm grant is per job and must not be divided by the study's -j value."""
    settings = config_mod.Settings(
        slurm={"modules": {"assembly": {"threads": 32, "memory_gb": 500}}}
    )
    step = run_cmd.Step("assembly:A", "assembly", [])

    threads, memory, _ = run_cmd.slurm_resources_for(
        step, settings, {}, default_threads=8, default_memory_gb=32
    )
    assert (threads, memory) == (32, 500)


# --- Single-step submission, monitoring, and dry run ----------------------------------------


def test_a_slurm_step_records_its_job_id_before_waiting(tmp_path, monkeypatch):
    """The job ID must be saved before monitoring so a restarted driver can recover the job."""
    events = []

    class _State:
        def check_slurm_resume(self, step):
            """There is no previous Slurm job in this test."""
            return

        def running_slurm_job(self, step):
            return

        def start(self, step, command, grant=None):
            events.append(("start", grant["executor"]))
            return 1.0

        def record_slurm_job(self, step, job_id):
            events.append(("record", job_id))

        def finish(self, step, started, status, log=None):
            events.append(("finish", status))

    monkeypatch.setattr(run_cmd, "submit_job", lambda *a, **k: "123")

    def fake_wait(job_id, config):
        events.append(("wait", job_id))
        return SlurmResult("123", "COMPLETED", 0)

    monkeypatch.setattr(run_cmd, "wait_for_job", fake_wait)
    step = run_cmd.Step("assembly", "assembly", ["-o", str(tmp_path / "assembly")])
    _, status, _ = run_cmd.run_slurm_step(
        step,
        run_cmd.Layout(str(tmp_path)),
        8,
        None,
        [],
        False,
        _State(),
        32,
        run_cmd.SlurmConfig(),
    )

    assert status == COMPLETED
    assert events == [
        ("start", "slurm"),
        ("record", "123"),
        ("wait", "123"),
        ("finish", COMPLETED),
    ]


def test_a_failed_slurm_job_is_not_marked_completed(tmp_path, monkeypatch):
    """A failed module must never report a successful result to the workflow driver."""
    monkeypatch.setattr(run_cmd, "submit_job", lambda *a, **k: "123")
    monkeypatch.setattr(
        run_cmd,
        "wait_for_job",
        lambda *a, **k: SlurmResult("123", "FAILED", 1, "Slurm job failed"),
    )
    step = run_cmd.Step("assembly", "assembly", [])
    _, status, _ = run_cmd.run_slurm_step(
        step,
        run_cmd.Layout(str(tmp_path)),
        8,
        None,
        [],
        False,
        memory_gb=32,
        slurm_config=run_cmd.SlurmConfig(),
    )
    assert status == FAILED


def test_a_slurm_dry_run_never_submits_a_scheduler_job(tmp_path, monkeypatch, capsys):
    """A dry run prints the step command without contacting Slurm or starting a real job."""

    def unexpected_submission(*args, **kwargs):
        raise AssertionError("dry run submitted a Slurm job")

    monkeypatch.setattr(run_cmd, "submit_job", unexpected_submission)
    step = run_cmd.Step("assembly", "assembly", [])
    _, status, _ = run_cmd.run_slurm_step(
        step,
        run_cmd.Layout(str(tmp_path)),
        4,
        None,
        [],
        True,
        memory_gb=16,
        slurm_config=run_cmd.SlurmConfig(),
    )
    assert status == COMPLETED
    assert "[slurm plan]" in capsys.readouterr().out


# --- Stage-level concurrency and barriers --------------------------------------------------


def test_all_parallel_slurm_steps_are_checked_before_advancing(tmp_path, monkeypatch):
    """A failed first sample must not disappear when a later sample finishes successfully.

    Previously `futures` was reset inside the submission loop, so only the final future was
    inspected. That could allow assembly to start after an earlier QC job failed.
    """

    class _State:
        resume = False

        def skip_reason(self, step):
            return None

        class manifest:
            @staticmethod
            def step(name):
                return None

    calls = []
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(run_cmd, "_after_group", lambda *a: calls.append("advanced"))
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "comm", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "warning", lambda *a, **k: None)
    monkeypatch.setattr(
        run_cmd, "error", lambda message: (_ for _ in ()).throw(RuntimeError(message))
    )

    def fake_step(step, *args):
        calls.append(step.name)
        return step, FAILED if step.name == "read_qc:A" else COMPLETED, "log"

    monkeypatch.setattr(run_cmd, "run_slurm_step", fake_step)
    settings = config_mod.Settings(slurm={"threads": 4, "memory_gb": 8})
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", [], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", [], "read_qc"),
    ]

    with pytest.raises(RuntimeError, match="failed"):
        run_cmd.execute(
            steps,
            run_cmd.Layout(str(tmp_path)),
            run_cmd.Study([]),
            8,
            2,
            None,
            [],
            _State(),
            False,
            memory_gb=16,
            executor="slurm",
            settings=settings,
            cli_overrides={},
        )

    assert "read_qc:A" in calls and "read_qc:B" in calls
    assert "advanced" not in calls


def test_slurm_stages_keep_their_order_and_receive_per_job_resources(tmp_path, monkeypatch):
    """Independent QC jobs finish before assembly, with full Slurm grants for each module."""

    class _State:
        resume = False

        def skip_reason(self, step):
            return None

        class manifest:
            @staticmethod
            def step(name):
                return None

    calls = []
    monkeypatch.setattr(run_cmd, "_before_group", lambda *a: None)
    monkeypatch.setattr(
        run_cmd, "_after_group", lambda step, *a: calls.append(("after", step.module))
    )
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "comm", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "audit_outputs", lambda *a: None)

    def fake_step(step, layout, threads, config, flags, dry_run, state, memory, slurm_config):
        calls.append(("run", step.name, threads, memory, slurm_config.partition))
        return step, COMPLETED, "log"

    monkeypatch.setattr(run_cmd, "run_slurm_step", fake_step)
    settings = config_mod.Settings(
        slurm={
            "threads": 8,
            "memory_gb": 32,
            "partition": "shared",
            "modules": {"assembly": {"threads": 32, "memory_gb": 500, "partition": "bigmem"}},
        }
    )
    steps = [
        run_cmd.Step("read_qc:A", "read_qc", [], "read_qc"),
        run_cmd.Step("read_qc:B", "read_qc", [], "read_qc"),
        run_cmd.Step("assembly", "assembly", []),
    ]

    assert (
        run_cmd.execute(
            steps,
            run_cmd.Layout(str(tmp_path)),
            run_cmd.Study([]),
            8,
            2,
            None,
            [],
            _State(),
            True,
            memory_gb=32,
            executor="slurm",
            settings=settings,
            cli_overrides={},
        )
        == 0
    )

    qc = [entry for entry in calls if entry[:2] in (("run", "read_qc:A"), ("run", "read_qc:B"))]
    assert len(qc) == 2
    assert all(entry[2:] == (8, 32, "shared") for entry in qc)
    assert ("after", "read_qc") in calls
    assert ("run", "assembly", 32, 500, "bigmem") in calls
    assert calls.index(("after", "read_qc")) < calls.index(("run", "assembly", 32, 500, "bigmem"))


# --- CLI backend selection -----------------------------------------------------------------


def test_the_slurm_cli_selects_the_slurm_executor(tmp_path, monkeypatch):
    """`--executor slurm` must reach execute instead of silently using its local default.

    Regression: the driver previously parsed `--executor` but did not forward it to execute().
    """
    sheet = tmp_path / "samples.toml"
    sheet.write_text(
        '[settings]\nsteps = ["read_qc"]\n'
        '[[samples]]\nname = "A"\nr1 = "A.fastq"\nlayout = "single"\n'
    )
    captured = {}

    def fake_execute(*args, **kwargs):
        captured["executor"] = kwargs.get("executor")
        return 0

    monkeypatch.setattr(run_cmd, "execute", fake_execute)
    monkeypatch.setattr(run_cmd, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "comm", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "record_strategy", lambda *a, **k: None)
    monkeypatch.setattr(run_cmd, "load_settings", lambda *a, **k: config_mod.Settings())

    assert (
        run_cmd.main(
            [
                str(sheet),
                "-o",
                str(tmp_path / "results"),
                "--executor",
                "slurm",
                "--dry-run",
                "-t",
                "32",
                "-m",
                "500",
            ]
        )
        == 0
    )
    assert captured["executor"] == "slurm"
