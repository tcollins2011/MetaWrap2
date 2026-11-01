"""Global flags given before `run`, `cohort` or `benchmark` reach every step's process."""

import itertools

import pytest

from metawrap2 import cli, command, limits


@pytest.fixture(autouse=True)
def _clean_state():
    yield
    limits.clear_overrides()
    for setter in (
        command.set_dry_run,
        command.set_force,
        command.set_resume,
        command.set_verbose_logs,
        command.set_skip_space_check,
        command.set_skip_validation,
    ):
        setter(False)
    command.set_heartbeat_seconds(None)


TYPED = [
    "--max-threads", "32",
    "--max-memory-gb", "1000",
    "--max-disk-fraction", "0.9",
    "--skip-space-check",
    "--skip-validation",
    "--verbose-logs",
    "--heartbeat-seconds", "15",
    "--gpu-strict",
    "--resume",
]  # fmt: skip


def test_every_typed_global_flag_is_passed_to_children():
    assert cli._apply_global_flags(TYPED + ["cohort", "inputs"]) == ["cohort", "inputs"]
    inherited = cli.inherited_global_flags()
    for flag in (
        "--skip-space-check",
        "--skip-validation",
        "--verbose-logs",
        "--gpu-strict",
        "--resume",
    ):
        assert flag in inherited
    pairs = dict(itertools.pairwise(inherited))
    assert pairs["--max-threads"] == "32"
    assert float(pairs["--max-memory-gb"]) == 1000
    assert float(pairs["--max-disk-fraction"]) == 0.9
    assert float(pairs["--heartbeat-seconds"]) == 15


def test_inherited_flags_parse_back_to_the_same_state():
    cli._apply_global_flags(TYPED + ["cohort"])
    first = (dict(limits.overrides()), command.runner.skip_space_check)
    inherited = cli.inherited_global_flags()
    limits.clear_overrides()
    command.set_skip_space_check(False)
    assert cli._apply_global_flags(inherited + ["read_qc"]) == ["read_qc"]
    assert (dict(limits.overrides()), command.runner.skip_space_check) == first


def test_no_flags_means_nothing_inherited():
    assert cli.inherited_global_flags() == []


def test_a_study_step_command_carries_the_global_flags():
    from metawrap2.commands import run as run_cmd

    cli._apply_global_flags(["--max-disk-fraction", "0.9", "--skip-space-check", "cohort"])
    flags = cli.inherited_global_flags()
    step = run_cmd.Step(name="read_qc:s1", module="read_qc", argv=["-o", "out"])
    argv = step.command(16, None, flags, memory_gb=120)
    head = argv[: argv.index("read_qc")]
    assert "--skip-space-check" in head and "--max-disk-fraction" in head
    assert "-t" in argv and "16" in argv


# --- from the tests of pull request #7 ------------------------------------------------------


@pytest.mark.parametrize(
    "key,value,flag",
    [
        ("max_threads", 32, "--max-threads"),
        ("max_thread_fraction", 0.8, "--max-thread-fraction"),
        ("max_memory_gb", 1000, "--max-memory-gb"),
        ("max_memory_fraction", 0.9, "--max-memory-fraction"),
        ("max_disk_fraction", 0.9, "--max-disk-fraction"),
        ("max_gpus", 0, "--max-gpus"),
        ("disk_check_interval", 30, "--disk-check-interval"),
    ],
)
def test_each_budget_override_is_inherited(key, value, flag):
    limits.set_override(key, value)
    inherited = cli.inherited_global_flags()
    assert inherited[:1] == [flag] and float(inherited[1]) == value and len(inherited) == 2


def test_gpu_strict_off_is_not_inherited():
    # --gpu-strict is a switch; it takes no value.
    limits.set_override("gpu_strict", False)
    assert "--gpu-strict" not in cli.inherited_global_flags()


def test_a_zero_heartbeat_is_inherited():
    # 0 turns heartbeats off, so it must not be dropped as if unset.
    command.set_heartbeat_seconds(0)
    inherited = cli.inherited_global_flags()
    assert inherited[:1] == ["--heartbeat-seconds"] and float(inherited[1]) == 0


def test_dry_run_and_force_are_inherited():
    command.set_dry_run(True)
    command.set_force(True)
    assert {"--dry-run", "--force"} <= set(cli.inherited_global_flags())


def test_a_step_command_keeps_its_module_arguments_after_the_global_flags():
    from metawrap2.commands import run as run_cmd

    step = run_cmd.Step(name="assembly:A", module="assembly", argv=["-o", "study/ASSEMBLY/A"])
    argv = step.command(16, "metawrap2.toml", ["--resume"], memory_gb=60)
    module = argv.index("assembly")
    assert argv.index("--resume") < module
    assert argv[argv.index("-t") + 1] == "16"
    assert argv[argv.index("-m", module) + 1] == "60"
    assert argv[argv.index("--config") + 1] == "metawrap2.toml"
