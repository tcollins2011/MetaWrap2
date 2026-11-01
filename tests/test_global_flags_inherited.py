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
