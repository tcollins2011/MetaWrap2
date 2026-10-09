import pytest

from metawrap2.command import CommandRunner, run


def test_success_and_stdout_redirect(tmp_path):
    out = tmp_path / "data.txt"
    rc = run(["python3", "-c", "print('hi')"], tool="echo", stdout_path=str(out))
    assert rc == 0
    assert out.read_text().strip() == "hi"  # stdout captured to the data file


def test_env_prefix_uses_mamba(monkeypatch):
    # --no-capture-output is conda-only; mamba 2.x mis-parses it and every command dies with
    # "exec: --: invalid option". So the flag is included only when `run --help` advertises it.
    from metawrap2 import command as command_mod

    r = CommandRunner(env_manager="mamba")
    # Whether mamba is installed on the machine running the tests is not what this asserts.
    monkeypatch.setattr(command_mod.shutil, "which", lambda name: "/opt/conda/bin/" + name)
    monkeypatch.setattr(command_mod, "_run_flags", lambda mgr: ["--no-capture-output"])
    assert r.env_prefix("metawrap2-bins") == [
        "mamba",
        "run",
        "--no-capture-output",
        "-n",
        "metawrap2-bins",
    ]

    monkeypatch.setattr(command_mod, "_run_flags", lambda mgr: [])
    assert r.env_prefix("metawrap2-bins") == ["mamba", "run", "-n", "metawrap2-bins"]
    assert r.env_prefix(None) == []


def test_skip_validation_flag_bypasses_input_checks():
    """--skip-validation makes validate_inputs a no-op, even for a would-fail check."""
    from metawrap2 import command, validate
    from metawrap2.modules import _common

    command.set_skip_validation(True)
    try:
        # A nonexistent file would normally abort the run; with the flag set it is skipped.
        _common.validate_inputs("test", [(validate.check_fastq, "/no/such/file.fastq", "reads")])
    finally:
        command.set_skip_validation(False)


def test_a_missing_env_manager_is_named_before_anything_is_spawned(monkeypatch):
    """An absent mamba must not reach Popen, where it reads as a missing *tool*.

    The bare ``FileNotFoundError: 'mamba'`` a subprocess raises points at the command being run
    rather than at the shell that has no conda on PATH, which is the actual cause - and in a
    benchmark, where a failed variant is recorded and the run continues, it produces a directory
    of identically failed runs instead of one message.
    """
    from metawrap2 import command as command_mod

    r = CommandRunner(env_manager="mamba")
    monkeypatch.setattr(command_mod.shutil, "which", lambda name: None)
    with pytest.raises(command_mod.EnvManagerMissing) as caught:
        r.env_prefix("metawrap2-reads")
    message = str(caught.value)
    assert "not on PATH" in message
    assert "miniforge3/bin" in message  # says how to fix it, not just what is wrong
    assert "use_conda_envs" in message  # and names the way to run without envs at all

    # A dry run is precisely for a machine where nothing is installed yet, so it still prints.
    dry = CommandRunner(env_manager="mamba", dry_run=True)
    assert dry.env_prefix("metawrap2-reads")[:2] == ["mamba", "run"]


def test_a_tool_in_a_conda_env_does_not_import_from_the_users_own_site_packages(tmp_path):
    """``~/.local`` comes before a conda env on sys.path and silently overrides what the env pinned.

    It cost a real debugging session: AMBER's env pins matplotlib below 3.9 for an API AMBER calls,
    and a matplotlib 3.11 in ~/.local was imported instead, so AMBER died after computing its metrics
    with ``Axes.boxplot() got an unexpected keyword argument 'labels'``.
    """
    r = CommandRunner(env_manager="mamba")
    assert r._child_env(True)["PYTHONNOUSERSITE"] == "1"


def test_a_tool_run_without_a_conda_env_keeps_the_environment_it_was_given():
    """No env, no pin to protect - and a tool the user pip-installed into ~/.local is then exactly
    the tool they meant to run."""
    assert CommandRunner()._child_env(False) is None


def test_an_interrupt_between_spawning_and_waiting_still_kills_the_tool(monkeypatch):
    """The window was two statements wide, and cost a whole machine when it was hit.

    The tool runs in its own process group, so only this process kills it. An interrupt arriving
    after Popen returned but before the wait it was guarded around left the tool running with every
    thread it had been granted, orphaned by the run that was exiting - and whatever the user started
    next then fought it for the machine.
    """
    import threading

    monkeypatch.setattr(
        threading.Thread, "start", lambda self: (_ for _ in ()).throw(KeyboardInterrupt())
    )
    killed = []
    monkeypatch.setattr(CommandRunner, "_terminate", staticmethod(lambda proc: killed.append(proc)))

    r = CommandRunner()
    with pytest.raises(KeyboardInterrupt):
        r.run("sleep 30", tool="sleep")

    assert len(killed) == 1, "the spawned tool must be terminated, not left behind"
    assert killed[0].pid > 0
    killed[0].kill()  # _terminate was stubbed out, so this test cleans up after itself
