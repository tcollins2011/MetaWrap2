import pytest

from metawrap2.command import CommandRunner, run


def test_success_and_stdout_redirect(tmp_path):
    out = tmp_path / "data.txt"
    rc = run(["python3", "-c", "print('hi')"], tool="echo", stdout_path=str(out))
    assert rc == 0
    assert out.read_text().strip() == "hi"  # stdout captured to the data file


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
