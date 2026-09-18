import pytest

from metawrap2.command import CommandRunner, ToolError, run


def test_success_and_stdout_redirect(tmp_path):
    out = tmp_path / "data.txt"
    rc = run(["python3", "-c", "print('hi')"], tool="echo", stdout_path=str(out))
    assert rc == 0
    assert out.read_text().strip() == "hi"  # stdout captured to the data file


def test_failure_raises_toolerror_with_stderr_tail():
    with pytest.raises(ToolError) as exc:
        run(["python3", "-c", "import sys; sys.stderr.write('boom\\n'); sys.exit(3)"],
            tool="thing", hint="check inputs")
    err = exc.value
    assert err.returncode == 3
    assert "boom" in err.log_tail          # stderr is captured for the error message
    assert "check inputs" in str(err) and "thing failed" in str(err)


def test_check_false_returns_code():
    assert run(["python3", "-c", "import sys; sys.exit(7)"], check=False) == 7


def test_env_prefix_uses_mamba():
    r = CommandRunner(env_manager="mamba")
    assert r.env_prefix("metawrap2-binning") == ["mamba", "run", "--no-capture-output", "-n", "metawrap2-binning"]
    assert r.env_prefix(None) == []


def test_dry_run_records_but_does_not_execute(tmp_path):
    r = CommandRunner(dry_run=True)

    class Rec:
        def __init__(self): self.cmds = []
        def record_command(self, c): self.cmds.append(c)

    rec = Rec()
    r.configure(recorder=rec)
    sentinel = tmp_path / "should_not_exist"
    rc = r.run(["python3", "-c", "open(%r,'w').close()" % str(sentinel)], tool="py")
    assert rc == 0
    assert not sentinel.exists()           # not executed
    assert rec.cmds and "python3" in rec.cmds[0]  # but recorded


def test_run_logs_tee_stdout_and_stderr(tmp_path):
    r = CommandRunner()
    so = tmp_path / "run.stdout"
    se = tmp_path / "run.stderr"
    r.configure(run_stdout_path=str(so), run_stderr_path=str(se))
    r.run(["python3", "-c", "import sys; print('to-out'); sys.stderr.write('to-err\\n')"], tool="py")
    assert "to-out" in so.read_text()
    assert "to-err" in se.read_text()


def test_recorder_notes_redirect(tmp_path):
    class Rec:
        def __init__(self): self.cmds = []
        def record_command(self, c): self.cmds.append(c)
    r = CommandRunner(dry_run=True)
    rec = Rec()
    r.configure(recorder=rec)
    r.run(["bwa", "mem", "a", "b"], env="metawrap2-binning", stdout_path="x.sam")
    assert "mamba run" in rec.cmds[0] and "> x.sam" in rec.cmds[0]
