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
    monkeypatch.setattr(command_mod, "_run_flags", lambda mgr: ["--no-capture-output"])
    assert r.env_prefix("metawrap2-binning") == [
        "mamba",
        "run",
        "--no-capture-output",
        "-n",
        "metawrap2-binning",
    ]

    monkeypatch.setattr(command_mod, "_run_flags", lambda mgr: [])
    assert r.env_prefix("metawrap2-binning") == ["mamba", "run", "-n", "metawrap2-binning"]
    assert r.env_prefix(None) == []
