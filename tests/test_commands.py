from metawrap2.commands import completion, config_cmd, doctor


def test_config_init_and_show(tmp_path, capsys):
    cfg = tmp_path / "c.toml"
    assert config_cmd.main(["init", "--path", str(cfg)]) == 0
    assert cfg.exists() and "[databases]" in cfg.read_text()
    # won't overwrite without --force
    assert config_cmd.main(["init", "--path", str(cfg)]) == 1
    assert config_cmd.main(["init", "--path", str(cfg), "--force"]) == 0
    # show reads it back
    assert config_cmd.main(["show", "--config", str(cfg)]) == 0
    out = capsys.readouterr().out
    assert "use_conda_envs" in out and str(cfg) in out


def test_completion_scripts(capsys):
    assert completion.main(["bash"]) == 0
    bash = capsys.readouterr().out
    assert "complete -F _metawrap2 metawrap2" in bash and "binning" in bash
    assert completion.main(["zsh"]) == 0
    assert "#compdef metawrap2" in capsys.readouterr().out


def test_doctor_reports_missing_when_no_envs(capsys, monkeypatch):
    # Force the "env absent" answer instead of depending on what is installed on the machine
    # running the tests: otherwise this passes on a fresh checkout and fails as soon as the
    # developer actually creates the envs (and probes real tools, making the suite slow).
    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rc = doctor.main(["binning", "--tools-only"])
    out = capsys.readouterr().out
    assert "MetaWrap2 module status" in out
    assert "metawrap2-binning" in out and "MISSING" in out
    assert rc == 1


def test_doctor_reports_ok_when_env_and_tools_are_healthy(capsys, monkeypatch):
    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: True)
    monkeypatch.setattr(doctor, "_probe_tool", lambda env, tool: (doctor.OK, "1.0"))
    rc = doctor.main(["binning", "--tools-only"])
    out = capsys.readouterr().out
    assert "tools OK" in out
    assert rc == 0


def test_doctor_parse_and_broken_signature():
    # broken-signature detection is what separates "broken" from "just exits nonzero"
    assert doctor._BROKEN_SIGNATURES.search("error while loading shared libraries: libfoo.so")
    assert not doctor._BROKEN_SIGNATURES.search("bwa 0.7.17-r1188")
