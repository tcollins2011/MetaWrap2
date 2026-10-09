from metawrap2.commands import config_cmd, doctor


def test_mamba_missing_from_path_is_a_message_not_a_traceback(monkeypatch, capsys):
    """The exception exists to explain a PATH problem in two sentences; a stack buries that.

    It was raised where the command is built precisely so the run stops early with something
    actionable, but nothing caught it, so the user got the message underneath a traceback - which
    reads as a crash, the exact confusion it was written to prevent.
    """
    import importlib

    from metawrap2 import cli
    from metawrap2.command import EnvManagerMissing

    module = importlib.import_module("metawrap2.commands.doctor")

    def boom(_argv):
        raise EnvManagerMissing("mamba", "metawrap2-bins")

    monkeypatch.setattr(module, "main", boom)
    monkeypatch.setattr(cli, "importlib", importlib)
    monkeypatch.setitem(cli.MODULES, "doctor", ("metawrap2.commands.doctor", "check the install"))

    assert cli.main(["doctor"]) == 2
    err = capsys.readouterr().err
    assert "Traceback" not in err
    assert "mamba is not on PATH" in err
    # Both ways out are named, since either one is a valid answer.
    assert "miniforge3/bin" in err and "use_conda_envs = false" in err


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


def test_doctor_reports_missing_when_no_envs(capsys, monkeypatch):
    # Force the "env absent" answer instead of depending on what is installed on the machine
    # running the tests: otherwise this passes on a fresh checkout and fails as soon as the
    # developer actually creates the envs (and probes real tools, making the suite slow).
    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rc = doctor.main(["binning", "--tools-only"])
    out = capsys.readouterr().out
    assert "MetaWrap2 module status" in out
    assert "metawrap2-bins" in out and "MISSING" in out
    assert rc == 1


def test_doctor_probes_an_opt_in_env_asked_for_by_its_yaml_name(monkeypatch):
    """`install-env bins` verifies by calling doctor with that name, which is an env, not a module.

    It is not a key of MODULE_ENVS or of MODULE_TOOLS, so it used to resolve to no known tools and
    report a freshly built env as "0/0 tools OK" - a pass that had probed nothing at all. An env
    several modules share is only healthy if every one of them can run, so the probe list is their
    union.
    """
    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: True)
    probed = []
    monkeypatch.setattr(
        doctor, "_probe_tool", lambda env, tool: probed.append((env, tool)) or (doctor.OK, "runs")
    )
    result = doctor._check_module("bins")
    assert result["env"] == "metawrap2-bins"
    probed_tools = [tool for env, tool in probed]
    assert all(env == "metawrap2-bins" for env, _ in probed)
    # binning's binners, coverage's aligners, dereplicate's skani and quant_bins' three quantifiers:
    # every module that resolves to this env, not just the first one.
    for expected in ("metabat2", "concoct", "minimap2", "fairy", "skani", "sylph", "coverm"):
        assert expected in probed_tools
    assert result["detail"] == "%d/%d tools OK" % (len(probed), len(probed))
