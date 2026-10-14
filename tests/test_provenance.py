import json

from metawrap2.provenance import RunRecorder


def _make_recorder(module, inputs):
    return RunRecorder(
        module=module,
        version="2.0.0",
        parameters={"threads": 4},
        config_file=None,
        conda_env="metawrap2-%s" % module,
        databases={},
        inputs=inputs,
        command_line=["metawrap2", module],
    )


def _finalize(rec, out, **kw):
    kw.setdefault("capture_versions", False)  # keep unit tests hermetic (no subprocess probes)
    rec.finalize(str(out), **kw)


def test_finalize_writes_three_files(tmp_path):
    rec = _make_recorder("binning", [])
    rec.record_command("metabat2 -i a.fa -t 4")
    out = tmp_path / "binning_out"
    _finalize(rec, out)
    for name in ("run_config.json", "run_commands.txt", "provenance.txt"):
        assert (out / name).exists(), name
    cfg = json.loads((out / "run_config.json").read_text())
    assert cfg["module"] == "binning"
    assert cfg["commands"] == ["metabat2 -i a.fa -t 4"]
    cmds = (out / "run_commands.txt").read_text()
    assert "metawrap2 binning" in cmds and "metabat2 -i a.fa -t 4" in cmds
