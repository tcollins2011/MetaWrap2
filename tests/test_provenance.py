import json

from metawrap2.provenance import CONFIG_NAME, RunRecorder


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


def test_lineage_propagates_from_mw2_input(tmp_path):
    up = _make_recorder("assembly", [])
    up.record_command("metaspades.py -o asm")
    asm_dir = tmp_path / "assembly_out"
    _finalize(up, asm_dir)

    down = _make_recorder("binning", [str(asm_dir / "final_assembly.fasta")])
    down.record_command("metabat2 ...")
    bin_dir = tmp_path / "binning_out"
    _finalize(down, bin_dir, input_paths=[str(asm_dir / "final_assembly.fasta")])

    cfg = json.loads((bin_dir / CONFIG_NAME).read_text())
    modules = [s["module"] for s in cfg["lineage"]]
    assert modules == ["assembly", "binning"]
    prov = (bin_dir / "provenance.txt").read_text()
    assert "assembly" in prov and "binning" in prov and "metaspades.py" in prov
