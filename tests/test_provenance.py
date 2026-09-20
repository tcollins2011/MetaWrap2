import json

from metawrap2 import provenance
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


def test_multiple_mw2_inputs_dedup_shared_ancestor(tmp_path):
    asm = _make_recorder("assembly", [])
    asm_dir = tmp_path / "asm"
    _finalize(asm, asm_dir)

    def child(mod, path):
        r = _make_recorder(mod, [str(asm_dir / "final_assembly.fasta")])
        r_dir = tmp_path / path
        _finalize(r, r_dir, input_paths=[str(asm_dir / "final_assembly.fasta")])
        return r_dir

    a = child("binning", "binsA_run")
    b = child("binning", "binsB_run")

    ref = _make_recorder("bin_refinement", [str(a), str(b)])
    ref_dir = tmp_path / "refined"
    _finalize(ref, ref_dir, input_paths=[str(a), str(b)])
    modules = [s["module"] for s in json.loads((ref_dir / CONFIG_NAME).read_text())["lineage"]]
    assert modules.count("assembly") == 1
    assert modules.count("binning") == 2
    assert modules[-1] == "bin_refinement"


def test_find_mw2_root_walks_up(tmp_path):
    root = tmp_path / "binning_out"
    (root / "metabat2_bins").mkdir(parents=True)
    (root / CONFIG_NAME).write_text("{}")
    assert provenance.find_mw2_root(str(root / "metabat2_bins" / "bin.1.fa")) == str(root)
    assert provenance.find_mw2_root(str(tmp_path / "nowhere.fa")) is None


def test_parse_env_tool():
    assert provenance._parse_env_tool(
        "mamba run --no-capture-output -n metawrap2-binning metabat2 -i a.fa"
    ) == ("metawrap2-binning", "metabat2")
    assert provenance._parse_env_tool("cd /x && mamba run -n e salmon quant") == ("e", "salmon")
    assert provenance._parse_env_tool("bash -c 'a | b'") is None
    assert provenance._parse_env_tool("prokka --outdir x") == (None, "prokka")


def test_capture_versions_graceful_when_tool_absent():
    # a tool that isn't installed -> no version recorded, no error
    versions = provenance.capture_tool_versions(
        ["mamba run -n metawrap2-binning definitely_not_a_real_tool_xyz --version"]
    )
    assert versions == {}
