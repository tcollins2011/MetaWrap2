"""Paths with spaces reach each tool as one argument.

``run`` splits a command string with ``shlex``, and the templates were filled with plain
``str.format``, so ``-o "o5b out"`` gave FastQC ``--outdir .../o5b`` and it failed with "does not
exist".
"""

from __future__ import annotations

import os
import shlex

import pytest

from metawrap2 import command
from metawrap2.modules import _common, read_qc


@pytest.fixture(autouse=True)
def _history(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    command.runner.configure(resume=False, force=False, dry_run=False)
    yield
    command.runner.configure(resume=False, force=False, dry_run=False)


def test_each_value_stays_one_argument():
    filled = _common.fill_command(
        "fastqc -t {threads} --outdir {outdir} {files}",
        threads=4,
        outdir="/data/o5b out/pre-QC_report",
        files=_common.quoted_words(["/data/a b_1.fastq", "/data/c_2.fastq"]),
    )
    assert shlex.split(filled) == [
        "fastqc",
        "-t",
        "4",
        "--outdir",
        "/data/o5b out/pre-QC_report",
        "/data/a b_1.fastq",
        "/data/c_2.fastq",
    ]


def test_ordinary_values_are_unchanged():
    """Paths without special characters are inserted as they are, so recorded commands match."""
    assert (
        _common.fill_command("bwa index {assembly}", assembly="/runs/asm/final.contigs.fa")
        == "bwa index /runs/asm/final.contigs.fa"
    )


def test_raw_fragments_are_inserted_as_they_are():
    assert _common.fill_command("x {extra}", extra=_common.Raw("-a -b")) == "x -a -b"


def test_a_filled_command_survives_bash_c():
    """The aligner pipeline runs under ``bash -c``; shell quoting keeps the path whole there too."""
    align = _common.fill_command("minimap2 -ax sr {index} {reads}", index="i", reads="/d x/r.fq")
    script = _common.fill_command(
        _common.ALIGN_SORT, align=_common.Raw(align), tmp="/o x/tmp", threads=2, bam="/o x/s.bam"
    )
    words = shlex.split(script)
    assert "/d x/r.fq" in words and "/o x/s.bam" in words and "/o x/tmp" in words


def test_read_qc_gives_fastqc_the_whole_output_path(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, env=None, tool=None, **kwargs):
        calls.append((tool, shlex.split(cmd)))
        parts = shlex.split(cmd)
        if tool == "fastp":
            for flag in ("-o", "-O"):
                with open(parts[parts.index(flag) + 1], "w") as fh:
                    fh.write("@r\nACGT\n+\nIIII\n")
        return 0

    monkeypatch.setattr(read_qc, "run", fake_run)
    monkeypatch.setattr(read_qc, "env_for", lambda module, settings: None)
    r1, r2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    r2.write_text("@r\nTGCA\n+\nIIII\n")
    out = tmp_path / "o5b out"
    argv = ["-1", str(r1), "-2", str(r2), "-o", str(out), "-t", "4", "--skip-bmtagger"]
    assert read_qc.main(argv) == 0
    fastqc = [parts for tool, parts in calls if tool == "fastqc"]
    assert fastqc
    for parts in fastqc:
        outdir = parts[parts.index("--outdir") + 1]
        assert outdir.startswith(str(out))
        assert os.path.isdir(outdir)
    (fastp,) = [parts for tool, parts in calls if tool == "fastp"]
    assert str(out / "trimmed_1.fastq") in fastp
    assert str(out / "fastp.json") in fastp


def test_the_odd_path_warning_says_what_happens(monkeypatch):
    warned = []
    monkeypatch.setattr(_common, "warning", warned.append)
    _common._warn_odd_paths(["/data/o5b out"])
    (message,) = warned
    assert "/data/o5b out" in message
    assert "one argument" in message
    assert "Consider renaming" not in message


def test_an_empty_value_adds_no_argument():
    """``{cuda}`` unset must not give VAMB an empty argument."""
    assert shlex.split(_common.fill_command("vamb -p {t} {cuda}", t=2, cuda="")) == [
        "vamb",
        "-p",
        "2",
    ]


#: Modules outside metawrap2.modules whose COMMANDS templates are filled with fill_command.
TEMPLATE_MODULES = ("aligners", "checkm", "gunc")


def _module_templates():
    import importlib
    import pkgutil
    import string

    from metawrap2 import modules

    names = ["modules." + info.name for info in pkgutil.iter_modules(modules.__path__)]
    for module_name in names + list(TEMPLATE_MODULES):
        mod = importlib.import_module("metawrap2." + module_name)
        for name, value in vars(mod).items():
            if name.isupper() and isinstance(value, str) and "{" in value:
                fields = {f for _, f, _, _ in string.Formatter().parse(value) if f}
                if fields:
                    yield "%s.%s" % (module_name, name), value, fields


def test_no_module_template_splits_a_path_with_a_space():
    """Every COMMANDS template keeps a filled-in path with a space as one piece."""
    checked = 0
    for label, template, fields in _module_templates():
        # Separators are a single character with no space (metawrap2.binsplit.SEPARATOR).
        values = {f: ("~" if f == "separator" else "/p q/f") for f in fields}
        try:
            words = shlex.split(_common.fill_command(template, **values))
        except (KeyError, IndexError, ValueError) as exc:  # pragma: no cover - reported below
            raise AssertionError("%s: %s" % (label, exc)) from exc
        assert "/p" not in words and not [w for w in words if w.startswith("q/f")], label
        checked += 1
    assert checked > 50
    labels = {label for label, _, _ in _module_templates()}
    assert {"gunc.GUNC_RUN", "checkm.CHECKM2_PREDICT", "checkm.CHECKM1_LINEAGE_WF"} <= labels


def _capture_run(monkeypatch, module):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(shlex.split(cmd))
        return 0

    monkeypatch.setattr(module, "run", fake_run)
    return calls


def test_checkm2_gets_whole_paths(tmp_path, monkeypatch):
    from metawrap2 import checkm

    calls = _capture_run(monkeypatch, checkm)
    command.runner.configure(resume=False, force=False, dry_run=True)
    bins = str(tmp_path / "my bins")
    checkm.run_checkm2(bins, threads=2, env=None, database="/db dir/x.dmnd", ttable=11)
    (argv,) = calls
    assert argv[argv.index("--input") + 1] == bins
    assert argv[argv.index("--output-directory") + 1] == bins + ".checkm"
    assert argv[argv.index("--database_path") + 1] == "/db dir/x.dmnd"
    assert argv[argv.index("--ttable") + 1] == "11"


def test_checkm1_gets_whole_paths(monkeypatch):
    from metawrap2 import checkm

    calls = _capture_run(monkeypatch, checkm)
    checkm._checkm1("/p q/bins", "/p q/out", "/t x", 4, 64, None, True, None, 1)
    (argv,) = calls
    assert argv[argv.index("-x") + 2 : argv.index("-x") + 4] == ["/p q/bins", "/p q/out"]
    assert argv[argv.index("--tmpdir") + 1] == "/t x"
    assert argv[-1] == "--reduced_tree"


def test_gunc_gets_whole_paths(monkeypatch):
    from metawrap2 import gunc

    calls = _capture_run(monkeypatch, gunc)
    command.runner.configure(resume=False, force=False, dry_run=True)
    gunc.run_gunc("/p q/bins", threads=2, db="/db x/gunc.dmnd", env=None)
    (argv,) = calls
    assert argv[argv.index("--input_dir") + 1] == "/p q/bins"
    assert argv[argv.index("--db_file") + 1] == "/db x/gunc.dmnd"
    assert argv[argv.index("--out_dir") + 1] == "/p q/bins.gunc"


def test_an_external_benchmark_arm_gets_whole_paths(monkeypatch):
    from metawrap2.benchmark import runner as bench
    from metawrap2.benchmark.registry import External, Variant

    monkeypatch.setattr(bench._runner, "env_prefix", lambda env: [])
    arm = Variant(
        name="mw1",
        external=External(
            command="metawrap binning -o {output}/out -t {threads} -a {assembly} {reads}",
            bins="out/final_bins",
            env="e",
            version="1",
        ),
    )
    inputs = bench.Inputs(
        pairs=[("s", "/r x/s_1.fastq", "/r x/s_2.fastq")], assembly="/a b/gold.fa"
    )
    argv = bench.external_argv(arm, inputs, "/o p/mw1", 2, 8)
    assert argv[argv.index("-a") + 1] == "/a b/gold.fa"
    assert argv[-2:] == ["/r x/s_1.fastq", "/r x/s_2.fastq"]
    assert argv[argv.index("-o") + 1] == "/o p/mw1/external/out"
