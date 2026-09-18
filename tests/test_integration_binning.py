"""End-to-end wiring test for the binning module using stub tools.

No real bioinformatics software is needed: we put tiny fake executables on PATH that
produce the files the module expects, disable conda envs, and run the module through its
CLI. This exercises arg parsing, command assembly, the command runner (stdout redirect),
output layout, and provenance writing - all the plumbing, without the science.
"""

import os
import stat

import pytest

from metawrap2.modules import binning

STUBS = {
    # bwa index -> touch; bwa mem -> print a minimal SAM header to stdout (redirected to .sam)
    "bwa": (
        "import sys\n"
        "if sys.argv[1] == 'mem':\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\nread1\\t0\\tc1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n')\n"
    ),
    # samtools sort -o OUT ...  -> create OUT ; samtools index -> noop
    "samtools": (
        "import sys\n"
        "if 'sort' in sys.argv and '-o' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-o')+1], 'w').write('BAM')\n"
    ),
    # jgi_summarize_bam_contig_depths --outputDepth PATH ... -> write a depth table
    "jgi_summarize_bam_contig_depths": (
        "import sys\n"
        "p = sys.argv[sys.argv.index('--outputDepth')+1]\n"
        "open(p,'w').write('contigName\\tcontigLen\\ttotalAvgDepth\\ts.bam\\n c1\\t1000\\t5\\t5\\n')\n"
    ),
    # metabat2 -o PREFIX ... -> write PREFIX.1.fa
    "metabat2": (
        "import os, sys\n"
        "pre = sys.argv[sys.argv.index('-o')+1]\n"
        "os.makedirs(os.path.dirname(pre), exist_ok=True)\n"
        "open(pre + '.1.fa', 'w').write('>c1\\nACGTACGT\\n')\n"
    ),
}


@pytest.fixture
def stub_path(tmp_path, monkeypatch):
    binp = tmp_path / "stubbin"
    binp.mkdir()
    for name, body in STUBS.items():
        f = binp / name
        f.write_text("#!/usr/bin/env python3\n" + body)
        f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)
    monkeypatch.setenv("PATH", str(binp) + os.pathsep + os.environ["PATH"])
    return binp


def test_binning_end_to_end_with_stubs(tmp_path, stub_path):
    # config with conda envs OFF, so tools run straight off PATH (our stubs)
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")

    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 300 + "\n")
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@read/1\nACGT\n+\nIIII\n")
    r2.write_text("@read/1\nACGT\n+\nIIII\n")
    out = tmp_path / "binning_out"

    rc = binning.main([
        "-a", str(asm), "-o", str(out), "--metabat2", "-t", "1",
        "--config", str(cfg), str(r1), str(r2),
    ])
    assert rc == 0

    # produced bins
    bins = out / "metabat2_bins"
    assert bins.is_dir() and any(f.endswith(".fa") for f in os.listdir(bins))
    # provenance written, with the real (unwrapped, since conda disabled) commands
    assert (out / "run_config.json").exists()
    cmds = (out / "run_commands.txt").read_text()
    assert "metabat2 -i" in cmds and "mamba run" not in cmds  # conda disabled -> no wrapping
