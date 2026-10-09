"""End-to-end wiring test for the binning module using stub tools.

No real bioinformatics software is needed: we put tiny fake executables on PATH that
produce the files the module expects, disable conda envs, and run the module through its
CLI. This exercises arg parsing, command assembly, the command runner (stdout redirect),
output layout, and provenance writing - all the plumbing, without the science.
"""

import os
import stat

import pytest

from metawrap2.command import ToolError
from metawrap2.constants import BIN_EXTENSION
from metawrap2.modules import binning

STUBS = {
    # minimap2 is the module's default aligner. -d writes the index; otherwise print a minimal SAM
    # to stdout, which the module pipes into samtools sort rather than to a file.
    "minimap2": (
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d') + 1], 'w').write('MMI')\n"
        "else:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\nread1\\t0\\tc1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n')\n"
    ),
    # --aligner bwa / bwa-mem2. Both index into files named after the assembly, so "index" needs
    # to write nothing the module then looks for; "mem" prints SAM, same as minimap2 above.
    "bwa": (
        "import sys\n"
        "if sys.argv[1] == 'mem':\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\nread1\\t0\\tc1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n')\n"
    ),
    "bwa-mem2": (
        "import sys\n"
        "if sys.argv[1] == 'mem':\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\nread1\\t0\\tc1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n')\n"
    ),
    # samtools sort -o OUT ...  -> create OUT ; samtools index -> noop
    #
    # The sort must drain stdin, as the real one does: binning pipes the aligner into it, and a
    # sort that exits while the aligner is still writing hands the aligner a SIGPIPE, which
    # pipefail then correctly reports as a failed run. That is a faithful stub's job, not a bug.
    # samtools wears four hats here, because the module uses it for four things:
    #
    #   sort -o OUT      create OUT, draining stdin first
    #   view -H BAM      print a header, which is where the depth table gets its contig lengths
    #   view -u BAM      print "alignments" for the pileup to read (the identity filter lives here)
    #   depth -a ... -   print per-base depth rows
    #
    # The last three are new: the contig depth table used to be built by shelling out to
    # jgi_summarize_bam_contig_depths, so a stub only had to write the finished file. It is computed
    # in-process now (see metawrap2.depth), which means this test exercises the real table assembly
    # and the stub has to supply believable pileup input rather than a believable answer.
    #
    # The sort must drain stdin, as the real one does: binning pipes the aligner into it, and a
    # sort that exits while the aligner is still writing hands the aligner a SIGPIPE, which
    # pipefail then correctly reports as a failed run. That is a faithful stub's job, not a bug.
    "samtools": (
        "import sys\n"
        "a = sys.argv\n"
        "if 'sort' in a:\n"
        "    sys.stdin.read()\n"
        "    if '-o' in a:\n"
        "        open(a[a.index('-o')+1], 'w').write('BAM')\n"
        "elif 'view' in a and '-H' in a:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n@SQ\\tSN:c1\\tLN:1000\\n')\n"
        "elif 'view' in a:\n"
        "    sys.stdout.write('ALIGNMENTS\\n')\n"
        "elif 'depth' in a:\n"
        "    sys.stdin.read()\n"
        "    for pos in range(1, 1001):\n"
        "        sys.stdout.write('c1\\t%d\\t5\\n' % pos)\n"
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

    rc = binning.main(
        [
            "-a",
            str(asm),
            "-o",
            str(out),
            "--metabat2",
            "-t",
            "1",
            "--config",
            str(cfg),
            str(r1),
            str(r2),
        ]
    )
    assert rc == 0

    # produced bins
    bins = out / "metabat2_bins"
    assert bins.is_dir() and any(f.endswith(BIN_EXTENSION) for f in os.listdir(bins))
    # provenance written, with the real (unwrapped, since conda disabled) commands
    assert (out / "run_config.json").exists()
    cmds = (out / "run_commands.txt").read_text()
    assert "metabat2 -i" in cmds and "mamba run" not in cmds  # conda disabled -> no wrapping
    # the default aligner, piped into the sort rather than routed through a SAM on disk
    assert "minimap2 -ax sr" in cmds and "| samtools sort" in cmds
    assert not any(f.endswith(".sam") for f in os.listdir(out / "work_files"))


def _binning_argv(tmp_path, cfg, out, extra=()):
    """The arg list for a one-sample, one-binner run, plus whatever *extra* flags are under test."""
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 300 + "\n")
    for suffix in ("_1", "_2"):
        (tmp_path / ("s%s.fastq" % suffix)).write_text("@read/1\nACGT\n+\nIIII\n")
    return (
        ["-a", str(asm), "-o", str(out), "--metabat2", "-t", "1", "--config", str(cfg)]
        + list(extra)
        + [str(tmp_path / "s_1.fastq"), str(tmp_path / "s_2.fastq")]
    )


@pytest.mark.parametrize("aligner", sorted(binning.ALIGNERS))
def test_every_aligner_choice_produces_coverage(tmp_path, stub_path, aligner):
    """--aligner has to reach a BAM by all three routes, not just the default one.

    The stub for each is keyed off that aligner's own index flag, so an aligner wired to the wrong
    index path - minimap2's separate .mmi against bwa's prefix-is-the-assembly - fails here rather
    than silently mapping against an index that was never built.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    out = tmp_path / "binning_out"

    rc = binning.main(_binning_argv(tmp_path, cfg, out, ["--aligner", aligner, "--keep-bams"]))

    assert rc == 0
    assert (out / "work_files" / "s.bam").exists()
    cmds = (out / "run_commands.txt").read_text()
    spec = binning.ALIGNERS[aligner]
    assert spec["index_cmd"].split()[0] in cmds
    assert spec["paired"].split()[0] in cmds and "| samtools sort" in cmds


def test_the_bams_are_gone_once_the_depth_tables_exist(tmp_path, stub_path):
    """The BAM is scratch, not output: what the binners needed out of it is the depth table.

    Also checks the resume marker went with it. Leaving align_reads marked done with no BAM on disk
    would let a later ``--resume --maxbin2`` skip alignment and bin on an empty depth table, which
    is the one way deleting them could do real harm.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    out = tmp_path / "binning_out"

    assert binning.main(_binning_argv(tmp_path, cfg, out)) == 0

    work = out / "work_files"
    assert not [f for f in os.listdir(work) if f.endswith(".bam")]
    assert (work / binning.DEPTH_TABLE).exists(), "the depth table is the output that matters"
    assert (work / "assembly.fa").exists(), "the binners are handed this path; it has to survive"
    assert not (out / ".metawrap2" / "steps" / "align_reads.done").exists()


def test_keep_bams_keeps_them(tmp_path, stub_path):
    """The escape hatch, for anyone who wants to look at an alignment afterwards."""
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    out = tmp_path / "binning_out"

    assert binning.main(_binning_argv(tmp_path, cfg, out, ["--keep-bams"])) == 0

    assert (out / "work_files" / "s.bam").exists()
    assert (out / ".metawrap2" / "steps" / "align_reads.done").exists()


def test_a_failed_alignment_is_not_hidden_by_a_successful_sort(tmp_path, stub_path):
    """pipefail, tested rather than assumed.

    ``samtools sort`` fed a truncated stream still exits 0 and still writes a BAM, so without
    ``set -o pipefail`` bash reports the sort's status and a crashed aligner reads as a clean run
    that produced an empty coverage table. This is the failure mode that makes the flag necessary.
    """
    (stub_path / "minimap2").write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d') + 1], 'w').write('MMI')\n"
        "    sys.exit(0)\n"
        "sys.stderr.write('out of memory\\n')\n"
        "sys.exit(1)\n"
    )
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    out = tmp_path / "binning_out"

    with pytest.raises(ToolError):
        binning.main(_binning_argv(tmp_path, cfg, out))
