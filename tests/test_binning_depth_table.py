"""The shared depth table, and the flag that used to hand MaxBin2 a column of zeros.

MaxBin2 got its own ``jgi_summarize_bam_contig_depths --noIntraDepthVariance`` pass, which reads as a
formatting flag - MaxBin2 has no use for the variance columns - but in jgi 2.18 that code path
reports values that are not depths. On a 60-contig, two-sample probe it zeroed 59 of the 60 contigs
and gave the 60th 7867 where the real per-sample depths were 0.267 and 0.453, and rerunning it
unchanged returned 7874, so the number came from uninitialized memory.

Nothing failed, which is why it needed finding on purpose: a file of zeros is a well-formed
abundance file and MaxBin2 still bins on composition. These tests pin the two halves of the fix -
that one table serves both binners, and that MaxBin2's abundance files carry the mean-depth columns
of that table rather than the variances next to them.

The table is no longer built by a third-party binary at all: `metawrap2.depth` computes it from
``samtools depth``, for the reason this file exists - the number every binner depends on should be
code with tests, not an opaque pass whose output could only be checked by eye.
"""

from __future__ import annotations

import os
import stat

import pytest

from metawrap2 import command
from metawrap2.modules import binning

#: A two-sample jgi depth table, in the layout the real one has: a mean and a variance per sample,
#: and the variances deliberately unlike the means so a test cannot pass by reading either.
DEPTH = (
    "contigName\tcontigLen\ttotalAvgDepth\tsample_0\tsample_0-var\tsample_1\tsample_1-var\n"
    "c1\t5000\t8.5\t3.5\t91.1\t5.0\t92.2\n"
    "c2\t2000\t1.25\t0.25\t93.3\t1.0\t94.4\n"
)


def test_maxbin2_abundance_files_carry_the_means_not_the_variances(tmp_path, monkeypatch):
    """One file per sample, each the contig and that sample's mean depth.

    The variances sit in the column next to each mean and are an order of magnitude larger here, so
    an off-by-one in the split shows up as a number rather than as a missing file.
    """
    work = tmp_path / "work_files"
    work.mkdir()
    depth = work / binning.DEPTH_TABLE
    depth.write_text(DEPTH)

    # MaxBin2 itself is not what is under test, and it is the only thing _run_maxbin2 shells out to.
    calls = []
    monkeypatch.setattr(binning, "run", lambda cmd, **kw: calls.append(cmd))
    monkeypatch.setattr(binning, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(binning, "comm", lambda *a, **k: None)
    monkeypatch.setattr(binning, "dry_run", lambda: False)
    # No bins to collect afterwards is an error in the module, which is correct and not this test's
    # subject; the error function is replaced so the split can be inspected on its own.
    monkeypatch.setattr(binning, "error", lambda msg: (_ for _ in ()).throw(SystemExit(msg)))

    out = tmp_path / "out"
    args = type(
        "Args",
        (),
        {
            "output": str(out),
            "threads": 4,
            "min_len": 1000,
            "universal": False,
        },
    )()
    with pytest.raises(SystemExit):  # the stub produced no bins
        binning._run_maxbin2(args, None, str(work), "assembly.fa", str(depth))

    assert (work / "mb2_sample_0.txt").read_text() == "c1\t3.5\nc2\t0.25\n"
    assert (work / "mb2_sample_1.txt").read_text() == "c1\t5.0\nc2\t1.0\n"
    assert not (work / "mb2_sample_0-var.txt").exists(), "a variance is not a sample"
    listed = (work / "mb2_abund_list.txt").read_text().split()
    assert [os.path.basename(p) for p in listed] == ["mb2_sample_0.txt", "mb2_sample_1.txt"]


# ─── the shared table, through main() ────────────────────────────────────────────────────

#: jgi, but counting its own invocations into a file so a test can assert there was only one.
#: Emits the two-sample layout above, named for whatever BAMs it was handed.
_OTHER_STUBS = {
    "minimap2": (
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d')+1], 'w').write('IDX')\n"
        "else:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\nr\\t0\\tc1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n')\n"
    ),
    # samtools now does the pileup the depth table is computed from, so the stub has to answer
    # three more questions than it used to: what contigs the BAM header names, what its alignments
    # are, and what the per-base depth is. Every `samtools depth` call is recorded, which is how the
    # "one pass for both binners" test counts passes.
    "samtools": (
        "import os, sys\n"
        "a = sys.argv\n"
        "if 'sort' in a:\n"
        "    sys.stdin.read()\n"
        "    if '-o' in a:\n"
        "        open(a[a.index('-o')+1], 'w').write('BAM')\n"
        "elif 'view' in a and '-H' in a:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n@SQ\\tSN:c1\\tLN:5000\\n')\n"
        "elif 'view' in a:\n"
        "    sys.stdout.write('ALIGNMENTS\\n')\n"
        "elif 'depth' in a:\n"
        "    sys.stdin.read()\n"
        "    with open(os.environ['DEPTH_CALLS'], 'a') as fh:\n"
        "        fh.write(' '.join(a[1:]) + chr(10))\n"
        "    for pos in range(1, 5001):\n"
        "        sys.stdout.write('c1\\t%d\\t7\\n' % pos)\n"
    ),
    "metabat2": (
        "import os, sys\n"
        "pre = sys.argv[sys.argv.index('-o')+1]\n"
        "os.makedirs(os.path.dirname(pre), exist_ok=True)\n"
        "open(pre + '.1.fa', 'w').write('>c1\\nACGTACGT\\n')\n"
    ),
    "run_MaxBin.pl": (
        "import os, sys\n"
        "pre = sys.argv[sys.argv.index('-out')+1]\n"
        "os.makedirs(os.path.dirname(pre), exist_ok=True)\n"
        "open(pre + '.001.fasta', 'w').write('>c1\\nACGTACGT\\n')\n"
    ),
}


@pytest.fixture
def stubs(tmp_path, monkeypatch):
    binp = tmp_path / "stubbin"
    binp.mkdir()
    bodies = dict(_OTHER_STUBS)
    for name, body in bodies.items():
        f = binp / name
        f.write_text("#!/usr/bin/env python3\n" + body)
        f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)
    calls = tmp_path / "depth_calls.txt"
    monkeypatch.setenv("DEPTH_CALLS", str(calls))
    monkeypatch.setenv("PATH", str(binp) + os.pathsep + os.environ["PATH"])
    return calls


def test_resuming_with_a_new_binner_rebuilds_the_table_it_needs(tmp_path, stubs):
    """The table is a stage of its own now, so --resume has to be able to get back to it.

    The sequence that has to work: metaBAT2 alone, which deletes the BAMs when it is done, then
    ``--resume --maxbin2``, which has no BAMs to read and no table that was built from the reads this
    binner is being asked about. Both markers are cleared with the BAMs, so it realigns and rebuilds;
    had only the alignment marker been cleared, MaxBin2 would have been handed a table from a
    previous read set, which is the kind of wrong answer that does not look like one.

    ``--resume`` is set through the runner rather than passed in argv: it is a global flag that
    ``cli.main`` strips before the module ever sees it, so a module parser handed it rejects it.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n")
    r1, r2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"
    r1.write_text("@r/1\nACGT\n+\nIIII\n")
    r2.write_text("@r/1\nACGT\n+\nIIII\n")
    out = tmp_path / "out"
    base = ["-a", str(asm), "-o", str(out), "-t", "1", "-m", "4", "--config", str(cfg)]

    assert binning.main(base + ["--metabat2", str(r1), str(r2)]) == 0
    steps = out / ".metawrap2" / "steps"
    assert not (steps / "depth_table.done").exists(), "cleared with the BAMs it was built from"
    assert not (steps / "align_reads.done").exists()
    first = len([ln for ln in stubs.read_text().splitlines() if ln.strip()])
    assert first == 1

    command.set_resume(True)
    try:
        assert binning.main(base + ["--metabat2", "--maxbin2", str(r1), str(r2)]) == 0
    finally:
        command.set_resume(False)
    built = [ln for ln in stubs.read_text().splitlines() if ln.strip()]
    assert len(built) == 2, "the resumed run has to build the table again: %r" % built
    assert (out / "work_files" / "mb2_s.txt").read_text() == "c1\t7\n"
