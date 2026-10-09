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


def test_the_depth_table_is_not_built_by_shelling_out_to_jgi_any_more():
    """No jgi command, and no second pass for MaxBin2 to get a broken one from."""
    assert not hasattr(binning, "JGI_DEPTH"), "the shell-out is gone, not just unused"
    assert not hasattr(binning, "JGI_DEPTH_MAXBIN"), "the second pass is gone, not just unused"
    # And the replacement keeps the two parameters the broken flag was about: the variance is
    # computed rather than suppressed, and the identity filter is still what decides what counts.
    # The filter's *value* is a measured default that has since moved off jgi's 97 - that belongs to
    # test_depth.py; what matters here is that writing our own table did not quietly drop either
    # parameter, which is the way this replacement could have gone wrong without failing anything.
    from metawrap2 import depth as depth_module

    assert 0 < depth_module.PERCENT_IDENTITY <= 100
    assert depth_module._identity_expression(depth_module.PERCENT_IDENTITY), "still a real filter"
    assert depth_module.EDGE_BASES == 75


@pytest.mark.parametrize(
    "header,expected",
    [
        (["contigName", "contigLen", "totalAvgDepth", "s0", "s0-var"], [3]),
        (["contigName", "contigLen", "totalAvgDepth", "s0", "s0-var", "s1", "s1-var"], [3, 5]),
        # No variance columns at all: every column from the fourth on is a mean. The function is
        # written not to care, which is the point of testing the name rather than the stride.
        (["contigName", "contigLen", "totalAvgDepth", "s0", "s1", "s2"], [3, 4, 5]),
        (["contigName", "contigLen", "totalAvgDepth"], []),
    ],
)
def test_depth_means_finds_the_mean_columns(header, expected):
    assert binning.depth_means(header) == expected


def test_a_sample_named_var_is_not_mistaken_for_a_variance_column():
    """``-var`` is a suffix, not a substring: a BAM called ``var.bam`` is a sample."""
    header = ["contigName", "contigLen", "totalAvgDepth", "var", "var-var", "variance", "x-var"]
    assert binning.depth_means(header) == [3, 5]


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


def test_two_binners_share_one_depth_pass(tmp_path, stubs):
    """The pileup reads every BAM end to end, so doing it twice over the same BAMs is the whole cost.

    It used to run once per binner because each binner built its own table. Both read the same
    numbers, so there was never a second table to build - only a second pass to pay for.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n")
    r1, r2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"
    r1.write_text("@r/1\nACGT\n+\nIIII\n")
    r2.write_text("@r/1\nACGT\n+\nIIII\n")
    out = tmp_path / "out"

    rc = binning.main(
        [
            "-a", str(asm),
            "-o", str(out),
            "--metabat2",
            "--maxbin2",
            "-t", "1",
            "-m", "4",
            "--config", str(cfg),
            str(r1), str(r2),
        ]
    )  # fmt: skip
    assert rc == 0

    # Only the pileup calls count; the stub records exactly those.
    invocations = [ln for ln in stubs.read_text().splitlines() if ln.strip()]
    # One BAM here, so one table means one pileup. Two would mean each binner built its own.
    assert len(invocations) == 1, "one table, built once: %r" % invocations
    assert "--noIntraDepthVariance" not in invocations[0]
    assert (out / "work_files" / binning.DEPTH_TABLE).exists()
    # And MaxBin2 was handed depths, not the zeros the old flag produced. The sample is named for
    # its BAM, which is what jgi puts in the header.
    assert (out / "work_files" / "mb2_s.txt").read_text() == "c1\t7\n"


def test_a_given_depth_table_replaces_the_alignment_and_the_jgi_pass(tmp_path, stubs):
    """--depth is the whole coverage stage: no aligner, no jgi, and the table used as given.

    This is what makes a cohort's coverage cost N mappings rather than N per assembly: the cohort
    stage maps each sample once against the concatenated catalogue and hands each sample's binning
    its own slice of the result.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n>c2\n" + "ACGT" * 500 + "\n")
    table = tmp_path / "cohort_depth.txt"
    table.write_text(DEPTH)
    out = tmp_path / "out"

    rc = binning.main(
        [
            "-a", str(asm),
            "-o", str(out),
            "--metabat2",
            "--maxbin2",
            "--depth", str(table),
            "-t", "1",
            "-m", "4",
            "--config", str(cfg),
        ]
    )  # fmt: skip
    assert rc == 0

    # jgi was never invoked at all, so its call log does not even exist.
    assert not stubs.exists() or not [ln for ln in stubs.read_text().splitlines() if ln.strip()]
    assert not (out / "work_files" / binning.DEPTH_TABLE).exists(), "nothing was built"
    assert not list((out / "work_files").glob("*.bam")), "nothing was aligned"
    # Both samples in the given table reached MaxBin2, with their own means.
    assert (out / "work_files" / "mb2_sample_0.txt").read_text() == "c1\t3.5\nc2\t0.25\n"
    assert (out / "work_files" / "mb2_sample_1.txt").read_text() == "c1\t5.0\nc2\t1.0\n"


def test_a_depth_table_for_a_different_assembly_is_refused(tmp_path):
    """The failure this catches: a well-formed table whose contigs are someone else's.

    A cohort maps every sample against one concatenated catalogue, so a wrong slice of that table
    is exactly as well-formed as the right one. Without the check metaBAT2 is handed a table that
    matches nothing and bins on composition alone, which is a different algorithm and not an error.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">sampleB_c1\n" + "ACGT" * 500 + "\n")
    table = tmp_path / "depth.txt"
    table.write_text(DEPTH)  # names c1 and c2, which this assembly does not contain

    with pytest.raises(SystemExit):
        binning.main(
            [
                "-a", str(asm), "-o", str(tmp_path / "out"),
                "--metabat2", "--depth", str(table),
                "-t", "1", "-m", "4", "--config", str(cfg),
            ]
        )  # fmt: skip


@pytest.mark.parametrize("binner", binning.BAM_BINNERS)
def test_depth_is_refused_for_the_binners_that_read_the_alignments(tmp_path, binner):
    """CONCOCT builds its own coverage table and the three learned binners read the BAM directory.

    Accepting --depth for them would align the reads anyway and ignore the table, which is the
    opposite of what was asked for and costs the whole alignment pass to find out.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n")
    table = tmp_path / "depth.txt"
    table.write_text(DEPTH)

    with pytest.raises(SystemExit):
        binning.main(
            [
                "-a", str(asm), "-o", str(tmp_path / "out"),
                "--" + binner, "--depth", str(table),
                "-t", "1", "-m", "4", "--config", str(cfg),
            ]
        )  # fmt: skip


def test_without_a_depth_table_reads_are_still_required(tmp_path):
    """The positional reads became optional, so their absence has to be caught here instead.

    argparse used to refuse this with "the following arguments are required: reads"; nargs="*"
    accepts it, and an empty read list would otherwise align nothing and hand metaBAT2 an empty
    depth table.
    """
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n")
    with pytest.raises(SystemExit):
        binning.main(
            [
                "-a", str(asm), "-o", str(tmp_path / "out"),
                "--metabat2", "-t", "1", "-m", "4", "--config", str(cfg),
            ]
        )  # fmt: skip


def test_the_contigs_of_a_depth_table_are_read_from_it(tmp_path):
    table = tmp_path / "d.txt"
    table.write_text(DEPTH)
    assert binning.depth_table_contigs(str(table)) == ["c1", "c2"]


def test_a_table_that_is_not_a_depth_table_is_rejected(tmp_path):
    """A TSV of three plausible columns is what a hand-rolled coverage file looks like."""
    table = tmp_path / "d.txt"
    table.write_text("contig\tlength\tcoverage\nc1\t5000\t3.5\n")
    with pytest.raises(ValueError, match="is not a contig depth table"):
        binning.depth_table_contigs(str(table))


def test_a_table_with_no_sample_column_is_rejected(tmp_path):
    """metaBAT2 accepts it and bins on composition alone, which is not what --depth asked for."""
    table = tmp_path / "d.txt"
    table.write_text("contigName\tcontigLen\ttotalAvgDepth\nc1\t5000\t3.5\n")
    with pytest.raises(ValueError, match="no per-sample depth column"):
        binning.depth_table_contigs(str(table))


def test_a_table_covering_more_than_this_assembly_is_accepted(tmp_path):
    """A cohort table sliced per sample is the normal case, not a problem: extra rows are ignored."""
    asm = tmp_path / "asm.fa"
    asm.write_text(">c1\n" + "ACGT" * 500 + "\n")
    table = tmp_path / "d.txt"
    table.write_text(DEPTH)  # c1 and c2; this assembly has only c1
    assert binning.check_depth_table(str(table), str(asm)) == 2


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
