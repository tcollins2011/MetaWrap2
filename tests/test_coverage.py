"""The cohort coverage stage: one catalogue, N alignments, one table, N slices.

What these tests are really about is the demultiplexing. Every piece of the module exists to make a
contig's sample of origin recoverable from the catalogue-wide depth table, and the failure mode if
that is wrong is not a crash: a slice of plausible depths attributed to the wrong assembly is a
well-formed depth table, and `binning --depth` would bin on it happily. So the label prefix, the
separator check, the restored contig names and the refusal to write an empty slice are each pinned
here.
"""

from __future__ import annotations

import os
import pathlib
import stat

import pytest

from metawrap2.modules import coverage

#: A catalogue-wide table in jgi's layout, two samples, contigs from two labelled assemblies.
CATALOGUE_DEPTH = (
    "contigName\tcontigLen\ttotalAvgDepth\ts0\ts0-var\ts1\ts1-var\n"
    "a~c1\t5000\t8.5\t3.5\t91.1\t5.0\t92.2\n"
    "a~c2\t2000\t1.25\t0.25\t93.3\t1.0\t94.4\n"
    "b~c1\t4000\t2.0\t1.5\t95.5\t0.5\t96.6\n"
)


# ─── the split ────────────────────────────────────────────────────────────────────────────


def _asm(tmp_path, name, contigs):
    path = tmp_path / name
    path.write_text("".join(">%s\n%s\n" % (c, "ACGT" * 250) for c in contigs))
    return str(path)


def test_splitting_restores_the_original_contig_names(tmp_path):
    """Each slice describes the assembly it will be binned against, not the catalogue.

    Prefixed names against an unprefixed assembly is exactly what `binning`'s own depth check
    refuses, so leaving them in would turn a correct run into a hard stop one command later.
    """
    table = tmp_path / "catalogue_depth.txt"
    table.write_text(CATALOGUE_DEPTH)

    paths = coverage.split_depth_table(str(table), str(tmp_path / "depth"), ["a", "b"])

    assert sorted(paths) == ["a", "b"]
    assert pathlib.Path(paths["a"]).read_text() == (
        "contigName\tcontigLen\ttotalAvgDepth\ts0\ts0-var\ts1\ts1-var\n"
        "c1\t5000\t8.5\t3.5\t91.1\t5.0\t92.2\n"
        "c2\t2000\t1.25\t0.25\t93.3\t1.0\t94.4\n"
    )
    assert pathlib.Path(paths["b"]).read_text() == (
        "contigName\tcontigLen\ttotalAvgDepth\ts0\ts0-var\ts1\ts1-var\n"
        "c1\t4000\t2.0\t1.5\t95.5\t0.5\t96.6\n"
    )
    # An empty slice is not an empty file to hand a binner: the catalogue is wrong.
    with pytest.raises(ValueError, match="no contig in"):
        coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a", "missing"])


# ─── through main() ───────────────────────────────────────────────────────────────────────

#: jgi, emitting a two-sample table over whatever contigs the catalogue it was pointed at holds.
#: Reading the catalogue is the point: a stub with a fixed contig list could not catch a module that
#: never prefixed anything.
_STUBS = {
    "minimap2": (
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d')+1], 'w').write('IDX')\n"
        "else:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
    ),
    # The depth table is computed in-process now (metawrap2.depth) rather than shelled out to
    # jgi_summarize_bam_contig_depths, so samtools has to answer what the header names and what the
    # per-base depth is. The contigs come from the catalogue this module wrote, which is where the
    # old jgi stub read them from too - a header that disagreed with the catalogue would be a BAM
    # that was never aligned against it.
    "samtools": (
        "import os, sys\n"
        "a = sys.argv\n"
        "def contigs(bam):\n"
        "    cat = os.path.join(os.path.dirname(os.path.dirname(bam)), 'catalogue.fa')\n"
        "    return [ln[1:].split()[0] for ln in open(cat) if ln.startswith('>')]\n"
        "if 'sort' in a:\n"
        "    sys.stdin.read()\n"
        "    if '-o' in a:\n"
        "        open(a[a.index('-o')+1], 'w').write('BAM')\n"
        "elif 'view' in a and '-H' in a:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
        "    for c in contigs(a[-1]):\n"
        "        sys.stdout.write('@SQ\\tSN:%s\\tLN:5000\\n' % c)\n"
        # view -u names the contigs it is streaming and depth reads them back off stdin: the pileup
        # is given "-" rather than a path, so this is the only way the two halves of the pipe agree
        # on which contigs exist - which is also true of the real pair.
        "elif 'view' in a:\n"
        "    for c in contigs(a[-1]):\n"
        "        sys.stdout.write('%s\\n' % c)\n"
        "elif 'depth' in a:\n"
        "    for line in sys.stdin:\n"
        "        c = line.strip()\n"
        "        if not c:\n"
        "            continue\n"
        "        for pos in range(1, 5001):\n"
        "            sys.stdout.write('%s\\t%d\\t7\\n' % (c, pos))\n"
    ),
}


@pytest.fixture
def stub_path(tmp_path, monkeypatch):
    binp = tmp_path / "stubbin"
    binp.mkdir()
    for name, body in _STUBS.items():
        f = binp / name
        f.write_text("#!/usr/bin/env python3\n" + body)
        f.chmod(f.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)
    monkeypatch.setenv("PATH", str(binp) + os.pathsep + os.environ["PATH"])
    return binp


def _reads(tmp_path, sample):
    r1 = tmp_path / ("%s_1.fastq" % sample)
    r2 = tmp_path / ("%s_2.fastq" % sample)
    r1.write_text("@r/1\nACGT\n+\nIIII\n")
    r2.write_text("@r/2\nACGT\n+\nIIII\n")
    return str(r1), str(r2)


def test_a_cohort_run_maps_once_per_sample_and_writes_one_slice_per_assembly(tmp_path, stub_path):
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    one = _asm(tmp_path, "one.fa", ["c1", "c2"])
    two = _asm(tmp_path, "two.fa", ["c1"])
    s0 = _reads(tmp_path, "s0")
    s1 = _reads(tmp_path, "s1")
    out = tmp_path / "cov"

    rc = coverage.main(
        [
            "-a", "a=" + one,
            "-a", "b=" + two,
            "-o", str(out),
            "-t", "1",
            "--config", str(cfg),
            s0[0], s0[1], s1[0], s1[1],
        ]
    )  # fmt: skip
    assert rc == 0

    # Two samples, two alignments - not two samples times two assemblies - and the BAMs the
    # alignment left are swept once the one table is built.
    assert not list((out / "work_files").glob("*.bam"))
    slices = sorted(p.name for p in (out / coverage.DEPTH_DIR).glob("*.txt"))
    assert slices == ["a.txt", "b.txt"]
    assert (out / coverage.DEPTH_DIR / "a.txt").read_text().splitlines()[1].startswith("c1\t")
    assert len((out / coverage.DEPTH_DIR / "b.txt").read_text().splitlines()) == 2
    # Both samples' columns are in both slices, which is the whole reason for the catalogue.
    header = (out / coverage.DEPTH_DIR / "b.txt").read_text().splitlines()[0].split("\t")
    assert [header[i] for i in coverage.depth_means(header)] == ["s0", "s1"]
    # And the catalogue itself is gone unless asked for.
    assert not (out / coverage.CATALOGUE).exists()
    assert (out / coverage.ASSEMBLIES_TABLE).read_text().splitlines()[1].startswith("a\t")
