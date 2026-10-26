"""The cohort coverage stage: one catalogue, N alignments, one table, N slices.

These tests cover demultiplexing: recovering each contig's sample of origin from the catalogue-wide
depth table. A slice attributed to the wrong assembly is still a valid depth table that `binning
--depth` would accept, so the label prefix, separator check, restored contig names and refusal to
write an empty slice are each tested.
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
    """Each slice uses the contig names of its own assembly, without the catalogue prefix.

    `binning`'s depth check rejects prefixed names against an unprefixed assembly.
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
    # An empty slice means the catalogue is wrong, so no file is written.
    with pytest.raises(ValueError, match="no contig in"):
        coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a", "missing"])


# ─── through main() ───────────────────────────────────────────────────────────────────────

#: Stub emitting a two-sample table over the contigs in the catalogue. Reading the catalogue lets
#: the test catch a module that never adds prefixes.
_STUBS = {
    "minimap2": (
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d')+1], 'w').write('IDX')\n"
        "else:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
    ),
    # The depth table is computed in-process (metawrap2.depth), so the samtools stub supplies the
    # header and per-base depth. Contigs come from the catalogue this module wrote, matching a BAM
    # aligned against it.
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
        # view -u passes the contig names on stdout and depth reads them from stdin ("-"), as in
        # the real pipe.
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

    # Two samples give two alignments (not samples x assemblies), and the BAMs are removed once
    # the table is built.
    assert not list((out / "work_files").glob("*.bam"))
    slices = sorted(p.name for p in (out / coverage.DEPTH_DIR).glob("*.txt"))
    assert slices == ["a.txt", "b.txt"]
    assert (out / coverage.DEPTH_DIR / "a.txt").read_text().splitlines()[1].startswith("c1\t")
    assert len((out / coverage.DEPTH_DIR / "b.txt").read_text().splitlines()) == 2
    # Both samples' columns are in both slices.
    header = (out / coverage.DEPTH_DIR / "b.txt").read_text().splitlines()[0].split("\t")
    assert [header[i] for i in coverage.depth_means(header)] == ["s0", "s1"]
    # The catalogue is removed unless requested.
    assert not (out / coverage.CATALOGUE).exists()
    assert (out / coverage.ASSEMBLIES_TABLE).read_text().splitlines()[1].startswith("a\t")
