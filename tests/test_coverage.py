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


# ─── labels ───────────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "spec,label,path",
    [
        ("s0=/data/asm.fa", "s0", "/data/asm.fa"),
        (" s0 = /data/asm.fa ", "s0", "/data/asm.fa"),
        ("/data/sample_7.fasta", "sample_7", "/data/sample_7.fasta"),
        ("/data/sample_7.fa.gz", "sample_7", "/data/sample_7.fa.gz"),
        ("asm.fa", "asm", "asm.fa"),
    ],
)
def test_an_assembly_is_labelled_explicitly_or_by_its_filename(spec, label, path):
    assert coverage.parse_assembly(spec) == coverage.Assembly(label, path)


@pytest.mark.parametrize("spec", ["=asm.fa", "s0=", "="])
def test_a_half_written_label_is_an_error(spec):
    with pytest.raises(ValueError, match="both sides"):
        coverage.parse_assembly(spec)


def test_two_assemblies_cannot_share_a_label():
    """The cohort layout that makes this likely: one filename, repeated per sample directory."""
    with pytest.raises(ValueError, match="both labelled"):
        coverage.parse_assemblies(
            ["sample_0/final_assembly.fasta", "sample_1/final_assembly.fasta"]
        )


def test_a_label_containing_the_separator_is_refused():
    with pytest.raises(ValueError, match="contains"):
        coverage.parse_assemblies(["a~b=asm.fa"])


def test_the_same_labels_are_fine_under_a_separator_they_do_not_use():
    got = coverage.parse_assemblies(["a~b=one.fa", "c=two.fa"], separator="|")
    assert [a.label for a in got] == ["a~b", "c"]


# ─── the catalogue ────────────────────────────────────────────────────────────────────────


def _asm(tmp_path, name, contigs):
    path = tmp_path / name
    path.write_text("".join(">%s\n%s\n" % (c, "ACGT" * 250) for c in contigs))
    return str(path)


def test_the_catalogue_prefixes_every_contig_with_its_label(tmp_path):
    one = _asm(tmp_path, "one.fa", ["c1", "c2"])
    two = _asm(tmp_path, "two.fa", ["c1"])
    out = tmp_path / "catalogue.fa"

    counts = coverage.write_catalogue(
        [coverage.Assembly("a", one), coverage.Assembly("b", two)], str(out)
    )

    assert counts == {"a": 2, "b": 1}
    names = [ln[1:] for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert names == ["a~c1", "a~c2", "b~c1"]
    # The sequences are untouched: this is a renaming, not a filter.
    assert out.read_text().count("ACGT" * 250) == 3


def test_a_single_assembly_is_still_prefixed(tmp_path):
    """One code path for one assembly and for twenty, so the N=1 case is the tested one too."""
    one = _asm(tmp_path, "one.fa", ["c1"])
    out = tmp_path / "catalogue.fa"
    coverage.write_catalogue([coverage.Assembly("only", one)], str(out))
    assert [ln[1:] for ln in out.read_text().splitlines() if ln.startswith(">")] == ["only~c1"]


def test_a_contig_name_containing_the_separator_stops_the_run(tmp_path):
    """Refused at the catalogue, where it can be reported, not at the split, where it cannot."""
    one = _asm(tmp_path, "one.fa", ["c1", "weird~name"])
    with pytest.raises(ValueError, match="--separator"):
        coverage.write_catalogue([coverage.Assembly("a", one)], str(tmp_path / "cat.fa"))


def test_another_separator_makes_those_names_usable(tmp_path):
    one = _asm(tmp_path, "one.fa", ["weird~name"])
    out = tmp_path / "cat.fa"
    coverage.write_catalogue([coverage.Assembly("a", one)], str(out), separator="|")
    assert [ln[1:] for ln in out.read_text().splitlines() if ln.startswith(">")] == ["a|weird~name"]


def test_an_empty_assembly_is_an_error(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    with pytest.raises(ValueError, match="no contigs"):
        coverage.write_catalogue([coverage.Assembly("a", str(empty))], str(tmp_path / "cat.fa"))


def test_the_catalogue_keeps_only_the_first_header_token(tmp_path):
    """The aligner reports the first token, so that is the name the depth table will carry."""
    one = tmp_path / "one.fa"
    one.write_text(">c1 length=5000 cov=3.5\nACGT\n")
    out = tmp_path / "cat.fa"
    coverage.write_catalogue([coverage.Assembly("a", str(one))], str(out))
    assert [ln[1:] for ln in out.read_text().splitlines() if ln.startswith(">")] == ["a~c1"]


# ─── the split ────────────────────────────────────────────────────────────────────────────


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


def test_every_slice_carries_every_sample_column(tmp_path):
    """The point of the exercise: sample b's contigs get sample a's coverage too."""
    table = tmp_path / "t.txt"
    table.write_text(CATALOGUE_DEPTH)
    paths = coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a", "b"])
    for path in paths.values():
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        assert coverage.depth_means(header) == [3, 5]


def test_a_label_with_no_contigs_in_the_catalogue_is_an_error(tmp_path):
    """An empty slice is not an empty file to hand a binner: the catalogue is wrong."""
    table = tmp_path / "t.txt"
    table.write_text(CATALOGUE_DEPTH)
    with pytest.raises(ValueError, match="no contig in"):
        coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a", "missing"])


def test_contigs_of_an_unknown_label_are_skipped_not_misfiled(tmp_path):
    table = tmp_path / "t.txt"
    table.write_text(CATALOGUE_DEPTH)
    paths = coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a"])
    written = pathlib.Path(paths["a"]).read_text()
    assert "b~c1" not in written
    assert "c1\t4000" not in written


def test_a_table_that_is_not_a_depth_table_is_refused(tmp_path):
    table = tmp_path / "t.txt"
    table.write_text("contig\tdepth\na~c1\t3.5\n")
    with pytest.raises(ValueError, match="is not a contig depth table"):
        coverage.split_depth_table(str(table), str(tmp_path / "d"), ["a"])


@pytest.mark.parametrize(
    "name,expected",
    [
        ("a~c1", ("a", "c1")),
        # Everything after the first separator is the contig: a name is only in the catalogue if it
        # had no separator of its own, so the first one is always the one this module wrote.
        ("a~c1~x", ("a", "c1~x")),
        ("unprefixed", ("", "unprefixed")),
    ],
)
def test_a_catalogue_name_splits_at_the_first_separator(name, expected):
    assert coverage.split_catalogue_name(name) == expected


def test_prefixing_and_splitting_are_inverses():
    for label, contig in (("s0", "k141_1"), ("a-b", "NODE_1_length_5_cov_3"), ("x", "c~1")):
        assert coverage.split_catalogue_name(coverage.catalogue_name(label, contig)) == (
            label,
            contig,
        )


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


def test_keeping_the_bams_and_the_catalogue(tmp_path, stub_path):
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    one = _asm(tmp_path, "one.fa", ["c1"])
    s0 = _reads(tmp_path, "s0")
    out = tmp_path / "cov"

    rc = coverage.main(
        [
            "-a", "a=" + one,
            "-o", str(out),
            "-t", "1",
            "--keep-bams",
            "--keep-catalogue",
            "--config", str(cfg),
            s0[0], s0[1],
        ]
    )  # fmt: skip
    assert rc == 0
    assert [p.name for p in (out / "work_files").glob("*.bam")] == ["s0.bam"]
    assert (out / coverage.CATALOGUE).exists()


def test_reads_are_required(tmp_path):
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    one = _asm(tmp_path, "one.fa", ["c1"])
    with pytest.raises(SystemExit):
        coverage.main(["-a", one, "-o", str(tmp_path / "cov"), "--config", str(cfg)])


# ─── the k-mer route (--coverage fairy) ───────────────────────────────────────────────────

#: fairy, standing in for both of its subcommands. The sketch pass writes one file per -S name so
#: the module's own naming is what the coverage pass is handed; the coverage pass reads the catalogue
#: it was pointed at, for the same reason the jgi stub does, and emits the MetaBAT2 layout fairy
#: really emits - which is what lets the slices, the split and `binning --depth` stay route-agnostic.
_FAIRY_STUB = """#!/usr/bin/env python3
import os, sys
argv = sys.argv[1:]
if argv[0] == "sketch":
    names = []
    i = argv.index("-S") + 1
    while i < len(argv) and not argv[i].startswith("-"):
        names.append(argv[i])
        i += 1
    out = argv[argv.index("-d") + 1]
    kind = "single" if "-r" in argv else "paired"
    os.makedirs(out, exist_ok=True)
    for n in names:
        open(os.path.join(out, "%s.%s.bcsp" % (n, kind)), "w").write("SKETCH")
    sys.exit(0)
sketches = [a for a in argv if a.endswith(".bcsp")]
names = [os.path.basename(s).split(".")[0] for s in sketches]
cat = [a for a in argv if a.endswith(".fa") or a.endswith(".fasta")][0]
contigs = [ln[1:].split()[0] for ln in open(cat) if ln.startswith(">")]
cols = ["contigName", "contigLen", "totalAvgDepth"]
for n in names:
    cols += [n, n + "-var"]
rows = [cols] + [[c, "5000", "8.5"] + ["3.5", "91.1"] * len(names) for c in contigs]
with open(argv[argv.index("-o") + 1], "w") as fh:
    fh.write("".join("\\t".join(r) + "\\n" for r in rows))
"""


@pytest.fixture
def fairy_path(tmp_path, monkeypatch):
    binp = tmp_path / "fairybin"
    binp.mkdir()
    fairy = binp / "fairy"
    fairy.write_text(_FAIRY_STUB)
    fairy.chmod(fairy.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)
    monkeypatch.setenv("PATH", str(binp) + os.pathsep + os.environ["PATH"])
    return binp


def _fairy_run(tmp_path, extra=()):
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
            "--coverage", "fairy",
            "--config", str(cfg),
            *extra,
            s0[0], s0[1], s1[0], s1[1],
        ]
    )  # fmt: skip
    return rc, out


def test_the_fairy_route_produces_the_same_slices_without_aligning(tmp_path, fairy_path):
    # No aligner, no samtools and no jgi on PATH in this fixture: if the module reached for one the
    # run would fail, which is the assertion. The point of the route is that it does not.
    rc, out = _fairy_run(tmp_path)
    assert rc == 0
    slices = sorted(p.name for p in (out / coverage.DEPTH_DIR).glob("*.txt"))
    assert slices == ["a.txt", "b.txt"]
    header = (out / coverage.DEPTH_DIR / "a.txt").read_text().splitlines()[0].split("\t")
    # The sample names come from the read filenames via -S, where the alignment route gets them
    # from the BAM names; downstream nothing can tell the two apart, which is the contract.
    assert [header[i] for i in coverage.depth_means(header)] == ["s0", "s1"]
    assert (out / coverage.DEPTH_DIR / "a.txt").read_text().splitlines()[1].startswith("c1\t")
    assert len((out / coverage.DEPTH_DIR / "b.txt").read_text().splitlines()) == 2


def test_the_fairy_route_sketches_each_sample_once_and_keeps_the_sketches(tmp_path, fairy_path):
    _rc, out = _fairy_run(tmp_path)
    sketches = sorted(p.name for p in (out / coverage.WORK_DIR / coverage.SKETCH_DIR).iterdir())
    assert sketches == ["s0.paired.bcsp", "s1.paired.bcsp"]


def test_the_fairy_route_sketches_interleaved_reads_as_single_end(tmp_path, fairy_path):
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    one = _asm(tmp_path, "one.fa", ["c1"])
    reads = tmp_path / "s0.fastq"
    reads.write_text("@r/1\nACGT\n+\nIIII\n@r/2\nACGT\n+\nIIII\n")
    out = tmp_path / "cov"
    rc = coverage.main(
        [
            "-a", "a=" + one,
            "-o", str(out),
            "-t", "1",
            "--coverage", "fairy",
            "--interleaved",
            "--config", str(cfg),
            str(reads),
        ]
    )  # fmt: skip
    assert rc == 0
    # Which mate a k-mer came from is not part of a sketch, so the pairing is not needed here - but
    # fairy names the file for the mode it was given, and the coverage call takes paths.
    sketches = [p.name for p in (out / coverage.WORK_DIR / coverage.SKETCH_DIR).iterdir()]
    assert sketches == ["s0.single.bcsp"]


def test_the_route_is_recorded_so_a_resumed_run_does_not_mix_the_two(tmp_path, fairy_path):
    _rc, out = _fairy_run(tmp_path)
    markers = sorted(p.name for p in pathlib.Path(out, ".metawrap2", "steps").iterdir())
    # Stage names differ by route, so resuming with --coverage align rebuilds the table rather than
    # inheriting one fairy wrote (and vice versa).
    assert markers == ["catalogue.done", "fairy_depth_table.done", "fairy_sketch.done"]


def test_the_alignment_route_is_still_the_default(tmp_path):
    assert coverage.DEFAULT_ROUTE == "align"
    # Not an opinion about which is better: the default is the one the binners' accuracy was
    # measured with, and B5's coverage arm is what would move it.
    assert coverage._parse_args(["-a", "x.fa", "-o", "o"]).coverage == "align"
