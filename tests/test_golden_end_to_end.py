"""A golden end-to-end test of the whole pipeline's plumbing, with stub tools.

This is the test that catches the mistakes unit tests structurally cannot see: a module that
writes its outputs where the *next* module does not look for them, a renamed suffix that only
half the codebase learned about, a manifest that records a path the run never produced.

It deliberately does not test the science. The wrapped tools are tiny deterministic stubs on
PATH, so what is being pinned is MetaWrap2's own contract:

* the **output tree** of each module - every file and directory it creates, by name;
* the **``.stats.tsv`` contract** - exact header, column order, row order, and cell formatting;
* that each module can **consume the previous module's output directory** unchanged;
* that the **manifest and provenance** describe what actually happened.

Golden values are literal in this file. When a change alters one, that shows up in review as a
diff of the expected output, which is the point - the failure is the notification.
"""

from __future__ import annotations

import json
import os
import stat
from typing import Dict, List

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2.constants import BIN_EXTENSION, STATS_SUFFIX
from metawrap2.modules import bin_refinement, binning

#: Contigs are 300 kb because bin_refinement discards bins outside 50 kb - 20 Mb, and
#: Binning_refiner discards refined bins under 0.5 Mb. Anything smaller makes this test pass by
#: producing nothing, which is the opposite of what it is for.
CONTIG_LEN = 300_000

#: Four contigs, named without any assembler-specific convention - nothing in the pipeline may
#: parse a length or a bin number out of a contig header.
CONTIGS = (1, 2, 3, 4)


def contig(i: int) -> str:
    return ">c%d\n%s\n" % (i, "ACGT" * (CONTIG_LEN // 4))


ASSEMBLY = "".join(contig(i) for i in CONTIGS)

# --- stub tools ---------------------------------------------------------------------------
#
# Each stub is the smallest program that produces what the real tool produces, keyed off the
# same command-line flags MetaWrap2 passes - so a change to a COMMANDS template that moves an
# output path makes the stub miss too, and the test fails.

#: Preamble shared by the binner stubs, so the contigs they emit are identical to ASSEMBLY.
_PREAMBLE = (
    "import os, sys\n"
    "CONTIG_LEN = %d\n"
    "def contig(i):\n"
    "    return '>c%%d\\n%%s\\n' %% (i, 'ACGT' * (CONTIG_LEN // 4))\n" % CONTIG_LEN
)

STUBS: Dict[str, str] = {
    "bwa": (
        "import sys\n"
        "if sys.argv[1] == 'mem':\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
        "    for i in (1, 2, 3, 4):\n"
        "        sys.stdout.write('r%d\\t0\\tc%d\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'\n"
        "                         % (i, i))\n"
    ),
    "samtools": (
        "import sys\n"
        "if 'sort' in sys.argv and '-o' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-o') + 1], 'w').write('BAM')\n"
    ),
    "jgi_summarize_bam_contig_depths": _PREAMBLE
    + (
        "p = sys.argv[sys.argv.index('--outputDepth') + 1]\n"
        "rows = ''.join('c%d\\t%d\\t5\\t5\\n' % (i, CONTIG_LEN) for i in (1, 2, 3, 4))\n"
        "open(p, 'w').write('contigName\\tcontigLen\\ttotalAvgDepth\\ts.bam\\n' + rows)\n"
    ),
    # The two binners agree that c1+c2 belong together and disagree about c3/c4. That gives
    # bin_refinement one bin worth keeping and one worth arguing over - a set where every binner
    # agrees, or none do, would not exercise consolidation at all.
    "metabat2": _PREAMBLE
    + (
        "pre = sys.argv[sys.argv.index('-o') + 1]\n"
        "os.makedirs(os.path.dirname(pre) or '.', exist_ok=True)\n"
        "open(pre + '.1.fa', 'w').write(contig(1) + contig(2))\n"
        "open(pre + '.2.fa', 'w').write(contig(3) + contig(4))\n"
    ),
    "run_MaxBin.pl": _PREAMBLE
    + (
        "pre = sys.argv[sys.argv.index('-out') + 1]\n"
        "os.makedirs(os.path.dirname(pre) or '.', exist_ok=True)\n"
        "open(pre + '.001.fasta', 'w').write(contig(1) + contig(2))\n"
        "open(pre + '.002.fasta', 'w').write(contig(3))\n"
        "open(pre + '.003.fasta', 'w').write(contig(4))\n"
    ),
}


@pytest.fixture
def stubs(tmp_path, monkeypatch):
    binp = tmp_path / "stubbin"
    binp.mkdir()
    for name, body in STUBS.items():
        path = binp / name
        path.write_text("#!/usr/bin/env python3\n" + body)
        path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)
    monkeypatch.setenv("PATH", str(binp) + os.pathsep + os.environ["PATH"])
    return binp


@pytest.fixture
def config(tmp_path):
    """Conda envs off, so the stubs on PATH are what runs."""
    cfg = tmp_path / "mw2.toml"
    cfg.write_text("[settings]\nuse_conda_envs = false\n")
    return str(cfg)


@pytest.fixture
def fake_checkm(monkeypatch):
    """Deterministic scores, so the .stats.tsv contract can be pinned exactly.

    Scores come from the number of contigs in the bin rather than a hash, because this test is
    about file format and wiring; the hash-based stand-in belongs in the refinement golden test.
    """

    def scores(bin_path):
        with open(bin_path) as fh:
            n = sum(1 for line in fh if line.startswith(">"))
        return 90.0 - 10.0 * (n - 1), 1.5 * (n - 1), 1000 * n, 1000 * n

    def run_checkm(bins_dir, **kwargs):
        rows: List[List[str]] = []
        for name in sorted(os.listdir(bins_dir)):
            if not name.endswith(BIN_EXTENSION):
                continue
            completeness, contamination, n50, size = scores(os.path.join(bins_dir, name))
            rows.append(
                [
                    os.path.splitext(name)[0],
                    "%.2f" % completeness,
                    "%.3f" % contamination,
                    "0.500",
                    "Bacteria",
                    str(n50),
                    str(size),
                ]
            )
        checkm_mod.write_stats(rows, bins_dir + STATS_SUFFIX)
        return bins_dir + STATS_SUFFIX

    monkeypatch.setattr(checkm_mod, "run_checkm", run_checkm)
    monkeypatch.setattr("metawrap2.modules.bin_refinement.run_checkm", run_checkm, raising=False)
    return run_checkm


def tree(root: str) -> List[str]:
    """Every path under *root*, relative and sorted - the output layout as a comparable value."""
    found = []
    for base, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if not d.startswith(".")]  # scratch is not part of the tree
        for name in sorted(dirs) + sorted(files):
            rel = os.path.relpath(os.path.join(base, name), root)
            if not any(part.startswith(".") for part in rel.split(os.sep)):
                found.append(rel)
    return sorted(found)


# --- binning ------------------------------------------------------------------------------


@pytest.fixture
def binned(tmp_path, stubs, config):
    """Run binning with two binners over four contigs. The rest of the file builds on this."""
    asm = tmp_path / "asm.fa"
    asm.write_text(ASSEMBLY)
    for suffix in ("_1", "_2"):
        (tmp_path / ("s%s.fastq" % suffix)).write_text(
            "".join("@r%d/1\nACGT\n+\nIIII\n" % i for i in CONTIGS)
        )
    out = tmp_path / "BINNING"
    rc = binning.main(
        [
            "-a",
            str(asm),
            "-o",
            str(out),
            "--metabat2",
            "--maxbin2",
            "-t",
            "1",
            "--config",
            config,
            str(tmp_path / "s_1.fastq"),
            str(tmp_path / "s_2.fastq"),
        ]
    )
    assert rc == 0, "binning failed"
    return out


def test_binning_produces_the_expected_bin_sets(binned):
    """Both binners' output directories exist under the names bin_refinement expects."""
    assert os.path.isdir(os.path.join(binned, "metabat2_bins"))
    assert os.path.isdir(os.path.join(binned, "maxbin2_bins"))
    assert sorted(os.listdir(os.path.join(binned, "metabat2_bins"))) == [
        "bin_001.fasta",
        "bin_002.fasta",
    ]

    # Every binner's output carries the same convention, whatever the tool called its files:
    # metaBAT2's "bin.1.fa" and MaxBin2's "bin.001.fasta" both become bin_001.fasta, counted
    # from 1 and zero-padded so filename order and bin order are the same thing.
    assert sorted(os.listdir(os.path.join(binned, "maxbin2_bins"))) == [
        "bin_001.fasta",
        "bin_002.fasta",
        "bin_003.fasta",
    ]


def test_binning_writes_provenance_and_a_log(binned):
    for name in ("run_config.json", "run_commands.txt", "provenance.txt", "metawrap2.log"):
        path = os.path.join(binned, name)
        assert os.path.isfile(path), "%s is missing" % name
        assert os.path.getsize(path) > 0, "%s is empty" % name


def test_binning_records_the_tool_versions_it_used(binned):
    with open(os.path.join(binned, "run_config.json")) as fh:
        record = json.load(fh)
    assert record["module"] == "binning"
    assert record["commands"], "no commands recorded"
    assert "tool_versions" in record


def test_binning_leaves_no_scratch_behind(binned):
    """Temporary space lives under the output directory and must not survive a clean run."""
    leftovers = [n for n in os.listdir(binned) if n.startswith(".metawrap2_tmp")]
    for name in leftovers:
        assert not os.listdir(os.path.join(binned, name)), "%s still holds files" % name


# --- the .stats.tsv contract --------------------------------------------------------------


def test_stats_header_is_exactly_the_documented_columns(tmp_path, fake_checkm):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin_001.fasta").write_text(">c1\nACGT\n")
    path = checkm_mod.run_checkm(str(bins))

    assert path.endswith(".stats.tsv"), "the suffix must say what the file is"
    with open(path) as fh:
        lines = fh.read().splitlines()
    assert lines[0].split("\t") == [
        "bin",
        "completeness",
        "contamination",
        "GC",
        "lineage",
        "N50",
        "size",
    ]


def test_stats_rows_are_sorted_by_completeness_descending(tmp_path, fake_checkm):
    """Downstream code and users both read the first row as the best bin."""
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin_001.fasta").write_text(">c1\nACGT\n>c2\nACGT\n>c3\nACGT\n")  # 3 contigs -> worst
    (bins / "bin_002.fasta").write_text(">c4\nACGT\n")  # 1 contig -> best
    with open(checkm_mod.run_checkm(str(bins))) as fh:
        rows = [line.split("\t") for line in fh.read().splitlines()[1:]]
    assert [r[0] for r in rows] == ["bin_002", "bin_001"]
    assert [float(r[1]) for r in rows] == [90.0, 70.0]


def test_the_stats_file_is_byte_for_byte_what_it_has_always_been(tmp_path, fake_checkm):
    """The golden value. A change here is a change to a file users parse."""
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin_001.fasta").write_text(">c1\n" + "ACGT" * 10 + "\n")
    (bins / "bin_002.fasta").write_text(">c2\nACGT\n>c3\nACGT\n")
    with open(checkm_mod.run_checkm(str(bins))) as fh:
        content = fh.read()
    assert content == (
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin_001\t90.00\t0.000\t0.500\tBacteria\t1000\t1000\n"
        "bin_002\t80.00\t1.500\t0.500\tBacteria\t2000\t2000\n"
    )


def test_stats_gains_a_binner_column_only_when_asked(tmp_path):
    """bin_refinement's summary names the source binner; plain CheckM output must not."""
    bins = tmp_path / "bins"
    bins.mkdir()
    rows = [["bin_001", "90.00", "0.000", "0.500", "Bacteria", "1000", "1000", "metabat2"]]
    checkm_mod.write_stats(rows, str(bins) + STATS_SUFFIX, binner="metabat2")
    with open(str(bins) + STATS_SUFFIX) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    assert header[-1] == "binner"


# --- binning -> bin_refinement ------------------------------------------------------------


def test_refinement_consumes_binning_output_and_names_its_results(
    tmp_path, binned, config, fake_checkm
):
    """The join between two modules: refinement must find bins where binning left them."""
    out = tmp_path / "BIN_REFINEMENT"
    rc = bin_refinement.main(
        [
            "-o",
            str(out),
            "-t",
            "1",
            "-c",
            "50",
            "-x",
            "10",
            "--config",
            config,
            "-A",
            os.path.join(binned, "metabat2_bins"),
            "-B",
            os.path.join(binned, "maxbin2_bins"),
        ]
    )
    assert rc == 0, "bin_refinement failed"

    produced = tree(str(out))
    assert "metawrap_50_10_bins" in produced
    assert "metawrap_50_10_bins.stats.tsv" in produced
    assert "metawrap_50_10_bins.contigs.tsv" in produced
    # No bare .stats / .contigs files may survive the rename anywhere in the tree.
    assert not [p for p in produced if p.endswith((".stats", ".contigs"))]


def test_refinement_stats_and_bins_agree(tmp_path, binned, config, fake_checkm):
    """Every bin on disk has a row, and every row has a bin - the invariant dereplication needs."""
    out = tmp_path / "BIN_REFINEMENT"
    assert (
        bin_refinement.main(
            [
                "-o",
                str(out),
                "-t",
                "1",
                "-c",
                "50",
                "-x",
                "10",
                "--config",
                config,
                "-A",
                os.path.join(binned, "metabat2_bins"),
                "-B",
                os.path.join(binned, "maxbin2_bins"),
            ]
        )
        == 0
    )
    bins_dir = os.path.join(out, "metawrap_50_10_bins")
    on_disk = sorted(
        os.path.splitext(n)[0] for n in os.listdir(bins_dir) if n.endswith(BIN_EXTENSION)
    )
    with open(os.path.join(out, "metawrap_50_10_bins" + STATS_SUFFIX)) as fh:
        in_stats = sorted(line.split("\t")[0] for line in fh.read().splitlines()[1:] if line)
    assert on_disk == in_stats


# --- refined bins -> quant_bins and annotate_bins -----------------------------------------
#
# These two are the other consumers of a bin directory, so they are where a change to how bins
# are named or laid out shows up. Their tools are stubbed the same way as the binners: salmon
# writes the quant.sf it really writes, prokka writes the .gff and .faa the module goes looking
# for. Anything that moves an expected path breaks the stub, which is the point.

QUANT_STUBS: Dict[str, str] = {
    "salmon": (
        "import os, sys\n"
        "if sys.argv[1] == 'index':\n"
        "    os.makedirs(sys.argv[sys.argv.index('-i') + 1], exist_ok=True)\n"
        "elif sys.argv[1] == 'quant':\n"
        "    out = sys.argv[sys.argv.index('-o') + 1]\n"
        "    os.makedirs(out, exist_ok=True)\n"
        "    rows = ''.join('c%d\\t1000\\t900\\t%d.5\\t%d\\n' % (i, i, i * 10)\n"
        "                   for i in (1, 2, 3, 4))\n"
        "    open(os.path.join(out, 'quant.sf'), 'w').write(\n"
        "        'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' + rows)\n"
    ),
}

PROKKA_STUB = {
    "prokka": (
        "import os, sys\n"
        "outdir = sys.argv[sys.argv.index('--outdir') + 1]\n"
        "prefix = sys.argv[sys.argv.index('--prefix') + 1]\n"
        "os.makedirs(outdir, exist_ok=True)\n"
        "open(os.path.join(outdir, prefix + '.gff'), 'w').write(\n"
        "    '##gff-version 3\\nc1\\tProdigal\\tCDS\\t1\\t99\\t.\\t+\\t0\\tID=%s_00001\\n' % prefix)\n"
        "open(os.path.join(outdir, prefix + '.faa'), 'w').write('>%s_00001\\nMKV\\n' % prefix)\n"
        "open(os.path.join(outdir, prefix + '.ffn'), 'w').write('>%s_00001\\nATGAAAGTT\\n' % prefix)\n"
    ),
}


def install_stubs(directory, stubs: Dict[str, str]) -> None:
    for name, body in stubs.items():
        path = directory / name
        path.write_text("#!/usr/bin/env python3\n" + body)
        path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IRWXU)


@pytest.fixture
def refined(tmp_path, binned, config, fake_checkm):
    """The refined bin set, as the later modules receive it."""
    out = tmp_path / "BIN_REFINEMENT"
    assert (
        bin_refinement.main(
            [
                "-o",
                str(out),
                "-t",
                "1",
                "-c",
                "50",
                "-x",
                "10",
                "--config",
                config,
                "-A",
                os.path.join(binned, "metabat2_bins"),
                "-B",
                os.path.join(binned, "maxbin2_bins"),
            ]
        )
        == 0
    )
    bins = out / "metawrap_50_10_bins"
    assert bins.is_dir() and list(bins.glob("*" + BIN_EXTENSION)), "refinement produced no bins"
    return bins


def test_quant_bins_consumes_the_refined_bins(tmp_path, refined, stubs, config):
    from metawrap2.modules import quant_bins

    install_stubs(stubs, QUANT_STUBS)
    asm = tmp_path / "asm.fa"
    asm.write_text(ASSEMBLY)
    out = tmp_path / "QUANT_BINS"
    rc = quant_bins.main(
        [
            "-b",
            str(refined),
            "-o",
            str(out),
            "-a",
            str(asm),
            "-t",
            "1",
            "--config",
            config,
            str(tmp_path / "s_1.fastq"),
            str(tmp_path / "s_2.fastq"),
        ]
    )
    assert rc == 0, "quant_bins failed"

    produced = tree(str(out))
    # The abundance table is the module's contract with the user, and its extension must say so.
    assert "bin_abundance_table.tsv" in produced
    assert any(p.startswith("quant_files") and p.endswith(".counts") for p in produced)
    assert not [p for p in produced if p.endswith((".stats", ".contigs"))]


def test_the_abundance_table_has_a_row_per_bin(tmp_path, refined, stubs, config):
    from metawrap2.modules import quant_bins

    install_stubs(stubs, QUANT_STUBS)
    asm = tmp_path / "asm.fa"
    asm.write_text(ASSEMBLY)
    out = tmp_path / "QUANT_BINS"
    assert (
        quant_bins.main(
            [
                "-b",
                str(refined),
                "-o",
                str(out),
                "-a",
                str(asm),
                "-t",
                "1",
                "--config",
                config,
                str(tmp_path / "s_1.fastq"),
                str(tmp_path / "s_2.fastq"),
            ]
        )
        == 0
    )
    with open(out / "bin_abundance_table.tsv") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    header, rows = lines[0], lines[1:]
    assert header.split("\t")[0] in ("Genomic bins", "bin", "Bin")
    on_disk = {p.stem for p in refined.glob("*" + BIN_EXTENSION)}
    named = {row.split("\t")[0] for row in rows}
    assert named == on_disk, "the table's bins must be exactly the bins it was given"


def test_annotate_bins_consumes_the_refined_bins(tmp_path, refined, stubs, config):
    from metawrap2.modules import annotate_bins

    install_stubs(stubs, PROKKA_STUB)
    out = tmp_path / "FUNCT_ANNOT"
    rc = annotate_bins.main(["-b", str(refined), "-o", str(out), "-t", "1", "--config", config])
    assert rc == 0, "annotate_bins failed"

    produced = tree(str(out))
    assert any(p.endswith(".gff") for p in produced), "no GFF was produced"
    assert any(p.endswith(".faa") for p in produced), "no protein FASTA was produced"


def test_annotation_covers_every_bin(tmp_path, refined, stubs, config):
    """One annotation per bin, in each of the three documented result folders.

    A bin silently skipped is the failure worth catching, and the place to see it is the
    reorganised output the module promises the user - not prokka's own working directory, where
    each bin also keeps a full copy.
    """
    from metawrap2.modules import annotate_bins

    install_stubs(stubs, PROKKA_STUB)
    out = tmp_path / "FUNCT_ANNOT"
    assert (
        annotate_bins.main(["-b", str(refined), "-o", str(out), "-t", "1", "--config", config]) == 0
    )
    bins = sorted(p.stem for p in refined.glob("*" + BIN_EXTENSION))
    for folder, suffix in (
        ("bin_funct_annotations", ".gff"),
        ("bin_translated_genes", ".faa"),
        ("bin_untranslated_genes", ".ffn"),
    ):
        directory = out / folder
        assert directory.is_dir(), "%s is missing" % folder
        found = sorted(p.name[: -len(suffix)] for p in directory.glob("*" + suffix))
        assert found == bins, "%s holds %s, expected %s" % (folder, found, bins)
