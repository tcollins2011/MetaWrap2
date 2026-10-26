"""A golden end-to-end test of the whole pipeline's plumbing, with stub tools.

It checks what unit tests cannot: that each module writes outputs where the next one looks for
them, that suffix renames are applied everywhere, and that manifests record paths the run
produced.

It does not test the science. The wrapped tools are small deterministic stubs on PATH, and the
test checks MetaWrap2's own contract:

* the **output tree** of each module: every file and directory it creates, by name;
* the **``.stats.tsv`` contract**: exact header, column order, row order and cell formatting;
* that each module can **consume the previous module's output directory** unchanged.

Golden values are literal in this file, so a change to one shows in review as a diff of the
expected output.

Covers one binning->refinement->consumer chain plus the .stats.tsv golden.
"""

from __future__ import annotations

import os
import stat
from typing import Dict, List

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2.constants import BIN_EXTENSION, STATS_SUFFIX
from metawrap2.modules import bin_refinement, binning

#: Contigs are 300 kb because bin_refinement discards bins outside 50 kb - 20 Mb, and
#: Binning_refiner discards refined bins under 0.5 Mb. Smaller contigs would let the test pass by
#: producing nothing.
CONTIG_LEN = 300_000

#: Four contigs, with no assembler-specific naming, since nothing may parse a length or bin number
#: from a contig header.
CONTIGS = (1, 2, 3, 4)


def contig(i: int) -> str:
    return ">c%d\n%s\n" % (i, "ACGT" * (CONTIG_LEN // 4))


ASSEMBLY = "".join(contig(i) for i in CONTIGS)

# --- stub tools ---------------------------------------------------------------------------
#
# Each stub is the smallest program that produces what the real tool produces, keyed off the
# flags MetaWrap2 passes, so a COMMANDS template change that moves an output path makes the test
# fail.

#: Preamble shared by the binner stubs, so the contigs they emit are identical to ASSEMBLY.
_PREAMBLE = (
    "import os, sys\n"
    "CONTIG_LEN = %d\n"
    "def contig(i):\n"
    "    return '>c%%d\\n%%s\\n' %% (i, 'ACGT' * (CONTIG_LEN // 4))\n" % CONTIG_LEN
)

STUBS: Dict[str, str] = {
    # binning's default aligner. "-d" builds the index; anything else maps, writing SAM to stdout
    # for the samtools sort on the other side of binning's pipe.
    "minimap2": (
        "import sys\n"
        "if '-d' in sys.argv:\n"
        "    open(sys.argv[sys.argv.index('-d') + 1], 'w').write('MMI')\n"
        "else:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
        "    for i in (1, 2, 3, 4):\n"
        "        sys.stdout.write('r%d\\t0\\tc%d\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'\n"
        "                         % (i, i))\n"
    ),
    # reassemble_bins still aligns with bwa, and binning offers it behind --aligner.
    "bwa": (
        "import sys\n"
        "if sys.argv[1] == 'mem':\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
        "    for i in (1, 2, 3, 4):\n"
        "        sys.stdout.write('r%d\\t0\\tc%d\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'\n"
        "                         % (i, i))\n"
    ),
    # Drains stdin like the real sort: binning pipes the aligner into it, and a sort that exits
    # early SIGPIPEs the aligner, which pipefail reports as a failed run.
    # samtools also supplies the pileup for the depth table, which metawrap2.depth builds (there is
    # no jgi_summarize_bam_contig_depths call), so the stub must emit alignments.
    "samtools": _PREAMBLE
    + (
        "a = sys.argv\n"
        "CONTIGS = ['c%d' % i for i in (1, 2, 3, 4)]\n"
        "if 'sort' in a:\n"
        "    sys.stdin.read()\n"
        "    if '-o' in a:\n"
        "        open(a[a.index('-o') + 1], 'w').write('BAM')\n"
        "elif 'view' in a and '-H' in a:\n"
        "    sys.stdout.write('@HD\\tVN:1.6\\n')\n"
        "    for c in CONTIGS:\n"
        "        sys.stdout.write('@SQ\\tSN:%s\\tLN:%d\\n' % (c, CONTIG_LEN))\n"
        "elif 'view' in a:\n"
        "    sys.stdout.write('ALIGNMENTS\\n')\n"
        "elif 'depth' in a:\n"
        "    sys.stdin.read()\n"
        "    for c in CONTIGS:\n"
        "        for pos in range(1, CONTIG_LEN + 1):\n"
        "            sys.stdout.write('%s\\t%d\\t5\\n' % (c, pos))\n"
    ),
    # The two binners agree that c1+c2 belong together and disagree about c3/c4, giving
    # bin_refinement one clear bin and one to resolve, so consolidation is exercised.
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
    """Deterministic scores, so the .stats.tsv contract can be checked literally.

    Scores come from the number of contigs in the bin, not a hash, since this test is about file
    format and wiring; the hash-based stand-in is in the refinement golden test.
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
    """Every path under *root*, relative and sorted, so the output layout can be compared."""
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

    # Every binner's output is renamed to one convention: metaBAT2's "bin.1.fa" and MaxBin2's
    # "bin.001.fasta" both become bin_001.fasta, counted from 1 and zero-padded so filename order
    # matches bin order.
    assert sorted(os.listdir(os.path.join(binned, "maxbin2_bins"))) == [
        "bin_001.fasta",
        "bin_002.fasta",
        "bin_003.fasta",
    ]


# --- the .stats.tsv contract --------------------------------------------------------------


def test_the_stats_file_is_byte_for_byte_what_it_has_always_been(tmp_path, fake_checkm):
    """The golden value. A change here changes a file users parse."""
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


# --- binning -> bin_refinement ------------------------------------------------------------


def test_refinement_consumes_binning_output_and_names_its_results(
    tmp_path, binned, config, fake_checkm
):
    """Refinement finds the bins where binning wrote them."""
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
            # GUNC needs a 6.7 GB reference database and is not part of the join under test.
            "--skip-gunc",
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


# --- refined bins -> annotate_bins --------------------------------------------------------
#
# annotate_bins consumes a bin directory, so changes to bin naming or layout show up here. prokka
# is stubbed to write the .gff and .faa the module looks for.

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
                "--skip-gunc",  # see above: no 6.7 GB reference database in a unit test
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


def test_annotation_covers_every_bin(tmp_path, refined, stubs, config):
    """One annotation per bin, in each of the three documented result folders.

    Checked in the reorganised output folders the module documents, not in prokka's working
    directory, where each bin also has a full copy.
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
