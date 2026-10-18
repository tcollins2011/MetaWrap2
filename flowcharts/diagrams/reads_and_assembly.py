"""READ_QC, ASSEMBLY, KRAKEN2 and NONPAREIL: from raw reads to an assembly and a community profile."""

from _grid import BLUE, GREEN, RED, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("reads_and_assembly", "Reads and assembly")
    d.box("raw", "Raw reads\n(paired, single or interleaved)", "data", x(0), y(0), w=W)
    d.box("fastqc", "FastQC report\non the raw reads", "step", x(1), y(0), w=W)
    d.box(
        "trim", "Trim adapters and quality\nfastp (default) or Trim Galore", "step", x(2), y(0), w=W
    )
    d.box("host", "Remove host reads\nbmtagger, -x host index", "step", x(3), y(0), w=W)
    d.box(
        "clean",
        "Clean reads\nfinal_pure_reads_1/2.fastq\n(FastQC report after)",
        "data",
        x(4),
        y(0),
        w=W,
    )
    d.panel_around("p_qc", "READ_QC", ["fastqc", "trim", "host"], *GREEN)

    d.box("long", "Long reads (-L)\nONT or PacBio", "data", x(4), y(1.6), w=W, dashed=True)
    d.box(
        "assemble",
        "Assemble\nshort reads: MEGAHIT (default) or\nmetaSPAdes; with both, MEGAHIT\nadds what metaSPAdes left out\nlong reads: metaFlye",
        "step",
        x(3),
        y(1.6),
        w=W,
    )
    d.box("sort", "Sort contigs and drop\nthose under -l (1,000 bp)", "step", x(2), y(1.6), w=W)
    d.box("quast", "Assembly statistics\nQUAST", "step", x(1), y(1.6), w=W)
    d.box("assembly", "final_assembly.fasta\nand assembly report", "final", x(0), y(1.6), w=W)
    d.panel_around("p_asm", "ASSEMBLY", ["assemble", "sort", "quast"], *RED)

    d.box(
        "kraken",
        "Taxonomic profile\nKraken2, reads or contigs;\n--bracken: abundance at a rank",
        "step",
        x(3),
        y(2.95),
        w=W,
    )
    d.box("krona", "Kronagram", "output", x(2), y(2.95))
    d.box("nonpareil", "Sequencing coverage\nNonpareil, per sample", "step", x(3), y(3.75), w=W)
    d.box(
        "diversity",
        "Diversity coverage and\ncomparability across samples",
        "output",
        x(2) - 40,
        y(3.75),
    )
    d.panel_around("p_prof", "KRAKEN2 and NONPAREIL", ["kraken", "nonpareil"], *BLUE)

    for a, b in [
        ("raw", "fastqc"),
        ("fastqc", "trim"),
        ("trim", "host"),
        ("host", "clean"),
        ("long", "assemble"),
        ("assemble", "sort"),
        ("sort", "quast"),
        ("quast", "assembly"),
        ("kraken", "krona"),
        ("nonpareil", "diversity"),
    ]:
        d.arrow(a, b)
    d.arrow(
        "clean",
        "assemble",
        exit=(0.25, 1),
        entry=(0.5, 0),
        points=[(d.nodes["clean"].x + d.nodes["clean"].w * 0.25, y(0.8)), (x(3), y(0.8))],
    )
    edge = d.nodes["clean"].x + d.nodes["clean"].w + 40
    for target in ("kraken", "nonpareil"):
        d.arrow(
            "clean",
            target,
            exit=(1, 0.5),
            entry=(1, 0.5),
            points=[(edge, y(0)), (edge, d.nodes[target].y + d.nodes[target].h / 2)],
        )
    return d
