"""DEREPLICATE, QUANT_BINS, CLASSIFY_BINS, APPRAISE and BLOBOLOGY: from bins to a described catalogue."""

from _grid import BLUE, GREEN, PURPLE, RED, YELLOW, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("genome_catalogue", "Genome catalogue")
    d.box("bins", "Refined or reassembled\nbins, per sample", "final", x(0), y(0), w=W)
    d.box("ani", "Pairwise ANI\nskani triangle", "step", x(1), y(0), w=W)
    d.box("cluster", "Cluster at 95% ANI,\n≥ 10% aligned fraction", "step", x(2), y(0), w=W)
    d.box("best", "Keep the best member\nof each cluster by CheckM score", "step", x(3), y(0), w=W)
    d.panel_around("p_derep", "DEREPLICATE (cohorts)", ["ani", "cluster", "best"], *RED)
    d.box("catalogue", "MAG catalogue\ndereplicated_genomes/", "final", x(4), y(0), w=W)

    row = 1.55
    d.box(
        "quant",
        "Abundance per sample\nsylph (default), salmon,\nCoverM or the depth table",
        "step",
        x(1),
        y(row),
        w=W,
    )
    d.box(
        "gtdbtk",
        "GTDB-Tk classify_wf\n(--gtdbtk, recommended), or\nMEGABLAST and taxator-tk\n(default, deprecated)",
        "step",
        x(2),
        y(row),
        w=W,
    )
    d.box(
        "appraise",
        "Share of the community\nin the genomes\nSingleM on reads and genomes",
        "step",
        x(3),
        y(row),
        w=W,
    )
    d.box(
        "blob",
        "GC against abundance\nper contig, with MEGABLAST\ntaxonomy and bin labels",
        "step",
        x(4),
        y(row),
        w=W,
    )
    d.panel_around("p_quant", "QUANT_BINS", ["quant"], *BLUE)
    d.panel_around("p_class", "CLASSIFY_BINS", ["gtdbtk"], *YELLOW)
    d.panel_around("p_appr", "APPRAISE", ["appraise"], *GREEN)
    d.panel_around("p_blob", "BLOBOLOGY", ["blob"], *PURPLE)

    out = row + 1.35
    d.box("o_quant", "bin_abundance_table.tsv\nand heatmap", "output", x(1), y(out))
    d.box("o_class", "bin_taxonomy.tsv", "output", x(2), y(out))
    d.box("o_appr", "Fraction represented,\nmissing lineages", "output", x(3), y(out))
    d.box("o_blob", "Blobplots", "output", x(4), y(out))
    d.box("reads", "Clean reads", "data", x(0), y(row), w=W)

    for a, b in [
        ("bins", "ani"),
        ("ani", "cluster"),
        ("cluster", "best"),
        ("best", "catalogue"),
        ("reads", "quant"),
        ("blob", "o_blob"),
        ("appraise", "o_appr"),
    ]:
        d.arrow(a, b)
    d.arrow("quant", "o_quant")
    d.arrow("gtdbtk", "o_class")
    bus = y(0.8)
    for target in ("quant", "gtdbtk", "appraise", "blob"):
        d.arrow(
            "catalogue",
            target,
            exit=(0.5, 1),
            entry=(0.5, 0),
            points=[(x(4), bus), (d.nodes[target].x + d.nodes[target].w / 2, bus)],
        )
    return d
