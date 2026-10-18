"""The function stack: genes, orthologs, pathways and resistance, per genome and per sample."""

from _grid import BLUE, GREEN, PURPLE, RED, YELLOW, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("function", "Functional annotation")
    d.box("bins", "Bins, or a\nwhole assembly", "final", x(0), y(0), w=W)
    d.box(
        "genes",
        "Call genes\nPyrodigal: trained per bin,\nmetagenomic mode for contigs",
        "step",
        x(1),
        y(0),
        w=W,
    )
    d.box(
        "annot",
        "Annotate the genes\nKOfam (--kofam) · Pfam (--pfam)\neggNOG (--eggnog) · dbCAN (--cazy)\nGECCO (--bgc) · Bakta (--bakta)",
        "step",
        x(2),
        y(0),
        w=W + 20,
    )
    d.panel_around("p_ann", "ANNOTATE", ["genes", "annot"], *RED)
    d.box(
        "tables",
        "kos.tsv, pfam_domains.tsv,\neggnog.tsv, cazymes.tsv,\ngene_clusters.tsv",
        "data",
        x(3) + 30,
        y(0),
        w=W + 20,
    )

    d.box(
        "weights",
        "Optional weighting\nbin abundance (-a), contig\ndepth (-C), gene abundance (-G)",
        "data",
        x(1),
        y(1.35),
        w=W,
    )
    d.box(
        "pathways",
        "KEGG module completeness,\npathwise and stepwise;\nfunctional hierarchy and Krona",
        "step",
        x(2),
        y(1.35),
        w=W + 20,
    )
    d.panel_around("p_path", "PATHWAYS", ["pathways"], *BLUE)
    d.box("o_path", "modules.tsv, hierarchy\ntables, Krona charts", "output", x(3) + 30, y(1.35))

    d.box("abund", "Bin abundance\nfrom QUANT_BINS", "data", x(1), y(2.7), w=W)
    d.box(
        "quantify",
        "KO and pathway abundance\nper sample, with the\nassignable fraction",
        "step",
        x(2),
        y(2.7),
        w=W + 20,
    )
    d.panel_around("p_qf", "QUANTIFY_FUNCTION", ["quantify"], *YELLOW)
    d.box("diff", "Test between groups\nANCOM-BC", "step", x(3) + 30, y(2.7), w=W)
    d.panel_around("p_diff", "DIFF_ABUNDANCE", ["diff"], *PURPLE)
    d.box("groups", "Sample groups", "data", x(4) + 30, y(2.7), w=W)
    d.box("o_diff", "Features that differ\nbetween groups", "output", x(3) + 30, y(3.8))

    d.box("cat_in", "ANNOTATE outputs from\nevery sample, and reads", "data", x(0), y(4.15), w=W)
    d.box(
        "catalogue",
        "Cluster genes across samples\nMMseqs2 linclust, 95% identity;\nquantify with salmon (--quantify)",
        "step",
        x(1),
        y(4.15),
        w=W + 20,
    )
    d.panel_around("p_cat", "CATALOGUE", ["catalogue"], *GREEN)
    d.box("o_cat", "Gene catalogue and\ngene abundance", "output", x(2) + 20, y(4.15))
    d.box(
        "amr_in",
        "Bins or assembly;\ngeNomad results\nfrom VIRUSES (optional)",
        "data",
        x(0),
        y(5.4),
        w=W,
    )
    d.box(
        "amr",
        "Resistance and virulence\nAMRFinderPlus, RGI (--rgi);\ncontext and evidence level",
        "step",
        x(1),
        y(5.4),
        w=W + 20,
    )
    d.panel_around("p_amr", "AMR", ["amr"], *RED)
    d.box("o_amr", "amr.tsv, virulence.tsv,\nharmonised tables", "output", x(2) + 20, y(5.4))

    for a, b in [
        ("bins", "genes"),
        ("genes", "annot"),
        ("annot", "tables"),
        ("weights", "pathways"),
        ("pathways", "o_path"),
        ("abund", "quantify"),
        ("quantify", "diff"),
        ("groups", "diff"),
        ("diff", "o_diff"),
        ("cat_in", "catalogue"),
        ("catalogue", "o_cat"),
        ("amr_in", "amr"),
        ("amr", "o_amr"),
    ]:
        d.arrow(a, b)
    d.arrow("pathways", "quantify", exit=(0.8, 1), entry=(0.8, 0))
    d.arrow(
        "tables",
        "pathways",
        exit=(0.5, 1),
        entry=(0.5, 0),
        points=[(x(3) + 30, y(0.68)), (x(2), y(0.68))],
    )
    return d
