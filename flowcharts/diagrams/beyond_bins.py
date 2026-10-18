"""VIRUSES, STRAINS, PANGENOME and PHYLOGENY: what the sample holds besides the bins themselves."""

from _grid import BLUE, GREEN, PURPLE, RED, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("beyond_bins", "Beyond the bins")
    d.box("v_in", "Assembly (free phage)\nor bins (proviruses)", "data", x(0), y(0), w=W)
    d.box("genomad", "Find viral contigs\ngeNomad", "step", x(1), y(0), w=W)
    d.box("checkv", "Estimate completeness\nCheckV", "step", x(2), y(0), w=W)
    d.box(
        "votus",
        "Cluster into vOTUs\nskani, 95% ANI, 85% coverage;\nCRISPR hosts (--crispr-hosts)",
        "step",
        x(3),
        y(0),
        w=W + 20,
    )
    d.panel_around("p_vir", "VIRUSES", ["genomad", "checkv", "votus"], *RED)
    d.box("o_vir", "Viral genomes, vOTUs\nand host links", "output", x(4) + 20, y(0))

    d.box("s_in", "Genomes, and each\nsample's reads", "data", x(0), y(1.45), w=W)
    d.box("map", "Map each sample\nminimap2", "step", x(1), y(1.45), w=W)
    d.box(
        "profile",
        "Profile each population\ninStrain (default) or TRACS",
        "step",
        x(2),
        y(1.45),
        w=W,
    )
    d.box(
        "compare",
        "Compare samples\npopANI, gated on\nshared breadth",
        "step",
        x(3),
        y(1.45),
        w=W + 20,
    )
    d.panel_around("p_str", "STRAINS", ["map", "profile", "compare"], *BLUE)
    d.box("o_str", "Strain sharing\nbetween samples", "output", x(4) + 20, y(1.45))

    d.box("p_in", "Dereplication clusters\nand the bin folders", "data", x(0), y(2.9), w=W)
    d.box("ppang", "One pangenome per species\nPPanGGOLiN, ≥ 3 genomes", "step", x(1), y(2.9), w=W)
    d.box("msa", "Core-gene alignment\n(--msa)", "step", x(2), y(2.9), w=W)
    d.panel_around("p_pan", "PANGENOME", ["ppang", "msa"], *GREEN)
    d.box(
        "tree",
        "Build a tree with support\nIQ-TREE, model chosen;\nFastTree with --fast",
        "step",
        x(3) + 40,
        y(2.9),
        w=W + 20,
    )
    d.panel_around("p_phy", "PHYLOGENY", ["tree"], *PURPLE)
    d.box("o_tree", "Trees with\nsupport values", "output", x(4) + 50, y(2.9))
    d.box(
        "gtdb_aln", "or a GTDB-Tk marker\nalignment", "data", x(3) + 40, y(4.05), w=W, dashed=True
    )
    d.box("o_pan", "pangenomes.tsv\nand partitions", "output", x(1), y(4.05))

    for a, b in [
        ("v_in", "genomad"),
        ("genomad", "checkv"),
        ("checkv", "votus"),
        ("votus", "o_vir"),
        ("s_in", "map"),
        ("map", "profile"),
        ("profile", "compare"),
        ("compare", "o_str"),
        ("p_in", "ppang"),
        ("ppang", "msa"),
        ("msa", "tree"),
        ("tree", "o_tree"),
        ("ppang", "o_pan"),
    ]:
        d.arrow(a, b)
    d.arrow("gtdb_aln", "tree", dashed=True)
    return d
