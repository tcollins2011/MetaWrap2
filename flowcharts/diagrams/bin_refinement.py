"""BIN_REFINEMENT: from up to six bin sets to one refined, chimera-checked set."""

from _grid import GREEN, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("bin_refinement", "BIN_REFINEMENT module")
    d.box("inputs", "Bin sets A–F\nfrom any binners", "bins", x(0), y(0), w=W)
    d.box(
        "prepare",
        "Prepare each set\nkeep bins of 50 kb–20 Mb,\ndrop unbinned contigs",
        "step",
        x(1),
        y(0),
        w=W,
    )
    d.box(
        "hybrids",
        "Hybrid bin sets\nBinning_refiner: AB, AC, BC, ABC\n(pairs and all, for more sets)",
        "step",
        x(2),
        y(0),
        w=W,
    )
    d.box(
        "score",
        "Score every set\nCheckM, or CheckM2 with --checkm2;\nidentical bins scored once",
        "step",
        x(3),
        y(0),
        w=W + 20,
    )
    d.box("checkm_db", "CheckM or\nCheckM2 data", "db", x(4), y(0))
    d.box(
        "consolidate",
        "Consolidate the sets\nbins overlapping ≥ 80% are one bin;\nkeep the higher of\ncompleteness − 5 × contamination",
        "step",
        x(3),
        y(1.25),
        w=W + 20,
    )
    d.box(
        "resolve", "Contigs in several bins\nkeep each in its best bin", "step", x(2), y(1.25), w=W
    )
    d.box(
        "rescore",
        "Re-score and filter\ncompleteness ≥ -c,\ncontamination ≤ -x",
        "step",
        x(1),
        y(1.25),
        w=W,
    )
    d.box(
        "gunc",
        "Flag chimeric bins\nGUNC, CSS > 0.45\n(--skip-gunc to turn off)",
        "step",
        x(1),
        y(2.5),
        w=W,
    )
    d.box("gunc_db", "GUNC data", "db", x(0), y(2.5))
    d.box(
        "dastool",
        "Alternative: --consolidator dastool\nDAS Tool scores the input sets once,\nreplacing hybrids, scoring\nand consolidation",
        "step",
        x(2.5),
        y(2.5),
        w=W + 120,
        dashed=True,
    )
    d.panel_around(
        "panel",
        "BIN_REFINEMENT",
        ["prepare", "hybrids", "score", "consolidate", "resolve", "rescore", "gunc", "dastool"],
        *GREEN,
    )
    d.box(
        "final",
        "metawrap_70_10_bins/\nand .stats.tsv\nwith MIMAG tags",
        "final",
        x(4),
        y(1.85),
        w=W,
    )
    d.box("plots", "binning_results.png", "output", x(4), y(3))
    for a, b in [
        ("inputs", "prepare"),
        ("prepare", "hybrids"),
        ("hybrids", "score"),
        ("score", "consolidate"),
        ("consolidate", "resolve"),
        ("resolve", "rescore"),
        ("rescore", "gunc"),
        ("final", "plots"),
    ]:
        d.arrow(a, b)
    d.arrow(
        "gunc",
        "final",
        exit=(1, 0.3),
        entry=(0, 0.5),
        points=[(x(1) + 220, d.nodes["gunc"].y + d.nodes["gunc"].h * 0.3), (x(1) + 220, y(1.85))],
    )
    d.arrow("checkm_db", "score", dashed=True)
    d.arrow("gunc_db", "gunc", dashed=True)
    return d
