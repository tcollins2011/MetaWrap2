"""REASSEMBLE_BINS: reassemble each bin from its own reads and keep the best version."""

from _grid import BLUE, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("reassemble_bins", "REASSEMBLE_BINS module")
    d.box("bins", "Refined bins", "final", x(0), y(0), w=W)
    d.box("reads", "Clean reads (paired)", "data", x(0), y(1), w=W)
    d.box("nanopore", "Nanopore reads\n(optional)", "data", x(0), y(2), w=W, dashed=True)
    d.box(
        "recruit",
        "Recruit reads to bins\nbwa mem against the binned\ncontigs; both mates in one bin",
        "step",
        x(1),
        y(0.5),
        w=W,
    )
    d.box(
        "split",
        "Split each bin's reads\nstrict: < 2 mismatches\npermissive: < 5",
        "step",
        x(2),
        y(0.5),
        w=W,
    )
    d.box(
        "thin",
        "Thin deep bins to\n--max-coverage\n(50× by default, seeded)",
        "step",
        x(3),
        y(0.5),
        w=W,
    )
    d.box(
        "spades",
        "Reassemble each bin twice\nSPAdes --careful, with the\noriginal as untrusted contigs",
        "step",
        x(3),
        y(1.75),
        w=W,
    )
    d.box(
        "score",
        "Score each version\noriginal, strict, permissive\nCheckM",
        "step",
        x(2),
        y(1.75),
        w=W,
    )
    d.box(
        "choose",
        "Keep the best version\ncompleteness +\n5 × (100 − contamination)",
        "step",
        x(1),
        y(1.75),
        w=W,
    )
    d.panel_around(
        "panel", "REASSEMBLE_BINS", ["recruit", "split", "thin", "spades", "score", "choose"], *BLUE
    )
    d.box("final", "reassembled_bins/\nand .stats.tsv", "final", x(1), y(3), w=W)
    d.box("plots", "Reassembly plots: N50,\ncompleteness, contamination", "output", x(2) + 40, y(3))
    for a, b in [
        ("recruit", "split"),
        ("split", "thin"),
        ("thin", "spades"),
        ("spades", "score"),
        ("score", "choose"),
        ("choose", "final"),
        ("final", "plots"),
    ]:
        d.arrow(a, b)
    d.arrow("bins", "recruit", exit=(1, 0.5), entry=(0, 0.3))
    d.arrow("reads", "recruit", exit=(1, 0.5), entry=(0, 0.75))
    d.arrow(
        "nanopore", "spades", dashed=True, exit=(1, 0.5), entry=(0.5, 1), points=[(x(3), y(2.55))]
    )
    return d
