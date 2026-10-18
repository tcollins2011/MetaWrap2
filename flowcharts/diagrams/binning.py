"""COVERAGE and BINNING: from assemblies and reads to one bin set per binner."""

from _grid import BLUE, RED, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("binning", "Coverage and binning")
    d.box("assemblies", "Every sample's\nassembly", "data", x(0), y(0), w=W)
    d.box(
        "concat",
        "Concatenate the assemblies,\nsample-prefixed contig names",
        "step",
        x(1),
        y(0),
        w=W,
    )
    d.box(
        "map_all",
        "Map every sample's reads once\nminimap2, or fairy k-mer sketches",
        "step",
        x(2),
        y(0),
        w=W,
    )
    d.box("split", "One depth table across\nsamples, split per assembly", "step", x(3), y(0), w=W)
    d.panel_around("p_cov", "COVERAGE (cohorts)", ["concat", "map_all", "split"], *BLUE)

    d.box("assembly", "Assembly", "data", x(0), y(1.8), w=W)
    d.box("reads", "Clean reads,\none set per sample", "data", x(0), y(3), w=W)
    d.box(
        "map",
        "Map each sample's reads\nminimap2 (default),\nbwa or bwa-mem2",
        "step",
        x(1),
        y(2.4),
        w=W,
    )
    d.box(
        "depth", "Depth table\nreads at ≥ 90% identity\n(--min-identity)", "step", x(2), y(1.8), w=W
    )
    d.box("classical", "metaBAT2 · MaxBin2\nbin from the depth table", "step", x(3), y(1.8), w=W)
    d.box(
        "bam",
        "CONCOCT · SemiBin2 · VAMB · COMEBin\nbin from the alignments;\n--device cpu, cuda or auto\nfor the learned three",
        "step",
        x(2.5),
        y(3),
        w=W + 110,
    )
    d.box("bins", "One bin set per binner\n<binner>_bins/", "bins", x(4.2), y(2.4), w=W)
    d.box(
        "presets",
        "--preset chooses the binners:\nvery-fast: metaBAT2\nfast: + MaxBin2\ndefault: metaBAT2 + SemiBin2\nslow: + COMEBin",
        "step",
        x(1),
        y(3.75),
        w=W,
        dashed=True,
    )
    d.panel_around("p_bin", "BINNING", ["map", "depth", "classical", "bam", "presets"], *RED)

    for a, b in [
        ("assemblies", "concat"),
        ("concat", "map_all"),
        ("map_all", "split"),
        ("depth", "classical"),
    ]:
        d.arrow(a, b)
    d.arrow("assembly", "map", exit=(1, 0.5), entry=(0, 0.3))
    d.arrow("reads", "map", exit=(1, 0.5), entry=(0, 0.7))
    d.arrow("map", "depth", exit=(1, 0.3), entry=(0, 0.5))
    d.arrow("map", "bam", exit=(1, 0.7), entry=(0, 0.5))
    d.arrow(
        "split",
        "depth",
        dashed=True,
        label="--depth",
        exit=(0.5, 1),
        entry=(0.5, 0),
        points=[(x(3), y(0.95)), (x(2), y(0.95))],
    )
    d.arrow("classical", "bins", exit=(1, 0.5), entry=(0, 0.3))
    d.arrow("bam", "bins", exit=(1, 0.5), entry=(0, 0.7))
    return d
