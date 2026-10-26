"""The module overview: the figure at the top of the README."""

from drawio import Diagram

C = [150, 480, 850, 1190, 1655, 2000]  # column centres
R = [62, 190, 330, 470, 610, 790, 925]  # row centres


def build() -> Diagram:
    d = Diagram("overview", "MetaWrap2 modules")
    # reads to assembly
    d.box("raw", "Raw reads from\nmetagenomic samples\n(short or long)", "data", C[0], R[1])
    d.box("read_qc", "READ_QC\nmodule", "module", C[0], R[2])
    d.box("clean", "Clean metagenomic\nreads", "data", C[0], R[3])
    d.box("kraken2", "KRAKEN2\nmodule", "module", C[0], R[4])
    d.box("krona", "Kronagrams", "output", C[0], R[5])
    d.box("coverage", "COVERAGE\nmodule", "module", C[1], R[0], dashed=True)
    d.box("binning", "BINNING\nmodule", "module", C[1], R[1])
    d.box("asm", "Metagenomic\nassembly", "data", C[1], R[2])
    d.box("assembly", "ASSEMBLY\nmodule", "module", C[1], R[3])
    d.box("nonpareil", "NONPAREIL\nmodule", "module", C[1], R[4])
    d.box("diversity", "Diversity\ncoverage", "output", C[1], R[5])
    # bins
    d.box("classical", "Classical bin sets\nmetaBAT2 · MaxBin2 · CONCOCT", "bins", C[2], 128, w=310)
    d.box("learned", "Learned bin sets\nSemiBin2 · VAMB · COMEBin", "bins", C[2], 252, w=310)
    d.box("refine_qc", "Bin quality and\nchimera report", "output", C[3], R[0])
    d.box("refine", "BIN_REFINEMENT\nmodule", "module", C[3], R[1])
    d.box("final", "Final bins", "final", C[3], R[2])
    d.box("blobology", "BLOBOLOGY\nmodule", "module", C[2], R[3])
    d.box("blobplots", "Blobplots", "output", C[2], R[4])
    d.box("quant", "QUANT_BINS\nmodule", "module", C[3], R[3])
    d.box("abundance", "Bin abundance\ntable", "output", C[3], R[4])
    # reassembly and what follows
    d.box("reassembly_qc", "QC plots", "output", C[4], R[0])
    d.box("reassemble", "REASSEMBLE_BINS\nmodule", "module", C[4], R[1])
    d.box("reassembled", "Reassembled bins", "final", C[4], R[2])
    d.box("classify", "CLASSIFY_BINS\nmodule", "module", 1500, R[3])
    d.box("taxonomy", "Bin taxonomy", "output", 1500, R[4])
    d.box("annotate", "ANNOTATE\nmodule", "module", 1810, R[3])
    d.box("pathways", "PATHWAYS\nmodule", "module", 1810, R[4])
    d.box("function", "Pathway completeness\nand abundance", "output", 1810, R[5])
    d.box("amr", "AMR\nmodule", "module", 2130, R[3])
    d.box("resistance", "Resistance genes\nin context", "output", 2130, R[4])
    # across a cohort
    d.box(
        "cohort", "STRAINS · PANGENOME · PHYLOGENY\nmodules", "module", C[5] + 40, R[0], dashed=True
    )
    d.box("catalogue", "Genome catalogue\nacross samples", "output", C[5], R[1])
    d.box("dereplicate", "DEREPLICATE\nmodule", "module", C[5], R[2], dashed=True)
    # straight from the assembly
    d.box("viruses", "VIRUSES\nmodule", "module", C[2], R[5])
    d.box("viral", "Viral genomes,\nvOTUs and hosts", "output", C[2], R[6])
    d.box("gene_catalogue", "CATALOGUE\nmodule", "module", C[3], R[5], dashed=True)
    d.box("genes", "Gene catalogue\nacross samples", "output", C[3], R[6])
    d.box("appraise", "APPRAISE\nmodule", "module", 1500, R[5])
    d.box("share", "Share of the\ncommunity binned", "output", 1500, R[6])
    d.box(
        "legend",
        "Dashed boxes: multi-sample steps, run across a whole cohort",
        "note",
        1950,
        R[6],
        w=360,
    )

    # Give every box the same minimum size, grown about its centre so positions stay fixed.
    for node in d.nodes.values():
        if node.kind == "note":
            continue
        width = max(node.w, 250 if node.kind != "output" else 240)
        height = max(node.h, 78)
        node.x -= (width - node.w) / 2
        node.y -= (height - node.h) / 2
        node.w, node.h = width, height

    n = d.nodes
    split_x = (n["clean"].x + n["clean"].w + n["assembly"].x) / 2
    for a, b in [
        ("raw", "read_qc"),
        ("read_qc", "clean"),
        ("clean", "kraken2"),
        ("kraken2", "krona"),
        ("clean", "assembly"),
        ("assembly", "asm"),
        ("asm", "binning"),
        ("nonpareil", "diversity"),
        ("refine", "refine_qc"),
        ("refine", "final"),
        ("blobology", "blobplots"),
        ("final", "quant"),
        ("quant", "abundance"),
        ("reassemble", "reassembly_qc"),
        ("reassemble", "reassembled"),
        ("classify", "taxonomy"),
        ("annotate", "pathways"),
        ("pathways", "function"),
        ("annotate", "amr"),
        ("amr", "resistance"),
        ("viruses", "viral"),
        ("gene_catalogue", "genes"),
        ("appraise", "share"),
    ]:
        d.arrow(a, b)
    d.arrow(
        "clean",
        "coverage",
        dashed=True,
        exit=(1, 0.5),
        entry=(0, 0.5),
        points=[(split_x, R[3]), (split_x, R[0])],
    )
    d.arrow(
        "clean",
        "nonpareil",
        exit=(1, 0.5),
        entry=(0, 0.5),
        points=[(split_x, R[3]), (split_x, R[4])],
    )
    d.arrow("coverage", "binning", dashed=True)
    d.arrow("binning", "classical", exit=(1, 0.25), entry=(0, 0.5))
    d.arrow("binning", "learned", exit=(1, 0.75), entry=(0, 0.5))
    d.arrow("classical", "refine", exit=(1, 0.5), entry=(0, 0.25))
    d.arrow("learned", "refine", exit=(1, 0.5), entry=(0, 0.75))
    d.arrow("final", "blobology", exit=(0, 0.5), entry=(0.5, 0))
    d.arrow("final", "reassemble", exit=(1, 0.5), entry=(0, 0.5))
    d.arrow("reassembled", "classify", exit=(0.25, 1), entry=(0.5, 0))
    d.arrow("reassembled", "annotate", exit=(0.75, 1), entry=(0.5, 0))
    d.arrow("reassembled", "dereplicate", dashed=True, exit=(1, 0.5), entry=(0, 0.5))
    d.arrow("dereplicate", "catalogue", dashed=True)
    d.arrow("catalogue", "cohort", dashed=True)
    bus_x = n["asm"].x + n["asm"].w + 34
    bus_y = (n["blobplots"].y + n["blobplots"].h + n["viruses"].y) / 2
    for target in ("viruses", "gene_catalogue", "appraise"):
        d.arrow(
            "asm",
            target,
            exit=(1, 0.5),
            entry=(0.5, 0),
            points=[(bus_x, R[2]), (bus_x, bus_y), (n[target].x + n[target].w / 2, bus_y)],
        )
    return d
