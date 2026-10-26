"""End-to-end diagram: read processing, the two genome recovery routes, and later steps."""

from _grid import BLUE, GREEN, PURPLE, RED, YELLOW, W, x, y
from drawio import Diagram


def build() -> Diagram:
    d = Diagram("end_to_end", "MetaWrap2 end to end")
    m = "module"
    d.box("sheet", "Sample sheet\nreads per sample", "data", x(0), y(0), w=W)
    d.box("qc", "READ_QC\n× each sample", m, x(1), y(0), w=W)
    d.box("clean", "Clean reads", "data", x(2), y(0), w=W)
    d.box("kraken", "KRAKEN2\ntaxonomic profile", m, x(3), y(0), w=W)
    d.box("nonpareil", "NONPAREIL\nsequencing coverage", m, x(4), y(0), w=W)
    d.panel_around("p_reads", "1  Reads", ["qc", "clean", "kraken", "nonpareil"], *GREEN)

    d.box("r_asm", "ASSEMBLY\none co-assembly", m, x(1), y(1.45), w=W)
    d.box("r_bin", "BINNING\nall samples' reads", m, x(2), y(1.45), w=W)
    d.box("r_ref", "BIN_REFINEMENT", m, x(3), y(1.45), w=W)
    d.box("r_rea", "REASSEMBLE_BINS", m, x(4), y(1.45), w=W)
    d.panel_around(
        "p_run",
        "2a  Genome recovery by co-assembly: metawrap2 run",
        ["r_asm", "r_bin", "r_ref", "r_rea"],
        *RED,
    )

    d.box("c_asm", "ASSEMBLY\n× each sample", m, x(1), y(2.75), w=W)
    d.box("c_cov", "COVERAGE\nall samples at once", m, x(2), y(2.75), w=W, dashed=True)
    d.box("c_bin", "BINNING and\nBIN_REFINEMENT\n× each sample", m, x(3), y(2.75), w=W)
    d.box("c_der", "DEREPLICATE\none MAG per species", m, x(4), y(2.75), w=W, dashed=True)
    d.panel_around(
        "p_coh",
        "2b  Genome recovery per sample: metawrap2 cohort",
        ["c_asm", "c_cov", "c_bin", "c_der"],
        *RED,
    )

    d.box("mags", "Final MAGs, or the MAG catalogue", "final", x(2.5), y(3.72), w=W + 120)

    d.box("quant", "QUANT_BINS\nabundance", m, x(1), y(4.6), w=W)
    d.box("classify", "CLASSIFY_BINS\ntaxonomy", m, x(2), y(4.6), w=W)
    d.box("appraise", "APPRAISE\nshare of community", m, x(3), y(4.6), w=W)
    d.box("blob", "BLOBOLOGY\nblobplots", m, x(4), y(4.6), w=W)
    d.panel_around("p_gen", "3  Genomes", ["quant", "classify", "appraise", "blob"], *YELLOW)

    d.box("annotate", "ANNOTATE\ngenes and orthologs", m, x(1), y(5.9), w=W)
    d.box("pathways", "PATHWAYS\nmodules and hierarchy", m, x(2), y(5.9), w=W)
    d.box("qfun", "QUANTIFY_FUNCTION\nper sample", m, x(3), y(5.9), w=W)
    d.box("diff", "DIFF_ABUNDANCE\nbetween groups", m, x(4), y(5.9), w=W)
    d.box("catalogue", "CATALOGUE\ngenes across samples", m, x(1), y(6.95), w=W, dashed=True)
    d.box("amr", "AMR\nresistance in context", m, x(2), y(6.95), w=W)
    d.panel_around(
        "p_fun", "4  Function", ["annotate", "pathways", "qfun", "diff", "catalogue", "amr"], *BLUE
    )

    d.box("viruses", "VIRUSES\nphage and hosts", m, x(1), y(8.35), w=W)
    d.box("strains", "STRAINS\nsharing between samples", m, x(2), y(8.35), w=W, dashed=True)
    d.box("pangenome", "PANGENOME\nper species", m, x(3), y(8.35), w=W, dashed=True)
    d.box("phylogeny", "PHYLOGENY\ntrees", m, x(4), y(8.35), w=W)
    d.panel_around(
        "p_bey", "5  Beyond the bins", ["viruses", "strains", "pangenome", "phylogeny"], *PURPLE
    )
    d.box("report", "REPORT\none PDF and mag_table.tsv", m, x(1), y(9.6), w=W)
    d.panel_around("p_rep", "6  Reading the result", ["report"], *GREEN)
    d.box("legend", "Dashed: multi-sample steps", "note", x(0), y(8.35), w=W)

    for a, b in [
        ("sheet", "qc"),
        ("qc", "clean"),
        ("clean", "kraken"),
        ("kraken", "nonpareil"),
        ("r_asm", "r_bin"),
        ("r_bin", "r_ref"),
        ("r_ref", "r_rea"),
        ("c_asm", "c_cov"),
        ("c_cov", "c_bin"),
        ("c_bin", "c_der"),
        ("annotate", "pathways"),
        ("pathways", "qfun"),
        ("qfun", "diff"),
        ("pangenome", "phylogeny"),
    ]:
        d.arrow(a, b)
    left = x(1) - W / 2 - 45
    d.arrow(
        "clean",
        "r_asm",
        exit=(0.5, 1),
        entry=(0, 0.5),
        points=[(x(2), y(0.72)), (left, y(0.72)), (left, y(1.45))],
    )
    d.arrow(
        "clean",
        "c_asm",
        exit=(0.5, 1),
        entry=(0, 0.5),
        points=[(x(2), y(0.72)), (left, y(0.72)), (left, y(2.75))],
    )
    right = x(4) + W / 2 + 45
    d.arrow(
        "r_rea", "mags", exit=(1, 0.5), entry=(1, 0.5), points=[(right, y(1.45)), (right, y(3.72))]
    )
    d.arrow("c_der", "mags", exit=(0.5, 1), entry=(1, 0.5), points=[(x(4), y(3.72))])
    for first, row in (
        ("quant", 4.6),
        ("annotate", 5.9),
        ("catalogue", 6.95),
        ("viruses", 8.35),
        ("report", 9.6),
    ):
        d.arrow(
            "mags", first, exit=(0, 0.5), entry=(0, 0.5), points=[(left, y(3.72)), (left, y(row))]
        )
    return d
