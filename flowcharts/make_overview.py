#!/usr/bin/env python
"""Draw the MetaWrap2 module overview: the figure at the top of the README.

One layout, two outputs, so they cannot drift apart:

* ``docs/figures/metawrap2_overview.png`` - rendered here with matplotlib, no draw.io needed;
* ``flowcharts/metawrap2_overview.drawio`` - the same boxes and arrows as an editable draw.io
  diagram, for anyone who would rather move a box by hand than edit coordinates.

Box sizes are measured from their text rather than guessed, and the script refuses to save a
figure in which any label would overflow its box. Colours and shapes follow the original metaWRAP
figure: red process boxes for modules, green for reads and assemblies, yellow for bins, blue
ellipses for what a user looks at. Dashed modules are the multi-sample (cohort) steps.

    python flowcharts/make_overview.py
"""

from __future__ import annotations

import html
import os
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PNG = os.path.join(ROOT, "docs", "figures", "metawrap2_overview.png")
DRAWIO = os.path.join(HERE, "metawrap2_overview.drawio")

# ── look ──────────────────────────────────────────────────────────────────────────────────────

#: (fill, stroke) per kind, taken from the original metaWRAP figure.
COLORS = {
    "module": ("#FF6666", "#B85450"),
    "data": ("#D5E8D4", "#82B366"),
    "bins": ("#FFF2CC", "#D6B656"),
    "final": ("#FFFF66", "#D6B656"),
    "output": ("#DAE8FC", "#6C8EBF"),
}
FONT = "Arial"
FONT_SIZE = 19  # in diagram pixels
STROKE = 2.0
#: Inner padding between a label and its box edge, in diagram pixels. Process boxes add their
#: side bars on top of this, so text never touches the bars either.
PAD_X, PAD_Y = 18, 14
BAR = 0.09  # width of a process box's side bars, as a fraction of its width (draw.io's default)
MIN_SIZE = {"module": (230, 92), "data": (230, 92), "bins": (300, 72), "final": (230, 80)}
MIN_ELLIPSE = (190, 76)

# ── layout ────────────────────────────────────────────────────────────────────────────────────

C = [150, 480, 850, 1190, 1655, 2000]  # column centres
R = [62, 190, 330, 470, 610, 790, 925]  # row centres; the last two are the assembly-level row


@dataclass
class Node:
    key: str
    label: str
    kind: str
    x: float  # centre
    y: float  # centre
    dashed: bool = False
    w: float = 0.0
    h: float = 0.0

    @property
    def left(self) -> float:
        return self.x - self.w / 2

    @property
    def right(self) -> float:
        return self.x + self.w / 2

    @property
    def top(self) -> float:
        return self.y - self.h / 2

    @property
    def bottom(self) -> float:
        return self.y + self.h / 2


@dataclass
class Edge:
    points: List[Tuple[float, float]]
    dashed: bool = False
    arrow: bool = True
    extra: Dict[str, str] = field(default_factory=dict)


def nodes() -> Dict[str, Node]:
    n = [
        # reads to assembly
        Node("raw", "Raw reads from\nmetagenomic samples\n(short or long)", "data", C[0], R[1]),
        Node("read_qc", "READ_QC\nmodule", "module", C[0], R[2]),
        Node("clean", "Clean metagenomic\nreads", "data", C[0], R[3]),
        Node("kraken2", "KRAKEN2\nmodule", "module", C[0], R[4]),
        Node("krona", "Kronagrams", "output", C[0], R[5]),
        Node("nonpareil", "NONPAREIL\nmodule", "module", C[1], R[4]),
        Node("diversity", "Diversity\ncoverage", "output", C[1], R[5]),
        Node("assembly", "ASSEMBLY\nmodule", "module", C[1], R[3]),
        Node("asm", "Metagenomic\nassembly", "data", C[1], R[2]),
        Node("binning", "BINNING\nmodule", "module", C[1], R[1]),
        Node("coverage", "COVERAGE\nmodule", "module", C[1], R[0], dashed=True),
        # bins
        Node("classical", "Classical bin sets\nmetaBAT2 · MaxBin2 · CONCOCT", "bins", C[2], 128),
        Node("learned", "Learned bin sets\nSemiBin2 · VAMB · COMEBin", "bins", C[2], 252),
        Node("refine", "BIN_REFINEMENT\nmodule", "module", C[3], R[1]),
        Node("refine_qc", "Bin quality and\nchimera report", "output", C[3], R[0]),
        Node("final", "Final bins", "final", C[3], R[2]),
        Node("blobology", "BLOBOLOGY\nmodule", "module", C[2], R[3]),
        Node("blobplots", "Blobplots", "output", C[2], R[4]),
        Node("quant", "QUANT_BINS\nmodule", "module", C[3], R[3]),
        Node("abundance", "Bin abundance\ntable", "output", C[3], R[4]),
        # reassembly, then what the genomes are and what they do
        Node("reassemble", "REASSEMBLE_BINS\nmodule", "module", C[4], R[1]),
        Node("reassembly_qc", "QC plots", "output", C[4], R[0]),
        Node("reassembled", "Reassembled bins", "final", C[4], R[2]),
        Node("classify", "CLASSIFY_BINS\nmodule", "module", 1500, R[3]),
        Node("taxonomy", "Bin taxonomy", "output", 1500, R[4]),
        Node("annotate", "ANNOTATE\nmodule", "module", 1810, R[3]),
        Node("pathways", "PATHWAYS\nmodule", "module", 1810, R[4]),
        Node("function", "Pathway completeness\nand abundance", "output", 1810, R[5]),
        Node("amr", "AMR\nmodule", "module", 2130, R[3]),
        Node("resistance", "Resistance genes\nin context", "output", 2130, R[4]),
        # across a cohort
        Node("dereplicate", "DEREPLICATE\nmodule", "module", C[5], R[2], dashed=True),
        Node("catalogue", "Genome catalogue\nacross samples", "output", C[5], R[1]),
        Node(
            "cohort_genomes",
            "STRAINS · PANGENOME · PHYLOGENY\nmodules",
            "module",
            C[5] + 40,
            R[0],
            dashed=True,
        ),
        # straight from the assembly, with or without bins
        Node("viruses", "VIRUSES\nmodule", "module", C[2], R[5]),
        Node("viral", "Viral genomes,\nvOTUs and hosts", "output", C[2], R[6]),
        Node("gene_catalogue", "CATALOGUE\nmodule", "module", C[3], R[5], dashed=True),
        Node("genes", "Gene catalogue\nacross samples", "output", C[3], R[6]),
        Node("appraise", "APPRAISE\nmodule", "module", 1500, R[5]),
        Node("share", "Share of the\ncommunity binned", "output", 1500, R[6]),
    ]
    return {node.key: node for node in n}


def edges(n: Dict[str, Node]) -> List[Edge]:
    """Orthogonal arrows, routed by hand so none crosses another or passes through a box."""

    def down(a: str, b: str, dx: float = 0.0) -> Edge:
        return Edge([(n[a].x + dx, n[a].bottom), (n[a].x + dx, n[b].top)])

    def up(a: str, b: str, dashed: bool = False) -> Edge:
        return Edge([(n[a].x, n[a].top), (n[a].x, n[b].bottom)], dashed=dashed)

    def right(a: str, b: str) -> Edge:
        return Edge([(n[a].right, n[a].y), (n[b].left, n[a].y)])

    # Clean reads branch once, in the gap between the first two columns: up to COVERAGE, along to
    # ASSEMBLY, and down to NONPAREIL.
    split_x = (n["clean"].right + n["assembly"].left) / 2
    elbow_x = (n["binning"].right + n["classical"].left) / 2
    refine_entry = n["refine"].left
    final_to_reassemble_x = (n["final"].right + n["reassembled"].left) / 2
    under_reassembled = (n["reassembled"].bottom + n["classify"].top) / 2
    # The assembly-level row hangs off one bus, fed from the assembly's right side down the gap
    # between the assembly column and the bin sets.
    bus_x = n["asm"].right + 34
    bus_y = (n["blobplots"].bottom + n["viruses"].top) / 2

    return [
        down("raw", "read_qc"),
        down("read_qc", "clean"),
        down("clean", "kraken2"),
        down("kraken2", "krona"),
        right("clean", "assembly"),
        Edge(
            [
                (split_x, n["clean"].y),
                (split_x, n["coverage"].y),
                (n["coverage"].left, n["coverage"].y),
            ],
            dashed=True,
        ),
        Edge(
            [
                (split_x, n["clean"].y),
                (split_x, n["nonpareil"].y),
                (n["nonpareil"].left, n["nonpareil"].y),
            ]
        ),
        down("nonpareil", "diversity"),
        up("assembly", "asm"),
        up("asm", "binning"),
        Edge(
            [(n["coverage"].x, n["coverage"].bottom), (n["coverage"].x, n["binning"].top)],
            dashed=True,
        ),
        # BINNING -> bin sets -> BIN_REFINEMENT
        Edge(
            [
                (n["binning"].right, n["binning"].y - 22),
                (elbow_x, n["binning"].y - 22),
                (elbow_x, n["classical"].y),
                (n["classical"].left, n["classical"].y),
            ]
        ),
        Edge(
            [
                (n["binning"].right, n["binning"].y + 22),
                (elbow_x, n["binning"].y + 22),
                (elbow_x, n["learned"].y),
                (n["learned"].left, n["learned"].y),
            ]
        ),
        Edge(
            [
                (n["classical"].right, n["classical"].y),
                ((n["classical"].right + refine_entry) / 2, n["classical"].y),
                ((n["classical"].right + refine_entry) / 2, n["refine"].y - 22),
                (refine_entry, n["refine"].y - 22),
            ]
        ),
        Edge(
            [
                (n["learned"].right, n["learned"].y),
                ((n["learned"].right + refine_entry) / 2, n["learned"].y),
                ((n["learned"].right + refine_entry) / 2, n["refine"].y + 22),
                (refine_entry, n["refine"].y + 22),
            ]
        ),
        up("refine", "refine_qc"),
        down("refine", "final"),
        # Final bins -> BLOBOLOGY, QUANT_BINS, REASSEMBLE_BINS
        Edge(
            [
                (n["final"].left, n["final"].y),
                (n["blobology"].x, n["final"].y),
                (n["blobology"].x, n["blobology"].top),
            ]
        ),
        down("blobology", "blobplots"),
        down("final", "quant"),
        down("quant", "abundance"),
        Edge(
            [
                (n["final"].right, n["final"].y),
                (final_to_reassemble_x, n["final"].y),
                (final_to_reassemble_x, n["reassemble"].y),
                (n["reassemble"].left, n["reassemble"].y),
            ]
        ),
        up("reassemble", "reassembly_qc"),
        down("reassemble", "reassembled"),
        # Reassembled bins -> CLASSIFY_BINS and ANNOTATE; ANNOTATE -> PATHWAYS and AMR
        Edge(
            [
                (n["reassembled"].x - 60, n["reassembled"].bottom),
                (n["reassembled"].x - 60, under_reassembled),
                (n["classify"].x, under_reassembled),
                (n["classify"].x, n["classify"].top),
            ]
        ),
        Edge(
            [
                (n["reassembled"].x + 60, n["reassembled"].bottom),
                (n["reassembled"].x + 60, under_reassembled),
                (n["annotate"].x, under_reassembled),
                (n["annotate"].x, n["annotate"].top),
            ]
        ),
        down("classify", "taxonomy"),
        down("annotate", "pathways"),
        down("pathways", "function"),
        right("annotate", "amr"),
        down("amr", "resistance"),
        # across a cohort: reassembled bins -> DEREPLICATE -> genome catalogue -> strains etc.
        Edge(
            [
                (n["reassembled"].right, n["reassembled"].y),
                (n["dereplicate"].left, n["reassembled"].y),
            ],
            dashed=True,
        ),
        up("dereplicate", "catalogue", dashed=True),
        Edge(
            [
                (n["catalogue"].x, n["catalogue"].top),
                (n["catalogue"].x, n["cohort_genomes"].bottom),
            ],
            dashed=True,
        ),
        # Metagenomic assembly -> the assembly-level row
        Edge([(n["asm"].right, n["asm"].y), (bus_x, n["asm"].y), (bus_x, bus_y)], arrow=False),
        Edge([(bus_x, bus_y), (n["appraise"].x, bus_y)], arrow=False),
        Edge([(n["viruses"].x, bus_y), (n["viruses"].x, n["viruses"].top)]),
        Edge([(n["gene_catalogue"].x, bus_y), (n["gene_catalogue"].x, n["gene_catalogue"].top)]),
        Edge([(n["appraise"].x, bus_y), (n["appraise"].x, n["appraise"].top)]),
        down("viruses", "viral"),
        down("gene_catalogue", "genes"),
        down("appraise", "share"),
    ]


# ── measuring and drawing ─────────────────────────────────────────────────────────────────────

WIDTH, HEIGHT = 2290, 990
DPI = 100
SCALE = 2  # render at 2x for a crisp README image


def _figure():
    fig = plt.figure(figsize=(WIDTH / DPI, HEIGHT / DPI), dpi=DPI)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(HEIGHT, 0)
    ax.axis("off")
    return fig, ax


def _points(px: float) -> float:
    """Diagram pixels to font points, at this figure's DPI."""
    return px * 72.0 / DPI


def _text_size(ax, label: str) -> Tuple[float, float]:
    """Width and height a label occupies, in diagram pixels."""
    t = ax.text(
        0,
        0,
        label,
        fontsize=_points(FONT_SIZE),
        family=FONT,
        ha="center",
        va="center",
        linespacing=1.25,
    )
    renderer = ax.figure.canvas.get_renderer()
    box = t.get_window_extent(renderer=renderer)
    t.remove()
    return box.width, box.height  # 1 figure pixel == 1 diagram pixel at this DPI


def size_nodes(ax, n: Dict[str, Node]) -> None:
    """Grow every box until its label fits inside it with padding to spare."""
    for node in n.values():
        tw, th = _text_size(ax, node.label)
        if node.kind == "output":
            # Text inside an ellipse fits when its corners lie inside: (x/a)^2 + (y/b)^2 <= 1. With
            # the semi-axes sqrt(2) times the padded half-sizes, the corners land exactly on it.
            w = max(MIN_ELLIPSE[0], (tw + 2 * PAD_X) * 1.4142)
            h = max(MIN_ELLIPSE[1], (th + 2 * PAD_Y) * 1.4142)
        else:
            min_w, min_h = MIN_SIZE[node.kind]
            needed_w = tw + 2 * PAD_X
            if node.kind == "module":
                needed_w = needed_w / (1 - 2 * BAR)  # the side bars are not text space
            w = max(min_w, needed_w)
            h = max(min_h, th + 2 * PAD_Y)
        node.w, node.h = round(w), round(h)
    # Boxes stacked in one column read as a set only if they line up, so each stack shares the
    # widest member's size.
    for group in (("classical", "learned"),):
        widest = max(n[k].w for k in group)
        for k in group:
            n[k].w = widest


def check_containment(ax, n: Dict[str, Node]) -> None:
    """Refuse to save a figure with a label outside its box."""
    problems = []
    for node in n.values():
        tw, th = _text_size(ax, node.label)
        if node.kind == "output":
            a, b = node.w / 2, node.h / 2
            if (tw / 2 / a) ** 2 + (th / 2 / b) ** 2 > 1.0:
                problems.append(node.key)
        else:
            inner = node.w * (1 - 2 * BAR) if node.kind == "module" else node.w
            if tw > inner - 2 or th > node.h - 2:
                problems.append(node.key)
    overlaps = [
        (a.key, b.key)
        for i, a in enumerate(n.values())
        for b in list(n.values())[i + 1 :]
        if a.left < b.right and b.left < a.right and a.top < b.bottom and b.top < a.bottom
    ]
    if problems or overlaps:
        raise SystemExit("Label overflow: %s; overlapping boxes: %s" % (problems, overlaps))


def draw(ax, n: Dict[str, Node], e: List[Edge]) -> None:
    for edge in e:
        xs, ys = zip(*edge.points)
        style = (0, (6, 4)) if edge.dashed else "solid"
        ax.plot(
            xs[:-1] + (xs[-1],),
            ys[:-1] + (ys[-1],),
            color="#000000",
            lw=STROKE * 0.75,
            linestyle=style,
            solid_capstyle="butt",
            zorder=1,
        )
        if not edge.arrow:
            continue  # a bus segment: the arrowheads are on the drops that leave it
        ax.add_patch(
            FancyArrowPatch(
                edge.points[-2],
                edge.points[-1],
                arrowstyle="-|>",
                mutation_scale=22,
                color="#000000",
                lw=0,
                shrinkA=0,
                shrinkB=0,
                zorder=2,
            )
        )
    for node in n.values():
        fill, stroke = COLORS[node.kind]
        ls = (0, (6, 3)) if node.dashed else "solid"
        if node.kind == "output":
            ax.add_patch(
                Ellipse(
                    (node.x, node.y),
                    node.w,
                    node.h,
                    facecolor=fill,
                    edgecolor=stroke,
                    lw=STROKE,
                    linestyle=ls,
                    zorder=3,
                )
            )
        else:
            ax.add_patch(
                Rectangle(
                    (node.left, node.top),
                    node.w,
                    node.h,
                    facecolor=fill,
                    edgecolor=stroke,
                    lw=STROKE,
                    linestyle=ls,
                    zorder=3,
                )
            )
            if node.kind == "module":
                for frac in (BAR, 1 - BAR):
                    x = node.left + node.w * frac
                    ax.plot(
                        [x, x],
                        [node.top, node.bottom],
                        color=stroke,
                        lw=STROKE,
                        linestyle=ls,
                        zorder=4,
                    )
        ax.text(
            node.x,
            node.y,
            node.label,
            fontsize=_points(FONT_SIZE),
            family=FONT,
            ha="center",
            va="center",
            linespacing=1.25,
            zorder=5,
            color="#000000",
        )
    _legend(ax)


def _legend(ax) -> None:
    """One line of key, bottom right: what a dashed box means. The colours are self-evident."""
    x, y = 1720, 950
    ax.add_patch(
        Rectangle(
            (x, y - 14),
            46,
            28,
            facecolor=COLORS["module"][0],
            edgecolor=COLORS["module"][1],
            lw=STROKE,
            linestyle=(0, (6, 3)),
            zorder=3,
        )
    )
    ax.text(
        x + 60,
        y,
        "multi-sample steps, run across a whole cohort",
        fontsize=_points(16),
        family=FONT,
        ha="left",
        va="center",
        color="#333333",
    )


# ── draw.io ───────────────────────────────────────────────────────────────────────────────────

_DRAWIO_STYLE = {
    "module": "shape=process;whiteSpace=wrap;html=1;size=%g;" % BAR,
    "data": "rounded=0;whiteSpace=wrap;html=1;",
    "bins": "rounded=0;whiteSpace=wrap;html=1;",
    "final": "rounded=0;whiteSpace=wrap;html=1;",
    "output": "ellipse;whiteSpace=wrap;html=1;",
}


def write_drawio(n: Dict[str, Node], e: List[Edge], path: str) -> None:
    cells = ['<mxCell id="0"/>', '<mxCell id="1" parent="0"/>']
    for node in n.values():
        fill, stroke = COLORS[node.kind]
        style = "%sfillColor=%s;strokeColor=%s;strokeWidth=2;fontFamily=%s;fontSize=%d;%s" % (
            _DRAWIO_STYLE[node.kind],
            fill,
            stroke,
            FONT,
            FONT_SIZE,
            "dashed=1;" if node.dashed else "",
        )
        value = html.escape(node.label).replace("\n", "&lt;br&gt;")
        cells.append(
            '<mxCell id="%s" value="%s" style="%s" vertex="1" parent="1">'
            '<mxGeometry x="%g" y="%g" width="%g" height="%g" as="geometry"/></mxCell>'
            % (node.key, value, style, node.left, node.top, node.w, node.h)
        )
    for i, edge in enumerate(e):
        style = "endArrow=%s;html=1;rounded=0;strokeWidth=1.5;%s" % (
            "classic" if edge.arrow else "none",
            "dashed=1;" if edge.dashed else "",
        )
        (sx, sy), (tx, ty) = edge.points[0], edge.points[-1]
        waypoints = "".join('<mxPoint x="%g" y="%g"/>' % p for p in edge.points[1:-1])
        cells.append(
            '<mxCell id="e%d" style="%s" edge="1" parent="1"><mxGeometry relative="1" as="geometry">'
            '<mxPoint x="%g" y="%g" as="sourcePoint"/><mxPoint x="%g" y="%g" as="targetPoint"/>'
            '<Array as="points">%s</Array></mxGeometry></mxCell>'
            % (i, style, sx, sy, tx, ty, waypoints)
        )
    xml = (
        '<mxfile host="metawrap2"><diagram id="overview" name="MetaWrap2 overview">'
        '<mxGraphModel dx="%d" dy="%d" grid="1" gridSize="10" page="1" pageWidth="%d" pageHeight="%d">'
        "<root>%s</root></mxGraphModel></diagram></mxfile>\n"
        % (WIDTH, HEIGHT, WIDTH, HEIGHT, "".join(cells))
    )
    with open(path, "w") as fh:
        fh.write(xml)


def _shrink(path: str) -> None:
    """Re-save as a 128-colour palette PNG: about a third of the size, no visible change.

    The diagram is a few flat colours plus anti-aliased edges, so a palette loses nothing a reader
    can see, and it keeps the image under the repository's 512 KB file limit.
    """
    from PIL import Image

    image = Image.open(path).convert("RGB")
    image.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(
        path, optimize=True
    )


def main() -> int:
    fig, ax = _figure()
    n = nodes()
    size_nodes(ax, n)
    check_containment(ax, n)
    e = edges(n)
    draw(ax, n, e)
    os.makedirs(os.path.dirname(PNG), exist_ok=True)
    fig.savefig(PNG, dpi=DPI * SCALE, facecolor="white")
    plt.close(fig)
    _shrink(PNG)
    write_drawio(n, e, DRAWIO)
    print("wrote %s\nwrote %s" % (os.path.relpath(PNG, ROOT), os.path.relpath(DRAWIO, ROOT)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
