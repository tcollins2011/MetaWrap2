"""A small layer for writing draw.io diagrams from Python and exporting them with draw.io itself.

Each diagram in ``flowcharts/diagrams/`` builds a :class:`Diagram`: boxes at grid positions, sized
from their text, and arrows attached to the boxes they connect so draw.io routes them. ``save()``
writes the ``.drawio`` file, which opens and edits in draw.io like any other, and ``export()`` renders
it with the draw.io desktop app's command-line exporter.

The colours and shapes follow the original metaWRAP figures: red process boxes for modules, green
for reads and assemblies, yellow for bins, blue ellipses for results, plus white boxes for the steps
inside a module and cylinders for reference databases.
"""

from __future__ import annotations

import html
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

FONT = "Helvetica"
FONT_SIZE = 16
STROKE = 2

#: kind -> (fill, stroke, extra style)
STYLES: Dict[str, Tuple[str, str, str]] = {
    "module": ("#FF6666", "#B85450", "shape=process;size=0.08;"),
    "data": ("#D5E8D4", "#82B366", "rounded=0;"),
    "bins": ("#FFF2CC", "#D6B656", "rounded=0;"),
    "final": ("#FFFF66", "#D6B656", "rounded=0;"),
    "output": ("#DAE8FC", "#6C8EBF", "ellipse;"),
    "step": ("#FFFFFF", "#666666", "rounded=1;arcSize=12;"),
    "tool": ("#F5F5F5", "#666666", "rounded=1;arcSize=12;fontStyle=2;"),
    "db": ("#E1D5E7", "#9673A6", "shape=cylinder3;boundedLbl=1;backgroundOutline=1;size=10;"),
    "decision": ("#FFFFFF", "#666666", "rhombus;"),
    "note": ("none", "none", "text;align=left;verticalAlign=top;"),
}

#: Container (module boundary) style: a rounded panel with its title at the top.
CONTAINER = (
    "rounded=1;arcSize=3;fillColor=%s;strokeColor=%s;strokeWidth=2;verticalAlign=top;"
    "align=left;spacingLeft=12;spacingTop=6;fontStyle=1;fontSize=16;fontFamily=%s;container=0;"
)

PAD_X, PAD_Y = 16, 12
MIN_W, MIN_H = 120, 48


def _text_size(label: str, size: int = FONT_SIZE, bold: bool = False) -> Tuple[float, float]:
    """Width and height of *label* in pixels, measured with matplotlib's Helvetica-like font."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    lines = label.split("\n")
    prop = FontProperties(
        family=["Helvetica", "Arial", "DejaVu Sans"], size=size, weight="bold" if bold else "normal"
    )
    width = 0.0
    for line in lines:
        if line:
            width = max(width, TextPath((0, 0), line, prop=prop).get_extents().width)
    # 1 pt == 1 px in draw.io's coordinate space; 1.25 line height.
    return width * 1.06, len(lines) * size * 1.25


@dataclass
class Node:
    id: str
    label: str
    kind: str
    x: float
    y: float
    w: float
    h: float
    dashed: bool = False
    parent: str = "1"


@dataclass
class Edge:
    source: str
    target: str
    label: str = ""
    dashed: bool = False
    points: List[Tuple[float, float]] = field(default_factory=list)
    exit: Optional[Tuple[float, float]] = None
    entry: Optional[Tuple[float, float]] = None


class Diagram:
    def __init__(self, name: str, title: str = ""):
        self.name = name
        self.title = title
        self.nodes: Dict[str, Node] = {}
        self.edges: List[Edge] = []
        self.containers: List[Tuple[str, str, float, float, float, float, str, str]] = []

    # -- building ----------------------------------------------------------------------------

    def box(
        self,
        id: str,
        label: str,
        kind: str,
        x: float,
        y: float,
        w: Optional[float] = None,
        h: Optional[float] = None,
        dashed: bool = False,
    ) -> Node:
        """A box centred on (x, y), at least as large as its text needs."""
        tw, th = _text_size(label)
        need_w, need_h = tw + 2 * PAD_X, th + 2 * PAD_Y
        if kind == "module":
            need_w /= 0.84  # the process shape's side bars are not text space
        if kind == "output":
            need_w, need_h = need_w * 1.38, need_h * 1.38
        if kind == "decision":
            need_w, need_h = need_w * 1.9, need_h * 1.9
        if kind == "db":
            need_h += 20
        width = max(w or 0, need_w, MIN_W)
        height = max(h or 0, need_h, MIN_H)
        node = Node(id, label, kind, x - width / 2, y - height / 2, width, height, dashed)
        self.nodes[id] = node
        return node

    def panel(
        self,
        id: str,
        title: str,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        fill: str,
        stroke: str,
    ) -> None:
        """A labelled module boundary from (x0, y0) to (x1, y1), drawn behind everything else."""
        self.containers.append((id, title, x0, y0, x1 - x0, y1 - y0, fill, stroke))

    def panel_around(
        self, id: str, title: str, members, fill: str, stroke: str, pad: float = 18, top: float = 40
    ) -> None:
        """A module boundary drawn tight around *members*, with room for its title on top."""
        nodes = [self.nodes[m] for m in members]
        x0 = min(n.x for n in nodes) - pad
        y0 = min(n.y for n in nodes) - top
        x1 = max(n.x + n.w for n in nodes) + pad
        y1 = max(n.y + n.h for n in nodes) + pad
        self.panel(id, title, x0, y0, x1, y1, fill, stroke)

    def arrow(self, source: str, target: str, **kw) -> None:
        self.edges.append(Edge(source, target, **kw))

    # -- writing -----------------------------------------------------------------------------

    def _check(self) -> None:
        """No two boxes may overlap; a diagram that breaks this is a layout bug, not a style."""
        nodes = list(self.nodes.values())
        for i, a in enumerate(nodes):
            for b in nodes[i + 1 :]:
                if a.x < b.x + b.w and b.x < a.x + a.w and a.y < b.y + b.h and b.y < a.y + a.h:
                    raise SystemExit("%s: boxes overlap: %s and %s" % (self.name, a.id, b.id))

    def xml(self) -> str:
        self._check()
        cells = ['<mxCell id="0"/>', '<mxCell id="1" parent="0"/>']
        for cid, title, x, y, w, h, fill, stroke in self.containers:
            cells.append(
                '<mxCell id="%s" value="%s" style="%s" vertex="1" parent="1">'
                '<mxGeometry x="%g" y="%g" width="%g" height="%g" as="geometry"/></mxCell>'
                % (cid, html.escape(title), CONTAINER % (fill, stroke, FONT), x, y, w, h)
            )
        for node in self.nodes.values():
            fill, stroke, extra = STYLES[node.kind]
            style = (
                "%swhiteSpace=wrap;html=1;fillColor=%s;strokeColor=%s;strokeWidth=%d;"
                "fontFamily=%s;fontSize=%d;fontColor=#000000;%s"
                % (extra, fill, stroke, STROKE, FONT, FONT_SIZE, "dashed=1;" if node.dashed else "")
            )
            # html=1 labels are parsed as HTML after the XML is read, so text is escaped twice:
            # once for HTML (so "<binner>" stays text) and once for the XML attribute.
            value = html.escape(html.escape(node.label, quote=False)).replace("\n", "&lt;br&gt;")
            cells.append(
                '<mxCell id="%s" value="%s" style="%s" vertex="1" parent="1">'
                '<mxGeometry x="%g" y="%g" width="%g" height="%g" as="geometry"/></mxCell>'
                % (node.id, value, style, node.x, node.y, node.w, node.h)
            )
        for i, edge in enumerate(self.edges):
            style = (
                "edgeStyle=orthogonalEdgeStyle;rounded=1;html=1;endArrow=block;endFill=1;"
                "strokeColor=#333333;strokeWidth=1.6;fontFamily=%s;fontSize=12;"
                "labelBackgroundColor=#FFFFFF;%s" % (FONT, "dashed=1;" if edge.dashed else "")
            )
            if edge.exit:
                style += "exitX=%g;exitY=%g;exitDx=0;exitDy=0;" % edge.exit
            if edge.entry:
                style += "entryX=%g;entryY=%g;entryDx=0;entryDy=0;" % edge.entry
            points = "".join('<mxPoint x="%g" y="%g"/>' % p for p in edge.points)
            geometry = '<mxGeometry relative="1" as="geometry">%s</mxGeometry>' % (
                '<Array as="points">%s</Array>' % points if points else ""
            )
            cells.append(
                '<mxCell id="e%d" value="%s" style="%s" edge="1" parent="1" source="%s" target="%s">'
                "%s</mxCell>"
                % (i, html.escape(edge.label), style, edge.source, edge.target, geometry)
            )
        return (
            '<mxfile host="metawrap2"><diagram id="%s" name="%s"><mxGraphModel grid="1" gridSize="10" '
            'page="0"><root>%s</root></mxGraphModel></diagram></mxfile>\n'
            % (self.name, html.escape(self.title or self.name), "".join(cells))
        )

    def save(self, directory: str) -> str:
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, self.name + ".drawio")
        with open(path, "w") as fh:
            fh.write(self.xml())
        return path


def drawio_binary() -> Optional[str]:
    for candidate in (
        os.environ.get("DRAWIO"),
        "/Applications/draw.io.app/Contents/MacOS/draw.io",
        shutil.which("drawio"),
        shutil.which("draw.io"),
    ):
        if candidate and os.path.exists(candidate):
            return candidate
    return None


def export(drawio_path: str, png_path: str, scale: float = 2.0) -> None:
    """Render with draw.io's own exporter, so the PNG is what draw.io shows."""
    binary = drawio_binary()
    if binary is None:
        raise SystemExit("draw.io desktop is needed to export (brew install --cask drawio)")
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    subprocess.run(
        [binary, "-x", "-f", "png", "-s", str(scale), "-b", "24", "-o", png_path, drawio_path],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    _shrink(png_path)


def _shrink(path: str) -> None:
    """Re-save as a 128-colour palette PNG: a fraction of the size, no visible change.

    The diagrams are flat colours plus anti-aliased edges, so a palette loses nothing a reader can
    see, and it keeps each image under the repository's 512 KB file limit.
    """
    from PIL import Image

    image = Image.open(path).convert("RGB")
    image.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(
        path, optimize=True
    )
