#!/usr/bin/env python
"""Render MetaWrap2 flowchart PNGs from their draw.io (.xml) sources, in pure Python.

Each ``<name>.xml`` in this folder is a draw.io diagram (compressed mxGraph). This script
decodes it and redraws the boxes, labels, and arrows with matplotlib into ``<name>.png`` -
no draw.io needed. Run ``python flowcharts/make_flowcharts.py`` to regenerate all figures.
"""

from __future__ import annotations

import base64
import glob
import html
import os
import re
import textwrap
import urllib.parse
import xml.etree.ElementTree as ET
import zlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch


def decode_diagram(path: str) -> str:
    """Return the decompressed mxGraphModel XML from a draw.io .xml file."""
    with open(path) as fh:
        text = fh.read()
    m = re.search(r"<diagram[^>]*>(.*?)</diagram>", text, re.DOTALL)
    if not m:  # already-plain mxGraphModel
        return text
    raw = base64.b64decode(m.group(1))
    return urllib.parse.unquote(zlib.decompress(raw, -15).decode("utf-8"))


def _clean_label(value: str) -> str:
    value = value.replace("<br>", "\n").replace("<br/>", "\n")
    value = re.sub(r"<[^>]+>", "", value)  # strip HTML tags
    value = value.replace("&nbsp;", " ")
    value = html.unescape(value)
    return value.strip()


def _style_dict(style: str) -> dict:
    out = {}
    for part in (style or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def parse_cells(xml: str):
    """Return (nodes, edges). nodes: id->dict(x,y,w,h,label,fill,rounded,ellipse)."""
    root = ET.fromstring(xml)
    nodes, edges = {}, []
    for cell in root.iter("mxCell"):
        cid = cell.get("id")
        if cell.get("edge") == "1":
            edges.append((cell.get("source"), cell.get("target")))
            continue
        if cell.get("vertex") != "1":
            continue
        geo = cell.find("mxGeometry")
        if geo is None:
            continue
        raw_style = cell.get("style", "") or ""
        style = _style_dict(raw_style)
        is_ellipse = "ellipse" in raw_style.split(";")
        nodes[cid] = {
            "x": float(geo.get("x", 0)),
            "y": float(geo.get("y", 0)),
            "w": float(geo.get("width", 80)),
            "h": float(geo.get("height", 40)),
            "label": _clean_label(cell.get("value", "") or ""),
            "fill": style.get("fillColor", "none"),
            "stroke": style.get("strokeColor", "#000000"),
            "rounded": style.get("rounded") == "1",
            "ellipse": is_ellipse,
        }
    return nodes, edges


def _wrap(label: str, width_px: float, fontsize: int = 10) -> str:
    """Wrap each explicit line to roughly fit the box width (draw.io auto-wraps text)."""
    chars = max(6, int(width_px / (0.78 * fontsize)))
    out = []
    for line in label.split("\n"):
        out.extend(textwrap.wrap(line, chars) or [""])
    return "\n".join(out)


def _draw_orthogonal(ax, src, tgt) -> None:
    """Draw a right-angle (Z) connector from src to tgt with an arrowhead, like draw.io."""
    scx, scy = src["x"] + src["w"] / 2, src["y"] + src["h"] / 2
    tcx, tcy = tgt["x"] + tgt["w"] / 2, tgt["y"] + tgt["h"] / 2
    dx, dy = tcx - scx, tcy - scy
    if abs(dy) >= abs(dx):  # vertical flow: exit top/bottom, enter opposite
        y0 = src["y"] + src["h"] if dy > 0 else src["y"]
        y1 = tgt["y"] if dy > 0 else tgt["y"] + tgt["h"]
        mid = (y0 + y1) / 2
        pts = [(scx, y0), (scx, mid), (tcx, mid), (tcx, y1)]
    else:  # horizontal flow: exit left/right, enter opposite
        x0 = src["x"] + src["w"] if dx > 0 else src["x"]
        x1 = tgt["x"] if dx > 0 else tgt["x"] + tgt["w"]
        mid = (x0 + x1) / 2
        pts = [(x0, scy), (mid, scy), (mid, tcy), (x1, tcy)]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color="#555555", lw=1.2, zorder=2, solid_capstyle="projecting")
    ax.add_patch(
        FancyArrowPatch(
            pts[-2],
            pts[-1],
            arrowstyle="-|>",
            mutation_scale=13,
            color="#555555",
            lw=1.2,
            shrinkA=0,
            shrinkB=0,
            zorder=2,
        )
    )


def render(path: str) -> str:
    nodes, edges = parse_cells(decode_diagram(path))
    if not nodes:
        return ""
    minx = min(n["x"] for n in nodes.values())
    miny = min(n["y"] for n in nodes.values())
    maxx = max(n["x"] + n["w"] for n in nodes.values())
    maxy = max(n["y"] + n["h"] for n in nodes.values())
    pad = 20
    W, H = (maxx - minx) + 2 * pad, (maxy - miny) + 2 * pad

    fig, ax = plt.subplots(figsize=(W / 100.0, H / 100.0), dpi=150)
    ax.set_xlim(minx - pad, maxx + pad)
    ax.set_ylim(miny - pad, maxy + pad)
    ax.invert_yaxis()  # draw.io y grows downward
    ax.set_aspect("equal")
    ax.axis("off")

    def center(nid):
        n = nodes[nid]
        return n["x"] + n["w"] / 2, n["y"] + n["h"] / 2

    # edges first (behind boxes), routed orthogonally like draw.io
    for src, tgt in edges:
        if src not in nodes or tgt not in nodes:
            continue
        _draw_orthogonal(ax, nodes[src], nodes[tgt])

    # nodes
    for n in nodes.values():
        fill = "none" if n["fill"] in ("none", "") else n["fill"]
        edge = n["stroke"] if n["stroke"] != "none" else "#333333"
        cx, cy = n["x"] + n["w"] / 2, n["y"] + n["h"] / 2
        if n["ellipse"]:
            ax.add_patch(
                Ellipse(
                    (cx, cy),
                    n["w"],
                    n["h"],
                    facecolor=fill,
                    edgecolor=edge,
                    linewidth=1.3,
                    zorder=3,
                )
            )
        else:
            style = "round,pad=0,rounding_size=8" if n["rounded"] else "square,pad=0"
            ax.add_patch(
                FancyBboxPatch(
                    (n["x"], n["y"]),
                    n["w"],
                    n["h"],
                    boxstyle=style,
                    facecolor=fill,
                    edgecolor=edge,
                    linewidth=1.3,
                    zorder=3,
                    mutation_aspect=1,
                )
            )
        if n["label"]:
            ax.text(
                cx, cy, _wrap(n["label"], n["w"]), ha="center", va="center", fontsize=10, zorder=4
            )

    out = os.path.splitext(path)[0] + ".png"
    fig.savefig(out, dpi=150, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    return out


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    for xmlf in sorted(glob.glob(os.path.join(here, "*.xml"))):
        out = render(xmlf)
        if out:
            print("wrote", os.path.basename(out))


if __name__ == "__main__":
    main()
