"""A static PNG sunburst of the functional hierarchy, drawn from the Krona input table.

``ktImportText`` writes only an interactive HTML file drawn by JavaScript, so there is no image
for reports or figures. This draws the same input file with matplotlib (already a dependency),
one PNG per Krona HTML.

* BRITE's non-pathway branches ("Brite Hierarchies", "Not Included in Pathway or Brite") are
  kept, since dropping them would inflate the other percentages. The pathway-only wheel is a
  separate file, as for Krona.
* Level A is drawn in a fixed order with a fixed colour per branch, so wheels from different
  samples can be compared side by side.
* Small wedges are drawn without labels, and the caption gives their count and total share.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional, Sequence, Tuple

#: KEGG BRITE's top-level branches in KEGG's order, each with a fixed colour used in every wheel.
#: Other names (another hierarchy, a later BRITE revision) use the overflow colours below.
LEVEL_A_COLOURS: Dict[str, str] = {
    "Metabolism": "#1f6fb4",
    "Genetic Information Processing": "#e07b39",
    "Environmental Information Processing": "#4c9f70",
    "Cellular Processes": "#b4436c",
    "Organismal Systems": "#7f63b8",
    "Human Diseases": "#a0522d",
    "Drug Development": "#c9a227",
    "Brite Hierarchies": "#8a8f98",
    "Not Included in Pathway or Brite": "#b9bec6",
}

#: Colours for branch names not in the palette above.
OVERFLOW_COLOURS: Tuple[str, ...] = ("#5b8ff9", "#f6bd16", "#5ad8a6", "#945fb9", "#e8684a")

#: Wedges below this share of the total are drawn without a label (2% of 3,000 genes is 60).
LABEL_THRESHOLD = 0.02

#: Label threshold for the outer (pathway) ring. BRITE has 570 pathways and a genome touches a
#: few hundred, so labels at 2% still overlap. On one CAMI2 bin, 228 of 251 pathway wedges were
#: below 2%.
OUTER_LABEL_THRESHOLD = 0.03

#: Number of rings, as in BRITE: branch -> chemical process -> pathway.
MAX_DEPTH = 3


def read_krona_text(path: str) -> List[Tuple[float, Tuple[str, ...]]]:
    """Read ``ktImportText``'s input: a weight, then one column per level.

    This is the file :func:`metawrap2.hierarchy.write_krona_input` writes. Rows whose weight does
    not parse are skipped.
    """
    rows: List[Tuple[float, Tuple[str, ...]]] = []
    with open(path) as fh:
        for fields in csv.reader(fh, delimiter="\t"):
            if not fields:
                continue
            try:
                weight = float(fields[0])
            except (TypeError, ValueError):
                continue
            if weight <= 0:
                continue
            labels = tuple(field.strip() for field in fields[1 : MAX_DEPTH + 1] if field.strip())
            if labels:
                rows.append((weight, labels))
    return rows


def collapse_ring(
    rows: Sequence[Tuple[float, Tuple[str, ...]]], depth: int
) -> List[Tuple[Tuple[str, ...], float]]:
    """Return the total weight per path prefix of length *depth*, in parent order.

    Sorted by parent so each child wedge is drawn inside its parent; sorting each ring by weight
    would misplace children.
    """
    totals: Dict[Tuple[str, ...], float] = {}
    for weight, labels in rows:
        if len(labels) < depth:
            # A path shorter than this ring is padded with blank labels, so every ring sums to
            # the same total.
            key = labels + ("",) * (depth - len(labels))
        else:
            key = labels[:depth]
        totals[key] = totals.get(key, 0.0) + weight
    return sorted(totals.items(), key=lambda item: (_branch_rank(item[0][0]), item[0]))


def _branch_rank(name: str) -> int:
    """Fixed position for a level-A branch, so ring order is stable across wheels."""
    order = list(LEVEL_A_COLOURS)
    return order.index(name) if name in order else len(order)


def colour_for(name: str, used: Optional[Dict[str, str]] = None) -> str:
    """Return a stable colour for a level-A branch name."""
    if name in LEVEL_A_COLOURS:
        return LEVEL_A_COLOURS[name]
    table = used if used is not None else {}
    if name not in table:
        table[name] = OVERFLOW_COLOURS[len(table) % len(OVERFLOW_COLOURS)]
    return table[name]


def _shade(colour: str, ring: int) -> str:
    """Lighten a branch's colour for each ring outwards."""
    colour = colour.lstrip("#")
    rgb = [int(colour[i : i + 2], 16) for i in (0, 2, 4)]
    fraction = 0.0 if ring == 0 else 0.22 * ring
    lightened = [int(value + (255 - value) * fraction) for value in rgb]
    return "#%02x%02x%02x" % tuple(lightened)


def draw(
    rows: Sequence[Tuple[float, Tuple[str, ...]]],
    out_path: str,
    title: str = "",
    unit_label: str = "weight",
) -> Optional[str]:
    """Draw the wheel to *out_path*. Returns the path written, or None for an empty distribution.

    An empty distribution (a genome with no KO) is a valid result, so this does not raise.
    """
    rows = [(weight, labels) for weight, labels in rows if weight > 0 and labels]
    if not rows:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    total = sum(weight for weight, _ in rows)
    figure, axes = plt.subplots(figsize=(8.4, 7.4))
    overflow: Dict[str, str] = {}
    legend_handles: List[object] = []
    outer_unlabelled = 0
    outer_unlabelled_weight = 0.0

    for ring in range(MAX_DEPTH):
        wedges = collapse_ring(rows, ring + 1)
        sizes = [weight for _, weight in wedges]
        colours = [_shade(colour_for(key[0], overflow), ring) for key, _ in wedges]
        labels = []
        for key, weight in wedges:
            name = key[-1]
            share = weight / total
            # The innermost ring's names would overlap at the centre, so they go in the legend.
            if ring == 0 or not name or share < _threshold(ring):
                if ring == MAX_DEPTH - 1 and name:
                    outer_unlabelled += 1
                    outer_unlabelled_weight += weight
                labels.append("")
            else:
                labels.append("%s\n%.1f%%" % (_wrap(name), 100.0 * share))
        axes.pie(
            sizes,
            radius=0.45 + 0.185 * (ring + 1),
            colors=colours,
            labels=labels,
            labeldistance=0.8,
            startangle=90,
            counterclock=False,
            textprops={"fontsize": 7.0},
            wedgeprops={"width": 0.185, "edgecolor": "white", "linewidth": 0.6},
        )
        if ring == 0:
            from matplotlib.patches import Patch

            legend_handles = [
                Patch(
                    facecolor=colour_for(key[0], overflow),
                    label="%s  %.1f%%" % (key[0], 100.0 * weight / total),
                )
                for key, weight in wedges
            ]

    axes.set(aspect="equal")
    if title:
        axes.set_title(title, fontsize=11)
    if legend_handles:
        axes.legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(1.0, 0.5),
            fontsize=7.5,
            frameon=False,
            title="branch (inner ring)",
            title_fontsize=8,
        )
    note = "%s %s; rings: branch, chemical process, pathway" % (_number(total), unit_label)
    if outer_unlabelled:
        # Counted on the outer ring only; summing across rings would count weight repeatedly.
        note += "\n%d pathway wedge(s) under %.0f%% are drawn unlabelled (%.1f%% of the total)" % (
            outer_unlabelled,
            100.0 * _threshold(MAX_DEPTH - 1),
            100.0 * outer_unlabelled_weight / total,
        )
    figure.text(0.5, 0.02, note, ha="center", fontsize=7.5, color="#4a5058")
    directory = os.path.dirname(os.path.abspath(out_path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    figure.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return out_path


def draw_from_krona_text(path: str, out_path: str = "", title: str = "") -> Optional[str]:
    """Draw the wheel for a Krona input file, defaulting to the same name with ``.png``."""
    out_path = out_path or os.path.splitext(path)[0] + ".png"
    return draw(read_krona_text(path), out_path, title=title)


def _threshold(ring: int) -> float:
    """Return the share a wedge needs to be labelled; higher on the outer ring."""
    return OUTER_LABEL_THRESHOLD if ring == MAX_DEPTH - 1 else LABEL_THRESHOLD


def _wrap(name: str, width: int = 18) -> str:
    """Break a long branch name over two lines; KEGG's are up to 44 characters."""
    if len(name) <= width:
        return name
    words = name.split()
    first: List[str] = []
    while words and len(" ".join(first + [words[0]])) <= width:
        first.append(words.pop(0))
    if not first:
        return name
    return " ".join(first) + "\n" + " ".join(words)


def _number(value: float) -> str:
    return "%.0f" % value if value >= 10 else "%.2f" % value
