"""A static wheel of the functional hierarchy, drawn from the same table Krona is given.

``ktImportText`` writes one thing: an interactive HTML file whose wheel is drawn by JavaScript in
the browser. Nothing on disk is an image. So a run's functional hierarchy could not appear in a
report, a release note, a paper figure or a Slack message - the only way to look at it was to open a
browser, and the only way to rasterise it would be to drive a headless browser, which is not a
dependency worth taking for a picture.

This draws the same data with matplotlib, which is already a hard dependency (the bundled plotting
helpers need it). One PNG per Krona HTML, from the identical input file, so the two can never
disagree. The interactive wheel stays: it is better for exploring, this is better for showing.

Three deliberate choices, each of which is a way this plot could lie:

* **Scaffolding branches are kept.** BRITE's non-pathway branches - "Brite Hierarchies", "Not
  Included in Pathway or Brite" - are the largest in the file, and dropping them to make a prettier
  wheel would inflate every remaining percentage. The pathway-only wheel exists for readability and
  is written as a separate file, exactly as Krona's is.
* **Wedges are not reordered between plots.** Level A is drawn in a fixed order with a fixed colour
  per branch name, so two samples' wheels can be put side by side and the same colour means the
  same branch. Sorting each wheel by its own weights would make two wheels of the same cohort
  incomparable at a glance, which is the main thing anybody wants to do with them.
* **Small wedges are drawn but not labelled**, and the caption says how many were left unlabelled
  and what share they hold. A wheel with 570 pathway labels is unreadable, and one that silently
  omits the long tail misrepresents the distribution.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, Optional, Sequence, Tuple

#: KEGG BRITE's eight top-level branches, in KEGG's own order, each with a fixed colour. Fixed so
#: that the same branch is the same colour in every wheel this project draws - a reader comparing
#: two samples should not have to re-read a legend. Anything outside this list (a different
#: hierarchy, a future BRITE revision) falls through to the overflow cycle below.
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

#: For a branch name the palette above does not know.
OVERFLOW_COLOURS: Tuple[str, ...] = ("#5b8ff9", "#f6bd16", "#5ad8a6", "#945fb9", "#e8684a")

#: A wedge smaller than this share of the total is drawn without a label. 2% of a 3,000-gene
#: distribution is 60 genes: worth drawing, too small to label legibly at this figure size.
LABEL_THRESHOLD = 0.02

#: The outer ring needs a higher bar than the middle one. BRITE has 570 pathways and a real genome
#: touches a few hundred of them, so labelling everything above 2% there still produces overlapping
#: text; 3% keeps the ring readable. Measured on one CAMI2 bin: 228 of 251 pathway wedges fall below
#: 2%, and they hold the long tail rather than the signal.
OUTER_LABEL_THRESHOLD = 0.03

#: Rings, outermost last. Three is what BRITE gives: branch -> chemical process -> pathway.
MAX_DEPTH = 3


def read_krona_text(path: str) -> List[Tuple[float, Tuple[str, ...]]]:
    """Read ``ktImportText``'s input: a weight, then one column per level.

    The same file :func:`metawrap2.hierarchy.write_krona_input` writes, read back rather than
    recomputed, so the picture cannot drift from the wheel beside it. Rows whose weight does not
    parse are skipped rather than treated as zero: a zero is a statement that a branch is absent,
    and a malformed line is not.
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
    """Total weight per distinct path prefix of length *depth*, in parent order.

    Parent order, not weight order: a child wedge has to sit inside the parent it belongs to, so the
    rings are built by walking the level-A order once and nesting beneath each branch. Sorting each
    ring independently by weight would draw children under the wrong parents - the classic way a
    hand-rolled sunburst ends up decorative rather than true.
    """
    totals: Dict[Tuple[str, ...], float] = {}
    for weight, labels in rows:
        if len(labels) < depth:
            # A path that stops short of this ring: counted at the depth it reaches, so the ring
            # totals still sum to the whole. Without this, every ring below the shallowest path
            # would silently total less than the one above it.
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
    """A stable colour per level-A branch name."""
    if name in LEVEL_A_COLOURS:
        return LEVEL_A_COLOURS[name]
    table = used if used is not None else {}
    if name not in table:
        table[name] = OVERFLOW_COLOURS[len(table) % len(OVERFLOW_COLOURS)]
    return table[name]


def _shade(colour: str, ring: int) -> str:
    """Lighten a branch's colour for each ring outwards, so nesting reads without a legend."""
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

    Returns rather than raises on an empty distribution: a genome with no KO at all is a real
    result, and the module that calls this should report it as such rather than treat a blank wheel
    as a failure.
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
            # The innermost ring is never labelled in place: all eight branch names would be drawn
            # at the centre on top of each other. They go in the legend, which also fixes the
            # colours in one place for a reader comparing two wheels.
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
        # Counted on the outer ring only. Summing unlabelled wedges across rings would count the
        # same weight up to three times and report a share above 100%.
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
    """The share a wedge needs before it is labelled, which is higher on the crowded outer ring."""
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
