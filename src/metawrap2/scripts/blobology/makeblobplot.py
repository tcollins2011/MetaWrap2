#!/usr/bin/env python
# Draw GC-vs-coverage blobplots from a blobplot data table, coloured by taxonomy or bin.
#
# Replaces makeblobplot.R / makeblobplot_with_bins.R / makeblobplot_with_colored_bins.R
# (the three only differed in the base/"catch-all" category label, which is now an
# optional argument). Same inputs and output filename: <input>.<taxlevel>.png
#
# Usage: makeblobplot.py <blobplot_table> <ignore_below_prop> <taxlevel_col> [base_category]
import sys

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Paul Tol qualitative scheme (http://www.sron.nl/~pault/), indexed by number of levels,
# with light grey (#DDDDDD) reserved for the base/catch-all category. Ported from the
# original R scripts so colours stay consistent.
PAULTOL = [
    ["#DDDDDD"],
    ["#DDDDDD", "#4477AA"],
    ["#DDDDDD", "#4477AA", "#CC6677"],
    ["#DDDDDD", "#4477AA", "#DDCC77", "#CC6677"],
    ["#DDDDDD", "#4477AA", "#117733", "#DDCC77", "#CC6677"],
    ["#DDDDDD", "#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"],
    ["#DDDDDD", "#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677", "#AA4499"],
    ["#DDDDDD", "#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677", "#AA4499"],
    ["#DDDDDD", "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#AA4499"],
    ["#DDDDDD", "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"],
    ["#DDDDDD", "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"],
    ["#DDDDDD", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"],
    ["#DDDDDD", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499"],
    ["#DDDDDD", "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499", "#777777"],
]

MAX_COLORS = 13


def main():
    input_file = sys.argv[1]
    ignore_below_prop = float(sys.argv[2])
    taxlevel = sys.argv[3]
    base = sys.argv[4] if len(sys.argv) > 4 else "Not annotated"

    df = pd.read_csv(input_file, sep="\t", header=0)
    df = df[df["len"] >= 200]

    cov_cols = [c for c in df.columns if c.startswith("cov_")]
    if len(cov_cols) > 1:
        df["cov_total"] = df[cov_cols].sum(axis=1)
        cov_cols = [c for c in df.columns if c.startswith("cov_")]

    df[taxlevel] = df[taxlevel].fillna(base).astype(str)

    # melt so each contig appears once per coverage column (read_set)
    m = df.melt(
        id_vars=["seqid", "len", "gc", taxlevel],
        value_vars=cov_cols,
        var_name="read_set",
        value_name="cov",
    )

    # collapse to at most ~13 distinct colours: rarest categories -> "other"
    counts = m[taxlevel].value_counts()
    if len(counts) > 14:
        keep = set(counts.nlargest(MAX_COLORS).index)
        m.loc[~m[taxlevel].isin(keep), taxlevel] = "other"

    # clean.blobs: categories rarer than threshold*total_annotated -> base category
    annotated_total = int((m[taxlevel] != base).sum())
    counts = m[taxlevel].value_counts()
    rare = counts[counts < ignore_below_prop * annotated_total].index
    m.loc[m[taxlevel].isin(rare), taxlevel] = base

    # order: base category first, then by descending frequency
    counts = m[taxlevel].value_counts()
    ordered = [base] + [c for c in counts.index if c != base]
    colors = PAULTOL[min(len(ordered), len(PAULTOL)) - 1]
    color_of = {cat: colors[i % len(colors)] for i, cat in enumerate(ordered)}

    read_sets = cov_cols
    ncols = len(read_sets)
    yticks = [10, 100, 1000, 10000]
    fig, axes = plt.subplots(1, ncols, figsize=(10 * ncols, 13), squeeze=False, sharey=True)
    for i, (ax, read_set) in enumerate(zip(axes[0], read_sets)):
        ax.set_axisbelow(True)
        ax.grid(True, which="major", color="#eaeaea", linewidth=1.0)
        sub = m[m["read_set"] == read_set]
        for cat in ordered:
            pts = sub[sub[taxlevel] == cat]
            if pts.empty:
                continue
            ax.scatter(pts["gc"], pts["cov"], s=14, alpha=1.0 / 3,
                       color=color_of[cat], label=cat, edgecolors="none")
        ax.set_yscale("log")
        ax.set_ylim(5, 10000)
        ax.set_xlim(0.2, 0.8)
        ax.set_xticks([0.2, 0.4, 0.6, 0.8])
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(v) for v in yticks])
        ax.minorticks_off()
        ax.tick_params(labelsize=20)
        for spine in ax.spines.values():
            spine.set_color("#808080")
            spine.set_linewidth(0.7)
        ax.set_xlabel("GC content", fontsize=26)
        if i == 0:
            ax.set_ylabel("Contig abundance", fontsize=26)

    fig.subplots_adjust(left=0.06, right=0.99, top=0.93, bottom=0.16, wspace=0.10)

    # ggplot-style grey facet strips above each panel
    from matplotlib.patches import Rectangle
    strip_h = 0.028
    for ax, read_set in zip(axes[0], read_sets):
        pos = ax.get_position()
        fig.patches.append(Rectangle(
            (pos.x0, pos.y1), pos.width, strip_h, transform=fig.transFigure,
            facecolor="#d9d9d9", edgecolor="#808080", linewidth=0.7, zorder=5, clip_on=False))
        fig.text(pos.x0 + pos.width / 2, pos.y1 + strip_h / 2, read_set,
                 ha="center", va="center", fontsize=20, zorder=6)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(len(labels), 4), markerscale=3,
               frameon=False, title="Annotation", fontsize=18, title_fontsize=20)
    fig.savefig("%s.%s.png" % (input_file, taxlevel), dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
