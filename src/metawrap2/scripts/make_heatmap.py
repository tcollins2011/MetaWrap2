#!/usr/bin/env python
"""Draw a clustered heatmap of genome (bin) abundances across samples with Seaborn.

Reads a tab-separated abundance table (bins x samples, as produced by
``split_salmon_out_into_bins``), standardizes each column to counts-per-million, log-scales,
and writes a clustered heatmap PNG.

Ported from metaWRAP's ``make_heatmap.py``. Kept as a standalone script because it pulls in
numpy/pandas/matplotlib/seaborn, which live in the module's conda environment; invoke it via
``python -m metawrap2.scripts.make_heatmap <table> <png>``.
"""

from __future__ import annotations

import sys
from typing import List


def load_data(filename: str):
    import pandas as pd

    print("loading abundance data...")
    df = pd.read_csv(filename, sep="\t", index_col=0)
    # remove all-zero rows
    df = df[(df.T != 0).any()]
    # standardize columns by total sum in each column, then to counts-per-million
    df = df.div(df.sum(axis=0), axis=1)
    df = 1000000 * df
    return df


def draw_clustermap(df, lut) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("drawing clustermap...")
    sns.set(font_scale=1)
    df = df.fillna(0)
    if lut is not False:
        g = sns.clustermap(df, figsize=(14, 8), col_colors=lut, col_cluster=True,
                           yticklabels=True, cmap="magma")
    else:
        g = sns.clustermap(df, figsize=(14, 8), col_cluster=True, yticklabels=True,
                           cmap="magma")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)


def main(argv: List[str]) -> int:
    print("loading libs...")
    import matplotlib
    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rc("font", family="arial")

    table = argv[0]
    df = load_data(table)

    # log standardize
    df += 0.01
    df = np.log(df)

    draw_clustermap(df, lut=False)

    out = argv[1] if len(argv) > 1 else ".".join(table.split(".")[:-1]) + ".png"
    plt.savefig(out, bbox_inches="tight", dpi=300)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
