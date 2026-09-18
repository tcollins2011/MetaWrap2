"""CheckM helpers shared by the binning and bin_refinement modules.

Runs CheckM (in its conda env) and summarizes its ``bin_stats_ext.tsv`` into the tab-
delimited ``.stats`` file the rest of MetaWrap2 consumes. Ports the old
``summarize_checkm.py`` logic, using ``ast.literal_eval`` instead of ``eval``.
"""

from __future__ import annotations

import ast
import os
from typing import Dict, List, Optional

from .command import run

STATS_HEADER = ["bin", "completeness", "contamination", "GC", "lineage", "N50", "size"]


def summarize(bin_stats_ext_tsv: str, binner: Optional[str] = None) -> List[List[str]]:
    """Parse CheckM's bin_stats_ext.tsv into rows matching the MetaWrap2 .stats format.

    Value formatting mirrors the original tool (numeric fields truncated to 5 chars).
    """
    rows: List[List[str]] = []
    with open(bin_stats_ext_tsv) as fh:
        for line in fh:
            name, _, payload = line.partition("\t")
            d = ast.literal_eval(payload.strip())
            lineage = d["marker lineage"]
            if "__" in lineage:
                lineage = lineage.split("__")[1]
            row = [
                name,
                str(d["Completeness"])[:5],
                str(d["Contamination"])[:5],
                str(d["GC"])[:5],
                lineage,
                str(d["N50 (contigs)"]),
                str(d["Genome size"]),
            ]
            if binner is not None:
                row.append(binner)
            rows.append(row)
    return rows


def write_stats(rows: List[List[str]], out_path: str, binner: Optional[str] = None) -> None:
    """Write a .stats file: header, then rows sorted by completeness (descending)."""
    header = list(STATS_HEADER) + (["binner"] if binner is not None else [])
    rows_sorted = sorted(rows, key=lambda r: float(r[1]), reverse=True)
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows_sorted:
            fh.write("\t".join(row) + "\n")


def pplacer_threads(threads: int, mem_gb: int) -> int:
    """pplacer needs ~40GB/thread; use min(threads, mem/40), at least 1."""
    return max(1, min(threads, mem_gb // 40))


def run_checkm(
    bins_dir: str,
    *,
    threads: int,
    mem_gb: int,
    env: Optional[str],
    binner: Optional[str] = None,
    quick: bool = False,
    log_path: Optional[str] = None,
) -> str:
    """Run CheckM lineage_wf on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    CheckM lives in the bin_refinement conda env, so callers pass that env even from other
    modules.
    """
    checkm_out = bins_dir + ".checkm"
    tmp = bins_dir + ".tmp"
    os.makedirs(tmp, exist_ok=True)
    cmd = (
        "checkm lineage_wf -x fa {bins} {out} -t {t} --tmpdir {tmp} --pplacer_threads {p}"
    ).format(bins=bins_dir, out=checkm_out, t=threads, tmp=tmp,
             p=pplacer_threads(threads, mem_gb))
    if quick:
        cmd += " --reduced_tree"
    run(cmd, env=env, tool="checkm", log_path=log_path,
        hint="CheckM needs substantial RAM (CheckM1 pplacer ~40GB/thread); try --quick or --checkm2.")

    ext = os.path.join(checkm_out, "storage", "bin_stats_ext.tsv")
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        raise RuntimeError("CheckM did not produce %s" % ext)
    stats_path = bins_dir + ".stats"
    write_stats(summarize(ext, binner=binner), stats_path, binner=binner)
    return stats_path
