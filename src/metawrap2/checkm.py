"""CheckM helpers shared by the binning and bin_refinement modules.

Runs CheckM (in its conda env) and summarizes its ``bin_stats_ext.tsv`` into the tab-
delimited ``.stats`` file the rest of MetaWrap2 consumes. Ports the old
``summarize_checkm.py`` logic, using ``ast.literal_eval`` instead of ``eval``.
"""

from __future__ import annotations

import ast
import os
from typing import List, Optional

from .command import run
from .constants import BIN_EXTENSION, STATS_SUFFIX

STATS_HEADER = ["bin", "completeness", "contamination", "GC", "lineage", "N50", "size"]

#: Re-exported so callers do not have to reach into metawrap2.constants for it.
__all__ = [
    "STATS_HEADER",
    "STATS_SUFFIX",
    "count_bins",
    "require_bins",
    "run_checkm",
    "run_checkm2",
    "summarize",
    "summarize_checkm2",
    "write_stats",
]


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
        fh.writelines("\t".join(row) + "\n" for row in rows_sorted)


def count_bins(bins_dir: str, extension: str = BIN_EXTENSION) -> int:
    """Number of bin fastas in *bins_dir*."""
    if not os.path.isdir(bins_dir):
        return 0
    return len([f for f in os.listdir(bins_dir) if f.endswith(extension)])


def require_bins(bins_dir: str, context: str = "") -> None:
    """Fail with an actionable message if *bins_dir* holds no bins.

    CheckM's own response to an empty directory is "No bins found. Check the extension (-x)
    used to identify bins", which sends users looking for a file-extension problem when the
    real cause is usually that no bin passed their -c/-x quality thresholds.
    """
    if count_bins(bins_dir) > 0:
        return
    message = ["No bins to analyse in %s." % bins_dir]
    if context:
        message.append(context)
    raise RuntimeError(" ".join(message))


def pplacer_threads(threads: int, mem_gb: int) -> int:
    """pplacer needs ~40GB/thread; use min(threads, mem/40), at least 1."""
    return max(1, min(threads, mem_gb // 40))


def summarize_checkm2(quality_report_tsv: str, binner: Optional[str] = None) -> List[List[str]]:
    """Parse CheckM2's ``quality_report.tsv`` into the same rows as :func:`summarize`.

    CheckM2 is a single classifier rather than a lineage-specific marker set, so there is no
    marker lineage to report; the lineage column is filled with ``checkm2`` to make it
    obvious in the output which estimator produced the numbers.
    """
    rows: List[List[str]] = []
    with open(quality_report_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = {
                name: header.index(name)
                for name in (
                    "Name",
                    "Completeness",
                    "Contamination",
                    "Genome_Size",
                    "GC_Content",
                    "Contig_N50",
                )
            }
        except ValueError as exc:
            raise RuntimeError(
                "Unexpected CheckM2 quality_report.tsv columns in %s: %s"
                % (quality_report_tsv, header)
            ) from exc
        for line in fh:
            if not line.strip():
                continue
            cut = line.rstrip("\n").split("\t")
            rows.append(
                [
                    cut[idx["Name"]],
                    str(cut[idx["Completeness"]])[:5],
                    str(cut[idx["Contamination"]])[:5],
                    str(cut[idx["GC_Content"]])[:5],
                    "checkm2",
                    cut[idx["Contig_N50"]],
                    cut[idx["Genome_Size"]],
                ]
                + ([binner] if binner is not None else [])
            )
    return rows


def run_checkm2(
    bins_dir: str,
    *,
    threads: int,
    env: Optional[str],
    binner: Optional[str] = None,
    log_path: Optional[str] = None,
) -> str:
    """Run CheckM2 on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    Produces the same ``.stats`` contract as :func:`run_checkm`, so every consumer
    (consolidation, dereplication, choose_best_bin, the plots) works unchanged. CheckM2
    lives in its own env (``metawrap2-bin_refinement-checkm2``) and needs no pplacer, so it
    has none of CheckM1's per-thread memory wall.
    """
    checkm_out = bins_dir + ".checkm"
    cmd = (
        f"checkm2 predict --threads {threads} --input {bins_dir} --output-directory {checkm_out} "
        f"-x {BIN_EXTENSION.lstrip('.')} --force"
    )
    run(
        cmd,
        env=env,
        tool="checkm2",
        log_path=log_path,
        hint="CheckM2 needs its DIAMOND database; install it with 'checkm2 database --download' "
        "and make sure CHECKM2DB points at it.",
    )

    stats_path = bins_dir + STATS_SUFFIX
    from .command import runner as _runner

    if _runner.dry_run:
        return stats_path
    report = os.path.join(checkm_out, "quality_report.tsv")
    if not os.path.isfile(report) or os.path.getsize(report) == 0:
        raise RuntimeError("CheckM2 did not produce %s" % report)
    write_stats(summarize_checkm2(report, binner=binner), stats_path, binner=binner)
    return stats_path


def run_checkm(
    bins_dir: str,
    *,
    threads: int,
    mem_gb: int,
    env: Optional[str],
    binner: Optional[str] = None,
    quick: bool = False,
    log_path: Optional[str] = None,
    checkm2: bool = False,
    checkm2_env: Optional[str] = None,
) -> str:
    """Run CheckM on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    CheckM1 lives in the bin_refinement conda env, so callers pass that env even from other
    modules. With ``checkm2=True`` this delegates to :func:`run_checkm2` instead, which runs
    in its own env; the ``.stats`` output is identical in shape either way.
    """
    if checkm2:
        return run_checkm2(
            bins_dir, threads=threads, env=checkm2_env or env, binner=binner, log_path=log_path
        )
    checkm_out = bins_dir + ".checkm"
    tmp = bins_dir + ".tmp"
    os.makedirs(tmp, exist_ok=True)  # harmless under --dry-run; checkm needs it to exist
    cmd = f"checkm lineage_wf -x {BIN_EXTENSION.lstrip('.')} {bins_dir} {checkm_out} -t {threads} --tmpdir {tmp} --pplacer_threads {pplacer_threads(threads, mem_gb)}"
    if quick:
        cmd += " --reduced_tree"
    run(
        cmd,
        env=env,
        tool="checkm",
        log_path=log_path,
        hint="CheckM needs substantial RAM (CheckM1 pplacer ~40GB/thread); try --quick or --checkm2.",
    )

    stats_path = bins_dir + STATS_SUFFIX
    from .command import runner as _runner

    if _runner.dry_run:
        # Nothing ran, so there is no bin_stats_ext.tsv to summarize. Returning the path lets
        # the caller keep printing the rest of the pipeline under --dry-run.
        return stats_path
    ext = os.path.join(checkm_out, "storage", "bin_stats_ext.tsv")
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        raise RuntimeError("CheckM did not produce %s" % ext)
    write_stats(summarize(ext, binner=binner), stats_path, binner=binner)
    return stats_path
