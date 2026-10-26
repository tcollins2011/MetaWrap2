"""CheckM helpers shared by the binning and bin_refinement modules.

Runs CheckM (in its conda env) and summarizes its ``bin_stats_ext.tsv`` into the tab-
delimited ``.stats`` file the rest of MetaWrap2 consumes. Ports the old
``summarize_checkm.py`` logic, using ``ast.literal_eval`` instead of ``eval``.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import shutil
import threading
from typing import Dict, List, Optional, Tuple

from . import mimag as _mimag
from .command import run
from .constants import BIN_EXTENSION, STATS_SUFFIX
from .logging import comm
from .utils import MAX_TMPDIR_LEN, short_tmpdir

STATS_HEADER = ["bin", "completeness", "contamination", "GC", "lineage", "N50", "size"]

#: Includes names re-exported from metawrap2.constants and metawrap2.utils.
__all__ = [
    "MAX_TMPDIR_LEN",
    "STATS_HEADER",
    "STATS_SUFFIX",
    "checkm2_database",
    "checkm_tmpdir",
    "count_bins",
    "require_bins",
    "run_checkm",
    "run_checkm2",
    "stats_scorer",
    "summarize",
    "summarize_checkm2",
    "tag_mimag",
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


def stats_scorer(stats_path: str) -> Optional[str]:
    """``checkm2`` or ``checkm1``: which CheckM wrote the ``.stats.tsv`` at *stats_path*.

    CheckM2 rows carry ``checkm2`` in the lineage column (see :func:`summarize_checkm2`); CheckM1
    rows carry a lineage. None if the file is missing or has no scored row.
    """
    if not os.path.isfile(stats_path):
        return None
    with open(stats_path) as fh:
        fh.readline()
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if len(cut) > 4 and not cut[4].startswith("skipped:"):
                return "checkm2" if cut[4] == "checkm2" else "checkm1"
    return None


def tag_mimag(stats_path: str, evidence_dir: Optional[str] = None) -> List[_mimag.Assessment]:
    """Append MIMAG quality columns to a ``.stats.tsv`` written here. Returns one verdict per row.

    Wraps :func:`metawrap2.mimag.annotate_stats` so callers of this module need no second import.

    Not called from :func:`run_checkm` or :func:`write_stats`: ``write_stats`` output is pinned by
    a golden test, and :func:`run_checkm_cached` feeds staged rows back into ``write_stats``, so
    appending there would break the test or make rows wider than the header. Like
    :func:`metawrap2.gunc.annotate_stats`, this is a pass over the finished file.

    *evidence_dir* is an ``annotate`` output directory to search for Bakta, Infernal and
    tRNAscan-SE output. If none is given or found, bins passing 90/5 are tagged
    ``high_quality_unverified_rrna``, because they were not searched for a 16S.
    """
    names: List[str] = []
    if os.path.isfile(stats_path):
        with open(stats_path) as fh:
            fh.readline()  # header; column 0 is the bin name by STATS_HEADER contract
            names = [line.split("\t")[0] for line in fh if line.strip()]
    evidence = _mimag.collect_evidence(evidence_dir, names) if evidence_dir else {}
    return _mimag.annotate_stats(stats_path, evidence)


def count_bins(bins_dir: str, extension: str = BIN_EXTENSION) -> int:
    """Number of bin fastas in *bins_dir*."""
    if not os.path.isdir(bins_dir):
        return 0
    return len([f for f in os.listdir(bins_dir) if f.endswith(extension)])


def require_bins(bins_dir: str, context: str = "") -> None:
    """Fail with an actionable message if *bins_dir* holds no bins.

    CheckM's own message for an empty directory points at the file extension (-x), while the usual
    cause is that no bin passed the -c/-x quality thresholds.
    """
    if count_bins(bins_dir) > 0:
        return
    message = ["No bins to analyse in %s." % bins_dir]
    if context:
        message.append(context)
    raise RuntimeError(" ".join(message))


#: pplacer gains little past 8 threads, and its placements are the same at any thread count.
PPLACER_USEFUL_THREADS = 8


#: CheckM1's memory guidance in GB: ~40 for the full genome tree, under 16 with ``--reduced_tree``
#: (``checkm lineage_wf -h``).
FULL_TREE_GB = 40
REDUCED_TREE_GB = 16


def pplacer_threads(threads: int, mem_gb: int, reduced_tree: bool = False) -> int:
    """Threads for pplacer: min(threads, mem/40, 8), at least 1 (mem/16 on the reduced tree).

    The memory term follows CheckM's guidance of ~40 GB per pplacer thread on the full tree and under
    16 GB on the reduced tree. pplacer's memory is per process (~34 GB for the full tree), so
    `bin_refinement` scores bin sets in parallel instead of giving one call more threads.
    """
    per_thread = REDUCED_TREE_GB if reduced_tree else FULL_TREE_GB
    return max(1, min(threads, mem_gb // per_thread, PPLACER_USEFUL_THREADS))


def summarize_checkm2(quality_report_tsv: str, binner: Optional[str] = None) -> List[List[str]]:
    """Parse CheckM2's ``quality_report.tsv`` into the same rows as :func:`summarize`.

    CheckM2 has no lineage-specific marker set, so the lineage column is set to ``checkm2`` to show
    which estimator produced the numbers.
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


def checkm2_database() -> str:
    """The configured ``CHECKM2DB``, or ``""`` if there is none.

    Resolved here for the four modules and the benchmark assessor that call CheckM2. Empty is not an
    error: CheckM2 also stores a path in its env, set by ``checkm2 database --download``, which a
    manually installed copy uses.
    """
    from .config import load_settings

    try:
        return load_settings(strict=False).db("CHECKM2DB")
    except Exception:  # noqa: BLE001 - an unreadable config is the config checker's business
        return ""


def checkm1_data_env() -> Dict[str, str]:
    """``CHECKM_DATA_PATH`` from the configured CHECKM_DB, or nothing if it is unset.

    CheckM1 stores its data path in a file inside its conda environment (``checkm data setRoot``).
    When the environment is rebuilt that path is lost, and CheckM fails in its tree stage with an
    EOFError from a worker even though ``metawrap2.toml`` and the data are correct. This happened
    when per-module environments were merged into per-workflow ones.

    CheckM1 prefers ``CHECKM_DATA_PATH`` over its stored path, so passing it on each call makes the
    configured value take effect.
    """
    from .config import load_settings

    try:
        path = load_settings().db("CHECKM_DB")
    except Exception:  # noqa: BLE001 - a missing or broken config must not stop a tool from running
        return {}
    return {"CHECKM_DATA_PATH": path} if path else {}


def run_checkm2(
    bins_dir: str,
    *,
    threads: int,
    env: Optional[str],
    binner: Optional[str] = None,
    log_path: Optional[str] = None,
    database: Optional[str] = None,
    ttable: Optional[int] = None,
) -> str:
    """Run CheckM2 on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    Writes the same ``.stats`` format as :func:`run_checkm`, so consolidation, dereplication,
    choose_best_bin and the plots work unchanged. CheckM2 runs in its own env
    (``metawrap2-checkm2``) and does not use pplacer, so it avoids CheckM1's per-thread memory cost.

    Its DIAMOND database is separate from CheckM1's. ``metawrap2 install-db checkm2`` fetches it and
    writes ``CHECKM2DB``, which is passed explicitly so the run does not depend on env state.

    *ttable* (4 or 11) is passed to ``--ttable``; by default CheckM2 calls genes under both tables
    and keeps the one with the higher coding density.
    """
    checkm_out = bins_dir + ".checkm"
    db = database if database is not None else checkm2_database()
    cmd = (
        f"checkm2 predict --threads {threads} --input {bins_dir} --output-directory {checkm_out} "
        f"-x {BIN_EXTENSION.lstrip('.')} --force"
    )
    if db:
        cmd += f" --database_path {db}"
    if ttable:
        cmd += f" --ttable {ttable}"
    run(
        cmd,
        env=env,
        tool="checkm2",
        log_path=log_path,
        hint="CheckM2 needs its own DIAMOND database, separate from CheckM1's data. Install it "
        "with 'metawrap2 install-db checkm2'.",
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


def checkm_tmpdir(bins_dir: str) -> str:
    """Create CheckM1's temp directory for *bins_dir* and return a path short enough to pass it.

    It is placed beside the bins, because CheckM's tree stage writes GBs that belong on the output
    filesystem. ``--tmpdir`` sets ``tempfile.tempdir`` process-wide and ``lineage_wf`` starts with
    ``mp.Manager().dict()``, so the path must be short enough for the Manager's listener socket
    (see :func:`metawrap2.utils.short_tmpdir`, shared with SemiBin2).
    """
    return short_tmpdir(bins_dir + ".tmp")


def _checkm1(
    bins_dir: str,
    checkm_out: str,
    tmp: str,
    threads: int,
    mem_gb: int,
    env: Optional[str],
    reduced: bool,
    log_path: Optional[str],
    pplacer: Optional[int],
) -> None:
    """One ``checkm lineage_wf`` invocation. *pplacer* overrides the computed thread count.

    Shared by the retries in :func:`run_checkm`. *reduced* adds ``--reduced_tree``.
    """
    placer = pplacer if pplacer is not None else pplacer_threads(threads, mem_gb, reduced)
    cmd = (
        f"checkm lineage_wf -x {BIN_EXTENSION.lstrip('.')} {bins_dir} {checkm_out} "
        f"-t {threads} --tmpdir {tmp} --pplacer_threads {placer}"
    )
    if reduced:
        cmd += " --reduced_tree"
    run(
        cmd,
        env=env,
        tool="checkm",
        log_path=log_path,
        hint="CheckM1 needs about %d GB of RAM per pplacer thread on the %s genome tree; check -m, "
        "or drop --checkm1 to score with CheckM2."
        % ((REDUCED_TREE_GB, "reduced") if reduced else (FULL_TREE_GB, "full")),
        extra_env=checkm1_data_env(),
        # On the first attempt a non-zero exit must return so the single-thread retry can run.
        check=pplacer is not None,
    )


def _clear_checkm_attempt(bins_dir: str, checkm_out: str, tmp: str) -> str:
    """Remove what a crashed CheckM left, and return a fresh tmpdir to retry in.

    Clears the tmpdir as well as the output, because pplacer resumes from the tmpdir and a
    half-written tree would make the retry crash again.
    """
    shutil.rmtree(checkm_out, ignore_errors=True)
    shutil.rmtree(tmp, ignore_errors=True)
    return checkm_tmpdir(bins_dir)


def run_checkm(
    bins_dir: str,
    *,
    threads: int,
    mem_gb: int,
    env: Optional[str],
    binner: Optional[str] = None,
    reduced_tree: bool = False,
    log_path: Optional[str] = None,
    checkm2: bool = False,
    checkm2_env: Optional[str] = None,
    ttable: Optional[int] = None,
) -> str:
    """Run CheckM on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    CheckM1 is in the bin_refinement conda env, so callers from other modules pass that env. With
    ``checkm2=True`` this calls :func:`run_checkm2`, which uses its own env; the ``.stats`` format
    is the same. *reduced_tree* (CheckM1's ``--reduced_tree``) applies to CheckM1 only, *ttable* to
    CheckM2 only.
    """
    if checkm2:
        return run_checkm2(
            bins_dir,
            threads=threads,
            env=checkm2_env or env,
            binner=binner,
            log_path=log_path,
            ttable=ttable,
        )
    checkm_out = bins_dir + ".checkm"
    tmp = checkm_tmpdir(bins_dir)
    _checkm1(bins_dir, checkm_out, tmp, threads, mem_gb, env, reduced_tree, log_path, pplacer=None)

    stats_path = bins_dir + STATS_SUFFIX
    from .command import runner as _runner

    if _runner.dry_run:
        # Nothing ran; return the path so --dry-run can continue printing the pipeline.
        return stats_path
    ext = os.path.join(checkm_out, "storage", "bin_stats_ext.tsv")
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        # Retry once on a single pplacer thread. Parallel placement intermittently fails with
        # `Dune__exe__Multiprocessing.Child_error` (an OCaml child crash with no cause): on three
        # CAMI2 skin samples, three of six refinement steps failed this way at seven pplacer
        # threads. The single-thread retry is slower but avoids rerunning the whole step.
        comm(
            "CheckM produced no marker statistics (parallel pplacer placement can fail without "
            "an error message); retrying once on a single pplacer thread (slower)"
        )
        tmp = _clear_checkm_attempt(bins_dir, checkm_out, tmp)
        _checkm1(bins_dir, checkm_out, tmp, threads, mem_gb, env, reduced_tree, log_path, pplacer=1)

    if (not os.path.isfile(ext) or os.path.getsize(ext) == 0) and not reduced_tree:
        # Second fallback: CheckM's documented --reduced_tree mode. In one case a single
        # 2,022-residue alignment out of 40 crashed pplacer at one thread in three of three runs,
        # with nothing unusual in the sequence, and placed without error against the reduced tree.
        # pplacer is unmaintained and reports only `Multiprocessing.Child_error`.
        #
        # A reduced tree can choose a different marker set, so completeness and contamination may
        # differ from a full-tree run; the message below says so.
        comm(
            "CheckM still produced nothing on a single pplacer thread. Retrying with "
            "--reduced_tree (CheckM's reduced reference tree). Completeness and contamination "
            "for this bin set may differ from a full-tree run; the flag is recorded in "
            "run_commands.txt"
        )
        tmp = _clear_checkm_attempt(bins_dir, checkm_out, tmp)
        _checkm1(
            bins_dir,
            checkm_out,
            tmp,
            threads,
            mem_gb,
            env,
            True,
            log_path,
            pplacer=1,
        )
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        raise RuntimeError("CheckM did not produce %s" % ext)
    write_stats(summarize(ext, binner=binner), stats_path, binner=binner)
    return stats_path


# --- Making consolidation cheaper -------------------------------------------------------------
#
# bin_refinement scores 2^n-1 hybrid bin sets, and CheckM is nearly all of the cost. Two measured
# redundancies are removed below without changing any reported number.

#: Size bounds outside which a bin cannot be a >=70% complete, <=10% contaminated genome, so it is
#: not sent to CheckM.
#:
#: Calibrated on the 17,616 benchmark bins with CheckM2 scores. The 4,002 that pass 70/10 range
#: from 0.99 Mb to 6.66 Mb; the bounds skip 1.4% of bins and none of the passing ones.
#:
#: CheckM's cost follows sequence length, not bin count: of 9,019 raw binner bins, 64 are above
#: 20 Mb and hold 11.8% of all sequence (metaBAT2's unbinned remainder reaches 159 Mb).
#:
#: Inside bin_refinement the filter has almost no effect, because
#: :func:`metawrap2.modules.bin_refinement.copy_bins` already limits bins to 50 kb - 20 Mb: of
#: 10,168 bins scored there, none is over the ceiling and 56 (0.02% of sequence) are under the
#: floor. The saving is in :func:`metawrap2.benchmark.assess.run_bin_quality`, which scores bin
#: sets as the binners wrote them. Both call sites apply it so neither depends on the other.
MAX_PLAUSIBLE_BIN_BP = 20_000_000
#: The largest known prokaryotic genomes are ~16 Mb, so the 20 Mb ceiling is loose.
MIN_PLAUSIBLE_BIN_BP = 200_000

# No N50 floor (an earlier version used 2 kb). In calibration, a floor low enough to keep every
# passing bin skips nothing the size bounds do not, and a higher one loses real genomes: the lowest
# N50 among passing bins is 2,025 bp, and a 2,500 bp floor loses 36 bins that pass 70/10. The 2 kb
# floor saved only 1.4 percentage points of sequence. Fragmentation reflects the assembly, not a
# mixed bin.


def bin_shape(path: str) -> Tuple[int, int, float]:
    """``(total_bp, n50, gc_fraction)`` for one bin, read once and reused.

    One pass over the file, cheap enough to run on every bin.
    """
    from .io.seqio import iter_fasta

    lengths: List[int] = []
    gc = at = 0
    for _header, seq in iter_fasta(path):
        lengths.append(len(seq))
        upper = seq.upper()
        gc += upper.count("G") + upper.count("C")
        at += upper.count("A") + upper.count("T")
    total = sum(lengths)
    half, running, n50 = total / 2.0, 0, 0
    for length in sorted(lengths, reverse=True):
        running += length
        n50 = length
        if running >= half:
            break
    return total, n50, (gc / (gc + at)) if (gc + at) else 0.0


def implausible_reason(total_bp: int, n50: int = 0) -> Optional[str]:
    """Why this bin cannot pass 70/10, or None if it might. See MAX_PLAUSIBLE_BIN_BP.

    ``n50`` is accepted for convenience (callers have it from bin_shape) and ignored, since no N50
    bound both saved work and kept every genome.
    """
    if total_bp > MAX_PLAUSIBLE_BIN_BP:
        return "%.1f Mb, larger than any prokaryotic genome" % (total_bp / 1e6)
    if total_bp < MIN_PLAUSIBLE_BIN_BP:
        return "%.0f kb, too small to be 70%% of a genome" % (total_bp / 1e3)
    return None


def bin_identity(path: str) -> str:
    """A content key for a bin, so the same set of contigs is scored only once.

    A hash of the sorted contig ids. Sorting handles hybrid sets that list the same contigs in a
    different order; contig ids identify the sequence within a run.

    On a three-binner refinement, 840 of the 1,230 bins across the seven hybrid sets were distinct,
    so 32% of CheckM calls re-scored a bin already scored.
    """
    from .io.seqio import contig_id, iter_fasta

    ids = sorted(contig_id(header) for header, _seq in iter_fasta(path))
    return hashlib.sha1("\n".join(ids).encode()).hexdigest()


class BinScoreCache:
    """Remembers what CheckM said about a bin, keyed by the bin's contents.

    Lives for one `bin_refinement` run, whose hybrid sets overlap heavily. Saved beside the output
    so ``--resume`` can reuse it. Guarded by a lock because bin sets can be scored in parallel.
    """

    def __init__(self, path: Optional[str] = None) -> None:
        self.path = path
        self.rows: Dict[str, List[str]] = {}
        self.lock = threading.RLock()
        self.hits = self.misses = self.skipped = 0
        if path and os.path.isfile(path):
            try:
                with open(path) as fh:
                    self.rows = {k: list(v) for k, v in json.load(fh).items()}
            except (OSError, ValueError):
                # A truncated cache is discarded and the bins are re-scored.
                self.rows = {}

    def save(self) -> None:
        if not self.path:
            return
        try:
            # The whole save holds the lock, so a save from an older snapshot cannot replace a
            # newer one. The temp file keeps a reader (the next --resume) from seeing a partial file.
            with self.lock:
                payload = json.dumps(dict(self.rows))
                tmp = "%s.%d.tmp" % (self.path, threading.get_ident())
                with open(tmp, "w") as fh:
                    fh.write(payload)
                os.replace(tmp, self.path)
        except OSError:
            pass  # the cache is an optimisation; failing to persist it must not fail the run

    def summary(self) -> str:
        with self.lock:
            return self._summary()

    def _summary(self) -> str:
        total = self.hits + self.misses + self.skipped
        if not total:
            return "no bins scored"
        return (
            "%d bin(s): %d scored, %d reused from an identical bin, %d skipped as implausible"
            % (
                total,
                self.misses,
                self.hits,
                self.skipped,
            )
        )


def run_checkm_cached(
    bins_dir: str,
    *,
    cache: Optional[BinScoreCache] = None,
    prefilter: bool = True,
    say=None,
    **kwargs,
) -> str:
    """:func:`run_checkm`, without re-scoring bins already scored or scoring hopeless ones.

    Writes the same ``<bins_dir>.stats`` as :func:`run_checkm`. Differences in what CheckM scores:

    * a bin whose contigs match one already scored in this run reuses that score (`bin_identity`);
    * a bin that cannot pass 70/10 on size is not scored (`implausible_reason`) and is recorded as
      failing so nothing downstream selects it;
    * the remaining bins are staged as symlinks in a directory and passed to CheckM.

    Calls :func:`run_checkm` on the whole directory when no FASTA bins are found or when no bin is
    reused or skipped.
    """
    from .io.seqio import fasta_files

    log = say or (lambda _msg: None)
    cache = cache if cache is not None else BinScoreCache()
    binner = kwargs.get("binner")
    stats_path = bins_dir + STATS_SUFFIX

    bins = fasta_files(bins_dir)
    if not bins:
        return run_checkm(bins_dir, **kwargs)

    reused: List[List[str]] = []
    skipped: List[List[str]] = []
    todo: List[str] = []
    for path in bins:
        name = os.path.splitext(os.path.basename(path))[0]
        total_bp, n50, gc = bin_shape(path)
        reason = implausible_reason(total_bp, n50) if prefilter else None
        if reason:
            with cache.lock:
                cache.skipped += 1
            # Recorded as 0/100 instead of omitted, so downstream code sees every bin and this one
            # cannot win a per-bin comparison.
            skipped.append(
                [name, "0.0", "100.0", "%.2f" % gc, "skipped:%s" % reason, str(n50), str(total_bp)]
            )
            continue
        key = bin_identity(path)
        with cache.lock:
            known = cache.rows.get(key)
            if known is not None:
                cache.hits += 1
        if known is not None:
            reused.append([name] + list(known[1:]))
            continue
        todo.append(path)

    if not todo:
        log("%s: nothing to score - %s" % (os.path.basename(bins_dir), cache.summary()))
        write_stats(reused + skipped, stats_path, binner=binner)
        return stats_path

    if not reused and not skipped:
        # Nothing reused or skipped, so no staging directory is needed.
        path = run_checkm(bins_dir, **kwargs)
        _remember(cache, bins, path)
        with cache.lock:
            cache.misses += len(todo)
        cache.save()
        return path

    with cache.lock:
        cache.misses += len(todo)
    log("%s: %s" % (os.path.basename(bins_dir), cache.summary()))

    fresh = _score_staged(bins_dir + ".to_score", todo, cache, kwargs)
    cache.save()

    write_stats(fresh + reused + skipped, stats_path, binner=binner)
    return stats_path


def plan_scoring(
    bin_dirs: List[str], cache: BinScoreCache, prefilter: bool = True
) -> Dict[str, List[str]]:
    """Which bins each bin set has to score, decided in set order before any scoring starts.

    A bin shared by several sets is scored by the first set that has it, as in a sequential run.
    Deciding this up front keeps the groups each CheckM call sees the same when sets are scored in
    parallel. This matters for CheckM2, whose DIAMOND hits for a protein depend on the rest of the
    query file. Bins already in the cache, and bins the size filter rules out, are left out.
    """
    from .io.seqio import fasta_files

    with cache.lock:
        claimed = set(cache.rows)
    plan: Dict[str, List[str]] = {}
    for directory in bin_dirs:
        todo: List[str] = []
        for path in fasta_files(directory):
            total_bp, n50, _gc = bin_shape(path)
            if prefilter and implausible_reason(total_bp, n50):
                continue
            key = bin_identity(path)
            if key in claimed:
                continue
            claimed.add(key)
            todo.append(path)
        plan[directory] = todo
    return plan


def score_bins(label: str, bins: List[str], cache: BinScoreCache, **kwargs) -> None:
    """Score *bins* with one CheckM call and store their rows in *cache*.

    Writes no ``.stats``; each bin set's file is written from the cache once all sets are scored.
    *label* is the bin set directory, and the bins are staged beside it.
    """
    if not bins:
        return
    _score_staged(label.rstrip("/") + ".to_score", bins, cache, kwargs)
    with cache.lock:
        cache.misses += len(bins)
    cache.save()


def _score_staged(staged: str, bins: List[str], cache: BinScoreCache, kwargs) -> List[List[str]]:
    """Score *bins* as one CheckM call over a directory of symlinks, and return their rows."""
    shutil.rmtree(staged, ignore_errors=True)
    os.makedirs(staged, exist_ok=True)
    for path in bins:
        link = os.path.join(staged, os.path.basename(path))
        try:
            os.symlink(os.path.abspath(path), link)
        except OSError:
            shutil.copy2(path, link)
    stats = run_checkm(staged, **kwargs)
    rows: List[List[str]] = []
    if os.path.isfile(stats):
        with open(stats) as fh:
            fh.readline()  # header; the columns are STATS_HEADER by contract
            for line in fh:
                if line.strip():
                    rows.append(line.rstrip("\n").split("\t"))
    _remember(cache, bins, stats)
    shutil.rmtree(staged, ignore_errors=True)
    shutil.rmtree(staged + ".checkm", ignore_errors=True)
    return rows


def _remember(cache: BinScoreCache, bins: List[str], stats_path: str) -> None:
    """Store each scored bin's row under its content key, for the next set that contains it."""
    if not os.path.isfile(stats_path):
        return
    by_name: Dict[str, List[str]] = {}
    with open(stats_path) as fh:
        fh.readline()
        for line in fh:
            if line.strip():
                row = line.rstrip("\n").split("\t")
                by_name[row[0]] = row
    for path in bins:
        name = os.path.splitext(os.path.basename(path))[0]
        found = by_name.get(name)
        if found is not None:
            with cache.lock:
                cache.rows[bin_identity(path)] = found
