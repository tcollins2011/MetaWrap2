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
from typing import Dict, List, Optional, Tuple

from . import mimag as _mimag
from .command import run
from .constants import BIN_EXTENSION, STATS_SUFFIX
from .logging import comm
from .utils import MAX_TMPDIR_LEN, short_tmpdir

STATS_HEADER = ["bin", "completeness", "contamination", "GC", "lineage", "N50", "size"]

#: Re-exported so callers do not have to reach into metawrap2.constants or metawrap2.utils for them.
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


def tag_mimag(stats_path: str, evidence_dir: Optional[str] = None) -> List[_mimag.Assessment]:
    """Append MIMAG quality columns to a ``.stats.tsv`` written here. Returns one row's verdict each.

    The checkm-side entry point for :func:`metawrap2.mimag.annotate_stats`, so a caller that already
    imports this module to produce the table does not need a second import to tag it, and so the
    append sits next to :data:`STATS_HEADER` - the contract the new columns have to come after.

    Not called from :func:`run_checkm` or :func:`write_stats`, deliberately. ``write_stats`` output
    is pinned byte for byte by a golden test (seven columns, three rows), and
    :func:`run_checkm_cached` reads a staged table back and re-feeds its rows into ``write_stats``,
    so appending inside either would either break that contract or write a table whose rows are
    three cells wider than its header. Tagging is therefore a pass over the finished file, the same
    shape as :func:`metawrap2.gunc.annotate_stats`, and the caller says when the file is finished.

    *evidence_dir* is an ``annotate`` output directory to look for Bakta/Infernal/tRNAscan-SE output
    in. With none, or with none found, bins that clear 90/5 are tagged
    ``high_quality_unverified_rrna``: nothing has searched them for a 16S, and that is not the same
    as a search that found none.
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


def checkm2_database() -> str:
    """The configured ``CHECKM2DB``, or ``""`` if there is none.

    Resolved here rather than at each call site because CheckM2 is reached from four modules and the
    benchmark assessor, all of which would otherwise have to know about a database that only the
    ``--checkm2`` path uses. Empty is not an error: CheckM2 also keeps a path of its own inside its
    env, set by ``checkm2 database --download``, which is what a hand-installed copy relies on.
    """
    from .config import load_settings

    try:
        return load_settings(strict=False).db("CHECKM2DB")
    except Exception:  # noqa: BLE001 - an unreadable config is the config checker's business
        return ""


def checkm1_data_env() -> Dict[str, str]:
    """``CHECKM_DATA_PATH`` from the configured CHECKM_DB, or nothing if it is unset.

    CheckM1 remembers its reference data in a file *inside its own conda environment*, written by
    ``checkm data setRoot``. That makes a correctly configured database silently stop working the
    moment the environment is rebuilt - the path in ``metawrap2.toml`` is right, the data is on disk,
    and CheckM looks somewhere else and dies in its tree stage with an EOFError from a worker.

    Found exactly that way: consolidating the per-module environments into per-workflow ones replaced
    the environment CheckM1 had been pointed at, and bin_refinement started failing on a machine where
    nothing about the database had changed.

    CheckM1 reads ``CHECKM_DATA_PATH`` in preference to its own stored path, so passing it per call
    makes the configured value the one that counts and removes the environment from the question.
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
) -> str:
    """Run CheckM2 on *bins_dir* and write ``<bins_dir>.stats``. Returns that path.

    Produces the same ``.stats`` contract as :func:`run_checkm`, so every consumer
    (consolidation, dereplication, choose_best_bin, the plots) works unchanged. CheckM2
    lives in its own env (``metawrap2-checkm2``) and needs no pplacer, so it
    has none of CheckM1's per-thread memory wall.

    Its DIAMOND reference is a different database from CheckM1's, and installing one does nothing
    for the other: ``metawrap2 install-db checkm2`` fetches it and writes ``CHECKM2DB``, which is
    passed explicitly below so the run does not depend on what state the env was left in.
    """
    checkm_out = bins_dir + ".checkm"
    db = database if database is not None else checkm2_database()
    cmd = (
        f"checkm2 predict --threads {threads} --input {bins_dir} --output-directory {checkm_out} "
        f"-x {BIN_EXTENSION.lstrip('.')} --force"
    )
    if db:
        cmd += f" --database_path {db}"
    run(
        cmd,
        env=env,
        tool="checkm2",
        log_path=log_path,
        hint="CheckM2 needs its own DIAMOND database, which is not CheckM1's: install it with "
        "'metawrap2 install-db checkm2'.",
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

    It sits beside the bins, which is the point of passing one: CheckM's tree stage writes GBs and
    they belong on the filesystem the caller chose for its output. ``--tmpdir`` makes CheckM set
    ``tempfile.tempdir`` process-wide, and the first thing ``lineage_wf`` does is
    ``mp.Manager().dict()``, so the path has to be short enough to hold that Manager's listener
    socket - see :func:`metawrap2.utils.short_tmpdir`, which is where the shortening lives because
    SemiBin2 has the same constraint for the same reason.
    """
    return short_tmpdir(bins_dir + ".tmp")


def _checkm1(
    bins_dir: str,
    checkm_out: str,
    tmp: str,
    threads: int,
    mem_gb: int,
    env: Optional[str],
    quick: bool,
    log_path: Optional[str],
    pplacer: Optional[int],
    reduced: bool = False,
) -> None:
    """One ``checkm lineage_wf`` invocation. *pplacer* overrides the computed thread count.

    Split out so the retries in :func:`run_checkm` run the same command with one thing changed,
    rather than a second copy of it that could drift from this one. *reduced* adds CheckM's own
    ``--reduced_tree``, which is also what ``quick`` asks for.
    """
    placer = pplacer if pplacer is not None else pplacer_threads(threads, mem_gb)
    cmd = (
        f"checkm lineage_wf -x {BIN_EXTENSION.lstrip('.')} {bins_dir} {checkm_out} "
        f"-t {threads} --tmpdir {tmp} --pplacer_threads {placer}"
    )
    if quick or reduced:
        cmd += " --reduced_tree"
    run(
        cmd,
        env=env,
        tool="checkm",
        log_path=log_path,
        hint="CheckM needs substantial RAM (CheckM1 pplacer ~40GB/thread); try --quick or --checkm2.",
        extra_env=checkm1_data_env(),
        # The retry is the error handling: a non-zero exit here must come back rather than end the
        # run, or the single-thread attempt never happens.
        check=pplacer is not None,
    )


def _clear_checkm_attempt(bins_dir: str, checkm_out: str, tmp: str) -> str:
    """Remove what a crashed CheckM left, and return a fresh tmpdir to retry in.

    Both directories, not just the output: pplacer resumes from what it finds in the tmpdir, so
    leaving a crashed run's half-written tree in place is how a retry reproduces the crash it is
    retrying.
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
    tmp = checkm_tmpdir(bins_dir)
    _checkm1(bins_dir, checkm_out, tmp, threads, mem_gb, env, quick, log_path, pplacer=None)

    stats_path = bins_dir + STATS_SUFFIX
    from .command import runner as _runner

    if _runner.dry_run:
        # Nothing ran, so there is no bin_stats_ext.tsv to summarize. Returning the path lets
        # the caller keep printing the rest of the pipeline under --dry-run.
        return stats_path
    ext = os.path.join(checkm_out, "storage", "bin_stats_ext.tsv")
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        # One retry on a single pplacer thread before giving up. pplacer's parallel tree placement
        # dies as `Dune__exe__Multiprocessing.Child_error` - an OCaml child crash with no cause
        # given - and it is not reproducible: on three CAMI2 skin samples, three of six refinement
        # steps crashed this way at seven pplacer threads and the other three, over the same kind
        # of bin set, finished. CheckM is 100% of what this step costs, so losing an hour of it to
        # a crash that one thread survives is worth the retry being slower.
        comm(
            "CheckM produced no marker statistics (pplacer's parallel placement can die without "
            "reporting a cause); retrying once on a single pplacer thread, which is slower"
        )
        tmp = _clear_checkm_attempt(bins_dir, checkm_out, tmp)
        _checkm1(bins_dir, checkm_out, tmp, threads, mem_gb, env, quick, log_path, pplacer=1)

    if (not os.path.isfile(ext) or os.path.getsize(ext) == 0) and not quick:
        # Second fallback, because the first one is not always the right diagnosis. Isolated here:
        # one 2,022-residue alignment out of 40 killed pplacer on its own, at one thread, three
        # times out of three, with no unusual character in it and nothing degenerate about its
        # length - and the same sequence placed fine against CheckM's *reduced* reference tree.
        # pplacer is an unmaintained OCaml tool whose failure says only
        # `Multiprocessing.Child_error`, so the sequence cannot be fixed from here; the reference
        # tree can be changed.
        #
        # This is CheckM's own documented low-memory mode, not an invention, but it is not free:
        # a reduced tree can choose a different lineage-specific marker set, so completeness and
        # contamination for this bin set are not identical to a full-tree run. Said out loud
        # because a mixed comparison is the kind of thing that has to be visible.
        comm(
            "CheckM still produced nothing on a single pplacer thread. Falling back to "
            "--reduced_tree, which places against CheckM's reduced reference tree: this bin set's "
            "completeness and contamination are then not strictly comparable with a full-tree run, "
            "and the flag is recorded in run_commands.txt"
        )
        tmp = _clear_checkm_attempt(bins_dir, checkm_out, tmp)
        _checkm1(
            bins_dir,
            checkm_out,
            tmp,
            threads,
            mem_gb,
            env,
            quick,
            log_path,
            pplacer=1,
            reduced=True,
        )
    if not os.path.isfile(ext) or os.path.getsize(ext) == 0:
        raise RuntimeError("CheckM did not produce %s" % ext)
    write_stats(summarize(ext, binner=binner), stats_path, binner=binner)
    return stats_path


# --- Making consolidation cheaper -------------------------------------------------------------
#
# bin_refinement scores 2^n-1 hybrid bin sets, and CheckM is the whole cost of it. Two measured
# redundancies in that work, both of which can be removed without changing a single reported number.

#: Sizes outside which a bin cannot be a >=70% complete, <=10% contaminated genome, so there is no
#: point asking CheckM about it.
#:
#: Calibrated against the 17,616 bins in this project's benchmarks that *do* have CheckM2 scores. Of the 4,002 that pass 70/10 the smallest is 0.99 Mb and the largest
#: 6.66 Mb, so both bounds sit clear of the real distribution: together they skip 1.4% of bins and lose
#: **none** of the passing ones.
#:
#: The saving is far larger than that 1.4% suggests, because CheckM's cost follows sequence rather than
#: bin count: of 9,019 raw binner bins, 64 are above 20 Mb and they hold **11.8% of all the sequence**.
#: Those are the bins that are obviously several organisms - metaBAT2's unbinnable remainder reaches
#: 159 Mb - and they are the most expensive thing in the run to gene-call.
#:
#: Where that saving lands is not where it was first claimed. bin_refinement scores bins that
#: :func:`metawrap2.modules.bin_refinement.copy_bins` has already held to 50 kb - 20 Mb, so inside
#: refinement this filter is very nearly a no-op: across 10,168 bins scored there, none is over the
#: ceiling and 56 - 0.02% of the sequence - fall under the floor. It pays in
#: :func:`metawrap2.benchmark.assess.run_bin_quality`, which scores bin sets exactly as a binner emitted
#: them. Keeping it in both places costs nothing and means neither call site depends on the other's
#: filtering, but the honest figure for refinement is "no measurable saving".
MAX_PLAUSIBLE_BIN_BP = 20_000_000
#: The largest known prokaryotic genomes are ~16 Mb, so 20 is a deliberately loose ceiling.
MIN_PLAUSIBLE_BIN_BP = 200_000

# There is deliberately no N50 floor here, though an earlier version of this had one at 2 kb. The
# calibration says it cannot be made to pay: a floor low enough to be safe skips nothing a size bound
# had not already skipped, and one high enough to skip anything starts discarding real genomes. The
# worst N50 among passing bins is 2,025 bp; at a 2,500 bp floor, 36 bins that pass 70/10 are lost, and
# all the floor had bought at 2 kb was 1.4 percentage points more sequence. Fragmentation is a property
# of the assembly, not evidence that a bin is several organisms.


def bin_shape(path: str) -> Tuple[int, int, float]:
    """``(total_bp, n50, gc_fraction)`` for one bin, read once and reused.

    Cheap enough to do for every bin before deciding anything: it is one pass over a file that CheckM
    would otherwise read several times.
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

    ``n50`` is accepted and ignored: it is what a caller naturally has to hand from bin_shape, and the
    calibration found no N50 bound that both saves work and keeps every genome.
    """
    if total_bp > MAX_PLAUSIBLE_BIN_BP:
        return "%.1f Mb, larger than any prokaryotic genome" % (total_bp / 1e6)
    if total_bp < MIN_PLAUSIBLE_BIN_BP:
        return "%.0f kb, too small to be 70%% of a genome" % (total_bp / 1e3)
    return None


def bin_identity(path: str) -> str:
    """A key for "this exact bin", so the same contigs are never scored twice.

    The sorted set of contig ids, hashed. Order is normalised because a bin is a set of contigs and two
    hybrid sets can emit the same one in a different order; sequence is not re-read because the contig
    ids already identify it within a run.

    Measured on a three-binner refinement: of 1,230 bins across the seven hybrid sets, only 840 are
    distinct, so **32% of the CheckM calls re-score a bin already scored**, and some bins appear in all
    seven sets. Nothing about the result changes - the same bin has the same completeness either way.
    """
    from .io.seqio import contig_id, iter_fasta

    ids = sorted(contig_id(header) for header, _seq in iter_fasta(path))
    return hashlib.sha1("\n".join(ids).encode()).hexdigest()


class BinScoreCache:
    """Remembers what CheckM said about a bin, keyed by the bin's contents.

    Lives for one `bin_refinement` run, which is where the duplication is: the hybrid sets it builds
    overlap heavily, and a bin that appears in four of them is the same bin each time. Persisted beside
    the output so a ``--resume`` does not re-pay either.
    """

    def __init__(self, path: Optional[str] = None) -> None:
        self.path = path
        self.rows: Dict[str, List[str]] = {}
        self.hits = self.misses = self.skipped = 0
        if path and os.path.isfile(path):
            try:
                with open(path) as fh:
                    self.rows = {k: list(v) for k, v in json.load(fh).items()}
            except (OSError, ValueError):
                # A truncated cache is a reason to re-measure, not to fail the run.
                self.rows = {}

    def save(self) -> None:
        if not self.path:
            return
        try:
            with open(self.path, "w") as fh:
                json.dump(self.rows, fh)
        except OSError:
            pass  # the cache is an optimisation; failing to persist it must not fail the run

    def summary(self) -> str:
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

    Produces exactly the ``<bins_dir>.stats`` :func:`run_checkm` would, so every consumer is unchanged.
    What differs is how much CheckM is asked to do:

    * a bin whose contigs match one already scored in this run takes that score (see `bin_identity`);
    * a bin that cannot pass 70/10 on its shape alone is not scored at all (see `implausible_reason`),
      and is recorded as failing so nothing downstream can select it;
    * whatever is left is staged into a directory of symlinks and handed to CheckM as usual.

    Falls back to scoring everything if anything about the staging looks wrong - the optimisation is
    never worth a wrong answer, and `run_checkm` is one call away.
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
            cache.skipped += 1
            # Recorded as a definite failure rather than omitted: a bin missing from the stats is a bin
            # downstream code has to guess about, and 0/100 can never win a per-bin comparison.
            skipped.append(
                [name, "0.0", "100.0", "%.2f" % gc, "skipped:%s" % reason, str(n50), str(total_bp)]
            )
            continue
        key = bin_identity(path)
        known = cache.rows.get(key)
        if known is not None:
            cache.hits += 1
            reused.append([name] + list(known[1:]))
            continue
        todo.append(path)

    if not todo:
        log("%s: nothing to score - %s" % (os.path.basename(bins_dir), cache.summary()))
        write_stats(reused + skipped, stats_path, binner=binner)
        return stats_path

    if not reused and not skipped:
        # Nothing to be saved here, so do not pay for a staging directory.
        path = run_checkm(bins_dir, **kwargs)
        _remember(cache, bins, path)
        cache.misses += len(todo)
        cache.save()
        return path

    staged = bins_dir + ".to_score"
    shutil.rmtree(staged, ignore_errors=True)
    os.makedirs(staged, exist_ok=True)
    for path in todo:
        link = os.path.join(staged, os.path.basename(path))
        try:
            os.symlink(os.path.abspath(path), link)
        except OSError:
            shutil.copy2(path, link)
    cache.misses += len(todo)
    log("%s: %s" % (os.path.basename(bins_dir), cache.summary()))

    staged_stats = run_checkm(staged, **kwargs)
    fresh: List[List[str]] = []
    with open(staged_stats) as fh:
        fh.readline()  # header; the columns are STATS_HEADER by contract
        for line in fh:
            if line.strip():
                fresh.append(line.rstrip("\n").split("\t"))
    _remember(cache, todo, staged_stats)
    cache.save()

    write_stats(fresh + reused + skipped, stats_path, binner=binner)
    shutil.rmtree(staged, ignore_errors=True)
    shutil.rmtree(staged + ".checkm", ignore_errors=True)
    return stats_path


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
            cache.rows[bin_identity(path)] = found
