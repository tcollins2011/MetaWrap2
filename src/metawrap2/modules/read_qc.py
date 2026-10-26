"""MetaWrap2 read_qc module: quality-trim reads and remove host/contaminant sequences.

Comprehensive QC of raw paired-end reads in preparation for assembly and downstream work:
optional FastQC report of the input, quality/adapter trimming with Trim Galore, removal of
host sequences with bmtagger, and an optional FastQC report of the cleaned reads. Gzipped
inputs are accepted and decompressed into the output folder first.

The commands each tool runs are in the COMMANDS block below, and the tunable numbers and
defaults in the CONSTANTS block.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shlex
import shutil
from typing import List, NamedTuple, Optional, Tuple

from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..io import reads as _reads
from ..io.seqio import smart_open
from ..scripts import select_human_reads, skip_human_reads
from ._common import (
    absolutize_paths,
    announcement,
    check_disk_space,
    check_fastq,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    expect_produced,
    file_size,
    finish_run,
    fs_move,
    fs_remove,
    make_checkpoint,
    read_layout,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["read_qc"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# --outdir, not -o: FastQC 0.13.0 moved to argparse and dropped -o, failing with
# `fastqc: error: unrecognized arguments: -o ...`. 0.12.x also accepts --outdir. 0.13.0 accepts
# -q, -t and -f fastq (its error named only -o).
FASTQC = "fastqc -q -t {threads} --outdir {outdir} -f fastq {files}"
# No --no_report_file: the trimming report is the only source of Trim Galore's retention numbers,
# which _report_retention reads.
TRIM_GALORE = "trim_galore --paired -j {trim_threads} -o {outdir} {r1} {r2}"
TRIM_GALORE_SINGLE = "trim_galore -j {trim_threads} -o {outdir} {reads}"
FASTP = "fastp -i {r1} -I {r2} -o {o1} -O {o2} -w {trim_threads} -j {json} -h {html}"
FASTP_SINGLE = "fastp -i {reads} -o {o1} -w {trim_threads} -j {json} -h {html}"
# Appended to either fastp command only by --fastp-quality-trim. fastp's defaults do no quality
# trimming: they discard a read when more than 40% of its bases are below Q15 and keep the rest at
# full length. These flags trim low-quality ends instead, as trim-galore does. Opt-in: changing the
# default needs a benchmark measurement, and the current benchmark set has no low-quality data where
# the two behaviours differ.
FASTP_QUALITY_TRIM = "--cut_front --cut_tail --cut_mean_quality 20"
BMTAGGER = "bmtagger.sh -b {bitmask} -x {srprism} -T {tmpdir} -q1" " -1 {r1} -2 {r2} -o {listfile}"
BMTAGGER_SINGLE = (
    "bmtagger.sh -b {bitmask} -x {srprism} -T {tmpdir} -q1" " -1 {reads} -o {listfile}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_HOST = "hg38"  # prefix of the host index inside the bmtagger database folder
# trim_galore's docs say -j above 4 gives no further speedup (cutadapt scaling), and it spawns
# ~4 processes per core, so -t is capped.
TRIM_GALORE_MAX_THREADS = 4
# fastp's docs put the useful thread ceiling at 16; beyond that it is I/O bound.
FASTP_MAX_THREADS = 16
# FastQC's -t is files processed in parallel, and above 4 it warns `[WARNING] You can't use more
# than 4 threads per file - setting this many is a waste`, so -t is capped.
FASTQC_MAX_THREADS = 4
#: Warn when the trimmer keeps less than this share of the input pairs. On the CAMI2 GI benchmark
#: (modern reads) both trimmers keep over 99% of pairs. On the MetaHIT 2010-era 75 bp cohort fastp's
#: defaults kept 41.5% (14,755,926 of 25,242,208 reads discarded whole, none shortened) and
#: trim-galore kept 95%. Below 70% usually means fastp is discarding reads it could trim; on that
#: cohort nonpareil then estimated 20x more sequencing, against 4.06x from the trimmed reads.
LOW_RETENTION_WARN = 0.70
#: Default trimmer. `benchmark run cami2_gi trimmer` over the 10 CAMI2 GI samples (333M read pairs):
#: fastp kept slightly more (+0.03% reads, +0.83% bases, +0.40% Q30 bases) with 5.8% fewer CPU-hours
#: and 29% less wall time, since cutadapt stops scaling at 4 threads. The quality difference is
#: negligible, so cost decided. --trimmer trim-galore is supported and peaks at 0.1 GB RSS against
#: fastp's 1.3 GB. Re-measure with `metawrap2 benchmark score cami2_gi trimmer`.
DEFAULT_TRIMMER = "fastp"
TRIMMERS = ("trim-galore", "fastp")
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 read_qc",
        usage="metawrap2 read_qc [options] -1 reads_1.fastq [-2 reads_2.fastq] -o output_dir",
        description="Quality-trim reads and remove host sequences. Paired files may use any "
        "of the usual mate conventions (_1/_2, _R1/_R2, _R1_001/_R2_001); .gz and "
        ".bz2 are accepted. Single-end and interleaved data are supported with "
        "--single-end / --interleaved.",
    )
    p.add_argument(
        "-1",
        dest="reads_1",
        required=True,
        help="forward fastq reads (or the only file, with --single-end/--interleaved)",
    )
    p.add_argument(
        "-2", dest="reads_2", help="reverse fastq reads (omit for --single-end/--interleaved)"
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "-x",
        "--host",
        default=DEFAULT_HOST,
        help="prefix of host index in bmtagger database folder (default hg38)",
    )
    p.add_argument(
        "--skip-bmtagger", action="store_true", help="do not remove host sequences with bmtagger"
    )
    p.add_argument(
        "--skip-trimming", action="store_true", help="do not trim sequences (skips the trimmer)"
    )
    p.add_argument(
        "--trimmer",
        choices=list(TRIMMERS),
        default=DEFAULT_TRIMMER,
        help="which trimmer to use (default %s). fastp makes one pass and uses up to %d threads; "
        "trim-galore is limited to %d threads by cutadapt. On CAMI2 fastp kept slightly more "
        "bases for less CPU. Both are supported."
        % (DEFAULT_TRIMMER, FASTP_MAX_THREADS, TRIM_GALORE_MAX_THREADS),
    )
    p.add_argument(
        "--fastp-quality-trim",
        action="store_true",
        help="make fastp trim low-quality read ends (%s) instead of only discarding whole "
        "low-quality reads. Use it on older or lower-quality libraries, where fastp's defaults "
        "discard most of the data. Off by default, because the two behaviours differ only on "
        "low-quality reads and the benchmark datasets contain none to measure the change on. "
        "No effect with --trimmer trim-galore, which already trims." % FASTP_QUALITY_TRIM,
    )
    p.add_argument(
        "--skip-pre-qc-report",
        action="store_true",
        help="do not make a FastQC report of the input sequences",
    )
    p.add_argument(
        "--skip-post-qc-report",
        action="store_true",
        help="do not make a FastQC report of the final sequences",
    )
    layout = p.add_mutually_exclusive_group()
    layout.add_argument(
        "--single-end", action="store_true", help="the input is single-end (one file, no mates)"
    )
    layout.add_argument(
        "--interleaved",
        action="store_true",
        help="the input is one file of interleaved paired reads; it is split "
        "into mates first, then processed as paired",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _sample_name(reads_1: str) -> str:
    """Sample name from a *_1.fastq(.gz) file, matching the original ${tmp%_*} behavior."""
    base = os.path.basename(reads_1)
    for token in ("_1.fastq", "_1.fq"):
        if token in base:
            return base.split(token)[0]
    return os.path.splitext(base)[0]


def _trim_stem(path: str) -> str:
    """The stem Trim Galore derives its output filename from: basename minus one extension.

    Trim Galore names its output after the input file, not the sample, so this follows the file it
    was given (which may be a decompressed or de-interleaved intermediate).
    """
    base = os.path.basename(path)
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def _is_ours(path: str, out_dir: str) -> bool:
    """True if *path* is an intermediate this run created inside *out_dir* (safe to move)."""
    out_abs = os.path.abspath(out_dir) + os.sep
    return os.path.abspath(path).startswith(out_abs)


def _decompress(path: str, out_dir: str) -> str:
    """If *path* is gzipped, decompress it into *out_dir* and return the plain path."""
    if not path.endswith(".gz"):
        return path
    comm("Decompressing %s" % path)
    dest = os.path.join(out_dir, os.path.basename(path)[:-3])
    try:
        with smart_open(path, "rb") as fin, open(dest, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    except OSError:
        error("Something went wrong with decompressing %s. Exiting." % path)
    return dest


def _fastqc_report(args, env, name: str, files: List[str]) -> None:
    outdir = os.path.join(args.output, name)
    ensure_dir(outdir)
    run(
        FASTQC.format(
            threads=min(args.threads, FASTQC_MAX_THREADS),
            outdir=outdir,
            files=" ".join(shlex.quote(f) for f in files),
        ),
        env=env,
        tool="fastqc",
    )
    for zipf in glob.glob(os.path.join(outdir, "*.zip")):
        os.remove(zipf)
    comm("%s saved to: %s" % (name, outdir))


def _run_trim_galore(args, env, streams: List[str], paired: bool) -> List[str]:
    """Trim with Trim Galore. Returns the files it produced, in mate order."""
    trim_threads = min(args.threads, TRIM_GALORE_MAX_THREADS)
    if paired:
        run(
            TRIM_GALORE.format(
                outdir=args.output, r1=streams[0], r2=streams[1], trim_threads=trim_threads
            ),
            env=env,
            tool="trim_galore",
        )
        # Trim Galore names paired output <stem>_val_1.fq / <stem>_val_2.fq, and single-end
        # output <stem>_trimmed.fq.
        return [
            os.path.join(args.output, "%s_val_%d.fq" % (_trim_stem(streams[i]), i + 1))
            for i in range(2)
        ]
    run(
        TRIM_GALORE_SINGLE.format(outdir=args.output, reads=streams[0], trim_threads=trim_threads),
        env=env,
        tool="trim_galore",
    )
    return [os.path.join(args.output, "%s_trimmed.fq" % _trim_stem(streams[0]))]


def _run_fastp(args, env, streams: List[str], trimmed: List[str], paired: bool) -> List[str]:
    """Trim with fastp. Returns the files it produced, in mate order.

    fastp writes to the given paths, so it is pointed at the final `trimmed_N.fastq` names and the
    later move is a no-op. Its JSON/HTML report is kept: the JSON has the before/after read counts
    and Q30 fractions that the trimmer benchmark scores.
    """
    trim_threads = min(args.threads, FASTP_MAX_THREADS)
    report_json = os.path.join(args.output, "fastp.json")
    report_html = os.path.join(args.output, "fastp.html")
    # Appended, so without --fastp-quality-trim the command is unchanged.
    extra = " " + FASTP_QUALITY_TRIM if args.fastp_quality_trim else ""
    if paired:
        run(
            FASTP.format(
                r1=streams[0],
                r2=streams[1],
                o1=trimmed[0],
                o2=trimmed[1],
                trim_threads=trim_threads,
                json=report_json,
                html=report_html,
            )
            + extra,
            env=env,
            tool="fastp",
        )
        return [trimmed[0], trimmed[1]]
    run(
        FASTP_SINGLE.format(
            reads=streams[0],
            o1=trimmed[0],
            trim_threads=trim_threads,
            json=report_json,
            html=report_html,
        )
        + extra,
        env=env,
        tool="fastp",
    )
    return [trimmed[0]]


class _Retention(NamedTuple):
    """How much of the input the trimmer kept, in pairs (paired) or reads (single-end).

    ``reasons`` is (what dropped them, how many), in the order the retention line prints them.
    """

    before: int
    after: int
    reasons: Tuple[Tuple[str, int], ...]

    @property
    def kept(self) -> float:
        """Surviving share of the input, 0.0-1.0."""
        return self.after / self.before if self.before else 0.0


#: fastp's ``filtering_result`` keys, with the labels the retention line uses.
#: ``low_quality_reads`` counts reads discarded whole because more than 40% of their bases were
#: under Q15. On the MetaHIT cohort it was 14,755,926, with too_short_reads at 0.
_FASTP_FILTER_REASONS = (
    ("low_quality_reads", "quality filtering"),
    ("too_short_reads", "too short"),
    ("too_many_N_reads", "too many Ns"),
    ("too_long_reads", "too long"),
)
#: Trim Galore has no machine-readable report, so these parse the text one. The wording varies
#: between versions ("Total number of sequences analysed:" vs "...analysed for the length cutoff of
#: 20 bp:"), hence the loose middles. If a pattern stops matching, only the retention line is lost.
_TG_PAIRS_IN = re.compile(r"Total number of sequences analysed[^:\n]*:\s*([\d,]+)")
_TG_PAIRS_DROPPED = re.compile(r"Number of sequence pairs removed[^:\n]*:\s*([\d,]+)")
_TG_READS_IN = re.compile(r"([\d,]+) sequences processed in total")
_TG_READS_DROPPED = re.compile(r"Sequences removed because they became shorter[^:\n]*:\s*([\d,]+)")


def _count(text: str) -> int:
    """A count out of a report, with the thousands separators both tools print stripped."""
    return int(text.replace(",", ""))


def _thousands(value: int) -> str:
    """Format with thousands separators: 12621104 -> "12,621,104"."""
    return f"{value:,}"


def _fastp_retention(path: str, paired: bool) -> Optional[_Retention]:
    """Retention out of fastp's JSON report, or None if it cannot be read.

    The keys are fastp's own: ``summary.before_filtering.total_reads``,
    ``summary.after_filtering.total_reads`` and the ``filtering_result`` block. fastp counts reads,
    so a paired run's numbers are halved here. A missing, truncated or differently keyed report
    returns None instead of raising, so it does not fail a finished QC step.
    """
    import json

    unit = 2 if paired else 1
    try:
        with open(path) as fh:
            doc = json.load(fh)
        summary = doc["summary"]
        before = int(summary["before_filtering"]["total_reads"])
        after = int(summary["after_filtering"]["total_reads"])
        result = doc.get("filtering_result") or {}
        reasons = []
        for key, label in _FASTP_FILTER_REASONS:
            dropped = int(result.get(key) or 0)
            if dropped:
                reasons.append((label, dropped // unit))
    except (AttributeError, KeyError, OSError, TypeError, ValueError):
        return None
    if before <= 0 or after < 0 or after > before:
        return None
    return _Retention(before // unit, after // unit, tuple(reasons))


def _trim_galore_retention(out_dir: str, paired: bool) -> Optional[_Retention]:
    """Retention out of Trim Galore's own ``*_trimming_report.txt``, or None if unreadable.

    All reports in the folder are searched: a paired run writes one per mate, and which one holds
    the pair-level totals has changed between Trim Galore releases.
    """
    text = ""
    for path in sorted(glob.glob(os.path.join(out_dir, "*_trimming_report.txt"))):
        try:
            with open(path, errors="replace") as fh:
                text += fh.read()
        except OSError:
            continue
    came_in = (_TG_PAIRS_IN if paired else _TG_READS_IN).search(text)
    lost = (_TG_PAIRS_DROPPED if paired else _TG_READS_DROPPED).search(text)
    if came_in is None:
        return None
    try:
        before = _count(came_in.group(1))
        dropped = _count(lost.group(1)) if lost else 0
    except ValueError:
        return None
    if before <= 0 or dropped > before:
        return None
    return _Retention(before, before - dropped, (("length filtering", dropped),) if dropped else ())


def _retention_advice(trimmer: str) -> str:
    """Advice for a low-retention run, depending on the trimmer."""
    if trimmer == "fastp":
        return (
            "By default fastp removes a read when more than 40% of its bases are below Q15 and "
            "does no quality trimming, so these reads were discarded, not shortened. On older or "
            "lower-quality data (for example 75 bp reads from around 2010), use --trimmer "
            "trim-galore, which trims low-quality tails and keeps the shorter read, or keep "
            "fastp and add --fastp-quality-trim. A loss this large affects every coverage, "
            "diversity and sequencing-effort estimate computed from these reads."
        )
    return (
        "Check the length cutoff and adapter detection in the *_trimming_report.txt files in the "
        "output folder before using these reads. On older or lower-quality data, short reads "
        "and a strict cutoff can remove most of a library."
    )


def _report_retention(args, paired: bool) -> None:
    """Say how much of the input survived trimming, and warn when it was too little.

    Always printed. Guards against an unreported loss, as when fastp discarded 58.5% of a library.
    The numbers come from the trimmer's own report; if it cannot be read, nothing is printed.
    """
    unit = "pairs" if paired else "reads"
    if args.trimmer == "fastp":
        retention = _fastp_retention(os.path.join(args.output, "fastp.json"), paired)
    else:
        retention = _trim_galore_retention(args.output, paired)
    if retention is None:
        return
    share = 100.0 * retention.kept
    detail = [
        ("%.1f%% removed by %s" if i == 0 else "%.1f%% %s")
        % (100.0 * dropped / retention.before, label)
        for i, (label, dropped) in enumerate(retention.reasons)
    ]
    comm(
        "Read retention: %s %s in -> %s out (%.1f%%)%s"
        % (
            _thousands(retention.before),
            unit,
            _thousands(retention.after),
            share,
            (": " + ", ".join(detail)) if detail else "",
        )
    )
    if retention.kept < LOW_RETENTION_WARN:
        warning(
            "%s kept only %.1f%% of the input %s (%s of %s). %s"
            % (
                args.trimmer,
                share,
                unit,
                _thousands(retention.after),
                _thousands(retention.before),
                _retention_advice(args.trimmer),
            )
        )


def _deinterleave_input(args, sample: str) -> List[str]:
    """Split interleaved input into two mate files inside the output dir; return both paths."""
    announcement("SPLITTING INTERLEAVED READS INTO MATES")
    r1 = os.path.join(args.output, "%s_deinterleaved_1.fastq" % sample)
    r2 = os.path.join(args.output, "%s_deinterleaved_2.fastq" % sample)
    if dry_run():
        comm("(dry run) would split %s into %s and %s" % (args.reads_1, r1, r2))
        return [r1, r2]
    try:
        pairs = _reads.deinterleave(args.reads_1, r1, r2)
    except ValueError as exc:
        error(str(exc))
    comm("split %d read pairs out of %s" % (pairs, args.reads_1))
    return [r1, r2]


def main(argv: List[str]) -> int:
    args = _parse_args(argv)

    if args.reads_2 and args.reads_1 == args.reads_2:
        error("The forward and reverse reads are the same file. Exiting pipeline.")

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    # Every read file must exist, be non-empty and not end mid-record; otherwise the failure shows
    # up several modules later as an empty assembly.
    check_disk_space(
        "read_qc", settings, [p for p in (args.reads_1, args.reads_2) if p], args.output
    )
    validate_inputs(
        "read_qc",
        [(check_fastq, args.reads_1, "forward reads (-1)")]
        + ([(check_fastq, args.reads_2, "reverse reads (-2)")] if args.reads_2 else []),
    )
    layout = read_layout(args)  # validates the layout and reports problems up front
    env = env_for("read_qc", settings)
    rec = start_run("read_qc", args, env, settings, inputs=list(layout.files))

    ckpt = make_checkpoint(args.output)
    try:
        do_bmtagger = not args.skip_bmtagger
        bmtagger_db = settings.db("BMTAGGER_DB")
        if do_bmtagger:
            bitmask = os.path.join(bmtagger_db, args.host + ".bitmask")
            if not (bmtagger_db and os.path.isfile(bitmask)):
                error("%s file doesnt exist. Please configure your bmtagger genome index" % bitmask)

        if os.path.isdir(args.output):
            warning("%s already exists." % args.output)
        ensure_dir(args.output)

        sample = _sample_name(args.reads_1)

        # Interleaved input is split into mates and then handled as paired: Trim Galore has no
        # interleaved mode, and single-end handling would lose the pairing.
        if layout.layout == "interleaved":
            streams = _deinterleave_input(args, sample)
            paired = True
        else:
            streams = [p for p in (args.reads_1, args.reads_2) if p]
            paired = layout.layout == "paired"

        # Decompress the tool inputs and remember the copies so they can be removed at the end
        # (several GB per sample).
        original = list(streams)
        streams = [_decompress(p, args.output) for p in streams]
        decompressed = [p for p, was in zip(streams, original) if p != was]
        if layout.layout == "interleaved":
            decompressed += [p for p in original if _is_ours(p, args.output)]

        if not args.skip_pre_qc_report and ckpt.todo("pre_qc_report"):
            announcement("MAKING PRE-QC REPORT")
            _fastqc_report(args, env, "pre-QC_report", streams)
            ckpt.done("pre_qc_report")

        trimmed = [
            os.path.join(args.output, "trimmed_%d.fastq" % (i + 1)) for i in range(len(streams))
        ]
        if not args.skip_trimming and ckpt.todo("trimming"):
            if args.trimmer == "fastp":
                announcement("RUNNING FASTP")
                produced = _run_fastp(args, env, streams, trimmed, paired)
            else:
                announcement("RUNNING TRIM-GALORE")
                produced = _run_trim_galore(args, env, streams, paired)
            for src, dest in zip(produced, trimmed):
                # fastp already wrote to the final names; only Trim Galore's output needs moving.
                if os.path.abspath(src) != os.path.abspath(dest):
                    fs_move(src, dest)
            if not file_size(trimmed[0]):
                error("Something went wrong with trimming the reads. Exiting.")
            for path in trimmed:
                expect_produced(
                    path, "trimmed reads", hint="Check %s's output in run.stderr." % args.trimmer
                )
            comm("Trimmed reads saved to: %s" % ", ".join(trimmed))
            _report_retention(args, paired)
            streams = list(trimmed)
            ckpt.done("trimming")
        elif all(os.path.isfile(t) for t in trimmed):
            # on --resume, pick up the trimmed reads a previous run left behind
            comm("skipping read trimming (already done; --resume)")
            streams = list(trimmed)

        clean = [
            os.path.join(args.output, "%s_%d.clean.fastq" % (sample, i + 1))
            for i in range(len(streams))
        ]
        if do_bmtagger and ckpt.todo("host_removal"):
            announcement("REMOVING HOST SEQUENCES WITH BMTAGGER")
            tmpdir = os.path.join(args.output, "bmtagger_tmp")
            ensure_dir(tmpdir)
            bitmask = os.path.join(bmtagger_db, args.host + ".bitmask")
            srprism = os.path.join(bmtagger_db, args.host + ".srprism")
            listfile = os.path.join(args.output, "%s.bmtagger.list" % sample)
            comm("running bmtagger with %s %s indexes..." % (bitmask, srprism))
            if paired:
                cmd = BMTAGGER.format(
                    bitmask=bitmask,
                    srprism=srprism,
                    tmpdir=tmpdir,
                    r1=streams[0],
                    r2=streams[1],
                    listfile=listfile,
                )
            else:
                cmd = BMTAGGER_SINGLE.format(
                    bitmask=bitmask,
                    srprism=srprism,
                    tmpdir=tmpdir,
                    reads=streams[0],
                    listfile=listfile,
                )
            run(cmd, env=env, tool="bmtagger.sh")
            if not dry_run() and not (os.path.isfile(listfile) and os.path.getsize(listfile)):
                warning("No host reads found. This is unusual; check the host index.")

            if not dry_run():
                comm("Now sorting out found host reads from the main fastq files...")
                for src, dest in zip(streams, clean):
                    with open(dest, "w") as out:
                        skip_human_reads.filter_reads(listfile, src, out)
                comm("Now writing the host reads to their own files...")
                for i, src in enumerate(streams):
                    host = os.path.join(args.output, "host_reads_%d.fastq" % (i + 1))
                    with open(host, "w") as out:
                        select_human_reads.select_reads(listfile, src, out)
                if not file_size(clean[0]):
                    error(
                        "Something went wrong with removing contaminant reads with bmtagger. "
                        "Exiting."
                    )

            shutil.rmtree(tmpdir, ignore_errors=True)
            fs_remove(listfile)
            streams = list(clean)
            for t in trimmed:
                fs_remove(t)
            ckpt.done("host_removal")
        elif do_bmtagger and all(os.path.isfile(c) for c in clean):
            comm("skipping host read removal (already done; --resume)")
            streams = list(clean)

        final = [
            os.path.join(args.output, "final_pure_reads_%d.fastq" % (i + 1))
            for i in range(len(streams))
        ]
        if ckpt.todo("finalize"):
            # `streams` are intermediates inside args.output, except with --skip-trimming
            # --skip-bmtagger on uncompressed input, where they are the user's raw files. Copy
            # those instead of moving them.
            for src, dest in zip(streams, final):
                if os.path.abspath(src) == os.path.abspath(dest):
                    continue
                if _is_ours(src, args.output):
                    fs_move(src, dest)
                elif dry_run():
                    continue
                else:
                    comm(
                        "Copying input %s instead of moving it, since it is outside the "
                        "output directory" % src
                    )
                    shutil.copyfile(src, dest)
            for path in final:
                expect_produced(path, "final QC'ed reads")
            comm("Contamination-free and trimmed reads are stored in: %s" % ", ".join(final))
            ckpt.done("finalize")

        if not args.skip_post_qc_report and ckpt.todo("post_qc_report"):
            announcement("MAKING POST-QC REPORT")
            _fastqc_report(args, env, "post-QC_report", final)
            ckpt.done("post_qc_report")

        # Remove the decompressed copies of compressed inputs now that the final reads exist.
        for path in decompressed if not dry_run() else []:
            if os.path.isfile(path) and os.path.abspath(path) not in {
                os.path.abspath(f) for f in final
            }:
                comm("removing intermediate file %s" % path)
                os.remove(path)

        announcement("READ QC PIPELINE COMPLETE!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
