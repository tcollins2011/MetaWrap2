"""MetaWrap2 read_qc module: quality-trim reads and remove host/contaminant sequences.

Comprehensive QC of raw paired-end reads in preparation for assembly and downstream work:
optional FastQC report of the input, quality/adapter trimming with Trim Galore, removal of
host sequences with bmtagger, and an optional FastQC report of the cleaned reads. Gzipped
inputs are accepted and transparently decompressed into the output folder first.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers/defaults in the CONSTANTS block. Change a flag there
and it takes effect — no need to follow the orchestration logic underneath.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import glob
import os
import shlex
import shutil
from typing import List

from ..config import load_settings
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

CONDA_ENV = "metawrap2-read_qc"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
FASTQC = "fastqc -q -t {threads} -o {outdir} -f fastq {files}"
TRIM_GALORE = "trim_galore --no_report_file --paired -j {trim_threads} -o {outdir} {r1} {r2}"
TRIM_GALORE_SINGLE = "trim_galore --no_report_file -j {trim_threads} -o {outdir} {reads}"
BMTAGGER = "bmtagger.sh -b {bitmask} -x {srprism} -T {tmpdir} -q1" " -1 {r1} -2 {r2} -o {listfile}"
BMTAGGER_SINGLE = (
    "bmtagger.sh -b {bitmask} -x {srprism} -T {tmpdir} -q1" " -1 {reads} -o {listfile}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_HOST = "hg38"  # prefix of the host index inside the bmtagger database folder
# trim_galore's own docs warn that -j above 4 gives no further speedup (cutadapt scaling),
# and it spawns ~4 processes per core, so we cap it rather than passing -t straight through.
TRIM_GALORE_MAX_THREADS = 4
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
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
        "--skip-bmtagger", action="store_true", help="dont remove host sequences with bmtagger"
    )
    p.add_argument(
        "--skip-trimming", action="store_true", help="dont trim sequences with trimgalore"
    )
    p.add_argument(
        "--skip-pre-qc-report",
        action="store_true",
        help="dont make FastQC report of input sequences",
    )
    p.add_argument(
        "--skip-post-qc-report",
        action="store_true",
        help="dont make FastQC report of final sequences",
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

    Trim Galore names its output after the *input file*, not after the sample, so this has to
    follow whatever we actually handed it (which may be a decompressed or de-interleaved
    intermediate rather than the user's original file).
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
            threads=args.threads, outdir=outdir, files=" ".join(shlex.quote(f) for f in files)
        ),
        env=env,
        tool="fastqc",
    )
    for zipf in glob.glob(os.path.join(outdir, "*.zip")):
        os.remove(zipf)
    comm("%s saved to: %s" % (name, outdir))


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
    # Every read file must exist, be non-empty, and not end mid-record: an empty or truncated
    # FASTQ otherwise produces an empty assembly, zero bins, and a failure several modules later.
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

        # Interleaved input is split into mates first and then follows the ordinary paired
        # path: Trim Galore has no interleaved mode, and treating it as single-end would throw
        # away the pairing the assembler and binners rely on.
        if layout.layout == "interleaved":
            streams = _deinterleave_input(args, sample)
            paired = True
        else:
            streams = [p for p in (args.reads_1, args.reads_2) if p]
            paired = layout.layout == "paired"

        # Decompress whatever we are about to feed the tools, remembering the copies so they
        # can be removed at the end (several GB per sample of pure duplication otherwise).
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
            announcement("RUNNING TRIM-GALORE")
            trim_threads = min(args.threads, TRIM_GALORE_MAX_THREADS)
            if paired:
                run(
                    TRIM_GALORE.format(
                        outdir=args.output, r1=streams[0], r2=streams[1], trim_threads=trim_threads
                    ),
                    env=env,
                    tool="trim_galore",
                )
                # Trim Galore names paired output <stem>_val_1.fq / <stem>_val_2.fq, and
                # single-end output <stem>_trimmed.fq.
                produced = [
                    os.path.join(args.output, "%s_val_%d.fq" % (_trim_stem(streams[i]), i + 1))
                    for i in range(2)
                ]
            else:
                run(
                    TRIM_GALORE_SINGLE.format(
                        outdir=args.output, reads=streams[0], trim_threads=trim_threads
                    ),
                    env=env,
                    tool="trim_galore",
                )
                produced = [os.path.join(args.output, "%s_trimmed.fq" % _trim_stem(streams[0]))]
            for src, dest in zip(produced, trimmed):
                fs_move(src, dest)
            if not file_size(trimmed[0]):
                error("Something went wrong with trimming the reads. Exiting.")
            for path in trimmed:
                expect_produced(
                    path, "trimmed reads", hint="Check Trim Galore's output in run.stderr."
                )
            comm("Trimmed reads saved to: %s" % ", ".join(trimmed))
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
                warning("No contamination reads found, which is very unlikely.")

            if not dry_run():
                comm("Now sorting out found host reads from the main fastq files...")
                for src, dest in zip(streams, clean):
                    with open(dest, "w") as out:
                        skip_human_reads.filter_reads(listfile, src, out)
                comm("Now keeping the host reads in their own files... for science...")
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
            # `streams` are intermediates we created inside args.output in every case except
            # "--skip-trimming --skip-bmtagger with uncompressed input", where they are still
            # the user's raw input. Moving those would destroy the user's data, so copy when
            # the source is not ours to move.
            for src, dest in zip(streams, final):
                if os.path.abspath(src) == os.path.abspath(dest):
                    continue
                if _is_ours(src, args.output):
                    fs_move(src, dest)
                elif dry_run():
                    continue
                else:
                    comm("Copying (not moving) input %s - it is outside the output directory" % src)
                    shutil.copyfile(src, dest)
            for path in final:
                expect_produced(path, "final QC'ed reads")
            comm("Contamination-free and trimmed reads are stored in: %s" % ", ".join(final))
            ckpt.done("finalize")

        if not args.skip_post_qc_report and ckpt.todo("post_qc_report"):
            announcement("MAKING POST-QC REPORT")
            _fastqc_report(args, env, "post-QC_report", final)
            ckpt.done("post_qc_report")

        # Drop the decompressed copies of the (compressed) inputs - the final reads are
        # written, so these are now just a duplicate of data the user already has.
        for path in (decompressed if not dry_run() else []):
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
