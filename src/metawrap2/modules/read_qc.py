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
import shutil
from typing import List

from ..config import load_settings
from ..io.seqio import smart_open
from ..scripts import select_human_reads, skip_human_reads
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run,
    start_run, warning,
)

CONDA_ENV = "metawrap2-read_qc"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
FASTQC       = "fastqc -q -t {threads} -o {outdir} -f fastq {r1} {r2}"
TRIM_GALORE  = "trim_galore --no_report_file --paired -o {outdir} {r1} {r2}"
BMTAGGER     = ("bmtagger.sh -b {bitmask} -x {srprism} -T {tmpdir} -q1"
                " -1 {r1} -2 {r2} -o {listfile}")
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_HOST = "hg38"   # prefix of the host index inside the bmtagger database folder
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 read_qc",
        usage="metawrap2 read_qc [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir",
        description="Quality-trim reads and remove host sequences. Read files must follow the "
                    "name_1.fastq / name_2.fastq convention (.gz accepted).",
    )
    p.add_argument("-1", dest="reads_1", required=True, help="forward fastq reads")
    p.add_argument("-2", dest="reads_2", required=True, help="reverse fastq reads")
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-t", "--threads", type=int, default=1, help="number of threads (default 1)")
    p.add_argument("-x", "--host", default=DEFAULT_HOST,
                   help="prefix of host index in bmtagger database folder (default hg38)")
    p.add_argument("--skip-bmtagger", action="store_true",
                   help="dont remove host sequences with bmtagger")
    p.add_argument("--skip-trimming", action="store_true",
                   help="dont trim sequences with trimgalore")
    p.add_argument("--skip-pre-qc-report", action="store_true",
                   help="dont make FastQC report of input sequences")
    p.add_argument("--skip-post-qc-report", action="store_true",
                   help="dont make FastQC report of final sequences")
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _sample_name(reads_1: str) -> str:
    """Sample name from a *_1.fastq(.gz) file, matching the original ${tmp%_*} behavior."""
    base = os.path.basename(reads_1)
    for token in ("_1.fastq", "_1.fq"):
        if token in base:
            return base.split(token)[0]
    return os.path.splitext(base)[0]


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


def _fastqc_report(args, env, name: str, r1: str, r2: str) -> None:
    outdir = os.path.join(args.output, name)
    ensure_dir(outdir)
    run(FASTQC.format(threads=args.threads, outdir=outdir, r1=r1, r2=r2),
        env=env, tool="fastqc")
    for zipf in glob.glob(os.path.join(outdir, "*.zip")):
        os.remove(zipf)
    comm("%s saved to: %s" % (name, outdir))


def main(argv: List[str]) -> int:
    args = _parse_args(argv)

    if args.reads_1 == args.reads_2:
        error("The forward and reverse reads are the same file. Exiting pipeline.")
    if not os.path.isfile(args.reads_1):
        error("%s file does not exist. Exiting..." % args.reads_1)
    if not os.path.isfile(args.reads_2):
        error("%s file does not exist. Exiting..." % args.reads_2)

    settings = load_settings(args.config)
    env = env_for("read_qc", settings)
    rec = start_run("read_qc", args, env, settings, inputs=[args.reads_1, args.reads_2])

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

        reads_1 = _decompress(args.reads_1, args.output)
        reads_2 = _decompress(args.reads_2, args.output)
        sample = _sample_name(reads_1)

        if not args.skip_pre_qc_report and ckpt.todo("pre_qc_report"):
            announcement("MAKING PRE-QC REPORT")
            _fastqc_report(args, env, "pre-QC_report", reads_1, reads_2)
            ckpt.done("pre_qc_report")

        if not args.skip_trimming and ckpt.todo("trimming"):
            announcement("RUNNING TRIM-GALORE")
            run(TRIM_GALORE.format(outdir=args.output, r1=reads_1, r2=reads_2),
                env=env, tool="trim_galore")
            trimmed_1 = os.path.join(args.output, "trimmed_1.fastq")
            trimmed_2 = os.path.join(args.output, "trimmed_2.fastq")
            os.replace(os.path.join(args.output, "%s_1_val_1.fq" % sample), trimmed_1)
            os.replace(os.path.join(args.output, "%s_2_val_2.fq" % sample), trimmed_2)
            if not os.path.getsize(trimmed_1):
                error("Something went wrong with trimming the reads. Exiting.")
            comm("Trimmed reads saved to: %s and %s" % (trimmed_1, trimmed_2))
            reads_1, reads_2 = trimmed_1, trimmed_2
            for leftover in ("%s_1_trimmed.fq" % sample, "%s_2_trimmed.fq" % sample):
                path = os.path.join(args.output, leftover)
                if os.path.isfile(path):
                    os.remove(path)
            ckpt.done("trimming")
        else:
            # on --resume, pick up the trimmed reads a previous run left behind
            trimmed_1 = os.path.join(args.output, "trimmed_1.fastq")
            trimmed_2 = os.path.join(args.output, "trimmed_2.fastq")
            if os.path.isfile(trimmed_1) and os.path.isfile(trimmed_2):
                comm("skipping read trimming (already done; --resume)")
                reads_1, reads_2 = trimmed_1, trimmed_2

        if do_bmtagger and ckpt.todo("host_removal"):
            announcement("REMOVING HOST SEQUENCES WITH BMTAGGER")
            tmpdir = os.path.join(args.output, "bmtagger_tmp")
            ensure_dir(tmpdir)
            bitmask = os.path.join(bmtagger_db, args.host + ".bitmask")
            srprism = os.path.join(bmtagger_db, args.host + ".srprism")
            listfile = os.path.join(args.output, "%s.bmtagger.list" % sample)
            comm("running bmtagger with %s %s indexes..." % (bitmask, srprism))
            run(BMTAGGER.format(bitmask=bitmask, srprism=srprism, tmpdir=tmpdir,
                                r1=reads_1, r2=reads_2, listfile=listfile),
                env=env, tool="bmtagger.sh")
            if not (os.path.isfile(listfile) and os.path.getsize(listfile)):
                warning("No contamination reads found, which is very unlikely.")

            clean_1 = os.path.join(args.output, "%s_1.clean.fastq" % sample)
            clean_2 = os.path.join(args.output, "%s_2.clean.fastq" % sample)
            comm("Now sorting out found human reads from the main fastq files...")
            with open(clean_1, "w") as out:
                skip_human_reads.filter_reads(listfile, reads_1, out)
            with open(clean_2, "w") as out:
                skip_human_reads.filter_reads(listfile, reads_2, out)

            comm("Now sorting out found human reads and putting them into a new file... for science...")
            with open(os.path.join(args.output, "host_reads_1.fastq"), "w") as out:
                select_human_reads.select_reads(listfile, reads_1, out)
            with open(os.path.join(args.output, "host_reads_2.fastq"), "w") as out:
                select_human_reads.select_reads(listfile, reads_2, out)
            if not os.path.getsize(clean_1):
                error("Something went wrong with removing contaminant reads with bmtagger. Exiting.")

            shutil.rmtree(tmpdir, ignore_errors=True)
            os.remove(listfile)
            reads_1, reads_2 = clean_1, clean_2
            for trimmed in ("trimmed_1.fastq", "trimmed_2.fastq"):
                path = os.path.join(args.output, trimmed)
                if os.path.isfile(path):
                    os.remove(path)
            ckpt.done("host_removal")
        elif do_bmtagger:
            # on --resume, pick up the host-filtered reads a previous run left behind
            clean_1 = os.path.join(args.output, "%s_1.clean.fastq" % sample)
            clean_2 = os.path.join(args.output, "%s_2.clean.fastq" % sample)
            if os.path.isfile(clean_1) and os.path.isfile(clean_2):
                comm("skipping host read removal (already done; --resume)")
                reads_1, reads_2 = clean_1, clean_2

        final_1 = os.path.join(args.output, "final_pure_reads_1.fastq")
        final_2 = os.path.join(args.output, "final_pure_reads_2.fastq")
        if ckpt.todo("finalize"):
            os.replace(reads_1, final_1)
            os.replace(reads_2, final_2)
            comm("Contamination-free and trimmed reads are stored in: %s and %s" % (final_1, final_2))
            ckpt.done("finalize")

        if not args.skip_post_qc_report and ckpt.todo("post_qc_report"):
            announcement("MAKING POST-QC REPORT")
            _fastqc_report(args, env, "post-QC_report", final_1, final_2)
            ckpt.done("post_qc_report")

        announcement("READ QC PIPELINE COMPLETE!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
