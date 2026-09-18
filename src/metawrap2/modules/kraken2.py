"""MetaWrap2 kraken2 module: taxonomic profiling of reads and/or contigs with KRAKEN2.

Runs KRAKEN2 on any number of assembly FASTA files and/or paired-end read sets, translates
the taxid output into full taxonomy lineages, collapses each into KRONA format, and draws a
single interactive kronagram of all samples with KronaTools.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers/defaults in the CONSTANTS block. Change a flag there
and it takes effect — no need to follow the orchestration logic underneath.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import glob
import os
import random
from typing import List

from ..config import load_settings
from ..io.seqio import smart_open
from ..scripts import kraken2_translate, kraken_to_krona
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run,
    start_run, warning,
)

CONDA_ENV = "metawrap2-kraken2"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
KRAKEN2_PAIRED    = "kraken2 --use-names --db {db} --paired --threads {threads} --output {out} {r1} {r2}"
KRAKEN2_SINGLE    = "kraken2 --use-names --db {db} --threads {threads} --output {out} {seq}"
MEMORY_MAPPING    = "--memory-mapping"   # appended when --no-preload is given (lower memory, slower)
KTIMPORTTEXT      = "ktImportText -o {out} {krona_files}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
FASTA_EXTS = (".fa", ".fasta", ".fa.gz", ".fasta.gz")
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 kraken2",
        usage="metawrap2 kraken2 [options] -o output_dir assembly.fasta reads_1.fastq reads_2.fastq ...",
        description="Taxonomic profiling of any number of assembly FASTA files (*.fa/*.fasta) "
                    "and/or paired-end reads (*_1.fastq/*_2.fastq, .gz accepted).",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-t", "--threads", type=int, default=1, help="number of threads (default 1)")
    p.add_argument("-s", "--subsample", type=int, default=None,
                   help="read subsampling number (default: all reads)")
    p.add_argument("--no-preload", action="store_true",
                   help="do not pre-load the kraken2 DB into memory (slower, lower memory)")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument("seqs", nargs="+", help="assembly FASTA files and/or paired-end read files")
    return p.parse_args(argv)


def _is_fasta(path: str) -> bool:
    return path.endswith(FASTA_EXTS)


def _fasta_sample(path: str) -> str:
    base = os.path.basename(path)
    if base.endswith(".gz"):
        base = base[:-3]
    return os.path.splitext(base)[0]


def _fastq_sample(reads_1: str) -> str:
    base = os.path.basename(reads_1)
    for token in ("_1.fastq", "_1.fq"):
        if token in base:
            return base.split(token)[0]
    return os.path.splitext(base)[0]


def _subsample_pairs(r1: str, r2: str, depth: int, out1: str, out2: str) -> None:
    """Randomly subsample *depth* read pairs from r1/r2 into out1/out2 (gz-transparent input)."""
    def _records(path):
        recs = []
        with smart_open(path, "rt") as fh:
            block = []
            for line in fh:
                block.append(line.rstrip("\n"))
                if len(block) == 4:
                    recs.append(block)
                    block = []
        return recs

    recs1 = _records(r1)
    recs2 = _records(r2)
    pairs = list(zip(recs1, recs2))
    random.shuffle(pairs)
    pairs = pairs[:depth]
    with open(out1, "w") as o1, open(out2, "w") as o2:
        for a, b in pairs:
            o1.write("\n".join(a) + "\n")
            o2.write("\n".join(b) + "\n")


def _kraken2(cmd: str, args) -> str:
    if args.no_preload:
        cmd = cmd + " " + MEMORY_MAPPING
    return cmd


def _run_kraken_on_inputs(args, env, db) -> None:
    for num in args.seqs:
        if num.endswith(("_1.fastq", "_1.fq", "_1.fastq.gz", "_1.fq.gz")):
            reads_1 = num
            reads_2 = num.replace("_1.fastq", "_2.fastq").replace("_1.fq", "_2.fq")
            if not os.path.isfile(reads_2):
                error("%s does not exist. Exiting..." % reads_2)
            sample = _fastq_sample(reads_1)
            comm("Now processing %s and %s with %d threads" % (reads_1, reads_2, args.threads))

            if args.subsample is not None:
                comm("subsampling down to %d reads..." % args.subsample)
                tmp_1 = os.path.join(args.output, "tmp_1.fastq")
                tmp_2 = os.path.join(args.output, "tmp_2.fastq")
                _subsample_pairs(reads_1, reads_2, args.subsample, tmp_1, tmp_2)
                reads_1, reads_2 = tmp_1, tmp_2
                comm("Subsampling done. Starting KRAKEN...")
                if not os.path.getsize(reads_1):
                    error("something went wrong with subsampling sequences. Exiting...")

            krak2 = os.path.join(args.output, "%s.krak2" % sample)
            if os.path.isfile(krak2) and os.path.getsize(krak2):
                comm("%s already exists - skipping running kraken2" % krak2)
            else:
                cmd = _kraken2(KRAKEN2_PAIRED.format(db=db, threads=args.threads, out=krak2,
                                                     r1=reads_1, r2=reads_2), args)
                run(cmd, env=env, tool="kraken2")
            if not (os.path.isfile(krak2) and os.path.getsize(krak2)):
                error("Something went wrong with running kraken2 on %s and %s . Exiting..."
                      % (reads_1, reads_2))
            if args.subsample is not None:
                os.remove(os.path.join(args.output, "tmp_1.fastq"))
                os.remove(os.path.join(args.output, "tmp_2.fastq"))

        elif _is_fasta(num):
            sample = _fasta_sample(num)
            comm("Now processing %s with %d threads" % (num, args.threads))
            krak2 = os.path.join(args.output, "%s.krak2" % sample)
            if os.path.isfile(krak2) and os.path.getsize(krak2):
                comm("%s already exists - skipping running kraken2" % krak2)
            else:
                cmd = _kraken2(KRAKEN2_SINGLE.format(db=db, threads=args.threads, out=krak2, seq=num),
                               args)
                run(cmd, env=env, tool="kraken2")
            if not (os.path.isfile(krak2) and os.path.getsize(krak2)):
                error("Something went wrong with running kraken2 on %s. Exiting..." % num)


def main(argv: List[str]) -> int:
    args = _parse_args(argv)

    settings = load_settings(args.config)
    env = env_for("kraken2", settings)
    rec = start_run("kraken2", args, env, settings, inputs=list(args.seqs))

    ckpt = make_checkpoint(args.output)
    try:
        db = settings.db("KRAKEN2_DB")
        if not (db and os.path.isdir(db)):
            error("The folder %s doesnt exist. Please consult the MetaWrap2 database guide to "
                  "download and build the KRAKEN2 database" % db)

        announcement("RUNNING KRAKEN ON ALL FILES")
        if os.path.isdir(args.output):
            warning("%s already exists." % args.output)
        ensure_dir(args.output)

        if ckpt.todo("kraken2"):
            _run_kraken_on_inputs(args, env, db)
            ckpt.done("kraken2")
        else:
            comm("skipping kraken2 on inputs (already done; --resume)")

        krak2_files = sorted(glob.glob(os.path.join(args.output, "*.krak2")))
        if not krak2_files:
            error("No fasta or fastq files detected! (must be in .fastq .fa .fasta or .fq format)")

        if ckpt.todo("translate"):
            announcement("RUNNING KRAKEN-TRANSLATE ON OUTPUT")
            for krak2 in krak2_files:
                comm("Translating %s" % krak2)
                translated = os.path.splitext(krak2)[0] + ".kraken2"
                kraken2_translate.translate_kraken2_annotations(annotation_file=krak2, kraken2_db=db,
                                                                output=translated)
            ckpt.done("translate")
        else:
            comm("skipping kraken-translate (already done; --resume)")

        if ckpt.todo("krona"):
            announcement("MAKING KRONAGRAM OF ALL FILES")
            for translated in sorted(glob.glob(os.path.join(args.output, "*.kraken2"))):
                krona = os.path.splitext(translated)[0] + ".krona"
                with open(krona, "w") as out:
                    kraken_to_krona.to_krona(translated, out)
                if not os.path.getsize(krona):
                    error("Something went wrong with making krona file from kraken file. Exiting...")

            krona_files = sorted(glob.glob(os.path.join(args.output, "*.krona")))
            kronagram = os.path.join(args.output, "kronagram.html")
            run(KTIMPORTTEXT.format(out=kronagram, krona_files=" ".join(krona_files)),
                env=env, tool="ktImportText")
            if not (os.path.isfile(kronagram) and os.path.getsize(kronagram)):
                error("Something went wrong with running KronaTools to make kronagram. Exiting...")
            ckpt.done("krona")
        else:
            comm("skipping kronagram (already done; --resume)")

        announcement("FINISHED RUNNING KRAKEN2 PIPELINE!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
