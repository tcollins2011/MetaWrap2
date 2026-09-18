"""MetaWrap2 annotate_bins module: quick functional annotation of bins with PROKKA.

Each bin's contig names are shortened (PROKKA rejects long ids), PROKKA annotates the bin,
and the results are reorganized into three per-run folders: the product-carrying gff lines
(``bin_funct_annotations``), translated genes (``bin_translated_genes``, .faa) and
untranslated genes (``bin_untranslated_genes``, .ffn).

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import os
import shutil
from typing import List

from ..config import load_settings
from ..scripts import shorten_contig_names
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run,
    start_run, warning,
)

CONDA_ENV = "metawrap2-annotate_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
PROKKA = "prokka --quiet --cpus {threads} --outdir {outdir} --prefix {prefix} {input}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 annotate_bins",
        usage="metawrap2 annotate_bins [options] -o output_dir -b bin_folder",
        description="Functionally annotate a set of bins with PROKKA.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-b", "--bins", required=True, help="folder with metagenomic bins in fasta format")
    p.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="number of threads (default 1)")
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _write_shortened(bin_file: str, tmp_bin: str) -> None:
    with open(tmp_bin, "w") as out:
        for line in shorten_contig_names.shorten_lines(bin_file):
            out.write(line + "\n")


def _grep_product(gff_in: str, gff_out: str) -> None:
    with open(gff_in) as fh, open(gff_out, "w") as out:
        for line in fh:
            if "product" in line:
                out.write(line)


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not os.path.isdir(args.bins):
        error("%s does not exist! Exiting." % args.bins)

    settings = load_settings(args.config)
    env = env_for("annotate_bins", settings)
    rec = start_run("annotate_bins", args, env, settings, inputs=[args.bins])

    ckpt = make_checkpoint(args.output)
    try:

        announcement("BEGIN ANNNOTAION PIPELINE!")
        comm("setting up output folder and copying over bins...")
        ensure_dir(args.output)
        prokka_out = os.path.join(args.output, "prokka_out")

        if ckpt.todo("prokka"):
            if os.path.isdir(prokka_out):
                shutil.rmtree(prokka_out)
            ensure_dir(prokka_out)

            tmp_bin = os.path.join(args.output, "tmp_bin.fa")
            for i in sorted(os.listdir(args.bins)):
                bin_name = os.path.splitext(i)[0]
                _write_shortened(os.path.join(args.bins, i), tmp_bin)
                comm("NOW ANNOTATING %s" % bin_name)
                rc = run(PROKKA.format(threads=args.threads, outdir=os.path.join(prokka_out, bin_name),
                                       prefix=bin_name, input=tmp_bin),
                         env=env, tool="prokka", check=False)
                if rc != 0:
                    warning("Something possibly went wrong with annotating %s. Proceeding anyways" % bin_name)
                gff = os.path.join(prokka_out, bin_name, bin_name + ".gff")
                if not (os.path.isfile(gff) and os.path.getsize(gff)):
                    error("Something went wrong with annotating %s. Exiting..." % bin_name)
                os.remove(tmp_bin)

            if not os.listdir(prokka_out):
                error("Something went wrong with running prokka on all the bins! Exiting...")
            comm("PROKKA finished annotating all the bins!")
            ckpt.done("prokka")
        else:
            comm("skipping PROKKA annotation (already done; --resume)")

        if ckpt.todo("format"):
            announcement("FORMATTING ANNNOTAIONS...")
            funct = os.path.join(args.output, "bin_funct_annotations")
            translated = os.path.join(args.output, "bin_translated_genes")
            untranslated = os.path.join(args.output, "bin_untranslated_genes")
            ensure_dir(funct)
            ensure_dir(translated)
            ensure_dir(untranslated)
            for i in sorted(os.listdir(prokka_out)):
                _grep_product(os.path.join(prokka_out, i, i + ".gff"),
                              os.path.join(funct, i + ".gff"))
                shutil.copy(os.path.join(prokka_out, i, i + ".faa"), translated)
                shutil.copy(os.path.join(prokka_out, i, i + ".ffn"), untranslated)

            comm("You will find the bin annotation gff files in %s." % funct)
            warning("Your contigs may be truncated in the annotation because they are too long for "
                    "PROKKA to take as input! Check the final annotation files to see what the contig "
                    "naming convention is.")
            ckpt.done("format")
        else:
            comm("skipping annotation formatting (already done; --resume)")

        announcement("ANNOTATE BINS PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
