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

from ..config import OPTIONAL_ENVS, load_settings
from ..constants import FASTA_EXTENSIONS
from ..progress import bar
from ..scripts import shorten_contig_names
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = "metawrap2-annotate_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
PROKKA = "prokka --quiet --cpus {threads} --outdir {outdir} --prefix {prefix} {input}"
# --bakta path. Bakta writes <prefix>.gff3/.faa/.ffn, so the gff extension differs from
# prokka's .gff; _annotator_outputs() accounts for that.
BAKTA = "bakta --threads {threads} --output {outdir} --prefix {prefix} --force {db} {input}"
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
    p.add_argument(
        "-b", "--bins", required=True, help="folder with metagenomic bins in fasta format"
    )
    threads_arg(p)
    p.add_argument(
        "--bakta",
        action="store_true",
        help="annotate with Bakta instead of PROKKA (needs the "
        "metawrap2-annotate_bins-bakta env and a Bakta database)",
    )
    p.add_argument(
        "--bakta-db",
        help="path to the Bakta database directory (default: "
        "BAKTA_DB from metawrap2.toml, else Bakta's own default)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _write_shortened(bin_file: str, tmp_bin: str) -> None:
    with open(tmp_bin, "w") as out:
        out.writelines(line + "\n" for line in shorten_contig_names.shorten_lines(bin_file))


def _annotator(args, settings):
    """Return (name, conda env, gff extension) for the selected annotator."""
    if args.bakta:
        env = OPTIONAL_ENVS["bakta"] if settings.use_conda_envs else None
        return "bakta", env, ".gff3"
    return "prokka", env_for("annotate_bins", settings), ".gff"


def _annotate_one(
    args, tool: str, env, outdir: str, prefix: str, tmp_bin: str, bakta_db: str
) -> int:
    if tool == "bakta":
        db = "--db %s" % bakta_db if bakta_db else ""
        cmd = BAKTA.format(threads=args.threads, outdir=outdir, prefix=prefix, db=db, input=tmp_bin)
    else:
        cmd = PROKKA.format(threads=args.threads, outdir=outdir, prefix=prefix, input=tmp_bin)
    return run(cmd, env=env, tool=tool, check=False)


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
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "annotate_bins",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                args.bins,
                "bin folder (-b)",
            )
        ],
    )
    tool, env, gff_ext = _annotator(args, settings)
    bakta_db = args.bakta_db or settings.db("BAKTA_DB")
    if args.bakta and not bakta_db:
        comm("No Bakta database given (--bakta-db / BAKTA_DB); relying on Bakta's own default.")
    rec = start_run("annotate_bins", args, env, settings, inputs=[args.bins])

    ckpt = make_checkpoint(args.output)
    try:

        announcement("BEGIN ANNNOTAION PIPELINE!")
        comm("setting up output folder and copying over bins...")
        ensure_dir(args.output)
        # Output folder keeps the historical "prokka_out" name whichever annotator ran, so
        # downstream paths and users' scripts do not change when --bakta is used.
        prokka_out = os.path.join(args.output, "prokka_out")

        if ckpt.todo("prokka"):
            if os.path.isdir(prokka_out):
                shutil.rmtree(prokka_out)
            ensure_dir(prokka_out)

            tmp_bin = os.path.join(args.output, "tmp_bin.fa")
            bin_files = [
                f
                for f in sorted(os.listdir(args.bins))
                if f.endswith(FASTA_EXTENSIONS) and os.path.isfile(os.path.join(args.bins, f))
            ]
            if not bin_files:
                error(
                    "No fasta bins (%s) found in %s. Exiting."
                    % ("/".join(FASTA_EXTENSIONS), args.bins)
                )
            for i in bar(bin_files, desc="annotating bins", unit=" bin", total=len(bin_files)):
                bin_name = os.path.splitext(i)[0]
                if not dry_run():
                    _write_shortened(os.path.join(args.bins, i), tmp_bin)
                comm("NOW ANNOTATING %s" % bin_name)
                rc = _annotate_one(
                    args, tool, env, os.path.join(prokka_out, bin_name), bin_name, tmp_bin, bakta_db
                )
                if rc != 0:
                    warning(
                        "Something possibly went wrong with annotating %s. Proceeding anyways"
                        % bin_name
                    )
                gff = os.path.join(prokka_out, bin_name, bin_name + gff_ext)
                if not dry_run():
                    if not (os.path.isfile(gff) and os.path.getsize(gff)):
                        error("Something went wrong with annotating %s. Exiting..." % bin_name)
                    os.remove(tmp_bin)

            if not dry_run() and not os.listdir(prokka_out):
                error("Something went wrong with running %s on all the bins! Exiting..." % tool)
            comm("%s finished annotating all the bins!" % tool.upper())
            ckpt.done("prokka")
        else:
            comm("skipping %s annotation (already done; --resume)" % tool.upper())

        if ckpt.todo("format") and not dry_run():
            announcement("FORMATTING ANNNOTAIONS...")
            funct = os.path.join(args.output, "bin_funct_annotations")
            translated = os.path.join(args.output, "bin_translated_genes")
            untranslated = os.path.join(args.output, "bin_untranslated_genes")
            ensure_dir(funct)
            ensure_dir(translated)
            ensure_dir(untranslated)
            for i in sorted(os.listdir(prokka_out)):
                _grep_product(
                    os.path.join(prokka_out, i, i + gff_ext), os.path.join(funct, i + ".gff")
                )
                shutil.copy(os.path.join(prokka_out, i, i + ".faa"), translated)
                shutil.copy(os.path.join(prokka_out, i, i + ".ffn"), untranslated)

            comm("You will find the bin annotation gff files in %s." % funct)
            warning(
                "Your contigs may be truncated in the annotation because they are too long for "
                "%s to take as input! Check the final annotation files to see what the contig "
                "naming convention is." % tool.upper()
            )
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
