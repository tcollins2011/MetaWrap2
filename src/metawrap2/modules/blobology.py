"""MetaWrap2 blobology module: GC-vs-coverage blobplots of contigs (and bins).

Contigs are taxonomically classified with MEGABLAST against NCBI nt, reads are mapped back
with minimap2 to estimate coverage, and the two are combined into a blobplot table that is
rendered into GC-vs-abundance plots coloured by taxonomy (and, optionally, by bin).

This is a MetaWrap2 port of a pipeline derived from BLOBOLOGY (Sujai Kumar,
github.com/blaxterlab/blobology); the heavy lifting lives in metawrap2.scripts.blobology.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import glob
import os
from typing import List, Optional

from .. import blastdb as _blastdb
from ..aligners import ALIGNERS, DEFAULT_ALIGNER, INDEX_SUFFIXES, index_path
from ..cliformat import make_parser
from ..command import tool_path_in_env
from ..config import MODULE_ENVS, load_settings
from ..io.seqio import iter_fasta
from ..progress import bar
from ..pyrun import PY
from ._common import (
    absolutize_paths,
    announcement,
    check_fasta,
    check_fastq,
    collect_read_pairs,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["blobology"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
FASTAQUAL_SELECT = (
    PY + " -m metawrap2.scripts.blobology.fastaqual_select -f {assembly} -s r -n {n_contigs}"
)
BLASTN = (
    "blastn -task megablast -query {query} -db {blastdb} -evalue 1e-5"
    " -num_threads {threads} -max_target_seqs 1 -outfmt '6 qseqid staxids'"
)
#: The aligner streams into the BAM conversion, so no SAM is ever written. pipefail matters here:
#: the conversion exits 0 on a truncated stream and leaves a short but perfectly valid BAM, so
#: without it an aligner that died halfway would show up as a blobplot built on partial coverage.
#:
#: No sort, deliberately. gc_cov_annotate reads this with a plain ``samtools view`` and sums aligned
#: span per contig, so record order is irrelevant and sorting it would be work spent on nothing.
ALIGN_BAM = "set -o pipefail; {align} | samtools view -O BAM -b -@ {threads} - > {bam}"
# Runs in the host interpreter (it imports metawrap2) but shells out to samtools, which only
# exists inside this module's conda env - hence --samtools with an absolute path.
GC_COV_ANNOTATE = (
    PY + " -m metawrap2.scripts.blobology.gc_cov_annotate --blasttaxid {megablast}"
    " --assembly {assembly} --bam {bams} --out {out} --taxdump {taxdump}"
    " --samtools {samtools}"
    " --taxlist species genus family order subclass phylum superkingdom"
)
ADD_BINS = PY + " -m metawrap2.scripts.blobology.add_bins_to_blobplot {blobplot} {bin_folder}"
MAKEBLOBPLOT = PY + " -m metawrap2.scripts.blobology.makeblobplot {table} {prop} {taxlevel}"
MAKEBLOBPLOT_BASE = (
    PY + " -m metawrap2.scripts.blobology.makeblobplot {table} {prop} {taxlevel} {base}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
EVALUE = "1e-5"  # MEGABLAST e-value cutoff (baked into the .megablast filename)
PLOT_PROP = 0.005  # collapse taxonomy categories rarer than this fraction
BIN_PLOT_PROP = 0.000001  # keep even tiny per-bin categories when colouring by bin
TAXLEVELS = ("taxlevel_order", "taxlevel_phylum", "taxlevel_superkingdom")
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 blobology",
        usage="metawrap2 blobology [options] -a assembly.fasta -o output_dir "
        "readsA_1.fastq readsA_2.fastq [readsB_1.fastq readsB_2.fastq ...]",
        description="GC-vs-coverage blobplots of contigs (and bins). Provide each "
        "replicate's paired reads (*_1.fastq/*_2.fastq, .gz accepted).",
    )
    p.add_argument("-a", "--assembly", required=True, help="assembly fasta file")
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "--subsample",
        type=int,
        default=None,
        help="number of contigs to run blobology on (randomized; default: all)",
    )
    p.add_argument(
        "--bins",
        default=None,
        help="folder of bins; contig names must match the assembly (default: none)",
    )
    p.add_argument(
        "--aligner",
        choices=sorted(ALIGNERS),
        default=DEFAULT_ALIGNER,
        help="aligner used to measure contig abundance (default %s; see metawrap2.aligners)"
        % DEFAULT_ALIGNER,
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument("reads", nargs="+", help="paired read files (*_1.fastq / *_2.fastq, .gz ok)")
    return p.parse_args(argv)


def _bams(out: str) -> str:
    bams = sorted(glob.glob(os.path.join(out, "*.bam")))
    if not bams and dry_run():
        return os.path.join(out, "<sample>.bam")
    return " ".join(bams)


def _annotate_with_bins(blobplot: str, binned_blobplot: str, bin_folder: str) -> bool:
    """Add per-contig bin annotations to *blobplot* in place. Returns True if any matched.

    Ports the shell logic: annotate to a temp file; if fewer than two contigs were binned
    the bins probably don't match the assembly, so we drop the annotation.
    """
    tmp = blobplot + ".binned.tmp"
    run(
        ADD_BINS.format(blobplot=blobplot, bin_folder=bin_folder),
        env=None,
        tool="add_bins_to_blobplot",
        log_path=tmp,
    )
    with open(tmp) as fh:
        binned = [ln for ln in fh if "Unbinned" not in ln]
    if len(binned) < 2:
        warning(
            "No contig matches were found in the bins provided. This may be because the "
            "contigs in the %s folder are not the same as the contigs in the assembly. "
            "The blobplot will not be annotated with bins." % bin_folder
        )
        os.remove(tmp)
        return False
    os.replace(tmp, blobplot)
    # binned-only blobplot: every line that is not "Unbinned"
    with open(blobplot) as fh, open(binned_blobplot, "w") as out:
        for ln in fh:
            if "Unbinned" not in ln:
                out.write(ln)
    return True


def _plot(
    table: str, prop: float, taxlevel: str, env: Optional[str], base: Optional[str] = None
) -> None:
    if base is None:
        cmd = MAKEBLOBPLOT.format(table=table, prop=prop, taxlevel=taxlevel)
    else:
        cmd = MAKEBLOBPLOT_BASE.format(table=table, prop=prop, taxlevel=taxlevel, base=base)
    run(cmd, env=None, tool="makeblobplot")


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "blobology",
        [(check_fasta, args.assembly, "assembly (-a)")]
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    env = env_for("blobology", settings)
    rec = start_run(
        "blobology", args, env, settings, inputs=[args.assembly, args.bins] + list(args.reads)
    )

    ckpt = make_checkpoint(args.output)
    try:

        if not os.path.isfile(args.assembly):
            error("Assembly file %s does not exist. Exiting..." % args.assembly)

        blastdb = settings.db("BLASTDB")
        taxdump = settings.db("TAXDUMP")
        ok, problem = _blastdb.find(settings)
        if not ok:
            error(problem or "the configured BLAST database cannot be used")
        blastdb = _blastdb.db_path(settings)
        if not os.path.isfile(os.path.join(taxdump, "citations.dmp")):
            error(
                "The file %s/citations.dmp does not exist, which likely means the NCBI taxonomy "
                "database path (TAXDUMP) is not set correctly or is not downloaded." % taxdump
            )

        try:
            pairs = collect_read_pairs(args.reads)
        except ValueError as e:
            error(str(e))
        comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

        out = args.output
        ensure_dir(out)
        assembly = os.path.basename(args.assembly)
        sample = os.path.splitext(assembly)[0]
        assembly_copy = os.path.join(out, assembly)

        megablast = os.path.join(out, "%s.nt.%s.megablast" % (sample, EVALUE))
        if ckpt.todo("megablast"):
            announcement("ASSIGN TAXONOMY TO CONTIGS WITH MEGABLAST")
            if args.subsample is None:
                comm("making copy of assembly file %s" % args.assembly)
                if not dry_run():
                    with open(assembly_copy, "w") as fh:
                        fh.writelines(
                            ">%s\n%s\n" % (header, seq) for header, seq in iter_fasta(args.assembly)
                        )
            else:
                comm("Choosing %d contigs at random from %s" % (args.subsample, args.assembly))
                run(
                    FASTAQUAL_SELECT.format(assembly=args.assembly, n_contigs=args.subsample),
                    env=None,
                    tool="fastaqual_select",
                    log_path=assembly_copy,
                )

            comm("Running MEGABLAST on %s" % assembly)
            if dry_run() or not (os.path.isfile(megablast) and os.path.getsize(megablast)):
                run(
                    BLASTN.format(query=assembly_copy, blastdb=blastdb, threads=args.threads),
                    env=env,
                    tool="blastn",
                    log_path=megablast,
                )
            else:
                comm("Looks like taxonomy assignment was already run. Skipping...")
            if not dry_run() and not (os.path.isfile(megablast) and os.path.getsize(megablast)):
                error(
                    "Something went wrong with assigning taxonomy to contigs with megablast. Exiting"
                )
            ckpt.done("megablast")
        else:
            comm("skipping taxonomy assignment with megablast (already done; --resume)")

        if ckpt.todo("map_reads"):
            announcement("MAP READS TO ASSEMBLY WITH %s" % args.aligner.upper())
            aligner = ALIGNERS[args.aligner]
            index = index_path(assembly_copy, args.aligner)
            comm("Indexing %s" % assembly_copy)
            run(
                aligner["index_cmd"].format(assembly=assembly_copy, index=index),
                env=env,
                tool="%s index" % args.aligner,
            )

            for sname, r1, r2 in bar(
                pairs, desc="mapping samples", unit=" sample", total=len(pairs)
            ):
                bam = os.path.join(out, sname + ".bam")
                if not dry_run() and os.path.isfile(bam) and os.path.getsize(bam):
                    comm("Looks like the alignment file for %s already exists. Skipping..." % sname)
                    continue
                comm(
                    "Now processing sample %s ... Aligning %s and %s to %s with %s"
                    % (sname, r1, r2, assembly_copy, args.aligner)
                )
                align = aligner["paired"].format(threads=args.threads, index=index, r1=r1, r2=r2)
                run(
                    [
                        "bash",
                        "-c",
                        ALIGN_BAM.format(align=align, threads=args.threads, bam=bam),
                    ],
                    env=env,
                    tool="%s | samtools view" % args.aligner,
                )
            ckpt.done("map_reads")
        else:
            comm("skipping read mapping (already done; --resume)")

        blobplot = os.path.join(out, sample + ".blobplot")
        binned_blobplot = os.path.join(out, sample + ".binned.blobplot")
        have_bins = False
        if ckpt.todo("blobfile"):
            announcement("MAKE BLOB FILE FROM BLAST AND ALIGNMENT OUTPUT")
            if dry_run() or not (os.path.isfile(blobplot) and os.path.getsize(blobplot)):
                run(
                    GC_COV_ANNOTATE.format(
                        megablast=megablast,
                        assembly=assembly_copy,
                        bams=_bams(out),
                        out=blobplot,
                        taxdump=taxdump,
                        samtools=tool_path_in_env(env, "samtools"),
                    ),
                    env=None,
                    tool="gc_cov_annotate",
                )
                if not dry_run() and not (os.path.isfile(blobplot) and os.path.getsize(blobplot)):
                    error(
                        "Something went wrong with making the blob file from the .bam and .megablast "
                        "files. Exiting."
                    )
                comm("blobplot text file saved to %s" % blobplot)
            else:
                comm("blobplot text file already exists")

            if args.bins:
                comm("adding bin annotations to blobfile %s" % blobplot)
                have_bins = (
                    True if dry_run() else _annotate_with_bins(blobplot, binned_blobplot, args.bins)
                )
            ckpt.done("blobfile")
        else:
            comm("skipping blob file creation (already done; --resume)")
            have_bins = (
                bool(args.bins)
                and os.path.isfile(binned_blobplot)
                and os.path.getsize(binned_blobplot) > 0
            )

        if ckpt.todo("plots"):
            announcement("MAKE FINAL BLOBPLOT IMAGES")
            comm("making blobplots with phylogeny annotations (order, phylum, and superkingdom).")
            for level in TAXLEVELS:
                _plot(blobplot, PLOT_PROP, level, env)
            if not dry_run() and not os.path.isfile(blobplot + ".taxlevel_phylum.png"):
                error("Something went wrong with making the plots from the blob file. Exiting.")

            if have_bins:
                comm("making blobplot images with bin annotations")
                _plot(blobplot, BIN_PLOT_PROP, "bin", env, base="Unbinned")
                _plot(blobplot, PLOT_PROP, "binned_yes_no", env, base="Unbinned")
                _plot(blobplot, PLOT_PROP, "binned_phylum", env, base="Unbinned")

                comm("making blobplot images of only the contigs that were binned")
                _plot(binned_blobplot, BIN_PLOT_PROP, "bin", env, base="Unbinned")
                comm(
                    "making blobplots of only binned contigs with phylogeny annotations "
                    "(order, phylum, and superkingdom)."
                )
                for level in TAXLEVELS:
                    _plot(binned_blobplot, PLOT_PROP, level, env)

                binned_dir = os.path.join(out, "blobplot_figures_only_binned_contigs")
                ensure_dir(binned_dir)
                for png in glob.glob(os.path.join(out, "*.binned.blobplot*png")):
                    os.replace(png, os.path.join(binned_dir, os.path.basename(png)))

            figures = os.path.join(out, "blobplot_figures")
            ensure_dir(figures)
            for png in glob.glob(os.path.join(out, "*png")):
                os.replace(png, os.path.join(figures, os.path.basename(png)))
            ckpt.done("plots")
        else:
            comm("skipping blobplot image generation (already done; --resume)")

        comm("cleaning up...")
        if not dry_run():
            for f in glob.glob(os.path.join(out, "*.bam")) + [
                assembly_copy + suffix for suffix in INDEX_SUFFIXES
            ]:
                if os.path.isfile(f):
                    os.remove(f)

        announcement("BLOBPLOT PIPELINE FINISHED SUCCESSFULLY!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
