"""MetaWrap2 quant_bins module: estimate the abundance of each bin in each sample.

Reads from every sample are quantified against the whole-assembly (or concatenated-bin)
transcriptome with salmon; the per-sample transcript counts are then collapsed to a
length-weighted median coverage per bin, giving a bins x samples abundance table. With more
than one sample a clustered Seaborn heatmap of those abundances is drawn as well.

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
from ..io.seqio import iter_fasta
from ..scripts import split_salmon_out_into_bins, summarize_salmon_files
from ._common import (
    announcement, collect_read_pairs, comm, ensure_dir, env_for, error, finish_run,
    make_checkpoint, run, start_run, warning,
)

CONDA_ENV = "metawrap2-quant_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
SALMON_INDEX = "salmon index -p {threads} -t {assembly} -i {index}"
SALMON_QUANT = "salmon quant -i {index} --libType IU -1 {r1} -2 {r2} -o {out} --meta -p {threads}"
MAKE_HEATMAP = ["python", "-m", "metawrap2.scripts.make_heatmap", "{table}", "{png}"]
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 quant_bins",
        usage="metawrap2 quant_bins [options] -b bins_folder -o output_dir -a assembly.fa "
              "readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]",
        description="Estimate the abundance of each bin across samples with salmon, and "
                    "draw a clustered abundance heatmap (.gz reads accepted).",
    )
    p.add_argument("-b", "--bins", required=True, help="folder containing draft genomes (bins) in fasta format")
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-a", "--assembly", default="",
                   help="fasta file with entire metagenomic assembly (strongly recommended!)")
    p.add_argument("-t", "--threads", type=int, default=DEFAULT_THREADS, help="number of threads (default 1)")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument("reads", nargs="+", help="paired read files (*_1.fastq/*_2.fastq, .gz accepted)")
    return p.parse_args(argv)


def _concat_bins_into_assembly(bin_folder: str, assembly: str) -> None:
    """Concatenate every bin's contigs into a single assembly fasta (supports .gz bins)."""
    if os.path.isfile(assembly):
        os.remove(assembly)
    with open(assembly, "w") as out:
        for f in sorted(os.listdir(bin_folder)):
            for header, seq in iter_fasta(os.path.join(bin_folder, f)):
                out.write(">%s\n%s\n" % (header, seq))


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not os.path.isdir(args.bins):
        error("%s does not exist. Exiting..." % args.bins)

    settings = load_settings(args.config)
    env = env_for("quant_bins", settings)
    rec = start_run("quant_bins", args, env, settings, inputs=[args.assembly, args.bins] + list(args.reads))

    ckpt = make_checkpoint(args.output)
    try:

        announcement("SETTING UP OUTPUT AND INDEXING ASSEMBLY")
        ensure_dir(args.output)

        assembly = args.assembly
        if not assembly:
            comm("Concatenating bins into a metagenomic assembly file.")
            assembly = os.path.join(args.output, "assembly.fa")
            _concat_bins_into_assembly(args.bins, assembly)
        elif not os.path.isfile(assembly):
            error("Assembly file %s does not exist. Exiting..." % assembly)

        index = os.path.join(args.output, "assembly_index")
        if ckpt.todo("index"):
            comm("Indexing assembly file with salmon. Ignore any warnings")
            run(SALMON_INDEX.format(threads=args.threads, assembly=assembly, index=index),
                env=env, tool="salmon index",
                hint="Something went wrong with indexing the assembly.")
            ckpt.done("index")
        else:
            comm("skipping salmon index (already done; --resume)")

        announcement("ALIGNING READS FROM ALL SAMPLES BACK TO BINS WITH SALMON")
        try:
            pairs = collect_read_pairs(args.reads)
        except ValueError as e:
            error(str(e))
        comm("%d forward and %d reverse read files detected" % (len(pairs), len(pairs)))

        align_dir = os.path.join(args.output, "alignment_files")
        ensure_dir(align_dir)
        if ckpt.todo("quant"):
            for sample, r1, r2 in pairs:
                out_quant = os.path.join(align_dir, sample + ".quant")
                comm("processing sample %s with reads %s and %s..." % (sample, r1, r2))
                run(SALMON_QUANT.format(index=index, r1=r1, r2=r2, out=out_quant, threads=args.threads),
                    env=env, tool="salmon quant",
                    hint="Something went wrong with aligning %s fastq files back to assembly!" % sample)
            ckpt.done("quant")
        else:
            comm("skipping salmon quant (already done; --resume)")

        quant_dir = os.path.join(args.output, "quant_files")
        if ckpt.todo("summarize"):
            comm("summarize salmon files...")
            summarize_salmon_files.summarize(align_dir)
            ensure_dir(quant_dir)
            for f in sorted(os.listdir(align_dir)):
                if f.endswith(".quant.counts"):
                    os.replace(os.path.join(align_dir, f), os.path.join(quant_dir, f))

            announcement("EXTRACTING AVERAGE ABUNDANCE OF EACH BIN")
            n = len([f for f in os.listdir(quant_dir) if "counts" in f])
            if n < 1:
                error("There were no files found in %s" % quant_dir)
            comm("There were %d samples detected. Making abundance table!" % n)
            table = os.path.join(args.output, "bin_abundance_table.tab")
            with open(table, "w") as out:
                split_salmon_out_into_bins.build_table(quant_dir, args.bins, assembly, out)
            comm("Average bin abundance table stored in %s" % table)

            if n > 1:
                announcement("MAKING GENOME ABUNDANCE HEATMAP WITH SEABORN")
                comm("making heatmap with Seaborn")
                heatmap = os.path.join(args.output, "bin_abundance_heatmap.png")
                run([a.format(table=table, png=heatmap) for a in MAKE_HEATMAP],
                    env=env, tool="make_heatmap", hint="Something went wrong with making the heatmap.")
                comm("cleaning up...")
                shutil.rmtree(align_dir, ignore_errors=True)
            else:
                warning("Cannot make clustered heatmap with just one sample... Skipping heatmap")
            ckpt.done("summarize")
        else:
            comm("skipping abundance summary (already done; --resume)")

        announcement("QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
