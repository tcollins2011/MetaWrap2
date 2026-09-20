"""MetaWrap2 assembly module: produce a metagenomic assembly from paired-end reads.

Assembles QC'd paired-end reads with MEGAHIT and/or metaSPAdes. When both are selected the
reads are first assembled with metaSPAdes, the reads that do not map back to the long
scaffolds are pulled out, and those leftovers are assembled with MEGAHIT (which does better
on low-coverage data); the two assemblies are then combined, sorted, and short contigs are
removed. The final assembly is QC'd with QUAST.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers/defaults in the CONSTANTS block. Change a flag there
and it takes effect — no need to follow the orchestration logic underneath.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import os
import shutil
from typing import List

from ..config import load_settings
from ..scripts import fix_megahit_contig_naming, rm_short_contigs, sam_to_fastq, sort_contigs
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
    expect_file,
    expect_produced,
    finish_run,
    fs_move,
    make_checkpoint,
    resolve_threads,
    run,
    scratch_dir,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = "metawrap2-assembly"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
METASPADES = "metaspades.py --tmp-dir {tmp} -t {threads} -m {mem} -o {out} -1 {r1} -2 {r2}"
METASPADES_RESTART = (
    "metaspades.py -o {out} --restart-from last -t {threads} -m {mem} --tmp-dir {tmp}"
)
BWA_INDEX = "bwa index {assembly}"
BWA_MEM = "bwa mem -t {threads} {assembly} {r1} {r2}"
MEGAHIT_LEFTOVERS = "megahit -r {reads} -o {out} -t {threads} -m {mem} --tmp-dir {tmp}"
MEGAHIT_PAIRED = "megahit -1 {r1} -2 {r2} -o {out} --tmp-dir {tmp} -t {threads} -m {mem} --continue"
QUAST = "quast -t {threads} -o {out} -m {quast_min} {assembly}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MEM_GB = 24
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are dropped
METASPADES_KEEP_LEN = 1500  # scaffolds kept for the leftover-read remapping step
QUAST_MIN_CONTIG = 500  # QUAST -m; ignore contigs shorter than this in the report
MEGAHIT_MEM_BYTES = 1_000_000_000  # MEGAHIT -m takes bytes, so memory_gb * this
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 assembly",
        usage="metawrap2 assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir",
        description="Assemble paired-end reads with MEGAHIT and/or metaSPAdes (.gz accepted).",
    )
    p.add_argument("-1", dest="reads_1", required=True, help="forward fastq reads")
    p.add_argument("-2", dest="reads_2", required=True, help="reverse fastq reads")
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-m", "--memory", type=int, default=DEFAULT_MEM_GB, help="memory in GB (default 24)"
    )
    threads_arg(p)
    p.add_argument(
        "-l",
        "--min-len",
        type=int,
        default=MIN_CONTIG_LEN,
        help="minimum length of assembled contigs (default 1000)",
    )
    p.add_argument("--megahit", action="store_true", help="assemble with MEGAHIT (default)")
    p.add_argument(
        "--metaspades",
        action="store_true",
        help="assemble with metaSPAdes (with --megahit too, runs the hybrid pipeline)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _run_metaspades(args, env, spades_out, settings):
    announcement("ASSEMBLING WITH METASPADES")
    comm("Using reads %s and %s for assembly," % (args.reads_1, args.reads_2))
    tmp = scratch_dir(settings, args.output, "metaspades")
    if not dry_run() and os.path.isfile(os.path.join(spades_out, "spades.log")):
        run(
            METASPADES_RESTART.format(
                out=spades_out, threads=args.threads, mem=args.memory, tmp=tmp
            ),
            env=env,
            tool="metaspades.py",
        )
    else:
        run(
            METASPADES.format(
                tmp=tmp,
                threads=args.threads,
                mem=args.memory,
                out=spades_out,
                r1=args.reads_1,
                r2=args.reads_2,
            ),
            env=env,
            tool="metaspades.py",
        )
    expect_file(
        os.path.join(spades_out, "scaffolds.fasta"),
        "Something went wrong with metaSPAdes assembly. Exiting.",
    )
    shutil.rmtree(tmp, ignore_errors=True)


def _extract_unmapped(args, env, spades_out):
    """Keep long metaSPAdes scaffolds, then pull out reads that do not map back to them."""
    announcement("SORTING OUT UNASSEMBLED READS")
    scaffolds = os.path.join(spades_out, "scaffolds.fasta")
    long_scaffolds = os.path.join(spades_out, "long_scaffolds.fasta")
    if not dry_run():
        with open(long_scaffolds, "w") as out:
            rm_short_contigs.remove_short_contigs(METASPADES_KEEP_LEN, scaffolds, out)
        if not os.path.getsize(long_scaffolds):
            error(
                "metaSPAdes produced no contigs over 1.5kb, so the rest of the pipeline wont "
                "work. Exiting."
            )

    run(BWA_INDEX.format(assembly=long_scaffolds), env=env, tool="bwa index")
    aln_sam = os.path.join(args.output, "aligned_to_metaspades.sam")
    unused_sam = os.path.join(args.output, "unused_by_metaspades.sam")
    run(
        BWA_MEM.format(
            threads=args.threads, assembly=long_scaffolds, r1=args.reads_1, r2=args.reads_2
        ),
        env=env,
        tool="bwa mem",
        log_path=aln_sam,
    )
    unused_fastq = os.path.join(args.output, "unused_by_metaspades.fastq")
    if dry_run():
        return
    # keep only reads that did NOT map (aligned records carry the NM:i: edit-distance tag)
    with open(aln_sam) as inp, open(unused_sam, "w") as out:
        for line in inp:
            if "NM:i:" not in line:
                out.write(line)
    with open(unused_fastq, "w") as out:
        sam_to_fastq.sam_to_fastq(unused_sam, out)
    os.remove(aln_sam)
    os.remove(unused_sam)
    if not os.path.getsize(unused_fastq):
        error("Something went wrong with pulling out unassembled reads. Exiting.")


def _run_megahit(args, env, megahit_out, spades_out, hybrid, settings):
    announcement("ASSEMBLING READS WITH MEGAHIT")
    tmp = scratch_dir(settings, args.output, "megahit")
    mem = args.memory * MEGAHIT_MEM_BYTES
    if hybrid:
        unused_fastq = os.path.join(args.output, "unused_by_metaspades.fastq")
        comm("assembling %s with megahit" % unused_fastq)
        run(
            MEGAHIT_LEFTOVERS.format(
                reads=unused_fastq, out=megahit_out, threads=args.threads, mem=mem, tmp=tmp
            ),
            env=env,
            tool="megahit",
        )
        fs_move(unused_fastq, os.path.join(spades_out, "unused_by_metaspades.fastq"))
    else:
        comm("assembling %s and %s with megahit" % (args.reads_1, args.reads_2))
        run(
            MEGAHIT_PAIRED.format(
                r1=args.reads_1,
                r2=args.reads_2,
                out=megahit_out,
                tmp=tmp,
                threads=args.threads,
                mem=mem,
            ),
            env=env,
            tool="megahit",
        )
    expect_file(
        os.path.join(megahit_out, "final.contigs.fa"),
        "Something went wrong with reassembling with Megahit. Exiting.",
    )
    shutil.rmtree(tmp, ignore_errors=True)


def _format_assembly(args, megahit_out, spades_out, do_megahit, do_metaspades):
    announcement("FORMAT THE ASSEMBLY")
    final_assembly = os.path.join(args.output, "final_assembly.fasta")
    if dry_run():
        comm("(dry run) would write the formatted assembly to %s" % final_assembly)
        return final_assembly
    if do_megahit and do_metaspades:
        comm("Reads were assembled with metaspades AND megahit")
        long_contigs = os.path.join(megahit_out, "long.contigs.fa")
        with open(long_contigs, "w") as out:
            fix_megahit_contig_naming.fix_naming(
                os.path.join(megahit_out, "final.contigs.fa"), args.min_len, out
            )
        combined = os.path.join(args.output, "combined_assembly.fasta")
        shutil.copyfile(os.path.join(spades_out, "long_scaffolds.fasta"), combined)
        with open(combined, "a") as out, open(long_contigs) as inp:
            shutil.copyfileobj(inp, out)
        with open(final_assembly, "w") as out:
            sort_contigs.sort_contigs(combined, out)
        os.remove(combined)
    elif do_metaspades:
        comm("Reads were assembled with metaspades")
        with open(final_assembly, "w") as out:
            rm_short_contigs.remove_short_contigs(
                args.min_len, os.path.join(spades_out, "scaffolds.fasta"), out
            )
        if not os.path.getsize(final_assembly):
            error("metaspades failed to produce long contigs (>500bp). Exiting.")
    elif do_megahit:
        comm("Reads were assembled with megahit. Formatting and sorting the assembly...")
        with open(final_assembly, "w") as out:
            fix_megahit_contig_naming.fix_naming(
                os.path.join(megahit_out, "final.contigs.fa"), args.min_len, out
            )
        if not os.path.getsize(final_assembly):
            error("megahit failed to produce long (>%d) contigs. Exiting." % args.min_len)
    else:
        error("Reads were NOT assembled with metaspades OR megahit...")

    if not os.path.getsize(final_assembly):
        error("Something went wrong with joining the assemblies. Exiting.")
    return final_assembly


def main(argv: List[str]) -> int:
    args = _parse_args(argv)

    do_metaspades = args.metaspades
    do_megahit = args.megahit or not args.metaspades  # default to MEGAHIT if nothing was picked
    hybrid = do_megahit and do_metaspades

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "assembly",
        [
            (check_fastq, args.reads_1, "forward reads (-1)"),
            (check_fastq, args.reads_2, "reverse reads (-2)"),
        ],
    )
    # Assembly is the step most likely to fill a filesystem, and the one where doing so wastes
    # the most time, so check before starting rather than several hours in.
    check_disk_space("assembly", settings, [args.reads_1, args.reads_2], args.output)
    env = env_for("assembly", settings)
    rec = start_run("assembly", args, env, settings, inputs=[args.reads_1, args.reads_2])

    ckpt = make_checkpoint(args.output)
    try:

        if os.path.isdir(args.output):
            warning("%s already exists." % args.output)
        ensure_dir(args.output)

        spades_out = os.path.join(args.output, "metaspades")
        megahit_out = os.path.join(args.output, "megahit")

        if do_metaspades and ckpt.todo("metaspades"):
            _run_metaspades(args, env, spades_out, settings)
            ckpt.done("metaspades")
        if hybrid and ckpt.todo("extract_leftovers"):
            _extract_unmapped(args, env, spades_out)
            ckpt.done("extract_leftovers")
        if do_megahit and ckpt.todo("megahit"):
            _run_megahit(args, env, megahit_out, spades_out, hybrid, settings)
            ckpt.done("megahit")

        if ckpt.todo("format"):
            _format_assembly(args, megahit_out, spades_out, do_megahit, do_metaspades)
            ckpt.done("format")
        else:
            comm("skipping assembly formatting (already done; --resume)")

        if ckpt.todo("quast"):
            announcement("RUNNING ASSEMBLY QC WITH QUAST")
            quast_out = os.path.join(args.output, "QUAST_out")
            run(
                QUAST.format(
                    threads=args.threads,
                    out=quast_out,
                    quast_min=QUAST_MIN_CONTIG,
                    assembly=os.path.join(args.output, "final_assembly.fasta"),
                ),
                env=env,
                tool="quast",
            )
            expect_produced(os.path.join(args.output, "final_assembly.fasta"), "final assembly")
            report = os.path.join(args.output, "assembly_report.html")
            if not dry_run():
                shutil.copyfile(os.path.join(quast_out, "report.html"), report)
                if not os.path.getsize(report):
                    error("Something went wrong with running QUAST. Exiting.")
            ckpt.done("quast")
        else:
            comm("skipping QUAST report (already done; --resume)")

        announcement("ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
