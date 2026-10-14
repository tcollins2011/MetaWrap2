"""MetaWrap2 assembly module: produce a metagenomic assembly from paired-end reads.

Assembles QC'd paired-end reads with MEGAHIT and/or metaSPAdes. When both are selected the
reads are first assembled with metaSPAdes, the reads that do not map back to the long
scaffolds are pulled out, and those leftovers are assembled with MEGAHIT (which does better
on low-coverage data); the two assemblies are then combined, sorted, and short contigs are
removed. The final assembly is QC'd with QUAST.

The remapping step in that hybrid path reads one bit out of the alignment - whether a read mapped -
so the aligner is a filter here rather than a measuring instrument, and `--aligner` picks it from the
same registry `binning` and `blobology` use (`metawrap2.aligners`). Its default is bwa, not the
registry's minimap2, and `HYBRID_ALIGNER` says why: the 5.7x saving was measured on depth agreement,
which is not the property this stage depends on.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers/defaults in the CONSTANTS block. Change a flag there
and it takes effect — no need to follow the orchestration logic underneath.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
from typing import List

from .. import longreads as _long
from ..aligners import ALIGNERS, INDEX_SUFFIXES, index_path
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..pyrun import PY
from ..scripts import fix_megahit_contig_naming, rm_short_contigs, sort_contigs
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
    memory_arg,
    require_file,
    resolve_memory_gb,
    resolve_threads,
    run,
    scratch_dir,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["assembly"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
METASPADES = "metaspades.py --tmp-dir {tmp} -t {threads} -m {mem} -o {out} -1 {r1} -2 {r2}"
METASPADES_RESTART = (
    "metaspades.py -o {out} --restart-from last -t {threads} -m {mem} --tmp-dir {tmp}"
)
# The aligner streams into the FASTQ conversion, so the alignment never lands on disk. pipefail is
# what makes that safe: the converter exits 0 on a truncated stream, so without it an aligner that
# ran out of memory would leave a short FASTQ and the hybrid assembly would quietly be built on a
# slice of the unassembled reads. PY is the host interpreter by absolute path, which resolves even
# though the left-hand side of the pipe runs inside this module's conda env.
#
# {align} comes from metawrap2.aligners rather than being spelled out here, so the three modules that
# map reads back to an assembly cannot drift apart on which aligners they offer.
EXTRACT_UNMAPPED = (
    "set -o pipefail; {align}"
    " | " + PY + " -m metawrap2.scripts.sam_to_fastq --unmapped - > {fastq}"
)
MEGAHIT_LEFTOVERS = "megahit -r {reads} -o {out} -t {threads} -m {mem} --tmp-dir {tmp}"
MEGAHIT_PAIRED = "megahit -1 {r1} -2 {r2} -o {out} --tmp-dir {tmp} -t {threads} -m {mem} --continue"
QUAST = "quast -t {threads} -o {out} -m {quast_min} {assembly}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are dropped
METASPADES_KEEP_LEN = 1500  # scaffolds kept for the leftover-read remapping step
QUAST_MIN_CONTIG = 500  # QUAST -m; ignore contigs shorter than this in the report
MEGAHIT_MEM_BYTES = 1_000_000_000  # MEGAHIT -m takes bytes, so memory_gb * this
# --aligner default, and deliberately not aligners.DEFAULT_ALIGNER, which is minimap2.
#
# The registry's evidence does not reach this stage, and measuring it here reversed the answer. There
# the aligner is an implementation detail of a coverage estimator and the output is a depth vector,
# which minimap2 and bwa agree on closely enough (Spearman 0.984) that the bins come out the same.
# Here it is a filter: the only thing read out of the alignment is the unmapped bit, and whichever
# reads carry it are what MEGAHIT assembles.
#
# Measured on the real thing - the 273k-scaffold metaSPAdes assembly from the `assembler` benchmark
# arm, filtered to METASPADES_KEEP_LEN, and the 3-sample co-assembly reads it was built from:
#
#   bwa       5.45 CPU-hours, 4.5 GB peak, recruits 2,867,680 reads as unassembled
#   minimap2  0.62 CPU-hours, 1.7 GB peak, recruits 4,113,486 - 44% more
#
# 8.8x cheaper, and not a drop-in. The disagreement is one-sided: minimap2 calls 1,249,798 reads
# unmapped that bwa placed on a scaffold, while bwa calls only 3,992 unmapped that minimap2 placed
# (Jaccard 0.695). So minimap2's `sr` preset is simply less sensitive on a fragmented real assembly,
# and the reads it misses are handed to MEGAHIT as "unassembled" - which means MEGAHIT reassembles
# sequence metaSPAdes already assembled, and _format_assembly concatenates both copies without being
# able to tell. The cheaper aligner changes what the assembly contains, in the bad direction, which
# is the one outcome a cost saving cannot buy. `--aligner minimap2` is there for anyone who wants the
# 8.8x knowing that.
#
# Note this is the opposite conclusion to the binning registry's, from the same two tools, because
# the property read out of the alignment is different. That is the argument for measuring per
# question rather than per tool.
HYBRID_ALIGNER = "bwa"
# metaFlye for long reads. --meta is the whole point: Flye's default assumes one genome at even
# coverage, and a metagenome is neither - --meta changes the repeat-resolution and coverage
# assumptions rather than just relaxing a threshold.
#
# The platform flag comes from metawrap2.longreads, because Flye, minimap2 and SemiBin2 each name
# the same four platforms differently and one place to get that wrong is enough.
#
# --keep-haplotypes is deliberately *not* passed. It keeps alternative haplotypes as separate
# contigs, which for a metagenome means one strain's variants arriving as extra contigs that
# binning will place in the same bin and CheckM2 will read as contamination.
FLYE = "flye {flag} {reads} --meta --out-dir {out} --threads {threads}"
#: Long-read assembly is the one path here where the assembler also polishes. Flye's default is one
#: round; more costs time per round and converges quickly, so this is left at the default and
#: recorded rather than tuned.
FLYE_ITERATIONS = 1
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 assembly",
        usage="metawrap2 assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir",
        description="Assemble paired-end reads with MEGAHIT and/or metaSPAdes (.gz accepted).",
    )
    # Not required any more: a long-read assembly has no short reads at all, and requiring them
    # would mean passing /dev/null to get metaFlye.
    p.add_argument("-1", dest="reads_1", help="forward fastq reads (short-read assembly)")
    p.add_argument("-2", dest="reads_2", help="reverse fastq reads (short-read assembly)")
    p.add_argument(
        "-L",
        "--long-reads",
        nargs="+",
        metavar="FILE",
        help="long-read fastq/fasta (gzip accepted), assembled with metaFlye. Several files are "
        "one dataset to Flye, so a co-assembly is this flag rather than a concatenation step",
    )
    p.add_argument(
        "--read-type",
        choices=sorted(_long.READ_TYPES),
        default=_long.DEFAULT_READ_TYPE,
        help="long-read platform (default: %(default)s). It sets Flye's error model and the "
        "minimap2 preset used downstream; the run checks the read lengths against it",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    memory_arg(p)
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
        help="assemble with metaSPAdes: longer contigs, but measured at 3x the CPU-hours and "
        "~140 GB peak memory on a 3-sample co-assembly (with --megahit too, runs the hybrid "
        "pipeline)",
    )
    p.add_argument(
        "--aligner",
        choices=sorted(ALIGNERS),
        default=HYBRID_ALIGNER,
        help="aligner used to find the reads metaSPAdes could not place, in the hybrid "
        "(--megahit --metaspades) path only (default %s)" % HYBRID_ALIGNER,
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

    aligner = ALIGNERS[args.aligner]
    index = index_path(long_scaffolds, args.aligner)
    run(
        aligner["index_cmd"].format(assembly=long_scaffolds, index=index),
        env=env,
        tool="%s index" % args.aligner,
    )
    unused_fastq = os.path.join(args.output, "unused_by_metaspades.fastq")
    run(
        [
            "bash",
            "-c",
            EXTRACT_UNMAPPED.format(
                align=aligner["paired"].format(
                    threads=args.threads, index=index, r1=args.reads_1, r2=args.reads_2
                ),
                fastq=unused_fastq,
            ),
        ],
        env=env,
        tool="%s | sam_to_fastq" % args.aligner,
    )
    if dry_run():
        return
    if not os.path.getsize(unused_fastq):
        error("Something went wrong with pulling out unassembled reads. Exiting.")
    # The index is spent the moment the FASTQ exists: nothing else in the module reads an alignment,
    # and the stage is checkpointed, so a --resume that reaches here rebuilds it anyway. Worth
    # sweeping rather than leaving, because it is not small next to what it was built from -
    # minimap2's .mmi for a 176 MB scaffold set is larger than the fasta.
    for suffix in INDEX_SUFFIXES:
        spent = long_scaffolds + suffix
        if os.path.exists(spent):
            os.remove(spent)


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


def _run_flye(args, env, flye_out, settings):
    """Assemble the long reads with metaFlye, after checking they are what the user said they are."""
    announcement("ASSEMBLING LONG READS WITH METAFLYE")
    spec = _long.resolve(args.read_type)
    comm("%s: %s" % (spec.name, spec.description))
    reads = _long.coassembly_inputs(args.long_reads)
    if not reads:
        error("none of the long-read files given exist")
    if len(reads) > 1:
        comm(
            "%d read files, assembled together as one dataset - Flye takes several and treats them "
            "as one, so a co-assembly needs no concatenation step" % len(reads)
        )

    rows = []
    for path in reads:
        stats = _long.scan_reads(path)
        notes = _long.check_type(stats, args.read_type)
        comm(
            "%s: %d reads scanned, N50 %d bp, median %d bp, %.0f%% over 10 kb"
            % (
                os.path.basename(path),
                stats.reads,
                stats.n50,
                stats.median,
                100.0 * stats.fraction_over_10kb,
            )
        )
        for note in notes:
            warning(note)
        rows.append(_long.read_report_row(os.path.basename(path), args.read_type, stats, notes))

    report = os.path.join(args.output, "long_read_report" + TSV_SUFFIX)
    if not dry_run():
        ensure_dir(args.output)
        with open(report, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_long.READ_REPORT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        comm("-> %s" % report)

    run(
        FLYE.format(
            flag=spec.flye_flag,
            reads=" ".join(reads),
            out=flye_out,
            threads=args.threads,
        ),
        env=env,
        tool="flye",
    )
    expect_file(
        os.path.join(flye_out, "assembly.fasta"),
        "metaFlye produced no assembly. Its own log is %s/flye.log" % flye_out,
    )
    comm(_long.LONG_READ_CAVEAT)
    if args.reads_1:
        comm(_long.HYBRID_CAVEAT)


def _format_assembly(
    args, megahit_out, spades_out, do_megahit, do_metaspades, flye_out=None, do_flye=False
):
    announcement("FORMAT THE ASSEMBLY")
    final_assembly = os.path.join(args.output, "final_assembly.fasta")
    if dry_run():
        comm("(dry run) would write the formatted assembly to %s" % final_assembly)
        return final_assembly
    if do_flye and not (do_megahit or do_metaspades):
        comm("Reads were assembled with metaFlye. Filtering and sorting the assembly...")
        # The same length filter and sort the other paths apply, so a long-read assembly is the
        # same shape of input to binning as a short-read one - including contig names, which Flye
        # writes as contig_1, contig_2 and which nothing downstream has to special-case.
        filtered = os.path.join(args.output, "flye_filtered.fasta")
        with open(filtered, "w") as out:
            rm_short_contigs.remove_short_contigs(
                args.min_len, os.path.join(flye_out, "assembly.fasta"), out
            )
        with open(final_assembly, "w") as out:
            sort_contigs.sort_contigs(filtered, out)
        os.remove(filtered)
        if not os.path.getsize(final_assembly):
            error("metaFlye produced no contigs above %d bp. Exiting." % args.min_len)
        return final_assembly
    if do_flye:
        error(
            "a long-read assembly and a short-read assembly are two assemblies, not one: this "
            "module will not merge them, because a merged metagenome assembly cannot tell a "
            "duplicate contig from a novel one (the same reason the hybrid short-read path is a "
            "measured special case). Run them separately and bin each"
        )
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

    do_flye = bool(args.long_reads)
    do_metaspades = args.metaspades
    # MEGAHIT when nothing was picked, and the benchmark says to keep it that way. On three CAMI2
    # gut samples co-assembled (benchmarks/results/cami2_gi.tsv, 2026-09-30) metaSPAdes returned a
    # 38% longer N50 (12.9 kb vs 9.3 kb) for 3.2x the CPU-hours (245 vs 77) and 9.7x the peak RSS
    # (140 GB vs 15 GB), while assembling 14% *fewer* total bases. The contiguity is real and worth
    # having when you can afford it - hence --metaspades - but 140 GB of resident memory is more
    # than most machines can give, so it cannot be what happens to someone who passed no flag.
    # Nothing short-read runs when only long reads were given: MEGAHIT is the fallback for "no
    # assembler chosen", not for "no short reads".
    do_megahit = (args.megahit or not args.metaspades) and bool(args.reads_1)
    hybrid = do_megahit and do_metaspades

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args)
    if not args.long_reads and not (args.reads_1 and args.reads_2):
        error(
            "nothing to assemble: pass -1/-2 for short reads, -L for long reads, or both. A "
            "short-read assembly needs both mates"
        )
    if args.reads_1 and not args.reads_2:
        error("short-read assembly needs both mates: -1 was given without -2")
    checks = []
    if args.reads_1:
        checks.append((check_fastq, args.reads_1, "forward reads (-1)"))
        checks.append((check_fastq, args.reads_2, "reverse reads (-2)"))
    for path in args.long_reads or []:
        checks.append((require_file, path, "long reads (-L)"))
    validate_inputs("assembly", checks)
    # Assembly is the step most likely to fill a filesystem, and the one where doing so wastes
    # the most time, so check before starting rather than several hours in.
    check_disk_space(
        "assembly",
        settings,
        [path for path in [args.reads_1, args.reads_2, *(args.long_reads or [])] if path],
        args.output,
    )
    env = env_for("assembly", settings)
    rec = start_run(
        "assembly",
        args,
        env,
        settings,
        inputs=[p for p in [args.reads_1, args.reads_2, *(args.long_reads or [])] if p],
    )

    ckpt = make_checkpoint(args.output)
    try:

        if os.path.isdir(args.output):
            warning("%s already exists." % args.output)
        ensure_dir(args.output)

        spades_out = os.path.join(args.output, "metaspades")
        megahit_out = os.path.join(args.output, "megahit")
        flye_out = os.path.join(args.output, "flye")

        if do_flye and ckpt.todo("flye"):
            _run_flye(args, env, flye_out, settings)
            ckpt.done("flye")

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
            _format_assembly(
                args, megahit_out, spades_out, do_megahit, do_metaspades, flye_out, do_flye
            )
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
