"""MetaWrap2 assembly module: produce a metagenomic assembly from paired-end reads.

Assembles QC'd paired-end reads with MEGAHIT and/or metaSPAdes. When both are selected the
reads are first assembled with metaSPAdes, the reads that do not map back to the long
scaffolds are pulled out, and those leftovers are assembled with MEGAHIT (which does better
on low-coverage data); the two assemblies are then combined, sorted, and short contigs are
removed. The final assembly is QC'd with QUAST.

The hybrid path's remapping step only reads whether each read mapped, so the aligner acts as a
filter. `--aligner` picks it from the registry `binning` and `blobology` use
(`metawrap2.aligners`). The default here is bwa, not the registry's minimap2; see
`HYBRID_ALIGNER`.

`--keep-graph` puts each assembler's assembly graph at a fixed, named path in the output
directory, for plasmid-aware or graph-aware binners. Without it the graphs stay inside the
assemblers' own directories under different names, and MEGAHIT writes none. It is off by default
and deletes nothing.

Commands are in the COMMANDS block below and tunable values in the CONSTANTS block.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from typing import Dict, List, Optional

from .. import longreads as _long
from ..aligners import ALIGNERS, INDEX_SUFFIXES, index_path
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..limits import human
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
    listdir,
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
# The aligner streams into the FASTQ conversion, so the alignment is never written to disk.
# pipefail is required because the converter exits 0 on a truncated stream: an aligner that ran out
# of memory would otherwise leave a short FASTQ for the hybrid assembly. PY is the host interpreter
# by absolute path, which resolves even though the left side of the pipe runs in the conda env.
#
# {align} comes from metawrap2.aligners, shared by every module that maps reads to an assembly.
EXTRACT_UNMAPPED = (
    "set -o pipefail; {align}"
    " | " + PY + " -m metawrap2.scripts.sam_to_fastq --unmapped - > {fastq}"
)
MEGAHIT_LEFTOVERS = "megahit -r {reads} -o {out} -t {threads} -m {mem} --tmp-dir {tmp}"
MEGAHIT_PAIRED = "megahit -1 {r1} -2 {r2} -o {out} --tmp-dir {tmp} -t {threads} -m {mem} --continue"
# MEGAHIT writes no graph: the de Bruijn graph exists only as the per-k contigs in
# intermediate_contigs/, which MEGAHIT's toolkit converts to FASTG. On an oral sample (117 MB
# final.contigs.fa), k141.contigs.fa converted to a 235 MB fastg in 1.8 s on one core, so it is done
# on demand (--keep-graph).
#
# The graph is written to stdout, so this command uses stdout_path rather than -o.
MEGAHIT_CONTIG2FASTG = "megahit_toolkit contig2fastg {k} {contigs}"
QUAST = "quast -t {threads} -o {out} -m {quast_min} {assembly}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are dropped
METASPADES_KEEP_LEN = 1500  # scaffolds kept for the leftover-read remapping step
QUAST_MIN_CONTIG = 500  # QUAST -m; ignore contigs shorter than this in the report
MEGAHIT_MEM_BYTES = 1_000_000_000  # MEGAHIT -m takes bytes, so memory_gb * this
# --aligner default. Differs from aligners.DEFAULT_ALIGNER (minimap2) because this stage reads only
# the unmapped bit, and the reads that carry it are what MEGAHIT assembles. For binning depth,
# minimap2 and bwa agree closely (Spearman 0.984).
#
# Measured on the 273k-scaffold metaSPAdes assembly from the `assembler` benchmark arm, filtered to
# METASPADES_KEEP_LEN, with the 3-sample co-assembly reads it was built from:
#
#   bwa       5.45 CPU-hours, 4.5 GB peak, recruits 2,867,680 reads as unassembled
#   minimap2  0.62 CPU-hours, 1.7 GB peak, recruits 4,113,486 - 44% more
#
# minimap2 calls 1,249,798 reads unmapped that bwa placed on a scaffold; bwa calls 3,992 unmapped
# that minimap2 placed (Jaccard 0.695). minimap2's `sr` preset is less sensitive on a fragmented
# assembly, so MEGAHIT reassembles sequence metaSPAdes already assembled and _format_assembly keeps
# both copies. `--aligner minimap2` is available for the 8.8x lower cost.
HYBRID_ALIGNER = "bwa"
# metaFlye for long reads. --meta changes Flye's repeat-resolution and coverage assumptions, which
# by default expect one genome at even coverage.
#
# The platform flag comes from metawrap2.longreads, which maps platform names for Flye, minimap2
# and SemiBin2.
#
# --keep-haplotypes is not passed: it keeps strain variants as separate contigs, which binning puts
# in the same bin and CheckM2 counts as contamination.
FLYE = "flye {flag} {reads} --meta --out-dir {out} --threads {threads}"
#: Number of Flye polishing rounds. Flye's default of one is kept; extra rounds cost time and
#: converge quickly.
FLYE_ITERATIONS = 1
# --keep-graph: the graph each assembler writes, and the name it is kept under.
#
# From metaSPAdes, assembly_graph_with_scaffolds.gfa: GFA with scaffold paths is what graph-aware
# binners (PlasMAAG, GraphMB) read, and it is smaller (121 MB against 231 MB for
# assembly_graph.fastg, on a gut co-assembly with a 54 MB final_assembly.fasta). metaFlye's
# assembly_graph.gfa was 2.3 MB for a 2.3 MB assembly; its assembly_graph.gv is a 21 KB Graphviz
# rendering that no binner reads, so it is not kept.
#
# Names are per assembler because `-L` with `-1/-2` runs two assemblers in one output directory,
# and metaSPAdes and metaFlye both name their graph assembly_graph.gfa.
SPADES_GRAPH = "assembly_graph_with_scaffolds.gfa"
FLYE_GRAPH = "assembly_graph.gfa"
MEGAHIT_INTERMEDIATES = "intermediate_contigs"
KEPT_GRAPHS: Dict[str, str] = {
    "megahit": "assembly_graph.megahit.fastg",
    "metaspades": "assembly_graph.metaspades.gfa",
    "flye": "assembly_graph.flye.gfa",
}
#: ``k<k>.contigs.fa`` in MEGAHIT's intermediates, the files the FASTG is built from.
MEGAHIT_K_CONTIGS = re.compile(r"^k(\d+)\.contigs\.fa$")
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 assembly",
        usage="metawrap2 assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir",
        description="Assemble paired-end reads with MEGAHIT and/or metaSPAdes (.gz accepted).",
    )
    # Optional: a long-read assembly has no short reads.
    p.add_argument("-1", dest="reads_1", help="forward fastq reads (short-read assembly)")
    p.add_argument("-2", dest="reads_2", help="reverse fastq reads (short-read assembly)")
    p.add_argument(
        "-L",
        "--long-reads",
        nargs="+",
        metavar="FILE",
        help="long-read fastq/fasta (gzip accepted), assembled with metaFlye. Several files are "
        "assembled together as one dataset (a co-assembly); there is no need to concatenate them",
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
    p.add_argument(
        "--keep-graph",
        action="store_true",
        help="keep each assembler's assembly graph as a named file in the output directory, for "
        "tools that need a graph (plasmid-aware and graph-aware binners, Hi-C placement). "
        "metaSPAdes' and metaFlye's GFA files are hard-linked and use no extra disk (121 MB and "
        "2.3 MB measured). MEGAHIT writes no graph, so one is built with megahit_toolkit "
        "contig2fastg from MEGAHIT's intermediate_contigs: 235 MB and 2 seconds of one core for a "
        "117 MB assembly. Off by default. Nothing is deleted either way, and --resume "
        "--keep-graph adds the graph to a finished assembly",
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
                "metaSPAdes produced no contigs over 1.5kb, so the rest of the pipeline cannot "
                "run. Exiting."
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
    # Remove the index once the FASTQ exists: nothing else reads the alignment, and a --resume
    # that reaches this checkpointed stage rebuilds it. minimap2's .mmi for a 176 MB scaffold set is
    # larger than the fasta.
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
        comm("%d read files, assembled together by Flye as one dataset" % len(reads))

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


def _link_or_copy(src: str, dest: str) -> None:
    """Put *src* at *dest* as a hard link, falling back to a copy.

    The graphs are large (121 MB for the metaSPAdes GFA of a 54 MB assembly), and a link avoids
    storing them twice while keeping the graph if the assembler's directory is cleared. ``os.link``
    fails across filesystems and on some network mounts, hence the fallback; *src* has already
    been checked to exist.
    """
    if os.path.lexists(dest):
        os.remove(dest)
    try:
        os.link(src, dest)
    except OSError:
        shutil.copyfile(src, dest)


def _keep_existing_graph(src: str, dest: str, assembler: str) -> List[str]:
    """Keep the graph *assembler* already wrote at *src* as *dest*. Returns what was kept."""
    if dry_run():
        comm("(dry run) would keep %s's assembly graph as %s" % (assembler, dest))
        return [dest]
    if not os.path.isfile(src) or not os.path.getsize(src):
        # A warning, not an error: the assembly is complete and the graph is optional, with no
        # downstream consumer yet.
        warning(
            "%s wrote no assembly graph at %s, so there is none to keep. The assembly itself is "
            "unaffected." % (assembler, src)
        )
        return []
    _link_or_copy(src, dest)
    comm("%s assembly graph -> %s (%s)" % (assembler, dest, human(os.path.getsize(dest))))
    return [dest]


def _megahit_graph_k(intermediates: str) -> Optional[int]:
    """The largest k MEGAHIT wrote contigs for, or None if it wrote none.

    ``options.json``'s ``k_max`` is not used: MEGAHIT's auto-k stops when a larger k adds nothing
    (a test leftover assembly ran k21..k59 with k_max 141). final.contigs.fa comes from the largest
    k present.
    """
    ks = [
        int(match.group(1))
        for match in (MEGAHIT_K_CONTIGS.match(name) for name in listdir(intermediates))
        if match and os.path.getsize(os.path.join(intermediates, match.group(0)))
    ]
    return max(ks) if ks else None


def _build_megahit_graph(
    args: argparse.Namespace, env: Optional[str], megahit_out: str
) -> List[str]:
    """Convert MEGAHIT's largest-k contigs into a FASTG. Returns what was kept."""
    dest = os.path.join(args.output, KEPT_GRAPHS["megahit"])
    intermediates = os.path.join(megahit_out, MEGAHIT_INTERMEDIATES)
    k = _megahit_graph_k(intermediates)
    if k is None and not dry_run():
        warning(
            "MEGAHIT left no k*.contigs.fa in %s, so its assembly graph cannot be built (MEGAHIT "
            "writes no graph itself). The assembly itself is unaffected." % intermediates
        )
        return []
    # Under --dry-run nothing has been assembled, so print a placeholder for k (as binning._bams
    # shows <sample>.bam).
    label = "<kmax>" if k is None else str(k)
    rc = run(
        MEGAHIT_CONTIG2FASTG.format(
            k=label, contigs=os.path.join(intermediates, "k%s.contigs.fa" % label)
        ),
        env=env,
        tool="megahit_toolkit contig2fastg",
        stdout_path=dest,
        check=False,
    )
    if dry_run():
        comm("(dry run) would build MEGAHIT's assembly graph as %s" % dest)
        return [dest]
    # stdout_path creates the file before the tool runs, so an empty file means the conversion
    # failed. Do not leave an empty graph where downstream tools would read it.
    if rc != 0 or not os.path.getsize(dest):
        os.remove(dest)
        warning(
            "megahit_toolkit contig2fastg failed on k%s, so no MEGAHIT graph was kept. The "
            "assembly itself is unaffected." % label
        )
        return []
    comm("MEGAHIT assembly graph (k%s) -> %s (%s)" % (label, dest, human(os.path.getsize(dest))))
    return [dest]


def _keep_graphs(
    args: argparse.Namespace,
    env: Optional[str],
    megahit_out: str,
    spades_out: str,
    flye_out: str,
    do_megahit: bool,
    do_metaspades: bool,
    do_flye: bool,
) -> List[str]:
    """Put every assembler's graph at a fixed path in the output directory, and report sizes.

    A separate stage, so ``--resume --keep-graph`` can add graphs to a finished assembly:
    metaSPAdes' and metaFlye's graphs are still in place and MEGAHIT's intermediates can still be
    converted.
    """
    announcement("KEEPING THE ASSEMBLY GRAPHS")
    kept: List[str] = []
    if do_metaspades:
        kept += _keep_existing_graph(
            os.path.join(spades_out, SPADES_GRAPH),
            os.path.join(args.output, KEPT_GRAPHS["metaspades"]),
            "metaSPAdes",
        )
    if do_flye:
        kept += _keep_existing_graph(
            os.path.join(flye_out, FLYE_GRAPH),
            os.path.join(args.output, KEPT_GRAPHS["flye"]),
            "metaFlye",
        )
    if do_megahit:
        kept += _build_megahit_graph(args, env, megahit_out)
    if not kept:
        warning("--keep-graph was given but no assembler left a graph to keep.")
    elif not dry_run():
        comm(
            "%d assembly graph(s) kept, %s in total"
            % (len(kept), human(sum(os.path.getsize(p) for p in kept)))
        )
    return kept


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
        # Apply the same length filter and sort as the other paths, so a long-read assembly has the
        # same shape as a short-read one. Flye's contig names (contig_1, contig_2) need no special
        # handling downstream.
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
            "long and short reads cannot be assembled together in one run: a merged metagenome "
            "assembly cannot tell a duplicate contig from a new one. Assemble them in separate "
            "runs and bin each assembly"
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
    # MEGAHIT is the default. On three CAMI2 gut samples co-assembled (2026-09-30), metaSPAdes gave
    # a 38% longer N50 (12.9 kb vs 9.3 kb) for 3.2x the CPU-hours (245 vs 77) and 9.7x the peak RSS
    # (140 GB vs 15 GB), and assembled 14% fewer total bases. --metaspades is available when the
    # memory is. With only long reads, no short-read assembler runs.
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
    # Assembly is the step most likely to fill a filesystem, so check space before starting.
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

        # Runs only with --keep-graph; without it the assemblers' output is left as is.
        if args.keep_graph and ckpt.todo("keep_graph"):
            _keep_graphs(
                args,
                env,
                megahit_out,
                spades_out,
                flye_out,
                do_megahit,
                do_metaspades,
                do_flye,
            )
            ckpt.done("keep_graph")

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
