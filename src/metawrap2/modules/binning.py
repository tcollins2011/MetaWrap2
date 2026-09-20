"""MetaWrap2 binning module: split assembly contigs into bins with metaBAT2/MaxBin2/CONCOCT.

Reads back to the assembly are aligned to produce coverage, then each selected binner runs.
Run the bin_refinement module afterwards to QC, consolidate, and reassemble the bins.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import os
import re
from typing import IO, Dict, List, Tuple

from .. import checkm as _checkm
from ..config import load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, UNBINNED_NAME, bin_filename
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar
from ._common import (
    absolutize_paths,
    announcement,
    check_disk_space,
    check_fasta,
    check_fastq,
    collect_read_pairs,
    collect_reads,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    listdir,
    make_checkpoint,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
)

CONDA_ENV = "metawrap2-binning"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BWA_INDEX = "bwa index {assembly}"
BWA_MEM_PAIRED = "bwa mem -v 1 -t {threads} {assembly} {r1} {r2}"
BWA_MEM_SINGLE = "bwa mem -t {threads} {assembly} {reads}"
BWA_MEM_INTERLEAVED = "bwa mem -v 1 -p -t {threads} {assembly} {reads}"
SAMTOOLS_SORT = "samtools sort -T {tmp} -@ {threads} -O BAM -o {bam} {sam}"
SAMTOOLS_INDEX = "samtools index -@ {threads} -b {bam}"
JGI_DEPTH = "jgi_summarize_bam_contig_depths --outputDepth {depth} {bams}"
JGI_DEPTH_MAXBIN = (
    "jgi_summarize_bam_contig_depths --outputDepth {depth} --noIntraDepthVariance {bams}"
)
METABAT2 = "metabat2 -i {assembly} -a {depth} -o {out}/metabat2_bins/bin -m {metabat_len} -t {threads} --unbinned"
MAXBIN2 = "run_MaxBin.pl -contig {assembly} -markerset {markers} -thread {threads} -min_contig_length {min_len} -out {maxbin_out}/bin -abund_list {abund_list}"
CONCOCT_CUTUP = "cut_up_fasta.py {assembly} -c 10000 --merge_last -b {bed} -o 0"
CONCOCT_COVERAGE = "concoct_coverage_table.py {bed} {bams}"
CONCOCT = "concoct -l {min_len} -t {threads} --coverage_file {cov} --composition_file {comp} -b {concoct_out}"
CONCOCT_MERGE = "merge_cutup_clustering.py {clustering}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
MIN_CONTIG_LEN = 1000  # -l default; contigs shorter than this are not binned
METABAT_MIN_LEN = 1500  # metaBAT2 floors the minimum contig length at 1500 bp
DEFAULT_MEM_GB = 4
MAXBIN_MARKERS_BACTERIAL = 107
#: Files metaBAT2 writes into its output prefix that are not genomes. Left in the bin set,
#: CheckM and bin_refinement score them as if they were bins.
METABAT2_NON_BINS = ("bin.unbinned.fa", "bin.lowDepth.fa", "bin.tooShort.fa")
MAXBIN_MARKERS_UNIVERSAL = 40  # --universal (better Archaea binning)
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 binning",
        usage="metawrap2 binning [options] -a assembly.fa -o output_dir reads_1.fastq reads_2.fastq ...",
        description="Bin assembly contigs with metaBAT2, MaxBin2 and/or CONCOCT. "
        "Provide each replicate's reads separately (.gz accepted).",
    )
    p.add_argument("-a", "--assembly", required=True, help="metagenomic assembly fasta")
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument("-m", "--memory", type=int, default=DEFAULT_MEM_GB, help="RAM in GB (default 4)")
    p.add_argument(
        "-l",
        "--min-len",
        type=int,
        default=MIN_CONTIG_LEN,
        help="minimum contig length to bin (default 1000; metaBAT2 floors at 1500)",
    )
    p.add_argument("--metabat2", action="store_true", help="bin with metaBAT2")
    p.add_argument("--maxbin2", action="store_true", help="bin with MaxBin2")
    p.add_argument("--concoct", action="store_true", help="bin with CONCOCT")
    p.add_argument(
        "--universal", action="store_true", help="use universal (not bacterial) markers in MaxBin2"
    )
    p.add_argument(
        "--run-checkm",
        action="store_true",
        help="run CheckM on the bins (uses the bin_refinement env; needs plenty of RAM)",
    )
    p.add_argument("--single-end", action="store_true", help="reads are single-end (*.fastq)")
    p.add_argument("--interleaved", action="store_true", help="reads are interleaved paired-end")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument(
        "reads", nargs="+", help="read files (paired *_1.fastq/*_2.fastq, or single/interleaved)"
    )
    return p.parse_args(argv)


def _split_concoct_bins(clustering_csv: str, assembly: str, out_dir: str) -> None:
    """Split contigs into per-cluster bin files from a CONCOCT merged clustering CSV.

    ``merge_cutup_clustering.py`` already folds the ``<contig>.<fragment>`` ids that
    ``cut_up_fasta.py`` produced back into whole-contig ids, so the keys here are plain
    contig names and must be used verbatim. The original split_concoct_bins.py truncated
    each name at the *first* '.', which silently mangled any assembler whose contig names
    contain a dot - metaSPAdes emits ``NODE_1_length_5000_cov_3.5`` - collapsing unrelated
    contigs onto one key and scattering them into the wrong bins.
    """
    contig_bin = {}
    with open(clustering_csv) as fh:
        for line in fh:
            if line.startswith("contig_id"):
                continue
            cid, _, cluster = line.strip().partition(",")
            if cid:
                contig_bin[cid] = cluster
    ensure_dir(out_dir)
    handles: Dict[str, IO[str]] = {}
    try:
        for header, seq in bar(iter_fasta(assembly), desc="splitting CONCOCT bins", unit=" contig"):
            contig = header.split()[0]
            contig_cluster = contig_bin.get(contig)
            # CONCOCT numbers its clusters from 0; bins are named from 1, like every other
            # binner's, so the cluster id is shifted rather than used as the filename.
            target = (
                bin_filename(int(contig_cluster) + 1)
                if contig_cluster is not None
                else UNBINNED_NAME
            )
            if target not in handles:
                # One handle per bin, all closed in the finally below. Reopening per contig
                # instead would make this O(contigs) in syscalls.
                handles[target] = open(os.path.join(out_dir, target), "w")  # noqa: SIM115
            handles[target].write(">%s\n%s\n" % (contig, seq))
    finally:
        for handle in handles.values():
            handle.close()


def _align_reads(args, env, work, assembly_copy) -> None:
    """Align each read set back to the assembly, producing sorted BAMs in work/."""
    run(BWA_INDEX.format(assembly=assembly_copy), env=env, tool="bwa index")
    tmp = os.path.join(work, "tmp-samtools")

    if args.single_end or args.interleaved:
        sets = collect_reads(args.reads)
    else:
        sets = [(s, r1, r2) for s, r1, r2 in collect_read_pairs(args.reads)]

    for entry in sets:
        sample = entry[0]
        bam = os.path.join(work, sample + ".bam")
        if not dry_run() and os.path.isfile(bam):
            comm("skipping aligning %s reads; %s already exists" % (sample, bam))
            continue
        sam = os.path.join(work, sample + ".sam")
        comm("Aligning %s reads back to assembly" % sample)
        if args.single_end:
            run(
                BWA_MEM_SINGLE.format(threads=args.threads, assembly=assembly_copy, reads=entry[1]),
                env=env,
                tool="bwa mem",
                log_path=sam,
            )
        elif args.interleaved:
            run(
                BWA_MEM_INTERLEAVED.format(
                    threads=args.threads, assembly=assembly_copy, reads=entry[1]
                ),
                env=env,
                tool="bwa mem",
                log_path=sam,
            )
        else:
            _s, r1, r2 = entry
            run(
                BWA_MEM_PAIRED.format(threads=args.threads, assembly=assembly_copy, r1=r1, r2=r2),
                env=env,
                tool="bwa mem",
                log_path=sam,
            )
        comm("Sorting the %s alignment file" % sample)
        run(
            SAMTOOLS_SORT.format(tmp=tmp, threads=args.threads, bam=bam, sam=sam),
            env=env,
            tool="samtools sort",
        )
        if not dry_run():
            os.remove(sam)


def _bams(work: str) -> str:
    bams = sorted(os.path.join(work, f) for f in listdir(work) if f.endswith(".bam"))
    # Under --dry-run no alignment ran, so show a readable placeholder rather than an empty
    # argument list that makes the printed command look malformed.
    if not bams and dry_run():
        return os.path.join(work, "<sample>.bam")
    return " ".join(bams)


def _run_metabat2(args, env, work, assembly_copy):
    announcement("RUNNING METABAT2")
    depth = os.path.join(work, "metabat_depth.txt")
    metabat_len = max(args.min_len, METABAT_MIN_LEN)
    comm("making contig depth file...")
    run(
        JGI_DEPTH.format(depth=depth, bams=_bams(work)),
        env=env,
        tool="jgi_summarize_bam_contig_depths",
    )
    comm("Starting binning with metaBAT2...")
    run(
        METABAT2.format(
            assembly=assembly_copy,
            depth=depth,
            out=args.output,
            metabat_len=metabat_len,
            threads=args.threads,
        ),
        env=env,
        tool="metabat2",
    )
    # metaBAT2 writes several things into the output prefix that are NOT genomes: the
    # leftovers from --unbinned, the contigs it rejected as too short or too low-depth, and
    # two plain-text summaries. Left in the bin directory, CheckM and bin_refinement score
    # them as if they were bins. Move them aside (keeping them - they are useful) so only
    # real bins remain.
    bins_dir = os.path.join(args.output, "metabat2_bins")
    if dry_run():
        return
    # metaBAT2 appends "total_depth=.. sample_depths=.." to every contig header. Strip it so
    # these bins carry the same contig ids as the assembly and as every other binner's bins -
    # otherwise nothing downstream can match them up.
    #
    # Note the extension test is FASTA_EXTENSIONS, not BIN_EXTENSION: metaBAT2 writes ".fa" and
    # is renumbered to the canonical naming further down, so at this point the files still carry
    # the tool's own names.
    for name in sorted(listdir(bins_dir)):
        if name.endswith(FASTA_EXTENSIONS) and name not in METABAT2_NON_BINS:
            path = os.path.join(bins_dir, name)
            records = list(iter_fasta(path))
            if any(" " in h or "\t" in h for h, _ in records):
                with open(path, "w") as out:
                    out.writelines(
                        ">%s\n%s\n" % (contig_id(header), seq) for header, seq in records
                    )

    extras = os.path.join(args.output, "metabat2_extra")
    for name in sorted(listdir(bins_dir)):
        is_non_bin = name.endswith((".txt", ".tsv")) or name in METABAT2_NON_BINS
        if not is_non_bin:
            continue
        ensure_dir(extras)
        os.replace(os.path.join(bins_dir, name), os.path.join(extras, name))
        comm("moved metaBAT2's %s out of the bin set (it is not a genome) into %s" % (name, extras))

    # Only now that the non-bins are out of the way: give the real bins the same names every
    # other binner's get. Done last so a rejected-contigs file can never consume a bin number.
    renumber_bins(bins_dir)


def renumber_bins(bins_dir: str) -> int:
    """Rename every FASTA in *bins_dir* to bin_001.fasta, bin_002.fasta, ... Returns the count.

    Bins are ordered by any number already in their filename, numerically - so metaBAT2's
    ``bin.2.fa`` stays ahead of its ``bin.10.fa`` instead of being reordered by the lexicographic
    accident that makes "10" sort before "2". Files without a number sort last, by name.

    Renaming goes via a temporary name because the source and target namespaces overlap: without
    it, renaming ``bin.1.fa`` to ``bin_001.fasta`` is fine, but a binner that already wrote
    ``bin_002.fasta`` could be overwritten by another bin being renamed onto it.
    """

    def order(name: str) -> Tuple[int, float, str]:
        digits = re.findall(r"\d+", os.path.splitext(name)[0])
        return (0, int(digits[-1]), name) if digits else (1, float("inf"), name)

    names = sorted(
        (n for n in listdir(bins_dir) if n.endswith(FASTA_EXTENSIONS)),
        key=order,
    )
    staged = []
    for i, name in enumerate(names, start=1):
        temporary = os.path.join(bins_dir, ".renumber.%d" % i)
        os.replace(os.path.join(bins_dir, name), temporary)
        staged.append((temporary, bin_filename(i)))
    for temporary, final in staged:
        os.replace(temporary, os.path.join(bins_dir, final))
    return len(staged)


def _run_maxbin2(args, env, work, assembly_copy):
    announcement("RUNNING MAXBIN2")
    depth = os.path.join(work, "mb2_master_depth.txt")
    comm("making contig depth file...")
    run(
        JGI_DEPTH_MAXBIN.format(depth=depth, bams=_bams(work)),
        env=env,
        tool="jgi_summarize_bam_contig_depths",
    )

    # split the master depth table into one abundance file per sample (cols 4..N)
    comm("splitting master depth file into per-sample abundance files")
    abund_list = os.path.join(work, "mb2_abund_list.txt")
    if dry_run():
        markers = MAXBIN_MARKERS_UNIVERSAL if args.universal else MAXBIN_MARKERS_BACTERIAL
        maxbin_out = os.path.join(work, "maxbin2_out")
        run(
            MAXBIN2.format(
                assembly=assembly_copy,
                markers=markers,
                threads=args.threads,
                min_len=args.min_len,
                maxbin_out=maxbin_out,
                abund_list=abund_list,
            ),
            env=env,
            tool="run_MaxBin.pl",
        )
        return
    with open(depth) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    with open(abund_list, "w") as listing:
        for i in range(3, len(header)):  # columns are 1-indexed 4..N -> 0-indexed 3..
            sample = header[i]
            per = os.path.join(work, "mb2_%s.txt" % os.path.splitext(sample)[0])
            with open(depth) as fh, open(per, "w") as out:
                for line in fh:
                    if "totalAvgDepth" in line:
                        continue
                    cols = line.rstrip("\n").split("\t")
                    out.write("%s\t%s\n" % (cols[0], cols[i]))
            listing.write(os.path.abspath(per) + "\n")

    markers = MAXBIN_MARKERS_UNIVERSAL if args.universal else MAXBIN_MARKERS_BACTERIAL
    maxbin_out = os.path.join(work, "maxbin2_out")
    ensure_dir(maxbin_out)
    comm("Starting binning with MaxBin2...")
    run(
        MAXBIN2.format(
            assembly=assembly_copy,
            markers=markers,
            threads=args.threads,
            min_len=args.min_len,
            maxbin_out=maxbin_out,
            abund_list=abund_list,
        ),
        env=env,
        tool="run_MaxBin.pl",
    )

    out_bins = os.path.join(args.output, "maxbin2_bins")
    ensure_dir(out_bins)
    if dry_run():
        return
    n = 0
    for f in sorted(os.listdir(maxbin_out)):
        if f.startswith("bin") and f.endswith(".fasta"):
            n += 1
            os.replace(os.path.join(maxbin_out, f), os.path.join(out_bins, bin_filename(n)))
    if n == 0:
        error("MaxBin2 did not produce a single bin. Something went wrong. Exiting.")
    comm("MaxBin2 finished successfully, and found %d bins!" % n)


def _run_concoct(args, env, work, assembly_copy):
    announcement("RUNNING CONCOCT")
    depth = os.path.join(work, "concoct_depth.txt")
    bed = os.path.join(work, "assembly_10K.bed")
    comp = os.path.join(work, "assembly_10K.fa")
    if dry_run() or not (os.path.isfile(depth) and os.path.getsize(depth)):
        comm("indexing .bam alignment files...")
        for f in sorted(listdir(work)):
            if f.endswith(".bam"):
                run(
                    SAMTOOLS_INDEX.format(threads=args.threads, bam=os.path.join(work, f)),
                    env=env,
                    tool="samtools index",
                )
        comm("cutting up contigs into 10kb fragments for CONCOCT...")
        run(
            CONCOCT_CUTUP.format(assembly=assembly_copy, bed=bed),
            env=env,
            tool="cut_up_fasta.py",
            log_path=comp,
        )
        comm("estimating contig fragment coverage...")
        run(
            CONCOCT_COVERAGE.format(bed=bed, bams=_bams(work)),
            env=env,
            tool="concoct_coverage_table.py",
            log_path=depth,
        )

    concoct_out = os.path.join(work, "concoct_out")
    ensure_dir(concoct_out)
    comm("Starting binning with CONCOCT...")
    run(
        CONCOCT.format(
            min_len=args.min_len,
            threads=args.threads,
            cov=depth,
            comp=comp,
            concoct_out=os.path.join(concoct_out, ""),
        ).rstrip(),
        env=env,
        tool="concoct",
    )
    clustering = os.path.join(concoct_out, "clustering_gt%d.csv" % args.min_len)
    merged = os.path.join(concoct_out, "clustering_gt%d_merged.csv" % args.min_len)
    comm("merging 10kb fragments back into contigs")
    run(
        CONCOCT_MERGE.format(clustering=clustering),
        env=env,
        tool="merge_cutup_clustering.py",
        log_path=merged,
    )
    comm("splitting contigs into bins")
    if not dry_run():
        _split_concoct_bins(merged, assembly_copy, os.path.join(args.output, "concoct_bins"))
    comm("CONCOCT finished successfully!")


def _count_bins(path: str) -> int:
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
    )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not (args.metabat2 or args.maxbin2 or args.concoct):
        error("You must select at least one binning method: --metabat2, --maxbin2, --concoct")

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    check_disk_space("binning", settings, [args.assembly] + args.reads, args.output)
    validate_inputs(
        "binning",
        [(check_fasta, args.assembly, "assembly (-a)")]
        + [(check_fastq, r, "read file") for r in args.reads],
    )
    env = env_for("binning", settings)
    rec = start_run("binning", args, env, settings, inputs=[args.assembly] + list(args.reads))

    ckpt = make_checkpoint(args.output)
    try:
        announcement("ALIGNING READS TO MAKE COVERAGE FILES")
        ensure_dir(args.output)
        work = os.path.join(args.output, "work_files")
        ensure_dir(work)
        assembly_copy = os.path.join(work, "assembly.fa")
        comm("making copy of assembly file %s" % args.assembly)
        if not dry_run():
            with open(assembly_copy, "w") as out:
                out.writelines(
                    ">%s\n%s\n" % (header, seq) for header, seq in iter_fasta(args.assembly)
                )

        try:
            if ckpt.todo("align_reads"):
                _align_reads(args, env, work, assembly_copy)
                ckpt.done("align_reads")
            else:
                comm("skipping read alignment (already done; --resume)")
            if args.metabat2 and ckpt.todo("metabat2"):
                _run_metabat2(args, env, work, assembly_copy)
                ckpt.done("metabat2")
            if args.maxbin2 and ckpt.todo("maxbin2"):
                _run_maxbin2(args, env, work, assembly_copy)
                ckpt.done("maxbin2")
            if args.concoct and ckpt.todo("concoct"):
                _run_concoct(args, env, work, assembly_copy)
                ckpt.done("concoct")
        except ValueError as e:  # read-collection / validation errors
            error(str(e))

        if args.run_checkm:
            refinement_env = env_for("bin_refinement", settings)
            for binner_dir in ("metabat2_bins", "maxbin2_bins", "concoct_bins"):
                path = os.path.join(args.output, binner_dir)
                if _count_bins(path):
                    _checkm.run_checkm(
                        path,
                        threads=args.threads,
                        mem_gb=args.memory,
                        env=refinement_env,
                        quick=False,
                    )

        announcement("BINNING PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
