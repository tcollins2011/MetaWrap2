"""MetaWrap2 reassemble_bins module: improve bins by recruiting reads and reassembling.

All reads are aligned to the combined bins; each read pair is routed to the bin it maps to
(at a strict and a permissive mismatch cutoff), and every bin is reassembled from its own
reads with SPAdes. CheckM then picks, per bin, whichever of the original / strict / permissive
version is best.

CheckM runs in the bin_refinement conda env (that is where CheckM lives); everything else
runs in this module's env. This file is meant to be read and edited: the tool commands live
in the COMMANDS block and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
from typing import List, Optional

from .. import checkm as _checkm
from ..config import load_settings
from ..io.seqio import iter_fasta
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run,
    start_run, warning,
)

CONDA_ENV = "metawrap2-reassemble_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BWA_INDEX          = "bwa index {assembly}"
# reads are streamed straight into the read-splitter (avoids writing a huge intermediate SAM)
BWA_MEM_FILTER     = ("bwa mem -t {threads} {assembly} {r1} {r2} | "
                      "python -m metawrap2.scripts.filter_reads_for_bin_reassembly "
                      "{bins} {reads_out} {strict} {permissive}")
MINIMAP2_FILTER    = ("minimap2 -t {threads} -ax map-ont {assembly} {nanopore} | "
                      "python -m metawrap2.scripts.filter_nanopore_reads_for_bin_reassembly "
                      "{bins} {reads_out}")
SPADES             = ("spades.py -t {threads} -m {mem} --tmp {tmp} --careful "
                      "--untrusted-contigs {contigs} -1 {r1} -2 {r2} -o {out}")
SPADES_NANOPORE    = ("spades.py -t {threads} -m {mem} --tmp {tmp} --careful "
                      "--untrusted-contigs {contigs} -1 {r1} -2 {r2} --nanopore {nano} -o {out}")
RM_SHORT_CONTIGS   = "python -m metawrap2.scripts.rm_short_contigs {min_len} {scaffolds}"
CHOOSE_BEST_BIN    = "python -m metawrap2.scripts.choose_best_bin {stats} {comp} {cont}"
CHECKM_QA_PLOT     = "checkm bin_qa_plot -x fa {checkm} {bins} {plot}"
PLOT_REASSEMBLY    = "python -m metawrap2.scripts.plot_reassembly {out} {comp} {cont} {stats} {orig_stats}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MEM_GB     = 40
DEFAULT_MIN_LEN    = 500   # -l; contigs shorter than this are dropped from reassemblies
DEFAULT_STRICT     = 2     # --strict-cut-off; max SNPs for "strict" read recruitment
DEFAULT_PERMISSIVE = 5     # --permissive-cut-off; max SNPs for "permissive" read recruitment
OPEN_FILE_LIMIT    = 10000 # raise NOFILE so many per-bin fastq handles can stay open at once
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 reassemble_bins",
        usage="metawrap2 reassemble_bins [options] -o output_dir -b bin_folder "
              "-1 reads_1.fastq -2 reads_2.fastq",
        description="Improve bins by recruiting reads to each bin and reassembling with SPAdes.",
    )
    p.add_argument("-b", "--bins", required=True, help="folder with metagenomic bins")
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-1", dest="reads_1", required=True, help="forward reads for reassembly")
    p.add_argument("-2", dest="reads_2", required=True, help="reverse reads for reassembly")
    p.add_argument("--nanopore", help="nanopore reads to use for reassembly")
    p.add_argument("-t", "--threads", type=int, default=1, help="number of threads (default 1)")
    p.add_argument("-m", "--memory", type=int, default=DEFAULT_MEM_GB, help="RAM in GB (default 40)")
    p.add_argument("-c", "--completeness", type=float, default=70,
                   help="minimum desired bin completion %% (default 70)")
    p.add_argument("-x", "--contamination", type=float, default=10,
                   help="maximum desired bin contamination %% (default 10)")
    p.add_argument("-l", "--min-len", type=int, default=DEFAULT_MIN_LEN,
                   help="minimum contig length to include in reassembly (default 500)")
    p.add_argument("--strict-cut-off", dest="strict", type=int, default=DEFAULT_STRICT,
                   help="maximum allowed SNPs for strict read mapping (default 2)")
    p.add_argument("--permissive-cut-off", dest="permissive", type=int, default=DEFAULT_PERMISSIVE,
                   help="maximum allowed SNPs for permissive read mapping (default 5)")
    p.add_argument("--skip-checkm", dest="run_checkm", action="store_false",
                   help="don't run CheckM to assess bins")
    p.add_argument("--parallel", action="store_true",
                   help="run SPAdes reassemblies in parallel, using 1 thread per bin")
    p.add_argument("--mdmcleaner", action="store_true",
                   help="the bin directory has results from MDMcleaner")
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _raise_open_file_limit() -> None:
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (min(OPEN_FILE_LIMIT, hard), hard))
    except (ImportError, ValueError, OSError):
        warning("Could not raise the open-file limit. If you are reassembling many bins you may "
                "hit a 'too many open files' error; try reassembling fewer bins at a time.")


def _prepare_original_bins(args, original_bins: str) -> None:
    """Populate *original_bins* with clean .fa bin files (handles the MDMcleaner layout)."""
    shutil.rmtree(original_bins, ignore_errors=True)
    if args.mdmcleaner:
        comm("Cleaning mdmcleaner contigs...")
        ensure_dir(original_bins)
        import glob
        for gz in glob.glob(os.path.join(args.bins, "*", "*_kept_contigs.fasta.gz")):
            base = os.path.basename(gz)[: -len(".gz")]
            for suffix in ("_filtered_kept_contigs.fasta", "_kept_contigs.fasta", ".fasta"):
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            dest = os.path.join(original_bins, base + ".fa")
            with open(dest, "w") as out:  # seqio reads the .gz transparently
                for header, seq in iter_fasta(gz):
                    out.write(">%s\n%s\n" % (header.split()[0], seq))
    else:
        shutil.copytree(args.bins, original_bins)


def _combine_bins(original_bins: str, assembly: str) -> None:
    ensure_dir(os.path.dirname(assembly))
    with open(assembly, "w") as out:
        for name in sorted(os.listdir(original_bins)):
            full = os.path.join(original_bins, name)
            if os.path.isfile(full):
                with open(full) as fh:
                    shutil.copyfileobj(fh, out)


def _assemble(fname: str, sp_threads: int, args, out: str, env: Optional[str]) -> None:
    """Reassemble one recruited-read set (``<bin>.strict_1.fastq`` etc.) with SPAdes."""
    stem = fname[: -len("_1.fastq")]           # e.g. bin.1.strict
    bin_name = stem
    orig_contigs = os.path.join(out, "original_bins", stem.rsplit(".", 1)[0] + ".fa")
    reads_dir = os.path.join(out, "reads_for_reassembly")
    r1 = os.path.join(reads_dir, stem + "_1.fastq")
    r2 = os.path.join(reads_dir, stem + "_2.fastq")
    asm_dir = os.path.join(out, "reassemblies", bin_name)
    scaffolds = os.path.join(asm_dir, "scaffolds.fasta")

    if os.path.isfile(scaffolds) and os.path.getsize(scaffolds):
        comm("Looks like %s was already re-assembled. Skipping..." % bin_name)
        return

    tmp_dir = asm_dir + ".tmp"
    ensure_dir(tmp_dir)
    comm("NOW REASSEMBLING %s" % bin_name)
    if args.nanopore:
        nano = os.path.join(reads_dir, stem.rsplit(".", 1)[0] + ".nanopore.fastq")
        cmd = SPADES_NANOPORE.format(threads=sp_threads, mem=args.memory, tmp=tmp_dir,
                                     contigs=orig_contigs, r1=r1, r2=r2, nano=nano, out=asm_dir)
    else:
        cmd = SPADES.format(threads=sp_threads, mem=args.memory, tmp=tmp_dir,
                            contigs=orig_contigs, r1=r1, r2=r2, out=asm_dir)
    run(cmd, env=env, tool="spades.py", check=False)
    if os.path.isfile(scaffolds) and os.path.getsize(scaffolds):
        comm("%s was reassembled successfully!" % bin_name)
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        warning("Something went wrong with reassembling %s" % bin_name)


def _run_checkm_stage(args, out: str, checkm_env: Optional[str]) -> None:
    announcement("RUN CHECKM ON REASSEMBLED BINS")
    reassembled = os.path.join(out, "reassembled_bins")
    original_bins = os.path.join(out, "original_bins")

    # keep the original bins alongside the reassemblies as the ".orig" candidate
    for base in sorted(os.listdir(original_bins)):
        if base.endswith(".fa"):
            shutil.copy(os.path.join(original_bins, base),
                        os.path.join(reassembled, os.path.splitext(base)[0] + ".orig.fa"))

    comm("Running CheckM on best bins (reassembled and original)")
    stats = _checkm.run_checkm(reassembled, threads=args.threads, mem_gb=args.memory,
                               env=checkm_env, quick=False)

    announcement("FINDING THE BEST VERSION OF EACH BIN")
    best_dir = os.path.join(out, "reassembled_best_bins")
    ensure_dir(best_dir)
    chosen = os.path.join(out, "chosen_bins.txt")
    run(CHOOSE_BEST_BIN.format(stats=stats, comp=int(args.completeness), cont=int(args.contamination)),
        env=checkm_env, tool="choose_best_bin", log_path=chosen)
    with open(chosen) as fh:
        names = [ln.strip() for ln in fh if ln.strip()]
    for name in names:
        comm("Copying best bin: %s" % name)
        shutil.copy(os.path.join(reassembled, name + ".fa"), best_dir)
    os.remove(chosen)

    best = os.listdir(best_dir)
    o = len([f for f in best if "orig" in f])
    s = len([f for f in best if "strict" in f])
    per = len([f for f in best if "permissive" in f])
    announcement("Reassembly results are in! %d bins were improved with 'strict' reassembly, %d "
                 "bins were improved with 'permissive' reassembly, and %d bins were not improved "
                 "by any reassembly, and thus will stay the same." % (s, per, o))
    if not best:
        error("there are no good bins found in %s - something went wrong with choosing the best "
              "bins between the reassemblies." % best_dir)

    # shuffle intermediates into work_files and promote the best bins to reassembled_bins
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    for item in ("reassembled_bins", "reassembled_bins.checkm", "reassembled_bins.stats",
                 "reads_for_reassembly", "binned_assembly", "reassemblies"):
        src = os.path.join(out, item)
        if os.path.exists(src):
            os.replace(src, os.path.join(work, item))
    os.replace(best_dir, reassembled)

    comm("Re-running CheckM on the best reassembled bins.")
    shutil.rmtree(reassembled + ".checkm", ignore_errors=True)
    stats = _checkm.run_checkm(reassembled, threads=args.threads, mem_gb=args.memory,
                               env=checkm_env, quick=False)

    comm("Making CheckM plot of %s bins" % reassembled)
    plot_dir = reassembled + ".plot"
    run(CHECKM_QA_PLOT.format(checkm=reassembled + ".checkm", bins=reassembled, plot=plot_dir),
        env=checkm_env, tool="checkm bin_qa_plot", check=False)
    qa_png = os.path.join(plot_dir, "bin_qa_plot.png")
    if os.path.isfile(qa_png):
        os.replace(qa_png, reassembled + ".png")
        shutil.rmtree(plot_dir, ignore_errors=True)
    else:
        warning("Something went wrong with making the CheckM plot.")

    comm("making reassembly N50, completion, and contamination summary plots.")
    orig_stats = os.path.join(out, "original_bins.stats")
    work_stats = os.path.join(work, "reassembled_bins.stats")
    with open(work_stats) as fh, open(orig_stats, "w") as out_fh:
        lines = fh.readlines()
        if lines:
            out_fh.write(lines[0])
        for line in lines[1:]:
            if "orig" in line:
                out_fh.write(line)
    run(PLOT_REASSEMBLY.format(out=out, comp=int(args.completeness), cont=int(args.contamination),
                               stats=stats, orig_stats=orig_stats),
        env=env_for("reassemble_bins", load_settings(args.config)), tool="plot_reassembly")


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    env = env_for("reassemble_bins", settings)
    checkm_env = env_for("bin_refinement", settings)  # CheckM lives in the bin_refinement env
    rec = start_run("reassemble_bins", args, env, settings, inputs=[args.bins, args.reads_1, args.reads_2, args.nanopore])

    ckpt = make_checkpoint(args.output)
    try:

        for path in (args.bins,):
            if not os.path.isdir(path):
                error("%s is not a valid directory. Exiting." % path)
        for path in (args.reads_1, args.reads_2):
            if not os.path.isfile(path):
                error("%s does not exist. Exiting." % path)

        announcement("BEGIN PIPELINE!")
        out = args.output
        comm("setting up output folder and copying over bins...")
        if os.path.isdir(out):
            warning("Warning: %s already exists!" % out)
        ensure_dir(out)

        original_bins = os.path.join(out, "original_bins")
        _prepare_original_bins(args, original_bins)

        assembly = os.path.join(out, "binned_assembly", "assembly.fa")
        _combine_bins(original_bins, assembly)

        reads_out = os.path.join(out, "reads_for_reassembly")
        if ckpt.todo("recruit_reads"):
            announcement("RECRUITING READS TO BINS FOR REASSEMBLY")
            _raise_open_file_limit()
            if not (os.path.isfile(assembly + ".amb")):
                comm("Indexing the assembly")
                run(BWA_INDEX.format(assembly=assembly), env=env, tool="bwa index")
                shutil.rmtree(reads_out, ignore_errors=True)
                ensure_dir(reads_out)
                comm("Aligning all reads back to entire assembly and splitting reads into individual "
                     "fastq files based on their bin membership")
                if args.nanopore:
                    run(["bash", "-c", MINIMAP2_FILTER.format(threads=args.threads, assembly=assembly,
                                                              nanopore=args.nanopore, bins=original_bins,
                                                              reads_out=reads_out)],
                        env=env, tool="minimap2 | filter_nanopore_reads_for_bin_reassembly")
                run(["bash", "-c", BWA_MEM_FILTER.format(threads=args.threads, assembly=assembly,
                                                         r1=args.reads_1, r2=args.reads_2,
                                                         bins=original_bins, reads_out=reads_out,
                                                         strict=args.strict, permissive=args.permissive)],
                    env=env, tool="bwa mem | filter_reads_for_bin_reassembly")
            else:
                comm("WARNING: Looks like the assembly was already indexed. Skipping indexing and read "
                     "splitting, and proceeding directly to assembly. If this is not intended, re-run "
                     "with a new/empty output folder.")
            ckpt.done("recruit_reads")
        else:
            comm("skipping read recruitment (already done; --resume)")

        reassembled = os.path.join(out, "reassembled_bins")
        if ckpt.todo("reassemble"):
            announcement("REASSEMBLING BINS WITH SPADES")
            ensure_dir(os.path.join(out, "reassemblies"))
            to_assemble = sorted(f for f in os.listdir(reads_out) if f.endswith("_1.fastq"))
            if args.parallel:
                with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
                    list(pool.map(lambda f: _assemble(f, 1, args, out, env), to_assemble))
                comm("all assemblies complete")
            else:
                for f in to_assemble:
                    _assemble(f, args.threads, args, out, env)
                comm("all assemblies complete")

            comm("Finalizing reassemblies")
            ensure_dir(reassembled)
            for name in sorted(os.listdir(os.path.join(out, "reassemblies"))):
                asm_dir = os.path.join(out, "reassemblies", name)
                scaffolds = os.path.join(asm_dir, "scaffolds.fasta")
                if not (os.path.isfile(scaffolds) and os.path.getsize(scaffolds)):
                    comm("%s was not successfully reassembled. It is possible there were not enough "
                         "reads." % name)
                    continue
                long_scaffolds = os.path.join(asm_dir, "long_scaffolds.fasta")
                run(RM_SHORT_CONTIGS.format(min_len=args.min_len, scaffolds=scaffolds),
                    env=env, tool="rm_short_contigs", log_path=long_scaffolds)
                if os.path.isfile(long_scaffolds) and os.path.getsize(long_scaffolds):
                    print("%s was reassembled! Processing..." % name)
                    os.replace(long_scaffolds, os.path.join(reassembled, name + ".fa"))
                else:
                    comm("%s was reassembled, but did not yield contigs %d bp. It is possible there were "
                         "not enough reads." % (name, args.min_len))

            if _count_fa(reassembled) < 1:
                error("None of the bins were successfully reassembled. %s is empty." % reassembled)
            comm("Looks like the reassemblies went well. Now to see if they made the bins better or "
                 "worse...")
            ckpt.done("reassemble")
        else:
            comm("skipping bin reassembly (already done; --resume)")

        if args.run_checkm and ckpt.todo("checkm"):
            _run_checkm_stage(args, out, checkm_env)
            ckpt.done("checkm")
        elif args.run_checkm:
            comm("skipping CheckM on reassembled bins (already done; --resume)")

        comm("you will find the final bins in %s" % reassembled)
        announcement("BIN REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


def _count_fa(path: str) -> int:
    return len([f for f in os.listdir(path) if f.endswith(".fa")]) if os.path.isdir(path) else 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
