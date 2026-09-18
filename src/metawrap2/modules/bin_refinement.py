"""MetaWrap2 bin_refinement module: consolidate several bin sets into one refined set.

Given 1-3 folders of bins produced from the same assembly (e.g. metaBAT2, MaxBin2, CONCOCT
outputs), this hybridizes them with Binning_refiner, scores every candidate set with CheckM,
consolidates the best version of each bin, and dereplicates shared contigs. The scientific
core lives in :mod:`metawrap2.refinement` (consolidate/dereplicate) and
:mod:`metawrap2.checkm` (CheckM); this module orchestrates them and preserves the original
output layout.

This file is meant to be read and edited: the commands each external tool runs live in the
COMMANDS block below and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import os
import shutil
from typing import List, Optional

from .. import checkm as _checkm
from .. import refinement as _refinement
from ..config import load_settings
from ..constants import MAX_BIN_SIZE, MIN_BIN_SIZE
from ..io.seqio import iter_fasta
from ..utils import strip_trailing_slash
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run,
    start_run, warning,
)

CONDA_ENV = "metawrap2-bin_refinement"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BINNING_REFINER    = "python -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -o {out}"
BINNING_REFINER_3  = "python -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -3 {b3} -o {out}"
PLOT_BINNING       = "python -m metawrap2.scripts.plot_binning_results {comp} {cont} {stats}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MEM_GB          = 40
DEFAULT_MIN_COMPLETION  = 70   # -c; bins must exceed this % completion to be "good"
DEFAULT_MAX_CONTAM      = 10   # -x; bins must stay below this % contamination
# bin_refinement ignores bins outside this size range (bytes); see constants.py
MIN_SIZE = MIN_BIN_SIZE        # 50 kb
MAX_SIZE = MAX_BIN_SIZE        # 20 Mb
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 bin_refinement",
        usage="metawrap2 bin_refinement [options] -o output_dir -A bin_folderA "
              "[-B bin_folderB -C bin_folderC]",
        description="Consolidate 1-3 bin sets (from the same assembly) into a refined set. "
                    "Contig names must be consistent across the bin folders.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-t", "--threads", type=int, default=1, help="number of threads (default 1)")
    p.add_argument("-m", "--memory", type=int, default=DEFAULT_MEM_GB, help="RAM in GB (default 40)")
    p.add_argument("-c", "--completeness", type=float, default=DEFAULT_MIN_COMPLETION,
                   help="minimum %% completion of bins [should be >50] (default 70)")
    p.add_argument("-x", "--contamination", type=float, default=DEFAULT_MAX_CONTAM,
                   help="maximum acceptable %% contamination of bins (default 10)")
    p.add_argument("-A", "--bins-a", required=True, help="folder with metagenomic bins (.fa/.fasta)")
    p.add_argument("-B", "--bins-b", help="another folder with metagenomic bins")
    p.add_argument("-C", "--bins-c", help="another folder with metagenomic bins")
    p.add_argument("--skip-refinement", dest="refine", action="store_false",
                   help="don't use Binning_refiner to make refined bins from binner combinations")
    p.add_argument("--skip-checkm", dest="run_checkm", action="store_false",
                   help="don't run CheckM to assess bins")
    p.add_argument("--skip-consolidation", dest="consolidate", action="store_false",
                   help="don't choose the best version of each bin; just pick the best folder")
    amb = p.add_mutually_exclusive_group()
    amb.add_argument("--keep-ambiguous", dest="dereplicate", action="store_const", const="keep",
                     help="for contigs in more than one bin, keep them in all bins")
    amb.add_argument("--remove-ambiguous", dest="dereplicate", action="store_const", const="remove",
                     help="for contigs in more than one bin, remove them from all bins")
    p.set_defaults(dereplicate="partial")  # default: keep only in the best bin
    p.add_argument("--quick", action="store_true",
                   help="add --reduced_tree to CheckM, reducing runtime (helps with low memory)")
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _fix_contig_naming(path: str) -> None:
    """Rewrite a bin fasta in place, replacing '=' with '_' in headers (ports fix_config_naming.py)."""
    records = [(h.replace("=", "_"), s) for h, s in iter_fasta(path)]
    with open(path, "w") as out:
        for header, seq in records:
            out.write(">%s\n%s\n" % (header, seq))


def _copy_bins(src: str, dest: str) -> int:
    """Copy size-filtered bins from *src* into *dest* as bin-name.fa. Returns count kept."""
    ensure_dir(dest)
    n = 0
    for name in sorted(os.listdir(src)):
        full = os.path.join(src, name)
        if not os.path.isfile(full):
            continue
        size = os.path.getsize(full)  # portable (no `stat -c`)
        if MIN_SIZE < size < MAX_SIZE:
            base = os.path.splitext(name)[0]
            shutil.copy(full, os.path.join(dest, base + ".fa"))
            n += 1
        else:
            comm("Skipping %s because the bin size is not between 50kb and 20Mb" % full)
    return n


def _count_bins(path: str) -> int:
    return len([f for f in os.listdir(path) if f.endswith(".fa")]) if os.path.isdir(path) else 0


def _refine(out: str, env: Optional[str], b1: str, b2: str, dest: str,
            b3: Optional[str] = None) -> None:
    """Run Binning_refiner on 2 or 3 bin folders (relative to *out*) into *dest*."""
    refined = "Refined_" + dest[len("bins"):]  # binsAB -> Refined_AB
    if b3 is None:
        cmd = BINNING_REFINER.format(b1=b1, b2=b2, out=refined)
    else:
        cmd = BINNING_REFINER_3.format(b1=b1, b2=b2, b3=b3, out=refined)
    run(cmd, env=env, tool="binning_refiner", cwd=out)
    produced = os.path.join(out, refined, "Refined")
    if not os.path.isdir(produced):
        error("Bin_refiner did not finish correctly for %s. Exiting..." % dest)
    os.replace(produced, os.path.join(out, dest))
    comm("there are %d refined bins in %s" % (_count_bins(os.path.join(out, dest)), dest))
    shutil.rmtree(os.path.join(out, refined))


def _rename_fasta_to_fa(bin_set_dir: str) -> None:
    """Rename any *.fasta in a bin-set directory to *.fa for consistency."""
    for f in os.listdir(bin_set_dir):
        if f.endswith(".fasta"):
            os.replace(os.path.join(bin_set_dir, f),
                       os.path.join(bin_set_dir, os.path.splitext(f)[0] + ".fa"))


def _good_bins(stats_path: str, comp: float, cont: float) -> List[str]:
    """Return names of bins in a .stats file meeting the completion/contamination thresholds."""
    good = []
    if not os.path.isfile(stats_path):
        return good
    with open(stats_path) as fh:
        for line in fh:
            if "completeness" in line or "compl" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            if comp <= float(cut[1]) <= 100 and 0 <= float(cut[2]) <= cont:
                good.append(cut[0])
    return good


def _run_checkm_set(bin_set_dir: str, args, env: Optional[str]) -> str:
    """Run CheckM on one bin set and report the number of 'good' bins. Returns the .stats path."""
    comm("Running CheckM on %s bins" % os.path.basename(bin_set_dir))
    stats = _checkm.run_checkm(bin_set_dir, threads=args.threads, mem_gb=args.memory,
                               env=env, quick=args.quick)
    n = len(_good_bins(stats, args.completeness, args.contamination))
    comm("There are %d 'good' bins found in %s! (>%s%% completion and <%s%% contamination)"
         % (n, os.path.basename(bin_set_dir), args.completeness, args.contamination))
    return stats


def _bin_set_dirs(out: str) -> List[str]:
    """Bin-set directories in *out* (binsA, binsB, ..., binsAB, ...), sorted, excluding binsM/binsO."""
    dirs = []
    for name in sorted(os.listdir(out)):
        if name.startswith("bins") and name not in ("binsM", "binsO") \
                and os.path.isdir(os.path.join(out, name)):
            dirs.append(name)
    return dirs


def main(argv: List[str]) -> int:  # noqa: C901 - orchestration mirrors the original pipeline
    args = _parse_args(argv)
    settings = load_settings(args.config)
    env = env_for("bin_refinement", settings)
    rec = start_run("bin_refinement", args, env, settings, inputs=[args.bins_a, args.bins_b, args.bins_c])

    ckpt = make_checkpoint(args.output)
    try:

        comp, cont = args.completeness, args.contamination
        folders = [(f, s) for f, s in (("binsA", args.bins_a), ("binsB", args.bins_b),
                                       ("binsC", args.bins_c)) if s]
        if not os.path.isdir(args.bins_a):
            error("%s is not a valid directory. Exiting." % args.bins_a)

        announcement("BEGIN PIPELINE!")
        out = args.output
        comm("setting up output folder and copying over bins...")
        if os.path.isdir(out):
            warning("Warning: %s already exists. Attempting to clean." % out)
            for name in ("binsA", "binsB", "binsC", "binsAB", "binsBC", "binsAC", "binsABC",
                         "binsM", "binsO"):
                shutil.rmtree(os.path.join(out, name), ignore_errors=True)
        ensure_dir(out)

        n_binnings = 0
        for dest, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)
            kept = _copy_bins(src, os.path.join(out, dest))
            comm("there are %d bins in %s" % (kept, dest))
            if kept == 0:
                error("Please provide valid input. Exiting...")
            n_binnings += 1
        comm("There are %d bin sets!" % n_binnings)

        comm("Fix contig naming by removing special characters...")
        for dest, _src in folders:
            for f in os.listdir(os.path.join(out, dest)):
                _fix_contig_naming(os.path.join(out, dest, f))

        if args.refine and n_binnings >= 2:
            if ckpt.todo("refine"):
                announcement("BEGIN BIN REFINEMENT")
                if n_binnings == 2:
                    comm("There are two bin folders, so we can consolidate them into a third, more "
                         "refined bin set.")
                    _refine(out, env, "binsA", "binsB", "binsAB")
                elif n_binnings == 3:
                    comm("There are three bin folders, so there are 4 ways to refine the bins "
                         "(A+B, B+C, A+C, A+B+C). Trying all four!")
                    _refine(out, env, "binsA", "binsB", "binsAB")
                    _refine(out, env, "binsC", "binsB", "binsBC")
                    _refine(out, env, "binsA", "binsC", "binsAC")
                    _refine(out, env, "binsA", "binsB", "binsABC", b3="binsC")
                comm("Bin refinement finished successfully!")
                ckpt.done("refine")
            else:
                comm("skipping bin refinement (already done; --resume)")
        elif not args.refine:
            comm("Skipping bin refinement. Will proceed with the %d bins specified." % n_binnings)
        else:
            comm("There is only one bin folder, so no refinement of bins possible. Moving on...")

        comm("fixing bin naming to .fa convention for consistency...")
        for name in _bin_set_dirs(out):
            _rename_fasta_to_fa(os.path.join(out, name))

        comm("making sure every refined bin set contains bins...")
        for name in _bin_set_dirs(out):
            if _count_bins(os.path.join(out, name)) == 0:
                comm("Removing bin set %s because it yielded 0 refined bins ..." % name)
                shutil.rmtree(os.path.join(out, name))

        # ── CheckM on all bin sets ──────────────────────────────────────────────────────────
        if args.run_checkm:
            if ckpt.todo("checkm"):
                announcement("RUNNING CHECKM ON ALL SETS OF BINS")
                for name in _bin_set_dirs(out):
                    _run_checkm_set(os.path.join(out, name), args, env)
                ckpt.done("checkm")
            else:
                comm("skipping CheckM on bin sets (already done; --resume)")
        else:
            comm("Skipping CheckM. Warning: bin consolidation will not be possible.")

        # ── consolidate ─────────────────────────────────────────────────────────────────────
        best_bin_set = None
        if args.consolidate and args.run_checkm:
            announcement("CONSOLIDATING ALL BIN SETS BY CHOOSING THE BEST VERSION OF EACH BIN")
            if n_binnings == 1:
                comm("There is only one original bin folder, so no consolidation possible. Moving on...")
                best_bin_set = "binsA"
            else:
                if ckpt.todo("consolidate"):
                    comm("There are %d original bin folders, plus the refined bins." % n_binnings)
                    binsM = os.path.join(out, "binsM")
                    shutil.copytree(os.path.join(out, "binsA"), binsM)
                    shutil.copy(os.path.join(out, "binsA.stats"), binsM + ".stats")
                    for name in _bin_set_dirs(out):
                        if name == "binsA":
                            continue
                        stats = os.path.join(out, name + ".stats")
                        if not os.path.isfile(stats):
                            continue
                        comm("merging %s and binsM" % name)
                        tmp = os.path.join(out, "binsM1")
                        _refinement.consolidate(binsM, os.path.join(out, name), binsM + ".stats",
                                                stats, tmp, comp, cont)
                        shutil.rmtree(binsM)
                        os.remove(binsM + ".stats")
                        os.replace(tmp, binsM)
                        os.replace(tmp + ".stats", binsM + ".stats")

                    binsO = os.path.join(out, "binsO")
                    if args.dereplicate == "keep":
                        comm("Skipping dereplication of contigs between bins...")
                        os.replace(binsM, binsO)
                        os.replace(binsM + ".stats", binsO + ".stats")
                    elif args.dereplicate == "partial":
                        comm("Scanning to find duplicate contigs between bins and only keep them in "
                             "the best bin...")
                        _refinement.dereplicate(binsM + ".stats", binsM, binsO, mode="best")
                        shutil.copy(binsM + ".stats", binsO + ".stats")
                    elif args.dereplicate == "remove":
                        comm("Scanning to find duplicate contigs between bins and deleting them in all "
                             "bins...")
                        _refinement.dereplicate(binsM + ".stats", binsM, binsO, mode="remove")
                        shutil.copy(binsM + ".stats", binsO + ".stats")
                    ckpt.done("consolidate")
                else:
                    comm("skipping bin consolidation (already done; --resume)")
                best_bin_set = "binsO"
        elif not args.consolidate:
            comm("Skipping bin consolidation. Will pick the best binning folder without mixing bins.")
            if not args.run_checkm:
                comm("cannot decide on best bin set because CheckM was not run. Assuming binsA.")
                best_bin_set = "binsA"
            else:
                best, max_good = None, -1
                for name in _bin_set_dirs(out):
                    n = len(_good_bins(os.path.join(out, name + ".stats"), comp, cont))
                    comm("There are %d 'good' bins found in %s!" % (n, name))
                    if n > max_good:
                        max_good, best = n, name
                best_bin_set = best
                comm("looks like the best bin set is %s" % best_bin_set)
        else:
            best_bin_set = "binsA"
        comm("You will find the best non-reassembled versions of the bins in %s" % best_bin_set)

        # ── finalize: re-QC the consolidated set and lay out the final outputs ───────────────
        if ckpt.todo("finalize"):
            if args.run_checkm and best_bin_set == "binsO" and args.dereplicate != "keep":
                comm("Re-running CheckM on binsO bins")
                binsO = os.path.join(out, "binsO")
                stats = _checkm.run_checkm(binsO, threads=args.threads, mem_gb=args.memory,
                                           env=env, quick=args.quick)
                comm("Removing bins that are inadequate quality...")
                good = set(_good_bins(stats, comp, cont))
                for f in list(os.listdir(binsO)):
                    if f.endswith(".fa") and os.path.splitext(f)[0] not in good:
                        print("%s will be removed because it fell below the quality threshold after "
                              "de-replication of contigs..." % os.path.splitext(f)[0])
                        os.remove(os.path.join(binsO, f))
                # rewrite stats to keep only the good bins
                with open(stats) as fh:
                    lines = fh.readlines()
                with open(stats, "w") as fh:
                    fh.write(lines[0])
                    for line in lines[1:]:
                        cut = line.rstrip("\n").split("\t")
                        if len(cut) >= 3 and comp <= float(cut[1]) <= 100 and 0 <= float(cut[2]) <= cont:
                            fh.write(line)
                comm("Re-evaluating bin quality after contig de-replication is complete! There are still "
                     "%d high quality bins." % len(good))

            if args.run_checkm:
                comm("making completion and contamination ranking plots of final outputs")
                stats_files = " ".join(sorted(os.path.join(out, f) for f in os.listdir(out)
                                              if f.endswith(".stats")))
                ensure_dir(os.path.join(out, "figures"))
                run(PLOT_BINNING.format(comp=int(comp), cont=int(cont), stats=stats_files),
                    env=env, tool="plot_binning_results", cwd=out)
                png = os.path.join(out, "binning_results.png")
                if os.path.isfile(png):
                    os.replace(png, os.path.join(out, "figures", "binning_results.png"))

            # ── move intermediate files aside and lay out the final outputs ──────────────────
            announcement("MOVING OVER TEMPORARY FILES")
            if n_binnings != 1:
                work = os.path.join(out, "work_files")
                ensure_dir(work)
                for f in list(os.listdir(out)):
                    if any(f.startswith(p) for p in ("binsA", "binsB", "binsC", "binsM", "binsO")):
                        os.replace(os.path.join(out, f), os.path.join(work, f))

            final = os.path.join(out, "metawrap_%d_%d_bins" % (int(comp), int(cont)))
            best_dir = os.path.join(out, "work_files", best_bin_set) if n_binnings != 1 \
                else os.path.join(out, best_bin_set)
            if os.path.isdir(best_dir):
                shutil.copytree(best_dir, final)
                if os.path.isfile(best_dir + ".stats"):
                    shutil.copy(best_dir + ".stats", final + ".stats")

            # restore the original input bin sets under their input folder names
            for dest, src in folders:
                base = os.path.basename(strip_trailing_slash(src))
                src_dir = os.path.join(out, "work_files", dest) if n_binnings != 1 \
                    else os.path.join(out, dest)
                if os.path.isdir(src_dir) and not os.path.exists(os.path.join(out, base)):
                    shutil.copytree(src_dir, os.path.join(out, base))
                    if os.path.isfile(src_dir + ".stats"):
                        shutil.copy(src_dir + ".stats", os.path.join(out, base + ".stats"))

            if args.run_checkm:
                comm("making contig membership files (for Anvio and other applications)")
                for name in os.listdir(out):
                    d = os.path.join(out, name)
                    if name.endswith("_bins") and os.path.isdir(d):
                        with open(d + ".contigs", "w") as out_fh:
                            for f in sorted(os.listdir(d)):
                                if f.endswith(".fa"):
                                    bin_name = os.path.splitext(f)[0]
                                    for header, _seq in iter_fasta(os.path.join(d, f)):
                                        out_fh.write("%s\t%s\n" % (header.split()[0], bin_name))
            ckpt.done("finalize")
        else:
            comm("skipping final output layout (already done; --resume)")

        announcement("BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
