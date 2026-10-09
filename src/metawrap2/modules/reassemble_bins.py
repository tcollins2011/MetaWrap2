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
from ..cliformat import bounded_float, bounded_int, make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, STATS_SUFFIX
from ..io.seqio import iter_fasta
from ..progress import bar
from ..pyrun import PY
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
    finish_run,
    force,
    fs_replace_dir,
    list_bins,
    listdir,
    make_checkpoint,
    memory_arg,
    require_nonempty_dir,
    resolve_memory_gb,
    resolve_threads,
    run,
    scratch_dir,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["reassemble_bins"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BWA_INDEX = "bwa index {assembly}"
# reads are streamed straight into the read-splitter (avoids writing a huge intermediate SAM)
# The right-hand side of these pipes is a MetaWrap2 helper, so it runs in the host
# interpreter (PY) - the module env has bwa/minimap2/spades, not python+metawrap2.
#
# pipefail, because bash otherwise reports only the last command's status: the filter reads
# whatever arrives on stdin and exits 0 on a truncated stream, so an aligner that died partway
# would look like a clean run that happened to select fewer reads, and the bins would be
# reassembled from a silently truncated subset.
_PIPEFAIL = "set -o pipefail; "
BWA_MEM_FILTER = (
    _PIPEFAIL
    + "bwa mem -t {threads} {assembly} {r1} {r2} | "
    + PY
    + " -m metawrap2.scripts.filter_reads_for_bin_reassembly "
    "{bins} {reads_out} {strict} {permissive}"
)
MINIMAP2_FILTER = (
    _PIPEFAIL
    + "minimap2 -t {threads} -ax map-ont {assembly} {nanopore} | "
    + PY
    + " -m metawrap2.scripts.filter_nanopore_reads_for_bin_reassembly "
    "{bins} {reads_out}"
)
SPADES = (
    "spades.py -t {threads} -m {mem} --tmp {tmp} --careful "
    "--untrusted-contigs {contigs} -1 {r1} -2 {r2} -o {out}"
)
SPADES_NANOPORE = (
    "spades.py -t {threads} -m {mem} --tmp {tmp} --careful "
    "--untrusted-contigs {contigs} -1 {r1} -2 {r2} --nanopore {nano} -o {out}"
)
RM_SHORT_CONTIGS = PY + " -m metawrap2.scripts.rm_short_contigs {min_len} {scaffolds}"
CHOOSE_BEST_BIN = PY + " -m metawrap2.scripts.choose_best_bin {stats} {comp} {cont}"
CHECKM_QA_PLOT = "checkm bin_qa_plot -x " + BIN_EXTENSION.lstrip(".") + " {checkm} {bins} {plot}"
PLOT_REASSEMBLY = (
    PY + " -m metawrap2.scripts.plot_reassembly {out} {comp} {cont} {stats} {orig_stats}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MIN_LEN = 500  # -l; contigs shorter than this are dropped from reassemblies
DEFAULT_STRICT = 2  # --strict-cut-off; max SNPs for "strict" read recruitment
DEFAULT_PERMISSIVE = 5  # --permissive-cut-off; max SNPs for "permissive" read recruitment
OPEN_FILE_LIMIT = 10000  # raise NOFILE so many per-bin fastq handles can stay open at once
#: Depth of coverage each bin's recruited reads are thinned to before SPAdes sees them, as a multiple
#: of the bin's own assembled size. ``--max-coverage 0`` turns the thinning off.
#:
#: This is the one knob in the module that is about cost rather than correctness. Read recruitment is
#: per bin, so an abundant organism hands SPAdes hundreds of times the coverage it can use, and SPAdes
#: time grows with reads, not with genome size - so a handful of high-abundance outliers decide how
#: long the whole step takes. 50x is comfortably above where short-read assembly of a prokaryotic
#: genome stops improving, which is why thinning to it is close to free in assembly quality and not
#: close to free in time.
#:
#: Expressed against the bin's assembled length rather than an assumed genome size: the bin is the
#: thing being reassembled, and a 1.2 Mbp bin and a 6 Mbp bin want very different numbers of reads
#: for the same depth.
DEFAULT_MAX_COVERAGE = 50.0
#: Seed for the thinning. Fixed, so a re-run recruits and keeps the same reads: a subsample that
#: changed between runs would make reassembly irreproducible for the same reason an unseeded binner is.
SUBSAMPLE_SEED = 1337
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
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
    threads_arg(p)
    memory_arg(p)
    p.add_argument(
        "-c",
        "--completeness",
        type=percent(),
        default=70,
        help="minimum desired bin completion %%",
    )
    p.add_argument(
        "-x",
        "--contamination",
        type=percent(),
        default=10,
        help="maximum desired bin contamination %%",
    )
    p.add_argument(
        "-l",
        "--min-len",
        type=bounded_int(1),
        default=DEFAULT_MIN_LEN,
        help="minimum contig length to include in reassembly",
    )
    p.add_argument(
        "--strict-cut-off",
        dest="strict",
        type=bounded_int(0),
        default=DEFAULT_STRICT,
        help="maximum allowed SNPs for strict read mapping",
    )
    p.add_argument(
        "--permissive-cut-off",
        dest="permissive",
        type=bounded_int(0),
        default=DEFAULT_PERMISSIVE,
        help="maximum allowed SNPs for permissive read mapping",
    )
    p.add_argument(
        "--max-coverage",
        type=bounded_float(0),
        default=DEFAULT_MAX_COVERAGE,
        metavar="X",
        help="thin each bin's recruited reads to this depth of coverage before reassembly, "
        "measured against that bin's own length; 0 keeps every read",
    )
    p.add_argument(
        "--skip-checkm",
        dest="run_checkm",
        action="store_false",
        help="don't run CheckM to assess bins",
    )
    p.add_argument(
        "--parallel",
        action="store_true",
        help="run SPAdes reassemblies in parallel, using 1 thread per bin",
    )
    p.add_argument(
        "--mdmcleaner", action="store_true", help="the bin directory has results from MDMcleaner"
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _bin_lengths(bins_dir: str) -> dict:
    """Assembled length of every bin in *bins_dir*, keyed by the bin's name without its extension."""
    lengths = {}
    for name in listdir(bins_dir):
        if not name.endswith(FASTA_EXTENSIONS):
            continue
        total = sum(len(seq) for _header, seq in iter_fasta(os.path.join(bins_dir, name)))
        if total:
            lengths[os.path.splitext(name)[0]] = total
    return lengths


def _estimate_bases(path: str, records: int = 2000) -> float:
    """Roughly how many sequence bases *path* holds, from the first *records* of it.

    Reads a sample, measures what fraction of those bytes were sequence, and scales by the file size.
    Exact for a file of uniform read length and header width, which recruited reads are, and close
    enough otherwise for a decision whose threshold is a round number.
    """
    size = os.path.getsize(path)
    if not size:
        return 0.0
    consumed = sequence = 0
    with open(path) as fh:
        for index, line in enumerate(fh):
            consumed += len(line)
            if index % 4 == 1:
                sequence += len(line.rstrip("\n"))
            if index >= records * 4 - 1:
                break
    if not consumed:
        return 0.0
    return size * (sequence / consumed)


def _thin_pair(r1: str, r2: str, keep: float, rng) -> int:
    """Keep *keep* of the read pairs in *r1*/*r2* in place. Returns how many pairs were kept.

    Pairs, not reads: SPAdes is being handed a paired library and a file where one mate survived and
    the other did not is worse input than a smaller library. Both files are walked in lockstep and the
    same coin decides both, which is what keeps them paired.

    Written to ``.thinned`` siblings and moved over the originals, so an interruption leaves the
    untouched originals rather than a half-written library that looks complete.
    """
    kept = 0
    out1, out2 = r1 + ".thinned", r2 + ".thinned"
    with open(r1) as f1, open(r2) as f2, open(out1, "w") as o1, open(out2, "w") as o2:
        while True:
            first = [f1.readline() for _ in range(4)]
            second = [f2.readline() for _ in range(4)]
            if not first[0] or not second[0]:
                break
            if rng.random() < keep:
                o1.writelines(first)
                o2.writelines(second)
                kept += 1
    os.replace(out1, r1)
    os.replace(out2, r2)
    return kept


def _cap_recruited_reads(args, reads_out: str, original_bins: str) -> None:
    """Thin each bin's recruited reads down to ``--max-coverage``, leaving thin bins alone.

    Runs after recruitment and before SPAdes, over the per-bin fastqs the filter wrote. A bin whose
    recruited reads already sit below the target is not touched at all, which is the common case -
    the thinning exists for the few abundant organisms that otherwise set the pace of the whole step.
    """
    if not args.max_coverage:
        comm("--max-coverage 0: every recruited read is kept")
        return
    import random

    lengths = _bin_lengths(original_bins)
    if not lengths:
        return
    rng = random.Random(SUBSAMPLE_SEED)
    thinned = 0
    for name in sorted(listdir(reads_out)):
        if not name.endswith("_1.fastq"):
            continue
        r1 = os.path.join(reads_out, name)
        r2 = os.path.join(reads_out, name[: -len("_1.fastq")] + "_2.fastq")
        if not os.path.isfile(r2):
            continue
        # "<bin>.strict_1.fastq" -> "<bin>"; the cut-off is not part of the bin's identity.
        bin_name = name[: -len("_1.fastq")].rsplit(".", 1)[0]
        length = lengths.get(bin_name)
        if not length:
            continue
        # Estimated rather than counted: deciding *whether* to thin does not need an exact base
        # count, and counting would mean reading every recruited read twice. The estimate is measured
        # instead of assumed - the ratio of sequence bytes to total bytes is read off the first
        # records of the file and applied to its size - because the obvious assumption is wrong. A
        # FASTQ record is a header, a sequence, a "+" and a quality string of the same length as the
        # sequence, so sequence is about *half* the bytes, less header overhead, and how much less
        # depends entirely on how long the read names are.
        bases = _estimate_bases(r1) + _estimate_bases(r2)
        coverage = bases / float(length)
        if coverage <= args.max_coverage:
            continue
        keep = args.max_coverage / coverage
        kept = _thin_pair(r1, r2, keep, rng)
        comm(
            "%s: ~%.0fx recruited, thinned to ~%.0fx (%d pairs kept)"
            % (name[: -len("_1.fastq")], coverage, args.max_coverage, kept)
        )
        thinned += 1
    if thinned:
        comm(
            "thinned %d recruited read set(s) to %gx (--max-coverage)"
            % (thinned, args.max_coverage)
        )
    else:
        comm("no recruited read set was above %gx, so none was thinned" % args.max_coverage)


def _raise_open_file_limit() -> None:
    try:
        import resource

        _soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        resource.setrlimit(resource.RLIMIT_NOFILE, (min(OPEN_FILE_LIMIT, hard), hard))
    except (ImportError, ValueError, OSError):
        warning(
            "Could not raise the open-file limit. If you are reassembling many bins you may "
            "hit a 'too many open files' error; try reassembling fewer bins at a time."
        )


def _prepare_original_bins(args, original_bins: str) -> None:
    """Populate *original_bins* with clean .fa bin files (handles the MDMcleaner layout)."""
    if dry_run():
        return
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
            dest = os.path.join(original_bins, base + BIN_EXTENSION)
            with open(dest, "w") as out:  # seqio reads the .gz transparently
                out.writelines(
                    ">%s\n%s\n" % (header.split()[0], seq) for header, seq in iter_fasta(gz)
                )
    else:
        shutil.copytree(args.bins, original_bins)


def _combine_bins(original_bins: str, assembly: str) -> None:
    """Concatenate the bins into the reference the reads get recruited against.

    Only the bins. A bin set may carry the unbinned catch-all, and this is the one place where
    including it would look defensible - unbinned contigs in the reference would hold on to reads that
    would otherwise be recruited to a bin they do not belong in. It still must not be here: the read
    filter assigns each pair to the bin of the contig it mapped to, so those reads would be recruited
    to a bin called "unbinned" and then reassembled as a genome, which is the opposite of leaving them
    out. Precision against misrecruitment is what the mate-agreement filter is for.
    """
    ensure_dir(os.path.dirname(assembly))
    if dry_run():
        return
    with open(assembly, "w") as out:
        for name in list_bins(original_bins):
            with open(os.path.join(original_bins, name)) as fh:
                shutil.copyfileobj(fh, out)


def _assemble(
    fname: str, sp_threads: int, args, out: str, env: Optional[str], settings=None
) -> None:
    """Reassemble one recruited-read set (``<bin>.strict_1.fastq`` etc.) with SPAdes."""
    stem = fname[: -len("_1.fastq")]  # e.g. bin.1.strict
    bin_name = stem
    orig_contigs = os.path.join(out, "original_bins", stem.rsplit(".", 1)[0] + BIN_EXTENSION)
    reads_dir = os.path.join(out, "reads_for_reassembly")
    r1 = os.path.join(reads_dir, stem + "_1.fastq")
    r2 = os.path.join(reads_dir, stem + "_2.fastq")
    asm_dir = os.path.join(out, "reassemblies", bin_name)
    scaffolds = os.path.join(asm_dir, "scaffolds.fasta")

    if not dry_run() and os.path.isfile(scaffolds) and os.path.getsize(scaffolds):
        comm("Looks like %s was already re-assembled. Skipping..." % bin_name)
        return

    tmp_dir = scratch_dir(settings, out, os.path.basename(asm_dir))
    comm("NOW REASSEMBLING %s" % bin_name)
    if args.nanopore:
        nano = os.path.join(reads_dir, stem.rsplit(".", 1)[0] + ".nanopore.fastq")
        cmd = SPADES_NANOPORE.format(
            threads=sp_threads,
            mem=args.memory,
            tmp=tmp_dir,
            contigs=orig_contigs,
            r1=r1,
            r2=r2,
            nano=nano,
            out=asm_dir,
        )
    else:
        cmd = SPADES.format(
            threads=sp_threads,
            mem=args.memory,
            tmp=tmp_dir,
            contigs=orig_contigs,
            r1=r1,
            r2=r2,
            out=asm_dir,
        )
    run(cmd, env=env, tool="spades.py", check=False)
    if dry_run():
        return
    if os.path.isfile(scaffolds) and os.path.getsize(scaffolds):
        comm("%s was reassembled successfully!" % bin_name)
        shutil.rmtree(tmp_dir, ignore_errors=True)
    else:
        # A per-bin SPAdes failure is tolerated (check=False) because one bad bin should not
        # kill the run - but the reason is in that bin's spades.log, and reporting a guess
        # ("not enough reads") sent users looking in the wrong place. Show what SPAdes said.
        warning("Something went wrong with reassembling %s.%s" % (bin_name, _spades_error(asm_dir)))


def _spades_error(asm_dir: str, lines: int = 6) -> str:
    """The ERROR lines from a failed SPAdes run, for inclusion in our own message."""
    log = os.path.join(asm_dir, "spades.log")
    if not os.path.isfile(log):
        return " (no spades.log was written - SPAdes may not have started)"
    try:
        with open(log, errors="replace") as fh:
            interesting = [
                ln.rstrip() for ln in fh if "ERROR" in ln or "== Error ==" in ln or "Cannot " in ln
            ]
    except OSError:
        return " (could not read %s)" % log
    if not interesting:
        return " See %s for details." % log
    return "\nSPAdes reported:\n  %s\nFull log: %s" % ("\n  ".join(interesting[-lines:]), log)


#: What the finished layout keeps out of the way: everything the chosen bins were chosen from. Named
#: here because both the shuffle that creates it and the --force clean that clears it need the list.
INTERMEDIATES = (
    "reassembled_bins",
    "reassembled_bins.checkm",
    "reassembled_bins" + STATS_SUFFIX,
    "reads_for_reassembly",
    "binned_assembly",
    "reassemblies",
)


def _shuffle_intermediates(out: str, best_dir: str, reassembled: str) -> None:
    """Move the candidates into ``work_files`` and promote the chosen bins to ``reassembled_bins``."""
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    for item in INTERMEDIATES:
        src = os.path.join(out, item)
        if not os.path.exists(src):
            continue
        # Destination-safe: on a second reassembly into the same folder the destination is the first
        # run's copy of the same intermediate. Unguarded, this raised at the end of the most expensive
        # module in the suite - after every reassembly was done - and left the rest half-moved.
        fs_replace_dir(src, os.path.join(work, item))
    os.replace(best_dir, reassembled)


def _clean_previous_run(out: str) -> None:
    """Clear a finished reassembly out of *out* so this run is not measured against it.

    The finished layout puts the chosen bins in ``reassembled_bins`` and everything they were chosen
    from in ``work_files``, and the reassembly step adds to ``reassembled_bins`` rather than replacing
    it. So without this, a second run's candidate pool contained the first run's promoted winners -
    including its ``.orig.fa`` copies of bins that may since have changed - and "the best version of
    each bin" was chosen across two runs. Nothing about the output would have said so.

    Callers gate this on ``--force``: ``--resume`` reaches the same already-existing folder and means
    the opposite thing.
    """
    for item in INTERMEDIATES + ("reassembled_best_bins", "work_files"):
        shutil.rmtree(os.path.join(out, item), ignore_errors=True)
    for name in listdir(out):
        if name.startswith("reassembled_bins.") or name == "chosen_bins.txt":
            path = os.path.join(out, name)
            if os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)
            else:
                os.remove(path)


def _run_checkm_stage(args, out: str, checkm_env: Optional[str]) -> None:
    announcement("RUN CHECKM ON REASSEMBLED BINS")
    reassembled = os.path.join(out, "reassembled_bins")
    original_bins = os.path.join(out, "original_bins")

    if dry_run():
        ensure_dir(reassembled)
        stats = _checkm.run_checkm(
            reassembled, threads=args.threads, mem_gb=args.memory, env=checkm_env, quick=False
        )
        run(
            CHOOSE_BEST_BIN.format(
                stats=stats, comp=int(args.completeness), cont=int(args.contamination)
            ),
            env=None,
            tool="choose_best_bin",
            stdout_path=os.path.join(out, "chosen_bins.txt"),
        )
        run(
            CHECKM_QA_PLOT.format(
                checkm=reassembled + ".checkm", bins=reassembled, plot=reassembled + ".plot"
            ),
            env=checkm_env,
            tool="checkm bin_qa_plot",
            check=False,
        )
        run(
            PLOT_REASSEMBLY.format(
                out=out,
                comp=int(args.completeness),
                cont=int(args.contamination),
                stats=stats,
                orig_stats=os.path.join(out, "original_bins" + STATS_SUFFIX),
            ),
            env=None,
            tool="plot_reassembly",
        )
        return

    # keep the original bins alongside the reassemblies as the ".orig" candidate
    for base in list_bins(original_bins, (BIN_EXTENSION,)):
        shutil.copy(
            os.path.join(original_bins, base),
            os.path.join(reassembled, os.path.splitext(base)[0] + ".orig.fa"),
        )

    comm("Running CheckM on best bins (reassembled and original)")
    stats = _checkm.run_checkm(
        reassembled, threads=args.threads, mem_gb=args.memory, env=checkm_env, quick=False
    )

    announcement("FINDING THE BEST VERSION OF EACH BIN")
    best_dir = os.path.join(out, "reassembled_best_bins")
    ensure_dir(best_dir)
    chosen = os.path.join(out, "chosen_bins.txt")
    run(
        CHOOSE_BEST_BIN.format(
            stats=stats, comp=int(args.completeness), cont=int(args.contamination)
        ),
        env=None,
        tool="choose_best_bin",
        log_path=chosen,
    )
    with open(chosen) as fh:
        names = [ln.strip() for ln in fh if ln.strip()]
    for name in names:
        comm("Copying best bin: %s" % name)
        shutil.copy(os.path.join(reassembled, name + BIN_EXTENSION), best_dir)
    os.remove(chosen)

    best = os.listdir(best_dir)
    o = len([f for f in best if "orig" in f])
    s = len([f for f in best if "strict" in f])
    per = len([f for f in best if "permissive" in f])
    announcement(
        "Reassembly results are in! %d bins were improved with 'strict' reassembly, %d "
        "bins were improved with 'permissive' reassembly, and %d bins were not improved "
        "by any reassembly, and thus will stay the same." % (s, per, o)
    )
    if not best:
        error(
            "there are no good bins found in %s - something went wrong with choosing the best "
            "bins between the reassemblies." % best_dir
        )

    _shuffle_intermediates(out, best_dir, reassembled)

    comm("Re-running CheckM on the best reassembled bins.")
    shutil.rmtree(reassembled + ".checkm", ignore_errors=True)
    stats = _checkm.run_checkm(
        reassembled, threads=args.threads, mem_gb=args.memory, env=checkm_env, quick=False
    )

    comm("Making CheckM plot of %s bins" % reassembled)
    plot_dir = reassembled + ".plot"
    run(
        CHECKM_QA_PLOT.format(checkm=reassembled + ".checkm", bins=reassembled, plot=plot_dir),
        env=checkm_env,
        tool="checkm bin_qa_plot",
        check=False,
    )
    qa_png = os.path.join(plot_dir, "bin_qa_plot.png")
    if os.path.isfile(qa_png):
        os.replace(qa_png, reassembled + ".png")
        shutil.rmtree(plot_dir, ignore_errors=True)
    else:
        warning("Something went wrong with making the CheckM plot.")

    comm("making reassembly N50, completion, and contamination summary plots.")
    orig_stats = os.path.join(out, "original_bins" + STATS_SUFFIX)
    work_stats = os.path.join(out, "work_files", "reassembled_bins" + STATS_SUFFIX)
    with open(work_stats) as fh, open(orig_stats, "w") as out_fh:
        lines = fh.readlines()
        if lines:
            out_fh.write(lines[0])
        for line in lines[1:]:
            if "orig" in line:
                out_fh.write(line)
    run(
        PLOT_REASSEMBLY.format(
            out=out,
            comp=int(args.completeness),
            cont=int(args.contamination),
            stats=stats,
            orig_stats=orig_stats,
        ),
        env=None,
        tool="plot_reassembly",
    )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args)
    check_disk_space("reassemble_bins", settings, [args.bins, args.reads_1], args.output)
    validate_inputs(
        "reassemble_bins",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                args.bins,
                "bin folder (-b)",
            ),
            (check_fastq, args.reads_1, "forward reads (-1)"),
            (check_fastq, args.reads_2, "reverse reads (-2)"),
        ],
    )
    env = env_for("reassemble_bins", settings)
    checkm_env = env_for("bin_refinement", settings)  # CheckM lives in the bin_refinement env
    rec = start_run(
        "reassemble_bins",
        args,
        env,
        settings,
        inputs=[args.bins, args.reads_1, args.reads_2, args.nanopore],
    )

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
            # Under --force only; --resume arrives here with the same folder and means the opposite.
            if force():
                comm("--force: clearing the previous reassembly out of %s" % out)
                _clean_previous_run(out)
        ensure_dir(out)

        original_bins = os.path.join(out, "original_bins")
        _prepare_original_bins(args, original_bins)

        assembly = os.path.join(out, "binned_assembly", "assembly.fa")
        _combine_bins(original_bins, assembly)

        reads_out = os.path.join(out, "reads_for_reassembly")
        if ckpt.todo("recruit_reads"):
            announcement("RECRUITING READS TO BINS FOR REASSEMBLY")
            _raise_open_file_limit()
            if dry_run() or not os.path.isfile(assembly + ".amb"):
                comm("Indexing the assembly")
                run(BWA_INDEX.format(assembly=assembly), env=env, tool="bwa index")
                shutil.rmtree(reads_out, ignore_errors=True)
                ensure_dir(reads_out)
                comm(
                    "Aligning all reads back to entire assembly and splitting reads into individual "
                    "fastq files based on their bin membership"
                )
                if args.nanopore:
                    run(
                        [
                            "bash",
                            "-c",
                            MINIMAP2_FILTER.format(
                                threads=args.threads,
                                assembly=assembly,
                                nanopore=args.nanopore,
                                bins=original_bins,
                                reads_out=reads_out,
                            ),
                        ],
                        env=env,
                        tool="minimap2 | filter_nanopore_reads_for_bin_reassembly",
                    )
                run(
                    [
                        "bash",
                        "-c",
                        BWA_MEM_FILTER.format(
                            threads=args.threads,
                            assembly=assembly,
                            r1=args.reads_1,
                            r2=args.reads_2,
                            bins=original_bins,
                            reads_out=reads_out,
                            strict=args.strict,
                            permissive=args.permissive,
                        ),
                    ],
                    env=env,
                    tool="bwa mem | filter_reads_for_bin_reassembly",
                )
            else:
                comm(
                    "WARNING: Looks like the assembly was already indexed. Skipping indexing and read "
                    "splitting, and proceeding directly to assembly. If this is not intended, re-run "
                    "with a new/empty output folder."
                )
            if not dry_run():
                _cap_recruited_reads(args, reads_out, original_bins)
            ckpt.done("recruit_reads")
        else:
            comm("skipping read recruitment (already done; --resume)")

        reassembled = os.path.join(out, "reassembled_bins")
        if ckpt.todo("reassemble"):
            announcement("REASSEMBLING BINS WITH SPADES")
            ensure_dir(os.path.join(out, "reassemblies"))
            to_assemble = sorted(f for f in listdir(reads_out) if f.endswith("_1.fastq"))
            if not to_assemble and dry_run():
                # No reads were recruited (nothing ran), so show the per-bin SPAdes calls for
                # the bins the user actually passed in.
                to_assemble = [
                    "%s.%s_1.fastq" % (os.path.splitext(b)[0], style)
                    for b in list_bins(args.bins, (".fa", ".fasta"))
                    for style in ("strict", "permissive")
                ]
            if args.parallel:
                with concurrent.futures.ThreadPoolExecutor(
                    max_workers=max(1, args.threads)
                ) as pool:
                    list(pool.map(lambda f: _assemble(f, 1, args, out, env, settings), to_assemble))
                comm("all assemblies complete")
            else:
                for f in bar(
                    to_assemble, desc="reassembling bins", unit=" bin", total=len(to_assemble)
                ):
                    _assemble(f, args.threads, args, out, env, settings)
                comm("all assemblies complete")

            comm("Finalizing reassemblies")
            ensure_dir(reassembled)
            # Skip the ".tmp" scratch directories _assemble creates next to each assembly;
            # they are not reassemblies and only survive when SPAdes failed for that bin.
            reassembly_dirs = sorted(
                d for d in listdir(os.path.join(out, "reassemblies")) if not d.endswith(".tmp")
            )
            for name in bar(
                reassembly_dirs,
                desc="collecting reassemblies",
                unit=" bin",
                total=len(reassembly_dirs),
            ):
                asm_dir = os.path.join(out, "reassemblies", name)
                scaffolds = os.path.join(asm_dir, "scaffolds.fasta")
                if not dry_run() and not (os.path.isfile(scaffolds) and os.path.getsize(scaffolds)):
                    comm(
                        "%s was not successfully reassembled (often too few recruited reads)."
                        "%s" % (name, _spades_error(asm_dir))
                    )
                    continue
                long_scaffolds = os.path.join(asm_dir, "long_scaffolds.fasta")
                run(
                    RM_SHORT_CONTIGS.format(min_len=args.min_len, scaffolds=scaffolds),
                    env=None,
                    tool="rm_short_contigs",
                    log_path=long_scaffolds,
                )
                if dry_run():
                    continue
                if os.path.isfile(long_scaffolds) and os.path.getsize(long_scaffolds):
                    print("%s was reassembled! Processing..." % name)
                    os.replace(long_scaffolds, os.path.join(reassembled, name + BIN_EXTENSION))
                else:
                    comm(
                        "%s was reassembled, but did not yield contigs %d bp. It is possible there were "
                        "not enough reads." % (name, args.min_len)
                    )

            if _count_fa(reassembled) < 1 and not dry_run():
                error(
                    "None of the bins were successfully reassembled, so %s is empty.\n"
                    "Every per-bin SPAdes run failed - look at the 'SPAdes reported' lines "
                    "above, or at %s/<bin>/spades.log. A single shared cause is likely; "
                    "common ones are too few recruited reads per bin, and read files SPAdes "
                    "cannot read (for example uniform quality strings, which make it fail "
                    "with 'Failed to determine offset')."
                    % (reassembled, os.path.join(out, "reassemblies"))
                )
            comm(
                "Looks like the reassemblies went well. Now to see if they made the bins better or "
                "worse..."
            )
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
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
    )


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
