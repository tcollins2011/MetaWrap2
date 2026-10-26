"""MetaWrap2 reassemble_bins module: improve bins by recruiting reads and reassembling.

All reads are aligned to the combined bins; each read pair is routed to the bin it maps to
(at a strict and a permissive mismatch cutoff), and every bin is reassembled from its own
reads with SPAdes. CheckM2 then scores every version, and the best of the original / strict /
permissive versions of each bin is kept.

CheckM2 runs in its own env; everything else runs in this module's env. This file is meant to be read and edited: the tool commands live
in the COMMANDS block and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
from typing import Dict, Iterable, List, Optional, Sequence

from .. import checkm as _checkm
from ..cliformat import bounded_float, bounded_int, make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, STATS_SUFFIX
from ..io.seqio import iter_fasta
from ..progress import bar
from ..pyrun import PY
from ._common import (
    absolutize_paths,
    add_checkm2_only_args,
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
    scorer_preflight,
    scratch_dir,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["reassemble_bins"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BWA_INDEX = "bwa index {assembly}"
# Reads are streamed into the read splitter, so no large intermediate SAM is written.
# The right-hand side of these pipes is a MetaWrap2 helper and runs in the host interpreter
# (PY); the module env has bwa/minimap2/spades but not python+metawrap2.
#
# pipefail, because bash otherwise reports only the last command's status. The filter exits 0
# on a truncated stream, so an aligner that died partway would otherwise look like a clean run.
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
PLOT_REASSEMBLY = (
    PY + " -m metawrap2.scripts.plot_reassembly {out} {comp} {cont} {stats} {orig_stats}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MIN_LEN = 500  # -l; contigs shorter than this are dropped from reassemblies
DEFAULT_STRICT = 2  # --strict-cut-off; max SNPs for "strict" read recruitment
DEFAULT_PERMISSIVE = 5  # --permissive-cut-off; max SNPs for "permissive" read recruitment
OPEN_FILE_LIMIT = 10000  # raise NOFILE so many per-bin fastq handles can stay open at once
#: Depth of coverage each bin's recruited reads are thinned to before SPAdes, as a multiple of the
#: bin's assembled length. ``--max-coverage 0`` turns thinning off.
#:
#: This setting affects run time. An abundant organism can recruit hundreds of times the coverage
#: SPAdes needs, and SPAdes time grows with read count, so a few abundant bins set the step's run
#: time. 50x is above the depth where short-read assembly of a prokaryotic genome stops improving.
#:
#: Measured against the bin's assembled length, not an assumed genome size, since a 1.2 Mbp bin
#: and a 6 Mbp bin need different read counts for the same depth.
DEFAULT_MAX_COVERAGE = 50.0
#: Seed for the thinning. Fixed, so a re-run keeps the same reads and reassembly is reproducible.
SUBSAMPLE_SEED = 1337
#: Marks the third candidate: the original bin, copied beside its two reassemblies so one CheckM
#: pass scores all three. ``choose_best_bin`` reads this suffix off the stats row as the
#: candidate's style, and the summary counts it.
ORIG_SUFFIX = ".orig"
#: The two read-recruitment cut-offs every bin is reassembled at.
REASSEMBLY_STYLES = ("strict", "permissive")
#: Every candidate style, so the summary reports zero for a style that won no bins.
STYLES = REASSEMBLY_STYLES + (ORIG_SUFFIX.lstrip("."),)
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
    add_checkm2_only_args(p, "reassemble_bins")
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
    """Assembled length of every bin in *bins_dir*, keyed by bin name without extension."""
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

    Measures the fraction of bytes that are sequence in the first records and scales by the file
    size. Exact for uniform read length and header width (as recruited reads are), approximate
    otherwise.
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

    Pairs, not reads: both files are read in lockstep and one random draw decides both mates, so
    SPAdes still gets a paired library.

    Written to ``.thinned`` files and moved over the originals, so an interruption leaves the
    originals intact.
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

    Runs after recruitment and before SPAdes, over the per-bin fastqs the filter wrote. Bins
    already below the target (most of them) are not touched.
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
        # Estimated, not counted, to avoid reading every recruited read twice. The ratio of
        # sequence bytes to total bytes is measured on the first records, because in FASTQ the
        # sequence is under half the bytes and the rest depends on read-name length.
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
    """Fill *original_bins* with clean bin files named BIN_EXTENSION (MDMcleaner layout too)."""
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

    Only the bins, not the unbinned catch-all. The read filter assigns each pair to the bin of
    the contig it mapped to, so including unbinned contigs would recruit reads to a bin called
    "unbinned" and reassemble it as a genome. Misrecruitment is handled by the mate-agreement
    filter.
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
        # A per-bin SPAdes failure is tolerated (check=False) so one bin does not stop the run.
        # Show the errors from that bin's spades.log.
        warning("Something went wrong with reassembling %s.%s" % (bin_name, _spades_error(asm_dir)))


def _spades_error(asm_dir: str, lines: int = 6) -> str:
    """The ERROR lines from a failed SPAdes run, for inclusion in our own message."""
    log = os.path.join(asm_dir, "spades.log")
    if not os.path.isfile(log):
        return " (no spades.log was written; SPAdes may not have started)"
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


#: Intermediates moved into ``work_files`` in the finished layout. Used by both the move and the
#: --force clean.
INTERMEDIATES = (
    "reassembled_bins",
    "reassembled_bins.checkm",
    "reassembled_bins" + STATS_SUFFIX,
    "reads_for_reassembly",
    "binned_assembly",
    "reassemblies",
)


def _shuffle_intermediates(out: str, best_dir: str, reassembled: str) -> None:
    """Move the candidates into ``work_files``; promote the chosen bins to ``reassembled_bins``."""
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    for item in INTERMEDIATES:
        src = os.path.join(out, item)
        if not os.path.exists(src):
            continue
        # Replaces an existing destination, which a second reassembly into the same folder has
        # from the first run. Without this, the move failed after every reassembly was done.
        fs_replace_dir(src, os.path.join(work, item))
    os.replace(best_dir, reassembled)


def _clean_previous_run(out: str) -> None:
    """Clear a finished reassembly out of *out* so this run is not measured against it.

    The reassembly step adds to ``reassembled_bins`` instead of replacing it, so without this a
    second run's candidates would include the first run's chosen bins and ``.orig`` copies.

    Callers run this only under ``--force``; ``--resume`` keeps the existing folder.
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


def _stage_original_candidates(original_bins: str, reassembled: str) -> int:
    """Copy each original bin beside its reassemblies as the ``.orig`` candidate. Returns the count.

    The candidate name must end in :data:`BIN_EXTENSION`, because one CheckM run scores the
    directory with ``-x fasta`` (see :func:`metawrap2.checkm.run_checkm`). A hardcoded ``.orig.fa``
    left the originals unscored, so ``choose_best_bin`` could only pick a reassembly; on a 26-bin
    run, 4 bins (15%) were replaced by a worse reassembly.
    """
    count = 0
    for base in list_bins(original_bins, (BIN_EXTENSION,)):
        shutil.copy(
            os.path.join(original_bins, base),
            os.path.join(reassembled, os.path.splitext(base)[0] + ORIG_SUFFIX + BIN_EXTENSION),
        )
        count += 1
    return count


def _style_counts(names: Iterable[str]) -> Dict[str, int]:
    """How many of the chosen ``<bin>.<style>`` names each style won, keyed by style.

    Read from the style suffix, not by substring match: a bin called ``bin_orig_site2`` matched all
    three substring tests.
    """
    counts: Dict[str, int] = {style: 0 for style in STYLES}
    for name in names:
        style = name.rpartition(".")[2]
        counts[style] = counts.get(style, 0) + 1
    return counts


def _summary_message(names: Sequence[str]) -> str:
    """The closing summary: how the chosen bins split across the three candidates.

    *names* is what ``choose_best_bin`` selected. The total can be smaller than the number of input
    bins: a bin whose every version fails ``-c``/``-x`` is not selected and is in none of the
    counts.
    """
    counts = _style_counts(names)
    return (
        "Reassembly finished. Of %d bins, %d were improved by the 'strict' reassembly and %d "
        "by the 'permissive' reassembly; %d kept their original bin, which scored better than "
        "either reassembly." % (len(names), counts["strict"], counts["permissive"], counts["orig"])
    )


def _score(args, bins_dir: str, checkm2_env: Optional[str]) -> str:
    """Score *bins_dir* with CheckM2 and return its ``.stats`` path."""
    return _checkm.run_checkm(
        bins_dir,
        threads=args.threads,
        mem_gb=args.memory,
        env=None,
        checkm2=True,
        checkm2_env=checkm2_env,
    )


def _run_checkm_stage(args, out: str, checkm2_env: Optional[str] = None) -> None:
    announcement("RUN CHECKM ON REASSEMBLED BINS")
    reassembled = os.path.join(out, "reassembled_bins")
    original_bins = os.path.join(out, "original_bins")

    if dry_run():
        ensure_dir(reassembled)
        stats = _score(args, reassembled, checkm2_env)
        run(
            CHOOSE_BEST_BIN.format(
                stats=stats, comp=int(args.completeness), cont=int(args.contamination)
            ),
            env=None,
            tool="choose_best_bin",
            stdout_path=os.path.join(out, "chosen_bins.txt"),
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
            check=False,
        )
        return

    staged = _stage_original_candidates(original_bins, reassembled)
    comm("staged %d original bin(s) for scoring as the '%s' candidate" % (staged, ORIG_SUFFIX))

    comm("Running CheckM2 on the original and reassembled bins")
    stats = _score(args, reassembled, checkm2_env)

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

    announcement(_summary_message(names))

    best = os.listdir(best_dir)
    if not best:
        error(
            "no good bins found in %s; choosing the best version of each bin "
            "from the reassemblies failed." % best_dir
        )

    _shuffle_intermediates(out, best_dir, reassembled)

    comm("Re-running CheckM on the best reassembled bins.")
    shutil.rmtree(reassembled + ".checkm", ignore_errors=True)
    stats = _score(args, reassembled, checkm2_env)

    comm("making reassembly N50, completion, and contamination summary plots.")
    orig_stats = os.path.join(out, "original_bins" + STATS_SUFFIX)
    work_stats = os.path.join(out, "work_files", "reassembled_bins" + STATS_SUFFIX)
    with open(work_stats) as fh, open(orig_stats, "w") as out_fh:
        lines = fh.readlines()
        if lines:
            out_fh.write(lines[0])
        for line in lines[1:]:
            # Matched on the bin name's style suffix, not on "orig" anywhere in the line, since the
            # lineage column and bin names are free text.
            if line.split("\t", 1)[0].endswith(ORIG_SUFFIX):
                out_fh.write(line)
    # check=False: the bins, both .stats.tsv tables and the CheckM output are complete by now, so a
    # plotting failure only warns.
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
        check=False,
    )
    if not os.path.isfile(os.path.join(out, "reassembly_results.png")):
        warning(
            "Could not make the reassembly summary plots. The bins in %s and "
            "the stats tables beside them are complete; only the figure is missing." % reassembled
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
    # CheckM2's env and database are checked now, before hours of SPAdes.
    checkm2_env = scorer_preflight(args, settings) if args.run_checkm else None
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
            # Under --force only; --resume keeps the existing folder.
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
                    "Aligning all reads to the assembly and splitting them into one "
                    "fastq file per bin"
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
                    "WARNING: the assembly was already indexed. Skipping indexing and read "
                    "splitting and going straight to assembly. To redo them, re-run "
                    "with a new or empty output folder."
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
                # the input bins. FASTA_EXTENSIONS matches the input validation, so ``.fna`` and
                # ``.fas`` bins are listed too.
                to_assemble = [
                    "%s.%s_1.fastq" % (os.path.splitext(b)[0], style)
                    for b in list_bins(args.bins, FASTA_EXTENSIONS)
                    for style in REASSEMBLY_STYLES
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
                        "%s was reassembled but gave no contigs of at least %d bp, possibly "
                        "because there were too few reads." % (name, args.min_len)
                    )

            if _count_fa(reassembled) < 1 and not dry_run():
                error(
                    "None of the bins were successfully reassembled, so %s is empty.\n"
                    "Every per-bin SPAdes run failed; see the 'SPAdes reported' lines "
                    "above, or %s/<bin>/spades.log. A shared cause is likely; "
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
            _run_checkm_stage(args, out, checkm2_env)
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
