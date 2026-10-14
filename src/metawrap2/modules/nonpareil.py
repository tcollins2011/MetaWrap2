"""MetaWrap2 nonpareil: how much of each community was actually sequenced, and can these samples
be compared?

Every cross-sample statement a metagenomics study makes - this taxon is enriched, this pathway is
absent, this sample is more diverse - assumes the samples were seen to a similar extent. Depth in
gigabases does not establish that: a 10 Gbp soil sample and a 10 Gbp gut sample sit at very
different coverage of their own diversity, and what looks absent in the first may simply be
unsampled. Nonpareil answers it from the reads alone, with no assembly and no reference, and it is
the cheapest step in the pipeline by a wide margin.

It runs per sample on **one** read file - the forward mate, not the pair. That is Nonpareil's own
requirement rather than a simplification here: redundancy between a read and its own mate is an
artefact of the library, and feeding both would measure the insert size as if it were community
structure.

Three outputs. A per-sample table of coverage, the Nd diversity index and the projected effort to
reach a target coverage; the rarefaction points behind those estimates, so the fit can be looked
at; and a short verdict on comparability across whatever samples were given, which is the reason
this belongs in a pipeline rather than in a notebook.

See :mod:`metawrap2.nonpareil` for the model, its port from the R package, and the caveat that
coverage is abundance-weighted.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Optional, Tuple

from .. import nonpareil as _np
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    require_file,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["nonpareil"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# -T kmer rather than alignment: the k-mer kernel is Nonpareil 3's own default recommendation and is
# orders of magnitude faster, which is what makes this affordable on every sample rather than on one.
# The alignment kernel remains reachable with --kernel for very short or heterogeneous reads, where
# the k-mer kernel's fixed k is the weaker assumption.
#
# -X and -n are left at the binary's defaults deliberately. They set how many query reads and how
# many subsamples per point, which is to say how precise each point of the curve is - and the model
# weights points by their standard deviation, so changing them changes the fit, not just the runtime.
NONPAREIL = (
    "nonpareil -s {reads} -T {kernel} -f {format} -b {prefix} -t {threads} -R {ram} -d {divide} "
    "-k {k}{seed}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: Nonpareil's own default RAM ceiling, in MiB. Raised by default here: the binary refuses to grow
#: past it and reports a truncated curve instead, and 1 GiB is small for a modern read set.
DEFAULT_RAM_MB = 8192
#: Passed as ``-r`` under the alignment kernel, which otherwise seeds from the clock. Measured: two
#: runs over the same read file returned 83.0% and 84.3% coverage, a 1.3-point difference that is
#: purely the subsampling draw. The same value the binners use, for the same reason - a number that
#: moves between identical runs cannot be compared with anything, including itself.
#:
#: The k-mer kernel ignores ``-r``: the binary's help says the seed is "currently only implemented
#: when -T alignment", so under that kernel this is passed and has no effect, and reproducibility
#: there is not something this module can promise.
DEFAULT_SEED = 1337
#: Logarithmic subsampling factor (-d). The binary's default, and the model's fallback weightings
#: are chosen by whether this is zero, so changing it changes which weighting sequence is tried.
DEFAULT_DIVIDE = 0.7
#: One row per sample: coverage, diversity, and what more sequencing would buy.
SUMMARY_TABLE = "diversity_coverage" + TSV_SUFFIX
#: Every rarefaction point, with the fitted value beside the observed one.
CURVES_TABLE = "diversity_curves" + TSV_SUFFIX
#: Forward-mate suffixes stripped off a filename to get a sample name, longest match first. The
#: first entry is what read_qc writes, so a QC'd read set names its sample the same way here as it
#: does in every other module's output.
READ_SUFFIXES = ("_final_pure_reads_1", "_1", "_R1", ".1", "_fwd")
#: Where the binary's own .npo files are kept. The summary is derived from them, and a derived table
#: should not be the only copy: re-fitting is free, re-running the binary is not.
NPO_DIR = "nonpareil"
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 nonpareil",
        usage="metawrap2 nonpareil [options] -o output_dir -r sample1_1.fastq [sample2_1.fastq ...]",
        description="Diversity coverage per sample, and whether the samples are comparable.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-r",
        "--reads",
        required=True,
        nargs="+",
        metavar="FILE",
        help="one read file per sample - the forward mate only, gzip accepted. Nonpareil takes "
        "single-end input: redundancy between a read and its own mate measures the library, not "
        "the community. QC'd reads, not raw ones: it assumes error well below 5%%",
    )
    p.add_argument(
        "--kernel",
        choices=("auto", "kmer", "alignment"),
        default="auto",
        help="redundancy kernel (default: %(default)s). auto measures the read lengths and picks "
        "one kernel for the whole cohort: kmer, which is far faster, unless a sample's mean read "
        "length is under twice k, where the kmer kernel underestimates coverage badly (50.5%% "
        "against the alignment kernel's 83.0%% on one real 34 bp sample)",
    )
    p.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=_np.DEFAULT_K,
        metavar="K",
        help="k-mer size for the kmer kernel (default: %(default)s, the binary's own). Reads "
        "shorter than this cannot hold a k-mer and are filtered out before it runs",
    )
    p.add_argument(
        "--no-filter",
        action="store_true",
        help="do not filter reads shorter than k. The kmer kernel then fails on any read set that "
        "contains one, which trimmed reads routinely do",
    )
    p.add_argument(
        "--format",
        choices=("fastq", "fasta"),
        default=None,
        help="input format (default: inferred from the file name)",
    )
    p.add_argument(
        "--target-coverage",
        type=float,
        default=_np.STAR,
        metavar="PERCENT",
        help="coverage to project the required effort for (default: %(default)s). Nonpareil's own "
        "default, so the projection means the same thing as one quoted in a paper",
    )
    p.add_argument(
        "--max-ram",
        type=int,
        default=DEFAULT_RAM_MB,
        metavar="MiB",
        help="RAM ceiling passed to nonpareil (default: %(default)s). It does not grow past this, "
        "it truncates the curve instead, so too low a value is reported as a short curve",
    )
    p.add_argument(
        "--subsample-factor",
        type=float,
        default=DEFAULT_DIVIDE,
        metavar="D",
        help="logarithmic subsampling factor -d (default: %(default)s). Raise it for more points "
        "on the curve when the model does not converge",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        metavar="N",
        help="random seed for the alignment kernel's subsampling (default: %(default)s). Without "
        "one nonpareil seeds from the clock, and two runs over the same reads returned 83.0%% and "
        "84.3%% coverage. The kmer kernel does not implement it",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _sample_name(path: str) -> str:
    """Strip the read-file conventions off a path to get a sample name.

    ``sample_1.fastq.gz`` is ``sample``. Done here rather than left as the basename because the
    name ends up as a row label in a table joined against other modules' tables, where
    ``sample_1.fastq.gz`` would not match anything.
    """
    name = os.path.basename(path)
    for extension in (".gz", ".bz2"):
        name = name.removesuffix(extension)
    name = os.path.splitext(name)[0]
    # Longest first. read_qc writes `<sample>_final_pure_reads_1.fastq`, and stripping "_1" from
    # that leaves `gut_final_pure_reads` as the sample name - which joins against nothing.
    for suffix in sorted(READ_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _infer_format(path: str, given: Optional[str]) -> str:
    if given:
        return given
    name = path.lower()
    if name.endswith((".gz", ".bz2")):
        name = os.path.splitext(name)[0]
    if name.endswith((".fastq", ".fq")):
        return "fastq"
    if name.endswith((".fasta", ".fa", ".fna", ".fsa")):
        return "fasta"
    error(
        "cannot tell whether %s is fastq or fasta from its name; pass --format. Guessing wrong "
        "makes nonpareil read quality lines as sequence" % path
    )
    raise AssertionError("unreachable")  # pragma: no cover - error() exits


def _choose_kernel(args, lengths: Dict[str, _np.ReadLengths]) -> str:
    """One kernel for the whole cohort, from the measured read lengths.

    Not per sample, which is the part that matters. The two kernels do not agree: on one real gut
    sample (mean read length 34 bp) the k-mer kernel reports 50.5% coverage and the alignment
    kernel 83.0%, from the same reads. Picking per sample would therefore make a cohort's coverage
    figures incomparable - which is the one thing this module exists to check - so the whole cohort
    gets whichever kernel suits its shortest-read member, and the choice is stated.
    """
    if args.kernel != "auto":
        return args.kernel
    floor = _np.SHORT_READ_K_MULTIPLE * args.kmer_size
    short = {name: ln for name, ln in lengths.items() if ln.mean and ln.mean < floor}
    if not short:
        comm(
            "kernel: kmer (every sample's mean read length is at least %.0f bp, %gx k)"
            % (floor, _np.SHORT_READ_K_MULTIPLE)
        )
        return "kmer"
    comm(
        "kernel: alignment, because %d of %d sample(s) have a mean read length under %.0f bp (%s). "
        "At those lengths a read carries few k-mers and the k-mer kernel underestimates coverage - "
        "measured at 50.5%% against the alignment kernel's 83.0%% on one 34 bp sample. One kernel "
        "for the whole cohort, since the two do not produce comparable numbers. Force either with "
        "--kernel"
        % (
            len(short),
            len(lengths),
            floor,
            ", ".join("%s %.0f bp" % (n, ln.mean) for n, ln in sorted(short.items())),
        )
    )
    return "alignment"


def _prepare(
    args, kernel: str, name: str, reads: str, workdir: str, lengths: _np.ReadLengths
) -> Tuple[str, Optional[str]]:
    """Return the file to hand the binary, and the format to declare for it.

    Usually the input itself. For the k-mer kernel it is a filtered copy when the read set holds
    reads shorter than k, because the binary refuses the whole file over one of them - see
    :func:`metawrap2.nonpareil.filter_short_reads`.
    """
    for line in _np.kernel_advice(lengths, kernel, k=args.kmer_size):
        warning("%s: %s" % (name, line))
    if kernel != "kmer" or args.no_filter or not lengths.shorter_than_k:
        return reads, None
    filtered = os.path.join(workdir, "%s.minlen%d.fasta.gz" % (name, args.kmer_size))
    report = _np.filter_short_reads(reads, filtered, minimum_length=args.kmer_size)
    comm(
        "%s: filtered out %d of %d reads shorter than %d bp (%.2f%% of reads, %.2f%% of bp). "
        "Effort and the projection below are against the %.4g bp that remain"
        % (
            name,
            report.dropped_reads,
            report.dropped_reads + report.kept_reads,
            args.kmer_size,
            report.fraction_reads_dropped * 100,
            report.fraction_bp_dropped * 100,
            report.kept_bp,
        )
    )
    return filtered, "fasta"


def _run_one(
    args, env, kernel: str, name: str, reads: str, workdir: str, lengths: _np.ReadLengths
) -> Optional[str]:
    """Run the binary for one sample; return the .npo path, or None if it failed."""
    prefix = os.path.join(workdir, name)
    npo = prefix + ".npo"
    prepared, forced_format = _prepare(args, kernel, name, reads, workdir, lengths)
    rc = run(
        NONPAREIL.format(
            reads=prepared,
            kernel=kernel,
            format=forced_format or _infer_format(reads, args.format),
            prefix=prefix,
            threads=args.threads,
            ram=args.max_ram,
            divide=args.subsample_factor,
            k=args.kmer_size,
            # Only under the kernel that implements it, so that a command line read back out of
            # run_commands.txt does not suggest the k-mer run was reproducible.
            seed=" -r %d" % args.seed if kernel == "alignment" else "",
        ),
        env=env,
        tool="nonpareil",
        check=False,
    )
    if rc != 0 or not os.path.isfile(npo):
        warning(
            "nonpareil failed on %s, so that sample has no coverage estimate. The run continues "
            "with the others; a missing sample is visible in %s as an absent row rather than as a "
            "zero" % (name, SUMMARY_TABLE)
        )
        return None
    return npo


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    absolutize_paths(args)
    validate_inputs(
        "nonpareil",
        [(require_file, path, "read file (-r)") for path in args.reads],
    )
    resolve_threads(args, settings)

    rec = start_run("nonpareil", args, None, settings, inputs=list(args.reads))
    try:
        announcement("BEGIN NONPAREIL PIPELINE!")
        ensure_dir(args.output)
        workdir = os.path.join(args.output, NPO_DIR)
        ensure_dir(workdir)
        env = env_for("nonpareil", settings)

        samples: List[Tuple[str, str]] = [(_sample_name(p), p) for p in args.reads]
        duplicates = {n for n, _ in samples if [x for x, _ in samples].count(n) > 1}
        if duplicates:
            error(
                "two read files reduce to the same sample name (%s). Rename them: the name is a "
                "row label, and two rows with one label cannot be told apart"
                % ", ".join(sorted(duplicates))
            )

        if dry_run():
            comm(
                "would estimate diversity coverage for %d sample(s) with the %s kernel"
                % (len(samples), args.kernel)
            )
            announcement("NONPAREIL PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        # Every sample is scanned before anything runs, because the kernel is chosen for the cohort
        # rather than per sample. The scan reads the head of each file and costs about a second.
        lengths: Dict[str, _np.ReadLengths] = {}
        for name, reads in samples:
            lengths[name] = _np.scan_read_lengths(reads, k=args.kmer_size)
            comm(
                "%s: %d reads scanned, length %d-%d bp (mean %.0f)"
                % (
                    name,
                    lengths[name].scanned,
                    lengths[name].minimum,
                    lengths[name].maximum,
                    lengths[name].mean,
                )
            )
        kernel = _choose_kernel(args, lengths)

        curves: List[_np.Curve] = []
        for name, reads in samples:
            comm("%s: measuring redundancy" % name)
            npo = _run_one(args, env, kernel, name, reads, workdir, lengths[name])
            if npo is None:
                continue
            try:
                curve = _np.curve_from_npo(npo, label=name, star=args.target_coverage)
            except _np.NonpareilUnavailable as problem:
                error(str(problem))
                return 1
            except ValueError as problem:
                warning("%s: %s" % (name, problem))
                continue
            curves.append(curve)
            if not curve.consistent:
                for note in curve.warnings:
                    warning("%s: %s" % (name, note))
            projection = (
                "%.3gx more sequencing would reach %.0f%%" % (fold, curve.star)
                if (fold := curve.fold_more_needed)
                else "no projection (the model did not converge)"
            )
            comm(
                "%s: coverage %.1f%%, Nd %.2f, %s"
                % (name, curve.C * 100, curve.diversity, projection)
            )

        if not curves:
            error("no sample produced a curve, so there is nothing to report")
            return 1

        summary = _np.write_summary(curves, os.path.join(args.output, SUMMARY_TABLE))
        points = _np.write_curve_points(curves, os.path.join(args.output, CURVES_TABLE))

        for line in _np.comparability(curves):
            comm(line)
        comm("-> %s" % summary)
        comm("-> %s" % points)
        announcement("NONPAREIL PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
