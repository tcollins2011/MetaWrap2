"""MetaWrap2 nonpareil: per-sample diversity coverage, and whether the samples are comparable.

Cross-sample comparisons assume the samples were sequenced to a similar coverage of their own
diversity, which depth in Gbp does not show: a 10 Gbp soil sample and a 10 Gbp gut sample have very
different coverage. Nonpareil estimates coverage from the reads alone, without assembly or a
reference, and is the cheapest step in the pipeline.

It runs per sample on one read file, the forward mate. Nonpareil requires single-end input:
redundancy between mates reflects the library, not the community.

Outputs: a per-sample table of coverage, the Nd diversity index and the projected effort for a
target coverage; the rarefaction points behind those estimates; and a comparability note across
the samples given.

See :mod:`metawrap2.nonpareil` for the model, its port from the R package, and the caveat that
coverage is abundance-weighted.
"""

from __future__ import annotations

import argparse
import os
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

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
# The k-mer kernel is Nonpareil 3's recommended default and is orders of magnitude faster than the
# alignment kernel. The alignment kernel is available with --kernel for very short or
# heterogeneous reads.
#
# -X and -n stay at the binary's defaults. They set the query reads and subsamples per point, and
# the model weights points by their standard deviation, so changing them changes the fit.
NONPAREIL = (
    "nonpareil -s {reads} -T {kernel} -f {format} -b {prefix} -t {threads} -R {ram} -d {divide} "
    "-k {k}{seed}"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: RAM ceiling in MiB. Raised from Nonpareil's 1 GiB default: the binary does not grow past it and
#: truncates the curve instead.
DEFAULT_RAM_MB = 8192
#: Passed as ``-r`` under the alignment kernel, which otherwise seeds from the clock (two unseeded
#: runs on one read file gave 83.0% and 84.3% coverage). The same seed the binners use.
#:
#: The k-mer kernel does not implement ``-r`` ("currently only implemented when -T alignment"),
#: so its results are not reproducible.
DEFAULT_SEED = 1337
#: Logarithmic subsampling factor (-d), the binary's default. Whether it is zero selects the
#: model's fallback weighting sequence.
DEFAULT_DIVIDE = 0.7
#: One row per sample: coverage, diversity, and projected effort for the target coverage.
SUMMARY_TABLE = "diversity_coverage" + TSV_SUFFIX
#: Every rarefaction point, with the fitted value beside the observed one.
CURVES_TABLE = "diversity_curves" + TSV_SUFFIX
#: Forward-mate suffixes stripped from a filename to get a sample name, longest match first. The
#: first is what read_qc writes, so sample names match other modules' output.
READ_SUFFIXES = ("_final_pure_reads_1", "_1", "_R1", ".1", "_fwd")
#: Directory for the binary's .npo files, kept so the model can be re-fitted without re-running.
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
        help="one read file per sample: the forward mates only, gzip accepted. Nonpareil takes "
        "single-end input, since redundancy between a read and its mate measures the library "
        "and not the community. Use QC'd reads; it assumes an error rate well below 5%%",
    )
    p.add_argument(
        "--kernel",
        choices=("auto", "kmer", "alignment"),
        default="auto",
        help="redundancy kernel (default: %(default)s). auto measures the read lengths and picks "
        "one kernel for the whole cohort: kmer (much faster), unless a sample's mean read "
        "length is under twice k, where the kmer kernel underestimates coverage (50.5%% "
        "against 83.0%% with the alignment kernel on one 34 bp sample)",
    )
    p.add_argument(
        "-k",
        "--kmer-size",
        type=int,
        default=_np.DEFAULT_K,
        metavar="K",
        help="k-mer size for the kmer kernel (default: %(default)s, the nonpareil default). "
        "Reads shorter than this are filtered out before it runs",
    )
    p.add_argument(
        "--no-filter",
        action="store_true",
        help="do not filter reads shorter than k. The kmer kernel then fails on any read set "
        "that contains one, which trimmed reads usually do",
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
        help="coverage to project the required sequencing effort for (default: %(default)s, "
        "the Nonpareil default, so projections match published ones)",
    )
    p.add_argument(
        "--max-ram",
        type=int,
        default=DEFAULT_RAM_MB,
        metavar="MiB",
        help="RAM limit passed to nonpareil (default: %(default)s). When it reaches the limit "
        "it truncates the curve, so a value too low gives a short curve",
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
        "a seed nonpareil uses the clock; two such runs on the same reads gave 83.0%% and "
        "84.3%% coverage. The kmer kernel ignores it",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _sample_name(path: str) -> str:
    """Return the sample name for a read file: ``sample_1.fastq.gz`` -> ``sample``.

    The name is a row label joined against other modules' tables.
    """
    name = os.path.basename(path)
    for extension in (".gz", ".bz2"):
        name = name.removesuffix(extension)
    name = os.path.splitext(name)[0]
    # Longest suffix first, so `<sample>_final_pure_reads_1` loses the whole suffix, not only "_1".
    for suffix in sorted(READ_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def sample_names(paths: Sequence[str]) -> List[str]:
    """Return one label per read file, using the parent directory when basenames collide.

    `read_qc` writes `<out>/<sample>/final_pure_reads_1.fastq`, so all files share one basename
    and the sample is the directory name (`READ_QC/ERR011347/...` is ERR011347). Only colliding
    names are replaced. The caller reports collisions that remain.
    """
    names = [_sample_name(path) for path in paths]
    counts = Counter(names)
    for index, (name, path) in enumerate(zip(names, paths)):
        if counts[name] < 2:
            continue
        parent = os.path.basename(os.path.dirname(os.path.abspath(path)))
        if parent:
            names[index] = parent
    return names


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
        "cannot tell whether %s is fastq or fasta from its name; pass --format (with the "
        "wrong format nonpareil reads quality lines as sequence)" % path
    )
    raise AssertionError("unreachable")  # pragma: no cover - error() exits


def _choose_kernel(args, lengths: Dict[str, _np.ReadLengths]) -> str:
    """Choose one kernel for the whole cohort from the measured read lengths.

    The kernels give different results (on one gut sample with 34 bp mean read length, 50.5%
    coverage with k-mer and 83.0% with alignment), so a per-sample choice would make samples
    incomparable. The kernel suits the sample with the shortest reads, and the choice is printed.
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
        "Such reads hold few k-mers and the k-mer kernel underestimates coverage (50.5%% "
        "against 83.0%% with the alignment kernel on one 34 bp sample). One kernel is used for "
        "the whole cohort because the two give different numbers. Set one with "
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
    """Return the file to pass to the binary, and the format to declare for it.

    Usually the input itself. For the k-mer kernel, a filtered copy when some reads are shorter
    than k, since the binary rejects the whole file; see
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
        "Effort and the projection below refer to the %.4g bp that remain"
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
            # Only for the alignment kernel, which implements it, so run_commands.txt does not
            # suggest a seeded k-mer run.
            seed=" -r %d" % args.seed if kernel == "alignment" else "",
        ),
        env=env,
        tool="nonpareil",
        check=False,
    )
    if rc != 0 or not os.path.isfile(npo):
        warning(
            "nonpareil failed on %s, so that sample has no coverage estimate. The run continues "
            "with the other samples; the failed sample has no row in %s" % (name, SUMMARY_TABLE)
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

        labels = sample_names(args.reads)
        samples: List[Tuple[str, str]] = list(zip(labels, args.reads))
        duplicates = {name for name, count in Counter(labels).items() if count > 1}
        if duplicates:
            error(
                "two read files give the same sample name (%s), and their parent directories "
                "do not distinguish them. Rename the files, or put each sample's reads in its "
                "own directory; sample names label the rows of the output"
                % ", ".join(sorted(duplicates))
            )

        if dry_run():
            comm(
                "would estimate diversity coverage for %d sample(s) with the %s kernel"
                % (len(samples), args.kernel)
            )
            announcement("NONPAREIL PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        # Scan every sample first, since one kernel is chosen for the cohort. The scan reads the
        # head of each file and takes about a second.
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
