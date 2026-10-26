"""Diversity coverage of a read set, estimated with Nonpareil.

Two samples at the same depth can have 90% and 40% coverage of their own diversity, which affects
any cross-sample comparison of richness, abundance or presence of a function. Nonpareil
(Rodriguez-R & Konstantinidis 2014; Rodriguez-R et al. 2018) estimates coverage from the reads,
without a reference or assembly, by measuring how often a read is redundant with reads already
seen and extrapolating the curve.

It gives three numbers:

* **C**, abundance-weighted average coverage at the depth sequenced. 0 to 1.
* **Nd**, the Nonpareil sequence-diversity index: log-scale diversity, comparable across samples
  (unlike richness estimated from unassembled reads).
* **LRstar**, the sequencing effort in base pairs projected to reach a target coverage (95% by
  default).

The redundancy curve comes from the `nonpareil` binary. The model fit, a weighted non-linear least
squares fit of a gamma CDF with a fallback sequence of weightings, is in Nonpareil's R package and
is ported here (:func:`fit_model`) so the pipeline can use the estimates directly. The port follows
Nonpareil 3.5.5's ``utils/Nonpareil/R/Nonpareil.R`` and keeps its parameter names.

On Nonpareil's three bundled datasets (HumanGut, IowaSoil, LakeLanier) this code reproduces R's C,
kappa, Nd, LRstar, model r and both gamma parameters to six significant figures. A test checks
this; the fixtures record the R invocation that produced the expected values.

Coverage is abundance-weighted. C=0.95 means 95% average coverage of the sequences weighted by
abundance, not 95% of the species; the missing part is the rare tail.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Target coverage in percent. Nonpareil's default, as used in its publications.
STAR = 95.0

#: Weight exponents applied to each point's standard deviation, tried in order until the fit
#: converges. Separate sequences for log and linear subsampling. Nonpareil's values and order.
WEIGHTS_LOG: Tuple[float, ...] = (0.0, 1.0, -1.0, 1.3, -1.1, 1.5, -1.5, 3.0, -3.0)
WEIGHTS_LINEAR: Tuple[float, ...] = (-1.1, -1.2, -0.9, -1.3, -1.0)

#: Printed with every coverage figure, since coverage is abundance-weighted.
COVERAGE_CAVEAT = (
    "Nonpareil coverage is the abundance-weighted average coverage of the community, not the "
    "fraction of its species seen. The uncovered part is the rare tail, so C=0.95 does not "
    "mean 95% of the lineages are present in the data"
)

#: Coverage below which a sample is flagged as not comparable. From Nonpareil's guidance.
LOW_COVERAGE = 0.5


#: Default k for the k-mer kernel (the binary's default). `nonpareil -T kmer` exits on reads shorter
#: than this ("Reads are required to have a minimum length of kmer size").
DEFAULT_K = 24

#: Mean read length, in multiples of k, below which the alignment kernel is advised (a 37 bp read
#: yields 14 24-mers). Nonpareil recommends the alignment kernel for short reads. Only a warning.
SHORT_READ_K_MULTIPLE = 2.0


class NonpareilUnavailable(RuntimeError):
    """Raised when scipy is not importable, so the gamma CDF fit cannot run."""


@dataclass(frozen=True)
class ReadLengths:
    """Read-length statistics from the first reads of a file, used to choose a kernel."""

    scanned: int
    minimum: int
    maximum: int
    mean: float
    shorter_than_k: int

    @property
    def fraction_short(self) -> float:
        return self.shorter_than_k / self.scanned if self.scanned else 0.0


def scan_read_lengths(path: str, k: int = DEFAULT_K, limit: int = 400_000) -> ReadLengths:
    """Return read-length statistics for the first *limit* reads.

    The head of the file is enough to choose a kernel and decide whether to filter; a full pass
    would double the module's runtime. The filtering pass itself reads every read.
    """
    lengths_seen = 0
    total = 0
    shortest: Optional[int] = None
    longest = 0
    short = 0
    for _header, sequence in _iter_reads(path):
        length = len(sequence)
        lengths_seen += 1
        total += length
        shortest = length if shortest is None else min(shortest, length)
        longest = max(longest, length)
        if length < k:
            short += 1
        if lengths_seen >= limit:
            break
    return ReadLengths(
        scanned=lengths_seen,
        minimum=shortest or 0,
        maximum=longest,
        mean=total / lengths_seen if lengths_seen else 0.0,
        shorter_than_k=short,
    )


def _open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        import gzip

        return gzip.open(path, "rt")
    return open(path)


def _iter_reads(path: str):
    """Yield (header, sequence) from fastq or fasta, gzipped or not.

    :mod:`metawrap2.io.seqio` reads only fasta. Fastq records are read in groups of four lines,
    since a quality line can also start with ``@``.
    """
    with _open_maybe_gzip(path) as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("@"):
            header = first.rstrip("\n")
            while True:
                sequence = handle.readline().rstrip("\n")
                handle.readline()  # the '+' separator
                handle.readline()  # the quality line
                if not sequence:
                    return
                yield header, sequence
                header = handle.readline().rstrip("\n")
                if not header:
                    return
        else:
            header = first.rstrip("\n")
            chunks: List[str] = []
            for line in handle:
                if line.startswith(">"):
                    yield header, "".join(chunks)
                    header, chunks = line.rstrip("\n"), []
                else:
                    chunks.append(line.strip())
            if chunks:
                yield header, "".join(chunks)


@dataclass(frozen=True)
class FilterReport:
    """Counts of reads and bases kept and removed by a minimum-length filter."""

    path: str
    kept_reads: int
    dropped_reads: int
    kept_bp: int
    dropped_bp: int
    minimum_length: int

    @property
    def fraction_reads_dropped(self) -> float:
        total = self.kept_reads + self.dropped_reads
        return self.dropped_reads / total if total else 0.0

    @property
    def fraction_bp_dropped(self) -> float:
        total = self.kept_bp + self.dropped_bp
        return self.dropped_bp / total if total else 0.0


def filter_short_reads(path: str, output: str, minimum_length: int = DEFAULT_K) -> FilterReport:
    """Write *path* to *output* as gzipped fasta, without reads shorter than *minimum_length*.

    The k-mer kernel rejects a file with any read shorter than k, and trimmed reads often have them
    (4.7% of one QC'd gut sample). Output is fasta because Nonpareil reads only the sequence, and
    gzipped to keep the temporary copy small.

    Effort is in base pairs, so effort and the projection are relative to the filtered total;
    the returned counts let the caller report this.
    """
    import gzip

    kept = dropped = 0
    kept_bp = dropped_bp = 0
    with gzip.open(output, "wt", compresslevel=1) as out:
        for header, sequence in _iter_reads(path):
            if len(sequence) < minimum_length:
                dropped += 1
                dropped_bp += len(sequence)
                continue
            kept += 1
            kept_bp += len(sequence)
            name = header[1:].split()[0] if header else "read%d" % kept
            out.write(">%s\n%s\n" % (name, sequence))
    return FilterReport(
        path=output,
        kept_reads=kept,
        dropped_reads=dropped,
        kept_bp=kept_bp,
        dropped_bp=dropped_bp,
        minimum_length=minimum_length,
    )


def kernel_advice(lengths: ReadLengths, kernel: str, k: int = DEFAULT_K) -> List[str]:
    """Return warnings about the chosen kernel given the measured read lengths, or none.

    The k-mer kernel suits typical read lengths but not heavily trimmed short reads.
    """
    lines: List[str] = []
    if kernel != "kmer":
        return lines
    if lengths.mean and lengths.mean < SHORT_READ_K_MULTIPLE * k:
        lines.append(
            "mean read length is %.0f bp against k=%d, so each read yields only ~%d k-mers. "
            "Nonpareil recommends the alignment kernel (--kernel alignment) for reads this short; "
            "the k-mer kernel still runs, on the reads long enough to hold a k-mer"
            % (lengths.mean, k, max(0, int(lengths.mean) - k + 1))
        )
    if lengths.shorter_than_k:
        lines.append(
            "%.1f%% of the first %d reads are shorter than k=%d. The k-mer kernel rejects a file "
            "that contains any such read, so they are filtered out before it runs"
            % (lengths.fraction_short * 100, lengths.scanned, k)
        )
    return lines


#: Values from the R package on Nonpareil's example data, matched to six significant figures.
R_AGREEMENT = {
    "HumanGut": {"C": 0.806352, "kappa": 0.790210, "diversity": 17.260378, "LRstar": 2.62214e09},
    "IowaSoil": {"C": 0.557296, "kappa": 0.527510, "diversity": 22.795958, "LRstar": 6.62396e11},
    "LakeLanier": {"C": 0.683087, "kappa": 0.659060, "diversity": 19.213483, "LRstar": 2.05296e10},
}


@dataclass
class Curve:
    """One sample's redundancy curve, its metadata, and the estimates fitted to it.

    Attribute names follow the R package's slots, so values can be checked against
    `Nonpareil.curve()` directly.
    """

    label: str
    #: Sequencing effort as reported in the .npo, in reads.
    x_obs: Tuple[float, ...]
    #: Observed redundancy, as measured by the binary.
    y_red: Tuple[float, ...]
    #: Redundancy transformed to coverage, and its spread.
    y_cov: Tuple[float, ...]
    y_sd: Tuple[float, ...]
    y_p50: Tuple[float, ...]
    #: Effort rescaled to base pairs, which is what LRstar is quoted in.
    x_adj: Tuple[float, ...]
    kernel: str
    #: Mean read length in bp, the .npo's `L`.
    read_length: float
    #: Length the effort axis is scaled by, the .npo's `AL`. With the alignment kernel it is the
    #: read length in bp; with the k-mer kernel it is k-mers per read (L - k + 1; 17.129 for
    #: L=40.129, k=24), so the effort axis and LRstar are in k-mers. Kept separate from
    #: `read_length`: mixing them up misreports the projected effort by about 2x on short reads.
    effort_length: float
    reads: float
    overlap: float
    log_sample: float
    version: str
    #: Redundancy at full effort.
    kappa: float
    #: Coverage at full effort.
    C: float
    star: float = STAR
    has_model: bool = False
    a: Optional[float] = None
    b: Optional[float] = None
    #: Nd. Zero when the fitted alpha is at most 1, where the index is undefined.
    diversity: float = 0.0
    #: Projected effort in bp for `star` coverage. None when it was not reached by extrapolation.
    LRstar: Optional[float] = None
    #: Pearson r between the fitted model and the rarefied data it was fitted to.
    model_r: Optional[float] = None
    consistent: bool = True
    warnings: List[str] = field(default_factory=list)

    @property
    def sequenced_bp(self) -> float:
        """Sequenced bases: reads x mean read length."""
        return self.reads * self.read_length

    @property
    def effort_unit(self) -> str:
        """The unit of the effort axis, and so of LRstar: "bp" or "kmers"."""
        return "kmers" if self.kernel == "kmer" else "bp"

    @property
    def effort_at_full_depth(self) -> float:
        """Return the end of the observed data on the effort axis, in :attr:`effort_unit`.

        Compare LRstar with this, not :attr:`sequenced_bp`; their units differ under the k-mer
        kernel.
        """
        return self.reads * self.effort_length

    @property
    def fold_more_needed(self) -> Optional[float]:
        """Return the multiple of current effort needed for `star` coverage, or None."""
        if self.LRstar is None or self.effort_at_full_depth <= 0:
            return None
        return self.LRstar / self.effort_at_full_depth


def _gamma_cdf(x: float, a: float, b: float) -> float:
    """`pgamma(log1p(x), shape=a, rate=b)` - Nonpareil's model function.

    Uses R's shape/rate parameterisation (scipy's scale = 1/rate), as in the published parameters.
    Using rate as scale gives a plausible curve with the wrong diversity index.
    """
    from scipy.stats import gamma  # local import: scipy is optional at install time

    return float(gamma.cdf(math.log1p(x), a, scale=1.0 / b))


def coverage_factor(overlap: float) -> float:
    """Return the exponent that converts redundancy to coverage, from Nonpareil's calibration."""
    return 1.0 - math.exp(2.23e-2 * overlap - 3.5698)


def read_npo(path: str, label: Optional[str] = None, correction: bool = True) -> Curve:
    """Parse a `.npo` file into a :class:`Curve`, with coverage and the adjusted effort axis.

    Equivalent to the R package's `read_metadata` and `read_data`. The model fit is separate, so
    a curve can still be reported when the fit fails, as it can on shallow data.
    """
    meta: Dict[str, str] = {}
    rows: List[Tuple[float, ...]] = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("# @"):
                key, _, value = line[3:].strip().partition(": ")
                meta[key] = value
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            rows.append(tuple(float(f) for f in fields[:6]))
    if not rows:
        raise ValueError("%s holds no curve: nonpareil produced no rarefaction points" % path)
    rows.sort(key=lambda r: r[0])

    # The .npo does not name the kernel; a k-mer size in range indicates the k-mer kernel, which
    # changes the overlap and the source of the effort length.
    ksize = float(meta.get("ksize", 0) or 0)
    kernel = "kmer" if 0 < ksize < 1001 else "alignment"
    log_sample = float(meta.get("logsampling", meta.get("divide", 0)) or 0)
    reads = float(meta["R"])
    length = float(meta["L"])
    if kernel == "kmer":
        overlap = 50.0
        # AL here is k-mers per read, not bp - see Curve.effort_length.
        effort_length = float(meta.get("AL", length))
    else:
        overlap = float(meta.get("overlap", 50))
        effort_length = length

    factor = coverage_factor(overlap) if correction else 1.0
    x_obs = tuple(r[0] for r in rows)
    y_red = tuple(r[1] for r in rows)
    y_cov = tuple(r[1] ** factor for r in rows)
    y_sd = tuple(r[2] ** factor for r in rows)
    y_p50 = tuple(r[4] ** factor for r in rows)
    kappa = y_red[-1]
    cov = y_cov[-1]

    # Effort is reported in reads but fitted in bp, with a coverage-dependent rescaling taken from
    # the R code (the 0.27 exponent is empirical). The first point is usually at zero reads; R
    # gives exp(log(0)) = 0 there, while Python raises, so zero is set directly.
    log_max = max(math.log(x) for x in x_obs if x > 0)
    x_adj = [
        math.exp(log_max + (cov**0.27) * (math.log(x) - log_max)) if x > 0 else 0.0 for x in x_obs
    ]
    scale = effort_length * reads / max(x_adj)
    x_adj = [x * scale for x in x_adj]

    curve = Curve(
        label=label or os.path.basename(path).replace(".npo", ""),
        x_obs=x_obs,
        y_red=y_red,
        y_cov=y_cov,
        y_sd=y_sd,
        y_p50=y_p50,
        x_adj=tuple(x_adj),
        kernel=kernel,
        read_length=length,
        effort_length=effort_length,
        reads=reads,
        overlap=overlap,
        log_sample=log_sample,
        version=meta.get("version", "unknown"),
        kappa=kappa,
        C=cov,
    )
    _check_consistency(curve)
    return curve


def _check_consistency(curve: Curve) -> None:
    """Record upstream's four checks for an unusable curve, each with the parameter to change.

    Failures are recorded as warnings, not raised.
    """
    half = 0.5 * curve.x_adj[-1]
    below = [i for i, x in enumerate(curve.x_adj) if x <= half]
    mid = max(below, key=lambda i: curve.x_adj[i]) if below else len(curve.x_adj) - 1
    if curve.y_p50[mid] == 0:
        curve.consistent = False
        curve.warnings.append(
            "the curve's median is zero at half the reads; re-run with a lower -L "
            "(alignment kernel) so reads overlap enough to register as redundant"
        )
    if curve.kappa <= 1e-5:
        curve.consistent = False
        curve.warnings.append(
            "redundancy is close to zero at full depth: the sample is too diverse for this "
            "sequencing depth, or -L is too high"
        )
    if len(curve.y_cov) > 1 and curve.y_cov[1] >= 1 - 1e-5:
        curve.consistent = False
        curve.warnings.append(
            "the curve saturates at the second point, so there is nothing to fit; re-run with a "
            "larger -d"
        )
    if sum(1 for y in curve.y_cov if 0 < y < 0.9) <= 10:
        curve.consistent = False
        curve.warnings.append(
            "fewer than 11 points are below 90% coverage, too few to fit; re-run with a "
            "larger -d"
        )


def fit_model(curve: Curve, weights_exp: Optional[Sequence[float]] = None) -> Curve:
    """Fit the gamma-CDF model and fill in Nd, LRstar and the model's r. Mutates and returns.

    As upstream, only points below 90% coverage are used: above that the curve is flat, and their
    small standard deviations would dominate the weighting.

    Each weighting in `weights_exp` is tried until one converges, as Nonpareil does; the first
    often fails on real data.
    """
    try:
        from scipy.optimize import curve_fit
        from scipy.stats import gamma
    except ImportError as exc:  # pragma: no cover - exercised by the module, not the unit tests
        raise NonpareilUnavailable(
            "scipy is needed to fit the Nonpareil model: pip install scipy, or install the "
            "metawrap2-reads environment, which includes it"
        ) from exc

    sel = [i for i, y in enumerate(curve.y_cov) if 0 < y < 0.9]
    if len(sel) < 3:
        curve.warnings.append(
            "too few points below 90% coverage to fit a model, so there is no diversity index or "
            "projected effort; the observed coverage above is measured, not fitted"
        )
        return curve
    xs = [curve.x_adj[i] for i in sel]
    ys = [curve.y_cov[i] for i in sel]
    sds = [curve.y_sd[i] for i in sel]

    if weights_exp is None:
        weights_exp = WEIGHTS_LOG if curve.log_sample != 0 else WEIGHTS_LINEAR

    def model(x: Any, a: float, b: float) -> Any:
        """The model for curve_fit: vectorised over x, with scipy's shape/scale."""
        import numpy as np

        return gamma.cdf(np.log1p(x), a, scale=1.0 / b)

    for exponent in weights_exp:
        # R's nls weights multiply the squared residual; curve_fit divides the residual by sigma,
        # so sigma = 1/sqrt(w). If any weight is zero or not finite, the fit is unweighted.
        try:
            weights = [sd**exponent if sd > 0 else 0.0 for sd in sds]
            if any(w <= 0 for w in weights) or any(not math.isfinite(w) for w in weights):
                sigma = None
            else:
                sigma = [1.0 / math.sqrt(w) for w in weights]
            popt, _ = curve_fit(
                model,
                xs,
                ys,
                p0=[1.0, 0.1],
                sigma=sigma,
                bounds=([0.0, 1e-12], [math.inf, math.inf]),
                maxfev=20000,
            )
        except (RuntimeError, ValueError, OverflowError):
            continue
        a, b = float(popt[0]), float(popt[1])
        if not (math.isfinite(a) and math.isfinite(b)) or b <= 0:
            continue
        curve.has_model = True
        curve.a, curve.b = a, b
        # Nd is undefined at alpha <= 1, so it is left at zero; (a-1)/b would be negative.
        curve.diversity = (a - 1.0) / b if a > 1.0 else 0.0
        try:
            curve.LRstar = float(math.exp(gamma.ppf(curve.star / 100.0, a, scale=1.0 / b)) - 1.0)
        except OverflowError:
            curve.LRstar = None
            curve.warnings.append(
                "the projected effort for %g%% coverage overflowed (the sample is far from that "
                "coverage); the diversity index is still usable" % curve.star
            )
        predicted = [_gamma_cdf(x, a, b) for x in xs]
        curve.model_r = _pearson(ys, predicted)
        return curve

    curve.warnings.append(
        "the model did not converge under any of Nonpareil's weightings, so there is no diversity "
        "index or projection; the observed coverage is still measured"
    )
    return curve


def _pearson(a: Sequence[float], b: Sequence[float]) -> Optional[float]:
    n = len(a)
    if n < 2:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((x - mb) ** 2 for x in b))
    if sa == 0 or sb == 0:
        return None
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


def curve_from_npo(path: str, label: Optional[str] = None, star: float = STAR) -> Curve:
    """Read a `.npo` file and fit the model."""
    curve = read_npo(path, label=label)
    curve.star = star
    return fit_model(curve)


#: Columns of the per-sample summary table, in order.
SUMMARY_COLUMNS: Tuple[str, ...] = (
    "sample",
    "kernel",
    "reads",
    "read_length_bp",
    "sequenced_bp",
    "effort_unit",
    "effort_at_full_depth",
    "redundancy",
    "coverage",
    "diversity_nd",
    "effort_for_target",
    "fold_more_for_target",
    "target_coverage",
    "model_r",
    "model_alpha",
    "model_beta",
    "consistent",
    "notes",
)


def summary_row(curve: Curve) -> Dict[str, str]:
    """Return the summary row for one sample. Values that were not fitted are empty."""

    def num(value: Optional[float], places: int = 4) -> str:
        if value is None:
            return ""
        return "%.*g" % (places + 2, round(value, places)) if abs(value) < 1e6 else "%.4g" % value

    return {
        "sample": curve.label,
        "kernel": curve.kernel,
        "reads": "%d" % curve.reads,
        "read_length_bp": num(curve.read_length, 1),
        "sequenced_bp": "%.4g" % curve.sequenced_bp,
        # With the k-mer kernel, effort and the projection are in k-mers, not base pairs.
        "effort_unit": curve.effort_unit,
        "effort_at_full_depth": "%.4g" % curve.effort_at_full_depth,
        "redundancy": num(curve.kappa),
        "coverage": num(curve.C),
        "diversity_nd": num(curve.diversity if curve.has_model and curve.diversity else None),
        "effort_for_target": num(curve.LRstar),
        "fold_more_for_target": num(curve.fold_more_needed, 2),
        "target_coverage": num(curve.star, 1),
        "model_r": num(curve.model_r),
        "model_alpha": num(curve.a),
        "model_beta": num(curve.b),
        "consistent": "yes" if curve.consistent else "no",
        "notes": "; ".join(curve.warnings),
    }


def write_summary(curves: Sequence[Curve], path: str) -> str:
    with open(path, "w") as out:
        out.write("\t".join(SUMMARY_COLUMNS) + "\n")
        for curve in curves:
            row = summary_row(curve)
            out.write("\t".join(row[column] for column in SUMMARY_COLUMNS) + "\n")
    return path


def write_curve_points(curves: Sequence[Curve], path: str) -> str:
    """Write every rarefaction point of every sample in long format, for plotting the curves.

    The last column holds the fitted coverage beside the observed value.
    """
    with open(path, "w") as out:
        out.write(
            "sample\teffort\teffort_unit\tredundancy\tcoverage\tcoverage_sd\tfitted_coverage\n"
        )
        for curve in curves:
            for i, effort in enumerate(curve.x_adj):
                fitted = (
                    "%.6g" % _gamma_cdf(effort, curve.a, curve.b)
                    if curve.has_model and curve.a is not None and curve.b is not None
                    else ""
                )
                out.write(
                    "%s\t%.6g\t%s\t%.6g\t%.6g\t%.6g\t%s\n"
                    % (
                        curve.label,
                        effort,
                        curve.effort_unit,
                        curve.y_red[i],
                        curve.y_cov[i],
                        curve.y_sd[i],
                        fitted,
                    )
                )
    return path


def comparability(curves: Sequence[Curve]) -> List[str]:
    """Return notes on comparing samples, given their coverage spread.

    Returns lines to print; low-coverage and spread warnings appear only when they apply.
    """
    lines: List[str] = []
    fitted = [c for c in curves if c.C > 0]
    if not fitted:
        return lines
    low = [c for c in fitted if c.C < LOW_COVERAGE]
    if low:
        lines.append(
            "%d of %d samples are below %.0f%% coverage (%s). Diversity and rare features cannot "
            "be compared between these and the other samples; a feature that looks absent may be "
            "unsampled."
            % (
                len(low),
                len(fitted),
                LOW_COVERAGE * 100,
                ", ".join(
                    "%s %.0f%%" % (c.label, c.C * 100) for c in sorted(low, key=lambda c: c.C)
                ),
            )
        )
    spread = max(c.C for c in fitted) - min(c.C for c in fitted)
    if spread >= 0.2:
        lines.append(
            "coverage ranges from %.0f%% to %.0f%% across samples, so differences in diversity "
            "between the best- and worst-covered samples partly reflect how much of each was seen."
            % (min(c.C for c in fitted) * 100, max(c.C for c in fitted) * 100)
        )
    lines.append(COVERAGE_CAVEAT)
    return lines
