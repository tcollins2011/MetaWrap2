"""Diversity coverage from a read set: can these samples be compared at all?

Sequencing depth is not the same thing as how much of a community was seen. Two samples at the
same depth can sit at 90% and 40% coverage of their own diversity, and every cross-sample
comparison drawn from them - richness, abundance, presence of a function - is then comparing a
well-described community with a glimpse of one. Nonpareil (Rodriguez-R & Konstantinidis 2014;
Rodriguez-R et al. 2018) estimates that coverage from the reads themselves, with no reference
database and no assembly, by measuring how often a read is redundant with the reads already seen
and extrapolating the curve.

Three numbers come out of it:

* **C**, abundance-weighted average coverage at the depth actually sequenced. 0 to 1.
* **Nd**, the Nonpareil sequence-diversity index: the log-scale diversity of the dataset, which is
  comparable across samples in a way that a richness estimate from unassembled reads is not.
* **LRstar**, the sequencing effort in base pairs projected to reach a target coverage (95% by
  default) - the answer to "how much more would I have to sequence".

The redundancy curve comes from the `nonpareil` binary. The model fit does not: upstream it lives
in an R package, and the fit is a weighted non-linear least squares of a gamma CDF with a documented
fallback sequence of weightings. That is ported here (see :func:`fit_model`) rather than shelled out
to, because the estimates are wanted *in* the pipeline - a coverage figure decides what the run can
say about cross-sample comparison, and reaching that conclusion should not involve parsing a PDF
report. The port follows Nonpareil 3.5.5's ``utils/Nonpareil/R/Nonpareil.R`` closely enough that
the parameter names match, deliberately: anyone checking a number here against the R package should
be able to read the two side by side.

Bioconda's `nonpareil` brings the R package with it, so the original is available as a cross-check,
and it was used as one: on Nonpareil's three bundled datasets (HumanGut, IowaSoil, LakeLanier) this
code reproduces R's C, kappa, Nd, LRstar, model r and both gamma parameters to six significant
figures. The comparison is pinned as a test; the fixtures carry the R invocation that produced the
expected values.

What this does not do: Nonpareil's coverage is of the *community*, abundance-weighted. A sample at
C=0.95 has 95% average coverage of the sequences present in proportion to their abundance, which is
not 95% of its species - the rare tail is systematically the part that is missing. That caveat
travels with the number.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Coverage considered "nearly complete", as a percentage. Nonpareil's own default, and the value
#: its published figures use, so a projected effort reported here means the same thing as one
#: reported in a paper.
STAR = 95.0

#: Weight exponents tried in order until the fit converges, applied to each point's standard
#: deviation. Two sequences because the curve's shape differs between log and linear subsampling.
#: These are Nonpareil's own, kept in its order: the first that converges wins, so changing the
#: order changes which fit a borderline dataset gets.
WEIGHTS_LOG: Tuple[float, ...] = (0.0, 1.0, -1.0, 1.3, -1.1, 1.5, -1.5, 3.0, -3.0)
WEIGHTS_LINEAR: Tuple[float, ...] = (-1.1, -1.2, -0.9, -1.3, -1.0)

#: Printed beside every coverage figure. The number is abundance-weighted, and read as a fraction
#: of species it overstates what was seen by the amount that matters most.
COVERAGE_CAVEAT = (
    "Nonpareil coverage is abundance-weighted average coverage of the community, not the fraction "
    "of its species that were seen: the uncovered remainder is the rare tail, so C=0.95 does not "
    "mean 95% of the lineages are present in the data"
)

#: Printed when coverage is low enough that cross-sample comparison is the thing at risk, rather
#: than any particular estimate. Nonpareil's own guidance for well-covered data.
LOW_COVERAGE = 0.5


#: Default k for the k-mer kernel, which is also the binary's. Reads shorter than this are fatal to
#: that kernel: `nonpareil -T kmer` exits with "Reads are required to have a minimum length of kmer
#: size" rather than skipping them.
DEFAULT_K = 24

#: Below this many multiples of k, the k-mer kernel has little to work with per read - a 37 bp read
#: yields 14 24-mers - and Nonpareil's own guidance is the alignment kernel for short reads. Not a
#: hard refusal: it is a warning, because the user may know their data better than a ratio does.
SHORT_READ_K_MULTIPLE = 2.0


class NonpareilUnavailable(RuntimeError):
    """Raised when scipy is not importable, so the gamma CDF fit cannot run."""


@dataclass(frozen=True)
class ReadLengths:
    """What a scan of the first reads of a file found. Used to choose a kernel and to warn."""

    scanned: int
    minimum: int
    maximum: int
    mean: float
    shorter_than_k: int

    @property
    def fraction_short(self) -> float:
        return self.shorter_than_k / self.scanned if self.scanned else 0.0


def scan_read_lengths(path: str, k: int = DEFAULT_K, limit: int = 400_000) -> ReadLengths:
    """Read lengths from the first *limit* reads: enough to decide, cheap enough to always do.

    A scan rather than a full pass. The decisions it feeds - which kernel suits this data, whether
    a filtering pass is needed at all - do not become different decisions on read 400,001, and a
    full pass over every read file would double the module's runtime to answer a question the head
    of the file already answers. The filtering pass itself, when it happens, sees every read.
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

    Written here rather than reused from :mod:`metawrap2.io.seqio` because that module's reader is
    a fasta reader, and this has to handle fastq - where a line starting with ``@`` can also be a
    quality line, so records are counted in fours rather than detected by their first character.
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
    """What a minimum-length filter removed. Every field is reported, none is inferred."""

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

    Needed because the k-mer kernel refuses a file containing a single read shorter than k, and
    trimmed reads routinely contain them - 4.7% of one real QC'd gut sample. Dropping them is the
    only way to run that kernel at all, so it is done here and counted, rather than left as a
    fatal error with a manual fix.

    Fasta, not fastq: Nonpareil reads only the sequence, so carrying the quality lines would double
    the size of a file that is already a copy of the input. Gzipped for the same reason - the input
    is gigabytes and this is a temporary copy of most of it.

    The counts matter for interpretation, not just for reassurance. Sequencing effort is in base
    pairs, so an effort axis built from the filtered file is shorter than the library really was,
    and the projected effort to reach a target coverage is quoted against the filtered total.
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
    """What the measured read lengths say about the chosen kernel. Lines to print, or none.

    The k-mer kernel is the right default for modern read lengths and the wrong one for heavily
    trimmed short reads, and the file itself says which it is - so this is measured rather than
    left to the user to know.
    """
    lines: List[str] = []
    if kernel != "kmer":
        return lines
    if lengths.mean and lengths.mean < SHORT_READ_K_MULTIPLE * k:
        lines.append(
            "mean read length is %.0f bp against k=%d, so each read yields only ~%d k-mers. The "
            "alignment kernel (--kernel alignment) is Nonpareil's recommendation for reads this "
            "short; the k-mer kernel will still run, on the reads long enough to hold a k-mer"
            % (lengths.mean, k, max(0, int(lengths.mean) - k + 1))
        )
    if lengths.shorter_than_k:
        lines.append(
            "%.1f%% of the first %d reads are shorter than k=%d. The k-mer kernel rejects the "
            "whole file over one such read, so they are filtered out before it runs"
            % (lengths.fraction_short * 100, lengths.scanned, k)
        )
    return lines


#: Verified against the R package on Nonpareil's own example data, to six significant figures. Kept
#: here rather than only in the test so that the claim and the numbers travel with the code.
R_AGREEMENT = {
    "HumanGut": {"C": 0.806352, "kappa": 0.790210, "diversity": 17.260378, "LRstar": 2.62214e09},
    "IowaSoil": {"C": 0.557296, "kappa": 0.527510, "diversity": 22.795958, "LRstar": 6.62396e11},
    "LakeLanier": {"C": 0.683087, "kappa": 0.659060, "diversity": 19.213483, "LRstar": 2.05296e10},
}


@dataclass
class Curve:
    """One sample's redundancy curve, its metadata, and the estimates fitted to it.

    Attribute names follow the R package's slots rather than being renamed to house style, so that
    a number here can be checked against a number from `Nonpareil.curve()` without a translation
    table.
    """

    label: str
    #: Sequencing effort as reported in the .npo, in reads.
    x_obs: Tuple[float, ...]
    #: Observed *redundancy* - what the binary measures.
    y_red: Tuple[float, ...]
    #: Redundancy transformed to coverage, and its spread.
    y_cov: Tuple[float, ...]
    y_sd: Tuple[float, ...]
    y_p50: Tuple[float, ...]
    #: Effort rescaled to base pairs, which is what LRstar is quoted in.
    x_adj: Tuple[float, ...]
    kernel: str
    #: Mean read length in bp, the .npo's `L`. What was sequenced.
    read_length: float
    #: The length the effort axis is scaled by, the .npo's `AL`, and **not always a length in bp**:
    #: with the alignment kernel it is the read length, but with the k-mer kernel it is k-mers per
    #: read (L - k + 1; measured as 17.129 for L=40.129, k=24). So the effort axis - and LRstar with
    #: it - is in k-mers there, not base pairs. Kept as its own field rather than folded into
    #: `read_length` because conflating them silently misreports the projected effort by the ratio
    #: between them, which on short reads is a factor of two.
    effort_length: float
    reads: float
    overlap: float
    log_sample: float
    version: str
    #: Redundancy at full effort.
    kappa: float
    #: Coverage at full effort - the headline number.
    C: float
    star: float = STAR
    has_model: bool = False
    a: Optional[float] = None
    b: Optional[float] = None
    #: Nd. Zero when the fit put alpha at or below 1, where the index is undefined rather than low.
    diversity: float = 0.0
    #: Projected effort in bp for `star` coverage. None when it was not reached by extrapolation.
    LRstar: Optional[float] = None
    #: Pearson r between the fitted model and the rarefied data it was fitted to.
    model_r: Optional[float] = None
    consistent: bool = True
    warnings: List[str] = field(default_factory=list)

    @property
    def sequenced_bp(self) -> float:
        """What was actually sequenced, in bp: reads x mean read length."""
        return self.reads * self.read_length

    @property
    def effort_unit(self) -> str:
        """The unit of the effort axis, and so of LRstar: "bp" or "kmers"."""
        return "kmers" if self.kernel == "kmer" else "bp"

    @property
    def effort_at_full_depth(self) -> float:
        """Where the observed data ends on the effort axis, in :attr:`effort_unit`.

        This, not :attr:`sequenced_bp`, is what LRstar must be compared against: the two are in
        different units under the k-mer kernel.
        """
        return self.reads * self.effort_length

    @property
    def fold_more_needed(self) -> Optional[float]:
        """How many times the current effort would reach `star` coverage. None without a model.

        The ratio rather than the absolute projection, because that is the form the answer gets
        used in: "4x more sequencing" is actionable where "1.2e11 bp" needs dividing first.
        """
        if self.LRstar is None or self.effort_at_full_depth <= 0:
            return None
        return self.LRstar / self.effort_at_full_depth


def _gamma_cdf(x: float, a: float, b: float) -> float:
    """`pgamma(log1p(x), shape=a, rate=b)` - Nonpareil's model function.

    Written against the shape/rate parameterisation R uses, not scipy's shape/scale, because the
    published parameters are in rate. Mixing them silently returns a plausible curve with the wrong
    diversity index.
    """
    from scipy.stats import gamma  # local import: scipy is optional at install time

    return float(gamma.cdf(math.log1p(x), a, scale=1.0 / b))


def coverage_factor(overlap: float) -> float:
    """The exponent that turns redundancy into coverage, from Nonpareil's calibration."""
    return 1.0 - math.exp(2.23e-2 * overlap - 3.5698)


def read_npo(path: str, label: Optional[str] = None, correction: bool = True) -> Curve:
    """Parse a `.npo` file into a :class:`Curve`, with coverage and the adjusted effort axis.

    Does everything the R package's `read_metadata` and `read_data` do and nothing else; the model
    fit is separate so that a curve can be reported even when the fit fails, which is a real case
    on shallow data and is more useful than an exception.
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

    # The kernel is not stated directly; a k-mer size in range is what distinguishes the kmer
    # kernel's output, and it changes both the overlap and where the read length comes from.
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

    # Effort is reported in reads but the model is fitted in bp, and the rescaling is coverage
    # dependent - a poorly covered dataset's curve is stretched differently from a well covered
    # one's. Straight from the R, including the 0.27 exponent, which is empirical.
    # The first rarefaction point is usually at zero reads, where R's log(0) is -Inf and the
    # rescaled effort comes out as exp(-Inf) = 0. Python raises instead, so the zero is written
    # out directly rather than letting the arithmetic produce it.
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
    """Record the four ways upstream considers a curve unusable, each with its remedy.

    Reported rather than raised. A curve that fails one of these has a coverage number that should
    not be quoted, but the run that produced it took hours and the user needs to know which
    parameter to change, not a traceback.
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
            "redundancy is essentially zero at full depth: either the sample is far too diverse "
            "for this depth to say anything, or -L is too high"
        )
    if len(curve.y_cov) > 1 and curve.y_cov[1] >= 1 - 1e-5:
        curve.consistent = False
        curve.warnings.append(
            "the curve saturates at the second point, so there is nothing to fit; re-run with a "
            "larger -d to spread the subsampling out"
        )
    if sum(1 for y in curve.y_cov if 0 < y < 0.9) <= 10:
        curve.consistent = False
        curve.warnings.append(
            "fewer than 11 points fall below 90% coverage, which is too few to fit; re-run with a "
            "larger -d"
        )


def fit_model(curve: Curve, weights_exp: Optional[Sequence[float]] = None) -> Curve:
    """Fit the gamma-CDF model and fill in Nd, LRstar and the model's r. Mutates and returns.

    Fitted only to the points below 90% coverage, as upstream does: above that the curve is flat
    and the points carry no information about the shape, while weighting by their small standard
    deviations would let them dominate.

    Each weighting in `weights_exp` is tried until one converges. That is not a heuristic added
    here - it is Nonpareil's own behaviour, and on real data the first choice does fail.
    """
    try:
        from scipy.optimize import curve_fit
        from scipy.stats import gamma
    except ImportError as exc:  # pragma: no cover - exercised by the module, not the unit tests
        raise NonpareilUnavailable(
            "scipy is needed to fit the Nonpareil model: pip install scipy, or install the "
            "metawrap2-reads environment, which has it"
        ) from exc

    sel = [i for i, y in enumerate(curve.y_cov) if 0 < y < 0.9]
    if len(sel) < 3:
        curve.warnings.append(
            "too few points below 90% coverage to fit a model, so there is no diversity index or "
            "projected effort; the observed coverage above is still measured, not fitted"
        )
        return curve
    xs = [curve.x_adj[i] for i in sel]
    ys = [curve.y_cov[i] for i in sel]
    sds = [curve.y_sd[i] for i in sel]

    if weights_exp is None:
        weights_exp = WEIGHTS_LOG if curve.log_sample != 0 else WEIGHTS_LINEAR

    def model(x: Any, a: float, b: float) -> Any:
        """The model as curve_fit wants it: vectorised over x, and in shape/scale for scipy."""
        import numpy as np

        return gamma.cdf(np.log1p(x), a, scale=1.0 / b)

    for exponent in weights_exp:
        # R's nls weights multiply the squared residual; curve_fit divides the residual by sigma.
        # So sigma = 1/sqrt(w), and a zero standard deviation would divide by zero - those points
        # are dropped from the weighting rather than from the fit.
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
        # Nd is undefined at alpha <= 1 rather than small, so it is left at zero and said to be
        # absent. Reporting (a-1)/b there would give a negative diversity.
        curve.diversity = (a - 1.0) / b if a > 1.0 else 0.0
        try:
            curve.LRstar = float(math.exp(gamma.ppf(curve.star / 100.0, a, scale=1.0 / b)) - 1.0)
        except OverflowError:
            curve.LRstar = None
            curve.warnings.append(
                "the projected effort for %g%% coverage overflowed, which means this sample is far "
                "from it; the diversity index is still usable" % curve.star
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
    """Read a `.npo` and fit it in one call, which is what every caller wants."""
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
    """One table row per sample. Everything unfitted is empty rather than zero."""

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
        # Stated in every row rather than left to the reader to infer from the kernel: with the
        # k-mer kernel the next column and the projection are in k-mers, not base pairs.
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
    """Every rarefaction point of every sample, long format - the input for plotting the curves.

    Written because the estimates are a model of these points and a reader who wants to see whether
    the model fits has to have them. The fitted value is in the last column beside the observed
    one.
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
    """What the coverage spread across samples means for comparing them.

    This is the reason to run Nonpareil in a pipeline rather than as a one-off: a single sample's
    coverage is interesting, and the spread across a cohort decides whether any cross-sample
    statement is safe. Returns lines to print, empty when there is nothing to warn about.
    """
    lines: List[str] = []
    fitted = [c for c in curves if c.C > 0]
    if not fitted:
        return lines
    low = [c for c in fitted if c.C < LOW_COVERAGE]
    if low:
        lines.append(
            "%d of %d samples are below %.0f%% coverage (%s). Diversity and the presence of rare "
            "features are not comparable between these and the rest: what looks absent may be "
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
            "coverage spans %.0f%% to %.0f%% across samples. Differences in diversity between the "
            "best- and worst-covered samples are partly differences in how much of them was seen."
            % (min(c.C for c in fitted) * 100, max(c.C for c in fitted) * 100)
        )
    lines.append(COVERAGE_CAVEAT)
    return lines
