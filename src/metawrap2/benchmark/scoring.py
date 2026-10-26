"""Collect benchmark tool output into one long-format results table and rank variants.

Each number is one row of ``dataset x comparison x variant x sample x metric``. The long format
lets new metrics be added without rewriting older files.

Cost is stored as ordinary metrics (``wall_seconds``, ``cpu_hours``, ``peak_rss_bytes``) so it can
be joined against quality metrics. :func:`value_against_cost` does that join and is used to choose
defaults. It ranks on the quality metric and reports cost beside it. An earlier version ranked on
the quality/cost ratio, which favours the cheapest arm and picked a default that recovered fewer
MAGs than metaWRAP 1. Where cost must drive the choice, it is exposed as a user-selected tier
(``binning --preset``).

The scorers parse the output of CheckM2, GUNC and AMBER and do not reimplement their scoring.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .. import gunc as _gunc

#: The TSV's columns, in order. Do not reorder: older results files must stay readable.
COLUMNS = (
    "dataset",
    "comparison",
    "variant",
    "sample",
    "metric",
    "value",
    "unit",
    "source",
)

#: Bin quality tiers: ``(name, min completeness, max contamination)``.
#:
#: * ``bins_medium_quality`` (>=50 / <10) is MIMAG medium quality (Bowers et al. 2017).
#: * ``bins_near_complete`` (>=90 / <5) is MIMAG high quality without its rRNA and tRNA requirement
#:   (23S, 16S and 5S plus at least 18 tRNAs), so it is not called high quality here.
#: * ``bins_high_quality`` (>=70 / <10) is not a MIMAG tier. It is this project's reporting
#:   threshold, matching bin_refinement's defaults and common MAG-recovery practice, and is the
#:   headline number in docs/benchmarks.md. The name is kept because it is the metric key in all
#:   recorded results rows.
QUALITY_TIERS = (
    ("bins_near_complete", 90.0, 5.0),
    ("bins_high_quality", 70.0, 10.0),
    ("bins_medium_quality", 50.0, 10.0),
)

#: Unit of each metric, written to the TSV. Metrics not listed get an empty unit.
UNITS = {
    "wall_seconds": "s",
    "cpu_hours": "cpu_h",
    "peak_rss_bytes": "B",
    "gpu_seconds": "s",
    "assembly_n50": "bp",
    "assembly_total_bp": "bp",
    "gold_recall_bp": "bp",
    "amber_f1_bp": "fraction",
    "amber_ari_bp": "index",
    "amber_completeness_bp": "fraction",
    "amber_purity_bp": "fraction",
    "gunc_pass_fraction": "fraction",
    "reads_kept": "reads",
    "bases_kept": "bp",
    "q30_bases_kept": "bp",
    "abundance_missing_bins": "bins",
}

#: Units whose metrics are averaged across samples instead of summed (fractions and indexes are
#: already normalised). A new metric gets this pooling by declaring one of these units in UNITS.
_AVERAGED_UNITS = frozenset(("fraction", "index"))

#: Byte lookup marking Phred+33 quality characters of 30 or better (ASCII 63, ``?``, and up). Used
#: with :meth:`bytes.translate` so the Q30 count runs in C instead of a per-base Python loop.
_Q30_BYTES = bytes(1 if i >= 63 else 0 for i in range(256))


@dataclass(frozen=True)
class Row:
    """One measurement. ``sample`` is empty for a figure pooled across samples."""

    dataset: str
    comparison: str
    variant: str
    metric: str
    value: float
    sample: str = ""
    unit: str = ""
    source: str = ""

    def as_dict(self) -> Dict[str, Any]:
        return {
            "dataset": self.dataset,
            "comparison": self.comparison,
            "variant": self.variant,
            "sample": self.sample,
            "metric": self.metric,
            # Not rounded, so later calculations keep full precision.
            "value": ("%g" % self.value) if isinstance(self.value, float) else str(self.value),
            "unit": self.unit or UNITS.get(self.metric, ""),
            "source": self.source,
        }


def write_tsv(rows: Iterable[Row], path: str, *, append: bool = False) -> int:
    """Write *rows* to *path*. Returns the number written.

    With *append*, no header is written if the file already has content.
    """
    rows = list(rows)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    exists = os.path.isfile(path) and os.path.getsize(path) > 0
    with open(path, "a" if append else "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(COLUMNS), delimiter="\t", lineterminator="\n")
        if not (append and exists):
            writer.writeheader()
        for row in rows:
            writer.writerow(row.as_dict())
    return len(rows)


def read_tsv(path: str) -> List[Row]:
    """Read a results TSV. Rows whose value is not a number are skipped."""
    out: List[Row] = []
    with open(path, newline="") as fh:
        for record in csv.DictReader(fh, delimiter="\t"):
            try:
                value = float(record["value"])
            except (TypeError, ValueError, KeyError):
                continue
            out.append(
                Row(
                    dataset=record.get("dataset", ""),
                    comparison=record.get("comparison", ""),
                    variant=record.get("variant", ""),
                    metric=record.get("metric", ""),
                    value=value,
                    sample=record.get("sample", ""),
                    unit=record.get("unit", ""),
                    source=record.get("source", ""),
                )
            )
    return out


# ── parsing scoring tool output ──────────────────────────────────────────────────────────────


def _read_table(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _first_column(record: Dict[str, str], *names: str) -> Optional[str]:
    """The first of *names* present in *record*, case- and space-insensitively.

    Used because the tools spell columns differently (CheckM1 ``Bin Id``, CheckM2 ``Name``).
    """
    normalised = {key.strip().lower(): key for key in record if key}
    for name in names:
        key = normalised.get(name.strip().lower())
        if key is not None:
            return record[key]
    return None


def bin_quality(path: str) -> List[Tuple[str, float, float]]:
    """``(bin, completeness, contamination)`` from a CheckM1 or CheckM2 report.

    The format is detected from the column names, not the filename.
    """
    out: List[Tuple[str, float, float]] = []
    for record in _read_table(path):
        name = _first_column(record, "Name", "Bin Id", "bin", "genome")
        completeness = _first_column(record, "Completeness", "Completeness_General", "completeness")
        contamination = _first_column(record, "Contamination", "contamination")
        if name is None or completeness is None or contamination is None:
            continue
        try:
            out.append((name, float(completeness), float(contamination)))
        except ValueError:
            # A repeated header or a bin CheckM could not assess. Skipped, not counted as zero.
            continue
    return out


def quality_counts(bins: Sequence[Tuple[str, float, float]]) -> Dict[str, int]:
    """How many bins fall in each :data:`QUALITY_TIERS` tier, plus the total.

    Tiers are cumulative: a near-complete bin is also counted as high and medium quality.
    """
    counts = {name: 0 for name, _c, _x in QUALITY_TIERS}
    for _name, completeness, contamination in bins:
        for tier, min_completeness, max_contamination in QUALITY_TIERS:
            if completeness >= min_completeness and contamination <= max_contamination:
                counts[tier] += 1
    counts["bins_total"] = len(bins)
    return counts


#: AMBER column headings mapped to this harness's metric names. Matched case-insensitively;
#: unrecognised headings are ignored.
_AMBER_METRICS = {
    "f1 score (bp)": "amber_f1_bp",
    "ari (bp)": "amber_ari_bp",
    "average completeness (bp)": "amber_completeness_bp",
    "average purity (bp)": "amber_purity_bp",
    "accuracy (bp)": "amber_accuracy_bp",
    "percentage of assigned bps": "amber_assigned_bp_fraction",
    # Column names in AMBER 2.x results.tsv (the headings above are from its HTML report and docs).
    # Recall is AMBER's average completeness and precision its average purity. The ``_per_bp``,
    # ``_cami1``, ``_sem`` and ``_var`` columns are omitted: variances are not scores, and CAMI 1's
    # F1 uses a different definition.
    "f1_score_bp": "amber_f1_bp",
    "adjusted_rand_index_bp": "amber_ari_bp",
    "recall_avg_bp": "amber_completeness_bp",
    "precision_avg_bp": "amber_purity_bp",
    "accuracy_bp": "amber_accuracy_bp",
    "percentage_of_assigned_bps": "amber_assigned_bp_fraction",
}


def amber_metrics_by_tool(path: str) -> Dict[str, Dict[str, float]]:
    """Metrics from an AMBER ``results.tsv``, grouped by the binning each row describes.

    AMBER writes all query binnings of one invocation into one table, tagged by label.
    :func:`amber_metrics` merges them into one dict, which is only correct for a single query.

    The label comes from whichever identifying column is present (``Tool``, ``tool`` or ``binning
    label``, depending on version). Rows without a label are grouped under the empty string.
    """
    out: Dict[str, Dict[str, float]] = {}
    for record in _read_table(path):
        label = _first_column(record, "Tool", "tool", "binning label", "Binning", "label") or ""
        bucket = out.setdefault(label.strip(), {})
        long_metric = _first_column(record, "metric")
        long_value = _first_column(record, "value")
        if long_metric is not None and long_value is not None:
            name = _AMBER_METRICS.get(long_metric.strip().lower())
            if name:
                try:
                    bucket[name] = float(long_value)
                except ValueError:
                    pass
            continue
        for heading, raw in record.items():
            if heading is None:
                continue
            name = _AMBER_METRICS.get(heading.strip().lower())
            if not name or raw in (None, ""):
                continue
            try:
                bucket[name] = float(raw)
            except ValueError:
                continue
    return {label: metrics for label, metrics in out.items() if metrics}


def amber_metrics(path: str) -> Dict[str, float]:
    """Metrics from an AMBER ``results.tsv``, keyed by this harness's metric names.

    Handles both the wide layout (one column per metric) and the long layout (``metric``/``value``
    columns). When several rows give the same metric, the last one wins; in AMBER's output that is
    the pooled row.
    """
    out: Dict[str, float] = {}
    for record in _read_table(path):
        long_metric = _first_column(record, "metric")
        long_value = _first_column(record, "value")
        if long_metric is not None and long_value is not None:
            name = _AMBER_METRICS.get(long_metric.strip().lower())
            if name:
                try:
                    out[name] = float(long_value)
                except ValueError:
                    pass
            continue
        for heading, raw in record.items():
            if heading is None:
                continue
            name = _AMBER_METRICS.get(heading.strip().lower())
            if not name or raw in (None, ""):
                continue
            try:
                out[name] = float(raw)
            except ValueError:
                continue
    return out


def gunc_pass_fraction(path: str) -> Optional[float]:
    """Fraction of bins GUNC passes, or None if the file says nothing about it.

    A fraction, so variants with different bin counts can be compared.
    """
    total = passed = 0
    for record in _read_table(path):
        verdict = _first_column(record, "pass.GUNC", "pass_GUNC", "pass.gunc")
        if verdict is None:
            continue
        total += 1
        # gunc.passed_value accepts both True/False and 1.0/0.0, which GUNC writes depending on
        # pandas typing. `nan` (too little reference representation to assess) is not a pass but
        # stays in the denominator, matching recorded rows; see gunc_unassessed_fraction.
        if _gunc.passed_value(verdict) is True:
            passed += 1
    if not total:
        return None
    return passed / total


def gunc_unassessed_fraction(path: str) -> Optional[float]:
    """Fraction of bins GUNC could not assess (``nan`` verdicts), or None if there are no verdicts.

    Kept separate from the pass fraction because an unassessed bin is a reference-coverage gap, not
    a chimera. On one 258-bin set it was 24 bins, shifting the pass fraction by about six points.
    """
    total = unassessed = 0
    for record in _read_table(path):
        verdict = _first_column(record, "pass.GUNC", "pass_GUNC", "pass.gunc")
        if verdict is None:
            continue
        total += 1
        if _gunc.passed_value(verdict) is None:
            unassessed += 1
    if not total:
        return None
    return unassessed / total


def assembly_stats(path: str) -> Dict[str, float]:
    """N50 and total length of a FASTA assembly.

    Computed here to avoid an extra dependency whose version could change the numbers between runs.
    Also returns the contig count.
    """
    lengths: List[int] = []
    current = 0
    opener: Any = open
    mode = "rt"
    if path.endswith(".gz"):
        import gzip

        opener = gzip.open
    with opener(path, mode) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    total = sum(lengths)
    n50 = 0
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            n50 = length
            break
    return {
        "assembly_total_bp": float(total),
        "assembly_n50": float(n50),
        "assembly_contigs": float(len(lengths)),
    }


def read_stats(path: str) -> Dict[str, float]:
    """Reads, bases and Q30 bases in a FASTQ.

    Counted from the FASTQ, not from trimmer reports, because fastp and Trim Galore report
    different things and this pipeline disables Trim Galore's summary.

    All three are counts so they can be summed across samples and mate files. Divide
    ``q30_bases_kept`` by ``bases_kept`` for a fraction.
    """
    reads = bases = q30 = 0
    opener: Any = open
    if path.endswith(".gz"):
        import gzip

        opener = gzip.open
    with opener(path, "rb") as fh:
        for index, line in enumerate(fh):
            phase = index % 4
            if phase == 1:
                reads += 1
                bases += len(line.strip())
            elif phase == 3:
                q30 += line.strip().translate(_Q30_BYTES).count(b"\x01")
    return {
        "reads_kept": float(reads),
        "bases_kept": float(bases),
        "q30_bases_kept": float(q30),
    }


# ── ranking variants ─────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Scored:
    """One variant's result on a metric, with what it cost to get there beside it.

    No quality/cost ratio is provided: two variants can share a ratio at very different costs, and
    whether the extra cost is acceptable is left to the reader.
    """

    comparison: str
    variant: str
    dataset: str
    metric: str
    value: float
    cost_metric: str
    cost: float
    is_default: bool = False
    #: Number of per-sample rows that went into ``value``.
    samples: int = 1

    def describe(self) -> str:
        return "%-28s %-16s %12.4g %12s%s" % (
            self.variant,
            self.metric,
            self.value,
            "unmeasured" if not self.cost else "%.4g" % self.cost,
            "  <- current default" if self.is_default else "",
        )


def value_against_cost(
    rows: Iterable[Row],
    metric: str,
    *,
    cost_metric: str = "cpu_hours",
    defaults: Sequence[str] = (),
) -> List[Scored]:
    """Pair each variant's *metric* against its *cost_metric*, best on the metric first.

    Rows are combined across samples, giving one result per variant. Counts and durations are
    summed. Metrics whose unit is a fraction or an index are averaged (unweighted), since summing
    them would favour variants measured on more samples. Cost is always summed. A figure that needs
    weighted pooling must be recorded already pooled, with an empty sample; :attr:`Scored.samples`
    shows how many rows were combined.

    A variant with no recorded cost is kept, with ``cost`` of zero, and shown as unmeasured.

    Rows for the same dataset, comparison, variant, sample and metric are collapsed to the last one,
    because the results TSV is appended to and re-scored comparisons repeat rows. ``source`` is not
    part of the key: a re-score that reads a different bin set records a new source, and keying on
    it once double-counted COMEBin's 178 high-quality MAGs as 356. Metrics that several files add to
    are suffixed per bin set before they reach here.
    """
    latest: Dict[Tuple[str, ...], Row] = {}
    for row in rows:
        if row.metric not in (metric, cost_metric):
            continue
        latest[(row.dataset, row.comparison, row.variant, row.sample, row.metric)] = row

    values: Dict[Tuple[str, str, str], float] = {}
    counts: Dict[Tuple[str, str, str], int] = {}
    costs: Dict[Tuple[str, str, str], float] = {}
    for row in latest.values():
        key = (row.dataset, row.comparison, row.variant)
        if row.metric == metric:
            values[key] = values.get(key, 0.0) + row.value
            counts[key] = counts.get(key, 0) + 1
        else:
            costs[key] = costs.get(key, 0.0) + row.value

    averaged = UNITS.get(metric, "") in _AVERAGED_UNITS
    out = [
        Scored(
            dataset=dataset,
            comparison=comparison,
            variant=variant,
            metric=metric,
            value=(value / counts[(dataset, comparison, variant)]) if averaged else value,
            cost_metric=cost_metric,
            cost=costs.get((dataset, comparison, variant), 0.0),
            is_default=variant in defaults,
            samples=counts[(dataset, comparison, variant)],
        )
        for (dataset, comparison, variant), value in values.items()
    ]
    # Ties on the metric are broken by lower cost.
    out.sort(key=lambda r: (-r.value, r.cost or float("inf")))
    return out


def format_scores(scored: Sequence[Scored]) -> str:
    """The table a reader uses to pick a default."""
    if not scored:
        return "No variant recorded a result on this metric, so there is nothing to compare."
    head = scored[0]
    lines = [
        "",
        "  %s on %s: %s, with %s beside it"
        % (head.comparison, head.dataset, head.metric, head.cost_metric),
        "  %-28s %-16s %12s %12s" % ("variant", "metric", "value", head.cost_metric),
        "  " + "-" * 74,
    ]
    lines.extend("  " + r.describe() for r in scored)
    lines.append("")
    lines.append(
        "  Ranked on %s. %s is shown for reference and is not used as a denominator:\n"
        "  a quality-per-cost ratio favours the cheapest variant. A ranking on that ratio once\n"
        "  chose a default that recovered fewer genomes than the tool it replaced."
        % (head.metric, head.cost_metric)
    )
    # State that the value is an unweighted mean over samples, so a reader can recompute it.
    if UNITS.get(head.metric, "") in _AVERAGED_UNITS and max(r.samples for r in scored) > 1:
        lines.append(
            "  %s is a %s, so it is the unweighted mean over each variant's samples: %s.\n"
            "  The cost column is a total."
            % (
                head.metric,
                UNITS[head.metric],
                ", ".join("%s over %d" % (r.variant, r.samples) for r in scored),
            )
        )
    lines.append("")
    return "\n".join(lines)
