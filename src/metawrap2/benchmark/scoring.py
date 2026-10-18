"""Turning tool output into one tidy table, and into the ratio the defaults are decided on.

Every number the benchmark produces ends up as one row of ``dataset x comparison x variant x
sample x metric``. Long format rather than a wide table with a column per metric, because the set
of metrics grows every time a scorer is added, and a wide table means rewriting every historical
file when it does. Long rows append.

Cost is stored as metrics too - ``wall_seconds``, ``cpu_hours``, ``peak_rss_bytes`` - rather than
as special columns beside the "real" results, so that cost is joinable against accuracy in the same
way one accuracy metric joins against another. :func:`value_against_cost` is that join, and it is
what settles which variant becomes a default.

It ranks on the quality metric and prints cost beside it. It used to rank on the *ratio* of the two,
and that is the one piece of machinery in this project that produced a wrong answer rather than a
missing one: the ratio is dominated by whichever arm is cheapest, and following it picked a default
that recovers fewer MAGs than metaWRAP 1 does. Quality decides; cost is for the reader to weigh, and
where it genuinely has to drive the choice it becomes a tier a user picks (``binning --preset``).

The scorers here parse output rather than reimplement scoring. CheckM2, GUNC and AMBER are the
tools that know how to judge a bin; this module's job is to get their answers into one file with
consistent names, so that two runs a year apart can be diffed.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from .. import gunc as _gunc

#: The TSV's columns, in order. Frozen on purpose: results tables accumulate over
#: time, and a column reorder would make every old file unreadable by the current code.
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
#: Two of the three are MIMAG's (Bowers et al. 2017) and one is ours, which is worth stating because
#: the file used to call all three "MIMAG-style":
#:
#: * ``bins_medium_quality`` (>=50 / <10) is MIMAG medium exactly.
#: * ``bins_near_complete`` (>=90 / <5) is MIMAG high-quality *minus* its rRNA and tRNA requirement -
#:   23S, 16S and 5S present plus at least 18 tRNAs - which is the half that is hard to satisfy and
#:   the reason this tier is not called high-quality here.
#: * ``bins_high_quality`` (>=70 / <10) is **not a MIMAG tier at all**. It is the working threshold
#:   this project reports on, chosen because it is where bin_refinement's defaults sit and where the
#:   MAG-recovery literature usually counts, and it is the headline number in docs/benchmarks.md.
#:
#: The name stays despite being a borrowed label: it is the metric key in every recorded results
#: row, so renaming it either invalidates those rows or needs a migration, and the
#: label was the inaccurate part rather than the number.
QUALITY_TIERS = (
    ("bins_near_complete", 90.0, 5.0),
    ("bins_high_quality", 70.0, 10.0),
    ("bins_medium_quality", 50.0, 10.0),
)

#: Units, so a reader of the TSV never has to guess whether a duration is seconds or hours. A
#: metric absent from here is recorded with an empty unit rather than a wrong one.
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

#: The units whose metrics must not be summed across samples. A fraction and an index are already
#: normalised - three samples at 0.9 purity are 0.9 purity, not 2.7 - so a variant's per-sample rows
#: are averaged instead. Read off UNITS rather than kept as a second list, so a new fraction metric
#: gets the right pooling by declaring its unit and nothing else.
_AVERAGED_UNITS = frozenset(("fraction", "index"))

#: Byte lookup marking Phred+33 quality characters of 30 or better (ASCII 63, ``?``, and up). Used
#: with :meth:`bytes.translate` so a Q30 count is one C-level pass per line rather than a Python
#: loop over every base - the difference between seconds and minutes on a real FASTQ.
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
            # Written with repr-ish precision rather than rounded: rounding at the storage layer
            # loses information that a later ratio needs.
            "value": ("%g" % self.value) if isinstance(self.value, float) else str(self.value),
            "unit": self.unit or UNITS.get(self.metric, ""),
            "source": self.source,
        }


def write_tsv(rows: Iterable[Row], path: str, *, append: bool = False) -> int:
    """Write *rows* to *path*. Returns the number written.

    Appending writes no second header, so a results file accumulated over several runs stays one
    readable table rather than a concatenation of files.
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
    """Read a results TSV back. Rows whose value is not a number are skipped, not guessed at."""
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


# ── parsing what the scoring tools produce ───────────────────────────────────────────────────


def _read_table(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _first_column(record: Dict[str, str], *names: str) -> Optional[str]:
    """The first of *names* present in *record*, case- and space-insensitively.

    The bin-quality tools disagree on spelling - CheckM1 writes ``Bin Id``, CheckM2 writes
    ``Name`` - and a benchmark that only understood one of them would quietly score nothing.
    """
    normalised = {key.strip().lower(): key for key in record if key}
    for name in names:
        key = normalised.get(name.strip().lower())
        if key is not None:
            return record[key]
    return None


def bin_quality(path: str) -> List[Tuple[str, float, float]]:
    """``(bin, completeness, contamination)`` from a CheckM1 or CheckM2 report.

    Both are tab tables with a completeness and a contamination column; which tool wrote it is
    decided by the columns present rather than by the filename, because the filename is whatever
    the caller passed.
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
            # A header repeated mid-file, or a bin CheckM could not judge. Skipping it is right;
            # counting it as zero would make a failed bin look like a bad one.
            continue
    return out


def quality_counts(bins: Sequence[Tuple[str, float, float]]) -> Dict[str, int]:
    """How many bins fall in each :data:`QUALITY_TIERS` tier, plus the total.

    Tiers are not exclusive - a near-complete bin is also high-quality - because that is how the
    thresholds are used in practice: "how many bins clear 90/5" is the question, not "how many
    clear 90/5 and not 95/2".
    """
    counts = {name: 0 for name, _c, _x in QUALITY_TIERS}
    for _name, completeness, contamination in bins:
        for tier, min_completeness, max_contamination in QUALITY_TIERS:
            if completeness >= min_completeness and contamination <= max_contamination:
                counts[tier] += 1
    counts["bins_total"] = len(bins)
    return counts


#: AMBER column headings mapped to the metric names this harness records. AMBER's own headings are
#: long and change punctuation between versions, so matching is done case-insensitively on these
#: exact strings and anything unrecognised is ignored rather than guessed at.
_AMBER_METRICS = {
    "f1 score (bp)": "amber_f1_bp",
    "ari (bp)": "amber_ari_bp",
    "average completeness (bp)": "amber_completeness_bp",
    "average purity (bp)": "amber_purity_bp",
    "accuracy (bp)": "amber_accuracy_bp",
    "percentage of assigned bps": "amber_assigned_bp_fraction",
    # What AMBER 2.x actually writes into results.tsv - the prose headings above are the ones on its
    # HTML page and in its docs. Completeness and purity appear under their information-retrieval
    # names: AMBER's "average completeness" is recall and its "average purity" is precision, the same
    # numbers under the vocabulary each field prefers. The ``_per_bp``, ``_cami1``, ``_sem`` and
    # ``_var`` columns are deliberately absent: a variance is not a score, and CAMI 1's F1 is a
    # different definition that must not be mixed into one column with this one.
    "f1_score_bp": "amber_f1_bp",
    "adjusted_rand_index_bp": "amber_ari_bp",
    "recall_avg_bp": "amber_completeness_bp",
    "precision_avg_bp": "amber_purity_bp",
    "accuracy_bp": "amber_accuracy_bp",
    "percentage_of_assigned_bps": "amber_assigned_bp_fraction",
}


def amber_metrics_by_tool(path: str) -> Dict[str, Dict[str, float]]:
    """Metrics from an AMBER ``results.tsv``, grouped by the binning each row describes.

    AMBER takes several query binnings in one invocation and writes them all into one table, tagged
    with the label they were submitted under. :func:`amber_metrics` folds that table into a single
    dict, which is right when there was one query and silently wrong when there were three: the last
    tool's numbers would be reported as though they were the run's.

    The label is taken from whichever of AMBER's identifying columns is present - it has been
    ``Tool``, ``tool`` and ``binning label`` across versions - and rows with no label at all are
    grouped under the empty string, which is how a single unlabelled run reads.
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

    Handles both shapes AMBER has shipped: a wide table with one column per metric, and a long one
    with ``metric``/``value`` columns. Where several rows describe the same metric - AMBER reports
    per-sample as well as pooled figures - the last one wins, which is the pooled row in the files
    AMBER writes.
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

    A fraction rather than a count so it is comparable between variants that produced different
    numbers of bins - which is the whole point of comparing binners.
    """
    total = passed = 0
    for record in _read_table(path):
        verdict = _first_column(record, "pass.GUNC", "pass_GUNC", "pass.gunc")
        if verdict is None:
            continue
        total += 1
        # Through gunc.passed_value because GUNC writes this column as True/False or as 1.0/0.0
        # depending on how pandas typed it, and accepting only the words reported a 166-bin set as
        # 0% non-chimeric. `nan` - a bin with too little reference representation to assess - is
        # neither a pass nor a fail, but it stays in the denominator, because that is what this
        # metric's committed rows mean; gunc_unassessed_fraction says how much of it that is.
        if _gunc.passed_value(verdict) is True:
            passed += 1
    if not total:
        return None
    return passed / total


def gunc_unassessed_fraction(path: str) -> Optional[float]:
    """Fraction of bins GUNC could not assess - the part of the denominator above that is not a
    verdict.

    Reported separately rather than folded in, because "chimeric" and "too little of this lineage
    in the reference to say" are different findings and only the first is about the bin set. On one
    real 258-bin set it is 24 bins, enough to move a pass fraction by six points.
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

    Computed here rather than shelled out to a stats tool because it is twenty lines, it removes a
    dependency from the benchmark env, and an assembly metric that disagrees with itself between
    runs because a helper changed version is worse than no metric.
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

    Counted off the FASTQ itself rather than read out of a trimmer's own report, because the
    trimmers do not agree on what to report: fastp writes a rich JSON, Trim Galore writes a text
    summary this pipeline switches off, and a comparison scored from two different tools'
    self-assessments is not a comparison. Reading the surviving reads asks both the same question.

    All three are counts, so they sum across samples and across mate files the way
    :func:`time_to_performance` expects. Q30 is recorded as a count of bases rather than a fraction
    for exactly that reason - a reader divides it by ``bases_kept`` themselves.
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


# ── the number the defaults are decided on ───────────────────────────────────────────────────


@dataclass(frozen=True)
class Scored:
    """One variant's result on a metric, with what it cost to get there beside it.

    Deliberately not a ratio. Two variants can share a ratio with one being twice as accurate for
    twice the cost, and which of those should be the default depends on whether the machine has the
    hours - a judgement a reader makes, not one this dataclass makes for them. There is no ``ratio``
    property to reach for, because when there was one it was what everything got ranked by.
    """

    comparison: str
    variant: str
    dataset: str
    metric: str
    value: float
    cost_metric: str
    cost: float
    is_default: bool = False
    #: How many per-sample rows went into ``value``. Carried so a reader can tell a figure measured
    #: once from one averaged over ten samples, which is the difference between a number and a mean.
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

    Rows are combined across samples for both sides, so a ten-sample dataset yields one ratio per
    variant rather than ten. Counts and durations are summed, because bins recovered and CPU-hours
    spent are both totals. Metrics whose unit is a fraction or an index are averaged instead: those
    are already normalised, and summing ten samples at 0.9 purity into 9.0 would not just inflate
    one row, it would rank the variant that was measured on the most samples first. The cost side is
    always summed, since ten samples really do cost ten samples' CPU-hours.

    The average is unweighted - a sample holding five bins counts as much as one holding two
    hundred. That is a pooling choice this function cannot make better, because the weights live in
    a different metric; a figure that needs bin-weighted pooling has to be recorded already pooled,
    with an empty sample, and :attr:`Ratio.samples` is carried so a reader can see which they have.

    A variant with no recorded cost is returned with ``cost`` of zero and reported as unmeasured
    rather than dropped, because a missing cost is a gap in the measurement worth seeing in the
    table, not a variant that was free. It still ranks on its quality, which is the point of ranking
    on quality: an arm whose cost was lost is still a comparable result.

    Rows that describe the same measurement - same dataset, comparison, variant, sample, metric and
    source - are collapsed to the last one rather than summed. The results TSV is appended to, so a
    comparison scored twice, or re-run after a fix, holds the same measurement more than once;
    summing those would double a variant's bins or its CPU-hours depending on which half was
    repeated. The latest value is the answer, and the older rows stay in the file as history.

    ``source`` is deliberately *not* part of that key, and leaving it in was a real defect rather than
    a harmless one. A re-score that changes which bin set is a variant's answer - which is exactly what
    fixing :func:`metawrap2.benchmark.assess.bin_dirs` did - records the new number against a new
    source file, so the stale row keyed differently and the two were added: COMEBin's 178 high-quality
    MAGs read as 356. Same dataset, comparison, variant, sample and metric is the same measurement, no
    matter which file it was read out of, and a metric that two files contribute to additively is
    suffixed per bin set long before it reaches here.
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
    # Ties on the metric are broken by the cheaper arm, which is the only place cost decides
    # anything here - and only between arms that scored identically.
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
        "  Ranked on %s, with %s reported beside it and never as a denominator: a\n"
        "  quality-per-cost ratio is dominated by whichever arm is cheapest, and ranking on one\n"
        "  picked a default here that recovered fewer genomes than the tool this replaces. Where\n"
        "  cost has to drive a choice it belongs in a tier the user picks, not in the ranking."
        % (head.metric, head.cost_metric)
    )
    # Said out loud rather than left to the units column, because an unweighted mean over samples of
    # different sizes is a figure a reader may want to recompute by hand, and they can only know to
    # want that if the table admits it took one.
    if UNITS.get(head.metric, "") in _AVERAGED_UNITS and max(r.samples for r in scored) > 1:
        lines.append(
            "  %s is a %s, so it is averaged over each variant's samples rather than summed,\n"
            "  unweighted: %s. The cost column is still a total."
            % (
                head.metric,
                UNITS[head.metric],
                ", ".join("%s over %d" % (r.variant, r.samples) for r in scored),
            )
        )
    lines.append("")
    return "\n".join(lines)
