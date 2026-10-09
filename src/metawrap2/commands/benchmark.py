"""`metawrap2 benchmark list|fetch|run|score` - measure the choices instead of arguing about them.

MetaWrap2's defaults are supposed to be chosen on measured quality, with cost reported beside it
rather than divided into it. This command is what makes
that claim checkable by someone other than the person who made it: the datasets are declared as
data, the comparisons are declared as data, and the results land in a committed TSV so a changed
default arrives in review next to the numbers that justify it.

Four verbs, in the order they are used:

- ``list``   what datasets and comparisons exist, and what a fetch would weigh.
- ``fetch``  download one dataset, checksum-verified, inside the ``[limits]`` disk ceiling.
- ``run``    run one comparison's variants over one dataset under one budget, recording cost.
- ``score``  read the tools' output into the tidy TSV, and rank the variants with their cost beside.

``fetch`` never runs implicitly. ``run`` on a dataset that is not present says how to fetch it and
stops; it does not helpfully start a 95 GB download because someone typed a dataset name.
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Sequence, Tuple

from .. import command as _command
from .. import gunc as _gunc
from .. import limits as _limits
from ..benchmark import abundance as _abundance
from ..benchmark import assess as _assess
from ..benchmark import fetch as _fetch
from ..benchmark import registry as _registry
from ..benchmark import runner as _runner
from ..benchmark import scoring as _scoring
from ..cliformat import make_parser
from ..config import load_settings
from ..constants import is_refined_bins_dirname


def default_results(dataset_name: str) -> str:
    """Where a dataset's results land by default. Committed, so a regression is a diff.

    Anchored to the registry - and therefore to the checkout - rather than to the current directory.
    This was ``os.path.join("benchmarks", "results")``, a relative path, so scoring a comparison from
    anywhere but the repository root silently wrote a *second* results file under wherever that happened
    to be, reported success, and left the committed one untouched. The numbers were real and had to be
    merged back in by hand; the run that produced them looked completely normal.
    """
    return os.path.join(_registry.registry_dir(), "results", "%s.tsv" % dataset_name)


def _dry_run() -> bool:
    """Is the global ``--dry-run`` in effect?

    Read from the shared runner rather than from an option of this command's own. ``--dry-run`` is
    a global flag that ``metawrap2``'s dispatcher strips before the subcommand's parser ever sees
    it, so a local ``--dry-run`` here would be unreachable - and, worse, would default to False
    and run for real while the user believed they had asked for a plan.
    """
    return bool(_command.runner.dry_run)


def _budget(args) -> _limits.Budget:
    """The budget for this run: the config's ``[limits]``, with the global ``--max-*`` flags applied.

    Resolved once here and handed down, rather than re-resolved per variant, so every variant of a
    comparison is measured against the same grant. Two variants given different grants are not a
    comparison.
    """
    settings = load_settings(getattr(args, "config", None))
    return _limits.resolve(settings.limits)


def _list(args) -> int:
    datasets = _registry.load_datasets()
    comparisons = _registry.load_configurations()
    print("\nDatasets (%s)" % _registry.datasets_path())
    for dataset in datasets.values():
        print("  %s" % dataset.describe())
        if args.verbose:
            if dataset.citation:
                print("      cite: %s" % dataset.citation)
            if dataset.doi:
                print("      doi:  %s" % dataset.doi)
            for gold in dataset.gold_standard:
                print("      gold: %-10s %s" % (gold.role, gold.path))

    print("\nComparisons (%s)" % _registry.configurations_path())
    for comparison in comparisons.values():
        incumbent = comparison.incumbent
        print("  %-24s %-16s %s" % (comparison.name, comparison.module, comparison.question or ""))
        for variant in comparison.variants:
            mark = " <- current default" if variant.default else ""
            print("      %-30s %s%s" % (variant.name, " ".join(variant.args), mark))
        if incumbent is None:
            print("      (no variant marked default; nothing to compare against)")
    print(
        "\nFetch one with:  metawrap2 benchmark fetch <dataset> --dest <dir>\n"
        "Then:            metawrap2 benchmark run <dataset> <comparison> --data <dir> -o <dir>\n"
    )
    return 0


def _fetch_cmd(args) -> int:
    dataset = _registry.get_dataset(args.dataset)
    dest = args.dest or os.path.join("benchmarks", "data", dataset.name)
    plan = _fetch.plan(dataset, dest, hash_it=args.verify)

    print()
    print(plan.describe())
    if plan.complete:
        return 0
    if dataset.citation:
        print("\n  Cite: %s" % dataset.citation)
    if dataset.license:
        print("  License: %s" % dataset.license)

    ok, message = _fetch.disk_verdict(plan, _budget(args))
    print("\nDisk: %s" % message)

    dry_run = _dry_run()
    if not ok:
        # Answered before the confirmation is asked for: there is no point confirming a download
        # that cannot land.
        print(
            "\nRefusing to start: this would cross the [limits] disk ceiling.\n"
            "Choose a --dest on another filesystem, or raise max_disk_fraction if you know the\n"
            "rest of the machine can spare it."
        )
        return 1
    # A benchmark dataset is ~95 GB from a public host. The size is printed and confirmed rather
    # than assumed acceptable, because "I typed the dataset name" is not the same as "I meant to
    # spend an afternoon and a hundred gigabytes".
    if not args.yes and not dry_run:
        print(
            "\nThis will download %s to %s.\nRe-run with --yes to start it, or --dry-run to see "
            "the plan without downloading." % (_limits.human(plan.bytes_needed), dest)
        )
        return 1

    return _fetch.fetch(
        dataset,
        dest,
        budget=_budget(args),
        keep_archives=args.keep_archives,
        dry_run=dry_run,
    )


def _run(args) -> int:
    dataset = _registry.get_dataset(args.dataset)
    comparison = _registry.get_comparison(args.comparison)
    data_dir = args.data or os.path.join("benchmarks", "data", dataset.name)
    work_dir = args.output or os.path.join("benchmarks", "work", dataset.name)
    dry_run = _dry_run()

    try:
        if not dry_run:
            _runner.require_env_manager()
        rows = _runner.run_comparison(
            dataset,
            comparison,
            data_dir,
            work_dir,
            variants=args.variant or (),
            samples=args.sample or (),
            max_samples=args.max_samples or 0,
            assembly=args.assembly,
            bins_from=args.bins_from,
            budget=_budget(args),
            dry_run=dry_run,
        )
    except _runner.RunnerError as exc:
        print("\n%s\n" % exc)
        return 1

    if dry_run:
        print("\n--dry-run: printed the commands without running them.\n")
        return 0
    if not rows:
        print("\nNo variant completed, so there is nothing to record.\n")
        return 1

    path = args.results or default_results(dataset.name)
    written = _scoring.write_tsv(rows, path, append=True)
    print("Recorded %d cost row(s) in %s" % (written, path))
    print(
        "\nScore the output when the runs are judged:\n"
        "  metawrap2 benchmark score %s %s --work %s\n" % (dataset.name, comparison.name, work_dir)
    )
    return 0


def _record_field(work_dir: str, comparison: str, field: str):
    """One field out of a run's record, or None when there is no readable record.

    Scoring has to work on a run someone else launched, possibly before a field existed, so a
    missing record or a missing key is an absence to fall back from rather than an error.
    """
    path = os.path.join(work_dir, "%s.run.json" % comparison)
    try:
        with open(path) as fh:
            return json.load(fh).get(field)
    except (OSError, ValueError):
        return None


def _say_sample_scope(work_dir: str, comparison: str) -> None:
    """Print which samples the run covered, if its record is still there.

    A cohort module's rows carry no sample label, so the table cannot show on its own whether a
    number came from three samples or ten - and `--max-samples` makes that a real possibility
    rather than a hypothetical one. Missing record, no claim: an older run predates this.
    """
    path = os.path.join(work_dir, "%s.run.json" % comparison)
    try:
        with open(path) as fh:
            samples = json.load(fh).get("samples") or []
    except (OSError, ValueError):
        return
    if samples:
        print("\n  Samples measured: %d  (%s)" % (len(samples), ", ".join(samples)))


def _score(args) -> int:
    dataset = _registry.get_dataset(args.dataset)
    comparison = _registry.get_comparison(args.comparison)
    work_dir = args.work or os.path.join("benchmarks", "work", dataset.name)
    path = args.results or default_results(dataset.name)

    metrics_wanted = args.metric or list(comparison.metrics)
    budget = _budget(args)
    # The run record names the data the run actually used, so scoring finds the same gold standard
    # without being told twice. --data overrides it, and covers records written before it was kept.
    data_dir = (
        args.data
        or _record_field(work_dir, comparison.name, "data_dir")
        or os.path.join("benchmarks", "data", dataset.name)
    )
    # The contigs the binners were given, for the sequence lengths AMBER weights its metrics by.
    # Taken from the record rather than from the dataset, because --assembly means the dataset's
    # gold assembly is not always what ran; it is only the fallback for records written without it.
    assembly = args.assembly or _record_field(work_dir, comparison.name, "assembly")
    if not assembly or not os.path.isfile(assembly):
        assembly = _runner.gold_path(dataset, "assembly", data_dir)

    # Built once, from the dataset rather than from any run, and only when something asks for it:
    # it reads every sample's read mapping, which is minutes of gunzip nobody should pay for
    # scoring a binner. All variants then score against the same object, which is the point - three
    # quantifiers measured against two truths would compare the truths.
    abundance_truth = None
    if any(m.startswith("abundance_") for m in metrics_wanted):
        abundance_truth = _abundance_truth(
            dataset, comparison, work_dir, data_dir=data_dir, assembly=assembly
        )

    # Which variants exited non-zero, as the runner recorded them. The runner already discards a
    # failed variant's cost rows; scoring read the record for data_dir and assembly but not for
    # this, so a variant that crashed after its first binner wrote a directory had those partial
    # bins scored, ranked, and written to the committed TSV - indistinguishable in the table from
    # an arm whose cost measurement was merely lost.
    failed = {
        name: record.get("returncode")
        for name, record in (_record_field(work_dir, comparison.name, "variants") or {}).items()
        if record.get("returncode")
    }

    rows: List[_scoring.Row] = []
    unranked: List[tuple] = []
    skipped: List[tuple] = []
    for variant in comparison.variants:
        root = os.path.join(work_dir, comparison.name, variant.name)
        if not os.path.isdir(root):
            continue
        if variant.name in failed and not getattr(args, "score_failed", False):
            skipped.append((variant.name, failed[variant.name]))
            continue
        bin_sets = _assess_variant(
            root,
            dataset=dataset,
            metrics=metrics_wanted,
            budget=budget,
            variant=variant.name,
            data_dir=data_dir,
            assembly=assembly,
            reassess=bool(getattr(args, "reassess", False)),
            checkm2=bool(getattr(args, "checkm2", False)),
        )
        found = _collect_scores(
            dataset.name,
            comparison.name,
            variant.name,
            root,
            per_sample=comparison.module in _runner.PER_SAMPLE_MODULES,
            ranked_bin_set=_ranked_bin_set(bin_sets),
            abundance_truth=abundance_truth,
        )
        if found:
            print("  %-30s %d metric(s)" % (variant.name, len(found)))
            rows.extend(found)
        else:
            print(
                "  %-30s nothing to score (no CheckM/GUNC/AMBER/read output under %s)"
                % (variant.name, root)
            )
        if len(bin_sets) > 1 and _ranked_bin_set(bin_sets) is None:
            unranked.append((variant.name, [os.path.basename(p) for p in bin_sets]))

    if skipped:
        print(
            "\n  Not scored - %d variant(s) whose run failed, so their output is partial:"
            % len(skipped)
        )
        for name, rc in skipped:
            print("    %-28s exited %s" % (name, rc))
        print("    Re-run them, or pass --score-failed to score what is there anyway.")

    if unranked:
        print(
            "\n  Not ranked - %d variant(s) left several raw bin sets and no consolidated one:"
            % len(unranked)
        )
        for name, labels in unranked:
            print("    %-28s %s" % (name, ", ".join(labels)))
        print(
            "  Their numbers are recorded per bin set (metric__<set>) but not summed into a\n"
            "  ranking: the sets overlap, so adding them counts the same genome more than once.\n"
            "  Consolidate them first - that is what bin_refinement is for - and score again."
        )

    if rows:
        written = _scoring.write_tsv(rows, path, append=True)
        print("\nRecorded %d accuracy row(s) in %s" % (written, path))

    if not os.path.isfile(path):
        print("\nNo results file at %s yet.\n" % path)
        return 1

    _say_sample_scope(work_dir, comparison.name)
    all_rows = _scoring.read_tsv(path)
    defaults = [v.name for v in comparison.variants if v.default]
    metrics = metrics_wanted
    printed = False
    for metric in metrics:
        scored = [
            r
            for r in _scoring.value_against_cost(
                all_rows, metric, cost_metric=args.cost_metric, defaults=defaults
            )
            if r.comparison == comparison.name and r.dataset == dataset.name
        ]
        if scored:
            print(_scoring.format_scores(scored))
            printed = True
    if not printed:
        print(
            "\nNothing to rank yet: %s holds no rows pairing %s with %s for this comparison.\n"
            % (path, ", ".join(metrics) or "<no metrics declared>", args.cost_metric)
        )
    return 0


def _abundance_truth(
    dataset,
    comparison,
    work_dir: str,
    *,
    data_dir: str,
    assembly: str,
) -> Optional[Dict[str, Dict[str, float]]]:
    """Build the expected fraction per bin per sample, or explain why it cannot be built.

    Returning None rather than raising: a missing abundance truth means the abundance metrics go
    unscored, and everything else in the comparison still scores. Refusing to print a CheckM number
    because a read mapping is absent would be the wrong trade.

    The bins come from the run record, not from the variants' own output directories. Every variant
    quantified the *same* bins - that is what makes the comparison a comparison - and the bins are
    the thing the truth is about, so reading them from one place keeps the three tables on one set of
    row labels.
    """
    bins_from = _record_field(work_dir, comparison.name, "bins_from")
    if not bins_from or not os.path.isdir(bins_from):
        print(
            "  Abundance metrics not scored: the run record for %s names no --bins-from directory\n"
            "  that still exists, and a quantifier's truth is a statement about specific bins."
            % comparison.name
        )
        return None
    contig_mapping = _runner.gold_path(dataset, "binning", data_dir)
    if not (contig_mapping and os.path.isfile(contig_mapping)) or not (
        assembly and os.path.isfile(assembly)
    ):
        print(
            "  Abundance metrics not scored: need both the gold assembly and its contig-to-genome\n"
            "  mapping to work out how much of each genome landed in each bin."
        )
        return None
    try:
        read_mappings = _runner.abundance_truth_files(dataset, data_dir)
        print("  Building the abundance truth from %d read mapping(s)..." % len(read_mappings))
        return _abundance.build_truth(
            _runner.one_bin_set(bins_from), assembly, contig_mapping, read_mappings
        )
    except (_runner.RunnerError, _abundance.AbundanceError, OSError) as exc:
        print("  Abundance metrics not scored: %s" % exc)
        return None


def _assess_variant(
    root: str,
    *,
    dataset,
    metrics: Sequence[str],
    budget: _limits.Budget,
    variant: str,
    data_dir: str,
    assembly: Optional[str] = None,
    reassess: bool = False,
    checkm2: bool = False,
) -> List[str]:
    """Run whichever assessors *metrics* needs over the bin sets under *root*. Returns those sets.

    Assessment happens here, at score time, rather than inside the measured run - see
    :mod:`metawrap2.benchmark.assess` for why the cost column depends on that choice.

    Already-assessed sets are skipped, because scoring is the cheap half of a benchmark and is meant
    to be re-runnable: re-reading a table to add a metric should not restart a CheckM run that took
    an hour. ``--reassess`` forces the work again for when the bins themselves have changed.
    """
    bin_sets = _assess.bin_dirs(root)
    if not bin_sets:
        return bin_sets
    if not (
        _assess.needs_amber(metrics) or _assess.needs_checkm(metrics) or _assess.needs_gunc(metrics)
    ):
        return bin_sets

    gold = _runner.gold_path(dataset, "binning", data_dir)
    if gold is not None and not os.path.isfile(gold):
        gold = None
    outstanding = (
        bin_sets if reassess else [p for p in bin_sets if not _already_assessed(p, metrics)]
    )
    # AMBER is asked about separately because its output is one table for the whole variant rather
    # than one file per bin set, so no per-set check can speak for it. It used to be covered by the
    # per-set check standing in for every assessor, which meant a variant whose sets already had
    # CheckM tables never got AMBER run on them at all.
    amber_outstanding = _assess.needs_amber(metrics) and not os.path.isfile(
        os.path.join(root, "amber", _assess.AMBER_RESULTS)
    )
    if not outstanding and not amber_outstanding:
        return bin_sets

    print(
        "  %-30s assessing %d bin set(s): %s"
        % (variant, len(outstanding), ", ".join(os.path.basename(p) for p in outstanding))
    )
    _assess.assess_bins(
        root,
        metrics=metrics,
        gold_binning=gold,
        assembly=assembly,
        threads=budget.threads,
        memory_gb=int(budget.memory_gb),
        checkm2=checkm2,
        reassess=reassess,
        say=print,
    )
    return bin_sets


def _already_assessed(bin_set: str, metrics: Sequence[str]) -> bool:
    """Has *bin_set* already got everything *metrics* asks of a per-set assessor?

    Per-set is the qualifier that matters: CheckM and GUNC write one file per bin set, so their
    output is evidence about this set, while AMBER writes one table for the variant and is left to
    :func:`metawrap2.benchmark.assess.assess_bins` to decide about. A set that has its CheckM table
    but not its GUNC report is not assessed, which is the case this function existed in a simpler
    form for and got wrong: GUNC was added as a second per-set assessor and every set that had been
    through CheckM would have been skipped before GUNC ever ran.
    """
    present = []
    if _assess.needs_checkm(metrics):
        present.append(os.path.isfile(bin_set + _STATS_SUFFIX))
    if _assess.needs_gunc(metrics):
        present.append(_gunc.find_report(bin_set + _GUNC_SUFFIX) is not None)
    # A set nothing per-set was asked of counts as done, which leaves the AMBER-only case to the
    # caller rather than answering it wrongly here.
    return all(present)


#: The suffix MetaWrap2 gives its own normalised bin-quality table, kept as a local name so the
#: filename rule below reads as one rule rather than a string literal repeated three times.
_STATS_SUFFIX = ".stats.tsv"

#: The bin set a consolidating module writes is recognised by its shape - metawrap_70_10_bins, named
#: for the thresholds it was cut at - not by one literal name. When it is present it is the variant's
#: answer, because it is the set a user would actually carry forward. See
#: :func:`metawrap2.constants.is_refined_bins_dirname`.


def _ranked_bin_set(bin_sets: Sequence[str]) -> Optional[str]:
    """Which of *bin_sets* holds the numbers that are the variant's own, if any.

    A variant that ran three binners leaves three bin sets, and the three overlap: the same genome
    is usually recovered by more than one of them. The ranker sums a metric across everything it
    finds for a variant, so recording all three as the variant's numbers would count that genome
    three times and report an F1 of 2.4 - a number with no meaning that nonetheless sorts first.

    So exactly one set can carry the variant's ranked metrics:

    - the consolidated set if the variant has one, since that is what a user keeps;
    - otherwise the only set, when there is only one;
    - otherwise none, and the caller says so rather than picking a winner arbitrarily. Three raw
      bin sets are not yet an answer to "which binner set is better" - they become one after
      consolidation, and until then the honest report is that the variant cannot be ranked.
    """
    labels = [os.path.basename(path) for path in bin_sets]
    consolidated = [label for label in labels if is_refined_bins_dirname(label)]
    if len(consolidated) == 1:
        return consolidated[0]
    if consolidated:
        # Several, which is a run resumed at different thresholds. Ranking an arbitrary one is the
        # thing this function exists to refuse, and it is not more defensible here than for raw sets.
        return None
    if len(labels) == 1:
        return labels[0]
    return None


def _for_bin_set(metric: str, label: str, ranked: Optional[str]) -> str:
    """*metric* as recorded for the bin set *label*: plain if it is the ranked set, suffixed if not.

    The suffix keeps a per-binner number in the TSV - it is genuinely interesting which binner found
    what - while making certain it is never summed into a ranking, because a ranking is over metric
    names and a suffixed name is a different metric.
    """
    if not label or (ranked is not None and label == ranked):
        return metric
    return "%s__%s" % (metric, label)


#: Where bin_refinement puts GUNC's report on the set it hands back (see its ``work_files/gunc``).
#: There is no bin set name in that path, and there does not need to be: the module runs GUNC once,
#: on the set it returns, so the report found there is the variant's own figure.
_GUNC_FINAL_DIR = os.path.join("work_files", "gunc")

#: What score time appends to a bin directory to hold its GUNC report, the way CheckM's goes to
#: ``<bin set>.checkm``. Here the bin set *is* in the path, which is what makes several sets per
#: variant tellable apart.
_GUNC_SUFFIX = ".gunc"


def _gunc_bin_set(path: str, root: str, sample: str) -> str:
    """Which bin set the GUNC report at *path* describes, as a label for :func:`_for_bin_set`.

    GUNC's filenames say which database and CSS level the report is for, never which bins, so unlike
    CheckM and AMBER there is no label to read out of the report itself. The directory holding it is
    the evidence, and there are two that mean something:

    * ``<bin set>.gunc/``, written at score time beside the bins it describes, gives ``<bin set>`` -
      the same label AMBER's rows carry, so the two tools' numbers line up on one set.
    * ``work_files/gunc/``, written by bin_refinement, gives the empty label, which
      :func:`_for_bin_set` reads as the variant's own figure. The module runs GUNC once and only on
      the set it hands back, so there is nothing else it could describe.

    Anything else keeps its path as the label, which is a number recorded without being ranked. That
    is the conservative end: a fraction is not summable, and two reports under one metric name would
    be added together by the ranker, so two sets that each passed every bin would rank at 2.0 -
    ahead of a variant that was perfect and only scored once.
    """
    directory = os.path.dirname(path)
    if directory.endswith(_GUNC_SUFFIX):
        return os.path.basename(directory)[: -len(_GUNC_SUFFIX)]
    relative = os.path.relpath(directory, root)
    if sample:
        relative = os.path.relpath(relative, sample)
    return "" if relative == _GUNC_FINAL_DIR else relative


def _has_sibling_stats(path: str) -> bool:
    """Is there a normalised ``.stats.tsv`` describing the same bins as the native report at *path*?

    Both CheckMs write their own report under ``<bin set>.checkm/`` and MetaWrap2 writes the
    normalised table next to it as ``<bin set>.stats.tsv``. Reading both would record every count
    twice under two sources, and the ranker sums sources.
    """
    checkm_dir = path
    while checkm_dir and not checkm_dir.endswith(".checkm"):
        parent = os.path.dirname(checkm_dir)
        if parent == checkm_dir:
            return False
        checkm_dir = parent
    return os.path.isfile(checkm_dir[: -len(".checkm")] + _STATS_SUFFIX)


def _collect_scores(
    dataset: str,
    comparison: str,
    variant: str,
    root: str,
    *,
    per_sample: bool = False,
    ranked_bin_set: Optional[str] = None,
    abundance_truth: Optional[Dict[str, Dict[str, float]]] = None,
) -> List[_scoring.Row]:
    """Read whatever scoring output a variant's directory happens to hold.

    Discovery by filename rather than by configuration: CheckM2, GUNC and AMBER each write a file
    with a stable name, and a benchmark that required the user to say where each one landed would
    be one more thing to get wrong between a run and its score.

    *per_sample* says the runner put one subdirectory per sample under *root*, so its name is the
    sample label. It is passed in rather than guessed at, because a cohort module's subdirectories
    are named after its own working files and labelling those as samples would invent samples that
    do not exist.

    *ranked_bin_set* names the one bin set whose numbers are the variant's own (see
    :func:`_ranked_bin_set`). Every other bin set's numbers are recorded under a suffixed metric
    name so they are kept as detail without being summed into the ranking.

    *abundance_truth* is the expected fraction per bin per sample, for a quantifier comparison. It
    has to be passed in rather than discovered here because, unlike every other file below, it is not
    written by the run: it is built from the dataset's own read mappings, once, and reused for every
    variant so the three quantifiers are held to one truth.
    """
    rows: List[_scoring.Row] = []
    # read_qc's two mate files contribute to one number per sample, so they are summed here rather
    # than recorded as two rows: collect_ratios keys a measurement on (dataset, comparison, variant,
    # sample, metric) and deliberately not on the source file, so two rows collapsed to whichever
    # mate was walked last and every read_qc figure was one mate's count. Keyed (sample, metric).
    read_totals: Dict[Tuple[str, str], float] = {}
    read_sources: Dict[Tuple[str, str], List[str]] = {}

    def add(metric: str, value: float, source: str, sample: str = "") -> None:
        rows.append(
            _scoring.Row(
                dataset=dataset,
                comparison=comparison,
                variant=variant,
                metric=metric,
                value=value,
                sample=sample,
                source=source,
            )
        )

    for base, _dirs, names in os.walk(root):
        relative = os.path.relpath(base, root).split(os.sep)[0]
        sample = relative if per_sample and relative != os.curdir else ""
        for name in sorted(names):
            path = os.path.join(base, name)
            if name.endswith(_STATS_SUFFIX):
                # MetaWrap2's own normalised table, written by either CheckM. Preferred over the
                # tool's native report below so one CheckM run contributes one set of rows: CheckM2
                # writes both, and reading both would double every count.
                label = name[: -len(_STATS_SUFFIX)]
                counts = _scoring.quality_counts(_scoring.bin_quality(path))
                for metric, count in sorted(counts.items()):
                    add(_for_bin_set(metric, label, ranked_bin_set), float(count), name, sample)
            elif name in ("quality_report.tsv",) or name.endswith(".checkm.tsv"):
                if _has_sibling_stats(path):
                    continue
                counts = _scoring.quality_counts(_scoring.bin_quality(path))
                for metric, count in sorted(counts.items()):
                    add(metric, float(count), name, sample)
            elif "maxCSS_level" in name and name.endswith(".tsv"):
                fraction = _scoring.gunc_pass_fraction(path)
                if fraction is not None:
                    label = _gunc_bin_set(path, root, sample)
                    add(
                        _for_bin_set("gunc_pass_fraction", label, ranked_bin_set),
                        fraction,
                        name,
                        sample,
                    )
            elif name == "results.tsv":
                # One AMBER run can hold several bin sets, tagged with the label they were submitted
                # under, so they are read per set rather than folded together.
                by_tool = _scoring.amber_metrics_by_tool(path)
                for label, metrics in sorted(by_tool.items()):
                    for metric, score in sorted(metrics.items()):
                        add(_for_bin_set(metric, label, ranked_bin_set), score, name, sample)
            elif name == "bin_abundance_table.tsv" and abundance_truth:
                # The one metric here that is not read off the run's own output: the table says what
                # the quantifier thought, and only the dataset's read mappings say what was true.
                # Scored whole rather than per sample - a correlation across two bins is noise, so
                # the module averages the samples itself.
                try:
                    scored = _abundance.score(path, abundance_truth)
                except _abundance.AbundanceError as exc:
                    print("    %s could not be scored: %s" % (name, exc))
                    continue
                for metric, value in sorted(scored.items()):
                    add(metric, value, name, sample)
            elif name in ("final_assembly.fasta", "final_assembly.fa"):
                for metric, stat in sorted(_scoring.assembly_stats(path).items()):
                    add(metric, float(stat), name, sample)
            elif name.startswith("final_pure_reads_") and name.endswith(
                (".fastq", ".fq", ".fastq.gz", ".fq.gz")
            ):
                # read_qc's product, one file per mate. Accumulated and emitted once below, so a
                # sample's reads_kept is both mates rather than whichever file was walked last.
                for metric, stat in sorted(_scoring.read_stats(path).items()):
                    key = (sample, metric)
                    read_totals[key] = read_totals.get(key, 0.0) + float(stat)
                    read_sources.setdefault(key, []).append(name)
    for (sample, metric), total in sorted(read_totals.items()):
        add(metric, total, "+".join(sorted(read_sources[(sample, metric)])), sample)
    return rows


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 benchmark",
        description="Measure MetaWrap2's options against known-truth datasets.",
    )
    ap.add_argument("--config", help="path to a metawrap2.toml")
    sub = ap.add_subparsers(dest="action", required=True)

    p_list = sub.add_parser("list", help="show the datasets and comparisons that are declared")
    p_list.add_argument(
        "-v", "--verbose", action="store_true", help="also show citations and gold-standard files"
    )

    p_fetch = sub.add_parser("fetch", help="download a dataset, checksum-verified")
    p_fetch.add_argument("dataset")
    p_fetch.add_argument("--dest", help="where to put it (default benchmarks/data/<dataset>)")
    p_fetch.add_argument(
        "--yes",
        action="store_true",
        help="confirm the download size and start (without this, the plan is printed and nothing "
        "is downloaded)",
    )
    p_fetch.add_argument(
        "--verify",
        action="store_true",
        help="re-checksum files already present instead of trusting their size",
    )
    p_fetch.add_argument(
        "--keep-archives",
        action="store_true",
        help="keep each .tar.gz after extracting it (roughly doubles the disk needed)",
    )

    p_run = sub.add_parser("run", help="run one comparison's variants over one dataset")
    p_run.add_argument("dataset")
    p_run.add_argument("comparison")
    p_run.add_argument(
        "--data", help="where the dataset was fetched (default benchmarks/data/<ds>)"
    )
    # --work is a synonym because `score` names the same directory that way, and a flag that only
    # works on one of two commands you type in sequence is a papercut with no upside.
    p_run.add_argument(
        "-o", "--output", "--work", help="working directory (default benchmarks/work/<dataset>)"
    )
    p_run.add_argument(
        "--variant",
        action="append",
        help="run only this variant (repeatable; default is every variant)",
    )
    p_run.add_argument(
        "--sample",
        action="append",
        help="run only this sample (repeatable; default is every sample in the dataset)",
    )
    p_run.add_argument(
        "--max-samples",
        type=int,
        metavar="N",
        help="use only the first N samples - for comparisons too expensive over the whole dataset",
    )
    p_run.add_argument(
        "--assembly",
        help="assembly to bin, instead of the dataset's gold-standard one",
    )
    p_run.add_argument(
        "--bins-from", help="a previous binning comparison's output, for bin_refinement"
    )
    p_run.add_argument("--results", help="results TSV to append to")

    p_score = sub.add_parser("score", help="read tool output into the results TSV and rank it")
    p_score.add_argument("dataset")
    p_score.add_argument("comparison")
    p_score.add_argument(
        "--work", "-o", "--output", help="the run's working directory (same one `run` was given)"
    )
    p_score.add_argument(
        "--data", help="where the dataset was fetched (default: whatever the run recorded)"
    )
    p_score.add_argument(
        "--assembly",
        help="the contigs the binners were given, for AMBER's sequence lengths (default: whatever "
        "the run recorded, else the dataset's gold-standard assembly)",
    )
    p_score.add_argument(
        "--reassess",
        action="store_true",
        help="re-run CheckM/AMBER on bins that already have results, for when the bins changed",
    )
    p_score.add_argument(
        "--checkm2",
        action="store_true",
        help="score bin quality with CheckM2 instead of CheckM1 (needs `metawrap2 install-env "
        "checkm2`). CheckM1's pplacer stage wants roughly 40 GB per thread; CheckM2 "
        "predicts completeness from a neural net and needs a fraction of that",
    )
    p_score.add_argument(
        "--score-failed",
        action="store_true",
        help="score variants whose run exited non-zero too. Their output is partial, so the number "
        "is not comparable with an arm that finished; skipped by default",
    )
    p_score.add_argument("--results", help="results TSV to read and append to")
    p_score.add_argument(
        "--metric",
        action="append",
        help="rank by this metric (repeatable; default is the comparison's declared metrics)",
    )
    p_score.add_argument(
        "--cost-metric",
        default="cpu_hours",
        help="the denominator of the ratio (default cpu_hours, which is what MetaWrap2 picks "
        "defaults on; wall_seconds rewards a tool for being given more cores)",
    )

    args = ap.parse_args(argv)
    try:
        return {"list": _list, "fetch": _fetch_cmd, "run": _run, "score": _score}[args.action](args)
    except _registry.RegistryError as exc:
        print("\n%s\n" % exc)
        return 2
