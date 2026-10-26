"""`metawrap2 benchmark list|fetch|run|score`: benchmark MetaWrap2's options on known-truth data.

Defaults are chosen on measured quality, with cost reported alongside. Datasets and comparisons
are declared as data, and results go to a committed TSV so a changed default is reviewed with its
numbers.

Subcommands, in order of use:

- ``list``   declared datasets and comparisons, and the download size of each dataset.
- ``fetch``  download one dataset, checksum-verified, within the ``[limits]`` disk ceiling.
- ``run``    run one comparison's variants over one dataset under one budget, recording cost.
- ``score``  read the tools' output into the TSV and rank the variants, with cost alongside.

``fetch`` never runs implicitly. ``run`` on a missing dataset prints how to fetch it and stops
(a dataset can be 95 GB).
"""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Sequence, Tuple

from .. import checkm as _checkm
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
    """Default results file for a dataset. It is committed, so a regression shows as a diff.

    Relative to the registry (the checkout), not the current directory. A relative path wrote a
    second results file when scoring was run outside the repository root.
    """
    return os.path.join(_registry.registry_dir(), "results", "%s.tsv" % dataset_name)


def _dry_run() -> bool:
    """True if the global ``--dry-run`` is in effect.

    Read from the shared runner: the dispatcher strips ``--dry-run`` before this parser sees it, so
    a local option would always be False.
    """
    return bool(_command.runner.dry_run)


def _budget(args) -> _limits.Budget:
    """The budget for this run: the config's ``[limits]`` with the global ``--max-*`` flags.

    Resolved once so every variant of a comparison runs under the same budget.
    """
    settings = load_settings(getattr(args, "config", None))
    return _limits.resolve(settings.limits)


def _variant_detail(variant) -> str:
    """What to show after an arm's name in ``benchmark list``.

    A native arm shows its flags. An external arm has no flags, so it shows its version and the
    env it runs in.
    """
    if not variant.is_external:
        return " ".join(variant.args)
    external = variant.external
    version = external.version or "version from `%s`" % external.version_command
    return "external: %s in env %s" % (version, external.env)


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
            print("      %-30s %s%s" % (variant.name, _variant_detail(variant), mark))
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
        # Checked before asking for confirmation.
        print(
            "\nNot started: the download would exceed the [limits] disk ceiling.\n"
            "Choose a --dest on another filesystem, or raise max_disk_fraction under [limits]."
        )
        return 1
    # A dataset is ~95 GB, so the download needs --yes.
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

    Records from older runs may lack a field, so a missing record or key returns None.
    """
    path = os.path.join(work_dir, "%s.run.json" % comparison)
    try:
        with open(path) as fh:
            return json.load(fh).get(field)
    except (OSError, ValueError):
        return None


def _say_sample_scope(work_dir: str, comparison: str) -> None:
    """Print which samples the run covered, if its record is still there.

    Cohort module rows have no sample label, and `--max-samples` can limit the samples used, so the
    count is printed here. Prints nothing if there is no record.
    """
    path = os.path.join(work_dir, "%s.run.json" % comparison)
    try:
        with open(path) as fh:
            samples = json.load(fh).get("samples") or []
    except (OSError, ValueError):
        return
    if samples:
        print("\n  Samples measured: %d  (%s)" % (len(samples), ", ".join(samples)))


def _say_external_arms(work_dir: str, comparison: str) -> None:
    """Print which arms were foreign pipelines, and which version of each ran.

    Read from the run record, which holds the version that ran; the registry only says what would
    run.
    """
    arms = _record_field(work_dir, comparison, "variants") or {}
    external = {
        name: record.get("external_version") or "an unrecorded version"
        for name, record in arms.items()
        if record.get("external")
    }
    if not external:
        return
    print("\n  External arm(s) - not MetaWrap2:")
    for name, version in sorted(external.items()):
        print("    %-28s %s" % (name, version))


def _score(args) -> int:
    dataset = _registry.get_dataset(args.dataset)
    comparison = _registry.get_comparison(args.comparison)
    work_dir = args.work or os.path.join("benchmarks", "work", dataset.name)
    path = args.results or default_results(dataset.name)

    metrics_wanted = args.metric or list(comparison.metrics)
    budget = _budget(args)
    # Use the data directory from the run record; --data overrides it and covers older records.
    data_dir = (
        args.data
        or _record_field(work_dir, comparison.name, "data_dir")
        or os.path.join("benchmarks", "data", dataset.name)
    )
    # Contigs the binners were given, for AMBER's sequence lengths. Taken from the record because
    # `run --assembly` may have replaced the gold assembly, which is the fallback.
    assembly = args.assembly or _record_field(work_dir, comparison.name, "assembly")
    if not assembly or not os.path.isfile(assembly):
        assembly = _runner.gold_path(dataset, "assembly", data_dir)

    # Built once from the dataset, and only for abundance metrics, since it decompresses every
    # sample's read mapping (minutes). All variants are scored against the same truth.
    abundance_truth = None
    if any(m.startswith("abundance_") for m in metrics_wanted):
        abundance_truth = _abundance_truth(
            dataset, comparison, work_dir, data_dir=data_dir, assembly=assembly
        )

    # Variants that exited non-zero, from the run record. Their output is partial, so they are
    # skipped unless --score-failed; otherwise partial bins were scored and ranked.
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
            checkm2=bool(getattr(args, "checkm2", True)),
        )
        found = _collect_scores(
            dataset.name,
            comparison.name,
            variant.name,
            root,
            per_sample=comparison.module in _runner.PER_SAMPLE_MODULES,
            ranked_bin_set=_ranked_bin_set(bin_sets),
            abundance_truth=abundance_truth,
            scorer="checkm2" if getattr(args, "checkm2", True) else "checkm1",
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
            "\n  Not scored: %d variant(s) whose run failed, so their output is partial:"
            % len(skipped)
        )
        for name, rc in skipped:
            print("    %-28s exited %s" % (name, rc))
        print("    Re-run them, or pass --score-failed to score what is there anyway.")

    if unranked:
        print(
            "\n  Not ranked: %d variant(s) left several raw bin sets and no consolidated one:"
            % len(unranked)
        )
        for name, labels in unranked:
            print("    %-28s %s" % (name, ", ".join(labels)))
        print(
            "  Their numbers are recorded per bin set (metric__<set>) and not summed, because\n"
            "  the sets overlap and a sum would count a genome more than once.\n"
            "  Consolidate them with bin_refinement and score again."
        )

    if rows:
        written = _scoring.write_tsv(rows, path, append=True)
        print("\nRecorded %d accuracy row(s) in %s" % (written, path))

    if not os.path.isfile(path):
        print("\nNo results file at %s yet.\n" % path)
        return 1

    _say_sample_scope(work_dir, comparison.name)
    _say_external_arms(work_dir, comparison.name)
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

    Returns None instead of raising, so the other metrics are still scored.

    The bins come from the run record's ``bins_from``. Every variant quantified the same bins, so
    the truth is built once on one set of bin labels.
    """
    bins_from = _record_field(work_dir, comparison.name, "bins_from")
    if not bins_from or not os.path.isdir(bins_from):
        print(
            "  Abundance metrics not scored: the run record for %s names no existing\n"
            "  --bins-from directory, and the abundance truth is computed per bin."
            % comparison.name
        )
        return None
    contig_mapping = _runner.gold_path(dataset, "binning", data_dir)
    if not (contig_mapping and os.path.isfile(contig_mapping)) or not (
        assembly and os.path.isfile(assembly)
    ):
        print(
            "  Abundance metrics not scored: this needs the gold assembly and its\n"
            "  contig-to-genome mapping, to compute how much of each genome is in each bin."
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
    checkm2: bool = True,
) -> List[str]:
    """Run whichever assessors *metrics* needs over the bin sets under *root*. Returns those sets.

    Assessment runs at score time, outside the measured run (see :mod:`metawrap2.benchmark.assess`).

    Already-assessed sets are skipped so scoring can be re-run cheaply. ``--reassess`` forces a
    rerun when the bins have changed.
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
    scorer = "checkm2" if checkm2 else "checkm1"
    outstanding = (
        bin_sets if reassess else [p for p in bin_sets if not _already_assessed(p, metrics, scorer)]
    )
    # AMBER writes one table per variant, not per bin set, so it is checked separately. Otherwise
    # sets with CheckM tables never got AMBER run on them.
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


def _already_assessed(bin_set: str, metrics: Sequence[str], scorer: str = "checkm2") -> bool:
    """True if *bin_set* has the output of every per-set assessor (CheckM, GUNC) *metrics* needs.

    AMBER writes one table per variant and is handled by the caller. A set with a CheckM table but
    no GUNC report is not assessed, and neither is one whose table *scorer* (``checkm1`` or
    ``checkm2``) did not write. A table with no scored bin counts for either scorer.
    """
    present = []
    if _assess.needs_checkm(metrics):
        stats = bin_set + _STATS_SUFFIX
        written_by = _checkm.stats_scorer(stats)
        present.append(os.path.isfile(stats) and written_by in (None, scorer))
    if _assess.needs_gunc(metrics):
        present.append(_gunc.find_report(bin_set + _GUNC_SUFFIX) is not None)
    # With no per-set assessor requested the set counts as done; the AMBER-only case is the
    # caller's.
    return all(present)


#: Suffix of MetaWrap2's normalised bin-quality table.
_STATS_SUFFIX = ".stats.tsv"

#: A consolidated bin set is recognised by its name pattern (metawrap_70_10_bins, from the
#: thresholds used). When present it is the variant's result. See
#: :func:`metawrap2.constants.is_refined_bins_dirname`.


def _ranked_bin_set(bin_sets: Sequence[str]) -> Optional[str]:
    """The bin set whose metrics are ranked as the variant's, or None.

    Bin sets from several binners overlap, and the ranker sums a metric over everything recorded
    for a variant, so ranking them all would count shared genomes more than once (an F1 of 2.4, for
    example). The ranked set is:

    - the consolidated set, if there is exactly one;
    - otherwise the only set, if there is one;
    - otherwise none, and the caller reports the variant as unranked.
    """
    labels = [os.path.basename(path) for path in bin_sets]
    consolidated = [label for label in labels if is_refined_bins_dirname(label)]
    if len(consolidated) == 1:
        return consolidated[0]
    if consolidated:
        # Several consolidated sets (a run resumed at different thresholds): do not pick one.
        return None
    if len(labels) == 1:
        return labels[0]
    return None


def _for_bin_set(metric: str, label: str, ranked: Optional[str]) -> str:
    """*metric* as recorded for the bin set *label*: plain if it is the ranked set, suffixed if not.

    The suffix keeps per-binner numbers in the TSV under a different metric name, so they are not
    summed into the ranking.
    """
    if not label or (ranked is not None and label == ranked):
        return metric
    return "%s__%s" % (metric, label)


#: Where bin_refinement writes GUNC's report on its output set. The module runs GUNC once, on that
#: set, so the report is the variant's own figure.
_GUNC_FINAL_DIR = os.path.join("work_files", "gunc")

#: Suffix of the directory that holds a bin set's score-time GUNC report (``<bin set>.gunc``, like
#: ``<bin set>.checkm``).
_GUNC_SUFFIX = ".gunc"


def _gunc_bin_set(path: str, root: str, sample: str) -> str:
    """Which bin set the GUNC report at *path* describes, as a label for :func:`_for_bin_set`.

    GUNC's file names do not identify the bins, so the label comes from the directory:

    * ``<bin set>.gunc/`` (score time) gives ``<bin set>``, the same label as AMBER's rows.
    * ``work_files/gunc/`` (bin_refinement) gives the empty label, which :func:`_for_bin_set`
      treats as the variant's own figure.

    Any other directory keeps its relative path as the label, so the number is recorded but not
    ranked. The ranker sums values under one metric name, and two pass fractions of 1.0 would rank
    as 2.0.
    """
    directory = os.path.dirname(path)
    if directory.endswith(_GUNC_SUFFIX):
        return os.path.basename(directory)[: -len(_GUNC_SUFFIX)]
    relative = os.path.relpath(directory, root)
    if sample:
        relative = os.path.relpath(relative, sample)
    return "" if relative == _GUNC_FINAL_DIR else relative


def _has_sibling_stats(path: str) -> bool:
    """True if a normalised ``.stats.tsv`` exists for the same bins as the native report at *path*.

    CheckM writes its report under ``<bin set>.checkm/`` and MetaWrap2 writes
    ``<bin set>.stats.tsv`` beside it. Reading both would count every bin twice.
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
    scorer: Optional[str] = None,
) -> List[_scoring.Row]:
    """Read the scoring output found in a variant's directory.

    Files are found by their fixed names (CheckM, GUNC, AMBER and module outputs).

    *per_sample* means the runner made one subdirectory per sample under *root*, named after the
    sample. Cohort modules' subdirectories are working files, not samples.

    *ranked_bin_set* is the bin set whose numbers are the variant's own (see
    :func:`_ranked_bin_set`). Other sets are recorded under suffixed metric names.

    *abundance_truth* is the expected fraction per bin per sample for a quantifier comparison. It
    is built from the dataset's read mappings, not the run, and shared by all variants.

    *scorer* (``checkm1`` or ``checkm2``) skips ``.stats.tsv`` tables the other CheckM wrote, such
    as an old table left when re-scoring failed.
    """
    rows: List[_scoring.Row] = []
    # read_qc's two mate files are summed into one value per (sample, metric). collect_ratios keys
    # rows without the source file, so two rows would collapse to one mate's count.
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
                # MetaWrap2's normalised table from either CheckM. Preferred over the native report
                # below so each CheckM run is counted once.
                if scorer and _checkm.stats_scorer(path) not in (None, scorer):
                    continue
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
                # One AMBER table can hold several bin sets; read them per label.
                by_tool = _scoring.amber_metrics_by_tool(path)
                for label, metrics in sorted(by_tool.items()):
                    for metric, score in sorted(metrics.items()):
                        add(_for_bin_set(metric, label, ranked_bin_set), score, name, sample)
            elif name == "bin_abundance_table.tsv" and abundance_truth:
                # Scored against the truth built from the dataset's read mappings. The scorer
                # averages over samples itself, so the whole table is scored at once.
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
                # read_qc output, one file per mate. Summed and emitted once below.
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
    # --work is accepted because `score` uses that name for the same directory.
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
        help="use only the first N samples, for comparisons too expensive to run on every sample",
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
    scorer = p_score.add_mutually_exclusive_group()
    scorer.add_argument(
        "--checkm1",
        dest="checkm2",
        action="store_false",
        help="deprecated: score bin quality with CheckM1 on its full genome tree instead of "
        "CheckM2 (about 40 GB of RAM per pplacer thread; needs `metawrap2 install-db checkm`)",
    )
    scorer.add_argument(
        "--checkm2",
        dest="checkm2",
        action="store_true",
        help="score bin quality with CheckM2 (the default)",
    )
    p_score.set_defaults(checkm2=True)
    p_score.add_argument(
        "--score-failed",
        action="store_true",
        help="also score variants whose run exited non-zero (default: skip them). Their output is "
        "partial, so the scores are not comparable with finished variants",
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
        help="cost metric shown beside each ranked metric (default cpu_hours, which MetaWrap2 "
        "uses when choosing defaults; wall_seconds favours runs given more cores)",
    )

    args = ap.parse_args(argv)
    try:
        return {"list": _list, "fetch": _fetch_cmd, "run": _run, "score": _score}[args.action](args)
    except _registry.RegistryError as exc:
        print("\n%s\n" % exc)
        return 2
