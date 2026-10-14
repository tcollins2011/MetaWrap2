"""Running one comparison's variants over one dataset, under one budget, and recording the cost.

The discipline this module exists to enforce is that only the thing being compared differs. Every
variant of a comparison gets the same inputs, the same thread and memory grant, and the same
machine; what changes is the flags under test. Without that, a results table is a table of
coincidences - a binner that happened to run when the node was quiet looks efficient.

Three specifics follow from it:

* **The budget is resolved once and handed to every variant.** Not re-resolved per variant, so a
  run that starts while something else is finishing does not give the first variant half a machine
  and the second one all of it.
* **Binning comparisons default to the gold-standard assembly.** CAMI2 ships one precisely so a
  binner can be measured without an assembler's variance in front of it. Comparing binners on
  their own assemblies measures the pair, which is a different question and a noisier one.
* **Cost is measured, not timed by hand.** :func:`metawrap2.usage.run_measured` reaps each child
  with ``wait4``, so peak RSS and CPU seconds are the step's own rather than a running maximum
  across everything the process has ever waited for.

What this module does not do is decide anything. It produces rows;
:func:`metawrap2.benchmark.scoring.value_against_cost` puts them beside the result, and a person
reads that.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence

from .. import config as _config
from .. import limits as _limits
from .. import usage as _usage
from ..command import runner as _runner
from ..constants import is_refined_bins_dirname
from ..io import reads as _reads
from ..utils import strip_trailing_slash
from .registry import Comparison, Dataset, RegistryError, Variant
from .scoring import Row

#: How many bin sets bin_refinement will consolidate. Its own -A..-F flags are the ceiling, and the
#: real limit is the 2^n-1 CheckM runs behind them: six sets is 63.
MAX_REFINEMENT_SETS = 6

#: Modules the runner knows how to feed, and what each needs beyond reads. Kept as an explicit
#: table rather than a generic argument builder: a benchmark that guessed wrong about a module's
#: arguments would produce numbers, and numbers from the wrong command are worse than an error.
SUPPORTED_MODULES = {
    "read_qc": "reads, one sample at a time",
    "assembly": "reads",
    "binning": "reads + an assembly (the gold standard one by default)",
    "bin_refinement": "bin sets from a previous binning run (--bins-from)",
    "quant_bins": "reads + one set of bins from a previous run (--bins-from)",
}

#: Modules that take one sample per invocation rather than a whole cohort. These are run once per
#: sample and their cost rows carry the sample label, so a ten-sample dataset yields ten measured
#: runs that sum into one ratio - the same arithmetic as a cohort module, from more evidence.
PER_SAMPLE_MODULES = frozenset({"read_qc"})


class RunnerError(Exception):
    """A comparison that cannot be run as asked."""


@dataclass
class Inputs:
    """The prepared inputs a dataset offers, discovered on disk rather than assumed."""

    pairs: List[Sequence[str]]  # (sample, r1, r2)
    assembly: Optional[str] = None
    bins_from: Optional[str] = None

    @property
    def read_files(self) -> List[str]:
        out: List[str] = []
        for _sample, r1, r2 in self.pairs:
            out.extend([r1, r2])
        return out


def _find(root: str, suffixes: Sequence[str]) -> List[str]:
    """Every file under *root* ending in one of *suffixes*, sorted.

    Sorted because the order reads are passed to a binner changes nothing about the result but
    everything about whether two runs are byte-comparable.
    """
    found: List[str] = []
    for base, _dirs, names in os.walk(root):
        for name in sorted(names):
            if name.endswith(tuple(suffixes)):
                found.append(os.path.join(base, name))
    return sorted(found)


def read_candidates(dataset: Dataset, source: str) -> List[str]:
    """The dataset's read files: its declared ``reads_glob`` if it has one, else every FASTQ.

    The glob matters more than it looks. CAMI2's ``reads/`` directory holds the 10 GB interleaved
    FASTQ *and* a 3 MB ``subsample_0.fq`` debugging extract, both under the same ``sample_N``
    directory - so a plain suffix scan finds two files that resolve to one sample label, and
    whichever sorts last silently wins. Benchmarking a trimmer on 3 MB instead of 10 GB produces
    numbers, and numbers off the wrong file are worse than no numbers.
    """
    if dataset.reads_glob:
        found = sorted(glob.glob(os.path.join(source, dataset.reads_glob)))
        if not found:
            raise RunnerError(
                "%s declares reads_glob = %r, which matches nothing under %s.\n"
                "Either the data is not extracted yet or the glob needs updating in "
                "benchmarks/datasets.toml." % (dataset.name, dataset.reads_glob, source)
            )
        return found
    return _find(source, (".fq.gz", ".fastq.gz", ".fq", ".fastq"))


def abundance_truth_files(dataset: Dataset, data_dir: str) -> Dict[str, str]:
    """The dataset's per-sample read-to-genome mappings, keyed by the same sample label as the reads.

    Keyed by :func:`_sample_label` rather than by directory name, because that is the label the
    quantifier's table columns carry - the join between an estimate and its truth is that label, and
    deriving both ends of it from one function is the only way to be sure they agree.
    """
    if not dataset.abundance_truth_glob:
        raise RunnerError(
            "%s declares no abundance_truth_glob, so there is no per-sample read-to-genome truth "
            "to score a quantifier against.\nAdd one in benchmarks/datasets.toml, or compare the "
            "quantifiers on a dataset that has one." % dataset.name
        )
    source = dataset.local_root if dataset.is_local else os.path.join(data_dir, "extracted")
    found = sorted(glob.glob(os.path.join(source, dataset.abundance_truth_glob)))
    if not found:
        raise RunnerError(
            "%s declares abundance_truth_glob = %r, which matches nothing under %s."
            % (dataset.name, dataset.abundance_truth_glob, source)
        )
    truth: Dict[str, str] = {}
    for path in found:
        sample = _sample_label(path)
        if sample in truth:
            raise RunnerError(
                "Two abundance truth files both resolve to sample %r:\n  %s\n  %s\n"
                "Narrow %s's abundance_truth_glob in benchmarks/datasets.toml."
                % (sample, truth[sample], path, dataset.name)
            )
        truth[sample] = path
    return truth


def _narrowing_label(path: str, layout: str) -> str:
    """The *sample* a read file belongs to, for deciding which files to prepare.

    Not :func:`_sample_label`, and the difference is the mate suffix. That function labels a file,
    which is what a results row and an abundance-truth join need; this one has to label a sample,
    because a paired dataset has two files per sample and narrowing on the file label selects one
    sample's two mates when asked for two samples.

    The suffix is stripped only for a paired layout. An interleaved dataset has one file per sample
    and CAMI2's labels are ``sample_0``, ``sample_1``, ``sample_2`` - stripping a trailing ``_1``
    there would merge two different samples into one.
    """
    label = _sample_label(path)
    if layout == "interleaved":
        return label
    for suffix in ("_R1", "_R2", "_r1", "_r2", "_1", "_2", ".1", ".2"):
        if label.endswith(suffix):
            return label[: -len(suffix)]
    return label


def _narrow_candidates(
    candidates: Sequence[str],
    wanted: Sequence[str],
    max_samples: int,
    say: Callable[[str], None],
    layout: str = "interleaved",
) -> List[str]:
    """The read files worth preparing, given the sample selection.

    Named samples that are absent are *not* an error here: :func:`select_samples` reports those
    against the full list of what the dataset has, and duplicating the check would report it
    against a list this function has already shortened.
    """
    if not wanted and max_samples <= 0:
        return list(candidates)
    # Paths per sample, not one path per sample. A paired dataset has two files under one label, so
    # keying a dict on the label keeps whichever came last - every mate file was silently dropped,
    # and the pairing downstream then read one file per sample and mislabelled the lot: a
    # three-sample run narrowed to two reported one. Grouping keeps both mates of every sample kept.
    by_sample: Dict[str, List[str]] = {}
    for path in candidates:
        by_sample.setdefault(_narrowing_label(path, layout), []).append(path)
    order = sorted(by_sample)
    if wanted:
        order = [name for name in order if name in set(wanted)]
    if max_samples > 0:
        order = order[:max_samples]
    if not order:
        return list(candidates)
    skipped = len(by_sample) - len(order)
    if skipped > 0:
        say(
            "  preparing reads for %d of %d samples (%s); the rest are not part of this run"
            % (len(order), len(by_sample), ", ".join(order))
        )
    kept = set(order)
    # Original order, so the files are prepared in the order they were found rather than in label
    # order - the downstream pairing matches mates by name and does not care, but a stable order
    # keeps a dry run's printed plan comparable between runs.
    return [path for path in candidates if _narrowing_label(path, layout) in kept]


def prepare_reads(
    dataset: Dataset,
    data_dir: str,
    work_dir: str,
    *,
    dry_run: bool = False,
    say: Callable[[str], None] = print,
    wanted: Sequence[str] = (),
    max_samples: int = 0,
) -> List[Sequence[str]]:
    """Locate the dataset's reads and, if interleaved, split them once into *work_dir*.

    *wanted* and *max_samples* are the same narrowing :func:`select_samples` applies, passed in
    here so the splitting is narrowed too. Without that, a comparison asked for three samples
    deinterleaved all ten first: 190 GB written and an hour spent to use 57 GB of it, on a
    filesystem the harness otherwise guards with a disk preflight. The selection still happens in
    one place - this only decides which files are worth writing.

    Splitting is cached: the second variant of a comparison finds the pairs already there. That is
    not just a speed-up - deinterleaving twice would put the second variant's reads on a colder
    page cache than the first's, which is exactly the kind of difference that turns into a
    fictitious performance gap.

    Under *dry_run* the paths are named but nothing is written. For CAMI2 that is the difference
    between printing six commands and spending an hour writing 100 GB of mate files - and the
    point of a dry run is to see the plan before committing to the work.
    """
    source = dataset.local_root if dataset.is_local else os.path.join(data_dir, "extracted")
    if not os.path.isdir(source):
        raise RunnerError(
            "No data for %s under %s.\nFetch it first: metawrap2 benchmark fetch %s"
            % (dataset.name, source, dataset.name)
        )
    candidates = read_candidates(dataset, source)
    if not candidates:
        raise RunnerError("Found no FASTQ under %s for %s." % (source, dataset.name))
    candidates = _narrow_candidates(
        candidates, wanted, max_samples, say, layout=dataset.reads_layout
    )

    if dataset.reads_layout == "interleaved":
        split_dir = os.path.join(work_dir, "reads")
        if not dry_run:
            os.makedirs(split_dir, exist_ok=True)
        pairs: List[Sequence[str]] = []
        claimed: Dict[str, str] = {}
        for path in candidates:
            sample = _sample_label(path)
            if sample in claimed:
                # A backstop for the glob above: two files under one label would mean the second
                # deinterleave overwrites the first, and the comparison would run on whichever
                # file happened to sort last rather than on the one meant.
                raise RunnerError(
                    "Two read files both resolve to sample %r:\n  %s\n  %s\n"
                    "Narrow %s's reads_glob in benchmarks/datasets.toml so it names one file per "
                    "sample." % (sample, claimed[sample], path, dataset.name)
                )
            claimed[sample] = path
            r1 = os.path.join(split_dir, "%s_1.fastq" % sample)
            r2 = os.path.join(split_dir, "%s_2.fastq" % sample)
            if os.path.isfile(r1) and os.path.isfile(r2) and os.path.getsize(r1):
                say("  reads for %s already split" % sample)
            elif dry_run:
                say("  (dry run) would split interleaved reads for %s" % sample)
            else:
                say("  splitting interleaved reads for %s" % sample)
                _reads.deinterleave(path, r1, r2)
            pairs.append((sample, r1, r2))
        return pairs

    try:
        triples = _reads.collect_pairs(candidates)
    except ValueError as exc:
        # collect_pairs already explains which conventions it recognises; re-raising as a
        # RunnerError keeps the benchmark's own errors in one class without losing that.
        raise RunnerError("Could not pair the reads under %s.\n%s" % (source, exc))
    return sorted(triples, key=lambda triple: triple[0])


def _sample_label(path: str) -> str:
    """A stable sample label from a CAMI2-style path.

    The published directories are timestamp-prefixed (``2017.12.04_18.45.54_sample_0``), and the
    timestamp differs per sample, so the label comes from the ``sample_N`` part rather than from
    the whole directory name - otherwise every results row would carry a publication timestamp
    nobody needs and no two body sites would sort together.
    """
    for part in os.path.abspath(path).split(os.sep):
        if "sample_" in part:
            tail = part.split("sample_", 1)[1]
            digits = "".join(c for c in tail if c.isdigit())
            if digits:
                return "sample_%s" % digits
    return os.path.splitext(os.path.basename(path))[0].replace(".fq", "").replace(".fastq", "")


def gather_inputs(
    dataset: Dataset,
    comparison: Comparison,
    data_dir: str,
    work_dir: str,
    *,
    assembly: Optional[str] = None,
    bins_from: Optional[str] = None,
    dry_run: bool = False,
    say: Callable[[str], None] = print,
    wanted: Sequence[str] = (),
    max_samples: int = 0,
) -> Inputs:
    """Everything *comparison* needs to run over *dataset*, or a clear error saying what is absent.

    *wanted* and *max_samples* are passed through to :func:`prepare_reads` so that preparing the
    reads costs what the run costs rather than what the dataset holds.
    """
    if comparison.module not in SUPPORTED_MODULES:
        raise RunnerError(
            "The benchmark runner cannot feed the %r module yet. It knows how to run:\n%s"
            % (
                comparison.module,
                "\n".join("  %-16s %s" % (k, v) for k, v in sorted(SUPPORTED_MODULES.items())),
            )
        )
    pairs = prepare_reads(
        dataset,
        data_dir,
        work_dir,
        dry_run=dry_run,
        say=say,
        wanted=wanted,
        max_samples=max_samples,
    )

    resolved_assembly = assembly
    # bin_refinement needs it too, for the same reason the command builder passes it: DAS Tool reads
    # the assembly contigs itself and refuses to start without them, so a consolidator comparison
    # could never run its DAS Tool arm. Resolved here rather than left to the caller, because the
    # command builder's `if inputs.assembly` was already right and still produced a command with no
    # -a - nothing had filled the field in.
    if comparison.module in ("binning", "bin_refinement") and resolved_assembly is None:
        candidate = gold_path(dataset, "assembly", data_dir)
        if candidate is None:
            if comparison.module == "bin_refinement":
                # Only DAS Tool needs it; the default consolidator works from the bins alone, so a
                # dataset without a gold assembly can still answer part of the question.
                candidate = ""
            else:
                raise RunnerError(
                    "%s has no gold-standard assembly, so a binning comparison needs one passed "
                    "with --assembly." % dataset.name
                )
        if candidate and not os.path.isfile(candidate):
            raise RunnerError(
                "Expected the gold-standard assembly at %s.\nFetch it: metawrap2 benchmark fetch %s"
                % (candidate, dataset.name)
            )
        resolved_assembly = candidate or None
        if not resolved_assembly:
            pass
        elif comparison.module == "binning":
            say(
                "  binning on the gold-standard assembly (%s), which isolates the binner from the\n"
                "  assembler's variance - pass --assembly to compare on a real assembly instead."
                % os.path.basename(candidate)
            )
        else:
            say(
                "  passing the gold-standard assembly (%s) to every arm; DAS Tool reads the contigs\n"
                "  itself and the arms have to be given the same inputs."
                % os.path.basename(candidate)
            )

    if comparison.module in ("bin_refinement", "quant_bins") and not bins_from:
        raise RunnerError(
            "A %s comparison needs the bins from a previous run. Run the binning comparison first, "
            "then pass its output with --bins-from." % comparison.module
        )
    return Inputs(pairs=pairs, assembly=resolved_assembly, bins_from=bins_from)


def one_bin_set(root: str) -> str:
    """The single bin directory under *root* a cohort module should be given.

    ``bin_refinement`` takes several bin sets and consolidates them; ``quant_bins`` takes one. *root*
    may already be that one - a folder of FASTA files - or it may be a variant's output directory
    holding several ``*_bins`` subdirectories, which is what a binning comparison leaves behind.

    Where there are several, the consolidated set wins and anything else is refused rather than
    picked. Quantifying one of three overlapping raw bin sets would produce a table that looks
    perfectly ordinary and answers a question nobody asked, and choosing between them arbitrarily is
    exactly the kind of silent decision a benchmark exists to avoid.
    """
    from ..io.seqio import fasta_files

    if fasta_files(root):
        return root
    sets = sorted(
        os.path.join(root, name)
        for name in os.listdir(root)
        if name.endswith("_bins") and os.path.isdir(os.path.join(root, name))
    )
    consolidated = [p for p in sets if is_refined_bins_dirname(p)]
    if len(consolidated) > 1:
        raise RunnerError(
            "%s holds %d consolidated bin sets (%s), so which one is the variant's answer is a "
            "guess. That is a refinement run resumed at different -c/-x thresholds; point "
            "--bins-from at the one you mean."
            % (root, len(consolidated), ", ".join(os.path.basename(p) for p in consolidated))
        )
    if consolidated:
        return consolidated[0]
    if len(sets) == 1:
        return sets[0]
    if not sets:
        raise RunnerError(
            "Found neither bins nor any *_bins directory under %s. --bins-from wants either a folder "
            "of bins or the output directory of a binning run." % root
        )
    raise RunnerError(
        "%s holds %d bin sets (%s) and no consolidated metawrap_<c>_<x>_bins. They overlap, so "
        "quantifying "
        "one of them would answer a different question from the one asked - refine them into one set "
        "first, or point --bins-from at the set you mean."
        % (root, len(sets), ", ".join(os.path.basename(p) for p in sets))
    )


def depth_table_for(root: str) -> str:
    """The contig depth table the binning run under *root* left behind.

    ``quant_bins --quantifier depth`` exists precisely so a cohort does not pay for a second round of
    alignment, so the benchmark arm has to be handed the table binning already wrote. Searched beside
    the bins and one level up, which covers both "--bins-from a binning output" and "--bins-from a
    refinement output inside it".
    """
    from ..modules.binning import DEPTH_TABLE

    candidates = [
        os.path.join(root, "work_files", DEPTH_TABLE),
        os.path.join(os.path.dirname(strip_trailing_slash(root)), "work_files", DEPTH_TABLE),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise RunnerError(
        "The depth quantifier needs the contig depth table from the binning run, and none of these "
        "exist:\n  %s\nThat table is written by `metawrap2 binning`, so point --bins-from at a "
        "binning output directory (or the refinement inside one) rather than at a bare folder of "
        "bins." % "\n  ".join(candidates)
    )


def _inherited_global_flags() -> List[str]:
    """``--force``/``--resume`` as typed on the benchmark's own command line.

    Each variant is a separate process, so the flags the user gave ``metawrap2`` have to be handed
    on explicitly or they apply to the harness and to nothing it runs - which is how a re-run of a
    failed comparison ended up refused by every variant in it.
    """
    flags = []
    if _runner.force:
        flags.append("--force")
    if _runner.resume:
        flags.append("--resume")
    return flags


def build_argv(
    comparison: Comparison,
    variant: Variant,
    inputs: Inputs,
    out_dir: str,
    threads: int,
    memory_gb: int,
    global_flags: Sequence[str] = (),
) -> List[str]:
    """The exact command a variant runs.

    Spawned as a process rather than imported and called, so that what the benchmark measures is a
    whole run - its env, its conda shim, its peak RSS - which is what a user's run is too.

    Addressed through ``python -m metawrap2.cli``, the same way ``metawrap2 run`` spawns a study's
    steps, because global flags are the CLI's and a module invoked as ``python -m
    metawrap2.modules.<name>`` has no way to receive them. Without ``--force`` reaching the variant,
    a comparison could not be re-run at all after a failure: every variant stopped at "already
    contains a MetaWrap2 run", and the only way forward was deleting five directories by hand.
    """
    from ..cli import MODULES

    if comparison.module not in MODULES:
        raise RunnerError(
            "%r is not a MetaWrap2 module. Known: %s"
            % (comparison.module, ", ".join(sorted(MODULES)))
        )
    argv = [sys.executable, "-m", "metawrap2.cli"] + list(global_flags) + [comparison.module]
    argv += list(comparison.common_args)
    argv += list(variant.args)
    argv += ["-o", out_dir, "-t", str(threads)]
    if comparison.module in ("assembly", "bin_refinement", "reassemble_bins"):
        argv += ["-m", str(memory_gb)]

    if comparison.module == "read_qc":
        # One sample per invocation, so the caller hands us a single pair.
        if len(inputs.pairs) != 1:
            raise RunnerError(
                "read_qc takes one sample at a time; got %d. This is a bug in the runner's "
                "per-sample loop, not something the registry can fix." % len(inputs.pairs)
            )
        _sample, r1, r2 = inputs.pairs[0]
        argv += ["-1", r1, "-2", r2]
    elif comparison.module == "assembly":
        # One pair, because `assembly` takes one file per mate - a comma-separated list is not a
        # form it accepts, and passing one made it exit immediately with "does not exist". A cohort
        # of samples is concatenated first, by prepare_coassembly, exactly as `metawrap2 run` does.
        if len(inputs.pairs) != 1:
            raise RunnerError(
                "assembly takes one pair of read files; got %d. The samples should have been "
                "concatenated first - this is a bug in the runner, not the registry."
                % len(inputs.pairs)
            )
        _sample, r1, r2 = inputs.pairs[0]
        argv += ["-1", r1, "-2", r2]
    elif comparison.module == "binning":
        argv += ["-a", str(inputs.assembly)]
        argv += inputs.read_files
    elif comparison.module == "quant_bins":
        # No -a: every quantifier under comparison works against the bins themselves, and giving
        # salmon alone a whole-assembly pseudo-transcriptome would compare three tools on two
        # different references.
        argv += ["-b", one_bin_set(str(inputs.bins_from))]
        if "depth" in variant.args:
            # The whole point of this quantifier is that the binning run already measured coverage,
            # so the arm has to be given that table rather than re-deriving one - which would be
            # measuring the aligner again under a different name.
            argv += ["--depth", depth_table_for(str(inputs.bins_from))]
        argv += inputs.read_files
    elif comparison.module == "bin_refinement":
        root = str(inputs.bins_from)
        sets = sorted(
            os.path.join(root, name)
            for name in os.listdir(root)
            if name.endswith("_bins") and os.path.isdir(os.path.join(root, name))
        )
        # A previous refinement into the same directory left its own answer behind as a
        # metawrap_<c>_<x>_bins folder. Refining a consolidated set together with the sets it was
        # consolidated from is not the comparison anyone declared, so it is dropped.
        sets = [path for path in sets if not is_refined_bins_dirname(path)]
        if not sets:
            raise RunnerError("Found no *_bins directories under %s to refine." % root)
        # --bins rather than -A/-B/-C: this used to zip the sets onto three flags, and zip
        # truncates. A variant that binned with six binners was silently refined from the three
        # whose directory names sort first - so metabat2_bins was dropped and the arm's published
        # numbers came from a different input set than the one it was named for.
        if len(sets) > MAX_REFINEMENT_SETS:
            raise RunnerError(
                "%s holds %d bin sets (%s) and bin_refinement consolidates at most %d, because it "
                "scores every combination of them with CheckM (2^n-1 runs). Point --bins-from at a "
                "variant with fewer binners, or refine in stages."
                % (
                    root,
                    len(sets),
                    ", ".join(os.path.basename(p) for p in sets),
                    MAX_REFINEMENT_SETS,
                )
            )
        for path in sets:
            argv += ["--bins", path]
        # The assembly too, when the dataset has one. MetaWrap2's own consolidation works from the
        # bins alone, but DAS Tool reads the assembly contigs itself and refuses to start without
        # them - so a consolidator comparison could never run its DAS Tool arm. Passed to every
        # variant rather than only the ones that ask, because an unused -a changes nothing and a
        # comparison whose arms are handed different inputs is not a comparison.
        if inputs.assembly:
            argv += ["-a", str(inputs.assembly)]
    return argv


def cost_rows(
    dataset: Dataset,
    comparison: Comparison,
    variant: Variant,
    measured: Dict[str, Any],
    sample: str = "",
) -> List[Row]:
    """The cost half of the table: what this variant spent.

    ``cpu_hours`` rather than wall time is the denominator of the ratio the defaults are chosen on,
    because wall time rewards a tool for being given more cores and CPU-hours do not. Wall time is
    recorded alongside it since it is what a user waits through.

    *sample* is set for a per-sample module, where one variant means several measured runs; it is
    empty for a cohort module, whose single run covers every sample at once.
    """

    def row(metric: str, value: float) -> Row:
        return Row(
            dataset=dataset.name,
            comparison=comparison.name,
            variant=variant.name,
            metric=metric,
            value=value,
            sample=sample,
            source="usage.run_measured",
        )

    rows = [row("wall_seconds", float(measured.get("wall_seconds") or 0.0))]
    user = measured.get("user_seconds")
    system = measured.get("system_seconds")
    if user is not None and system is not None:
        rows.append(row("cpu_hours", (float(user) + float(system)) / 3600.0))
    peak = measured.get("peak_rss_bytes")
    if peak is not None:
        rows.append(row("peak_rss_bytes", float(peak)))
    return rows


def _units(comparison: Comparison, inputs: Inputs) -> List[Sequence[Any]]:
    """The invocations one variant needs: ``(sample_label, inputs)`` pairs.

    A cohort module gets one invocation with every sample in it and an empty label; a per-sample
    module gets one per sample, labelled, so its rows sum the same way a cohort module's single
    measurement does.
    """
    if comparison.module not in PER_SAMPLE_MODULES:
        return [("", inputs)]
    units: List[Sequence[Any]] = []
    for pair in inputs.pairs:
        one = Inputs(pairs=[pair], assembly=inputs.assembly, bins_from=inputs.bins_from)
        units.append((str(pair[0]), one))
    return units


def gold_path(dataset, role: str, data_dir: str) -> Optional[str]:
    """Where *dataset*'s gold standard for *role* lives, or None if it declares none.

    The path is not checked here. Callers differ in what a missing file means - a binning run cannot
    start without the gold assembly, while scoring can still report cost without the gold binning -
    so each decides, and neither has to know the layout rule twice.
    """
    gold = dataset.gold(role)
    if gold is None:
        return None
    if dataset.is_local:
        return dataset.local_path_for(gold)
    return os.path.join(data_dir, "gold", gold.name)


def require_env_manager() -> None:
    """Refuse to start a comparison whose every variant would fail for the same reason.

    Each variant is a subprocess, so an absent ``mamba`` is not discovered until the first tool
    inside the first variant tries to start - and the runner's policy of recording a failure and
    carrying on then turns one environment problem into a directory of failed runs and an empty
    results table. A comparison over ten samples and two variants can spend a long time doing that.

    Called from the ``benchmark run`` command rather than from :func:`run_comparison`, which is
    driven in tests with the executor faked out and has no business requiring conda to be installed.
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return
    manager = _runner.env_manager
    if shutil.which(manager) is None:
        raise RunnerError(
            "%s is not on PATH, so no variant of this comparison could run its tools.\n"
            'Add it (usually PATH="$HOME/miniforge3/bin:$PATH") or set use_conda_envs = false '
            "under [settings] in ~/.metawrap2/config.toml." % manager
        )


def prepare_coassembly(
    inputs: Inputs,
    work_dir: str,
    *,
    dry_run: bool = False,
    say: Callable[[str], None] = print,
) -> Inputs:
    """One pair of read files holding every selected sample, for a cohort assembly.

    ``metawrap2 assembly`` takes one file per mate, so a co-assembly is fed a concatenation - which
    is what ``metawrap2 run`` builds for a study too, so the benchmark measures the same thing a
    user's run does rather than a form only the harness uses.

    Built once and reused by every variant, for the same reason the deinterleaved reads are: doing
    it per variant would hand the second one a colder page cache than the first, and that difference
    turns into a performance gap nobody meant to measure.
    """
    if len(inputs.pairs) == 1:
        return inputs
    dest = [os.path.join(work_dir, "coassembly_%d.fastq" % index) for index in (1, 2)]
    sources = [[pair[1] for pair in inputs.pairs], [pair[2] for pair in inputs.pairs]]
    for target, group in zip(dest, sources):
        if _is_current(target, group):
            say("  %s is already current" % os.path.basename(target))
            continue
        if dry_run:
            say("  (dry run) would concatenate %d file(s) into %s" % (len(group), target))
            continue
        say("  concatenating %d file(s) into %s" % (len(group), os.path.basename(target)))
        os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
        with open(target, "wb") as out:
            for src in group:
                with open(src, "rb") as fh:
                    shutil.copyfileobj(fh, out, length=8 << 20)
    return Inputs(
        pairs=[("coassembly", dest[0], dest[1])],
        assembly=inputs.assembly,
        bins_from=inputs.bins_from,
    )


def _is_current(dest: str, sources: Sequence[str]) -> bool:
    """True if *dest* exists, is non-empty, and is newer than every source."""
    if not os.path.isfile(dest) or not os.path.getsize(dest):
        return False
    mtime = os.path.getmtime(dest)
    return all(os.path.isfile(s) and os.path.getmtime(s) <= mtime for s in sources)


def select_samples(inputs: Inputs, *, samples: Sequence[str] = (), max_samples: int = 0) -> Inputs:
    """*inputs* narrowed to the named samples, or to the first *max_samples* of them.

    Naming a sample that is not there is an error listing what is, rather than a quiet run over
    nine samples when ten were meant: a typo in a sample label would otherwise change what the
    comparison measured without changing anything a reader of the results could see.
    """
    if not samples and max_samples <= 0:
        return inputs
    available = [str(pair[0]) for pair in inputs.pairs]
    if samples:
        missing = [name for name in samples if name not in available]
        if missing:
            raise RunnerError(
                "No such sample(s): %s\nThis dataset has: %s"
                % (", ".join(missing), ", ".join(available))
            )
        wanted = [pair for pair in inputs.pairs if str(pair[0]) in set(samples)]
    else:
        wanted = list(inputs.pairs)
    if max_samples > 0:
        wanted = wanted[:max_samples]
    if not wanted:
        raise RunnerError("The sample selection left nothing to run.")
    return Inputs(pairs=wanted, assembly=inputs.assembly, bins_from=inputs.bins_from)


def run_comparison(
    dataset: Dataset,
    comparison: Comparison,
    data_dir: str,
    work_dir: str,
    *,
    variants: Sequence[str] = (),
    samples: Sequence[str] = (),
    max_samples: int = 0,
    assembly: Optional[str] = None,
    bins_from: Optional[str] = None,
    budget: Optional[_limits.Budget] = None,
    dry_run: bool = False,
    say: Callable[[str], None] = print,
) -> List[Row]:
    """Run every (or the named) variant and return the cost rows.

    Accuracy rows come later, from :mod:`metawrap2.benchmark.scoring` reading the tool output that
    lands under *work_dir*; keeping the two apart means a scorer can be improved and re-applied to
    a run that already happened rather than requiring the run again.

    *samples* and *max_samples* narrow which of the dataset's samples take part, because some
    comparisons are not affordable over all of them - a 10-sample CAMI2 co-assembly through
    metaSPAdes is days and hundreds of GB of RAM, and the question of which assembler is more
    contiguous does not need every sample to answer. The subset applies to every variant equally,
    so the comparison stays internally valid; what it is not is comparable with a *different*
    subset, which is why the chosen samples are written into the run record and printed.
    """
    budget = budget or _limits.resolve()
    chosen = (
        [comparison.variant(name) for name in variants] if variants else list(comparison.variants)
    )
    inputs = gather_inputs(
        dataset,
        comparison,
        data_dir,
        work_dir,
        assembly=assembly,
        bins_from=bins_from,
        dry_run=dry_run,
        say=say,
        wanted=samples,
        max_samples=max_samples,
    )
    # Still narrowed here, not only in gather_inputs: that one decides which files to write, this
    # one is the single place that decides which samples the comparison measured, and it is what
    # errors on a sample name the dataset does not have.
    inputs = select_samples(inputs, samples=samples, max_samples=max_samples)
    measured_samples = [str(pair[0]) for pair in inputs.pairs]
    if comparison.module == "assembly":
        inputs = prepare_coassembly(inputs, work_dir, dry_run=dry_run, say=say)

    say("")
    say("Comparison: %s (%s)" % (comparison.name, comparison.question or comparison.module))
    say("Budget:     %s" % budget.describe())
    say("Variants:   %s" % ", ".join(v.name for v in chosen))
    say("Samples:    %d  (%s)" % (len(measured_samples), ", ".join(measured_samples)))

    rows: List[Row] = []
    record: Dict[str, Any] = {
        "dataset": dataset.name,
        "comparison": comparison.name,
        "budget": budget.as_dict(),
        # What the numbers cover. A cohort module's rows carry no sample label, so without this a
        # three-sample co-assembly and a ten-sample one are indistinguishable in the results TSV.
        "samples": measured_samples,
        # Where the dataset was read from, so `score` can find the same gold standard the run used
        # without being told again - and so a record says which copy of the data produced it.
        "data_dir": os.path.abspath(data_dir),
        # The contigs the binners were actually given. `score` needs them for AMBER's sequence
        # lengths, and --assembly means the dataset's default is not always the right guess.
        "assembly": os.path.abspath(inputs.assembly) if inputs.assembly else None,
        # The bins the variants were given, for the comparisons that take a previous run's output.
        # A quantifier's truth is built from these, so scoring has to read the same set the run did
        # rather than whichever bin directory happens to be lying around at scoring time.
        "bins_from": os.path.abspath(inputs.bins_from) if inputs.bins_from else None,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "variants": {},
    }
    for variant in chosen:
        variant_record: Dict[str, Any] = {"returncode": 0, "runs": []}
        record["variants"][variant.name] = variant_record
        variant_rows: List[Row] = []
        say("")
        say("--- %s ---" % variant.name)
        for label, unit in _units(comparison, inputs):
            out_dir = os.path.join(
                work_dir, comparison.name, variant.name, *([label] if label else [])
            )
            argv = build_argv(
                comparison,
                variant,
                unit,
                out_dir,
                budget.threads,
                int(budget.memory_gb),
                global_flags=_inherited_global_flags(),
            )
            if label:
                say("    [%s]" % label)
            say("    %s" % " ".join(argv))
            if dry_run:
                continue
            os.makedirs(out_dir, exist_ok=True)
            rc, measured = _usage.run_measured(argv)
            say("    %s  (exit %d)" % (_usage.summarise(measured), rc))
            variant_record["runs"].append(
                {
                    "sample": label,
                    "command": argv,
                    "returncode": rc,
                    "usage": measured,
                    "output": out_dir,
                }
            )
            if rc != 0:
                # Recorded, not raised: one variant failing is a result about that variant, and
                # the others are still worth having. A failed variant contributes no rows at all -
                # not even the samples that did finish - because a ratio summed over nine samples
                # is not comparable with one summed over ten, and the partial number would look
                # like a cheap win.
                variant_record["returncode"] = rc
                say("    FAILED - recorded, abandoning this variant and continuing with the rest")
                variant_rows = []
                break
            variant_rows.extend(cost_rows(dataset, comparison, variant, measured, sample=label))
        rows.extend(variant_rows)

    if not dry_run:
        os.makedirs(work_dir, exist_ok=True)
        path = os.path.join(work_dir, "%s.run.json" % comparison.name)
        with open(path, "w") as fh:
            json.dump(record, fh, indent=2, sort_keys=True)
        say("")
        say("Run record: %s" % path)
    return rows


def load_registry_error_hint(exc: RegistryError) -> str:
    """A registry problem, phrased for someone editing ``benchmarks/datasets.toml``."""
    return (
        "%s\n\nThe registry is benchmarks/datasets.toml; it is data, so this is an edit away." % exc
    )
