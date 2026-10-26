"""Run one comparison's variants over one dataset under one budget, and record the cost.

Every variant of a comparison gets the same inputs, the same thread and memory grant and the same
machine; only the flags under test differ.

* The budget is resolved once and given to every variant, so a variant that starts while other
  work is finishing gets the same grant as the rest.
* Binning comparisons default to the gold-standard assembly, which CAMI2 provides so binners can
  be measured without assembler variance.
* :func:`metawrap2.usage.run_measured` reaps each child with ``wait4``, so peak RSS and CPU seconds
  are per step, not a running maximum over every child the process has waited for.
* A variant can be another pipeline (:class:`metawrap2.benchmark.registry.External`). It gets the
  same inputs, budget and measurement. Only its argv and the publishing of its bins under the name
  the scorer reads are different.

This module produces cost rows; :func:`metawrap2.benchmark.scoring.value_against_cost` sets them
beside the accuracy results.
"""

from __future__ import annotations

import glob
import json
import os
import shlex
import shutil
import sys
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Sequence

from .. import config as _config
from .. import limits as _limits
from .. import usage as _usage
from ..command import capture_in_env as _capture_in_env
from ..command import runner as _runner
from ..constants import BIN_EXTENSION, bin_label, is_refined_bins_dirname
from ..io import reads as _reads
from ..utils import strip_trailing_slash
from .registry import PLACEHOLDER_RE, Comparison, Dataset, RegistryError, Variant
from .scoring import Row

#: How many bin sets bin_refinement will consolidate: its -A..-F flags, and the 2^n-1 CheckM runs
#: behind them (63 for six sets).
MAX_REFINEMENT_SETS = 6

#: Modules the runner can feed, and what each needs beyond reads. An explicit table, so an
#: unsupported module is an error instead of a guessed command line.
SUPPORTED_MODULES = {
    "read_qc": "reads, one sample at a time",
    "assembly": "reads",
    "binning": "reads + an assembly (the gold standard one by default)",
    "bin_refinement": "bin sets from a previous binning run (--bins-from)",
    "quant_bins": "reads + one set of bins from a previous run (--bins-from)",
}

#: Where an external arm writes its own output, inside the arm's directory. The scorer looks for
#: bin sets by name directly under the arm, so the foreign pipeline's directories are kept one
#: level down, away from the bin set the runner publishes.
EXTERNAL_RUN_DIRNAME = "external"

#: Modules that take one sample per invocation instead of a whole cohort. They run once per sample
#: and their cost rows carry the sample label; the rows sum into one ratio, as for a cohort module.
PER_SAMPLE_MODULES = frozenset({"read_qc"})


class RunnerError(Exception):
    """A comparison that cannot be run as asked."""


@dataclass
class Inputs:
    """The prepared inputs a dataset offers, as found on disk."""

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

    Sorted so two runs pass reads in the same order and their outputs are byte-comparable.
    """
    found: List[str] = []
    for base, _dirs, names in os.walk(root):
        for name in sorted(names):
            if name.endswith(tuple(suffixes)):
                found.append(os.path.join(base, name))
    return sorted(found)


def read_candidates(dataset: Dataset, source: str) -> List[str]:
    """The dataset's read files: its declared ``reads_glob`` if it has one, else every FASTQ.

    CAMI2's ``reads/`` directory holds the 10 GB interleaved FASTQ and a 3 MB ``subsample_0.fq``
    debugging extract under the same ``sample_N`` directory. A plain suffix scan finds both under
    one sample label, so datasets like this declare a glob.
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

    Keyed by :func:`_sample_label`, the label the quantifier's table columns carry, so estimate
    and truth are joined on labels from the same function.
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

    Unlike :func:`_sample_label`, which labels a file, this strips the mate suffix so both files
    of a paired sample get one label.

    The suffix is stripped only for a paired layout. Interleaved CAMI2 labels are ``sample_0``,
    ``sample_1``, ``sample_2``, and stripping a trailing ``_1`` would merge two samples.
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

    Named samples that are absent are not an error here; :func:`select_samples` reports them
    against the dataset's full sample list.
    """
    if not wanted and max_samples <= 0:
        return list(candidates)
    # A list of paths per sample: a paired dataset has two files under one label, and both mates
    # must be kept.
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
    # Original order, so a dry run's printed plan is stable between runs. The downstream pairing
    # matches mates by name.
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

    *wanted* and *max_samples* are the narrowing :func:`select_samples` applies, passed here so
    only the selected samples are split. (Splitting all ten CAMI2 samples for a three-sample
    comparison wrote 190 GB to use 57 GB.) :func:`select_samples` still makes the final selection.

    Splitting is cached and shared by all variants, so every variant reads from the same files
    with a similar page cache.

    Under *dry_run* the paths are named but nothing is written (for CAMI2, about 100 GB of mate
    files).
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
                # Backstop for the glob: two files under one label would deinterleave to the same
                # output, and the second would overwrite the first.
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
        # collect_pairs' message lists the naming conventions it recognises; re-raised as a
        # RunnerError so the benchmark's errors share one class.
        raise RunnerError("Could not pair the reads under %s.\n%s" % (source, exc))
    return sorted(triples, key=lambda triple: triple[0])


def _sample_label(path: str) -> str:
    """A stable sample label from a CAMI2-style path.

    The published directories are timestamp-prefixed (``2017.12.04_18.45.54_sample_0``) with a
    different timestamp per sample, so the label is taken from the ``sample_N`` part.
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

    *wanted* and *max_samples* are passed to :func:`prepare_reads` so only the selected samples
    are prepared.
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
    # bin_refinement needs the assembly too: DAS Tool reads the contigs itself and will not start
    # without them. build_argv passes -a only if this field is filled in here.
    if comparison.module in ("binning", "bin_refinement") and resolved_assembly is None:
        candidate = gold_path(dataset, "assembly", data_dir)
        if candidate is None:
            if comparison.module == "bin_refinement":
                # Only DAS Tool needs it; the default consolidator works from the bins alone, so
                # the other arms can still run.
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
                "  binning on the gold-standard assembly (%s), so assembler variation does not\n"
                "  affect the comparison. Pass --assembly to compare on a real assembly instead."
                % os.path.basename(candidate)
            )
        else:
            say(
                "  passing the gold-standard assembly (%s) to every arm, since DAS Tool reads the\n"
                "  contigs itself and every arm needs the same inputs."
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

    ``bin_refinement`` takes several bin sets; ``quant_bins`` takes one. *root* may be a folder of
    FASTA files, or a variant's output directory holding several ``*_bins`` subdirectories (what
    a binning comparison leaves).

    With several sets, the consolidated set is used. Otherwise RunnerError is raised instead of
    picking one of several overlapping raw bin sets.
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
            "%s holds %d consolidated bin sets (%s), usually from a refinement run resumed at "
            "different -c/-x thresholds. Point --bins-from at the one to use."
            % (root, len(consolidated), ", ".join(os.path.basename(p) for p in consolidated))
        )
    if consolidated:
        return consolidated[0]
    if len(sets) == 1:
        return sets[0]
    if not sets:
        raise RunnerError(
            "Found no bins and no *_bins directory under %s. --bins-from takes a folder of bins "
            "or the output directory of a binning run." % root
        )
    raise RunnerError(
        "%s holds %d bin sets (%s) and no consolidated metawrap_<c>_<x>_bins. The sets overlap. "
        "Refine them into one set first, "
        "or point --bins-from at the set to quantify."
        % (root, len(sets), ", ".join(os.path.basename(p) for p in sets))
    )


def depth_table_for(root: str) -> str:
    """The contig depth table the binning run under *root* left behind.

    ``quant_bins --quantifier depth`` reuses binning's table instead of aligning again. Searched
    beside the bins and one level up, covering ``--bins-from`` pointing at a binning output or at
    a refinement output inside it.
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
        "exist:\n  %s\n`metawrap2 binning` writes that table. Point --bins-from at a binning "
        "output directory (or the refinement inside one), not at a folder of bins "
        "alone." % "\n  ".join(candidates)
    )


def _inherited_global_flags() -> List[str]:
    """``--force``, ``--resume`` and ``--gpu-strict`` as typed on the benchmark's command line.

    Each variant is a separate process, so these flags are passed on explicitly.
    """
    flags = []
    if _runner.force:
        flags.append("--force")
    if _runner.resume:
        flags.append("--resume")
    if _limits.overrides().get("gpu_strict"):
        flags.append("--gpu-strict")
    return flags


def external_output_dir(out_dir: str) -> str:
    """Where an external arm's own output goes: ``<arm>/external``."""
    return os.path.join(out_dir, EXTERNAL_RUN_DIRNAME)


def external_bins_dir(variant: Variant, out_dir: str) -> str:
    """The directory the arm declared its bins land in, resolved under *out_dir*."""
    if variant.external is None:
        raise RunnerError("%s is not an external arm." % variant.name)
    return os.path.join(external_output_dir(out_dir), variant.external.bins)


def published_bins_dir(variant: Variant, out_dir: str) -> str:
    """Where the arm's bins are presented for scoring: ``<arm>/<arm name>_bins``.

    :func:`metawrap2.benchmark.assess.bin_dirs` finds bin sets by the ``_bins`` suffix, so an
    external arm is scored by the same code as a native one.
    """
    return os.path.join(out_dir, "%s_bins" % variant.name)


def external_env(variant: Variant) -> Optional[str]:
    """The conda env an external arm runs in, or None when the config says not to use conda envs.

    Resolved at call time, as :func:`assess.amber_env` is, so an edit to the config file takes
    effect.
    """
    if variant.external is None:
        return None
    if not _config.load_settings().use_conda_envs:
        return None
    return variant.external.env


def external_values(
    inputs: Inputs,
    out_dir: str,
    threads: int,
    memory_gb: int,
    sample: str = "",
) -> Dict[str, str]:
    """What each placeholder resolves to for this invocation.

    Only the values this run has. :func:`external_argv` refuses a placeholder missing here by
    name, so a pipeline is never given the string "None" as a path.
    """
    values = {
        "output": external_output_dir(out_dir),
        "threads": str(threads),
        "memory": str(memory_gb),
    }
    if inputs.assembly:
        values["assembly"] = inputs.assembly
    if inputs.read_files:
        values["reads"] = " ".join(inputs.read_files)
    if len(inputs.pairs) == 1:
        # Only for one pair. A per-sample module is called once per sample, so {r1} is defined
        # there. Filling it from the first of several pairs would run one sample and label the row
        # with the whole dataset.
        _sample, r1, r2 = inputs.pairs[0]
        values["r1"] = r1
        values["r2"] = r2
    if inputs.bins_from:
        values["bins"] = one_bin_set(str(inputs.bins_from))
    if sample:
        values["sample"] = sample
    return values


def external_argv(
    variant: Variant,
    inputs: Inputs,
    out_dir: str,
    threads: int,
    memory_gb: int,
    sample: str = "",
) -> List[str]:
    """The exact command an external arm runs, inside its own conda env.

    The template is substituted first and then split with :mod:`shlex`, because ``{reads}`` is a
    list of files and becomes several argv entries. A pipeline with two steps declares them as its
    own ``bash -c '... && ...'``, which is still one measured process tree.
    """
    external = variant.external
    if external is None:
        raise RunnerError("%s is not an external arm." % variant.name)
    values = external_values(inputs, out_dir, threads, memory_gb, sample=sample)
    missing = [name for name in external.placeholders if name not in values]
    if missing:
        raise RunnerError(
            "External arm %r uses %s, which this run has nothing to fill from.\n"
            "The registry accepted the placeholder, so this is about the run rather than the "
            "declaration: an {assembly} needs a gold-standard assembly or --assembly, {bins} needs "
            "--bins-from, and {r1}/{r2} need a module that takes one sample at a time."
            % (variant.name, ", ".join("{%s}" % name for name in missing))
        )
    resolved = PLACEHOLDER_RE.sub(lambda m: values[m.group(1)], external.command)
    return _runner.env_prefix(external_env(variant)) + shlex.split(resolved)


def external_version(variant: Variant, *, probe: bool = True) -> str:
    """What version an external arm ran, for the run record.

    A declared ``version`` is used as given. Otherwise the arm's ``version_command`` is run in its
    env, so the version comes from the installed binary.

    A version that cannot be read is an error naming the arm, so every results row can be
    reproduced from the table alone. *probe* is False under a dry run, where the env may not exist
    yet.
    """
    external = variant.external
    if external is None:
        raise RunnerError("%s is not an external arm." % variant.name)
    if external.version:
        return external.version
    if not probe:
        return ""
    output = _capture_in_env(external_env(variant), shlex.split(external.version_command))
    for line in output.splitlines():
        if line.strip():
            return line.strip()
    raise RunnerError(
        "External arm %r: `%s` in env %r said nothing, so the version it would run is unknown.\n"
        "A results row claiming a pipeline without a version is not reproducible. Either fix that "
        'command or pin `version = "..."` in benchmarks/configurations.toml.'
        % (variant.name, external.version_command, external.env)
    )


def publish_external_bins(variant: Variant, out_dir: str) -> str:
    """Present an external arm's bins as ``<arm name>_bins``, and return that directory.

    One symlink per file, named with the bin extension the scorer looks for:
    :func:`metawrap2.checkm.count_bins` counts ``.fasta``, and metaWRAP 1.3.0 writes ``bin.1.fa``.
    A link to the whole directory would show no bins.

    The pipeline's own output is not changed, and nothing is copied.
    """
    from ..io.seqio import fasta_files

    source = external_bins_dir(variant, out_dir)
    if not os.path.isdir(source):
        raise RunnerError(
            "External arm %r exited 0 but wrote no bins at %s, which is where its declaration says "
            "they land.\nEither the pipeline failed quietly or `bins` in "
            "benchmarks/configurations.toml names the wrong path; both would otherwise be an empty "
            "row in the results table." % (variant.name, source)
        )
    bins = fasta_files(source)
    if not bins:
        raise RunnerError(
            "External arm %r exited 0 and %s holds no FASTA, so there is nothing to score."
            % (variant.name, source)
        )
    gzipped = [path for path in bins if path.endswith(".gz")]
    if gzipped:
        raise RunnerError(
            "External arm %r wrote %d gzipped bin(s) (%s). CheckM1 and prokka read uncompressed "
            "bins only, so the arm has to decompress them itself - have its command end with a "
            "`gunzip` over %s." % (variant.name, len(gzipped), os.path.basename(gzipped[0]), source)
        )
    published = published_bins_dir(variant, out_dir)
    if os.path.islink(published):
        os.unlink(published)
    elif os.path.isdir(published):
        # A previous publish is a directory of symlinks and is rebuilt. Real files under this name
        # are other data and are not deleted.
        real = [
            name
            for name in os.listdir(published)
            if not os.path.islink(os.path.join(published, name))
        ]
        if real:
            raise RunnerError(
                "%s already exists and holds %d file(s) that are not links to an external arm's "
                "bins (%s). Refusing to replace it - rename it, or rename the arm."
                % (published, len(real), ", ".join(sorted(real)[:3]))
            )
        shutil.rmtree(published)
    os.makedirs(published)
    for path in bins:
        link = os.path.join(published, bin_label(os.path.basename(path)) + BIN_EXTENSION)
        # Relative, so the published set still resolves after the work directory is moved or
        # copied.
        os.symlink(os.path.relpath(path, published), link)
    return published


def build_argv(
    comparison: Comparison,
    variant: Variant,
    inputs: Inputs,
    out_dir: str,
    threads: int,
    memory_gb: int,
    global_flags: Sequence[str] = (),
    sample: str = "",
) -> List[str]:
    """The exact command a variant runs.

    Spawned as a process, so the measurement covers a whole run (env, conda shim, peak RSS) as a
    user would run it.

    Invoked through ``python -m metawrap2.cli``, as ``metawrap2 run`` spawns a study's steps,
    because global flags such as ``--force`` belong to the CLI and ``python -m
    metawrap2.modules.<name>`` cannot receive them.

    An external arm is built by :func:`external_argv`, so this one function answers what any arm
    runs, and ``--dry-run`` prints foreign arms beside native ones.
    """
    from ..cli import MODULES

    if variant.is_external:
        return external_argv(variant, inputs, out_dir, threads, memory_gb, sample=sample)
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
        # One sample per invocation, so the caller passes a single pair.
        if len(inputs.pairs) != 1:
            raise RunnerError(
                "read_qc takes one sample at a time; got %d. This is a bug in the runner's "
                "per-sample loop." % len(inputs.pairs)
            )
        _sample, r1, r2 = inputs.pairs[0]
        argv += ["-1", r1, "-2", r2]
    elif comparison.module == "assembly":
        # One pair: `assembly` takes one file per mate, not a comma-separated list. A cohort is
        # concatenated first by prepare_coassembly, as `metawrap2 run` does.
        if len(inputs.pairs) != 1:
            raise RunnerError(
                "assembly takes one pair of read files; got %d. The samples should have been "
                "concatenated first. This is a bug in the runner." % len(inputs.pairs)
            )
        _sample, r1, r2 = inputs.pairs[0]
        argv += ["-1", r1, "-2", r2]
    elif comparison.module == "binning":
        argv += ["-a", str(inputs.assembly)]
        argv += inputs.read_files
    elif comparison.module == "quant_bins":
        # No -a: every quantifier under comparison works against the bins, so all use the same
        # reference (salmon would otherwise index the whole assembly).
        argv += ["-b", one_bin_set(str(inputs.bins_from))]
        if "depth" in variant.args:
            # This quantifier reuses the coverage the binning run measured, so pass that table.
            argv += ["--depth", depth_table_for(str(inputs.bins_from))]
        argv += inputs.read_files
    elif comparison.module == "bin_refinement":
        root = str(inputs.bins_from)
        sets = sorted(
            os.path.join(root, name)
            for name in os.listdir(root)
            if name.endswith("_bins") and os.path.isdir(os.path.join(root, name))
        )
        # Drop a metawrap_<c>_<x>_bins folder left by an earlier refinement into this directory, so
        # a consolidated set is not refined together with its own inputs.
        sets = [path for path in sets if not is_refined_bins_dirname(path)]
        if not sets:
            raise RunnerError("Found no *_bins directories under %s to refine." % root)
        # --bins instead of -A/-B/-C. Zipping the sets onto three flags truncated a six-binner
        # variant to the three sets that sorted first.
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
        # The assembly too, when there is one: DAS Tool reads the contigs itself and will not start
        # without them. Passed to every variant so all arms get the same inputs; MetaWrap2's own
        # consolidation ignores it.
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
    """The cost rows for this variant.

    ``cpu_hours`` is the denominator of the ratio defaults are chosen on, because wall time
    favours tools given more cores. Wall time is recorded too.

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

    A cohort module gets one invocation with every sample and an empty label; a per-sample module
    gets one labelled invocation per sample.
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

    The path is not checked here; callers decide what a missing file means (binning cannot start
    without the gold assembly, while scoring can still report cost without the gold binning).
    """
    gold = dataset.gold(role)
    if gold is None:
        return None
    if dataset.is_local:
        return dataset.local_path_for(gold)
    return os.path.join(data_dir, "gold", gold.name)


def require_env_manager() -> None:
    """Raise RunnerError if the conda env manager is missing, before any variant starts.

    Each variant is a subprocess, so a missing ``mamba`` would otherwise be found only when the
    first tool starts, and every variant would be recorded as failed.

    Called from the ``benchmark run`` command, not :func:`run_comparison`, which tests drive with
    a fake executor and no conda.
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return
    manager = _runner.env_manager
    if shutil.which(manager) is None:
        raise RunnerError(
            "%s is not on PATH, so the tools for this comparison cannot run.\n"
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

    ``metawrap2 assembly`` takes one file per mate, so a co-assembly gets a concatenation, as
    ``metawrap2 run`` builds for a study.

    Built once and reused by every variant, like the deinterleaved reads.
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

    Naming a sample that does not exist raises RunnerError listing the available samples, so a
    typo does not change what the comparison measures.
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

    Accuracy rows come later, from :mod:`metawrap2.benchmark.scoring` reading the tool output
    under *work_dir*, so a scorer can be changed and re-applied without re-running.

    *samples* and *max_samples* narrow which samples take part, since some comparisons are too
    expensive over all of them (a 10-sample CAMI2 co-assembly through metaSPAdes takes days and
    hundreds of GB of RAM). The subset applies to every variant. Results are not comparable across
    different subsets, so the chosen samples are printed and written into the run record.
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
    # gather_inputs only decides which files to write; select_samples decides which samples are
    # measured and rejects unknown sample names.
    inputs = select_samples(inputs, samples=samples, max_samples=max_samples)
    measured_samples = [str(pair[0]) for pair in inputs.pairs]
    if comparison.module == "assembly":
        inputs = prepare_coassembly(inputs, work_dir, dry_run=dry_run, say=say)

    say("")
    say("Comparison: %s (%s)" % (comparison.name, comparison.question or comparison.module))
    say("Budget:     %s" % budget.describe())
    say("Variants:   %s" % ", ".join(v.name for v in chosen))
    say("Samples:    %d  (%s)" % (len(measured_samples), ", ".join(measured_samples)))

    # Every external arm's version, resolved before anything runs (as require_env_manager is), so
    # a comparison that cannot record a version stops before any compute is spent.
    external_versions = {
        variant.name: external_version(variant, probe=not dry_run)
        for variant in chosen
        if variant.is_external
    }
    for name, version in sorted(external_versions.items()):
        say("External:   %s is %s" % (name, version or "an unprobed version (dry run)"))

    rows: List[Row] = []
    record: Dict[str, Any] = {
        "dataset": dataset.name,
        "comparison": comparison.name,
        "budget": budget.as_dict(),
        # The samples measured. A cohort module's rows carry no sample label, so this is the only
        # record of how many samples a co-assembly used.
        "samples": measured_samples,
        # Where the dataset was read from, so `score` finds the same gold standard the run used.
        "data_dir": os.path.abspath(data_dir),
        # The contigs the binners were given (possibly from --assembly). `score` needs them for
        # AMBER's sequence lengths.
        "assembly": os.path.abspath(inputs.assembly) if inputs.assembly else None,
        # The bins the variants were given, for comparisons that take a previous run's output.
        # A quantifier's truth is built from these, so scoring must read the same set.
        "bins_from": os.path.abspath(inputs.bins_from) if inputs.bins_from else None,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "variants": {},
    }
    for variant in chosen:
        variant_record: Dict[str, Any] = {"returncode": 0, "runs": []}
        external = variant.external
        if external is not None:
            # Provenance for the arm, in the record `score` reads: the version that ran, its env
            # and the command template, so the row can be reproduced.
            variant_record["external"] = True
            variant_record["external_version"] = external_versions[variant.name]
            variant_record["external_env"] = external.env
            variant_record["external_command_template"] = external.command
            variant_record["external_bins"] = external.bins
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
                sample=label,
            )
            if label:
                say("    [%s]" % label)
            # Shell-quoted, so the printed command can be copied and re-run. An external arm's
            # `bash -c '...'` is one argument.
            say("    %s" % shlex.join(argv))
            if dry_run:
                continue
            os.makedirs(out_dir, exist_ok=True)
            if variant.is_external:
                os.makedirs(external_output_dir(out_dir), exist_ok=True)
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
                # Recorded, not raised, so the other variants still run. A failed variant
                # contributes no rows, not even for finished samples, because a ratio summed over
                # fewer samples is not comparable with the others.
                variant_record["returncode"] = rc
                say("    FAILED - recorded; skipping the rest of this variant and continuing")
                variant_rows = []
                break
            if variant.is_external:
                try:
                    published = publish_external_bins(variant, out_dir)
                except RunnerError as exc:
                    # Recorded as a failed arm, like a non-zero exit, so the run record keeps the
                    # arms that already finished. The arm contributes no rows.
                    variant_record["returncode"] = 1
                    variant_record["error"] = str(exc)
                    say("    %s" % exc)
                    say(
                        "    FAILED - recorded, abandoning this variant and continuing with the rest"
                    )
                    variant_rows = []
                    break
                variant_record["runs"][-1]["published_bins"] = published
                say("    bins published as %s" % os.path.basename(published))
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
    return "%s\n\nFix the entry in the registry, benchmarks/datasets.toml." % exc
