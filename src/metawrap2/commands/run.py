"""`metawrap2 run <samples.yaml|toml> -o study/` - drive the whole pipeline from one file.

Without it, a study means chaining the modules by hand: creating intermediate directories,
moving ``final_pure_reads_*`` into ``CLEAN_READS``, concatenating reads for a co-assembly and
passing each bin directory to the next module. None of that is recorded in provenance or
resumable.

This command does that work. It takes the inputs, sets the layout, runs the steps in order with
independent samples in parallel, and skips steps that are already done.

    metawrap2 run samples.toml -o study/ -t 24
    metawrap2 run samples.toml -o study/ --dry-run     # print the plan, run nothing
    metawrap2 run samples.toml -o study/ --resume      # continue where it stopped

Each step is an ordinary ``metawrap2 <module>`` invocation built by :func:`plan_steps`, so every
step can also be run by hand; the commands are in the study's ``run_commands.txt``. Step state is
kept in ``<output>/.metawrap2/run_steps``, and ``--resume`` uses it to tell which steps finished.

The sample sheet is TOML, or YAML if PyYAML is installed:

    [settings]
    coassemble = true
    steps = ["read_qc", "assembly", "binning", "bin_refinement"]

    [[samples]]
    name = "ERR011347"
    r1 = "raw/ERR011347_1.fastq.gz"
    r2 = "raw/ERR011347_2.fastq.gz"

    [[samples]]
    name = "ERR011348"
    r1 = "raw/ERR011348_1.fastq.gz"
    r2 = "raw/ERR011348_2.fastq.gz"

``r2`` may be omitted for single-end data (set ``layout = "single"``) or interleaved data
(``layout = "interleaved"``). If you give a directory instead of a sample sheet, the reads in
it are paired up automatically using the conventions in :mod:`metawrap2.io.reads`.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import limits as _limits
from .. import usage as _usage
from ..cliformat import make_parser
from ..config import load_settings
from ..io import reads as _reads
from ..logging import announcement, comm, error, warning
from ..manifest import COMPLETED, FAILED, INTERRUPTED, MANIFEST_NAME, Manifest
from ..modules import _common
from ..modules import binning as _binning
from ..modules import coverage as _coverage
from ..modules import dereplicate as _dereplicate

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml

#: Steps in pipeline order. Each is a module name; the driver knows how to wire each one.
ALL_STEPS: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "kraken2",
    "coverage",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "dereplicate",
    "quant_bins",
    "classify_bins",
    "annotate_bins",
    "blobology",
)

#: Steps that apply only when each sample is assembled separately: `coverage` builds one depth
#: table across the assemblies and `dereplicate` merges the per-sample bin sets into one catalogue.
COHORT_ONLY_STEPS = frozenset({"coverage", "dereplicate"})

#: What a cohort loses when each of those steps is left out of a pinned `steps` list.
COHORT_ONLY_LOSSES: Dict[str, str] = {
    "coverage": "one depth table across every sample",
    "dereplicate": "a non-redundant MAG catalogue",
}

#: Default steps: the genome-recovery path. kraken2, classify_bins and blobology need the largest
#: databases (~300 GB BLAST, ~90 GB Kraken2) and are opt-in via `steps`.
DEFAULT_STEPS: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "quant_bins",
)

#: Default cohort steps: the same path plus the two cohort-only stages. `reassemble_bins` is left
#: out because it is the most expensive step and would run once per sample; add it to `steps` if
#: wanted.
COHORT_STEPS: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "coverage",
    "binning",
    "bin_refinement",
    "dereplicate",
    "quant_bins",
)

#: Binner flags used when the sheet names none, taken from the binning module's default tier so the
#: two cannot diverge.
#:
#: CONCOCT was removed from the default on 2026-10-02 for cost. On cami2_gi it used 192 of the 213
#: CPU-hours of a three-binner run, for 10 extra GUNC-clean MAGs (74 against 64); metaBAT2+SemiBin2
#: reach 105 for much less. It is still available as --concoct. This pair is also the most chimeric
#: of the three consolidated sets measured (18.9%); see benchmarks/configurations.toml.
DEFAULT_BINNERS: Tuple[str, ...] = tuple(
    "--" + binner for binner in _binning.PRESETS[_binning.DEFAULT_PRESET]
)

#: Modules that accept -m/--memory in GB. The driver divides its budget between concurrent
#: steps and passes each one its share.
MEMORY_MODULES = frozenset({"assembly", "binning", "bin_refinement", "reassemble_bins"})

#: Per-sample steps; everything else operates on the study as a whole.
PER_SAMPLE_STEPS = frozenset({"read_qc"})

#: Fraction of its thread grant below which a step's use suggests a larger -j. It is the same
#: threshold as the mid-step warning in metawrap2.usage, so the two never disagree about what counts
#: as idle. read_qc on the three-sample tutorial cohort kept 1.2-2.0 cores busy of a large grant,
#: because its parts are largely serial (bmtagger is single-threaded, FastQC parallelises per file,
#: fastp stops scaling around 16 threads).
IDLE_GRANT_FRACTION = _usage.LOW_UTILISATION_FRACTION

#: The two study strategies, recorded in the manifest. `run` co-assembles all reads and bins the
#: result; `cohort` assembles each sample and dereplicates the per-sample bin sets. Their catalogues
#: are not comparable (genome counts, redundancy and per-genome coverage differ).
CO_ASSEMBLY = "co-assembly"
PER_SAMPLE = "per-sample"

#: Top-level key the strategy is recorded under in the study manifest.
STRATEGY_KEY = "study"

_STATE_DIR = os.path.join(".metawrap2", "run_steps")


# --- the sample sheet ---------------------------------------------------------------------


@dataclass
class Sample:
    name: str
    r1: str
    r2: Optional[str] = None
    layout: str = "paired"  # "paired" | "single" | "interleaved"

    @property
    def files(self) -> List[str]:
        return [f for f in (self.r1, self.r2) if f]

    def layout_flags(self) -> List[str]:
        if self.layout == "single":
            return ["--single-end"]
        if self.layout == "interleaved":
            return ["--interleaved"]
        return []


@dataclass
class Study:
    samples: List[Sample]
    coassemble: bool = True
    steps: Tuple[str, ...] = DEFAULT_STEPS
    module_options: Dict[str, List[str]] = field(default_factory=dict)
    #: True when the sheet set `coassemble` either way, as opposed to taking the default.
    coassemble_chosen: bool = False
    #: True when the sheet named the steps. The default depends on the assembly mode, which can
    #: still change after the sheet is read (`--no-coassembly`).
    steps_chosen: bool = False

    def extra(self, module: str) -> List[str]:
        """Extra command-line options the sheet specified for *module*."""
        return list(self.module_options.get(module, []))

    def default_steps(self) -> Tuple[str, ...]:
        """Return the default steps for the assembly mode."""
        return DEFAULT_STEPS if self.coassemble else COHORT_STEPS


def resolve_steps(study: Study, mode_final: bool = True) -> None:
    """Set *study*'s step list for its assembly mode.

    Called again after command-line parsing, since `--no-coassembly` changes the default. It is
    idempotent for explicitly named steps.

    *mode_final* is False when called from the sheet reader, where the mode may still change (for
    example under `metawrap2 cohort`). Then only the default is set and no steps are dropped.
    """
    if not study.steps_chosen:
        study.steps = study.default_steps()
    if study.coassemble and mode_final:
        dropped = [s for s in study.steps if s in COHORT_ONLY_STEPS]
        if dropped:
            # With one assembly and one bin set these steps do nothing useful at full cost.
            warning(
                "Co-assembling, so %s %s nothing to do (one assembly, one bin set); skipping. "
                "Set coassemble = false to run a cohort per sample."
                % (
                    ", ".join(dropped),
                    "has" if len(dropped) == 1 else "have",
                )
            )
            study.steps = tuple(s for s in study.steps if s not in COHORT_ONLY_STEPS)


def warn_if_cohort_steps_pinned(study: Study) -> None:
    """Warn when a cohort sheet's pinned `steps` omits the cohort-only stages.

    A sheet written for `metawrap2 run` usually omits `coverage` and `dereplicate`; used with
    `cohort`, it would produce per-sample bins and no catalogue.
    """
    if study.coassemble or not study.steps_chosen:
        return
    missing = [step for step in COHORT_ONLY_LOSSES if step not in study.steps]
    if not missing:
        return
    warning(
        "the sheet pins steps and omits %s, so this cohort will not produce %s. Remove `steps` "
        "from the sheet to use the cohort defaults (%s)."
        % (
            " and ".join(missing),
            " or ".join(COHORT_ONLY_LOSSES[step] for step in missing),
            " -> ".join(COHORT_STEPS),
        )
    )


def _as_option_list(value: object, module: str) -> List[str]:
    """Normalise a sheet's per-module options into an argv list."""
    if value is None:
        return []
    if isinstance(value, str):
        return shlex.split(value)
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    raise ValueError(
        "[modules.%s] options must be a string or a list, got %r" % (module, type(value).__name__)
    )


def load_sheet(path: str) -> Study:
    """Parse a sample sheet (TOML, or YAML when PyYAML is available)."""
    if path.lower().endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise ValueError(
                "%s looks like YAML but PyYAML is not installed. Either `pip install pyyaml` "
                "or write the sheet as TOML (see `metawrap2 run --example`)." % path
            ) from None
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    else:
        with open(path, "rb") as fh:
            data = _toml.load(fh)

    raw_samples = data.get("samples") or []
    if not raw_samples:
        raise ValueError("%s defines no samples. See `metawrap2 run --example`." % path)

    samples: List[Sample] = []
    seen = set()
    for i, entry in enumerate(raw_samples, start=1):
        if not isinstance(entry, dict):
            raise TypeError("sample %d in %s is not a table/mapping" % (i, path))
        r1 = entry.get("r1") or entry.get("reads_1")
        if not r1:
            raise ValueError("sample %d in %s has no 'r1'" % (i, path))
        r2 = entry.get("r2") or entry.get("reads_2")
        layout = entry.get("layout") or ("paired" if r2 else "single")
        if layout not in ("paired", "single", "interleaved"):
            raise ValueError(
                "sample %d in %s: layout must be paired/single/interleaved, "
                "got %r" % (i, path, layout)
            )
        if layout == "paired" and not r2:
            raise ValueError(
                "sample %d in %s: layout is 'paired' but no 'r2' was given" % (i, path)
            )
        name = entry.get("name") or _reads.sample_name(r1)
        if name in seen:
            raise ValueError("%s uses the sample name %r more than once" % (path, name))
        seen.add(name)
        samples.append(Sample(name=name, r1=r1, r2=r2, layout=layout))

    settings = data.get("settings") or {}
    chosen = bool(settings.get("steps"))
    coassemble = bool(settings.get("coassemble", True))
    coassemble_chosen = "coassemble" in settings
    steps = tuple(settings.get("steps") or (DEFAULT_STEPS if coassemble else COHORT_STEPS))
    unknown = [s for s in steps if s not in ALL_STEPS]
    if unknown:
        raise ValueError(
            "unknown step(s) %s in %s. Known steps: %s"
            % (", ".join(unknown), path, ", ".join(ALL_STEPS))
        )

    modules = data.get("modules") or {}
    module_options = {m: _as_option_list(v, m) for m, v in modules.items()}

    study = Study(
        samples=samples,
        coassemble=coassemble,
        coassemble_chosen=coassemble_chosen,
        steps=steps,
        module_options=module_options,
        steps_chosen=chosen,
    )
    resolve_steps(study, mode_final=False)
    return study


def study_from_directory(directory: str) -> Study:
    """Build a study by pairing up every read file in *directory*."""
    files = [os.path.join(directory, f) for f in sorted(os.listdir(directory))]
    pairs = _reads.collect_pairs(files)  # raises with an actionable message if none
    return Study(samples=[Sample(name=name, r1=r1, r2=r2) for name, r1, r2 in pairs])


#: The sheet `--example` prints. The binning tier comes from the binning module's default (the
#: sheet once listed the `fast` tier's binners as the default), and `coassemble` matches the
#: command that printed it.
EXAMPLE_SHEET = """\
# MetaWrap2 sample sheet. Run it with:  metawrap2 %(command)s samples.toml -o study/ -t 24

[settings]
# Co-assemble all samples together (true) or assemble each sample separately (false).
# `metawrap2 run` honours this line; `metawrap2 cohort` always assembles each sample separately.
coassemble = %(coassemble)s
# There is no `steps` line. The defaults are the genome-recovery path, and which steps that
# means depends on the mode: `metawrap2 cohort` adds the coverage and dereplicate stages that
# only exist when each sample is assembled on its own, and pinning a list here would switch them
# off without saying so. Use --steps/--skip to change the plan for one run. kraken2,
# classify_bins and blobology are opt-in either way, because of the databases they need.

# Extra options passed through to a module, exactly as you would type them. `binning` names the
# default tier instead of spelling out its binners, so it still means the default if the tiers
# change; `metawrap2 binning -h` lists them (--preset fast is cheaper, --preset gpu wants a
# GPU). Naming a binner (--metabat2 --semibin2 --concoct) replaces the tier rather than adding
# to it. The refinement thresholds are bin_refinement's own, written out because every MAG count
# this project reports uses them.
[modules]
binning = "--preset %(preset)s"
bin_refinement = "-c 70 -x 10"

[[samples]]
name = "ERR011347"
r1 = "RAW_READS/ERR011347_1.fastq.gz"
r2 = "RAW_READS/ERR011347_2.fastq.gz"

[[samples]]
name = "ERR011348"
r1 = "RAW_READS/ERR011348_1.fastq.gz"
r2 = "RAW_READS/ERR011348_2.fastq.gz"

# Single-end or interleaved samples set `layout` and omit r2:
# [[samples]]
# name = "SRR999"
# r1 = "RAW_READS/SRR999.fastq.gz"
# layout = "single"
"""


def example_sheet(cohort: bool = False) -> str:
    """Return the example sheet for `run` or `cohort`.

    The cohort sheet sets `coassemble = false`, so `metawrap2 cohort` on it prints no override
    warning.
    """
    return EXAMPLE_SHEET % {
        "command": "cohort" if cohort else "run",
        "coassemble": "false" if cohort else "true",
        "preset": _binning.DEFAULT_PRESET,
    }


# --- the layout the driver owns -----------------------------------------------------------


class Layout:
    """Paths inside the study directory.

    Every path the driver passes to or expects from a module is defined here.
    """

    def __init__(self, root: str):
        self.root = os.path.abspath(root)

    def _p(self, *parts: str) -> str:
        return os.path.join(self.root, *parts)

    @property
    def clean_reads(self) -> str:
        return self._p("CLEAN_READS")

    def read_qc(self, sample: str) -> str:
        return self._p("READ_QC", sample)

    def clean_pair(self, sample: str) -> Tuple[str, str]:
        return (
            os.path.join(self.clean_reads, "%s_1.fastq" % sample),
            os.path.join(self.clean_reads, "%s_2.fastq" % sample),
        )

    def clean_single(self, sample: str) -> str:
        return os.path.join(self.clean_reads, "%s.fastq" % sample)

    @property
    def assembly_dir(self) -> str:
        return self._p("ASSEMBLY")

    @property
    def assembly(self) -> str:
        return os.path.join(self.assembly_dir, "final_assembly.fasta")

    def sample_assembly_dir(self, sample: str) -> str:
        return self._p("ASSEMBLY", sample)

    def sample_assembly(self, sample: str) -> str:
        return os.path.join(self.sample_assembly_dir(sample), "final_assembly.fasta")

    @property
    def coverage_dir(self) -> str:
        return self._p("COVERAGE")

    def coverage_depth(self, sample: str) -> str:
        """Return the per-sample depth table `coverage` writes for *sample*."""
        return os.path.join(self.coverage_dir, _coverage.DEPTH_DIR, "%s.txt" % sample)

    def coverage_bams(self, sample: str) -> str:
        """Return the directory where `coverage --split-bams` writes one assembly's BAMs.

        Named for the assembly (one per sample in a cohort), it holds a BAM of every sample's reads
        against that assembly, for the binners that read alignments.
        """
        return os.path.join(self.coverage_dir, _coverage.BAMS_DIR, sample)

    @property
    def binning_dir(self) -> str:
        return self._p("INITIAL_BINNING")

    def binner_bins(self, binner: str) -> str:
        return os.path.join(self.binning_dir, "%s_bins" % binner)

    def sample_binning_dir(self, sample: str) -> str:
        return self._p("INITIAL_BINNING", sample)

    def sample_binner_bins(self, sample: str, binner: str) -> str:
        return os.path.join(self.sample_binning_dir(sample), "%s_bins" % binner)

    @property
    def refinement_dir(self) -> str:
        return self._p("BIN_REFINEMENT")

    def refined_bins(self, completeness: int, contamination: int) -> str:
        return os.path.join(
            self.refinement_dir, "metawrap_%d_%d_bins" % (completeness, contamination)
        )

    def sample_refinement_dir(self, sample: str) -> str:
        return self._p("BIN_REFINEMENT", sample)

    def sample_refined_bins(self, sample: str, completeness: int, contamination: int) -> str:
        return os.path.join(
            self.sample_refinement_dir(sample),
            "metawrap_%d_%d_bins" % (completeness, contamination),
        )

    @property
    def reassembly_dir(self) -> str:
        return self._p("BIN_REASSEMBLY")

    @property
    def reassembled_bins(self) -> str:
        return os.path.join(self.reassembly_dir, "reassembled_bins")

    def sample_reassembly_dir(self, sample: str) -> str:
        return self._p("BIN_REASSEMBLY", sample)

    def sample_reassembled_bins(self, sample: str) -> str:
        return os.path.join(self.sample_reassembly_dir(sample), "reassembled_bins")

    @property
    def catalogue_dir(self) -> str:
        return self._p("MAG_CATALOGUE")

    @property
    def catalogue_genomes(self) -> str:
        return os.path.join(self.catalogue_dir, _dereplicate.CATALOGUE_DIR)

    @property
    def kraken_dir(self) -> str:
        return self._p("KRAKEN")

    def sample_kraken_dir(self, sample: str) -> str:
        return self._p("KRAKEN", sample)

    def sample_blobology_dir(self, sample: str) -> str:
        return self._p("BLOBOLOGY", sample)

    @property
    def quant_dir(self) -> str:
        return self._p("QUANT_BINS")

    @property
    def classify_dir(self) -> str:
        return self._p("BIN_CLASSIFICATION")

    @property
    def annotate_dir(self) -> str:
        return self._p("FUNCT_ANNOT")

    @property
    def blobology_dir(self) -> str:
        return self._p("BLOBOLOGY")

    @property
    def state_dir(self) -> str:
        return self._p(_STATE_DIR)

    @property
    def log_dir(self) -> str:
        return self._p("LOGS")

    #: Top-level entries the driver writes itself, excluded from the end-of-run audit of
    #: undeclared outputs.
    DRIVER_OWNED = ("LOGS", "CLEAN_READS", "manifest.json", "run_config.json", "provenance.txt")


# --- steps -------------------------------------------------------------------------------


@dataclass
class Step:
    """One module invocation the driver will make."""

    name: str  # unique step id, e.g. "read_qc:ERR011347"
    module: str
    argv: List[str]
    parallel_group: Optional[str] = None  # steps sharing a group may run concurrently
    #: What this step reads and must produce, so --resume can check that outputs still exist and
    #: inputs have not changed (see metawrap2.manifest).
    inputs: List[str] = field(default_factory=list)
    outputs: List[str] = field(default_factory=list)

    def command(
        self,
        threads: int,
        config: Optional[str],
        global_flags: Sequence[str],
        memory_gb: Optional[int] = None,
    ) -> List[str]:
        argv = [sys.executable, "-m", "metawrap2.cli"] + list(global_flags) + [self.module]
        argv += self.argv + ["-t", str(threads)]
        # Modules that take -m get their share of the budget; N concurrent steps each given the
        # full -m could request N times the machine's RAM.
        if memory_gb is not None and self.module in MEMORY_MODULES and "-m" not in self.argv:
            argv += ["-m", str(memory_gb)]
        if config:
            argv += ["--config", config]
        return argv


def step_log_path(layout: Layout, step: Step) -> str:
    """Return the path of a step's live log.

    Used by the runner and by the opening and failure banners.
    """
    return os.path.join(layout.log_dir, "%s.log" % step.name.replace(":", "_"))


def _qc_output_reads(layout: Layout, sample: Sample) -> List[str]:
    """Return read_qc's cleaned reads for *sample*, before the driver links them."""
    out = layout.read_qc(sample.name)
    if sample.layout == "single":
        return [os.path.join(out, "final_pure_reads_1.fastq")]
    return [
        os.path.join(out, "final_pure_reads_1.fastq"),
        os.path.join(out, "final_pure_reads_2.fastq"),
    ]


def clean_reads_for(layout: Layout, study: Study, sample: Sample, ran_read_qc: bool) -> List[str]:
    """Return the reads later steps use for *sample*: QC'd if read_qc ran, else the input."""
    if not ran_read_qc:
        return list(sample.files)
    if sample.layout == "single":
        return [layout.clean_single(sample.name)]
    return list(layout.clean_pair(sample.name))


def plan_steps(study: Study, layout: Layout) -> List[Step]:
    """Turn a study into the ordered list of module invocations to make.

    It touches no files, so ``--dry-run`` shows the plan a real run would execute.
    """
    steps: List[Step] = []
    ran_read_qc = "read_qc" in study.steps
    all_reads: List[str] = []
    for sample in study.samples:
        all_reads += clean_reads_for(layout, study, sample, ran_read_qc)

    if ran_read_qc:
        for sample in study.samples:
            argv = ["-1", sample.r1]
            if sample.r2:
                argv += ["-2", sample.r2]
            argv += sample.layout_flags()
            argv += ["-o", layout.read_qc(sample.name)]
            steps.append(
                Step(
                    name="read_qc:%s" % sample.name,
                    module="read_qc",
                    argv=argv + study.extra("read_qc"),
                    parallel_group="read_qc",
                    inputs=sample.files,
                    outputs=_qc_output_reads(layout, sample),
                )
            )

    if "assembly" in study.steps:
        if study.coassemble:
            # Whether the co-assembly is paired depends on the samples' layout; single-end
            # samples produce only an ALL_READS_1.
            per_sample = [clean_reads_for(layout, study, s, ran_read_qc) for s in study.samples]
            paired_coassembly = all(len(r) > 1 for r in per_sample)
            argv = ["-1", os.path.join(layout.clean_reads, "ALL_READS_1.fastq")]
            if paired_coassembly:
                argv += ["-2", os.path.join(layout.clean_reads, "ALL_READS_2.fastq")]
            argv += ["-o", layout.assembly_dir]
            steps.append(
                Step(
                    name="assembly",
                    module="assembly",
                    argv=argv + study.extra("assembly"),
                    inputs=[os.path.join(layout.clean_reads, "ALL_READS_1.fastq")],
                    outputs=[layout.assembly],
                )
            )
        else:
            for sample in study.samples:
                reads = clean_reads_for(layout, study, sample, ran_read_qc)
                argv = ["-1", reads[0]]
                if len(reads) > 1:
                    argv += ["-2", reads[1]]
                argv += ["-o", layout.sample_assembly_dir(sample.name)]
                steps.append(
                    Step(
                        name="assembly:%s" % sample.name,
                        module="assembly",
                        argv=argv + study.extra("assembly"),
                        parallel_group="assembly",
                        inputs=reads,
                        outputs=[
                            os.path.join(
                                layout.sample_assembly_dir(sample.name), "final_assembly.fasta"
                            )
                        ],
                    )
                )

    binning_opts = study.extra("binning") or list(DEFAULT_BINNERS)
    comp, cont = _quality_thresholds(study)
    if study.coassemble:
        steps += _coassembly_steps(study, layout, all_reads, binning_opts, comp, cont)
    else:
        steps += _cohort_steps(study, layout, ran_read_qc, all_reads, binning_opts, comp, cont)
    return steps


def _final_bins(study: Study, layout: Layout, binning_opts: Sequence[str]) -> str:
    """Return the bin folder later steps read in a co-assembly study: the last one produced."""
    comp, cont = _quality_thresholds(study)
    if "reassemble_bins" in study.steps:
        return layout.reassembled_bins
    if "bin_refinement" in study.steps:
        return layout.refined_bins(comp, cont)
    return layout.binner_bins(selected_binners(binning_opts)[0])


def _sample_final_bins(
    study: Study, layout: Layout, sample: str, binning_opts: Sequence[str]
) -> str:
    """Return the bin folder later steps read for one sample of a cohort."""
    comp, cont = _quality_thresholds(study)
    if "reassemble_bins" in study.steps:
        return layout.sample_reassembled_bins(sample)
    if "bin_refinement" in study.steps:
        return layout.sample_refined_bins(sample, comp, cont)
    return layout.sample_binner_bins(sample, selected_binners(binning_opts)[0])


def _coassembly_steps(
    study: Study,
    layout: Layout,
    all_reads: List[str],
    binning_opts: List[str],
    comp: int,
    cont: int,
) -> List[Step]:
    """Plan the steps after assembly when all samples are co-assembled.

    One assembly and one bin set, so one of each downstream step.
    """
    steps: List[Step] = []
    if "kraken2" in study.steps:
        steps.append(
            Step(
                name="kraken2",
                module="kraken2",
                argv=["-o", layout.kraken_dir]
                + study.extra("kraken2")
                + all_reads
                + [layout.assembly],
                inputs=all_reads + [layout.assembly],
                outputs=[os.path.join(layout.kraken_dir, "kronagram.html")],
            )
        )

    if "binning" in study.steps:
        steps.append(
            Step(
                name="binning",
                module="binning",
                argv=["-a", layout.assembly, "-o", layout.binning_dir] + binning_opts + all_reads,
                inputs=[layout.assembly] + all_reads,
                outputs=[layout.binner_bins(b) for b in selected_binners(binning_opts)],
            )
        )

    if "bin_refinement" in study.steps:
        # Pass only the bin sets binning produces, as repeated --bins (documented for generated
        # arguments). A fixed -A/-B/-C failed when other binners were selected, and the letters
        # stop at F.
        selected = selected_binners(binning_opts)
        argv = ["-o", layout.refinement_dir]
        for binner in selected:
            argv += ["--bins", layout.binner_bins(binner)]
        steps.append(
            Step(
                name="bin_refinement",
                module="bin_refinement",
                argv=argv + study.extra("bin_refinement"),
                inputs=[layout.binner_bins(b) for b in selected],
                outputs=[layout.refined_bins(comp, cont)],
            )
        )

    # Without refinement, use the first binner that was run.
    best_bins = (
        layout.refined_bins(comp, cont)
        if "bin_refinement" in study.steps
        else layout.binner_bins(selected_binners(binning_opts)[0])
    )

    if "reassemble_bins" in study.steps:
        reads = clean_reads_for(layout, study, study.samples[0], "read_qc" in study.steps)
        argv = [
            "-b",
            best_bins,
            "-o",
            layout.reassembly_dir,
            "-1",
            os.path.join(layout.clean_reads, "ALL_READS_1.fastq"),
        ]
        if len(reads) > 1:
            argv += ["-2", os.path.join(layout.clean_reads, "ALL_READS_2.fastq")]
        steps.append(
            Step(
                name="reassemble_bins",
                module="reassemble_bins",
                argv=argv + study.extra("reassemble_bins"),
                inputs=[best_bins, os.path.join(layout.clean_reads, "ALL_READS_1.fastq")],
                outputs=[layout.reassembled_bins],
            )
        )
    final_bins = _final_bins(study, layout, binning_opts)

    if "quant_bins" in study.steps:
        steps.append(
            Step(
                name="quant_bins",
                module="quant_bins",
                argv=["-b", final_bins, "-a", layout.assembly, "-o", layout.quant_dir]
                + study.extra("quant_bins")
                + all_reads,
                inputs=[final_bins, layout.assembly] + all_reads,
                outputs=[os.path.join(layout.quant_dir, "bin_abundance_table.tsv")],
            )
        )

    steps += _bin_annotation_steps(study, layout, final_bins)

    if "blobology" in study.steps:
        steps.append(
            Step(
                name="blobology",
                module="blobology",
                argv=["-a", layout.assembly, "-o", layout.blobology_dir, "--bins", final_bins]
                + study.extra("blobology")
                + all_reads,
                inputs=[layout.assembly, final_bins] + all_reads,
                outputs=[os.path.join(layout.blobology_dir, "blobplot_figures")],
            )
        )
    return steps


def _bin_annotation_steps(study: Study, layout: Layout, bins: str) -> List[Step]:
    """Plan classify_bins and annotate_bins, which take only a folder of bins.

    Shared by both modes, which differ only in the folder: the catalogue or the single bin set.
    """
    steps: List[Step] = []
    if "classify_bins" in study.steps:
        steps.append(
            Step(
                name="classify_bins",
                module="classify_bins",
                argv=["-b", bins, "-o", layout.classify_dir] + study.extra("classify_bins"),
                inputs=[bins],
                outputs=[os.path.join(layout.classify_dir, "bin_taxonomy.tsv")],
            )
        )
    if "annotate_bins" in study.steps:
        steps.append(
            Step(
                name="annotate_bins",
                module="annotate_bins",
                argv=["-b", bins, "-o", layout.annotate_dir] + study.extra("annotate_bins"),
                inputs=[bins],
                outputs=[os.path.join(layout.annotate_dir, "bin_funct_annotations")],
            )
        )
    return steps


def _cohort_steps(
    study: Study,
    layout: Layout,
    ran_read_qc: bool,
    all_reads: List[str],
    binning_opts: List[str],
    comp: int,
    cont: int,
) -> List[Step]:
    """Plan the steps after assembly when each sample is assembled separately.

    One coverage pass across all assemblies, per-sample binning on that depth table, per-sample
    refinement, then one dereplication into a catalogue used by quantification and annotation.

    Each per-sample stage is a parallel group, so samples run through a stage together and the
    driver divides the thread budget between them.
    """
    steps: List[Step] = []
    samples = [s.name for s in study.samples]
    reads_of = {s.name: clean_reads_for(layout, study, s, ran_read_qc) for s in study.samples}

    if "kraken2" in study.steps:
        for sample in samples:
            steps.append(
                Step(
                    name="kraken2:%s" % sample,
                    module="kraken2",
                    argv=["-o", layout.sample_kraken_dir(sample)]
                    + study.extra("kraken2")
                    + reads_of[sample]
                    + [layout.sample_assembly(sample)],
                    parallel_group="kraken2",
                    inputs=reads_of[sample] + [layout.sample_assembly(sample)],
                    outputs=[os.path.join(layout.sample_kraken_dir(sample), "kronagram.html")],
                )
            )

    # One coverage pass for the cohort: every sample's reads against the concatenation of all
    # assemblies. That is N alignments instead of N per assembly.
    cohort_coverage = "coverage" in study.steps
    # metaBAT2 and MaxBin2 use a depth table; SemiBin2, VAMB, COMEBin and CONCOCT read alignments
    # and refuse `--depth`. If any of those is selected (the default includes SemiBin2), coverage
    # must also write per-assembly BAMs.
    needs_bams = [b for b in selected_binners(binning_opts) if b in _binning.BAM_BINNERS]
    if cohort_coverage:
        argv: List[str] = []
        for sample in samples:
            argv += ["-a", "%s=%s" % (sample, layout.sample_assembly(sample))]
        argv += ["-o", layout.coverage_dir]
        outputs = [layout.coverage_depth(s) for s in samples]
        if needs_bams and "binning" in study.steps:
            argv += ["--split-bams"]
            comm(
                "%s %s the alignments, not a coverage table, so the coverage stage also splits the "
                "catalogue BAMs into one set per assembly"
                % (
                    ", ".join("--" + b for b in needs_bams),
                    "read" if len(needs_bams) > 1 else "reads",
                )
            )
        steps.append(
            Step(
                name="coverage",
                module="coverage",
                argv=argv + study.extra("coverage") + all_reads,
                inputs=[layout.sample_assembly(s) for s in samples] + all_reads,
                outputs=outputs,
            )
        )

    if "binning" in study.steps:
        for sample in samples:
            argv = ["-a", layout.sample_assembly(sample), "-o", layout.sample_binning_dir(sample)]
            if cohort_coverage:
                # This sample's contigs with every sample's depths; no alignment in this step.
                argv += ["--depth", layout.coverage_depth(sample)]
                inputs = [layout.sample_assembly(sample), layout.coverage_depth(sample)]
                if needs_bams:
                    # BAMs for the binners that need alignments; the table still serves metaBAT2
                    # and MaxBin2.
                    argv += ["--bams", layout.coverage_bams(sample)]
            else:
                # Without the coverage stage, binning uses single-sample coverage.
                argv += reads_of[sample]
                inputs = [layout.sample_assembly(sample)] + reads_of[sample]
            steps.append(
                Step(
                    name="binning:%s" % sample,
                    module="binning",
                    argv=argv + binning_opts,
                    parallel_group="binning",
                    inputs=inputs,
                    outputs=[
                        layout.sample_binner_bins(sample, b) for b in selected_binners(binning_opts)
                    ],
                )
            )

    if "bin_refinement" in study.steps:
        selected = selected_binners(binning_opts)
        for sample in samples:
            argv = ["-o", layout.sample_refinement_dir(sample)]
            for binner in selected:
                argv += ["--bins", layout.sample_binner_bins(sample, binner)]
            steps.append(
                Step(
                    name="bin_refinement:%s" % sample,
                    module="bin_refinement",
                    argv=argv + study.extra("bin_refinement"),
                    parallel_group="bin_refinement",
                    inputs=[layout.sample_binner_bins(sample, b) for b in selected],
                    outputs=[layout.sample_refined_bins(sample, comp, cont)],
                )
            )

    if "reassemble_bins" in study.steps:
        for sample in samples:
            bins = (
                layout.sample_refined_bins(sample, comp, cont)
                if "bin_refinement" in study.steps
                else layout.sample_binner_bins(sample, selected_binners(binning_opts)[0])
            )
            reads = reads_of[sample]
            argv = ["-b", bins, "-o", layout.sample_reassembly_dir(sample), "-1", reads[0]]
            if len(reads) > 1:
                argv += ["-2", reads[1]]
            steps.append(
                Step(
                    name="reassemble_bins:%s" % sample,
                    module="reassemble_bins",
                    argv=argv + study.extra("reassemble_bins"),
                    parallel_group="reassemble_bins",
                    inputs=[bins] + reads,
                    outputs=[layout.sample_reassembled_bins(sample)],
                )
            )

    per_sample_bins = [_sample_final_bins(study, layout, s, binning_opts) for s in samples]
    catalogue = layout.catalogue_genomes
    if "dereplicate" in study.steps:
        argv = []
        for bins in per_sample_bins:
            argv += ["-b", bins]
        argv += ["-o", layout.catalogue_dir]
        steps.append(
            Step(
                name="dereplicate",
                module="dereplicate",
                argv=argv + study.extra("dereplicate"),
                inputs=per_sample_bins,
                outputs=[catalogue],
            )
        )
    # Without dereplication there is no single cohort bin folder, so the first sample's is used.
    final_bins = catalogue if "dereplicate" in study.steps else per_sample_bins[0]

    if "quant_bins" in study.steps:
        # Against the catalogue, with no -a: near-identical MAGs from different samples would
        # split reads, and there is no single assembly.
        steps.append(
            Step(
                name="quant_bins",
                module="quant_bins",
                argv=["-b", final_bins, "-o", layout.quant_dir]
                + study.extra("quant_bins")
                + all_reads,
                inputs=[final_bins] + all_reads,
                outputs=[os.path.join(layout.quant_dir, "bin_abundance_table.tsv")],
            )
        )

    steps += _bin_annotation_steps(study, layout, final_bins)

    if "blobology" in study.steps:
        for sample, bins in zip(samples, per_sample_bins):
            steps.append(
                Step(
                    name="blobology:%s" % sample,
                    module="blobology",
                    argv=[
                        "-a",
                        layout.sample_assembly(sample),
                        "-o",
                        layout.sample_blobology_dir(sample),
                        "--bins",
                        bins,
                    ]
                    + study.extra("blobology")
                    + reads_of[sample],
                    parallel_group="blobology",
                    inputs=[layout.sample_assembly(sample), bins] + reads_of[sample],
                    outputs=[os.path.join(layout.sample_blobology_dir(sample), "blobplot_figures")],
                )
            )
    return steps


def selected_binners(binning_options: Sequence[str]) -> List[str]:
    """Return the binners a set of `binning` options selects, in the order refinement gets them.

    The order comes from the binning module's ``BINNERS``, so new binners need no edit here.
    ``--preset`` is resolved from the binning module's table with the same precedence: named
    binners win, and the preset is used only when none is named.
    """
    selected = [b for b in _binning.BINNERS if "--" + b in binning_options]
    if selected:
        return selected
    if "--preset" in binning_options:
        index = binning_options.index("--preset") + 1
        if index < len(binning_options):
            preset = _binning.PRESETS.get(_binning.resolve_preset(binning_options[index]) or "")
            if preset:
                return [b for b in _binning.BINNERS if b in preset]
    return [b.lstrip("-") for b in DEFAULT_BINNERS]


def _quality_thresholds(study: Study) -> Tuple[int, int]:
    """Return the -c/-x the refinement step will use, which name its output directory.

    bin_refinement writes ``metawrap_<c>_<x>_bins``, and the driver passes that path on.
    """
    comp, cont = 70, 10  # bin_refinement's own defaults
    opts = study.extra("bin_refinement")
    for flag, index in (("-c", 0), ("-x", 1)):
        if flag in opts:
            try:
                value = int(float(opts[opts.index(flag) + 1]))
            except (IndexError, ValueError):
                continue
            if index == 0:
                comp = value
            else:
                cont = value
    return comp, cont


# --- step state ---------------------------------------------------------------------------


class StepState:
    """Decide what a resume may skip, using the study's manifest.

    The manifest records each step's command, timing, declared inputs and outputs with content
    fingerprints, and its outcome. A step re-runs when its outputs are missing or its inputs
    have changed.
    """

    def __init__(self, layout: Layout, resume: bool, strict_fingerprints: bool = False):
        self.manifest = Manifest(os.path.join(layout.state_dir, MANIFEST_NAME))
        self.resume = resume
        self.strict = strict_fingerprints

    def skip_reason(self, step: Step) -> Optional[str]:
        """Return why *step* may be skipped, or None if it must run."""
        if not self.resume:
            return None
        ok, reason = self.manifest.can_skip(step.name, whole_file=self.strict)
        return reason if ok else None

    def rerun_reason(self, step: Step) -> str:
        """Return why a step is re-run despite --resume, for a message."""
        _ok, reason = self.manifest.can_skip(step.name, whole_file=self.strict)
        return reason

    def start(
        self, step: Step, command: Sequence[str], grant: Optional[Dict[str, Any]] = None
    ) -> float:
        return self.manifest.start(step.name, command, step.inputs, grant=grant)

    def finish(
        self,
        step: Step,
        started: float,
        status: str,
        log: Optional[str] = None,
        usage: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.manifest.finish(step.name, started, status, outputs=step.outputs, log=log, usage=usage)


def strategy_record(study: Study) -> Dict[str, Any]:
    """Return how this study's catalogue was built: strategy, sample count, and dereplication.

    A per-sample study without `dereplicate` keeps redundant per-sample bin sets, so its genome
    count is not comparable with a co-assembly's.
    """
    return {
        "strategy": CO_ASSEMBLY if study.coassemble else PER_SAMPLE,
        "samples": len(study.samples),
        "dereplicated": "dereplicate" in study.steps,
    }


def record_strategy(state: StepState, study: Study) -> Dict[str, Any]:
    """Record the strategy in the study manifest under a top-level `study` key.

    `metawrap2 status <study> --json` prints it. No schema change is needed: older manifests have
    no `study` block, and readers use `data["steps"]`.
    """
    record = strategy_record(study)
    state.manifest.data[STRATEGY_KEY] = record
    state.manifest.save()
    return record


# --- the driver ---------------------------------------------------------------------------


def _link_or_copy(src: str, dest: str) -> None:
    """Symlink *dest* to *src*, or copy if symlinks are not supported."""
    if os.path.lexists(dest):
        os.remove(dest)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        os.symlink(os.path.abspath(src), dest)
    except OSError:  # e.g. a filesystem with no symlink support
        import shutil

        shutil.copyfile(src, dest)


def collect_clean_reads(layout: Layout, study: Study, dry_run: bool) -> None:
    """Expose read_qc's output as CLEAN_READS/<sample>_1.fastq etc. for the later steps.

    Uses symlinks (the tutorial uses ``mv``), so the QC output directory stays intact for its
    provenance record and for --resume.
    """
    for sample in study.samples:
        produced = _qc_output_reads(layout, sample)
        if sample.layout == "single":
            targets = [layout.clean_single(sample.name)]
        else:
            targets = list(layout.clean_pair(sample.name))
        for src, dest in zip(produced, targets):
            if dry_run:
                comm("(dry run) would link %s -> %s" % (dest, src))
                continue
            if not os.path.isfile(src):
                error(
                    "read_qc did not produce %s.\nThe %s step must have failed; look at "
                    "%s."
                    % (
                        src,
                        sample.name,
                        os.path.join(layout.log_dir, "read_qc_%s.log" % sample.name),
                    )
                )
            _link_or_copy(src, dest)


def concatenate_for_coassembly(layout: Layout, study: Study, dry_run: bool) -> None:
    """Build the ALL_READS files a co-assembly needs, skipping the work if they are current."""
    groups: List[Tuple[str, List[str]]] = []
    first = clean_reads_for(layout, study, study.samples[0], "read_qc" in study.steps)
    for index in range(len(first)):
        sources = []
        for sample in study.samples:
            reads = clean_reads_for(layout, study, sample, "read_qc" in study.steps)
            if index < len(reads):
                sources.append(reads[index])
        groups.append(
            (os.path.join(layout.clean_reads, "ALL_READS_%d.fastq" % (index + 1)), sources)
        )

    for dest, sources in groups:
        if dry_run:
            comm("(dry run) would concatenate %d file(s) into %s" % (len(sources), dest))
            continue
        if _is_current(dest, sources):
            comm("%s is already up to date" % dest)
            continue
        comm("concatenating %d file(s) into %s" % (len(sources), dest))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        with open(dest, "wb") as out:
            for src in sources:
                with open(src, "rb") as fh:
                    # copyfileobj streams, so the reads are not held in memory
                    import shutil

                    shutil.copyfileobj(fh, out, length=1024 * 1024)


def _is_current(dest: str, sources: List[str]) -> bool:
    """Return True if *dest* exists, is non-empty, and is newer than every source."""
    if not os.path.isfile(dest) or not os.path.getsize(dest):
        return False
    dest_mtime = os.path.getmtime(dest)
    return all(os.path.isfile(s) and os.path.getmtime(s) <= dest_mtime for s in sources)


def run_step(
    step: Step,
    layout: Layout,
    threads: int,
    config: Optional[str],
    global_flags: Sequence[str],
    dry_run: bool,
    state: Optional[StepState] = None,
    memory_gb: Optional[int] = None,
) -> Tuple[Step, str, str]:
    """Run one step as its own ``metawrap2 <module>`` process.

    Returns ``(step, status, log path)`` where status is one of the manifest's
    ``completed`` / ``failed`` / ``interrupted``. An interrupted step's outputs are in an unknown
    state, so ``--resume`` always re-runs it; a failed step's log may show the reason.
    """
    os.makedirs(layout.log_dir, exist_ok=True)
    log_path = step_log_path(layout, step)
    argv = step.command(threads, config, global_flags, memory_gb=memory_gb)
    if dry_run:
        print("[plan] %s" % " ".join(shlex.quote(a) for a in argv))
        return step, COMPLETED, log_path

    grant: Dict[str, Any] = {"threads": threads}
    if memory_gb is not None:
        grant["memory_gb"] = memory_gb
    started = state.start(step, argv, grant=grant) if state is not None else time.time()
    status = FAILED
    usage: Dict[str, Any] = {}
    try:
        with open(log_path, "w") as log:
            log.write("$ %s\n\n" % " ".join(shlex.quote(a) for a in argv))
            log.flush()
            # Measure peak memory and CPU time (see metawrap2.usage), which `calibrate` uses.
            # The driver writes the cost line, so the child is told it runs under `run` and does
            # not print its own. The variable is set only in the child's environment.
            step_env = dict(os.environ)
            step_env[_common.DRIVEN_BY_RUN_ENV] = "1"
            # Pass threads and step so the utilisation warning names the step, which matters when
            # several samples run at once.
            rc, usage = _usage.run_measured(
                argv,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=step_env,
                threads=threads,
                step=step.name,
            )
            log.write("\n[metawrap2] %s\n" % _usage.summarise(usage))
        status = COMPLETED if rc == 0 else FAILED
    except KeyboardInterrupt:
        # Ctrl-C leaves outputs half-written; record it so --resume redoes the step.
        status = INTERRUPTED
        raise
    finally:
        if state is not None:
            state.finish(step, started, status, log=log_path, usage=usage)
    return step, status, log_path


def cores_busy(usage: Dict[str, Any]) -> Optional[float]:
    """Return the mean busy cores for a step: CPU seconds over wall seconds.

    The same figure ``metawrap2.usage.summarise`` prints. None when CPU time was not measured.
    """
    if "user_seconds" not in usage:
        return None
    wall = float(usage.get("wall_seconds") or 0)
    if wall <= 0:
        return None
    return (float(usage.get("user_seconds", 0)) + float(usage.get("system_seconds", 0))) / wall


def parallelism_advice(
    module: str, busy: float, threads: int, samples: int, workers: int
) -> Optional[str]:
    """Return a ``-j`` suggestion for a per-sample step that left most of its threads idle, or None.

    For example, read_qc using 1.7 cores of a large grant at ``-j 1`` over three samples takes
    three times the wall time it needs. This is only advice: a larger -j also multiplies peak
    memory, which has not been measured for concurrent runs.
    """
    if busy > threads * IDLE_GRANT_FRACTION:
        return None
    # Do not suggest more concurrent steps than the grant can supply at the measured usage.
    fits = max(1, int(threads // max(1.0, busy)))
    suggest = min(samples, fits)
    if suggest <= workers:
        return None
    return (
        "%s used %.1f of its %d granted cores; with %d sample(s) to run, -j %d would spend the "
        "same CPU in less wall time." % (module, busy, threads, samples, suggest)
    )


def execute(
    steps: List[Step],
    layout: Layout,
    study: Study,
    threads: int,
    jobs: int,
    config: Optional[str],
    global_flags: Sequence[str],
    state: StepState,
    dry_run: bool,
    memory_gb: Optional[int] = None,
) -> int:
    """Run the planned steps in order, grouping the parallelisable ones. Returns an exit code."""
    index = 0
    completed = 0
    advised = False  # whether the -j suggestion has already been made for this run
    while index < len(steps):
        step = steps[index]
        # Gather a run of consecutive steps that share a parallel group.
        group = [step]
        if step.parallel_group:
            while (
                index + len(group) < len(steps)
                and steps[index + len(group)].parallel_group == step.parallel_group
            ):
                group.append(steps[index + len(group)])
        index += len(group)

        pending = []
        for s in group:
            reason = state.skip_reason(s)
            if reason:
                comm("skipping %s - %s" % (s.name, reason))
                completed += 1
            else:
                if state.resume and state.manifest.step(s.name) is not None:
                    # Say why --resume is re-running this step.
                    comm("re-running %s - %s" % (s.name, state.rerun_reason(s)))
                pending.append(s)

        if not pending:
            _after_group(group[0], layout, study, dry_run)
            continue

        _before_group(group[0], layout, study, dry_run)

        workers = max(1, min(jobs, len(pending)))
        per_step_threads = max(1, threads // workers)
        # Divide the memory budget like the threads, so N concurrent steps (e.g. SPAdes) do not
        # each request the whole budget.
        per_step_memory = max(1, memory_gb // workers) if memory_gb else None
        if len(pending) > 1:
            announcement(
                "RUNNING %d x %s (%d at a time, %d threads%s each)"
                % (
                    len(pending),
                    group[0].module,
                    workers,
                    per_step_threads,
                    " and %dGB" % per_step_memory if per_step_memory else "",
                )
            )
        else:
            announcement("RUNNING %s" % pending[0].name.upper())
        # Print the live log path before the step starts, since the console is quiet during long
        # steps.
        if len(pending) == 1:
            comm("live log: %s" % step_log_path(layout, pending[0]))
        else:
            comm("live logs: %s" % os.path.join(layout.log_dir, "<step>.log"))

        failures = []
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [
                    pool.submit(
                        run_step,
                        s,
                        layout,
                        per_step_threads,
                        config,
                        global_flags,
                        dry_run,
                        state,
                        per_step_memory,
                    )
                    for s in pending
                ]
                for future in concurrent.futures.as_completed(futures):
                    s, status, log_path = future.result()
                    if status == COMPLETED:
                        completed += 1
                        record = state.manifest.step(s.name)
                        cost = (record.data.get("usage") or {}) if record else {}
                        comm(
                            "finished %s (%s)" % (s.name, _usage.summarise(cost))
                            if cost
                            else "finished %s" % s.name
                        )
                        # Advise once per run; later steps of a stage would repeat it.
                        busy = cores_busy(cost) if not advised and len(pending) > 1 else None
                        if busy is not None:
                            advice = parallelism_advice(
                                s.module, busy, per_step_threads, len(pending), workers
                            )
                            if advice:
                                comm(advice)
                                advised = True
                    else:
                        failures.append((s, log_path))
        except KeyboardInterrupt:
            # run_step's finally has recorded the running steps as interrupted, so --resume
            # redoes only those.
            warning(
                "Interrupted. The running steps are recorded as interrupted. Re-run the same "
                "command with --resume to continue: those steps are redone and finished steps "
                "are skipped."
            )
            return 130

        if failures:
            for s, log_path in failures:
                warning("%s failed. Its log is %s" % (s.name, log_path))
            error(
                "%d of %d step(s) failed; stopping.\n"
                "Fix the cause and re-run the same command with --resume to continue from "
                "this point." % (len(failures), len(pending))
            )

        _after_group(group[0], layout, study, dry_run)

    announcement("STUDY COMPLETE - %d step(s) finished" % completed)
    if not dry_run:
        audit_outputs(state, layout)
    comm("results are under %s" % layout.root)
    return 0


def audit_outputs(state: StepState, layout: Layout) -> None:
    """Check at the end of a run that the manifest matches what is on disk.

    ``--resume`` relies on the manifest, so mismatches are reported now. They are warnings, not
    errors, since every step reported success.
    """
    problems = state.manifest.audit(layout.root, ignore=Layout.DRIVER_OWNED)
    for path in problems["missing"]:
        warning("a step reported success but its output is missing: %s" % path)
    for path in problems["empty"]:
        warning("a step reported success but its output is empty: %s" % path)
    for name in problems["no_outputs"]:
        warning(
            "%s declared no outputs, so --resume cannot verify it and will always re-run it" % name
        )
    for path in problems["undeclared"]:
        comm("note: %s was written but no step declared it as an output" % path)
    if not any(problems.values()):
        comm("output audit: every recorded output is present and non-empty")


def _before_group(step: Step, layout: Layout, study: Study, dry_run: bool) -> None:
    """Prepare inputs a step needs before it runs.

    A co-assembly reads ``ALL_READS_1.fastq``. It is built here so that a study without read_qc
    also gets it (from the raw reads); after read_qc the files are already current and this is a
    no-op.
    """
    if step.module == "assembly" and study.coassemble:
        announcement("CONCATENATING READS FOR THE CO-ASSEMBLY")
        concatenate_for_coassembly(layout, study, dry_run)


def _after_group(step: Step, layout: Layout, study: Study, dry_run: bool) -> None:
    """Do the file handling between steps that the tutorial does with ``mv`` and ``cat``."""
    if step.module == "read_qc":
        announcement("COLLECTING THE QC'ED READS")
        collect_clean_reads(layout, study, dry_run)
    if step.module == "read_qc" and study.coassemble:
        concatenate_for_coassembly(layout, study, dry_run)


# --- entry point --------------------------------------------------------------------------


def main(argv: List[str], prog: str = "metawrap2 run", cohort: bool = False) -> int:
    """Plan and execute a study.

    `metawrap2 cohort` passes *cohort* to force per-sample assembly with the same planner.
    """
    ap = make_parser(
        prog=prog,
        description=(
            "Assemble and bin each sample of a cohort separately, then collapse the results "
            "into one non-redundant genome catalogue."
            if cohort
            else "Run the whole MetaWrap2 pipeline for a study described by a sample sheet."
        ),
    )
    ap.add_argument(
        "sheet", nargs="?", help="sample sheet (.toml/.yaml), or a directory of read files"
    )
    ap.add_argument("-o", "--output", help="study output directory")
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=None,
        help="total threads to use, divided between concurrent steps "
        "(default: [settings] threads in metawrap2.toml, else the [limits] budget)",
    )
    ap.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="how many independent samples to process at once; the thread budget "
        "is divided between them (default 1)",
    )
    ap.add_argument(
        "--steps", nargs="+", choices=ALL_STEPS, help="override the steps from the sheet"
    )
    ap.add_argument(
        "--skip", nargs="+", choices=ALL_STEPS, default=(), help="drop these steps from the plan"
    )
    ap.add_argument(
        "--no-coassembly",
        action="store_true",
        help=(
            argparse.SUPPRESS
            if cohort
            else "assemble each sample separately instead of co-assembling"
        ),
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print the plan (every module command) and run nothing",
    )
    ap.add_argument(
        "--resume",
        action="store_true",
        help="skip steps that already finished in this output directory",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="pass --force to each module, overwriting existing output",
    )
    ap.add_argument(
        "-m",
        "--memory",
        type=int,
        default=None,
        help="total memory budget in GB, divided between concurrent steps "
        "(default: the [limits] budget)",
    )
    ap.add_argument(
        "--strict-fingerprints",
        action="store_true",
        help="hash whole files when deciding whether a --resume may skip a step, "
        "instead of sampling their ends (slower, exact)",
    )
    ap.add_argument("--config", help="path to metawrap2.toml")
    ap.add_argument("--example", action="store_true", help="print an example sample sheet and exit")
    args = ap.parse_args(argv)

    if args.example:
        print(example_sheet(cohort=cohort))
        return 0
    if not args.sheet or not args.output:
        ap.error(
            "a sample sheet (or directory) and -o/--output are required; "
            "see `metawrap2 run --example`"
        )

    try:
        study = (
            study_from_directory(args.sheet)
            if os.path.isdir(args.sheet)
            else load_sheet(args.sheet)
        )
    except (OSError, ValueError, TypeError) as exc:
        # Report a malformed sheet as a user error, not a traceback (a non-mapping sample entry
        # raises TypeError).
        error(str(exc))

    if args.steps:
        study.steps = tuple(args.steps)
        study.steps_chosen = True
    if args.no_coassembly or cohort:
        if cohort and study.coassemble_chosen and study.coassemble:
            # The sheet asks for a co-assembly but the command is `cohort`; the command wins,
            # with a warning.
            warning(
                "the sheet sets coassemble = true, but `metawrap2 cohort` assembles each sample "
                "separately; running per-sample. Use `metawrap2 run` to honour the sheet."
            )
        study.coassemble = False
    # Resolve after the mode is final, since --no-coassembly changes the default steps and --skip
    # applies to the result.
    resolve_steps(study)
    if not args.steps:
        # Only for steps pinned in the sheet, not for `--steps` on the command line.
        warn_if_cohort_steps_pinned(study)
    if args.skip:
        study.steps = tuple(s for s in study.steps if s not in args.skip)
    if not study.steps:
        error("no steps left to run.")

    # metawrap2.cli strips the global --dry-run/--force/--resume and applies them to the shared
    # runner, so accept them from either place: `metawrap2 --dry-run run ...` and
    # `metawrap2 run ... --dry-run` behave the same.
    from .. import command as _command

    dry_run = bool(args.dry_run or _command.runner.dry_run)
    resume = bool(args.resume or _command.runner.resume)
    force = bool(args.force or _command.runner.force)

    settings = load_settings(args.config)
    # The driver knows how many steps run at once, so it resolves the budget and divides it;
    # otherwise each concurrent module would claim the whole budget. See metawrap2.limits.
    budget = _limits.resolve(settings.limits)
    threads = args.threads or settings.threads or budget.threads
    if threads > budget.threads:
        if args.threads is not None:
            warning(
                "Asked for %d threads but the budget allows %d of this machine's %d cores; "
                "using %d. Raise [limits] max_thread_fraction or max_threads to use more."
                % (threads, budget.threads, budget.machine_cores, budget.threads)
            )
        threads = budget.threads
    memory_gb = args.memory or max(1, int(budget.memory_gb))
    if memory_gb > budget.memory_gb:
        if args.memory is not None:
            warning(
                "Asked for %d GB but the budget allows %d GB of this machine's %s; using %d. "
                "Raise [limits] max_memory_fraction or max_memory_gb to use more."
                % (
                    memory_gb,
                    int(budget.memory_gb),
                    _limits.human(budget.machine_memory_bytes),
                    int(budget.memory_gb),
                )
            )
        memory_gb = max(1, int(budget.memory_gb))
    layout = Layout(args.output)
    os.makedirs(layout.root, exist_ok=True)

    global_flags: List[str] = []
    if force:
        global_flags.append("--force")
    if resume:
        global_flags.append("--resume")
    if dry_run:
        global_flags.append("--dry-run")

    steps = plan_steps(study, layout)

    announcement("METAWRAP2 STUDY: %d SAMPLE(S), %d STEP(S)" % (len(study.samples), len(steps)))
    comm("samples: %s" % ", ".join(s.name for s in study.samples))
    comm("steps: %s" % " -> ".join(study.steps))
    comm(
        "%s assembly; %d threads total, %d job(s) at a time"
        % ("co-" if study.coassemble else "per-sample", threads, args.jobs)
    )
    comm("output: %s" % layout.root)
    strategy = strategy_record(study)
    comm(
        "catalogue strategy: %s, %d sample(s), dereplication %s"
        % (
            strategy["strategy"],
            strategy["samples"],
            "on" if strategy["dereplicated"] else "off",
        )
    )
    comm("Resource budget: %s" % budget.describe())
    overrides = _limits.overrides()
    if overrides:
        comm(
            "Budget overridden on the command line: %s"
            % ", ".join("--%s %s" % (k.replace("_", "-"), v) for k, v in sorted(overrides.items()))
        )

    state = StepState(layout, resume=resume, strict_fingerprints=args.strict_fingerprints)
    # Written before the first step, so a failed or interrupted study still records its strategy.
    # Not written under --dry-run.
    if not dry_run:
        record_strategy(state, study)
    return execute(
        steps,
        layout,
        study,
        threads,
        args.jobs,
        args.config,
        global_flags,
        state,
        dry_run,
        memory_gb=memory_gb,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
