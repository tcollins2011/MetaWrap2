"""`metawrap2 run <samples.yaml|toml> -o study/` - drive the whole pipeline from one file.

Until now a study meant hand-chaining ten modules and doing the plumbing between them
yourself: making the intermediate directories, moving ``final_pure_reads_*`` into a
``CLEAN_READS`` folder, concatenating them for a co-assembly, passing the right bin directory
to the next module. The project's own tutorial contains those steps as shell one-liners. Every
one of them is a place to make a mistake, none of it is recorded in provenance, and none of it
is resumable.

This command owns that plumbing. You describe the inputs; it decides the layout, runs the
steps in order, runs independent samples concurrently, and skips whatever is already done.

    metawrap2 run samples.toml -o study/ -t 24
    metawrap2 run samples.toml -o study/ --dry-run     # print the plan, run nothing
    metawrap2 run samples.toml -o study/ --resume      # continue where it stopped

**Design.** The driver never reimplements a module: each step is one ordinary
``metawrap2 <module>`` invocation, built by :func:`plan_steps`, so anything the driver can do
you can also do by hand, and the commands it ran are in the study's ``run_commands.txt``.
Step-level state lives in ``<output>/.metawrap2/run_steps`` so ``--resume`` is decided by
whether a step finished, not by guessing from output files.

The sample sheet is TOML (or YAML, if PyYAML happens to be installed - TOML needs no
dependency, since the config loader already parses it):

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

#: Steps that only mean anything when each sample is assembled separately. `coverage` builds one
#: depth table across a cohort's assemblies, and `dereplicate` collapses the per-sample bin sets into
#: one catalogue - a co-assembly has one assembly and one bin set, so neither has anything to do.
COHORT_ONLY_STEPS = frozenset({"coverage", "dereplicate"})

#: The steps run by default: the genome-recovery path. The three that need the largest
#: databases (kraken2, classify_bins, blobology) are opt-in via `steps`, so a default run
#: does not demand a ~300 GB BLAST database and a ~90 GB Kraken2 database.
DEFAULT_STEPS: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "binning",
    "bin_refinement",
    "reassemble_bins",
    "quant_bins",
)

#: The steps a cohort runs when the sheet does not say: the same path, plus the two stages that only
#: exist once each sample has its own assembly. `reassemble_bins` is left out rather than forgotten -
#: it is the most expensive step in the pipeline and it runs per sample here, so a twenty-sample
#: cohort would pay for it twenty times; add it to `steps` deliberately.
COHORT_STEPS: Tuple[str, ...] = (
    "read_qc",
    "assembly",
    "coverage",
    "binning",
    "bin_refinement",
    "dereplicate",
    "quant_bins",
)

#: The binners `binning` runs when the sheet does not say. Two, not three: CONCOCT came out on
#: 2026-10-02 on cost. On cami2_gi it was 192 of the 213 cpu-hours a three-binner run spent binning,
#: and the 10 extra clean MAGs it bought - 74 against 64, where clean means passing GUNC - are reached
#: for a fraction of that by the learned binners (105 for metaBAT2+SemiBin2). It is CONCOCT's price
#: against those, not its accuracy, that retired it. It remains available as --concoct, and
#: bin_refinement still takes three
#: bin sets; what changed is what runs when nobody says. See benchmarks/configurations.toml for the
#: arms, and note the honest caveat recorded there: this pair is also the dirtiest of the three
#: consolidated sets measured (18.9% chimeric), so it is a throughput choice made knowingly.
#: The binning flags a sheet gets when it names none, derived from the binning module's own
#: default tier rather than restated - a literal pair here is how `--preset` and the driver's idea
#: of "the default" drifted apart once before.
DEFAULT_BINNERS: Tuple[str, ...] = tuple(
    "--" + binner for binner in _binning.PRESETS[_binning.DEFAULT_PRESET]
)

#: Modules that accept -m/--memory in GB. The driver divides its budget between concurrent
#: steps and passes each one its share.
MEMORY_MODULES = frozenset({"assembly", "binning", "bin_refinement", "reassemble_bins"})

#: Per-sample steps; everything else operates on the study as a whole.
PER_SAMPLE_STEPS = frozenset({"read_qc"})

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
    #: True when the steps were named rather than defaulted. The default depends on the assembly
    #: mode, and the mode can still change after the sheet is read (`--no-coassembly`), so the
    #: driver has to be able to tell "the user asked for exactly these" from "nobody said".
    steps_chosen: bool = False

    def extra(self, module: str) -> List[str]:
        """Extra command-line options the sheet specified for *module*."""
        return list(self.module_options.get(module, []))

    def default_steps(self) -> Tuple[str, ...]:
        """The steps to run when the sheet did not say, which differs by assembly mode."""
        return DEFAULT_STEPS if self.coassemble else COHORT_STEPS


def resolve_steps(study: Study) -> None:
    """Settle *study*'s step list now that its assembly mode is known.

    Called again after the command line has had its say, because `--no-coassembly` changes which
    default applies; it is idempotent for a study whose steps were named explicitly.
    """
    if not study.steps_chosen:
        study.steps = study.default_steps()
    if study.coassemble:
        dropped = [s for s in study.steps if s in COHORT_ONLY_STEPS]
        if dropped:
            # Running them would be worse than skipping them: `coverage` would build a catalogue
            # of one assembly and `dereplicate` would compare a bin set against itself, both at
            # full cost. Say so rather than doing it.
            warning(
                "Co-assembling, so %s %s nothing to do (one assembly, one bin set); skipping. "
                "Set coassemble = false to run a cohort per sample."
                % (
                    ", ".join(dropped),
                    "has" if len(dropped) == 1 else "have",
                )
            )
            study.steps = tuple(s for s in study.steps if s not in COHORT_ONLY_STEPS)


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
    resolve_steps(study)
    return study


def study_from_directory(directory: str) -> Study:
    """Build a study by pairing up every read file in *directory*."""
    files = [os.path.join(directory, f) for f in sorted(os.listdir(directory))]
    pairs = _reads.collect_pairs(files)  # raises with an actionable message if none
    return Study(samples=[Sample(name=name, r1=r1, r2=r2) for name, r1, r2 in pairs])


EXAMPLE_SHEET = """\
# MetaWrap2 sample sheet. Run it with:  metawrap2 run samples.toml -o study/ -t 24

[settings]
# Co-assemble all samples together (true) or assemble each sample separately (false).
coassemble = true
# Which steps to run, in this order. The default omits kraken2/classify_bins/blobology
# because they need the very large Kraken2 and BLAST databases.
steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]

# Extra options passed through to a module, exactly as you would type them. The binners here are
# the defaults, spelled out so they are easy to change; adding --concoct roughly quadruples what
# binning costs and is worth it when the catalogue matters more than the clock.
[modules]
binning = "--metabat2 --maxbin2"
bin_refinement = "-c 50 -x 10"

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


# --- the layout the driver owns -----------------------------------------------------------


class Layout:
    """Where everything lives inside the study directory.

    One object so the paths are stated once. Anything the driver hands to a module, or expects
    back from one, is named here rather than rebuilt at each call site.
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
        """The depth-table slice `coverage` writes for one sample, named as that module names it."""
        return os.path.join(self.coverage_dir, _coverage.DEPTH_DIR, "%s.txt" % sample)

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

    #: Top-level entries the driver itself owns, so the end-of-run audit does not report them as
    #: written by nobody. Anything else appearing here really was produced by a step that did not
    #: declare it.
    DRIVER_OWNED = ("LOGS", "CLEAN_READS", "manifest.json", "run_config.json", "provenance.txt")


# --- steps -------------------------------------------------------------------------------


@dataclass
class Step:
    """One module invocation the driver will make."""

    name: str  # unique step id, e.g. "read_qc:ERR011347"
    module: str
    argv: List[str]
    parallel_group: Optional[str] = None  # steps sharing a group may run concurrently
    #: What this step reads and what it must produce. Declared here so --resume can verify a
    #: skipped step's outputs still exist and its inputs have not changed, rather than trusting
    #: a marker file (see metawrap2.manifest).
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
        # Modules that take -m must be told their *share* of the budget, not all of it: N
        # concurrent steps each given the full -m can request N times the machine's RAM.
        if memory_gb is not None and self.module in MEMORY_MODULES and "-m" not in self.argv:
            argv += ["-m", str(memory_gb)]
        if config:
            argv += ["--config", config]
        return argv


def _qc_output_reads(layout: Layout, sample: Sample) -> List[str]:
    """Where read_qc's cleaned reads end up for *sample* (before the driver renames them)."""
    out = layout.read_qc(sample.name)
    if sample.layout == "single":
        return [os.path.join(out, "final_pure_reads_1.fastq")]
    return [
        os.path.join(out, "final_pure_reads_1.fastq"),
        os.path.join(out, "final_pure_reads_2.fastq"),
    ]


def clean_reads_for(layout: Layout, study: Study, sample: Sample, ran_read_qc: bool) -> List[str]:
    """The reads later steps should use for *sample*: QC'd if read_qc ran, else the input."""
    if not ran_read_qc:
        return list(sample.files)
    if sample.layout == "single":
        return [layout.clean_single(sample.name)]
    return list(layout.clean_pair(sample.name))


def plan_steps(study: Study, layout: Layout) -> List[Step]:
    """Turn a study into the ordered list of module invocations to make.

    Deliberately pure: it touches no files, so ``--dry-run`` shows exactly the plan that a real
    run would execute, and the planning is straightforward to test.
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
    """The bin folder later steps should read in a co-assembly study: the latest one produced."""
    comp, cont = _quality_thresholds(study)
    if "reassemble_bins" in study.steps:
        return layout.reassembled_bins
    if "bin_refinement" in study.steps:
        return layout.refined_bins(comp, cont)
    return layout.binner_bins(selected_binners(binning_opts)[0])


def _sample_final_bins(
    study: Study, layout: Layout, sample: str, binning_opts: Sequence[str]
) -> str:
    """The same question for one sample of a cohort, answered from the same list of steps."""
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
    """Everything after the assembly when all the samples were assembled together.

    One assembly, one bin set, one of each downstream step - which is the whole reason this is
    separate from the cohort plan rather than a flag inside one.
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
        # Hand refinement only the bin sets binning was actually asked to produce. Passing a
        # fixed -A/-B/-C meant a study that selected, say, "--metabat2 --concoct" failed here
        # with "INITIAL_BINNING/maxbin2_bins is not a valid directory". Repeated --bins rather
        # than -A/-B/-C because the letters stop at F and zipping against three of them silently
        # dropped the fourth bin set; --bins is the form bin_refinement documents for a caller
        # generating its arguments, which is what this is.
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

    # Without refinement, fall back to the first binner that was actually run.
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
    """classify_bins and annotate_bins, which care only about a folder of bins.

    Shared between the two modes because they are the two steps where a cohort and a co-assembly
    differ only in which folder is handed over - the catalogue, or the one bin set there was.
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
    """Everything after the assemblies when each sample was assembled on its own.

    This is the cohort pipeline of the plan: one coverage pass across every assembly, per-sample
    binning on a slice of that one depth table, per-sample consolidation, then one dereplication into
    a species catalogue that the quantifier and the annotators see instead of any single sample's
    bins.

    Each per-sample stage is its own parallel group, so a cohort runs N samples at a time through a
    stage rather than one sample end to end - which is what keeps the thread budget divided by the
    driver rather than by luck.
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

    # One coverage pass for the whole cohort: every sample's reads against the concatenated
    # catalogue of every sample's assembly. N alignments, not N per assembly - which is the only
    # reason multi-sample binning is affordable at cohort scale.
    cohort_coverage = "coverage" in study.steps
    if cohort_coverage:
        argv: List[str] = []
        for sample in samples:
            argv += ["-a", "%s=%s" % (sample, layout.sample_assembly(sample))]
        argv += ["-o", layout.coverage_dir]
        steps.append(
            Step(
                name="coverage",
                module="coverage",
                argv=argv + study.extra("coverage") + all_reads,
                inputs=[layout.sample_assembly(s) for s in samples] + all_reads,
                outputs=[layout.coverage_depth(s) for s in samples],
            )
        )

    if "binning" in study.steps:
        for sample in samples:
            argv = ["-a", layout.sample_assembly(sample), "-o", layout.sample_binning_dir(sample)]
            if cohort_coverage:
                # The point of the coverage stage: this sample's contigs, every sample's depths,
                # and nothing aligned here at all.
                argv += ["--depth", layout.coverage_depth(sample)]
                inputs = [layout.sample_assembly(sample), layout.coverage_depth(sample)]
            else:
                # No coverage stage means single-sample coverage, which is a weaker binning and is
                # what the user asked for by dropping the step.
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
    # Without dereplication there is no one folder of bins for the cohort, so the steps that want
    # one get the first sample's - which is wrong for a cohort and is why dereplicate is a default.
    final_bins = catalogue if "dereplicate" in study.steps else per_sample_bins[0]

    if "quant_bins" in study.steps:
        # Against the catalogue, with no -a: near-identical MAGs from different samples would
        # otherwise steal each other's reads, and there is no single assembly to point at anyway.
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
    """Which binners a set of `binning` options turns on, in the order refinement is given them.

    The order comes from the binning module's own ``BINNERS``, so a binner added there is wired
    here without a second edit. It used to be a literal list of the three classical binners, which
    meant a sheet saying ``--metabat2 --semibin2`` had SemiBin2's bins dropped before refinement,
    and a sheet saying ``--semibin2`` alone fell through to the default pair and pointed refinement
    at two directories that were never created.

    ``--preset`` is resolved from the binning module's own table, with the same precedence that module
    applies: an explicitly named binner wins, and the preset is only consulted when none was. Without
    this, a sheet saying ``--preset very-fast`` matched no ``--<binner>`` flag, fell through to the
    default pair, and failed at refinement with "maxbin2_bins is not a directory" - after paying for
    the binning. Every preset but ``fast`` was unusable from a study sheet.
    """
    selected = [b for b in _binning.BINNERS if "--" + b in binning_options]
    if selected:
        return selected
    if "--preset" in binning_options:
        index = binning_options.index("--preset") + 1
        if index < len(binning_options):
            preset = _binning.PRESETS.get(binning_options[index])
            if preset:
                return [b for b in _binning.BINNERS if b in preset]
    return [b.lstrip("-") for b in DEFAULT_BINNERS]


def _quality_thresholds(study: Study) -> Tuple[int, int]:
    """The -c/-x the refinement step will use, so we can name its output directory.

    bin_refinement names its output ``metawrap_<c>_<x>_bins``, so the driver has to know the
    thresholds to hand that directory to the next step.
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
    """Decides what a resume may skip, backed by the study's manifest.

    The manifest records each step's command, timing, declared inputs and outputs with content
    fingerprints, and its outcome. ``--resume`` consults it rather than a marker file, so a
    step re-runs when its outputs have gone missing or its inputs have changed - both of which
    a marker file cannot tell you.
    """

    def __init__(self, layout: Layout, resume: bool, strict_fingerprints: bool = False):
        self.manifest = Manifest(os.path.join(layout.state_dir, MANIFEST_NAME))
        self.resume = resume
        self.strict = strict_fingerprints

    def skip_reason(self, step: Step) -> Optional[str]:
        """Why *step* may be skipped, or None if it must run."""
        if not self.resume:
            return None
        ok, reason = self.manifest.can_skip(step.name, whole_file=self.strict)
        return reason if ok else None

    def rerun_reason(self, step: Step) -> str:
        """Why a step is being re-run despite --resume (for an explanatory message)."""
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


# --- the driver ---------------------------------------------------------------------------


def _link_or_copy(src: str, dest: str) -> None:
    """Point *dest* at *src* without duplicating the data when we can help it."""
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

    The tutorial does this with ``mv``; doing it with symlinks instead means the QC output
    directory stays intact, so its provenance record still describes files that exist and a
    resumed run can still find them.
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
                    # copyfileobj streams, so this never holds a library in memory
                    import shutil

                    shutil.copyfileobj(fh, out, length=1024 * 1024)


def _is_current(dest: str, sources: List[str]) -> bool:
    """True if *dest* exists and is newer than every source (so it need not be rebuilt)."""
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
    ``completed`` / ``failed`` / ``interrupted``. The distinction matters for ``--resume``: an
    interrupted step has outputs in an unknown state and must always re-run, whereas a failed
    step may have a diagnosable reason in its log.
    """
    os.makedirs(layout.log_dir, exist_ok=True)
    log_path = os.path.join(layout.log_dir, "%s.log" % step.name.replace(":", "_"))
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
            # Measured rather than merely timed: the step's real peak memory is the one number
            # that would let the memory budgets and space multipliers be calibrated instead of
            # guessed. See metawrap2.usage.
            rc, usage = _usage.run_measured(argv, stdout=log, stderr=subprocess.STDOUT)
            log.write("\n[metawrap2] %s\n" % _usage.summarise(usage))
        status = COMPLETED if rc == 0 else FAILED
    except KeyboardInterrupt:
        # Ctrl-C leaves this step's outputs half-written. Recording that - rather than letting
        # the process die with no trace - is what lets --resume know to redo it instead of
        # trusting whatever happens to be on disk.
        status = INTERRUPTED
        raise
    finally:
        if state is not None:
            state.finish(step, started, status, log=log_path, usage=usage)
    return step, status, log_path


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
                    # Explain why a resume is redoing something, so "it ran this again" is
                    # never a mystery.
                    comm("re-running %s - %s" % (s.name, state.rerun_reason(s)))
                pending.append(s)

        if not pending:
            _after_group(group[0], layout, study, dry_run)
            continue

        _before_group(group[0], layout, study, dry_run)

        workers = max(1, min(jobs, len(pending)))
        per_step_threads = max(1, threads // workers)
        # Divide the memory budget the same way as the threads: N concurrent SPAdes runs each
        # told they may use the whole budget will between them ask for N times the machine's RAM.
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
                    else:
                        failures.append((s, log_path))
        except KeyboardInterrupt:
            # Every step still running has been recorded as interrupted by run_step's finally,
            # so a later --resume redoes exactly those and keeps what genuinely finished.
            warning(
                "Interrupted. The steps that were running are recorded as interrupted; "
                "re-run the same command with --resume to continue - they will be redone, "
                "and finished steps will be skipped."
            )
            return 130

        if failures:
            for s, log_path in failures:
                warning("%s failed. Its log is %s" % (s.name, log_path))
            error(
                "%d of %d step(s) failed; stopping here so the failure is not compounded.\n"
                "Fix the cause and re-run the same command with --resume to carry on from "
                "this point." % (len(failures), len(pending))
            )

        _after_group(group[0], layout, study, dry_run)

    announcement("STUDY COMPLETE - %d step(s) finished" % completed)
    if not dry_run:
        audit_outputs(state, layout)
    comm("results are under %s" % layout.root)
    return 0


def audit_outputs(state: StepState, layout: Layout) -> None:
    """Check at the end of a run that the manifest describes what is actually on disk.

    The manifest is what ``--resume`` trusts, so a claim in it that does not match reality is a
    problem now rather than on the next run - when the cause will be a run in the past. This
    warns instead of failing: every step reported success, and refusing to finish a completed
    study over a bookkeeping discrepancy would be worse than saying so clearly.
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
    """The plumbing a step needs in place before it runs.

    One case, and it was a real hole: a co-assembly reads ``ALL_READS_1.fastq``, which is built by
    concatenating the samples, and that only happened after a ``read_qc`` step. A study that
    co-assembles *without* read QC - every arm of a benchmark that holds QC constant, for instance -
    therefore planned an assembly whose inputs nothing created, and failed on "does not exist" naming
    the files but not the cause.

    Building it here instead of after read_qc covers both: with QC the files are already current by the
    time assembly starts and this is a no-op, and without it the raw reads get concatenated.
    ``clean_reads_for`` already knew how to answer "which reads, QC'ed or not" - it was simply never
    asked in this case.
    """
    if step.module == "assembly" and study.coassemble:
        announcement("CONCATENATING READS FOR THE CO-ASSEMBLY")
        concatenate_for_coassembly(layout, study, dry_run)


def _after_group(step: Step, layout: Layout, study: Study, dry_run: bool) -> None:
    """The plumbing between steps: exactly the ``mv``/``cat`` the tutorial asks users to do."""
    if step.module == "read_qc":
        announcement("COLLECTING THE QC'ED READS")
        collect_clean_reads(layout, study, dry_run)
    if step.module == "read_qc" and study.coassemble:
        concatenate_for_coassembly(layout, study, dry_run)


# --- entry point --------------------------------------------------------------------------


def main(argv: List[str], prog: str = "metawrap2 run", cohort: bool = False) -> int:
    """Plan and execute a study.

    *cohort* is what `metawrap2 cohort` passes: the same planner with per-sample assembly forced on,
    so there is one driver to keep correct rather than two.
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
        print(EXAMPLE_SHEET)
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
        # A malformed sample sheet is a user error, so report it as one rather than as a
        # traceback (TypeError included: a non-mapping sample entry raises it).
        error(str(exc))

    if args.steps:
        study.steps = tuple(args.steps)
        study.steps_chosen = True
    if args.no_coassembly or cohort:
        if cohort and study.coassemble_chosen and study.coassemble:
            # The sheet asked for a co-assembly and the command asked for a cohort. The command
            # wins, because `metawrap2 cohort` is nothing but per-sample assembly - but say so,
            # rather than quietly ignoring a line the user wrote.
            warning(
                "the sheet sets coassemble = true, but `metawrap2 cohort` assembles each sample "
                "separately; running per-sample. Use `metawrap2 run` to honour the sheet."
            )
        study.coassemble = False
    # After the mode is final, not before: --no-coassembly changes which default step list applies,
    # and --skip has to subtract from whichever one won.
    resolve_steps(study)
    if args.skip:
        study.steps = tuple(s for s in study.steps if s not in args.skip)
    if not study.steps:
        error("no steps left to run.")

    # metawrap2.cli strips the global --dry-run/--force/--resume out of argv and applies them
    # to the shared runner before this subcommand parses anything, so accept them from either
    # place: `metawrap2 --dry-run run ...` and `metawrap2 run ... --dry-run` must behave the
    # same. Without this the driver silently *executed* a run the user asked to only plan.
    from .. import command as _command

    dry_run = bool(args.dry_run or _command.runner.dry_run)
    resume = bool(args.resume or _command.runner.resume)
    force = bool(args.force or _command.runner.force)

    settings = load_settings(args.config)
    # The driver is the one process that knows how many steps will run at once, so it is the one
    # that has to divide the budget. Resolving it here (rather than letting each module resolve
    # the whole machine for itself) is what keeps N concurrent steps from each claiming 80% of
    # the RAM between them. See metawrap2.limits.
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
    comm("Resource budget: %s" % budget.describe())
    overrides = _limits.overrides()
    if overrides:
        comm(
            "Budget overridden on the command line: %s"
            % ", ".join("--%s %s" % (k.replace("_", "-"), v) for k, v in sorted(overrides.items()))
        )

    state = StepState(layout, resume=resume, strict_fingerprints=args.strict_fingerprints)
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
