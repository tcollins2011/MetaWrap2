"""The benchmark dataset and comparison registry, read from ``benchmarks/*.toml``.

Datasets are declared as data, not as a fetch script, so the URL, pinned checksum, measured
download size and gold-standard layout sit together and can be reviewed. The loader rejects a
dataset missing a checksum or a measured size at load time, before any download starts.

* **No network access.** Importing, listing or planning a fetch never touches the network; only
  :func:`metawrap2.benchmark.fetch.fetch` does, for a dataset named by the user.
* **No size estimates.** Every ``bytes`` value in the file was measured with a HEAD request, so
  the disk preflight uses real figures.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml

#: Gold-standard roles the harness understands. ``assembly`` is the reference contigs, ``binning``
#: the contig-to-genome truth a binning scorer needs, ``profile`` a taxonomic abundance truth.
#: Unknown roles are rejected, since a mistyped role would score against the wrong file.
GOLD_STANDARD_ROLES = {
    "assembly": "reference contigs the gold standard bins were derived from",
    "binning": "contig -> genome truth that binning is scored against",
    "binning_raw": "uncorrected variant of the binning truth, kept for comparison",
    "profile": "taxonomic abundance truth for profiling comparisons",
}

#: How reads are laid out inside a dataset, which decides whether the harness has to split them.
READS_LAYOUTS = {
    "paired": "separate _1/_2 files per sample",
    "interleaved": "one file per sample with mates adjacent; deinterleaved on input",
    "single": "unpaired reads",
}


#: The placeholders an external arm's command template may use, and what the runner fills each
#: from. An arm cannot name inputs of its own, so it runs on the same inputs and samples as the
#: native arms.
EXTERNAL_PLACEHOLDERS = {
    "output": "the arm's own output directory, created before the command starts",
    "threads": "the budget's thread grant - the same number every arm of the comparison gets",
    "memory": "the budget's memory grant, in GB",
    "assembly": "the contigs the comparison bins (the gold-standard one unless --assembly said so)",
    "reads": "every selected sample's read files, space separated",
    "r1": "the forward reads, for a module that takes one sample per invocation",
    "r2": "the reverse reads",
    "bins": "the one bin set --bins-from names",
    "sample": "the sample label, for a per-sample module",
}

#: The placeholders an external arm must use for each module. ``output`` and ``threads`` are in
#: every entry, so the runner controls where output goes and every arm gets the same thread grant.
#: The rest are the module's inputs. Requiring them rejects a command that reads data from a path
#: typed into it.
EXTERNAL_REQUIRED_PLACEHOLDERS = {
    "read_qc": ("output", "threads", "r1", "r2"),
    "assembly": ("output", "threads", "r1", "r2"),
    "binning": ("output", "threads", "assembly", "reads"),
    "bin_refinement": ("output", "threads", "bins"),
    "quant_bins": ("output", "threads", "bins", "reads"),
}

#: The required placeholders for a module not listed above: the output directory and the thread
#: grant.
EXTERNAL_REQUIRED_FALLBACK = ("output", "threads")

#: Commands that return before the work they start has finished.
#: :func:`metawrap2.usage.run_measured` reaps the process it started, so a submitted or backgrounded
#: arm would be measured as the submission: a few seconds of wall time and almost no CPU. They are
#: rejected at load time.
UNMEASURABLE_WORDS = (
    "nohup",
    "setsid",
    "disown",
    "sbatch",
    "qsub",
    "bsub",
    "screen",
    "tmux",
    "systemd-run",
)

#: A placeholder in a command template: a lowercase name with no spaces, so shell and awk braces
#: (``${HOME}``, ``{print $1}``) are not taken for one. A command that needs a literal ``{word}``
#: has to use a wrapper script.
PLACEHOLDER_RE = re.compile(r"\{([a-z_][a-z0-9_]*)\}")


class RegistryError(Exception):
    """A dataset or comparison declaration that cannot be used as written."""


@dataclass(frozen=True)
class Remote:
    """One file belonging to a dataset.

    *md5* and *size* are optional only for local fixtures. For anything with a ``base_url``
    :func:`_check_dataset` requires both.
    """

    path: str
    md5: Optional[str] = None
    size: Optional[int] = None
    sample: Optional[int] = None
    role: Optional[str] = None
    note: str = ""

    @property
    def name(self) -> str:
        """The basename, which is the file's name after download."""
        return os.path.basename(self.path)

    @property
    def is_archive(self) -> bool:
        return self.path.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar"))


@dataclass(frozen=True)
class Dataset:
    """A benchmark dataset: where to get it, what it costs, and what the truth is."""

    name: str
    title: str
    description: str = ""
    citation: str = ""
    doi: str = ""
    license: str = ""
    base_url: str = ""
    checksum: str = "md5"
    checksum_manifest: str = ""
    #: Set for a dataset that ships with the repository, as a path relative to the registry file.
    #: Mutually exclusive with *base_url*: a dataset is either fetched or committed.
    local: str = ""
    local_root: str = ""
    samples: Sequence[int] = ()
    keep: Sequence[str] = ()
    reads_glob: str = ""
    reads_layout: str = "paired"
    #: Glob, relative to the extracted data, matching one read-to-genome mapping per sample. The
    #: ``binning`` gold standard maps contigs to genomes; this maps reads, which is what measures an
    #: organism's share of a sample. Needed only for quantifier comparisons.
    abundance_truth_glob: str = ""
    disk_note: str = ""
    files: Sequence[Remote] = ()
    gold_standard: Sequence[Remote] = ()

    @property
    def is_local(self) -> bool:
        return bool(self.local)

    @property
    def download_bytes(self) -> int:
        """Total bytes on the wire for a complete fetch.

        Also used as the disk preflight figure. It is an upper bound: the fetcher deletes each
        archive after extraction, and the extracted content is compressed FASTQ, so it cannot
        exceed the download.
        """
        return sum(f.size or 0 for f in list(self.files) + list(self.gold_standard))

    def url_for(self, remote: Remote) -> str:
        if self.is_local:
            raise RegistryError("%s is a local dataset; it has no URLs." % self.name)
        return self.base_url.rstrip("/") + "/" + remote.path.lstrip("/")

    def local_path_for(self, remote: Remote) -> str:
        """Where *remote* lives for a committed dataset."""
        if not self.is_local:
            raise RegistryError("%s is not a local dataset." % self.name)
        return os.path.join(self.local_root, remote.path)

    def gold(self, role: str) -> Optional[Remote]:
        """The gold-standard file for *role*, or None when this dataset has no such truth."""
        for remote in self.gold_standard:
            if remote.role == role:
                return remote
        return None

    def files_for_sample(self, sample: int) -> List[Remote]:
        return [f for f in self.files if f.sample == sample]

    def describe(self) -> str:
        """One line for ``metawrap2 benchmark list``."""
        from ..limits import human

        if self.is_local:
            where = "committed fixture"
        else:
            where = "%s to download" % human(self.download_bytes)
        return "%-12s %-4s %s  (%s)" % (
            self.name,
            "%dx" % len(self.samples) if self.samples else "-",
            self.title,
            where,
        )


@dataclass(frozen=True)
class External:
    """Another pipeline run in place of MetaWrap2 in one arm of a comparison, such as metaWRAP 1.

    The arm is declared in the registry, so it gets the same inputs, samples, budget and cost
    measurement as the native arms.

    * *command*: the command template. :data:`EXTERNAL_PLACEHOLDERS` stand for the inputs the
      runner resolved; the arm cannot name inputs of its own.
    * *bins*: where its bins are written under its output directory. The runner publishes them as
      ``<arm name>_bins`` so the scorer finds them.
    * *env*: the conda env its tools are in. metaWRAP 1.3.0 needs CheckM 1.0.18 and Python 2.7.
    * *version* or *version_command*: one is required, and the runner records the version.
    """

    command: str
    bins: str
    env: str
    version: str = ""
    version_command: str = ""

    @property
    def placeholders(self) -> Tuple[str, ...]:
        """The placeholder names *command* uses, in the order they first appear."""
        seen: List[str] = []
        for name in PLACEHOLDER_RE.findall(self.command):
            if name not in seen:
                seen.append(name)
        return tuple(seen)


@dataclass(frozen=True)
class Variant:
    """One option in a comparison: a name, the flags that select it, and whether it ships.

    *external* is set instead of *args* for an arm that runs another pipeline. The loader rejects
    a variant that sets both.
    """

    name: str
    args: Sequence[str] = ()
    default: bool = False
    note: str = ""
    external: Optional[External] = None

    @property
    def is_external(self) -> bool:
        return self.external is not None


@dataclass(frozen=True)
class Comparison:
    """One decision the benchmark settles, and the options for it."""

    name: str
    module: str
    question: str = ""
    metrics: Sequence[str] = ()
    #: Flags every variant gets, for the parts of the command not under test (`--skip-bmtagger`
    #: on synthetic data with no host, for example). Declared once so all variants get the same.
    common_args: Sequence[str] = ()
    variants: Sequence[Variant] = field(default_factory=tuple)

    @property
    def incumbent(self) -> Optional[Variant]:
        """The variant MetaWrap2 currently ships, for showing beside the alternatives."""
        for variant in self.variants:
            if variant.default:
                return variant
        return None

    def variant(self, name: str) -> Variant:
        for variant in self.variants:
            if variant.name == name:
                return variant
        raise RegistryError(
            "%s has no variant %r. Known: %s"
            % (self.name, name, ", ".join(v.name for v in self.variants))
        )


def registry_dir() -> str:
    """The ``benchmarks/`` directory.

    Resolved from the source tree, since benchmarks run only from a checkout (the registry
    describes 100 GB downloads and committed fixtures, which do not belong in a wheel).
    ``METAWRAP2_BENCHMARKS`` overrides it; the tests use that.
    """
    override = os.environ.get("METAWRAP2_BENCHMARKS")
    if override:
        return override
    # src/metawrap2/benchmark/registry.py -> the repository root -> benchmarks/
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(here))), "benchmarks")


def datasets_path() -> str:
    return os.path.join(registry_dir(), "datasets.toml")


def configurations_path() -> str:
    return os.path.join(registry_dir(), "configurations.toml")


def _read_toml(path: str) -> Dict[str, Any]:
    if not os.path.isfile(path):
        raise RegistryError(
            "No benchmark registry at %s.\n"
            "Benchmarks are run from a source checkout; set METAWRAP2_BENCHMARKS to point at a "
            "benchmarks/ directory if yours is elsewhere." % path
        )
    try:
        with open(path, "rb") as fh:
            return _toml.load(fh)
    except (OSError, _toml.TOMLDecodeError) as exc:
        raise RegistryError("%s is not readable as TOML: %s" % (path, exc))


def _remote(entry: Dict[str, Any], where: str, *, gold: bool = False) -> Remote:
    if not isinstance(entry, dict):
        raise RegistryError("%s: expected a table, got %r" % (where, entry))
    path = entry.get("path")
    if not path or not isinstance(path, str):
        raise RegistryError("%s: every file needs a 'path'." % where)
    if path.startswith("/") or ".." in path.split("/"):
        # Registry paths decide where files are written; absolute or escaping paths would write
        # outside the user's destination.
        raise RegistryError(
            "%s: path %r must be relative and must not contain '..'." % (where, path)
        )
    role = entry.get("role")
    if gold:
        if not role:
            raise RegistryError("%s: a gold_standard entry needs a 'role'." % where)
        if role not in GOLD_STANDARD_ROLES:
            raise RegistryError(
                "%s: unknown gold-standard role %r. Known roles: %s"
                % (where, role, ", ".join(sorted(GOLD_STANDARD_ROLES)))
            )
    size = entry.get("bytes")
    if size is not None and (not isinstance(size, int) or size <= 0):
        raise RegistryError("%s: 'bytes' must be a positive integer, got %r" % (where, size))
    return Remote(
        path=path,
        md5=entry.get("md5"),
        size=size,
        sample=entry.get("sample"),
        role=role,
        note=entry.get("note", ""),
    )


def _external(entry: Dict[str, Any], where: str, module: str) -> External:
    """One external arm, or a refusal naming the arm and the problem.

    These checks need nothing to run. Without them the mistakes would show up as wrong results
    rows: no bins path gives "nothing to score", a backgrounded command gives a few seconds of CPU,
    and an input path typed into the command gives a row measured on different data.
    """
    if not isinstance(entry, dict):
        raise RegistryError("%s: 'external' should be a table, got %r" % (where, entry))
    command = entry.get("command") or ""
    if not isinstance(command, str) or not command.strip():
        raise RegistryError(
            "%s: an external arm needs a 'command', the template the pipeline is run from." % where
        )
    bins = entry.get("bins") or ""
    if not isinstance(bins, str) or not bins.strip():
        raise RegistryError(
            "%s: an external arm needs 'bins', the path under {output} where its bins are "
            "written. Without it the scorer finds no bins." % where
        )
    bins = bins.strip()
    if bins.startswith("/") or ".." in bins.split("/"):
        raise RegistryError(
            "%s: 'bins' must be relative to {output} and must not contain '..'; got %r."
            % (where, bins)
        )
    from ..constants import is_refined_bins_dirname

    if is_refined_bins_dirname(bins.split("/")[0]):
        # The scorer finds a refined bin set one level under an arm's output by itself, and the
        # runner also publishes the declared directory under the arm's name. With both, the same
        # bins would be assessed twice.
        raise RegistryError(
            "%s: 'bins' starts with %r, a name the scorer finds by itself, so these bins would be "
            "scored twice (directly and through the <arm>_bins link). Put the pipeline's output "
            "one level down." % (where, bins.split("/")[0])
        )
    env = entry.get("env") or ""
    if not isinstance(env, str) or not env.strip():
        raise RegistryError(
            "%s: an external arm needs 'env', the conda environment with the pipeline's tools "
            "(metaWRAP 1.3.0, for example, needs Python 2.7)." % where
        )
    if not (entry.get("version") or entry.get("version_command")):
        raise RegistryError(
            "%s: an external arm needs 'version' or 'version_command', so its results record "
            "which version ran." % where
        )
    for word in command.split():
        bare = word.strip("\"'();")
        if bare in UNMEASURABLE_WORDS or bare == "&" or (bare.endswith("&") and bare != "&&"):
            raise RegistryError(
                "%s: 'command' contains %r, which returns before the work finishes, so the runner "
                "would measure the submission instead of the run. Run the pipeline in the "
                "foreground." % (where, bare)
            )
    external = External(
        command=command.strip(),
        bins=bins,
        env=env.strip(),
        version=entry.get("version", ""),
        version_command=entry.get("version_command", ""),
    )
    unknown = [name for name in external.placeholders if name not in EXTERNAL_PLACEHOLDERS]
    if unknown:
        raise RegistryError(
            "%s: 'command' uses %s, which the runner cannot fill. It knows:\n%s"
            % (
                where,
                ", ".join("{%s}" % name for name in unknown),
                "\n".join(
                    "  {%-10s %s" % (name + "}", EXTERNAL_PLACEHOLDERS[name])
                    for name in sorted(EXTERNAL_PLACEHOLDERS)
                ),
            )
        )
    required = EXTERNAL_REQUIRED_PLACEHOLDERS.get(module, EXTERNAL_REQUIRED_FALLBACK)
    missing = [name for name in required if name not in external.placeholders]
    if missing:
        raise RegistryError(
            "%s: an external %s arm's command must use %s, so that it gets the same inputs, "
            "samples and thread grant as the other arms."
            % (where, module, ", ".join("{%s}" % name for name in missing))
        )
    return external


def _check_dataset(dataset: Dataset) -> None:
    """Reject declarations that would otherwise fail later."""
    where = "dataset %r" % dataset.name
    if dataset.is_local and dataset.base_url:
        raise RegistryError("%s: set either 'local' or 'base_url', not both." % where)
    if not dataset.is_local and not dataset.base_url:
        raise RegistryError("%s: needs a 'base_url' (or 'local' for a committed fixture)." % where)
    if not dataset.files:
        raise RegistryError("%s: declares no files." % where)
    if dataset.reads_layout not in READS_LAYOUTS:
        raise RegistryError(
            "%s: unknown reads_layout %r. Known: %s"
            % (where, dataset.reads_layout, ", ".join(sorted(READS_LAYOUTS)))
        )
    if not dataset.is_local:
        if dataset.checksum != "md5":
            raise RegistryError(
                "%s: checksum %r is not supported; only 'md5' is (the checksum the dataset "
                "publishers provide)." % (where, dataset.checksum)
            )
        for remote in list(dataset.files) + list(dataset.gold_standard):
            if not remote.md5:
                raise RegistryError(
                    "%s: %s has no md5. Remote files need a checksum so that a truncated "
                    "download is detected." % (where, remote.path)
                )
            if remote.size is None:
                raise RegistryError(
                    "%s: %s has no 'bytes'. The fetch needs the size to check the [limits] "
                    "disk ceiling before it starts; get it from a HEAD request on the URL."
                    % (where, remote.path)
                )
    declared = {f.sample for f in dataset.files if f.sample is not None}
    unknown = declared - set(dataset.samples) if dataset.samples else set()
    if unknown:
        raise RegistryError(
            "%s: files reference sample(s) %s that 'samples' does not list."
            % (where, ", ".join(str(s) for s in sorted(unknown)))
        )
    if dataset.samples:
        missing = set(dataset.samples) - declared
        if missing:
            raise RegistryError(
                "%s: 'samples' lists %s but no file declares them."
                % (where, ", ".join(str(s) for s in sorted(missing)))
            )


def load_datasets(path: Optional[str] = None) -> Dict[str, Dataset]:
    """Every declared dataset, keyed by name. Never touches the network."""
    path = path or datasets_path()
    data = _read_toml(path)
    root = os.path.dirname(os.path.abspath(path))
    schema = data.pop("schema", None)
    if schema != 1:
        raise RegistryError("%s: expected 'schema = 1' at the top, found %r." % (path, schema))
    out: Dict[str, Dataset] = {}
    for name, table in data.items():
        if not isinstance(table, dict):
            raise RegistryError("%s: top-level %r should be a [table]." % (path, name))
        where = "%s: dataset %r" % (path, name)
        local = table.get("local", "")
        dataset = Dataset(
            name=name,
            title=table.get("title") or name,
            description=table.get("description", ""),
            citation=table.get("citation", ""),
            doi=table.get("doi", ""),
            license=table.get("license", ""),
            base_url=table.get("base_url", ""),
            checksum=table.get("checksum", "md5"),
            checksum_manifest=table.get("checksum_manifest", ""),
            local=local,
            local_root=os.path.join(root, local) if local else "",
            samples=tuple(table.get("samples", ())),
            keep=tuple(table.get("keep", ())),
            reads_glob=table.get("reads_glob", ""),
            reads_layout=table.get("reads_layout", "paired"),
            abundance_truth_glob=table.get("abundance_truth_glob", ""),
            disk_note=table.get("disk_note", ""),
            files=tuple(
                _remote(entry, "%s file %d" % (where, i))
                for i, entry in enumerate(table.get("files", ()))
            ),
            gold_standard=tuple(
                _remote(entry, "%s gold_standard %d" % (where, i), gold=True)
                for i, entry in enumerate(table.get("gold_standard", ()))
            ),
        )
        _check_dataset(dataset)
        out[name] = dataset
    if not out:
        raise RegistryError("%s declares no datasets." % path)
    return out


def load_configurations(path: Optional[str] = None) -> Dict[str, Comparison]:
    """Every declared comparison, keyed by name."""
    path = path or configurations_path()
    data = _read_toml(path)
    schema = data.pop("schema", None)
    if schema != 1:
        raise RegistryError("%s: expected 'schema = 1' at the top, found %r." % (path, schema))
    out: Dict[str, Comparison] = {}
    for name, table in data.items():
        if not isinstance(table, dict):
            raise RegistryError("%s: top-level %r should be a [table]." % (path, name))
        module = table.get("module")
        if not module:
            raise RegistryError("%s: comparison %r needs a 'module'." % (path, name))
        variants = []
        seen = set()
        for entry in table.get("variants", ()):
            vname = entry.get("name")
            if not vname:
                raise RegistryError(
                    "%s: comparison %r has a variant with no 'name'." % (path, name)
                )
            if vname in seen:
                raise RegistryError(
                    "%s: comparison %r declares variant %r twice; its results would overwrite "
                    "each other in the TSV." % (path, name, vname)
                )
            seen.add(vname)
            where = "%s: comparison %r variant %r" % (path, name, vname)
            external = None
            if "external" in entry:
                if entry.get("args"):
                    raise RegistryError(
                        "%s: declares both 'args' and an 'external' command. An arm is MetaWrap2 "
                        "with flags or another pipeline, never both." % where
                    )
                if entry.get("default"):
                    raise RegistryError(
                        "%s: an external arm cannot be 'default = true'; the default marks what "
                        "MetaWrap2 ships." % where
                    )
                external = _external(entry["external"], where, module)
            variants.append(
                Variant(
                    name=vname,
                    args=tuple(entry.get("args", ())),
                    default=bool(entry.get("default", False)),
                    note=entry.get("note", ""),
                    external=external,
                )
            )
        if len(variants) < 2:
            raise RegistryError(
                "%s: comparison %r has %d variant(s); it needs at least two."
                % (path, name, len(variants))
            )
        defaults = [v.name for v in variants if v.default]
        if len(defaults) > 1:
            raise RegistryError(
                "%s: comparison %r marks %s as default; only one variant can be the default."
                % (path, name, " and ".join(defaults))
            )
        out[name] = Comparison(
            name=name,
            module=module,
            question=table.get("question", ""),
            metrics=tuple(table.get("metrics", ())),
            common_args=tuple(table.get("common_args", ())),
            variants=tuple(variants),
        )
    if not out:
        raise RegistryError("%s declares no comparisons." % path)
    return out


def get_dataset(name: str, path: Optional[str] = None) -> Dataset:
    datasets = load_datasets(path)
    if name not in datasets:
        raise RegistryError(
            "Unknown dataset %r. Declared: %s" % (name, ", ".join(sorted(datasets)))
        )
    return datasets[name]


def get_comparison(name: str, path: Optional[str] = None) -> Comparison:
    comparisons = load_configurations(path)
    if name not in comparisons:
        raise RegistryError(
            "Unknown comparison %r. Declared: %s" % (name, ", ".join(sorted(comparisons)))
        )
    return comparisons[name]
