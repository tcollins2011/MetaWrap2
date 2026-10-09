"""The benchmark dataset and comparison registry, read from ``benchmarks/*.toml``.

Datasets are declared as data rather than as a fetch script so that adding one is a change a
reviewer can read: the URL, the pinned checksum, the measured size on the wire, and the
gold-standard layout all sit next to each other. The loader's job is to turn that into objects and
to be unforgiving about it - a dataset missing a checksum, or with a size nobody measured, is
rejected at load time rather than three quarters of the way through a 100 GB download.

Two things this module deliberately does not do:

* **It does not fetch anything.** Importing the registry, listing it, or planning a fetch never
  touches the network. Only :func:`metawrap2.benchmark.fetch.fetch` does, and only when a user
  asked for it by name.
* **It does not estimate sizes.** Every ``bytes`` value in the file was measured with a HEAD
  request. An estimated size in a disk pre-flight is worse than no size, because it is believed.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml

#: Gold-standard roles the harness understands. ``assembly`` is the reference contigs, ``binning``
#: the contig-to-genome truth a binning scorer needs, ``profile`` a taxonomic abundance truth.
#: Anything else is a typo, and a typo in a gold-standard role means scoring silently against the
#: wrong file, so unknown roles are refused rather than ignored.
GOLD_STANDARD_ROLES = {
    "assembly": "reference contigs the gold standard bins were derived from",
    "binning": "contig -> genome truth, the input a binning scorer compares against",
    "binning_raw": "an unfixed variant of the binning truth, kept only for comparison",
    "profile": "taxonomic abundance truth for profiling comparisons",
}

#: How reads are laid out inside a dataset, which decides whether the harness has to split them.
READS_LAYOUTS = {
    "paired": "separate _1/_2 files per sample",
    "interleaved": "one file per sample with mates adjacent; deinterleaved on the way in",
    "single": "unpaired reads",
}


class RegistryError(Exception):
    """A dataset or comparison declaration that cannot be used as written."""


@dataclass(frozen=True)
class Remote:
    """One file belonging to a dataset.

    *md5* and *size* are optional only because a local fixture has neither a checksum to pin nor a
    wire size to budget for. For anything with a ``base_url`` they are required, and
    :func:`_check_dataset` enforces that - a remote file without a checksum is a file nobody can
    tell is intact.
    """

    path: str
    md5: Optional[str] = None
    size: Optional[int] = None
    sample: Optional[int] = None
    role: Optional[str] = None
    note: str = ""

    @property
    def name(self) -> str:
        """The basename, which is what the file is called once it has landed."""
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
    #: Mutually exclusive with *base_url*: a dataset is either fetched or committed, never both.
    local: str = ""
    local_root: str = ""
    samples: Sequence[int] = ()
    keep: Sequence[str] = ()
    reads_glob: str = ""
    reads_layout: str = "paired"
    #: Glob, relative to the extracted data, matching one read-to-genome mapping per sample. This is
    #: a different truth from the ``binning`` gold standard: that one says which genome a *contig*
    #: came from, this one which genome a *read* came from, and only the second can say how much of
    #: a sample an organism was. A quantifier comparison needs it; nothing else does, so a dataset
    #: without it simply cannot be used for one.
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

        This doubles as the disk pre-flight figure, and it is a true upper bound rather than a
        guess: the fetcher deletes each archive once extracted, and what it keeps out of a
        ``.tar.gz`` is already-compressed FASTQ, so the extracted footprint cannot exceed the
        download it came from.
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
class Variant:
    """One option in a comparison: a name, the flags that select it, and whether it ships."""

    name: str
    args: Sequence[str] = ()
    default: bool = False
    note: str = ""


@dataclass(frozen=True)
class Comparison:
    """One decision the benchmark exists to settle, and the options for it."""

    name: str
    module: str
    question: str = ""
    metrics: Sequence[str] = ()
    #: Flags every variant gets, for the parts of the command that are not what is under test -
    #: `--skip-bmtagger` on synthetic data with no host, say. Declared once at the comparison level
    #: rather than repeated per variant, because a flag copy-pasted into three variants and
    #: forgotten in the fourth is a difference nobody meant to measure.
    common_args: Sequence[str] = ()
    variants: Sequence[Variant] = field(default_factory=tuple)

    @property
    def incumbent(self) -> Optional[Variant]:
        """The variant MetaWrap2 currently ships, so a table can show it beside its challengers."""
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

    Resolved from the source tree when running from a checkout, which is the only place a
    benchmark is ever run from - the registry describes 100 GB downloads and committed fixtures,
    neither of which belongs in an installed wheel. ``METAWRAP2_BENCHMARKS`` overrides it, which is
    what the tests use.
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
        # A registry is data, but it is data that names where bytes land. An absolute or escaping
        # path would write outside the destination the user chose.
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


def _check_dataset(dataset: Dataset) -> None:
    """Refuse a declaration that would fail late rather than early."""
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
                "%s: checksum %r is not supported; only 'md5' is, because that is what the "
                "publishers of these datasets provide." % (where, dataset.checksum)
            )
        for remote in list(dataset.files) + list(dataset.gold_standard):
            if not remote.md5:
                raise RegistryError(
                    "%s: %s has no md5. A remote file without a checksum is a file nobody can "
                    "tell is intact, and a truncated 10 GB download looks exactly like a good "
                    "one." % (where, remote.path)
                )
            if remote.size is None:
                raise RegistryError(
                    "%s: %s has no 'bytes'. The size is what lets the fetch check the [limits] "
                    "disk ceiling before it starts; HEAD the URL and record it."
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
        raise RegistryError(
            "%s: expected 'schema = 1' at the top, found %r. The schema number is there so a "
            "future format change is a clear error rather than a misread field." % (path, schema)
        )
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
            variants.append(
                Variant(
                    name=vname,
                    args=tuple(entry.get("args", ())),
                    default=bool(entry.get("default", False)),
                    note=entry.get("note", ""),
                )
            )
        if len(variants) < 2:
            raise RegistryError(
                "%s: comparison %r has %d variant(s). A comparison with fewer than two options "
                "is not comparing anything." % (path, name, len(variants))
            )
        defaults = [v.name for v in variants if v.default]
        if len(defaults) > 1:
            raise RegistryError(
                "%s: comparison %r marks %s as default. Only one option can be the incumbent."
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
