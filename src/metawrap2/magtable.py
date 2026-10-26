"""One row per genome, joined from the tables a finished run already wrote.

Completeness and contamination are in ``<bins>.stats.tsv``, the chimerism verdict in GUNC's
report, the lineage in ``classify_bins``' two-column file and per-sample abundance in
``quant_bins``' matrix. This module joins them.

No values are computed. Every cell is copied from a pipeline output, so the GC column stays the
fraction CheckM writes and is not converted to a percentage.

Rules:

* An absent measurement is blank, never zero. :data:`BLANK` is written for every measurement that
  was not made, since a 0 would read as a result.
* Every input is optional and no input drops a genome. A run without ``classify_bins`` gets an
  empty taxonomy column. A genome found only in the abundance matrix still gets a row.

Genome names are joined through :func:`metawrap2.constants.bin_label`, which maps
``bin_001.fasta``, ``bin.1.fa`` and ``bin_001`` to the same genome. A real run uses three of these
forms across the inputs (``bin.1`` in the stats table, ``bin.1.fa`` in the taxonomy file,
``bin_001`` in GUNC's report).
"""

from __future__ import annotations

import csv
import os
from typing import Dict, Iterator, List, Mapping, NamedTuple, Optional, Sequence, Tuple

from . import gunc as _gunc
from . import mimag as _mimag
from .benchmark.scoring import QUALITY_TIERS
from .constants import (
    CONTIGS_SUFFIX,
    STATS_SUFFIX,
    TSV_SUFFIX,
    bin_label,
    is_refined_bins_dirname,
    is_unbinned_name,
)

#: Value of an unmeasured cell: an empty string, which TSV readers treat as missing.
BLANK = ""

#: File name of the table this module writes.
TABLE_NAME = "mag_table" + TSV_SUFFIX

#: Label for a genome scored below every tier, distinct from :data:`BLANK` ("not measured").
BELOW_TIERS = "below_medium_quality"

#: Prefix for ``quant_bins``' abundance columns. It marks the values as abundances (coverage for
#: salmon and the depth route, relative abundance for sylph) and stops a sample named e.g.
#: "genome" or "n50" from colliding with a base column.
ABUNDANCE_PREFIX = "abundance."

#: Leftmost column when the table covers more than one bin set: the sample a row's bin set came
#: from, taken from the stats table's directory as ``metawrap2 run`` lays it out
#: (``BIN_REFINEMENT/sample_7/metawrap_70_10_bins.stats.tsv`` -> ``sample_7``). Omitted for a
#: single bin set.
SAMPLE_COLUMN = "sample"

#: Columns before the per-sample ones: name, size, quality, then taxonomy.
#:
#: The ``mimag_tag`` name comes from :data:`metawrap2.mimag.STATS_COLUMNS`, so it matches the
#: column ``mimag.annotate_stats`` appends to a ``.stats.tsv``.
BASE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "size",
    "contigs",
    "n50",
    "gc",
    "completeness",
    "contamination",
    "quality_tier",
    _mimag.STATS_COLUMNS[0],
    "gunc_pass",
    "gunc_css",
    "taxonomy",
)

#: Tier order for a distribution plot: best first, with :data:`BLANK` last under its own label,
#: so unscored genomes are shown.
TIER_ORDER: Tuple[str, ...] = tuple(name[len("bins_") :] for name, _c, _x in QUALITY_TIERS) + (
    BELOW_TIERS,
    BLANK,
)

#: Label for the :data:`BLANK` tier in a plot or a count.
NOT_MEASURED = "not measured"

#: Current suffixes mapped to the ones earlier runs wrote. Both are accepted so older output can be
#: read; the repo's reference run has ``metawrap_50_10_bins.stats`` and ``bin_abundance_table.tab``.
LEGACY_SUFFIX = {STATS_SUFFIX: ".stats", CONTIGS_SUFFIX: ".contigs", TSV_SUFFIX: ".tab"}

#: ``classify_bins``' per-bin consensus taxonomy, and ``quant_bins``' abundance matrix.
TAXONOMY_NAME = "bin_taxonomy" + TSV_SUFFIX
ABUNDANCE_NAME = "bin_abundance_table" + TSV_SUFFIX

#: Directory levels to descend below a run directory. Files in a directory at this depth are still
#: read, so 2 reaches ``study/BIN_REFINEMENT/sample_7/metawrap_70_10_bins.stats.tsv`` from the
#: study root. Bounded because an assembly directory holds a metaSPAdes tree thousands of entries
#: deep.
MAX_DEPTH = 2

#: ``bin_refinement``'s scratch directory, skipped when looking for a bin set.
#:
#: It holds stats tables for intermediate bin sets (``binsA``, ``binsAB``, ``binsO.to_score``;
#: seven per sample in the reference study). It also holds the GUNC report, so it is only skipped
#: when looking for stats tables.
WORK_DIR = "work_files"

#: First-column values treated as a header in two-column files that are specified without one but
#: sometimes have one (``classify_bins --gtdbtk`` writes ``user_genome``).
HEADER_WORDS = frozenset({"bin", "bins", "genome", "genomic bins", "user_genome", "name"})

__all__ = [
    "ABUNDANCE_PREFIX",
    "BASE_COLUMNS",
    "BELOW_TIERS",
    "BLANK",
    "INPUT_KINDS",
    "NOT_MEASURED",
    "SAMPLE_COLUMN",
    "TABLE_NAME",
    "TIER_ORDER",
    "BinSet",
    "Inputs",
    "Table",
    "build",
    "discover",
    "find_gunc",
    "find_stats",
    "quality_tier",
    "read_abundance",
    "read_contig_counts",
    "read_stats",
    "read_taxonomy",
    "read_tsv",
    "sample_of",
    "tier_counts",
    "write",
]


class BinSet(NamedTuple):
    """One set of genomes and the three per-bin-set files that describe it.

    Grouped because the three files must describe the same bins. Every sample numbers its bins from
    ``bin_001``, so files from different samples would join without error; in the reference study
    that gave 23 wrong chimerism verdicts and one genome present in neither bin set.
    """

    #: Which sample's bins these are, for the :data:`SAMPLE_COLUMN` cell.
    sample: str
    stats: Optional[str] = None
    contigs: Optional[str] = None
    gunc: Optional[str] = None


class Inputs(NamedTuple):
    """Inputs for a table. Every part is optional and there may be many bin sets.

    *taxonomy* and *abundance* are study-wide: ``classify_bins`` and ``quant_bins`` each write one
    table per study. In a cohort run they are keyed on the dereplicated catalogue's names, which
    need not match any sample's bin names.
    """

    bin_sets: Tuple[BinSet, ...] = ()
    taxonomy: Optional[str] = None
    abundance: Optional[str] = None


#: Input kinds, as named in :attr:`Table.absent` and :attr:`Table.unjoined`. Listed explicitly
#: because three of them are fields of :class:`BinSet`, not of ``Inputs``.
INPUT_KINDS: Tuple[str, ...] = ("stats", "contigs", "gunc", "taxonomy", "abundance")


class Table(NamedTuple):
    """One row per genome, plus the source files, missing inputs and inputs that did not join."""

    columns: List[str]
    rows: List[List[str]]
    #: The abundance matrix's sample names - the ``abundance.*`` columns, not
    #: :data:`SAMPLE_COLUMN`, which names the bin set a row came from.
    abundance_samples: List[str]
    #: input kind -> the file it was read from.
    sources: Dict[str, str]
    #: Input kinds with no file anywhere; their columns are blank in every row.
    absent: List[str]
    #: ``(kind, sample)`` for a file that exists but matched none of that bin set's genomes. The
    #: cells are left blank and the caller reports it.
    unjoined: List[Tuple[str, str]]
    #: How many bin sets the rows came from. More than one means :data:`SAMPLE_COLUMN` is present.
    bin_sets: int = 0


# ── finding the inputs ───────────────────────────────────────────────────────────────────────


def _files(dirs: Sequence[str]) -> Iterator[str]:
    """Every file at most :data:`MAX_DEPTH` below each of *dirs*, in a stable order."""
    for root in dirs:
        if not os.path.isdir(root):
            continue
        base_depth = os.path.abspath(root).count(os.sep)
        for parent, subdirs, names in os.walk(os.path.abspath(root)):
            if parent.count(os.sep) - base_depth >= MAX_DEPTH:
                subdirs[:] = []
            subdirs.sort()
            for name in sorted(names):
                yield os.path.join(parent, name)


def _endings(suffix: str) -> Tuple[str, ...]:
    """A suffix and the one earlier runs wrote, for matching either."""
    legacy = LEGACY_SUFFIX.get(suffix)
    return (suffix,) if legacy is None else (suffix, legacy)


def _named(dirs: Sequence[str], name: str) -> Optional[str]:
    """The first file called *name* (or its legacy spelling) under *dirs*."""
    stem, suffix = os.path.splitext(name)
    wanted = {stem + ending for ending in _endings(suffix)}
    for path in _files(dirs):
        if os.path.basename(path) in wanted:
            return path
    return None


def find_stats(dirs: Sequence[str]) -> List[str]:
    """Every bin set's CheckM table under *dirs*: one per directory, preferring a consolidated set.

    One table per directory, from every directory:

    * A ``metawrap2 run`` study has one bin set per sample (``BIN_REFINEMENT/sample_6``,
      ``sample_7``, ``sample_13``); all are returned (in the reference study, 62 genomes, against 24
      from the first sample alone).
    * Tables in one directory are alternative bin sets of one assembly (``bin_refinement`` leaves
      one per input binner beside the consolidated set), so only one is taken to avoid counting
      genomes several times. The consolidated set is preferred; two consolidated sets (a re-run at
      different thresholds, ``metawrap_50_10_bins`` and ``metawrap_70_10_bins``) still give one.

    A directory with no consolidated set contributes its alphabetically first table (the case for
    a ``binning``-only run).
    """
    by_directory: Dict[str, List[str]] = {}
    for path in _files(dirs):
        if not any(path.endswith(ending) for ending in _endings(STATS_SUFFIX)):
            continue
        if WORK_DIR in os.path.dirname(path).split(os.sep):
            continue
        by_directory.setdefault(os.path.dirname(path), []).append(path)

    refined = {
        directory: [p for p in paths if is_refined_bins_dirname(_stats_stem(p))]
        for directory, paths in by_directory.items()
    }
    # If any consolidated set exists, use only consolidated sets. Mixing in a raw binner's table
    # from another directory would list bins alongside the bins they were consolidated from.
    if any(refined.values()):
        return [min(paths) for _directory, paths in sorted(refined.items()) if paths]
    return [min(paths) for _directory, paths in sorted(by_directory.items())]


def _stats_stem(stats_path: str) -> str:
    """``.../metawrap_50_10_bins.stats.tsv`` -> ``metawrap_50_10_bins``."""
    base = os.path.basename(stats_path)
    for ending in _endings(STATS_SUFFIX):
        if base.endswith(ending):
            return base[: -len(ending)]
    return base


def find_contigs(stats_path: str) -> Optional[str]:
    """The contig-membership table for the same bin set as *stats_path*.

    Derived from the stats path so the contig counts are for the same bins as the rest of the row.
    """
    root = os.path.join(os.path.dirname(stats_path), _stats_stem(stats_path))
    for ending in _endings(CONTIGS_SUFFIX):
        if os.path.isfile(root + ending):
            return root + ending
    return None


def sample_of(stats_path: str) -> str:
    """The sample a bin set belongs to: the name of its stats table's directory."""
    return os.path.basename(os.path.dirname(os.path.abspath(stats_path)))


def discover(run_dirs: Sequence[str]) -> Inputs:
    """Find every bin set under *run_dirs*, plus the study-wide taxonomy and abundance tables."""
    bin_sets = tuple(
        BinSet(
            sample=sample_of(stats),
            stats=stats,
            contigs=find_contigs(stats),
            gunc=find_gunc(stats),
        )
        for stats in find_stats(run_dirs)
    )
    if not bin_sets:
        # No stats table, but GUNC may have run. With no bin set to mismatch, a broad search for
        # its report is safe here. Rows then come from any taxonomy or abundance table.
        found = _find_gunc_under(run_dirs)
        if found:
            bin_sets = (BinSet(sample=sample_of(found), gunc=found),)
    return Inputs(
        bin_sets=bin_sets,
        taxonomy=_named(run_dirs, TAXONOMY_NAME),
        abundance=_named(run_dirs, ABUNDANCE_NAME),
    )


def find_gunc(stats_path: str) -> Optional[str]:
    """GUNC's report for the same bins as *stats_path*, or None.

    Searched only in the stats table's directory (``bin_refinement`` leaves the report in that
    sample's ``work_files/gunc/``). A study-wide search would return one report for all samples,
    and since every sample numbers bins from ``bin_001`` the wrong verdicts would join.

    Uses :func:`metawrap2.gunc.find_report`, since the report filename contains the reference
    database version.
    """
    return _gunc.find_report(os.path.dirname(os.path.abspath(stats_path)))


def _find_gunc_under(dirs: Sequence[str]) -> Optional[str]:
    """Any GUNC report under *dirs*, for a run with no stats table to pair one with."""
    for root in dirs:
        if not os.path.isdir(root):
            continue
        # abspath so all paths printed in a report have the same form.
        found = _gunc.find_report(os.path.abspath(root))
        if found:
            return found
    return None


# ── reading each input ───────────────────────────────────────────────────────────────────────


def read_tsv(path: str) -> Iterator[List[str]]:
    """Non-empty tab-separated rows of *path*.

    Public so :mod:`metawrap2.report` reads the same tables the same way.
    """
    with open(path, newline="") as fh:
        for fields in csv.reader(fh, delimiter="\t"):
            if fields and any(field.strip() for field in fields):
                yield fields


def read_stats(path: str) -> Dict[str, Dict[str, str]]:
    """A ``.stats.tsv`` as ``genome -> column -> value``, keyed on :func:`bin_label`.

    Column names are lowercased (``GC`` and ``gc`` are one column). Every column is kept, including
    those ``gunc.annotate_stats`` and ``mimag.annotate_stats`` append, so this table reports the
    values already in the file.

    File order is preserved (``checkm.write_stats`` sorts by completeness, descending).
    """
    out: Dict[str, Dict[str, str]] = {}
    rows = read_tsv(path)
    try:
        header = [cell.strip().lower() for cell in next(rows)]
    except StopIteration:
        return out
    for fields in rows:
        name = bin_label(fields[0].strip())
        if not name or is_unbinned_name(fields[0].strip()):
            # Skip the unbinned catch-all (e.g. `binning` writes concoct_bins/unbinned.fasta); it
            # is not a genome and would be counted in per-genome totals.
            continue
        out[name] = {
            key: fields[i].strip() if i < len(fields) else BLANK
            for i, key in enumerate(header)
            if key
        }
    return out


def read_contig_counts(path: str) -> Dict[str, int]:
    """Contig count per genome, from the ``contig -> bin`` membership table."""
    counts: Dict[str, int] = {}
    for fields in read_tsv(path):
        if len(fields) < 2:
            continue
        name = bin_label(fields[1].strip())
        if not name or name.lower() in HEADER_WORDS:
            continue
        counts[name] = counts.get(name, 0) + 1
    return counts


def read_taxonomy(path: str) -> Dict[str, str]:
    """``classify_bins``' per-bin consensus lineage, keyed on :func:`bin_label`.

    An empty lineage is kept as :data:`BLANK`: ``classify_bins`` writes a row with no lineage for a
    bin it could not place, which differs from a bin that was not classified.
    """
    out: Dict[str, str] = {}
    for fields in read_tsv(path):
        raw = fields[0].strip()
        if raw.lower() in HEADER_WORDS:
            continue
        name = bin_label(raw)
        if name:
            out[name] = fields[1].strip() if len(fields) > 1 else BLANK
    return out


def read_abundance(path: str) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    """``quant_bins``' matrix as ``(sample names, genome -> sample -> value)``.

    Values are kept as text so they match the matrix exactly (the four quantifiers write different
    precisions).
    """
    rows = read_tsv(path)
    try:
        header = [cell.strip() for cell in next(rows)]
    except StopIteration:
        return [], {}
    samples = [cell for cell in header[1:] if cell]
    values: Dict[str, Dict[str, str]] = {}
    for fields in rows:
        name = bin_label(fields[0].strip())
        if not name or name.lower() in HEADER_WORDS:
            continue
        values[name] = {
            sample: fields[i + 1].strip() if i + 1 < len(fields) else BLANK
            for i, sample in enumerate(samples)
        }
    return samples, values


# ── the table ────────────────────────────────────────────────────────────────────────────────


def _number(raw: object) -> Optional[float]:
    """*raw* as a float, or None if it does not state one. ``""`` and ``NA`` are None, not 0."""
    if raw is None:
        return None
    text = str(raw).strip()
    if text.lower() in ("", "na", "n/a", "nan", "none", "-"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def quality_tier(completeness: Optional[float], contamination: Optional[float]) -> str:
    """The strongest :data:`metawrap2.benchmark.scoring.QUALITY_TIERS` tier a genome clears.

    Uses the same thresholds as the counts in ``docs/benchmarks.md``. ``bins_high_quality``
    (70/10) is this project's threshold, not a MIMAG tier; MIMAG is the separate ``mimag_tag``
    column, which also uses rRNA and tRNA evidence.

    One tier per genome, strongest first. ``scoring.quality_counts`` counts non-exclusively, which
    suits benchmarks but would put one genome in several bars of a chart.

    Either measurement missing gives :data:`BLANK`.
    """
    if completeness is None or contamination is None:
        return BLANK
    for name, min_completeness, max_contamination in QUALITY_TIERS:
        if completeness >= min_completeness and contamination <= max_contamination:
            return name[len("bins_") :]
    return BELOW_TIERS


def build(inputs: Inputs) -> Table:
    """Join whatever of *inputs* exists into one row per genome, across every bin set.

    Rows are built per bin set and appended, not keyed by genome name across bin sets, because
    every sample numbers its bins from ``bin_001`` (in the reference study, keying by name would
    collapse 62 genomes into 24).
    """
    taxonomy = read_taxonomy(inputs.taxonomy) if inputs.taxonomy else {}
    abundance_samples, abundance = (
        read_abundance(inputs.abundance) if inputs.abundance else ([], {})
    )

    many = len(inputs.bin_sets) > 1
    columns = (
        ([SAMPLE_COLUMN] if many else [])
        + list(BASE_COLUMNS)
        + [ABUNDANCE_PREFIX + sample for sample in abundance_samples]
    )
    rows: List[List[str]] = []
    unjoined: List[Tuple[str, str]] = []
    sources: Dict[str, str] = {}
    placed: set = set()

    for bin_set in inputs.bin_sets:
        stats = read_stats(bin_set.stats) if bin_set.stats else {}
        contigs = read_contig_counts(bin_set.contigs) if bin_set.contigs else {}
        # GUNC keys its report on the bin filenames it was given, which may differ from the stats
        # table's names (gzipped or other suffixes), so re-key through bin_label.
        gunc = {
            bin_label(name): verdict
            for name, verdict in (_gunc.parse_report(bin_set.gunc) if bin_set.gunc else {}).items()
        }
        for kind, path in (
            ("stats", bin_set.stats),
            ("contigs", bin_set.contigs),
            ("gunc", bin_set.gunc),
        ):
            if path:
                sources.setdefault(kind, path)

        names = _genomes(stats, contigs, gunc)
        unjoined.extend(
            (kind, bin_set.sample)
            for kind, table in (
                ("gunc", gunc),
                ("taxonomy", taxonomy),
                ("abundance", abundance),
            )
            if table and names and not any(name in table for name in names)
        )
        prefix = [bin_set.sample] if many else []
        for name in names:
            rows.append(
                prefix
                + _row(
                    name,
                    stats.get(name, {}),
                    contigs,
                    gunc,
                    taxonomy,
                    abundance_samples,
                    abundance,
                )
            )
        placed.update(names)

    # Genomes only in a study-wide table (e.g. every row of the abundance matrix after cohort
    # dereplication) are kept, with a blank sample cell since their bin set is unknown.
    for name in sorted((set(taxonomy) | set(abundance)) - placed):
        rows.append(
            ([BLANK] if many else [])
            + _row(name, {}, {}, {}, taxonomy, abundance_samples, abundance)
        )

    for kind, path in (("taxonomy", inputs.taxonomy), ("abundance", inputs.abundance)):
        if path:
            sources[kind] = path
    return Table(
        columns=columns,
        rows=rows,
        abundance_samples=abundance_samples,
        sources=sources,
        absent=[kind for kind in INPUT_KINDS if kind not in sources],
        unjoined=unjoined,
        bin_sets=len(inputs.bin_sets),
    )


def _genomes(*tables: Mapping[str, object]) -> List[str]:
    """Every genome named by any of one bin set's files, in the stats table's order then the rest.

    The stats table's order comes first. Genomes only in a later file (e.g. the contig membership
    table) are appended, sorted.
    """
    first = tables[0] if tables else {}
    order = list(first)
    extra = {name for table in tables[1:] for name in table} - set(order)
    return order + sorted(extra)


def _row(
    name: str,
    stat: Dict[str, str],
    contigs: Dict[str, int],
    gunc: Dict[str, _gunc.Verdict],
    taxonomy: Dict[str, str],
    samples: Sequence[str],
    abundance: Dict[str, Dict[str, str]],
) -> List[str]:
    """One genome's row. Every cell is either copied from an input or :data:`BLANK`."""
    completeness = _number(stat.get("completeness"))
    contamination = _number(stat.get("contamination"))
    verdict = gunc.get(name)
    # Tags already in the stats table (appended by `gunc` and `mimag`) take precedence over values
    # derived here, so the two files agree.
    gunc_pass = stat.get("gunc_pass") or (
        BLANK if verdict is None else ("yes" if verdict.passed else "no")
    )
    gunc_css = stat.get("gunc_css") or (BLANK if verdict is None else "%.4f" % verdict.css)
    tag = stat.get(_mimag.STATS_COLUMNS[0]) or _mimag_tag(name, completeness, contamination)
    counted = contigs.get(name)
    row = [
        name,
        stat.get("size", BLANK),
        BLANK if counted is None else str(counted),
        stat.get("n50", BLANK),
        stat.get("gc", BLANK),
        stat.get("completeness", BLANK),
        stat.get("contamination", BLANK),
        quality_tier(completeness, contamination),
        tag,
        gunc_pass,
        gunc_css,
        taxonomy.get(name, BLANK),
    ]
    per_sample = abundance.get(name, {})
    return row + [per_sample.get(sample, BLANK) for sample in samples]


def _mimag_tag(name: str, completeness: Optional[float], contamination: Optional[float]) -> str:
    """:mod:`metawrap2.mimag`'s tag for a genome whose stats table was never tagged.

    Called without rRNA/tRNA evidence, since this module does not search ``annotate`` output for a
    16S. A genome that clears 90/5 is therefore tagged high quality unverified.
    """
    if completeness is None and contamination is None:
        # Blank instead of mimag's "unscored" tag, so unmeasured cells look the same everywhere.
        return BLANK
    return _mimag.assess(name, completeness, contamination).tag.value


def tier_counts(table: Table) -> List[Tuple[str, int]]:
    """How many genomes are in each tier, in :data:`TIER_ORDER`, for a distribution plot.

    Empty tiers are kept so charts of different runs have the same bars. The :data:`BLANK` tier is
    labelled :data:`NOT_MEASURED`.
    """
    at = table.columns.index("quality_tier")
    counts = {tier: 0 for tier in TIER_ORDER}
    for row in table.rows:
        tier = row[at] if row[at] in counts else BLANK
        counts[tier] += 1
    return [(tier or NOT_MEASURED, counts[tier]) for tier in TIER_ORDER]


def write(table: Table, path: str) -> str:
    """Write *table* as a TSV and return the path."""
    directory = os.path.dirname(os.path.abspath(path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(table.columns) + "\n")
        fh.writelines("\t".join(row) + "\n" for row in table.rows)
    return path
