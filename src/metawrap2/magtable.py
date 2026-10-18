"""One row per genome, joined from the tables a finished run already wrote.

A finished run leaves twenty-odd correct tables and not one row per genome. Completeness and
contamination are in ``<bins>.stats.tsv``, the chimerism verdict is in GUNC's own report, the
lineage is in ``classify_bins``' two-column file, and the per-sample abundance is in
``quant_bins``' matrix - so the simplest question anyone asks of a MAG ("how good is it, what is
it, and how much of it is there") needs four files opened side by side and joined by hand. This
module does that join once.

Nothing here measures anything. Every cell is copied from a file the pipeline already wrote, which
is why the GC column is left as the fraction CheckM writes rather than rescaled to a percentage:
a value copied verbatim cannot disagree with the file beside it, and a value converted twice is a
silent factor of 100.

Two things decide whether this table is trustworthy:

* **An absent measurement is blank, never zero.** A contamination of 0 is a finding; an unmeasured
  contamination is not. Three other tools' tables in this project have conflated the two, and the
  resulting numbers were believed, because 0 looks like an answer. :data:`BLANK` is the only thing
  this module ever writes for a measurement that was not made.
* **Every input is optional and no input can drop a genome.** A run that skipped ``classify_bins``
  gets the whole table with an empty taxonomy column, not a shorter table and not a failure. A
  genome that appears only in the abundance matrix still gets a row, because the alternative is a
  table that silently answers a different question than the one it was asked.

Genome names are joined through :func:`metawrap2.constants.bin_label`, which is the one place in
this project that knows ``bin_001.fasta``, ``bin.1.fa`` and ``bin_001`` are the same genome. The
four files being joined here use three of those conventions between them - a real run has
``bin.1`` in the stats table, ``bin.1.fa`` in the taxonomy file and ``bin_001`` in GUNC's report -
so a fourth copy of the stripping rule is the obvious way for this join to lose rows.
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

#: What an unmeasured cell holds. Not "0", not "0.0", not "NA" - an empty cell, which every reader
#: of a TSV already treats as "no value" and no reader averages, sums or plots.
BLANK = ""

#: The table this module writes, beside whatever else a caller is producing.
TABLE_NAME = "mag_table" + TSV_SUFFIX

#: The one label for a genome scored but below every tier, so "measured and poor" is distinct from
#: :data:`BLANK` ("not measured"). Those two were the same cell in the tables this replaces.
BELOW_TIERS = "below_medium_quality"

#: ``quant_bins``' abundance columns are prefixed rather than used bare. The unit depends on which
#: quantifier ran (salmon and the depth route give coverage, sylph a relative abundance), so the
#: prefix says what kind of number it is; and a sample legitimately named "genome" or "n50" would
#: otherwise collide with a column of this table and overwrite it.
ABUNDANCE_PREFIX = "abundance."

#: The leftmost column when a report covers more than one bin set: which sample's bin set a row
#: came from. Taken from the directory the stats table sits in, because that is where ``metawrap2
#: run`` puts it (``BIN_REFINEMENT/sample_7/metawrap_70_10_bins.stats.tsv`` -> ``sample_7``), and
#: omitted entirely for a single bin set - a per-run report should not grow an empty field.
SAMPLE_COLUMN = "sample"

#: The columns before the per-sample ones, in the order a reader wants them: what the genome is
#: called, how big it is, how good it is, and what it is.
#:
#: ``mimag_tag`` is named from :data:`metawrap2.mimag.STATS_COLUMNS` rather than spelled out, so the
#: column here and the column ``mimag.annotate_stats`` appends to a ``.stats.tsv`` cannot drift
#: apart and be read as two different things.
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

#: Tier order for a distribution plot: best first, with :data:`BLANK` last under its own label.
#: An unscored genome is shown rather than dropped - a bar chart that omits them makes a run where
#: CheckM failed look like a run with fewer genomes.
TIER_ORDER: Tuple[str, ...] = tuple(name[len("bins_") :] for name, _c, _x in QUALITY_TIERS) + (
    BELOW_TIERS,
    BLANK,
)

#: Label for the :data:`BLANK` tier in a plot or a count, where an empty string reads as a bug.
NOT_MEASURED = "not measured"

#: Suffixes the current code writes, each with the one earlier runs wrote. Both are accepted
#: because a report is most wanted over output that already exists: this repo's own reference run
#: has ``metawrap_50_10_bins.stats`` and ``bin_abundance_table.tab``, written before the suffixes
#: were normalised, and a reader that only knew the new names would call that whole run absent.
LEGACY_SUFFIX = {STATS_SUFFIX: ".stats", CONTIGS_SUFFIX: ".contigs", TSV_SUFFIX: ".tab"}

#: ``classify_bins``' per-bin consensus taxonomy, and ``quant_bins``' abundance matrix.
TAXONOMY_NAME = "bin_taxonomy" + TSV_SUFFIX
ABUNDANCE_NAME = "bin_abundance_table" + TSV_SUFFIX

#: How many directory levels below a run directory to descend. Files *in* a directory at this depth
#: are still read, so 2 reaches a per-sample study's stats table from the study root
#: (``study/BIN_REFINEMENT/sample_7/metawrap_70_10_bins.stats.tsv``). Bounded because an assembly
#: directory holds a metaSPAdes tree thousands of entries deep, and unbounded recursion here turned
#: a report over a real run directory into a walk of the whole assembly.
MAX_DEPTH = 2

#: ``bin_refinement``'s scratch directory, skipped when looking for a bin set.
#:
#: It holds a stats table per *intermediate* bin set - ``binsA``, ``binsAB``, ``binsO.to_score`` -
#: which are the arithmetic of consolidation rather than anybody's answer: the reference study has
#: seven of them per sample against one final table. It also holds the GUNC report, which *is*
#: wanted, so this is skipped for stats tables specifically and not walked around entirely.
WORK_DIR = "work_files"

#: First-column values that mean "this is a header", for the two-column files that have no header
#: by contract but have gained one twice (``classify_bins --gtdbtk`` writes ``user_genome``).
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
    """One set of genomes, and the three per-bin-set files that describe it.

    Grouped rather than looked up independently because all three have to be *the same bins*. Every
    sample in a study names its bins ``bin_001`` onwards, so a stats table from one sample and a
    GUNC report from another join perfectly and silently: on the reference study that produced 23
    wrong chimerism verdicts plus one phantom genome that existed in neither sample's bin set.
    """

    #: Which sample's bins these are, for the :data:`SAMPLE_COLUMN` cell.
    sample: str
    stats: Optional[str] = None
    contigs: Optional[str] = None
    gunc: Optional[str] = None


class Inputs(NamedTuple):
    """What a table is built from. Every part is optional and there may be many bin sets.

    *taxonomy* and *abundance* are study-wide rather than per-bin-set: ``classify_bins`` and
    ``quant_bins`` each write one table for the whole study, and in a cohort run they are keyed on
    the *dereplicated* catalogue's names, which need not match any sample's bin names.
    """

    bin_sets: Tuple[BinSet, ...] = ()
    taxonomy: Optional[str] = None
    abundance: Optional[str] = None


#: The kinds of input a cell can come from, which are the names used in :attr:`Table.absent` and
#: :attr:`Table.unjoined`. Spelled out rather than taken from ``Inputs._fields``, because three of
#: them live inside a :class:`BinSet` and a reader of a report wants to hear "no taxonomy", not
#: "no bin_sets".
INPUT_KINDS: Tuple[str, ...] = ("stats", "contigs", "gunc", "taxonomy", "abundance")


class Table(NamedTuple):
    """One row per genome, plus what it was built from, what was missing, and what would not join."""

    columns: List[str]
    rows: List[List[str]]
    #: The abundance matrix's sample names - the ``abundance.*`` columns, not
    #: :data:`SAMPLE_COLUMN`, which names the bin set a row came from.
    abundance_samples: List[str]
    #: input kind -> the file(s) it was read from, so a reader can check a surprising cell.
    sources: Dict[str, str]
    #: input kinds with no file anywhere. These are the columns that are blank for every row, and a
    #: reader has to be told the difference between that and a measurement that came back empty.
    absent: List[str]
    #: ``(kind, sample)`` for a file that exists but matched none of that bin set's genomes - a
    #: join that could not be made rather than one that came back empty. The cells are left blank
    #: and the caller says so, because the alternative is guessing which genome is which.
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

    **One per directory, and all of the directories.** Two different mistakes are possible here and
    the rule has to avoid both:

    * *Reporting a subset as the whole.* A ``metawrap2 run`` study puts one bin set per sample
      (``BIN_REFINEMENT/sample_6``, ``sample_7``, ``sample_13``), and returning only the first of
      them described 24 genomes in a study that recovered 62, with nothing on the page saying which
      sample the 24 came from.
    * *Counting the same genomes twice.* Within one directory the several tables are competing
      answers over one assembly - ``bin_refinement`` leaves one per input binner beside the
      consolidated set, so reporting all four would report every genome up to four times and
      quadruple the study's apparent yield.

    So: group by directory, and within a directory take the consolidated set if there is one. Two
    consolidated sets in one directory (a re-run at different thresholds, ``metawrap_50_10_bins``
    beside ``metawrap_70_10_bins``) are also the same bins cut twice, so that is still one table.

    A directory with no consolidated set contributes its alphabetically first table, which is what
    a ``binning``-only run looks like.
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
    # Consolidated sets win outright where any exist. Mixing a consolidated set from one directory
    # with a raw binner's table from another would put bins and the bins they were consolidated
    # from in one table as if they were different genomes.
    if any(refined.values()):
        return [min(paths) for _directory, paths in sorted(refined.items()) if paths]
    return [min(paths) for _directory, paths in sorted(by_directory.items())]


def _stats_stem(stats_path: str) -> str:
    """``.../metawrap_50_10_bins.stats.tsv`` -> ``metawrap_50_10_bins``: the bin set it describes."""
    base = os.path.basename(stats_path)
    for ending in _endings(STATS_SUFFIX):
        if base.endswith(ending):
            return base[: -len(ending)]
    return base


def find_contigs(stats_path: str) -> Optional[str]:
    """The contig-membership table for the same bin set as *stats_path*.

    Looked up from the stats path rather than searched for, because the contig count has to be of
    the bins the rest of the row describes. A refined row with CONCOCT's contig count would be
    wrong in a way nothing on the page could reveal.
    """
    root = os.path.join(os.path.dirname(stats_path), _stats_stem(stats_path))
    for ending in _endings(CONTIGS_SUFFIX):
        if os.path.isfile(root + ending):
            return root + ending
    return None


def sample_of(stats_path: str) -> str:
    """Which sample a bin set belongs to: the name of the directory its stats table sits in."""
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
        # No stats table anywhere, but GUNC may still have run - and with no bin set to pair it
        # with there is only one thing it can be about, so the broad search is safe here and only
        # here. The rows then come from whatever taxonomy or abundance table exists.
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

    Searched inside the stats table's own directory, not across the run directories. GUNC's output
    is a sibling of the bin set it scored (``bin_refinement`` leaves it in that sample's
    ``work_files/gunc/``), and a study-wide search returned *one* report for *all* samples: on the
    reference study, sample_6's verdicts were joined onto sample_13's genomes, which matched
    because both sets are numbered from ``bin_001``.

    Through :func:`metawrap2.gunc.find_report` rather than a second pattern of our own, because the
    middle of that filename carries the reference database version and a literal name here would
    find nothing the first time the database is updated.
    """
    return _gunc.find_report(os.path.dirname(os.path.abspath(stats_path)))


def _find_gunc_under(dirs: Sequence[str]) -> Optional[str]:
    """Any GUNC report under *dirs*, for a run with no stats table to pair one with."""
    for root in dirs:
        if not os.path.isdir(root):
            continue
        # abspath because these paths are printed next to each other on a report page, and a
        # mixture of absolute and relative ones reads as two different runs.
        found = _gunc.find_report(os.path.abspath(root))
        if found:
            return found
    return None


# ── reading each input ───────────────────────────────────────────────────────────────────────


def read_tsv(path: str) -> Iterator[List[str]]:
    """Non-empty tab-separated rows of *path*.

    Public because :mod:`metawrap2.report` reads the same tables for the same reason, and two
    copies of "skip the blank lines" is how two readers of one file come to disagree about how many
    rows it has.
    """
    with open(path, newline="") as fh:
        for fields in csv.reader(fh, delimiter="\t"):
            if fields and any(field.strip() for field in fields):
                yield fields


def read_stats(path: str) -> Dict[str, Dict[str, str]]:
    """A ``.stats.tsv`` as ``genome -> column -> value``, keyed on :func:`bin_label`.

    Column names are lowercased so ``GC`` and ``gc`` are one column, and every column is kept -
    including the ones ``gunc.annotate_stats`` and ``mimag.annotate_stats`` append - because a
    value already written into this file is the value this table must report. Re-deriving a tag
    that is sitting in the file is how two tables of the same run come to disagree.

    File order is preserved: ``checkm.write_stats`` sorts by completeness descending, which is the
    order a reader of that file already knows, and re-sorting here would make the two tables of the
    same bin set look like tables of different ones.
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
            # The unbinned catch-all is not a genome. It is a legitimate member of a bin set
            # (`binning` writes concoct_bins/unbinned.fasta on purpose) and a row for it in a MAG
            # table would be counted in every per-genome total on the page.
            continue
        out[name] = {
            key: fields[i].strip() if i < len(fields) else BLANK
            for i, key in enumerate(header)
            if key
        }
    return out


def read_contig_counts(path: str) -> Dict[str, int]:
    """How many contigs each genome holds, from the ``contig -> bin`` membership table."""
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

    An empty lineage in the file is kept as :data:`BLANK` rather than skipped: ``classify_bins``
    writes a row with no lineage for a bin it could not place, and that is a different statement
    from the bin not having been classified at all.
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

    Values are carried as text, not parsed to float and reformatted. The matrix is written by one
    of four quantifiers with different precision, and rounding it here would make the report and
    the table it was built from disagree in the last digit for no gain.
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

    The thresholds are the project's own, reused rather than restated, so this column means the
    same thing as the counts in ``docs/benchmarks.md``. Note that
    ``bins_high_quality`` (70/10) is this project's working threshold and not a MIMAG tier - the
    MIMAG answer is the separate ``mimag_tag`` column, which also accounts for the rRNA and tRNA
    evidence that completeness and contamination cannot speak to.

    One tier per genome, strongest first, because a table needs one label per row; the
    non-exclusive counting ``scoring.quality_counts`` does is right for a benchmark and would put
    the same genome in three bars of one chart here.

    Either measurement missing gives :data:`BLANK`. A genome with no completeness has not scored 0.
    """
    if completeness is None or contamination is None:
        return BLANK
    for name, min_completeness, max_contamination in QUALITY_TIERS:
        if completeness >= min_completeness and contamination <= max_contamination:
            return name[len("bins_") :]
    return BELOW_TIERS


def build(inputs: Inputs) -> Table:
    """Join whatever of *inputs* exists into one row per genome, across every bin set.

    Rows are built per bin set and appended, never keyed by genome name across bin sets. Every
    sample in a study numbers its bins from ``bin_001``, so a table keyed on the name alone would
    have each sample's ``bin_001`` overwrite the last - 62 genomes collapsing into the 24 distinct
    names of the largest bin set, with no duplicate to notice.
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
        # GUNC keys its report on the bin filenames it was given, which for a gzipped or
        # differently suffixed bin set is not the spelling the stats table uses. Re-keyed once per
        # bin set rather than at each lookup, so the join is on one convention throughout.
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

    # A genome named only by a study-wide table - in the abundance matrix but in no bin set, which
    # is every row of it once a cohort has been dereplicated. Kept, because dropping it would make
    # this table answer a narrower question than the files it came from; and with a blank sample
    # cell, because which bin set it came from is exactly what is not known.
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

    The stats table comes first and keeps its order, because that is the ordering a reader of the
    ``.stats.tsv`` already has in their head. A genome named only by a later file - in the contig
    membership table but not in the stats table - is appended rather than dropped, because dropping
    it would make this table quietly answer a narrower question than the files it came from.
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
    # The tags already written into the stats table win over anything derived here: `gunc` and
    # `mimag` both append their columns to that file, and a report that re-derived them could
    # disagree with the file a user opens next to it.
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

    Called with no rRNA/tRNA evidence, which is the honest state here: this module joins tables and
    does not go looking for an ``annotate`` output to search for a 16S. A genome that clears 90/5
    therefore comes back as high quality *unverified*, which is what MIMAG allows to be said when
    nothing has looked.
    """
    if completeness is None and contamination is None:
        # Left blank rather than passed to mimag, which would correctly return its "unscored" tag -
        # a word in a cell that reads as a verdict. A cell nothing measured should look unmeasured.
        return BLANK
    return _mimag.assess(name, completeness, contamination).tag.value


def tier_counts(table: Table) -> List[Tuple[str, int]]:
    """How many genomes are in each tier, in :data:`TIER_ORDER`, for a distribution plot.

    Empty tiers are kept so two runs' charts have the same bars in the same places, and the
    :data:`BLANK` tier is labelled :data:`NOT_MEASURED` rather than left as an empty axis label.
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
