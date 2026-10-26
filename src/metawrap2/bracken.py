"""Bracken: abundance estimates at one rank from Kraken2 read assignments.

A Kraken2 report's third column counts reads assigned at a node, and Kraken2 assigns a read to the
lowest node all its k-mers agree on. Reads from sequence conserved across a genus stop at the genus
node, so a species at 10% of the community can have only a small share of reads at its own node.
Reading the species rows as composition understates every species and adds a large "genus-level"
share that is not an organism.

Bracken (Lu et al., *PeerJ Comput Sci* 2017, 10.7717/peerj-cs.104) redistributes those reads down
to a chosen rank. For each node above the rank it uses the k-mer content of the reference genomes
to estimate how likely a read of a given length from each candidate genome is to stop at that node,
and splits the node's reads accordingly. Two constraints follow:

* The k-mer distribution belongs to one Kraken2 database: it is built by classifying that
  database's reference genomes against it. ``bracken-build`` writes its files into the Kraken2
  database directory, so ``--bracken-db`` defaults there.
* It also belongs to one read length, because the chance of a read stopping at a node depends on
  how many k-mers it has. Each length is stored as ``database<N>mers.kmer_distrib``, and Bracken
  cannot estimate for a length that was not built.

On the fixture under ``tests/fixtures/bracken/`` (a hand-built 1,000-read report run through
Bracken 3.0.1), *Streptococcus pneumoniae* has 150 reads at its own node (15%) and 450 reads sit at
the genus and phylum nodes above it. Bracken estimates 301 reads (31.7%), twice the species count.

``fraction_total_reads`` is a fraction in 0-1, not a percentage: Bracken's ``est_abundance.py``
writes ``"%0.5f" % (new_est_reads / sum_all_reads)``. The same fraction/percentage confusion has
caused errors with SingleM's ``read_fraction`` (a percentage), skani's alignment fraction and
CheckM's completeness. :data:`FRACTION_UNIT` records the unit in every row of the output table.

Bracken's exit status cannot be trusted. Its wrapper ends every error path with ``echo "ERROR:
..."`` and a bare ``exit``, which returns the ``echo``'s status of zero. A missing database, a
missing ``kmer_distrib`` for the read length and a missing input file all exit 0. This module
therefore runs its guards before the call and checks the file Bracken should have written.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: Lu J, Breitwieser FP, Thielen P, Salzberg SL. "Bracken: estimating species abundance in
#: metagenomics data." *PeerJ Computer Science* 3:e104 (2017). Checked against Crossref metadata
#: for 10.7717/peerj-cs.104.
CITATION = (
    "Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in "
    "metagenomics data. PeerJ Comput Sci 3:e104 (2017). doi:10.7717/peerj-cs.104"
)

#: The ranks offered, in descending order. Bracken's ``-l`` also accepts strain sub-levels (``S1``,
#: ``S2``, ...), which are not offered: a standard Kraken2 database has almost no nodes below
#: species, so ``-l S1`` gives a near-empty table that is easily misread as a negative result.
RANKS: Tuple[str, ...] = ("D", "P", "C", "O", "F", "G", "S")

#: Bracken's rank letters spelled out, for messages.
RANK_NAMES: Mapping[str, str] = {
    "D": "domain",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}

#: The columns Bracken writes, in order, taken from the ``o_file.write`` call in
#: ``est_abundance.py`` that writes the header (the README does not name them).
BRACKEN_COLUMNS: Tuple[str, ...] = (
    "name",
    "taxonomy_id",
    "taxonomy_lvl",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
    "fraction_total_reads",
)

#: The file ``bracken -r N`` reads inside the database directory. The read-length guard uses it to
#: list the lengths that are available.
KMER_DISTRIB = "database%dmers.kmer_distrib"
_KMER_DISTRIB_RE = re.compile(r"^database(\d+)mers\.kmer_distrib$")

#: The database's k, from the header of ``kraken2-inspect`` output (``inspect.txt``). A read
#: shorter than k has no k-mers and cannot be classified.
_KMER_LENGTH_RE = re.compile(r"\bk\s*=\s*(\d+)")

#: If less than this share of the movable reads is moved, the database and the reads do not match.
#: Above zero because Bracken moves nothing for a node whose children are all below the ``-t``
#: threshold; low because a single-organism sample has little to move. On the fixture under
#: tests/fixtures/bracken/ a matched database moves 100% of the 450 movable reads, a mismatched 0.
MIN_REDISTRIBUTED_FRACTION = 0.01

#: The unit written in the long abundance table, named as in Bracken. A fraction of the rank's
#: estimated reads, in 0-1.
FRACTION_UNIT = "fraction_total_reads"
#: Bracken's re-estimated read counts. Not comparable across samples without a depth correction.
READS_UNIT = "est_reads"

#: Columns of MetaWrap2's long-format abundance table, as written by `catalogue --quantify` and read
#: by `diff_abundance`, so a Bracken table needs no conversion.
LONG_COLUMNS: Tuple[str, ...] = ("sample", "feature", "value", "unit")


class BrackenError(ValueError):
    """A Bracken table, rank or database that cannot be used, with what was seen and expected."""


@dataclass(frozen=True)
class BrackenRow:
    """One taxon's row of a ``.bracken`` table.

    ``added_reads`` equals ``int(new_est_reads) - int(kraken_assigned_reads)``. It is parsed from
    the file because it shows whether Bracken redistributed anything.
    """

    name: str
    taxonomy_id: int
    taxonomy_lvl: str
    #: Reads Kraken2 put at this taxon's own node.
    kraken_assigned_reads: int
    #: Reads Bracken moved down to it from nodes above.
    added_reads: int
    #: The estimate: assigned plus added. This is the abundance number.
    new_est_reads: int
    #: Share of the rank's total estimated reads, a fraction in 0-1 (see the module docstring).
    fraction_total_reads: float

    @property
    def added_share(self) -> float:
        """Share of this taxon's estimate that Bracken added on top of Kraken2's assignments.

        A taxon whose estimate is mostly redistributed reads had few reads assigned directly.
        """
        if self.new_est_reads <= 0:
            return 0.0
        return self.added_reads / float(self.new_est_reads)


def check_rank(rank: str) -> str:
    """Normalise and validate a rank letter, or raise saying which ones exist."""
    normalised = (rank or "").strip().upper()
    if normalised not in RANKS:
        raise BrackenError(
            "rank %r is not supported; expected one of %s (%s)"
            % (
                rank,
                ",".join(RANKS),
                ", ".join("%s=%s" % (letter, RANK_NAMES[letter]) for letter in RANKS),
            )
        )
    return normalised


def read_table(path: str) -> List[BrackenRow]:
    """Parse a ``.bracken`` / ``*_bracken_species.txt`` table.

    The header must equal :data:`BRACKEN_COLUMNS`. Bracken's ``-w`` output is a six-column
    Kraken2-style report with no header, and parsing it positionally as the ``-o`` table gives
    plausible but wrong values.
    """
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if fields != BRACKEN_COLUMNS:
            raise BrackenError(
                "%s is not a Bracken abundance table: expected the columns %s, found %s. The "
                "six-column file with no header is Bracken's -w report, not its -o table"
                % (path, list(BRACKEN_COLUMNS), list(fields))
            )
        rows: List[BrackenRow] = []
        for number, record in enumerate(reader, start=2):
            rows.append(_row(path, number, record))

    _check_fractions(path, rows)
    return rows


def _row(path: str, line_number: int, record: Mapping[str, str]) -> BrackenRow:
    try:
        return BrackenRow(
            name=(record["name"] or "").strip(),
            taxonomy_id=int(record["taxonomy_id"]),
            taxonomy_lvl=(record["taxonomy_lvl"] or "").strip(),
            kraken_assigned_reads=int(record["kraken_assigned_reads"]),
            added_reads=int(record["added_reads"]),
            new_est_reads=int(record["new_est_reads"]),
            fraction_total_reads=float(record["fraction_total_reads"]),
        )
    except (KeyError, TypeError, ValueError) as problem:
        raise BrackenError("%s line %d is not a Bracken row: %s" % (path, line_number, problem))


def _check_fractions(path: str, rows: Sequence[BrackenRow]) -> None:
    """Check that the last column is a fraction in 0-1.

    A value above 1 means a percentage (100x too large). The fractions are shares of one rank's
    total, so they must also sum to about 1. The tolerance allows for rounding to five decimals.
    """
    if not rows:
        return
    largest = max(row.fraction_total_reads for row in rows)
    if largest > 1.0:
        raise BrackenError(
            "%s has fraction_total_reads=%g, above 1. Bracken writes a fraction in 0-1 "
            "(new_est_reads / total estimated reads at the rank); this column is a percentage "
            "or the file is not a Bracken table" % (path, largest)
        )
    total = sum(row.fraction_total_reads for row in rows)
    if not 0.95 <= total <= 1.05:
        raise BrackenError(
            "%s: fraction_total_reads sums to %.4f over %d row(s), expected 1. The table has "
            "been filtered, concatenated from several samples, "
            "or is not a Bracken table" % (path, total, len(rows))
        )


def check_table_rank(rows: Sequence[BrackenRow], rank: str, path: str = "table") -> None:
    """Raise unless every row was estimated at *rank*.

    Bracken writes the requested rank into ``taxonomy_lvl`` on every row. A mismatch usually comes
    from ``--resume`` over tables written with a different ``--bracken-rank``.
    """
    wanted = check_rank(rank)
    found = sorted({row.taxonomy_lvl for row in rows})
    if not found:
        return
    if found != [wanted]:
        raise BrackenError(
            "%s was estimated at rank %s; requested %s (%s). The rank comes from each row's "
            "taxonomy_lvl. Re-run with --force, or ask for rank %s"
            % (path, "/".join(found), wanted, RANK_NAMES[wanted], found[0])
        )


# ── database read lengths and k ───────────────────────────────────────────────────────────────


def available_read_lengths(db: str) -> List[int]:
    """The read lengths *db* has a Bracken k-mer distribution for, ascending.

    Empty means ``bracken-build`` has not been run on this Kraken2 database (a missing install
    step), which gets a different message from a wrong read length.
    """
    if not os.path.isdir(db):
        return []
    lengths: List[int] = []
    for entry in os.listdir(db):
        match = _KMER_DISTRIB_RE.match(entry)
        if match:
            lengths.append(int(match.group(1)))
    return sorted(lengths)


def database_kmer_length(db: str) -> Optional[int]:
    """The Kraken2 database's k, from the header ``kraken2-inspect`` writes. None if unavailable.

    Used to reject a read length below k: such a read has no k-mers and cannot be classified.
    """
    inspect = os.path.join(db, "inspect.txt")
    if not os.path.isfile(inspect):
        return None
    try:
        with open(inspect) as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                match = _KMER_LENGTH_RE.search(line)
                if match:
                    return int(match.group(1))
    except OSError:
        return None
    return None


def check_read_length(db: str, read_length: int) -> Optional[str]:
    """What is wrong with asking *db* for *read_length*, or None if nothing is.

    Returns a message instead of raising, so the calling module can report it with the sample
    through its own ``error``/``warning`` channel. Run it before Bracken, which exits zero on this
    error.
    """
    if not os.path.isdir(db):
        return (
            "the Bracken database directory %s does not exist. bracken-build writes its files "
            "into the Kraken2 database directory; pass that directory as --bracken-db" % db
        )

    available = available_read_lengths(db)
    if not available:
        return (
            "%s holds no Bracken k-mer distribution (no database<N>mers.kmer_distrib). Build "
            "one from this Kraken2 database with `bracken-build -d %s -l %d`"
            " (needed in addition to the Kraken2 database)" % (db, db, read_length)
        )

    kmer_length = database_kmer_length(db)
    if kmer_length is not None and read_length < kmer_length:
        return (
            "read length %d is shorter than this Kraken2 database's k of %d; reads that "
            "short contain no complete k-mer and cannot be classified" % (read_length, kmer_length)
        )

    if read_length not in available:
        return (
            "%s has no k-mer distribution for %d bp reads: it has %s (looked for %s). "
            "Distributions are built per read length. Pass --bracken-read-length with one of "
            "those, or build the missing one with "
            "`bracken-build -d %s -l %d`"
            % (
                db,
                read_length,
                ", ".join("%d bp" % length for length in available),
                KMER_DISTRIB % read_length,
                db,
                read_length,
            )
        )
    return None


# ── reads above the rank in the Kraken2 report ───────────────────────────────────────────────

#: Kraken2's report columns, by position (the file has no header).
_REPORT_RANK_COLUMN = 3
_REPORT_TAXON_READS_COLUMN = 2
_REPORT_COLUMNS = 6

#: Rank codes, coarsest first, copied from ``main_lvls`` in Bracken's ``est_abundance.py`` so "above
#: the rank" matches Bracken. Bracken puts ``K`` (NCBI kingdom) before ``D`` (NCBI superkingdom);
#: this follows Bracken's order, since the aim is to predict what Bracken will move.
_REPORT_RANK_ORDER: Tuple[str, ...] = ("R", "K", "D", "P", "C", "O", "F", "G", "S")


def reads_above_rank(report: str, rank: str) -> Optional[int]:
    """Reads Kraken2 assigned at nodes coarser than *rank*, which Bracken can redistribute.

    This is the denominator for the redistribution guard. The ``.bracken`` table has no rows for
    reads left above the rank, so without the report "Bracken added almost nothing" cannot be told
    apart from "there was almost nothing to add".

    Unranked nodes (code ``-``, such as "cellular organisms") count as movable: Bracken gives them a
    derived sub-level of their parent and redistributes them like other nodes above the rank.
    Unclassified reads (``U``) are not counted, since Bracken cannot place them. A sub-level code
    such as ``S1`` is below its letter and never counts.

    None when the report cannot be read; the guard then falls back to a weaker denominator.
    """
    target = _REPORT_RANK_ORDER.index(check_rank(rank))
    total = 0
    try:
        with open(report) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < _REPORT_COLUMNS:
                    continue
                code = fields[_REPORT_RANK_COLUMN].strip()
                letter = code[:1]
                if not letter or letter == "U":
                    continue
                if letter in _REPORT_RANK_ORDER and _REPORT_RANK_ORDER.index(letter) >= target:
                    continue
                try:
                    total += int(fields[_REPORT_TAXON_READS_COLUMN])
                except ValueError:
                    continue
    except OSError:
        return None
    return total


# ── redistribution check ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Redistribution:
    """How much work Bracken did on one sample, and whether that is a believable amount."""

    rank: str
    #: Reads Kraken2 had already placed at the rank.
    assigned: int
    #: Reads Bracken moved down to the rank.
    added: int
    #: The estimate at the rank: assigned + added.
    estimated: int
    #: Reads Kraken2 left above the rank, from the Kraken2 report. When None, :attr:`moved_share` is
    #: relative to the estimate instead, which is a weaker measure.
    movable: Optional[int] = None

    @property
    def moved_share(self) -> float:
        """Share of the movable reads that moved, or of the estimate when *movable* is unknown."""
        denominator = self.movable if self.movable is not None else self.estimated
        if not denominator:
            return 0.0
        return self.added / float(denominator)

    @property
    def suspicious(self) -> bool:
        """True when reads were available to move and nearly none of them moved.

        With no reads above the rank a share of zero is correct, so the guard does not fire.
        """
        if self.movable is not None and self.movable <= 0:
            return False
        if not self.estimated:
            return True
        return self.moved_share < MIN_REDISTRIBUTED_FRACTION

    def describe(self) -> str:
        basis = (
            "%d read(s) sat above %s" % (self.movable, self.rank)
            if self.movable is not None
            else "estimate of %d read(s)" % self.estimated
        )
        return "%s: Kraken2 assigned %d read(s) at rank %s, Bracken added %d (%.1f%% of the %s)" % (
            RANK_NAMES.get(self.rank, self.rank),
            self.assigned,
            self.rank,
            self.added,
            100.0 * self.moved_share,
            basis,
        )


def redistribution(
    rows: Sequence[BrackenRow], rank: str, movable: Optional[int] = None
) -> Redistribution:
    """Summarise what Bracken moved, for the guard below and for the run's log."""
    return Redistribution(
        rank=check_rank(rank),
        assigned=sum(row.kraken_assigned_reads for row in rows),
        added=sum(row.added_reads for row in rows),
        estimated=sum(row.new_est_reads for row in rows),
        movable=movable,
    )


def check_redistribution(summary: Redistribution) -> Optional[str]:
    """What is wrong with this redistribution, or None if nothing is.

    Catches a Bracken database built from a different Kraken2 database. Its k-mer distributions use
    taxids the report does not contain, so no reads move, yet Bracken writes a full table whose
    fractions sum to 1. The result is Kraken2's assigned counts renormalised, with no error.
    """
    if not summary.suspicious:
        return None
    if not summary.estimated:
        return (
            "Bracken estimated no reads at rank %s (%s). The Kraken2 report has no nodes at "
            "that rank, or all of them are below Bracken's read threshold"
            % (summary.rank, RANK_NAMES.get(summary.rank, summary.rank))
        )
    return (
        "%s. Almost no reads were redistributed, so these numbers are close to Kraken2's "
        "assigned counts renormalised. The usual cause is a Bracken "
        "database built from a different Kraken2 database than the one that produced the report: "
        "its k-mer distributions use taxids that are not in this report, and Bracken "
        "still writes a full table" % summary.describe()
    )


# ── MetaWrap2 long-format abundance table ─────────────────────────────────────────────────────


def long_rows(
    sample: str, rows: Sequence[BrackenRow], unit: str = FRACTION_UNIT
) -> List[List[str]]:
    """One sample's rows of the long abundance table: sample, feature, value, unit.

    The feature is the taxon name, not the taxid, as in the pipeline's other long tables, because
    `diff_abundance` reports this string to the user. Names are checked to be unique per taxid:
    two rows sharing a name would be summed into one feature by consumers that key on it.
    """
    if unit not in (FRACTION_UNIT, READS_UNIT):
        raise BrackenError(
            "unit %r is not supported; expected %r (a fraction in 0-1) or %r "
            "(Bracken's re-estimated read counts)" % (unit, FRACTION_UNIT, READS_UNIT)
        )
    seen: Dict[str, int] = {}
    out: List[List[str]] = []
    for row in rows:
        if row.name in seen and seen[row.name] != row.taxonomy_id:
            raise BrackenError(
                "sample %s has two taxa named %r (taxids %d and %d); the name must identify one "
                "taxon, or both would be summed into one row"
                % (sample, row.name, seen[row.name], row.taxonomy_id)
            )
        seen[row.name] = row.taxonomy_id
        value = row.fraction_total_reads if unit == FRACTION_UNIT else float(row.new_est_reads)
        out.append([sample, row.name, "%.6g" % value, unit])
    return out


def write_long_table(
    path: str, per_sample: Mapping[str, Sequence[BrackenRow]], unit: str = FRACTION_UNIT
) -> int:
    """Write the combined long-format abundance table across samples. Returns the row count.

    Samples and taxa are sorted so the file is byte-identical between runs on the same inputs and
    can be diffed.
    """
    written = 0
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(list(LONG_COLUMNS))
        for sample in sorted(per_sample):
            rows = sorted(per_sample[sample], key=lambda row: (-row.new_est_reads, row.name))
            for record in long_rows(sample, rows, unit):
                writer.writerow(record)
                written += 1
    return written
