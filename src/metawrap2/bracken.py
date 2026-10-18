"""Bracken: Kraken2's read assignments turned into an abundance estimate at one rank.

**A Kraken2 report is not an abundance table.** Its third column is the reads assigned *at* a node,
and Kraken2 assigns a read to the lowest node all of its k-mers agree on. A read drawn from sequence
conserved across a genus therefore stops at the genus node, so a species genuinely at 10% of the
community can show a small share of reads at its own node with the rest sitting at genus, family or
higher. Reading the species rows of a report as composition understates every species and invents a
large "genus-level" population that is not an organism.

Bracken (Lu et al., *PeerJ Comput Sci* 2017, 10.7717/peerj-cs.104) fixes this by redistributing the
stalled reads down to a chosen rank. For each node above the rank it asks, from the k-mer content of
the reference genomes themselves, how probable it is that a read of a given length from each
candidate genome would have stopped at that node, and splits the node's reads accordingly.

Two consequences, and they are the whole reason this file has guards rather than just a parser:

* **The k-mer distribution is a property of one Kraken2 database.** It is computed by mapping that
  database's own reference genomes back through it, so a distribution built from a different database
  describes different genomes. ``bracken-build`` writes its files *into* the Kraken2 database
  directory for exactly this reason, which is why ``--bracken-db`` defaults there.
* **It is also a property of one read length.** The probability that a read stops at the genus node
  depends on how many k-mers the read has, so the distribution is built per read length and stored as
  ``database<N>mers.kmer_distrib``. Asking for a length the database was not built for is not a
  degraded answer, it is no answer.

The size of the correction, on a hand-built 1,000-read report run through real Bracken 3.0.1 (the
fixture under ``tests/fixtures/bracken/``): *Streptococcus pneumoniae* has 150 reads at its own node,
15% of the sample, and 450 of the 1,000 reads sit at the genus and phylum nodes above it. Bracken
estimates 301 reads for it, 31.7%. The species-level read count was wrong by a factor of two.

**``fraction_total_reads`` is a fraction in 0-1, not a percentage.** Verified in Bracken's own
``est_abundance.py``, which writes ``"%0.5f" % (new_est_reads / sum_all_reads)``. This is stated
rather than inferred because the project has been bitten three times by the same confusion in other
tools' tables - SingleM's ``read_fraction`` (a percentage), skani's alignment fraction and CheckM's
completeness - and in each case the magnitude of a typical value made the wrong guess look right.
:data:`FRACTION_UNIT` carries the answer into every row of the output table.

One thing this file will not do: Bracken has no exit status worth trusting. Its wrapper script ends
every error path with ``echo "ERROR: ..."`` followed by a bare ``exit``, which exits with the status
of the ``echo`` - zero. A missing database, a missing ``kmer_distrib`` for the requested read length
and a missing input file all look like success to the caller. Everything here therefore checks the
file Bracken was supposed to write, and the guards below run before the call rather than after it.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: Lu J, Breitwieser FP, Thielen P, Salzberg SL. "Bracken: estimating species abundance in
#: metagenomics data." *PeerJ Computer Science* 3:e104 (2017). Verified against Crossref metadata
#: for 10.7717/peerj-cs.104 rather than quoted from memory.
CITATION = (
    "Lu J, Breitwieser FP, Thielen P, Salzberg SL. Bracken: estimating species abundance in "
    "metagenomics data. PeerJ Comput Sci 3:e104 (2017). doi:10.7717/peerj-cs.104"
)

#: The ranks offered, in descending order. Bracken's ``-l`` also takes ``S1``, ``S2`` and so on for
#: strain sub-levels, which are deliberately not offered: a Kraken2 standard database has almost no
#: nodes below species, so ``-l S1`` on one produces a near-empty table rather than strain calls, and
#: a flag whose usual answer is "nothing found" is a flag that gets misread as a negative result.
RANKS: Tuple[str, ...] = ("D", "P", "C", "O", "F", "G", "S")

#: Bracken's rank letters spelled out, for messages. A guard that says "rank K" helps nobody.
RANK_NAMES: Mapping[str, str] = {
    "D": "domain",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}

#: The columns Bracken writes, in order. Taken from the single ``o_file.write`` call in
#: ``est_abundance.py`` that emits the header, not from documentation - Bracken's README describes
#: the columns in prose and does not name them.
BRACKEN_COLUMNS: Tuple[str, ...] = (
    "name",
    "taxonomy_id",
    "taxonomy_lvl",
    "kraken_assigned_reads",
    "added_reads",
    "new_est_reads",
    "fraction_total_reads",
)

#: ``bracken -r N`` resolves to this file inside the database directory, and nothing else. Knowing
#: the name is what lets the read-length guard say which lengths *are* available instead of only
#: that the requested one is not.
KMER_DISTRIB = "database%dmers.kmer_distrib"
_KMER_DISTRIB_RE = re.compile(r"^database(\d+)mers\.kmer_distrib$")

#: ``kraken2-inspect``/``inspect.txt`` opens with this, and it carries the database's k. A read
#: shorter than k yields no k-mers at all, so such a request is meaningless rather than imprecise.
_KMER_LENGTH_RE = re.compile(r"\bk\s*=\s*(\d+)")

#: Below this share of the movable reads actually being moved, the database and the reads disagree.
#: Not zero, because Bracken legitimately moves nothing for a node whose children are all below the
#: ``-t`` threshold; and not high, because a clean single-organism sample has little to move. On the
#: fixture under tests/fixtures/bracken/ a matched database moves 100% of the 450 movable reads and
#: a mismatched one moves 0.
MIN_REDISTRIBUTED_FRACTION = 0.01

#: The unit written in the long abundance table, spelled the way Bracken spells it so that nothing
#: downstream has to guess. A fraction of the rank's estimated reads, in 0-1.
FRACTION_UNIT = "fraction_total_reads"
#: The alternative: Bracken's re-estimated read counts. Not compositional, but also not comparable
#: across samples without a depth correction, since a deeper library has more reads of everything.
READS_UNIT = "est_reads"

#: Columns of MetaWrap2's long-format abundance table. The same four `catalogue --quantify` writes
#: and `diff_abundance` reads, in the same order, so a Bracken table can be fed straight into the
#: statistics module without a conversion step.
LONG_COLUMNS: Tuple[str, ...] = ("sample", "feature", "value", "unit")


class BrackenError(ValueError):
    """A Bracken table, rank or database that cannot be used, with what was seen and expected."""


@dataclass(frozen=True)
class BrackenRow:
    """One taxon's row of a ``.bracken`` table.

    ``added_reads`` is redundant - Bracken writes ``int(new_est_reads) - int(kraken_assigned_reads)``
    - but it is parsed rather than recomputed because it is the one column that says whether Bracken
    did anything at all, and a column read from the file can be checked against the arithmetic.
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
    #: Share of the rank's total estimated reads. A fraction in 0-1 - see the module docstring.
    fraction_total_reads: float

    @property
    def added_share(self) -> float:
        """What share of this taxon's estimate Bracken supplied rather than Kraken2.

        Worth having per row and not only in aggregate: a taxon whose estimate is 95% redistributed
        reads was barely observed directly, and that is a different kind of claim from one Kraken2
        assigned outright, even though the two look identical in the abundance column.
        """
        if self.new_est_reads <= 0:
            return 0.0
        return self.added_reads / float(self.new_est_reads)


def check_rank(rank: str) -> str:
    """Normalise and validate a rank letter, or raise saying which ones exist."""
    normalised = (rank or "").strip().upper()
    if normalised not in RANKS:
        raise BrackenError(
            "rank %r is not one this module estimates at; expected one of %s (%s)"
            % (
                rank,
                ",".join(RANKS),
                ", ".join("%s=%s" % (letter, RANK_NAMES[letter]) for letter in RANKS),
            )
        )
    return normalised


def read_table(path: str) -> List[BrackenRow]:
    """Parse a ``.bracken`` / ``*_bracken_species.txt`` table.

    The header is required to be exactly :data:`BRACKEN_COLUMNS`. Checked rather than trusted
    because the two files Bracken writes look alike and are not: ``-o`` is this table and ``-w`` is a
    six-column Kraken2-style report with no header, so a swapped pair parses into plausible nonsense
    if columns are taken positionally.
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
    """Confirm the last column really is a fraction, with the file in front of us.

    A value above 1 means the column is a percentage, which would make every abundance derived from
    it 100x too large; and since Bracken's fractions are shares of one rank's total they must also
    sum to about 1, so a sum near 100 catches the same mistake from the other side. The tolerance is
    loose because the column is rounded to five decimals, so a table with many rows loses a little.
    """
    if not rows:
        return
    largest = max(row.fraction_total_reads for row in rows)
    if largest > 1.0:
        raise BrackenError(
            "%s has fraction_total_reads=%g, which is above 1. Bracken writes a fraction in 0-1 "
            "(new_est_reads / total estimated reads at the rank), so this column is a percentage "
            "or this is not a Bracken table" % (path, largest)
        )
    total = sum(row.fraction_total_reads for row in rows)
    if not 0.95 <= total <= 1.05:
        raise BrackenError(
            "%s: fraction_total_reads sums to %.4f over %d row(s), not to 1. Bracken's fractions "
            "are shares of one rank's estimated reads, so a table whose fractions do not sum to 1 "
            "has been filtered, concatenated from several samples, or is not a Bracken table"
            % (path, total, len(rows))
        )


def check_table_rank(rows: Sequence[BrackenRow], rank: str, path: str = "table") -> None:
    """Raise unless every row was estimated at *rank*.

    Bracken stamps the rank it was asked for into ``taxonomy_lvl`` on every row, which makes this
    cheap and makes it worth doing: the common way to end up with the wrong one is ``--resume`` over
    an output directory whose tables were written at a different ``--bracken-rank``, and a genus
    table read as a species table is wrong without looking wrong.
    """
    wanted = check_rank(rank)
    found = sorted({row.taxonomy_lvl for row in rows})
    if not found:
        return
    if found != [wanted]:
        raise BrackenError(
            "%s was estimated at rank %s, not at %s (%s) as requested. Bracken writes the rank it "
            "ran at into every row's taxonomy_lvl; re-run with --force, or ask for rank %s"
            % (path, "/".join(found), wanted, RANK_NAMES[wanted], found[0])
        )


# ── the database: which read lengths it was built for, and what k it used ─────────────────────


def available_read_lengths(db: str) -> List[int]:
    """The read lengths *db* has a Bracken k-mer distribution for, ascending.

    Empty means ``bracken-build`` has never run against this Kraken2 database. That is a different
    failure from a wrong read length and gets a different message: one is a missing install step,
    the other is a wrong flag.
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

    Only used to reject a read length below k, which is not a tuning question: Kraken2 classifies on
    k-mers, so a read shorter than k produces none and cannot be classified at all.
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

    Returns a message rather than raising, because the caller is a module that wants to report this
    through its own ``error``/``warning`` channel alongside the sample it happened on. Bracken
    itself prints the same problem and then exits zero, so this has to run before the call.
    """
    if not os.path.isdir(db):
        return (
            "the Bracken database directory %s does not exist. bracken-build writes its files into "
            "the Kraken2 database directory, so --bracken-db should name that directory" % db
        )

    available = available_read_lengths(db)
    if not available:
        return (
            "%s holds no Bracken k-mer distribution (no database<N>mers.kmer_distrib). Bracken "
            "needs one built from this Kraken2 database with `bracken-build -d %s -l %d`; the "
            "Kraken2 database alone is not enough" % (db, db, read_length)
        )

    kmer_length = database_kmer_length(db)
    if kmer_length is not None and read_length < kmer_length:
        return (
            "read length %d is shorter than this Kraken2 database's k of %d, so a read of that "
            "length contains no complete k-mer and cannot be classified at all"
            % (read_length, kmer_length)
        )

    if read_length not in available:
        return (
            "%s has no k-mer distribution for %d bp reads: it expects %s (looked for %s). The "
            "distribution is built per read length, so this is not an approximation Bracken can "
            "make - pass --bracken-read-length with one of those, or build the missing one with "
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


# ── the Kraken2 report: how many reads there were to move ────────────────────────────────────

#: Kraken2's report columns. No header in the file, so they are positional by necessity.
_REPORT_RANK_COLUMN = 3
_REPORT_TAXON_READS_COLUMN = 2
_REPORT_COLUMNS = 6

#: Rank codes ordered coarsest first, copied from ``main_lvls`` in Bracken's ``est_abundance.py`` so
#: that this module's notion of "above the rank" is the same as the tool's. Note ``K`` before ``D``:
#: NCBI's superkingdom is Kraken2's ``D`` and its kingdom is ``K``, and Bracken orders them this way.
#: Agreeing with Bracken matters more here than being right about taxonomy, since the point of the
#: number is to predict what Bracken will try to move.
_REPORT_RANK_ORDER: Tuple[str, ...] = ("R", "K", "D", "P", "C", "O", "F", "G", "S")


def reads_above_rank(report: str, rank: str) -> Optional[int]:
    """Reads Kraken2 assigned at nodes coarser than *rank*: the reads Bracken exists to move.

    This is the denominator the redistribution guard needs and the ``.bracken`` table cannot supply.
    A table shows the reads that arrived at the rank; it has no row for one that stayed at the family
    node, so "Bracken added almost nothing" is indistinguishable from "there was almost nothing to
    add" without reading the report as well.

    Unranked nodes (code ``-``, such as "cellular organisms") count as movable: Bracken gives them a
    derived sub-level of their parent and redistributes their reads like any other node above the
    rank, and a read parked there is as stuck as one parked on a phylum. Unclassified reads (``U``)
    do not - Bracken cannot place them, and counting them would make every sample look like a
    failure. A sub-level code such as ``S1`` is below its letter, so it is never movable upward.

    None when the report cannot be read, which leaves the guard with a weaker denominator rather
    than a wrong one.
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


# ── did Bracken actually redistribute anything? ──────────────────────────────────────────────


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
    #: Reads Kraken2 left above the rank, from the Kraken2 report. None when it was not read, and
    #: then :attr:`moved_share` is of the estimate instead - a weaker number, flagged as such.
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

        Both halves matter. With nothing above the rank there is nothing to redistribute and a share
        of zero is correct, so the guard must not fire on a sample that was simply easy to classify.
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

    The failure this catches is the quiet one. A Bracken database built from a *different* Kraken2
    database has k-mer distributions keyed by taxids the report never mentions, so every lookup
    misses, every node's reads stay where they were, and Bracken still writes a full table whose
    fractions still sum to 1. Nothing in the output says the estimate is just Kraken2's assigned
    counts renormalised - which is the number the user ran Bracken to avoid.
    """
    if not summary.suspicious:
        return None
    if not summary.estimated:
        return (
            "Bracken estimated no reads at rank %s (%s). Either the Kraken2 report has no nodes at "
            "that rank, or every one of them fell below Bracken's read threshold"
            % (summary.rank, RANK_NAMES.get(summary.rank, summary.rank))
        )
    return (
        "%s. Almost nothing was redistributed, so these numbers are close to Kraken2's assigned "
        "counts renormalised rather than an abundance estimate. The usual cause is a Bracken "
        "database built from a different Kraken2 database than the one that produced the report: "
        "its k-mer distributions are keyed by taxids this report never mentions, so every lookup "
        "misses and Bracken still writes a complete-looking table" % summary.describe()
    )


# ── MetaWrap2's own long-format abundance table ───────────────────────────────────────────────


def long_rows(
    sample: str, rows: Sequence[BrackenRow], unit: str = FRACTION_UNIT
) -> List[List[str]]:
    """One sample's rows of the long abundance table: sample, feature, value, unit.

    The feature is the taxon name rather than its taxid, matching every other long table in the
    pipeline in being readable - `diff_abundance` reports whichever string it finds here straight
    back to the user, and a list of significant taxids is not a result anyone can read. Names from
    one Kraken2 database are unique per taxid, but that is checked rather than assumed: two rows
    sharing a name would be silently summed into one feature by any consumer that keys on it.
    """
    if unit not in (FRACTION_UNIT, READS_UNIT):
        raise BrackenError(
            "unit %r is not one this writer produces; expected %r (a fraction in 0-1) or %r "
            "(Bracken's re-estimated read counts)" % (unit, FRACTION_UNIT, READS_UNIT)
        )
    seen: Dict[str, int] = {}
    out: List[List[str]] = []
    for row in rows:
        if row.name in seen and seen[row.name] != row.taxonomy_id:
            raise BrackenError(
                "sample %s has two taxa named %r (taxids %d and %d), so the name cannot identify a "
                "feature here; they would be summed into one row"
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

    Samples and taxa are both sorted so the file is byte-identical between runs over the same
    inputs: a table that reorders itself cannot be diffed, which is how a changed number gets
    missed in review.
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
