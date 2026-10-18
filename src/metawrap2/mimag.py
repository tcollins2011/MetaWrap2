"""MIMAG quality tags for recovered genomes (Bowers et al., Nat Biotechnol 2017;35:725-731).

CheckM and CheckM2 answer two of MIMAG's three questions: how complete a bin is and how contaminated.
The third is whether the genome carries the ribosomal and transfer RNA genes a submitted MAG is
expected to have - 5S, 16S and 23S rRNA, and tRNAs for at least 18 of the 20 amino acids. That
evidence is not in a CheckM report at all, so a tier assigned from completeness and contamination
alone is MIMAG high quality *minus its hardest criterion*, which is why
:data:`metawrap2.benchmark.scoring.QUALITY_TIERS` calls its 90/5 tier ``bins_near_complete`` rather
than high quality.

This module adds the missing half. ``annotate --bakta`` already runs tRNAscan-SE, Aragorn and
Infernal/Rfam over every genome, so where Bakta has run the counts are a parse of a table the
pipeline has already written; where it has not, there is no rRNA evidence either way.

That last case is the whole reason this module is written carefully rather than as three
comparisons. A bin at 95% complete and 2% contaminated with no Bakta output has not been shown to
lack a 16S - nobody looked. Calling it medium quality would be a measurement this pipeline never
made, reported as a negative result, and it would be indistinguishable in the output from a genome
where Infernal searched and found nothing. The first is a gap in the run; the second is a fact about
the genome, and only the second is a reason to downgrade a bin. So an unmeasured genome that clears
90/5 is tagged :data:`MimagTag.HIGH_UNVERIFIED` with ``mimag_rrna=unknown`` beside it, and a measured
absence is tagged :data:`MimagTag.MEDIUM` with the genes that were found spelled out. A missing
measurement is not a missing gene.

Nothing here filters bins, for the same reason GUNC does not (see :mod:`metawrap2.gunc`): the tag is
a column, and what to submit is the reader's decision.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from enum import Enum
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .constants import MIMAG_STATS_COLUMNS

#: MIMAG high quality: >90% complete, <5% contaminated, 5S + 16S + 23S present, >=18 tRNAs.
#:
#: The comparisons are strict where Bowers et al. Table 1 is strict and inclusive where it is
#: inclusive, which is asymmetric on purpose: high quality is ">90% completion" and "<5%
#: contamination", medium is ">=50% completion". So 90.0 exactly is not MIMAG high quality and 50.0
#: exactly is MIMAG medium.
#:
#: That differs by one boundary from :func:`metawrap2.benchmark.scoring.quality_counts`, which
#: compares ``>= 90`` and ``<= 5`` on the same two numbers. The difference only ever shows on a bin
#: scored exactly 90.00 or exactly 5.00, and the two are not interchangeable: the benchmark counts
#: tiers for comparing bin sets, while this tag is what a reader would put in a submission, so it
#: follows the paper's own inequalities rather than the benchmark's.
HIGH_MIN_COMPLETENESS = 90.0
HIGH_MAX_CONTAMINATION = 5.0

#: MIMAG medium quality: >=50% complete, <10% contaminated. Low quality is <50% complete with the
#: same contamination ceiling.
MEDIUM_MIN_COMPLETENESS = 50.0

#: Contamination at or above this is outside MIMAG entirely - there is no tier for it. Bowers et al.
#: put a <10% contamination ceiling on *every* tier including low quality, so a 30%-contaminated bin
#: is not "low quality", it is unclassifiable, and saying "low" would promote it.
MAX_CONTAMINATION = 10.0

#: tRNAs for at least 18 of the 20 amino acids. Counted as distinct amino acids rather than as tRNA
#: features, because that is what the criterion says: a genome with 40 tRNA genes covering 15 amino
#: acids does not satisfy it, and a feature count would pass it comfortably.
MIN_TRNA_TYPES = 18

#: The three rRNA genes MIMAG high quality requires, in the order they are reported.
RRNA_GENES = ("5S", "16S", "23S")

#: The 20 amino acids MIN_TRNA_TYPES counts against. Selenocysteine and pyrrolysine are deliberately
#: absent: they are the 21st and 22nd, not part of "18 of the 20", and tRNAscan-SE reports SeC often
#: enough that counting it would hand two free amino acids to any genome that has one.
AMINO_ACIDS = frozenset(
    (
        "Ala",
        "Arg",
        "Asn",
        "Asp",
        "Cys",
        "Gln",
        "Glu",
        "Gly",
        "His",
        "Ile",
        "Leu",
        "Lys",
        "Met",
        "Phe",
        "Pro",
        "Ser",
        "Thr",
        "Trp",
        "Tyr",
        "Val",
    )
)

#: The columns this module appends to a ``.stats.tsv``, in order. Appended after whatever is already
#: there and never reordered, because benchmark tables and users' own scripts read the
#: base columns by position (see :data:`metawrap2.checkm.STATS_HEADER`). Declared in
#: :mod:`metawrap2.constants` beside the suffix of the file they go into, so a reader has one place
#: to look for the table's full shape, and re-exported here because this is the module that writes
#: them.
STATS_COLUMNS = MIMAG_STATS_COLUMNS

#: Written in ``mimag_rrna`` when no rRNA source was found, so that "we did not look" and "we looked
#: and found none" are different cells rather than both being blank.
UNKNOWN = "unknown"

#: Said when a bin set has no rRNA evidence at all, which is the normal state in ``bin_refinement``:
#: refinement runs before ``annotate``, so at that point nothing has searched for a 16S.
NO_EVIDENCE_HINT = (
    "No rRNA/tRNA evidence was found, so bins that clear 90/5 are tagged %s rather than "
    "high quality. Run 'metawrap2 annotate --bakta' over the final bins and re-run this step to "
    "resolve them." % "high_quality_unverified_rrna"
)

#: Bakta names a tRNA by amino acid in its product column: ``tRNA-Ala(tgc)``. The trailing digit in
#: ``tRNA-Ile2`` is an isotype, not a different amino acid, so the capture stops at three letters.
_BAKTA_TRNA = re.compile(r"tRNA-([A-Za-z]{3})")

#: Bakta's rRNA products, e.g. ``16S ribosomal RNA``. Matched on the subunit prefix because the rest
#: of the string varies ("16S ribosomal RNA" vs "16S ribosomal RNA (partial)").
_BAKTA_RRNA = re.compile(r"\b(5S|16S|23S)\b")

#: Rfam family names for the three subunits, as they appear in an Infernal ``--tblout`` target
#: column. ``5_8S_rRNA`` is a eukaryotic family that must not satisfy the 5S criterion, so it is
#: excluded explicitly rather than left to a substring match on "5S".
_RFAM_FAMILIES = {
    "5S": ("5S_rRNA",),
    "16S": ("SSU_rRNA_bacteria", "SSU_rRNA_archaea", "SSU_rRNA_eukarya", "SSU_rRNA_microsporidia"),
    "23S": ("LSU_rRNA_bacteria", "LSU_rRNA_archaea", "LSU_rRNA_eukarya"),
}

__all__ = [
    "AMINO_ACIDS",
    "HIGH_MAX_CONTAMINATION",
    "HIGH_MIN_COMPLETENESS",
    "MAX_CONTAMINATION",
    "MEDIUM_MIN_COMPLETENESS",
    "MIN_TRNA_TYPES",
    "NO_EVIDENCE_HINT",
    "RRNA_GENES",
    "STATS_COLUMNS",
    "UNKNOWN",
    "Assessment",
    "MimagTag",
    "RnaEvidence",
    "annotate_stats",
    "assess",
    "collect_evidence",
    "evidence_for_genome",
    "parse_bakta_gff",
    "parse_bakta_tsv",
    "parse_infernal_tblout",
    "parse_trnascan",
    "summarise",
    "tally",
]


class MimagTag(str, Enum):
    """What MIMAG tier a genome is in, plus the two states MIMAG itself does not name.

    A ``str`` enum so the value can go straight into a TSV cell and be compared to a parsed one
    without a conversion at each site.
    """

    #: >90 / <5 with 5S, 16S, 23S and >=18 tRNA types all measured present.
    HIGH = "high_quality"
    #: >90 / <5, but no rRNA/tRNA measurement exists to confirm or deny the rest. Not a MIMAG tier -
    #: it is this pipeline saying which half of the criteria it checked.
    HIGH_UNVERIFIED = "high_quality_unverified_rrna"
    #: >=50 / <10. Also where a genome that clears 90/5 lands once a *measured* rRNA or tRNA absence
    #: rules out high quality.
    MEDIUM = "medium_quality"
    #: <50 / <10.
    LOW = "low_quality"
    #: Contamination at or above 10%: no MIMAG tier applies. See :data:`MAX_CONTAMINATION`.
    NOT_MIMAG = "not_mimag"
    #: No completeness or contamination to judge by. An absence, not a zero.
    UNSCORED = "unscored"


@dataclass(frozen=True)
class RnaEvidence:
    """What was found when something looked for a genome's rRNA and tRNA genes.

    Every field is ``Optional`` and ``None`` means "not measured", never "absent". A run with
    Infernal output but no tRNAscan-SE output has the three rRNA fields set and ``trna_types`` at
    ``None``, and that genome's tRNA criterion is unverified rather than failed.
    """

    rrna_5s: Optional[bool] = None
    rrna_16s: Optional[bool] = None
    rrna_23s: Optional[bool] = None
    trna_types: Optional[int] = None
    #: Which file the evidence came from, for the stats table and for a reader asking why.
    source: str = ""

    @property
    def rrna(self) -> Tuple[Optional[bool], Optional[bool], Optional[bool]]:
        """The three rRNA flags in :data:`RRNA_GENES` order."""
        return (self.rrna_5s, self.rrna_16s, self.rrna_23s)

    def describe_rrna(self) -> str:
        """The ``mimag_rrna`` cell: ``5S+16S+23S``, ``16S``, ``none``, or ``unknown``.

        ``none`` and ``unknown`` are the two cells this whole module exists to keep apart - the
        first is a search that found nothing, the second is no search.
        """
        if all(flag is None for flag in self.rrna):
            return UNKNOWN
        found = [gene for gene, flag in zip(RRNA_GENES, self.rrna) if flag]
        if not found:
            return "none" if all(flag is not None for flag in self.rrna) else UNKNOWN
        missing_unmeasured = any(flag is None for flag in self.rrna)
        return "+".join(found) + ("+%s" % UNKNOWN if missing_unmeasured else "")


@dataclass(frozen=True)
class Assessment:
    """One genome's MIMAG verdict and everything it was decided from."""

    bin: str
    completeness: Optional[float]
    contamination: Optional[float]
    evidence: RnaEvidence
    tag: MimagTag

    @property
    def rrna_measured(self) -> bool:
        """True if something actually searched this genome for all three rRNA genes."""
        return all(flag is not None for flag in self.evidence.rrna)

    @property
    def trna_measured(self) -> bool:
        return self.evidence.trna_types is not None


def assess(
    name: str,
    completeness: Optional[float],
    contamination: Optional[float],
    evidence: Optional[RnaEvidence] = None,
) -> Assessment:
    """The MIMAG tier for one genome.

    The order of the checks is the order of the guards that matter:

    1. no completeness or contamination at all -> :attr:`MimagTag.UNSCORED`. A bin CheckM could not
       score has not scored zero, and 0/100 would read as a complete failure of the genome rather
       than of the measurement.
    2. contamination >= 10% -> :attr:`MimagTag.NOT_MIMAG`, because every MIMAG tier including low
       quality has that ceiling.
    3. >90 / <5 -> high quality, unverified high quality, or medium, depending on what the rRNA and
       tRNA evidence says and on whether it exists. A *measured* absence demotes to medium; a
       *missing* measurement leaves the genome at :attr:`MimagTag.HIGH_UNVERIFIED`.
    4. >=50 -> medium, else low.
    """
    found = evidence if evidence is not None else RnaEvidence()
    if completeness is None or contamination is None:
        return Assessment(name, completeness, contamination, found, MimagTag.UNSCORED)
    if contamination >= MAX_CONTAMINATION:
        return Assessment(name, completeness, contamination, found, MimagTag.NOT_MIMAG)

    if completeness > HIGH_MIN_COMPLETENESS and contamination < HIGH_MAX_CONTAMINATION:
        # A criterion only fails on a measurement. ``flag is False`` and ``flag is None`` are
        # deliberately not folded into ``not flag``: that one-character shortcut is exactly how a
        # genome nobody annotated gets demoted to medium and becomes indistinguishable from one
        # where Infernal searched and found no 16S.
        failed = any(flag is False for flag in found.rrna) or (
            found.trna_types is not None and found.trna_types < MIN_TRNA_TYPES
        )
        unmeasured = any(flag is None for flag in found.rrna) or found.trna_types is None
        if failed:
            tag = MimagTag.MEDIUM
        elif unmeasured:
            tag = MimagTag.HIGH_UNVERIFIED
        else:
            tag = MimagTag.HIGH
        return Assessment(name, completeness, contamination, found, tag)

    tag = MimagTag.MEDIUM if completeness >= MEDIUM_MIN_COMPLETENESS else MimagTag.LOW
    return Assessment(name, completeness, contamination, found, tag)


# --- Parsing the evidence the pipeline already writes -----------------------------------------


def _trna_count(amino_acids: Iterable[str]) -> int:
    """How many of the 20 amino acids are covered, from raw tRNA labels.

    Normalised to ``Xxx`` capitalisation first, because the three sources spell them three ways:
    Bakta writes ``tRNA-Ala``, tRNAscan-SE writes ``Ala`` in its Type column, and some versions
    upper-case it. Anything not in :data:`AMINO_ACIDS` - ``Undet``, ``Pseudo``, ``SeC``, ``Sup`` -
    is dropped rather than counted, so a genome whose tRNAs are mostly undetermined does not clear
    18 on labels that name no amino acid.
    """
    covered = {
        label.strip().capitalize() for label in amino_acids if label and label.strip()
    } & AMINO_ACIDS
    return len(covered)


def parse_bakta_tsv(path: str) -> RnaEvidence:
    """rRNA and tRNA evidence from a Bakta ``<genome>.tsv``.

    Bakta's table is one row per feature with a ``Type`` column (``cds``, ``tRNA``, ``rRNA``,
    ``ncRNA``, ...) and a ``Product``. The header line starts with ``#Sequence Id`` and is preceded
    by several ``#``-comment lines carrying the Bakta and database versions, so the columns are
    located by name from the first ``#``-line that has them rather than by counting lines - the
    number of preamble lines differs between Bakta versions.

    Both rRNA and tRNA come from this one file, so the returned evidence has all four fields set:
    Bakta having run is itself the measurement, and a genome with no rRNA row has a *measured*
    absence of rRNA, which is the case that may legitimately demote a 95/2 bin to medium.
    """
    rrna: set[str] = set()
    trnas: List[str] = []
    columns: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                lowered = [cell.lstrip("#").strip().lower() for cell in cut]
                if "type" in lowered and "product" in lowered:
                    columns = {name: i for i, name in enumerate(lowered)}
                continue
            if not columns or len(cut) <= columns["product"]:
                continue
            kind = cut[columns["type"]].strip().lower()
            product = cut[columns["product"]]
            if kind == "rrna":
                rrna.update(_BAKTA_RRNA.findall(product))
            elif kind == "trna":
                trnas.extend(_BAKTA_TRNA.findall(product))
                gene = columns.get("gene")
                if not trnas and gene is not None and len(cut) > gene:
                    trnas.extend(_BAKTA_TRNA.findall(cut[gene]))
    return RnaEvidence(
        rrna_5s="5S" in rrna,
        rrna_16s="16S" in rrna,
        rrna_23s="23S" in rrna,
        trna_types=_trna_count(trnas),
        source=os.path.basename(path),
    )


def parse_bakta_gff(path: str) -> RnaEvidence:
    """The same evidence from a Bakta ``<genome>.gff3``, for a run that kept only the GFF.

    Column 3 is the feature type and column 9 the attributes, where the subunit and amino acid are
    in ``product=``/``Name=``. Read as a fallback and not as the primary source because the TSV's
    columns are fixed while GFF attribute keys vary with Bakta's version.
    """
    rrna: set[str] = set()
    trnas: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 9:
                continue
            kind = cut[2].strip().lower()
            attributes = cut[8]
            if kind == "rrna":
                rrna.update(_BAKTA_RRNA.findall(attributes))
            elif kind == "trna":
                trnas.extend(_BAKTA_TRNA.findall(attributes))
    return RnaEvidence(
        rrna_5s="5S" in rrna,
        rrna_16s="16S" in rrna,
        rrna_23s="23S" in rrna,
        trna_types=_trna_count(trnas),
        source=os.path.basename(path),
    )


def parse_trnascan(path: str) -> Optional[int]:
    """How many amino acids a tRNAscan-SE table covers, or None if it names none.

    tRNAscan-SE's default output is three header lines then whitespace-separated columns, of which
    the fifth is the tRNA type (``Ala``, ``Undet``, ``Pseudo``). Rows are taken by position because
    the header is wrapped across two lines and cannot be parsed as one.

    None rather than 0 when nothing in the file names an amino acid: a file that is present but
    holds no usable rows is more likely a truncated or unexpected format than a genome with no
    tRNAs, and reporting 0 would demote a bin on it.
    """
    types: List[str] = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("Sequence", "Name", "---")):
                continue
            cut = line.split()
            if len(cut) < 5:
                continue
            types.append(cut[4])
    count = _trna_count(types)
    return count if count else None


def parse_infernal_tblout(path: str) -> RnaEvidence:
    """rRNA presence from an Infernal ``cmscan --tblout`` table.

    The target name is the Rfam family, so the three subunits are recognised by family rather than
    by a product string (see :data:`_RFAM_FAMILIES`). Only the rRNA fields are set: an Infernal run
    says nothing about tRNAs, and leaving ``trna_types`` at None keeps that criterion unverified
    instead of failed.

    ``--tblout`` is space-padded to column widths rather than delimited, which is why the line is
    split on whitespace; the target name is the first field, or the second in the ``--fmt 2`` layout
    where an index comes first.
    """
    found: set[str] = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cut = line.split()
            candidates = cut[:2] if len(cut) > 1 else cut[:1]
            for gene, families in _RFAM_FAMILIES.items():
                if any(field in families for field in candidates):
                    found.add(gene)
    return RnaEvidence(
        rrna_5s="5S" in found,
        rrna_16s="16S" in found,
        rrna_23s="23S" in found,
        source=os.path.basename(path),
    )


#: Where the fallback rRNA/tRNA files are looked for, relative to the search directory and with
#: ``%s`` standing in for the genome name. These are the names a standalone Infernal/tRNAscan-SE
#: pass over bins is expected to use; nothing in the pipeline writes them today, so they only ever
#: match for a run that produced them deliberately.
_TRNASCAN_NAMES = ("rna/%s.trnascan.tsv", "rna/%s.trnascan.txt", "%s.trnascan.tsv")
_INFERNAL_NAMES = ("rna/%s.infernal.tblout", "rna/%s.cmscan.tblout", "%s.infernal.tblout")


def evidence_for_genome(search_dir: str, name: str) -> Optional[RnaEvidence]:
    """Whatever rRNA/tRNA evidence exists for genome *name* under *search_dir*, or None.

    *search_dir* is an ``annotate`` output directory. Bakta is preferred because one file answers
    both criteria; failing that, a separate Infernal table and tRNAscan-SE table are merged, each
    contributing only what it measured. None means no source was found at all, which is the state
    :attr:`MimagTag.HIGH_UNVERIFIED` exists for - not an error, and not an absence of genes.
    """
    bakta = os.path.join(search_dir, "bakta", name)
    for parse, path in (
        (parse_bakta_tsv, os.path.join(bakta, name + ".tsv")),
        (parse_bakta_gff, os.path.join(bakta, name + ".gff3")),
    ):
        if os.path.isfile(path):
            return parse(path)

    rrna: Optional[RnaEvidence] = None
    for template in _INFERNAL_NAMES:
        path = os.path.join(search_dir, template % name)
        if os.path.isfile(path):
            rrna = parse_infernal_tblout(path)
            break
    trnas: Optional[int] = None
    trna_source = ""
    for template in _TRNASCAN_NAMES:
        path = os.path.join(search_dir, template % name)
        if os.path.isfile(path):
            trnas = parse_trnascan(path)
            trna_source = os.path.basename(path)
            break
    if rrna is None and trnas is None:
        return None
    base = rrna if rrna is not None else RnaEvidence()
    sources = [part for part in (base.source, trna_source) if part]
    return RnaEvidence(
        rrna_5s=base.rrna_5s,
        rrna_16s=base.rrna_16s,
        rrna_23s=base.rrna_23s,
        trna_types=trnas,
        source=",".join(sources),
    )


def collect_evidence(search_dir: str, names: Iterable[str]) -> Dict[str, RnaEvidence]:
    """:func:`evidence_for_genome` for several genomes, omitting the ones with no source.

    Omission rather than an empty :class:`RnaEvidence` per genome, so a caller can tell "no
    annotation run exists" from "the annotation ran and found nothing" by whether the key is there.
    """
    if not os.path.isdir(search_dir):
        return {}
    out: Dict[str, RnaEvidence] = {}
    for name in names:
        found = evidence_for_genome(search_dir, name)
        if found is not None:
            out[name] = found
    return out


# --- The stats table -----------------------------------------------------------------------------


def _number(raw: str) -> Optional[float]:
    """A stats cell as a float, or None if it does not hold one. Empty is absent, not zero."""
    text = raw.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def annotate_stats(
    stats_path: str, evidence: Optional[Dict[str, RnaEvidence]] = None
) -> List[Assessment]:
    """Append :data:`STATS_COLUMNS` to a ``.stats.tsv``; return the assessment of every row.

    Appended and never reordered, and the base columns are not touched, because
    :data:`metawrap2.checkm.STATS_HEADER` is a published contract read by position - a golden test
    pins it byte for byte.

    A genome with no entry in *evidence* gets the unknown rRNA state rather than a failed table:
    tagging is worth having on a run that never called ``annotate``, which is most runs. A row whose
    completeness or contamination will not parse is tagged :attr:`MimagTag.UNSCORED` for the same
    reason - a table that raises here would cost a finished run its whole stats file over a cell.

    Re-running is safe: the columns are reused when already present, which a ``--resume`` of a step
    that annotates twice would otherwise turn into six columns and a table nothing can parse.
    """
    if not os.path.isfile(stats_path):
        return []
    with open(stats_path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        return []

    header = rows[0]
    at = {name: i for i, name in enumerate(header)}
    if STATS_COLUMNS[0] in at:
        positions = [at[name] for name in STATS_COLUMNS]
    else:
        positions = [len(header) + i for i in range(len(STATS_COLUMNS))]
        header = header + list(STATS_COLUMNS)
    completeness_at = at.get("completeness", 1)
    contamination_at = at.get("contamination", 2)

    width = len(header)
    out: List[List[str]] = [header]
    assessments: List[Assessment] = []
    for row in rows[1:]:
        row = row + [""] * (width - len(row))
        found = (evidence or {}).get(row[0])
        verdict = assess(
            row[0],
            _number(row[completeness_at]) if completeness_at < len(row) else None,
            _number(row[contamination_at]) if contamination_at < len(row) else None,
            found,
        )
        tag_at, rrna_at, trna_at = positions
        row[tag_at] = verdict.tag.value
        row[rrna_at] = verdict.evidence.describe_rrna()
        row[trna_at] = (
            "" if verdict.evidence.trna_types is None else str(verdict.evidence.trna_types)
        )
        out.append(row)
        assessments.append(verdict)
    with open(stats_path, "w") as fh:
        fh.writelines("\t".join(row) + "\n" for row in out)
    return assessments


#: The tiers a summary reports, in the order they are reported: best first, then the two states that
#: are not tiers. ``UNSCORED`` is last because it is a fact about the run rather than about a bin.
_SUMMARY_ORDER = (
    (MimagTag.HIGH, "high quality"),
    (MimagTag.HIGH_UNVERIFIED, "high quality with rRNA/tRNA unmeasured"),
    (MimagTag.MEDIUM, "medium quality"),
    (MimagTag.LOW, "low quality"),
    (MimagTag.NOT_MIMAG, "too contaminated for any MIMAG tier"),
    (MimagTag.UNSCORED, "unscored"),
)


def tally(assessments: Sequence[Assessment]) -> Dict[MimagTag, int]:
    """How many genomes landed in each tag. Every tag is present, at zero if nothing is in it."""
    counts = {tag: 0 for tag in MimagTag}
    for verdict in assessments:
        counts[verdict.tag] += 1
    return counts


def summarise(assessments: Sequence[Assessment]) -> str:
    """One line for the end of a run: how many genomes in each MIMAG tier.

    Tiers with no genomes in them are left out, so the line stays short on a small bin set, except
    that high quality and its unverified form are both shown whenever either is non-zero - the point
    of the pair is the comparison between them, and printing only the one that happens to be
    occupied hides it.
    """
    if not assessments:
        return "no bins to tag"
    counts = tally(assessments)
    show_high = counts[MimagTag.HIGH] or counts[MimagTag.HIGH_UNVERIFIED]
    parts = [
        "%d %s" % (counts[tag], label)
        for tag, label in _SUMMARY_ORDER
        if counts[tag] or (show_high and tag in (MimagTag.HIGH, MimagTag.HIGH_UNVERIFIED))
    ]
    return "MIMAG tiers over %d bin(s): %s" % (len(assessments), ", ".join(parts))
