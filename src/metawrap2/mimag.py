"""MIMAG quality tags for recovered genomes (Bowers et al., Nat Biotechnol 2017;35:725-731).

CheckM and CheckM2 give completeness and contamination. MIMAG high quality also requires 5S, 16S
and 23S rRNA and tRNAs for at least 18 of the 20 amino acids, which CheckM does not report. For
this reason :data:`metawrap2.benchmark.scoring.QUALITY_TIERS` calls its 90/5 tier
``bins_near_complete``.

``annotate --bakta`` runs tRNAscan-SE, Aragorn and Infernal/Rfam over every genome, so where Bakta
has run, this module parses its table. Where it has not, there is no rRNA/tRNA evidence.

Missing evidence and a measured absence are tagged differently. A genome that passes 90/5 with no
evidence is tagged :data:`MimagTag.HIGH_UNVERIFIED` with ``mimag_rrna=unknown``. A genome where a
search found genes missing is tagged :data:`MimagTag.MEDIUM`, with the genes found listed.

Bins are not filtered, as with GUNC (see :mod:`metawrap2.gunc`); the tag is a column.
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
#: The inequalities follow Bowers et al. Table 1: high quality is ">90% completion" and "<5%
#: contamination", medium is ">=50% completion". So exactly 90.0 is not high quality and exactly
#: 50.0 is medium.
#:
#: :func:`metawrap2.benchmark.scoring.quality_counts` uses ``>= 90`` and ``<= 5`` instead. The two
#: differ only for bins scored exactly 90.00 or 5.00; this tag follows the paper because it is
#: meant for submissions.
HIGH_MIN_COMPLETENESS = 90.0
HIGH_MAX_CONTAMINATION = 5.0

#: MIMAG medium quality: >=50% complete, <10% contaminated. Low quality is <50% complete with the
#: same contamination ceiling.
MEDIUM_MIN_COMPLETENESS = 50.0

#: Contamination at or above this fits no MIMAG tier: Bowers et al. require <10% contamination for
#: every tier, including low quality.
MAX_CONTAMINATION = 10.0

#: tRNAs for at least 18 of the 20 amino acids. Counted as distinct amino acids, not tRNA features
#: (40 tRNA genes covering 15 amino acids do not pass).
MIN_TRNA_TYPES = 18

#: The three rRNA genes MIMAG high quality requires, in the order they are reported.
RRNA_GENES = ("5S", "16S", "23S")

#: The 20 amino acids MIN_TRNA_TYPES counts against. Selenocysteine and pyrrolysine (the 21st and
#: 22nd) are excluded; tRNAscan-SE often reports SeC.
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

#: Columns this module appends to a ``.stats.tsv``, in order, after the existing columns. The base
#: columns are read by position (see :data:`metawrap2.checkm.STATS_HEADER`). Defined in
#: :mod:`metawrap2.constants` and re-exported here.
STATS_COLUMNS = MIMAG_STATS_COLUMNS

#: Written in ``mimag_rrna`` when no rRNA source was found, as distinct from ``none`` (searched,
#: none found).
UNKNOWN = "unknown"

#: Printed when a bin set has no rRNA evidence, which is usual in ``bin_refinement`` because it
#: runs before ``annotate``.
NO_EVIDENCE_HINT = (
    "No rRNA/tRNA evidence was found, so bins that clear 90/5 are tagged %s rather than "
    "high quality. Run 'metawrap2 annotate --bakta' over the final bins and re-run this step to "
    "resolve them." % "high_quality_unverified_rrna"
)

#: Bakta names a tRNA by amino acid in its product column: ``tRNA-Ala(tgc)``. The trailing digit in
#: ``tRNA-Ile2`` is an isotype, not a different amino acid, so the capture stops at three letters.
_BAKTA_TRNA = re.compile(r"tRNA-([A-Za-z]{3})")

#: Bakta's rRNA products, e.g. ``16S ribosomal RNA``. Matched on the subunit, since the rest of the
#: string varies ("16S ribosomal RNA (partial)").
_BAKTA_RRNA = re.compile(r"\b(5S|16S|23S)\b")

#: Rfam family names for the three subunits, as in an Infernal ``--tblout`` target column. Matched
#: exactly, so the eukaryotic ``5_8S_rRNA`` family does not count as 5S.
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
    """A genome's MIMAG tier, plus states MIMAG does not define.

    A ``str`` enum, so values can be written to and compared with TSV cells directly.
    """

    #: >90 / <5 with 5S, 16S, 23S and >=18 tRNA types all measured present.
    HIGH = "high_quality"
    #: >90 / <5, with rRNA/tRNA not measured. Not a MIMAG tier.
    HIGH_UNVERIFIED = "high_quality_unverified_rrna"
    #: >=50 / <10. Also used for a genome passing 90/5 with a measured rRNA or tRNA absence.
    MEDIUM = "medium_quality"
    #: <50 / <10.
    LOW = "low_quality"
    #: Contamination at or above 10%: no MIMAG tier applies. See :data:`MAX_CONTAMINATION`.
    NOT_MIMAG = "not_mimag"
    #: No completeness or contamination available.
    UNSCORED = "unscored"


@dataclass(frozen=True)
class RnaEvidence:
    """rRNA and tRNA genes found for a genome.

    ``None`` in any field means "not measured", not "absent". For example, Infernal output without
    tRNAscan-SE output sets the rRNA fields and leaves ``trna_types`` as ``None``.
    """

    rrna_5s: Optional[bool] = None
    rrna_16s: Optional[bool] = None
    rrna_23s: Optional[bool] = None
    trna_types: Optional[int] = None
    #: File(s) the evidence came from.
    source: str = ""

    @property
    def rrna(self) -> Tuple[Optional[bool], Optional[bool], Optional[bool]]:
        """The three rRNA flags in :data:`RRNA_GENES` order."""
        return (self.rrna_5s, self.rrna_16s, self.rrna_23s)

    def describe_rrna(self) -> str:
        """Return the ``mimag_rrna`` cell: ``5S+16S+23S``, ``16S``, ``none``, or ``unknown``.

        ``none`` means a search found nothing; ``unknown`` means no search.
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
    """One genome's MIMAG tag and the values it was based on."""

    bin: str
    completeness: Optional[float]
    contamination: Optional[float]
    evidence: RnaEvidence
    tag: MimagTag

    @property
    def rrna_measured(self) -> bool:
        """True if this genome was searched for all three rRNA genes."""
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
    """Return the MIMAG assessment for one genome.

    Checks, in order:

    1. no completeness or contamination -> :attr:`MimagTag.UNSCORED`.
    2. contamination >= 10% -> :attr:`MimagTag.NOT_MIMAG`.
    3. >90 / <5 -> high quality if rRNA/tRNA are measured present, medium if a measured absence
       fails a criterion, else :attr:`MimagTag.HIGH_UNVERIFIED`.
    4. >=50 -> medium, else low.
    """
    found = evidence if evidence is not None else RnaEvidence()
    if completeness is None or contamination is None:
        return Assessment(name, completeness, contamination, found, MimagTag.UNSCORED)
    if contamination >= MAX_CONTAMINATION:
        return Assessment(name, completeness, contamination, found, MimagTag.NOT_MIMAG)

    if completeness > HIGH_MIN_COMPLETENESS and contamination < HIGH_MAX_CONTAMINATION:
        # A criterion fails only on a measurement, so False and None are tested separately; using
        # ``not flag`` would demote unannotated genomes to medium.
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
    """Return how many of the 20 amino acids the tRNA labels cover.

    Labels are normalised to ``Xxx`` capitalisation, since tools and versions differ in case.
    Labels not in :data:`AMINO_ACIDS` (``Undet``, ``Pseudo``, ``SeC``, ``Sup``) are ignored.
    """
    covered = {
        label.strip().capitalize() for label in amino_acids if label and label.strip()
    } & AMINO_ACIDS
    return len(covered)


def parse_bakta_tsv(path: str) -> RnaEvidence:
    """Return rRNA and tRNA evidence from a Bakta ``<genome>.tsv``.

    The table has one row per feature with ``Type`` and ``Product`` columns. The header line
    (``#Sequence Id``) follows a variable number of ``#`` comment lines, so columns are found by
    name in the ``#`` lines.

    All four fields are set: Bakta searched for both, so a genome with no rRNA row has a measured
    absence.
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
    """Return the same evidence from a Bakta ``<genome>.gff3``, for runs without the TSV.

    Column 3 is the feature type and column 9 the attributes (``product=``/``Name=``). Used only
    as a fallback, since GFF attribute keys vary between Bakta versions.
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
    """Return how many amino acids a tRNAscan-SE table covers, or None if it names none.

    The default output has three header lines, then whitespace-separated columns; the fifth is the
    tRNA type (``Ala``, ``Undet``, ``Pseudo``). Columns are read by position because the header
    wraps over two lines.

    Returns None, not 0, when no row names an amino acid, since that more likely means a truncated
    or unexpected file than a genome without tRNAs.
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
    """Return rRNA presence from an Infernal ``cmscan --tblout`` table.

    Subunits are recognised by Rfam family name (see :data:`_RFAM_FAMILIES`). Only the rRNA fields
    are set; ``trna_types`` stays None.

    ``--tblout`` is space-padded, so lines are split on whitespace. The target name is the first
    field, or the second in the ``--fmt 2`` layout.
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


#: Fallback rRNA/tRNA file names, relative to the search directory, with ``%s`` for the genome
#: name. For a separate Infernal/tRNAscan-SE run; the pipeline does not write these files.
_TRNASCAN_NAMES = ("rna/%s.trnascan.tsv", "rna/%s.trnascan.txt", "%s.trnascan.tsv")
_INFERNAL_NAMES = ("rna/%s.infernal.tblout", "rna/%s.cmscan.tblout", "%s.infernal.tblout")


def evidence_for_genome(search_dir: str, name: str) -> Optional[RnaEvidence]:
    """Return the rRNA/tRNA evidence for genome *name* under *search_dir*, or None if there is none.

    *search_dir* is an ``annotate`` output directory. Bakta output is preferred because it covers
    both criteria; otherwise Infernal and tRNAscan-SE tables are merged, each setting only what it
    measured.
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
    """Run :func:`evidence_for_genome` for several genomes, omitting those with no source."""
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
    """Return a stats cell as a float, or None if empty or not a number."""
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

    The base columns are unchanged, since :data:`metawrap2.checkm.STATS_HEADER` is read by
    position (a golden test checks it).

    A genome with no entry in *evidence* gets the unknown rRNA state. A row whose completeness or
    contamination does not parse is tagged :attr:`MimagTag.UNSCORED`.

    Re-running reuses existing columns, so a ``--resume`` does not add them twice.
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


#: Tiers in summary order: best first, then the two states that are not tiers.
_SUMMARY_ORDER = (
    (MimagTag.HIGH, "high quality"),
    (MimagTag.HIGH_UNVERIFIED, "high quality with rRNA/tRNA unmeasured"),
    (MimagTag.MEDIUM, "medium quality"),
    (MimagTag.LOW, "low quality"),
    (MimagTag.NOT_MIMAG, "too contaminated for any MIMAG tier"),
    (MimagTag.UNSCORED, "unscored"),
)


def tally(assessments: Sequence[Assessment]) -> Dict[MimagTag, int]:
    """Return the number of genomes per tag, including tags with zero."""
    counts = {tag: 0 for tag in MimagTag}
    for verdict in assessments:
        counts[verdict.tag] += 1
    return counts


def summarise(assessments: Sequence[Assessment]) -> str:
    """Return a one-line count of genomes per MIMAG tier, for the end of a run.

    Empty tiers are omitted, except that high quality and its unverified form are both shown when
    either is non-zero, so they can be compared.
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
