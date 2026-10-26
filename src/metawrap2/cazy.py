"""CAZyme families from dbCAN, with the agreement between its three methods.

KEGG places a gene in a pathway; CAZy says which glycosidic bond it acts on, which matters for
studies of polysaccharide degradation (gut, soil carbon, plant biomass). The two annotations
overlap little for these enzymes.

dbCAN runs three methods (a DIAMOND search of CAZy, the dbCAN family HMMs and the dbCAN-sub
subfamily HMMs) and reports each one's call per gene. On one bin, 65 of 1,772 proteins got a
CAZyme call and a third of those were supported by one method only. dbCAN recommends trusting a
call supported by at least two, so every row carries the number of agreeing methods and a
confident call needs two.

A CBM is a carbohydrate-binding module, not a catalytic domain: a protein with one binds a
polysaccharide and may or may not cut it. The class summary counts CBMs separately from the
catalytic classes.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

#: dbCAN's recommendation: a family call supported by at least two of its three methods. One-method
#: calls are kept and counted separately, because the method is informative (DIAMOND finds close
#: homologues of characterised enzymes; the HMMs find divergent ones).
MIN_TOOLS = 2

#: The CAZy class of a family, by prefix, in CAZy's conventional order. :func:`family_class`
#: matches the longest prefix first.
CAZY_CLASSES: Tuple[Tuple[str, str], ...] = (
    ("GH", "glycoside hydrolase"),
    ("GT", "glycosyltransferase"),
    ("PL", "polysaccharide lyase"),
    ("CE", "carbohydrate esterase"),
    ("AA", "auxiliary activity (redox)"),
    ("CBM", "carbohydrate-binding module"),
)

#: Classes whose members catalyse a reaction on a carbohydrate. CBM is excluded because it only
#: binds one.
CATALYTIC_CLASSES = ("GH", "GT", "PL", "CE", "AA")

#: dbCAN's columns in overview.tsv, as of 5.2. A schema change raises KeyError at parse time.
OVERVIEW_COLUMNS = {
    "gene": "Gene ID",
    "ec": "EC#",
    "hmm": "dbCAN_hmm",
    "sub": "dbCAN_sub",
    "diamond": "DIAMOND",
    "tools": "#ofTools",
    "recommended": "Recommend Results",
    "substrate": "Substrate",
}

#: dbCAN writes a lone "-" for every field it has nothing for.
EMPTY = "-"

#: ``GH13(45-389)`` -> ``GH13``; ``GH13_e82`` keeps its dbCAN-sub subfamily suffix. Coordinates
#: are dropped because they differ between methods.
_COORDINATES = re.compile(r"\((\d+)-(\d+)\)")


def _families(field: Optional[str]) -> Tuple[str, ...]:
    """Split one method's call into families, without its coordinates.

    A multi-domain protein gets several families from one method, separated by ``+``, and each may
    carry coordinates. Both are stripped so methods can be compared by family.
    """
    if not field or field.strip() in ("", EMPTY):
        return ()
    out: List[str] = []
    for piece in re.split(r"[+|]", field.strip()):
        name = _COORDINATES.sub("", piece).strip()
        if name and name != EMPTY and name not in out:
            out.append(name)
    return tuple(out)


def family_class(family: str) -> Optional[str]:
    """``GH13_e82`` -> ``GH``. None for anything that is not a CAZy family name.

    Checked longest prefix first. dbCAN's DIAMOND column can hold a CAZy accession that is not a
    family; that returns None.
    """
    for prefix, _description in sorted(CAZY_CLASSES, key=lambda item: -len(item[0])):
        if family.startswith(prefix) and (
            len(family) == len(prefix) or family[len(prefix)].isdigit()
        ):
            return prefix
    return None


@dataclass(frozen=True)
class CazymeCall:
    """One gene's CAZyme annotation, with which of dbCAN's three methods found it."""

    genome: str
    gene: str
    #: Every family any method proposed, deduplicated. The consensus is in `families_agreed`.
    families: Tuple[str, ...]
    #: Families proposed by at least `MIN_TOOLS` methods (dbCAN's bar for a call).
    families_agreed: Tuple[str, ...]
    #: dbCAN's pick, empty when its methods disagree. Kept as written so a row can be checked
    #: against dbCAN's output.
    recommended: str
    diamond: Tuple[str, ...]
    hmm: Tuple[str, ...]
    subfamily: Tuple[str, ...]
    #: dbCAN's count of agreeing methods. The tests check it against a recomputed count.
    tools: int
    ec: str
    substrate: str

    @property
    def is_confident(self) -> bool:
        return self.tools >= MIN_TOOLS

    @property
    def classes(self) -> Tuple[str, ...]:
        found: List[str] = []
        for family in self.families_agreed or self.families:
            klass = family_class(family)
            if klass and klass not in found:
                found.append(klass)
        return tuple(found)

    @property
    def is_catalytic(self) -> bool:
        """Whether any agreed family is an enzyme class rather than a binding module."""
        return any(klass in CATALYTIC_CLASSES for klass in self.classes)


def parse_overview(path: str, genome: str) -> List[CazymeCall]:
    """Read dbCAN's ``overview.tsv``.

    Only genes with at least one call appear in that file, so the CAZyme fraction of a genome has
    to be computed against the gene table.
    """
    calls: List[CazymeCall] = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and OVERVIEW_COLUMNS["gene"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like a dbCAN overview: no %r column"
                % (path, OVERVIEW_COLUMNS["gene"])
            )
        for row in reader:
            diamond = _families(row.get(OVERVIEW_COLUMNS["diamond"]))
            hmm = _families(row.get(OVERVIEW_COLUMNS["hmm"]))
            subfamily = _families(row.get(OVERVIEW_COLUMNS["sub"]))
            families: List[str] = []
            for group in (diamond, hmm, subfamily):
                for family in group:
                    if family not in families:
                        families.append(family)
            agreed = _agreed(diamond, hmm, subfamily)
            try:
                tools = int(row.get(OVERVIEW_COLUMNS["tools"], 0) or 0)
            except ValueError:
                tools = 0
            calls.append(
                CazymeCall(
                    genome=genome,
                    gene=row[OVERVIEW_COLUMNS["gene"]],
                    families=tuple(families),
                    families_agreed=agreed,
                    recommended=_clean(row.get(OVERVIEW_COLUMNS["recommended"])),
                    diamond=diamond,
                    hmm=hmm,
                    subfamily=subfamily,
                    tools=tools,
                    ec=_clean_list(row.get(OVERVIEW_COLUMNS["ec"])),
                    substrate=_clean_list(row.get(OVERVIEW_COLUMNS["substrate"])),
                )
            )
    return calls


def _clean(value: Optional[str]) -> str:
    text = (value or "").strip()
    return "" if text in ("", EMPTY) else text


def _clean_list(value: Optional[str]) -> str:
    """Clean a field that can be a list with empty members, e.g. ``-;peptidoglycan`` or ``-|-``.

    dbCAN joins a multi-domain protein's per-domain annotations and writes ``-`` for an unmapped
    domain, e.g. ``-;peptidoglycan``. Those members are dropped so "-" is never counted as a
    substrate.

    ``|`` is also a separator: dbCAN writes ``-|-`` in the EC column for a two-domain gene with no
    EC number (2 of 11 rows in the vendored fixture). Whole members are dropped, not dashes, so a
    partial EC number such as ``3.2.1.-`` is kept.
    """
    text = _clean(value)
    if not text:
        return ""
    parts = [part.strip() for part in re.split(r"[;,|]", text)]
    kept = [part for part in parts if part and part != EMPTY]
    return ";".join(dict.fromkeys(kept))


def _agreed(
    diamond: Sequence[str], hmm: Sequence[str], subfamily: Sequence[str], minimum: int = MIN_TOOLS
) -> Tuple[str, ...]:
    """Families proposed by at least *minimum* of the three methods.

    Compared at family level with the subfamily suffix removed: ``GH13`` from the HMMs and
    ``GH13_e82`` from dbCAN-sub count as the same call.
    """
    votes: Dict[str, int] = {}
    for group in (diamond, hmm, subfamily):
        for family in {_base(name) for name in group}:
            votes[family] = votes.get(family, 0) + 1
    return tuple(sorted(family for family, count in votes.items() if count >= minimum))


def _base(family: str) -> str:
    """``GH13_e82`` -> ``GH13``; ``GH13`` unchanged."""
    return family.split("_", 1)[0]


#: Columns of the per-gene CAZyme table.
CAZYME_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "cazy_families",
    "cazy_families_agreed",
    "cazy_class",
    "dbcan_recommended",
    "methods_agreeing",
    "confident",
    "catalytic",
    "ec",
    "substrate",
    "diamond",
    "dbcan_hmm",
    "dbcan_sub",
)


def cazyme_row(call: CazymeCall) -> Dict[str, str]:
    return {
        "genome": call.genome,
        "gene": call.gene,
        "cazy_families": ";".join(call.families),
        "cazy_families_agreed": ";".join(call.families_agreed),
        "cazy_class": ";".join(call.classes),
        "dbcan_recommended": call.recommended,
        "methods_agreeing": str(call.tools),
        "confident": "yes" if call.is_confident else "no",
        # Per row, so CBM-only proteins (binding, not necessarily catalytic) can be told apart.
        "catalytic": "yes" if call.is_catalytic else "no",
        "ec": call.ec,
        "substrate": call.substrate,
        "diamond": ";".join(call.diamond),
        "dbcan_hmm": ";".join(call.hmm),
        "dbcan_sub": ";".join(call.subfamily),
    }


@dataclass(frozen=True)
class ClassSummary:
    """One CAZy class in one genome: how many genes, and how confident."""

    genome: str
    cazy_class: str
    description: str
    genes: int
    confident_genes: int
    families: int


def summarise_classes(calls: Sequence[CazymeCall], genome: str) -> List[ClassSummary]:
    """Per-class counts for one genome, in CAZy's conventional class order."""
    descriptions = dict(CAZY_CLASSES)
    genes: Dict[str, int] = {}
    confident: Dict[str, int] = {}
    families: Dict[str, set] = {}
    for call in calls:
        for klass in call.classes:
            genes[klass] = genes.get(klass, 0) + 1
            if call.is_confident:
                confident[klass] = confident.get(klass, 0) + 1
            families.setdefault(klass, set()).update(
                family
                for family in (call.families_agreed or call.families)
                if family_class(family) == klass
            )
    return [
        ClassSummary(
            genome=genome,
            cazy_class=klass,
            description=descriptions[klass],
            genes=genes[klass],
            confident_genes=confident.get(klass, 0),
            families=len(families.get(klass, ())),
        )
        for klass, _description in CAZY_CLASSES
        if klass in genes
    ]


#: Columns of the per-class summary.
CLASS_COLUMNS: Tuple[str, ...] = (
    "genome",
    "cazy_class",
    "description",
    "genes",
    "genes_with_two_or_more_methods",
    "distinct_families",
)


def class_row(summary: ClassSummary) -> Dict[str, str]:
    return {
        "genome": summary.genome,
        "cazy_class": summary.cazy_class,
        "description": summary.description,
        "genes": str(summary.genes),
        "genes_with_two_or_more_methods": str(summary.confident_genes),
        "distinct_families": str(summary.families),
    }


def substrate_counts(calls: Sequence[CazymeCall]) -> List[Tuple[str, int]]:
    """substrate -> how many genes dbCAN assigned to it, commonest first.

    dbCAN predicts a substrate only for families it has a mapping for, so most genes have none.
    """
    counts: Dict[str, int] = {}
    for call in calls:
        if not call.substrate:
            continue
        for substrate in (part.strip() for part in re.split(r"[;,]", call.substrate)):
            if substrate:
                counts[substrate] = counts.get(substrate, 0) + 1
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))
