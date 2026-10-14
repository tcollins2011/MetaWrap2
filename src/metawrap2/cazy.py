"""CAZyme families from dbCAN, and the agreement between the three ways it calls them.

Carbohydrate-active enzymes are the part of a genome's function that a KO table answers worst. KEGG
places a gene in a pathway; CAZy says which glycosidic bond it acts on, which is the question in any
study about polysaccharide degradation - the gut, soil carbon, plant biomass. The two annotations
overlap hardly at all for the enzymes that matter, so this is not a second opinion on the same
calls.

dbCAN runs three independent methods - a DIAMOND search of CAZy, the dbCAN family HMMs, and the
dbCAN-sub subfamily HMMs - and writes what each one said per gene. That disagreement is the useful
part of its output and is the thing a single "CAZy family" column throws away: on one real bin, 65
of 1,772 proteins got a CAZyme call and a third of them were supported by one method only. dbCAN's
own guidance is to trust a call supported by at least two, so the number of agreeing methods travels
with every row here and the default for a confident call is two.

One more distinction the output keeps. A **CBM** is a carbohydrate-binding module, not a catalytic
domain: a protein with a CBM binds a polysaccharide, and may or may not cut it. Counting CBMs
alongside glycoside hydrolases inflates a count of degradative enzymes with proteins that do not
degrade anything, so the class summary separates them and says so.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

#: dbCAN's own recommendation: a family call supported by at least two of its three methods. Not a
#: threshold invented here - a one-method call is kept, reported and counted separately rather than
#: dropped, because which method found it is itself informative (DIAMOND finds close homologues of
#: characterised enzymes; the HMMs find divergent ones).
MIN_TOOLS = 2

#: The CAZy class of a family, by prefix. Order matters: "GT" and "GH" share no prefix but "CBM"
#: must be checked before "CE" would ever match, so these are matched whole rather than by
#: startswith over a dict in arbitrary order.
CAZY_CLASSES: Tuple[Tuple[str, str], ...] = (
    ("GH", "glycoside hydrolase"),
    ("GT", "glycosyltransferase"),
    ("PL", "polysaccharide lyase"),
    ("CE", "carbohydrate esterase"),
    ("AA", "auxiliary activity (redox)"),
    ("CBM", "carbohydrate-binding module"),
)

#: Classes whose members catalyse a reaction on a carbohydrate. CBM is deliberately absent: it binds
#: one. Summing it into "degradative capacity" is the error this constant exists to prevent.
CATALYTIC_CLASSES = ("GH", "GT", "PL", "CE", "AA")

#: dbCAN's columns in overview.tsv, as of 5.2. Named so a schema change is a KeyError at parse time
#: rather than an empty column downstream.
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

#: ``GH13(45-389)`` -> ``GH13``; ``GH13_e82`` keeps its subfamily suffix, which is the point of
#: dbCAN-sub. Coordinates are dropped because they are per-method and the family is not.
_COORDINATES = re.compile(r"\((\d+)-(\d+)\)")


def _families(field: Optional[str]) -> Tuple[str, ...]:
    """Split one method's call into families, without its coordinates.

    A multi-domain protein gets several families from one method, separated by ``+``, and each may
    carry coordinates. Both are stripped here so that agreement between methods can be compared on
    the family rather than on the string.
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

    Checked longest-prefix first so ``CBM`` is never read as ``CE``-something, and returning None
    rather than guessing: dbCAN's DIAMOND column can hold a CAZy accession that is not a family at
    all, and a wrong class is worse than an unclassified one in a summary table.
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
    #: Families at least `MIN_TOOLS` methods proposed, which is dbCAN's own bar for a call.
    families_agreed: Tuple[str, ...]
    #: dbCAN's own pick, which can be empty when its methods disagree. Kept as it is written rather
    #: than recomputed, so a row can be checked against dbCAN's output directly.
    recommended: str
    diamond: Tuple[str, ...]
    hmm: Tuple[str, ...]
    subfamily: Tuple[str, ...]
    #: As dbCAN counts it. Recomputed here as well, and the two are asserted equal in the tests,
    #: because its number is the one a reader will compare against.
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

    Only genes with at least one call appear in that file, so its length is the CAZyme count and not
    the gene count - which is why the fraction of a genome that is CAZyme has to be computed against
    the gene table rather than from this file alone.
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
                    ec=_clean(row.get(OVERVIEW_COLUMNS["ec"])),
                    substrate=_clean_list(row.get(OVERVIEW_COLUMNS["substrate"])),
                )
            )
    return calls


def _clean(value: Optional[str]) -> str:
    text = (value or "").strip()
    return "" if text in ("", EMPTY) else text


def _clean_list(value: Optional[str]) -> str:
    """Clean a field that can be a list with empty members, e.g. ``-;peptidoglycan``.

    dbCAN joins a multi-domain protein's per-domain substrates, and writes ``-`` for a domain it has
    no mapping for - so the lone-dash check is not enough: real output carries ``-;peptidoglycan``
    for a protein with two domains, one of them unmapped. Left in, that dash becomes a substrate
    named "-" in any count.
    """
    text = _clean(value)
    if not text:
        return ""
    parts = [part.strip() for part in re.split(r"[;,]", text)]
    kept = [part for part in parts if part and part != EMPTY]
    return ";".join(dict.fromkeys(kept))


def _agreed(
    diamond: Sequence[str], hmm: Sequence[str], subfamily: Sequence[str], minimum: int = MIN_TOOLS
) -> Tuple[str, ...]:
    """Families proposed by at least *minimum* of the three methods.

    Compared at family level with the subfamily suffix removed, because dbCAN-sub exists to add that
    suffix: ``GH13`` from the HMMs and ``GH13_e82`` from dbCAN-sub are the same call, and treating
    them as different would make the subfamily method's every agreement look like a disagreement.
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
        # Stated per row because the distinction is lost the moment these are counted: a CBM-only
        # protein binds a polysaccharide and need not act on one.
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

    dbCAN predicts a substrate only for the families it has a mapping for, so the genes with no
    substrate are the majority and are not an error. Counted here rather than left as a column so
    that "what can this genome eat" has an answer at the genome level.
    """
    counts: Dict[str, int] = {}
    for call in calls:
        if not call.substrate:
            continue
        for substrate in (part.strip() for part in re.split(r"[;,]", call.substrate)):
            if substrate:
                counts[substrate] = counts.get(substrate, 0) + 1
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))
