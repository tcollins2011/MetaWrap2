"""Resistance and virulence calls, their genetic context, and how strong a claim each one supports.

A list of resistance genes is the easiest output in this pipeline to over-read. Two findings shape
everything here:

**Assemblies break exactly where the interesting things are.** Direct read mapping found 85 ARGs in
a dataset where individual assemblers found 20-44; only **12%** of ARG genomic contexts were captured
correctly; MEGAHIT's ARG-containing contigs averaged **284 bp**, so a 300 bp length filter would have
discarded 52% of them. Two rules follow, and they cost nothing: annotate before any length filter,
and record how much sequence sits between a gene and the end of its contig. A statement that there
is no mobile element within 10 kb of an ARG is unsupported when the ARG is 2 kb from a contig end -
the absence is of sequence, not of mobile elements. That is why every call here carries
:attr:`ArgCall.distance_to_contig_end` and why the context verdict has a ``context_unknown`` state
rather than defaulting to "no MGE".

**There is no field-wide threshold for "on a mobile element".** Published practice runs from 1 kb
(insertion-sequence adjacency) through 10 kb flanking windows to mere same-contig co-location, and
one paper describes its own 1 kb choice as "somehow arbitrary". So the window is a parameter, its
value goes into provenance, and the distance itself is reported rather than only a boolean.

The evidence ladder is Cao et al., *J Bacteriol* 2026 (10.1128/jb.00090-26). It exists because
"resistance gene detected" and "confers resistance" are four levels apart, and metagenomics reaches
the first two:

===========  =====================================================  =============================================
level        evidence                                               language it licenses
===========  =====================================================  =============================================
L1 presence  the gene is detected, with controls                    "resistance determinants detected"
L2 context   the gene is linked to a host and a neighbourhood        "resistance potential"
L3 activity  expression or protein evidence                         "active under the conditions tested"
L4 phenotype isolate AST, functional metagenomics, genetics         "confers resistance"
===========  =====================================================  =============================================

Every row carries its level, and :data:`LEVEL_LANGUAGE` is in the output so a reader does not have to
remember which phrasing a level licenses. Nothing this module produces can reach L3 or L4.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: Default flanking window for "is there a mobile element beside this gene", in bases. 10 kb is the
#: most common choice in the literature rather than a measured optimum, which is exactly why it is a
#: parameter recorded in provenance instead of a constant buried in a function.
MOBILITY_WINDOW = 10_000


class EvidenceLevel(Enum):
    """Cao et al.'s four levels. This module assigns L1 or L2 and can never assign more."""

    PRESENCE = "L1_presence"
    CONTEXT = "L2_context"


#: What each level licenses a reader to say. Emitted in the table so the constraint travels with the
#: data rather than living only in documentation nobody re-reads.
LEVEL_LANGUAGE: Mapping[str, str] = {
    EvidenceLevel.PRESENCE.value: "resistance determinants detected in the sample",
    EvidenceLevel.CONTEXT.value: "resistance potential, in this genomic context",
}


class Context(Enum):
    """Whether a mobile element was found beside a gene, or whether the question is answerable.

    The third state is the point. A gene closer to the end of its contig than the mobility window is
    one whose neighbourhood was not observed, and reporting that as "no mobile element nearby" turns
    missing sequence into a negative finding.
    """

    MGE_NEARBY = "mge_nearby"
    NO_MGE_IN_WINDOW = "no_mge_in_window"
    UNKNOWN = "context_unknown"


@dataclass(frozen=True)
class ArgCall:
    """One resistance, stress or virulence call, with where it is and what is next to it."""

    genome: str
    contig: str
    start: int
    end: int
    strand: str
    #: The gene symbol, e.g. ``aph(3')-IIIa``.
    symbol: str
    name: str
    #: AMRFinderPlus' Type: AMR, STRESS or VIRULENCE. Kept distinct because they are not
    #: interchangeable units of risk and must not be summed into one score.
    element_type: str
    subtype: str
    drug_class: str
    subclass: str
    #: AMRFinderPlus' Method: EXACTX, BLASTX, PARTIALX, HMM and so on. A PARTIALX call is a fragment
    #: of a reference gene, which is a weaker observation than EXACTX and is usually a contig end.
    method: str
    identity: Optional[float]
    coverage: Optional[float]
    contig_length: int
    distance_to_contig_end: int
    tool: str
    database_version: str
    #: The ARO accession, when the backend supplies one. RGI does - its calls are made against
    #: CARD, so the ontology term is the identity of the hit rather than something inferred
    #: afterwards. AMRFinderPlus does not, and argNorm fills that in; an empty string here means
    #: "this backend did not say", not "no term exists".
    aro: str = ""

    @property
    def is_partial_method(self) -> bool:
        return "PARTIAL" in self.method.upper()

    @property
    def is_virulence(self) -> bool:
        return self.element_type.upper() == "VIRULENCE"


def _as_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _as_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


#: AMRFinderPlus' column names, as of 4.x. Named here so a schema change is a KeyError at parse time
#: rather than silently empty columns in the output.
AMRFINDER_COLUMNS = {
    "contig": "Contig id",
    "start": "Start",
    "stop": "Stop",
    "strand": "Strand",
    "symbol": "Element symbol",
    "name": "Element name",
    "type": "Type",
    "subtype": "Subtype",
    "class": "Class",
    "subclass": "Subclass",
    "method": "Method",
    "identity": "% Identity to reference",
    "coverage": "% Coverage of reference",
}


def parse_amrfinder(
    path: str,
    genome: str,
    contig_lengths: Mapping[str, int],
    database_version: str,
    tool: str = "amrfinderplus",
) -> List[ArgCall]:
    """Read an AMRFinderPlus report into :class:`ArgCall` objects.

    *contig_lengths* is what turns a coordinate into a context claim: without it there is no way to
    know whether a gene 2 kb from a contig's start is 2 kb from the end of the observed sequence or
    in the middle of a megabase. AMRFinderPlus does not report contig lengths, so they come from the
    assembly the caller already has.
    """
    calls: List[ArgCall] = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames and AMRFINDER_COLUMNS["contig"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like an AMRFinderPlus report: no %r column"
                % (path, AMRFINDER_COLUMNS["contig"])
            )
        for row in reader:
            contig = row[AMRFINDER_COLUMNS["contig"]]
            start = _as_int(row[AMRFINDER_COLUMNS["start"]])
            end = _as_int(row[AMRFINDER_COLUMNS["stop"]])
            if start > end:
                start, end = end, start
            length = contig_lengths.get(contig, 0)
            # min() of the two ends: a neighbourhood is only as observed as its shorter side.
            distance = min(start - 1, length - end) if length else 0
            calls.append(
                ArgCall(
                    genome=genome,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=row.get(AMRFINDER_COLUMNS["strand"], "") or "",
                    symbol=row.get(AMRFINDER_COLUMNS["symbol"], "") or "",
                    name=row.get(AMRFINDER_COLUMNS["name"], "") or "",
                    element_type=row.get(AMRFINDER_COLUMNS["type"], "") or "",
                    subtype=row.get(AMRFINDER_COLUMNS["subtype"], "") or "",
                    drug_class=row.get(AMRFINDER_COLUMNS["class"], "") or "",
                    subclass=row.get(AMRFINDER_COLUMNS["subclass"], "") or "",
                    method=row.get(AMRFINDER_COLUMNS["method"], "") or "",
                    identity=_as_float(row.get(AMRFINDER_COLUMNS["identity"], "")),
                    coverage=_as_float(row.get(AMRFINDER_COLUMNS["coverage"], "")),
                    contig_length=length,
                    distance_to_contig_end=max(distance, 0),
                    tool=tool,
                    database_version=database_version,
                )
            )
    return calls


@dataclass(frozen=True)
class MobileElement:
    """A mobile-element feature on a contig, from whichever caller found it."""

    contig: str
    start: int
    end: int
    kind: str
    source: str


def classify_context(
    call: ArgCall,
    elements: Optional[Sequence[MobileElement]],
    window: int = MOBILITY_WINDOW,
    whole_contig: Optional[Mapping[str, str]] = None,
) -> Tuple[Context, Optional[int], Tuple[str, ...]]:
    """(verdict, distance to the nearest element, the kinds found) for one call.

    *elements* is ``None`` when no mobile-element caller was run, and an empty sequence when one ran
    and found nothing. Those are different claims and collapsing them is the mistake this whole
    module is arranged to avoid: an empty list means "looked, found none", while ``None`` means the
    question was never asked. Only the first can ever support ``no_mge_in_window``.

    The order of the checks is the rest of the substance. A mobile element inside the window is a
    positive finding regardless of contig ends. Only when nothing is found does the contig end
    matter, and then it decides between "nothing there" and "could not tell" - a gene with less than
    *window* bases of observed sequence beside it cannot support either answer.
    """
    if elements is None and not whole_contig:
        return Context.UNKNOWN, None, ()
    # A contig that *is* a plasmid or a free viral genome is the strongest context statement
    # available, and it is not a distance: the gene sits on the element rather than near it. So it
    # is reported at distance 0 with the kind named, rather than squeezed into the interval search
    # by inventing coordinates for a whole contig.
    if whole_contig and call.contig in whole_contig:
        return Context.MGE_NEARBY, 0, (whole_contig[call.contig],)
    same_contig = [e for e in (elements or ()) if e.contig == call.contig]
    nearest: Optional[int] = None
    kinds: List[str] = []
    for element in same_contig:
        if element.end < call.start:
            gap = call.start - element.end
        elif element.start > call.end:
            gap = element.start - call.end
        else:
            gap = 0
        if gap <= window:
            kinds.append(element.kind)
            nearest = gap if nearest is None else min(nearest, gap)
    if kinds:
        return Context.MGE_NEARBY, nearest, tuple(sorted(set(kinds)))
    if elements is None:
        # Whole-contig calls were supplied and this contig was not among them, which says nothing
        # about what sits beside the gene - only that the contig as a whole was not called mobile.
        return Context.UNKNOWN, None, ()
    if call.distance_to_contig_end < window:
        return Context.UNKNOWN, None, ()
    return Context.NO_MGE_IN_WINDOW, None, ()


def evidence_level(context: Context) -> EvidenceLevel:
    """L2 only when the neighbourhood was actually observed; L1 otherwise.

    A call whose context is unknown has not been linked to a genomic context, which is precisely
    what L2 requires, so it stays at L1 however confident the detection itself is.
    """
    if context in (Context.MGE_NEARBY, Context.NO_MGE_IN_WINDOW):
        return EvidenceLevel.CONTEXT
    return EvidenceLevel.PRESENCE


#: Columns of the AMR table. Long and explicit on purpose: Cao's L1 validity requirements are the
#: database name, version and date, the tool version, the identity and coverage thresholds, the
#: multi-mapping rule and the normalisation denominator. They are requirements, not metadata.
ARG_COLUMNS: Tuple[str, ...] = (
    "genome",
    "contig",
    "start",
    "end",
    "strand",
    "element_symbol",
    "element_name",
    "element_type",
    "subtype",
    "class",
    "subclass",
    "method",
    "percent_identity",
    "percent_coverage",
    "contig_length",
    "distance_to_contig_end",
    "aro",
    "context",
    "distance_to_nearest_mge",
    "mge_kinds",
    "evidence_level",
    "permitted_language",
    "mobility_window",
    "tool",
    "database_version",
)


def arg_row(
    call: ArgCall,
    context: Context,
    distance_to_mge: Optional[int],
    kinds: Sequence[str],
    window: int,
) -> Dict[str, object]:
    level = evidence_level(context)
    return {
        "genome": call.genome,
        "contig": call.contig,
        "start": call.start,
        "end": call.end,
        "strand": call.strand,
        "element_symbol": call.symbol,
        "element_name": call.name,
        "element_type": call.element_type,
        "subtype": call.subtype,
        "class": call.drug_class,
        "subclass": call.subclass,
        "method": call.method,
        "percent_identity": "" if call.identity is None else "%.2f" % call.identity,
        "percent_coverage": "" if call.coverage is None else "%.2f" % call.coverage,
        "contig_length": call.contig_length,
        "distance_to_contig_end": call.distance_to_contig_end,
        # Empty for a backend that does not supply one; argNorm fills those in
        # amr_normalised.tsv rather than this table guessing.
        "aro": call.aro,
        "context": context.value,
        "distance_to_nearest_mge": "" if distance_to_mge is None else distance_to_mge,
        "mge_kinds": ",".join(kinds),
        "evidence_level": level.value,
        "permitted_language": LEVEL_LANGUAGE[level.value],
        "mobility_window": window,
        "tool": call.tool,
        "database_version": call.database_version,
    }


def cluster_calls(calls: Sequence[ArgCall], max_gap: int = 10_000) -> List[List[ArgCall]]:
    """Group calls that sit together on one contig, nearest-neighbour within *max_gap*.

    Co-located resistance genes are the unit that actually moves: the three in this project's dev
    fixture - ant(6)-Ia, aph(3')-IIIa and sat4 within 2.3 kb - are a known aminoglycoside and
    streptothricin cassette, and reporting them as three unrelated findings loses the fact that they
    are one. Clustering is reported alongside the per-gene rows rather than instead of them.
    """
    by_contig: Dict[str, List[ArgCall]] = {}
    for call in calls:
        by_contig.setdefault(call.contig, []).append(call)
    clusters: List[List[ArgCall]] = []
    for contig_calls in by_contig.values():
        ordered = sorted(contig_calls, key=lambda c: (c.start, c.end))
        current = [ordered[0]]
        for call in ordered[1:]:
            if call.start - current[-1].end <= max_gap:
                current.append(call)
            else:
                clusters.append(current)
                current = [call]
        clusters.append(current)
    return sorted(clusters, key=lambda group: (-len(group), group[0].contig, group[0].start))


def summarise(calls: Sequence[ArgCall]) -> Dict[str, int]:
    """Counts worth printing at the end of a run, split by the things that must not be summed."""
    out = {
        "amr": 0,
        "stress": 0,
        "virulence": 0,
        "partial_method": 0,
        "within_window_of_contig_end": 0,
    }
    for call in calls:
        kind = call.element_type.upper()
        if kind == "AMR":
            out["amr"] += 1
        elif kind == "STRESS":
            out["stress"] += 1
        elif kind == "VIRULENCE":
            out["virulence"] += 1
        if call.is_partial_method:
            out["partial_method"] += 1
        if call.distance_to_contig_end < MOBILITY_WINDOW:
            out["within_window_of_contig_end"] += 1
    return out


# ── argNorm: one name per gene, and the drugs it implicates ──────────────────────────────────


@dataclass(frozen=True)
class NormalisedCall:
    """One resistance call mapped onto the Antibiotic Resistance Ontology.

    Worth having for two separate reasons. The gene symbols differ between databases, so an ARO
    accession is what makes a call found by two backends one finding rather than two; and the drugs
    an ARO term confers resistance to turn a gene list into a statement a reader can use - "cepA"
    says little, "confers resistance to cephalosporin" is the result.
    """

    genome: str
    gene: str
    #: Empty when argNorm could not place the gene. Absence from the ontology is not evidence that
    #: a call is wrong - a recently described gene may simply not be in it yet.
    aro: str
    aro_name: str
    #: The drugs the ARO term confers resistance to, by name. One gene can implicate many: APH(3')-IIIa
    #: covers nine aminoglycosides.
    drugs: Tuple[str, ...]


def parse_argnorm(path: str) -> List[NormalisedCall]:
    """Parse argNorm's output over a hAMRonization table.

    The leading ``# argNorm version: ...`` line is why this is not a bare DictReader: handed the
    file directly, csv takes that comment as the header and every column comes back under the key
    ``# argNorm version: 1.1.0``, which fails as an empty result rather than as an error.
    """
    with open(path) as handle:
        lines = [line for line in handle if not line.startswith("#")]
    out: List[NormalisedCall] = []
    for row in csv.DictReader(lines, delimiter="\t"):
        drugs = tuple(
            drug.strip()
            for drug in (row.get("confers_resistance_to_names") or "").split(",")
            if drug.strip()
        )
        out.append(
            NormalisedCall(
                genome=row.get("input_file_name", ""),
                gene=row.get("gene_symbol", ""),
                aro=(row.get("ARO") or "").strip(),
                aro_name=(row.get("ARO_name") or "").strip(),
                drugs=drugs,
            )
        )
    return out


def drugs_implicated(rows: Sequence[NormalisedCall]) -> List[Tuple[str, List[str]]]:
    """drug -> the genes implicating it, most-implicated drug first.

    Grouped by drug rather than by gene because that is the direction the question runs in: a
    reader wants to know what this sample is resistant to, and one gene covering nine
    aminoglycosides should not be read as nine findings.
    """
    by_drug: Dict[str, List[str]] = {}
    for row in rows:
        for drug in row.drugs:
            genes = by_drug.setdefault(drug, [])
            if row.gene not in genes:
                genes.append(row.gene)
    return sorted(by_drug.items(), key=lambda item: (-len(item[1]), item[0]))


# ── mobile-element context from geNomad ──────────────────────────────────────────────────────

#: geNomad's provirus table, which is the one product here with coordinates: a prophage integrated
#: at a stated position in a host contig. Everything else geNomad calls is a whole contig.
PROVIRUS_SUFFIX = "_provirus.tsv"
#: What a whole-contig call is recorded as in the `mge_kinds` column.
PLASMID_KIND = "plasmid_contig"
VIRAL_KIND = "viral_contig"
PROVIRUS_KIND = "provirus"
#: A provirus carrying an integrase is a different claim from one without: the integrase is the
#: machinery that moved it, so it is named separately rather than folded into `provirus`.
PROVIRUS_WITH_INTEGRASE = "provirus_with_integrase"


def read_proviruses(path: str, source: str = "genomad") -> List[MobileElement]:
    """Parse geNomad's ``*_provirus.tsv`` into mobile elements on their host contigs.

    The coordinates are on ``source_seq``, the host contig, which is what makes these usable for
    the distance question - unlike geNomad's other outputs, which classify whole contigs.
    """
    elements: List[MobileElement] = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig = row.get("source_seq") or ""
            if not contig:
                continue
            integrases = (row.get("integrases") or "").strip()
            kind = (
                PROVIRUS_WITH_INTEGRASE
                if integrases and integrases not in ("NA", "")
                else PROVIRUS_KIND
            )
            try:
                start, end = int(float(row["start"])), int(float(row["end"]))
            except (KeyError, TypeError, ValueError):
                continue
            elements.append(
                MobileElement(contig=contig, start=start, end=end, kind=kind, source=source)
            )
    return elements


def read_contig_calls(path: str, kind: str, column: str = "seq_name") -> Dict[str, str]:
    """contig -> *kind*, from a geNomad summary or the `viruses` module's table of one.

    Rows whose name contains ``|provirus_`` are skipped. Those are the *trimmed viral part* of a
    host contig, not the contig itself, so treating one as a whole-contig call would mark the host
    contig mobile under a name no ARG call ever carries - and the provirus table already covers
    that case with real coordinates.
    """
    out: Dict[str, str] = {}
    if not os.path.isfile(path):
        return out
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get(column) or "").strip()
            if not name or "|provirus_" in name:
                continue
            out[name] = kind
    return out


def mobile_elements_from_genomad(directory: str) -> Tuple[List[MobileElement], Dict[str, str]]:
    """(elements with coordinates, whole-contig calls) from a geNomad or `viruses` output tree.

    Accepts either, because both are things a user has: a geNomad run directory, or the output of
    `metawrap2 viruses`, which keeps geNomad's own output under ``genomad/<input>/`` and writes its
    own ``plasmid_contigs.tsv``. Searched rather than demanded by exact path - geNomad nests its
    results under ``<name>_summary/`` and ``<name>_find_proviruses/``, and asking a user for those
    paths is asking them to know the layout.

    Returns empty containers when it finds nothing, which the caller must distinguish from "not
    asked for": an empty result still means a caller looked, and that is what licenses
    ``no_mge_in_window``.
    """
    elements: List[MobileElement] = []
    whole: Dict[str, str] = {}
    for root, _dirs, files in os.walk(directory):
        for name in files:
            path = os.path.join(root, name)
            if name.endswith(PROVIRUS_SUFFIX):
                elements.extend(read_proviruses(path))
            elif name.endswith("_plasmid_summary.tsv") or name == "plasmid_contigs.tsv":
                whole.update(read_contig_calls(path, PLASMID_KIND))
            elif name.endswith("_virus_summary.tsv"):
                whole.update(read_contig_calls(path, VIRAL_KIND))
    return elements, whole


def context_summary(
    rows: Sequence[Tuple[Context, Optional[int], Sequence[str]]],
) -> Dict[str, int]:
    """Counts per context verdict and per element kind, for the end-of-run report."""
    out: Dict[str, int] = {context.value: 0 for context in Context}
    for context, _distance, kinds in rows:
        out[context.value] += 1
        for kind in kinds:
            out[kind] = out.get(kind, 0) + 1
    return out


# ── RGI, the second backend ──────────────────────────────────────────────────────────────────

#: RGI's columns, named so a schema change is a KeyError at parse time rather than empty output.
RGI_COLUMNS = {
    "orf": "ORF_ID",
    "contig": "Contig",
    "start": "Start",
    "end": "Stop",
    "strand": "Orientation",
    "cut_off": "Cut_Off",
    "aro_name": "Best_Hit_ARO",
    "aro": "ARO",
    "identity": "Best_Identities",
    "model": "Model_type",
    "drug_class": "Drug Class",
    "mechanism": "Resistance Mechanism",
    "family": "AMR Gene Family",
    "coverage": "Percentage Length of Reference Sequence",
    "nudged": "Nudged",
    "snps": "SNPs_in_Best_Hit_ARO",
}

#: RGI's confidence tiers, strongest first. `Perfect` is an exact match to a curated reference,
#: `Strict` a hit above the curated bit-score cutoff but not identical, `Loose` below it. RGI
#: excludes Loose unless asked; this records the tier per call because the tiers are not
#: interchangeable - on one real bin RGI's four non-Perfect calls sit at 27.8-46.1% identity.
RGI_TIERS: Tuple[str, ...] = ("Perfect", "Strict", "Loose")

#: What RGI writes in `Nudged` when it promoted a Loose hit to Strict because the gene looked
#: partial. Worth carrying: a nudged call is RGI overriding its own threshold.
NUDGED = "TRUE"


def parse_rgi(
    path: str,
    genome: str,
    contig_lengths: Optional[Mapping[str, int]] = None,
    database_version: str = "unknown",
) -> List[ArgCall]:
    """Read RGI's ``*.txt`` into the same calls :func:`parse_amrfinder` produces.

    One dataclass for both backends on purpose: the context classification, the evidence ladder and
    the clustering are properties of a *call on a contig*, not of the tool that made it, and a
    second type would mean a second copy of all three.

    Two RGI-specific details are kept rather than flattened. ``Cut_Off`` goes into ``method``,
    because Perfect and Strict are different claims - RGI's Strict hits on one real bin sit at
    27.8-46.1% identity to their reference, which is a homolog rather than the gene. And
    ``Model_type`` goes into ``subtype``: a protein *variant* model call means a point mutation
    against a wild-type sequence, which is a different kind of finding from an acquired gene.
    """
    lengths = contig_lengths or {}
    calls: List[ArgCall] = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and RGI_COLUMNS["aro_name"] not in reader.fieldnames:
            raise KeyError(
                "%s does not look like an RGI table: no %r column" % (path, RGI_COLUMNS["aro_name"])
            )
        for row in reader:
            contig = (row.get(RGI_COLUMNS["contig"]) or "").strip()
            symbol = (row.get(RGI_COLUMNS["aro_name"]) or "").strip()
            if not contig or not symbol:
                continue
            start = _as_int(row.get(RGI_COLUMNS["start"]) or "")
            end = _as_int(row.get(RGI_COLUMNS["end"]) or "")
            if start > end:
                start, end = end, start
            tier = (row.get(RGI_COLUMNS["cut_off"]) or "").strip() or "unknown"
            if str(row.get(RGI_COLUMNS["nudged"], "")).strip().upper() == NUDGED:
                # RGI promoted a Loose hit over its own cutoff because the gene looked partial. The
                # call is kept and the override is visible in the method rather than hidden in a
                # column nothing downstream reads.
                tier = "%s_nudged" % tier
            contig_length = lengths.get(contig, 0)
            distance = max(0, min(start - 1, contig_length - end)) if contig_length else 0
            calls.append(
                ArgCall(
                    genome=genome,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=(row.get(RGI_COLUMNS["strand"]) or "").strip() or "?",
                    symbol=symbol,
                    # The gene family rather than a repeat of the symbol: RGI's Best_Hit_ARO is
                    # already the name, and the family is the extra information.
                    name=(row.get(RGI_COLUMNS["family"]) or "").strip() or symbol,
                    element_type="AMR",
                    subtype=(row.get(RGI_COLUMNS["model"]) or "").strip() or "AMR",
                    drug_class=(row.get(RGI_COLUMNS["drug_class"]) or "").strip(),
                    subclass=(row.get(RGI_COLUMNS["mechanism"]) or "").strip(),
                    method=tier,
                    identity=_as_float(row.get(RGI_COLUMNS["identity"]) or ""),
                    coverage=_as_float(row.get(RGI_COLUMNS["coverage"]) or ""),
                    contig_length=contig_length,
                    distance_to_contig_end=distance,
                    tool="rgi",
                    database_version=database_version,
                    aro=_aro_accession(row.get(RGI_COLUMNS["aro"])),
                )
            )
    return calls


def _aro_accession(raw: object) -> str:
    """RGI's ARO column as a prefixed accession.

    RGI writes the bare number (``3002647``) where argNorm, hAMRonization and CARD's own site all
    write ``ARO:3002647``. Prefixing here means a row from either backend can be joined on this
    column without a per-tool rule at the join.
    """
    text = str(raw or "").strip()
    if not text:
        return ""
    return text if text.upper().startswith("ARO:") else "ARO:%s" % text


def card_version(card_json: str) -> str:
    """CARD's own version string, from the data directory's ``card.json``.

    Read rather than taken from a directory name: an ARG call is only as datable as its database,
    and a path can say anything. Returns "unknown" when the file cannot be read, which is reported
    rather than guessed.
    """
    import json

    path = card_json
    if os.path.isdir(path):
        path = os.path.join(path, "card.json")
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, ValueError):
        return "unknown"
    version = data.get("_version") if isinstance(data, dict) else None
    return str(version) if version else "unknown"


@dataclass(frozen=True)
class BackendAgreement:
    """How two resistance backends' calls on the same genomes overlap."""

    #: tool -> how many (genome, gene) pairs it called.
    per_tool: Mapping[str, int]
    #: (genome, gene) pairs every tool called.
    shared: int
    #: tool -> the genes only it called, sorted.
    unique: Mapping[str, Tuple[str, ...]]

    @property
    def tools(self) -> Tuple[str, ...]:
        return tuple(sorted(self.per_tool))


def backend_agreement(
    calls: Sequence[ArgCall],
    aro_for: Optional[Mapping[Tuple[str, str], str]] = None,
) -> BackendAgreement:
    """How two backends' calls on the same genomes overlap, on ARO terms where they are known.

    The join is per genome, on the ARO accession when the call has one and on the case-folded gene
    symbol otherwise. Symbols alone undercount agreement, measurably: on five CAMI2 bins
    AMRFinderPlus calls ``cepA`` and RGI calls ``CepA-29`` - the same beta-lactamase at the same
    99.67% identity, one gene by ARO and two by name. A symbol join reported one shared gene where
    two were shared.

    *aro_for* supplies terms for the backend that does not carry its own: RGI's calls are made
    against CARD and arrive with ``aro`` filled, while AMRFinderPlus' get theirs from argNorm, which
    runs after the tables are written. Keyed by (genome, lowercased symbol), which is how argNorm's
    output can be read back.
    """
    lookup = aro_for or {}
    by_tool: Dict[str, set] = {}
    for call in calls:
        symbol = call.symbol.lower()
        aro = call.aro or lookup.get((call.genome, symbol), "")
        by_tool.setdefault(call.tool, set()).add((call.genome, aro or symbol))
    shared: set = set.intersection(*by_tool.values()) if len(by_tool) > 1 else set()
    # Unique calls are reported by *name*, not by the key they were joined on: "only rgi:
    # ARO:3000777" is not something a reader can act on.
    label = {
        (
            call.genome,
            call.aro or lookup.get((call.genome, call.symbol.lower()), "") or call.symbol.lower(),
        ): call.symbol
        for call in calls
    }
    return BackendAgreement(
        per_tool={tool: len(found) for tool, found in by_tool.items()},
        shared=len(shared),
        unique={
            tool: tuple(sorted(label.get(key, key[1]) for key in (found - shared)))
            for tool, found in by_tool.items()
        },
    )
