"""Resistance and virulence calls, their genetic context, and the evidence level of each call.

Assemblies often break around resistance genes. In one dataset direct read mapping found 85 ARGs
where individual assemblers found 20-44, only 12% of ARG genomic contexts were assembled correctly,
and MEGAHIT's ARG-containing contigs averaged 284 bp (a 300 bp length filter would have removed 52%
of them). So genes are annotated before any length filter, and every call records
:attr:`ArgCall.distance_to_contig_end`. When a gene is closer to a contig end than the mobility
window and no mobile element is found, the verdict is ``context_unknown``, not "no MGE".

Published windows for "near a mobile element" range from 1 kb (insertion-sequence adjacency)
through 10 kb flanks to same-contig co-location. The window is therefore a parameter recorded in
provenance, and the distance itself is reported.

The evidence ladder is from Cao et al., *J Bacteriol* 2026 (10.1128/jb.00090-26). Metagenomics
reaches the first two levels:

===========  =====================================================  =============================================
level        evidence                                               language it licenses
===========  =====================================================  =============================================
L1 presence  the gene is detected, with controls                    "resistance determinants detected"
L2 context   the gene is linked to a host and a neighbourhood        "resistance potential"
L3 activity  expression or protein evidence                         "active under the conditions tested"
L4 phenotype isolate AST, functional metagenomics, genetics         "confers resistance"
===========  =====================================================  =============================================

Every row carries its level and the matching phrasing from :data:`LEVEL_LANGUAGE`. This module
assigns only L1 or L2.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: Default flanking window, in bases, for finding a mobile element beside a gene. 10 kb is the most
#: common choice in the literature, not a measured optimum, so it is a parameter recorded in
#: provenance.
MOBILITY_WINDOW = 10_000


class EvidenceLevel(Enum):
    """The two of Cao et al.'s four levels that this module can assign."""

    PRESENCE = "L1_presence"
    CONTEXT = "L2_context"


#: The phrasing each level supports, written into the table.
LEVEL_LANGUAGE: Mapping[str, str] = {
    EvidenceLevel.PRESENCE.value: "resistance determinants detected in the sample",
    EvidenceLevel.CONTEXT.value: "resistance potential, in this genomic context",
}


class Context(Enum):
    """Whether a mobile element was found beside a gene, or whether that cannot be determined.

    ``UNKNOWN`` covers genes closer to a contig end than the mobility window, whose neighbourhood
    was not fully observed.
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
    #: AMRFinderPlus' Type: AMR, STRESS or VIRULENCE. These are counted separately, never summed.
    element_type: str
    subtype: str
    drug_class: str
    subclass: str
    #: AMRFinderPlus' Method: EXACTX, BLASTX, PARTIALX, HMM and so on. PARTIALX is a fragment of a
    #: reference gene, usually at a contig end, and weaker evidence than EXACTX.
    method: str
    identity: Optional[float]
    coverage: Optional[float]
    contig_length: int
    distance_to_contig_end: int
    tool: str
    database_version: str
    #: The ARO accession, when the backend supplies one. RGI does, since it searches CARD;
    #: AMRFinderPlus does not, and argNorm fills it in later. Empty means the backend gave none.
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


#: AMRFinderPlus' column names, as of 4.x. A missing contig column raises KeyError at parse time.
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

    *contig_lengths* comes from the assembly, since AMRFinderPlus does not report contig lengths;
    it is needed to compute each call's distance to the contig end.
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
            # Distance to the nearer contig end.
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

    *elements* is ``None`` when no mobile-element caller was run, and empty when one ran and found
    nothing. Only the second can give ``no_mge_in_window``.

    A mobile element within the window gives ``mge_nearby`` regardless of contig ends. If none is
    found, a gene with less than *window* bases to the contig end gives ``context_unknown``.
    """
    if elements is None and not whole_contig:
        return Context.UNKNOWN, None, ()
    # A gene on a contig called plasmid or viral as a whole is on the element: distance 0, with
    # the kind named.
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
        # Only whole-contig calls were supplied and this contig was not one; nothing is known about
        # the gene's neighbourhood.
        return Context.UNKNOWN, None, ()
    if call.distance_to_contig_end < window:
        return Context.UNKNOWN, None, ()
    return Context.NO_MGE_IN_WINDOW, None, ()


def evidence_level(context: Context) -> EvidenceLevel:
    """L2 when the neighbourhood was observed; L1 when the context is unknown."""
    if context in (Context.MGE_NEARBY, Context.NO_MGE_IN_WINDOW):
        return EvidenceLevel.CONTEXT
    return EvidenceLevel.PRESENCE


#: Columns of the AMR table. Cao et al.'s L1 requirements include the database name, version and
#: date, the tool version, the identity and coverage thresholds, the multi-mapping rule and the
#: normalisation denominator.
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
        # Empty for a backend that does not supply one; argNorm fills these in amr_normalised.tsv.
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

    Co-located resistance genes often move together: ant(6)-Ia, aph(3')-IIIa and sat4 within 2.3 kb
    in the dev fixture are a known aminoglycoside and streptothricin cassette. Clusters are
    reported in addition to the per-gene rows.
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
    """End-of-run counts, by element type, partial calls and calls near a contig end."""
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

    Gene symbols differ between databases, so the ARO accession joins calls from two backends. The
    ARO term also gives the drugs the gene confers resistance to.
    """

    genome: str
    gene: str
    #: Empty when argNorm could not map the gene; a recently described gene may not be in the
    #: ontology yet.
    aro: str
    aro_name: str
    #: The drugs the ARO term confers resistance to, by name; APH(3')-IIIa covers nine
    #: aminoglycosides.
    drugs: Tuple[str, ...]


def parse_argnorm(path: str) -> List[NormalisedCall]:
    """Parse argNorm's output over a hAMRonization table.

    Lines starting with ``#`` (the leading ``# argNorm version: ...``) are dropped before
    DictReader, which would otherwise take that line as the header.
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

    Grouped by drug to show what the sample may be resistant to.
    """
    by_drug: Dict[str, List[str]] = {}
    for row in rows:
        for drug in row.drugs:
            genes = by_drug.setdefault(drug, [])
            if row.gene not in genes:
                genes.append(row.gene)
    return sorted(by_drug.items(), key=lambda item: (-len(item[1]), item[0]))


# ── mobile-element context from geNomad ──────────────────────────────────────────────────────

#: geNomad's provirus table, the only output with coordinates (a prophage's position in its host
#: contig). Other geNomad calls are whole contigs.
PROVIRUS_SUFFIX = "_provirus.tsv"
#: Labels used in the `mge_kinds` column.
PLASMID_KIND = "plasmid_contig"
VIRAL_KIND = "viral_contig"
PROVIRUS_KIND = "provirus"
#: A provirus with an integrase is labelled separately from one without.
PROVIRUS_WITH_INTEGRASE = "provirus_with_integrase"


def read_proviruses(path: str, source: str = "genomad") -> List[MobileElement]:
    """Parse geNomad's ``*_provirus.tsv`` into mobile elements on their host contigs.

    Coordinates are on ``source_seq``, the host contig.
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

    Rows whose name contains ``|provirus_`` are skipped: they are the viral part of a host contig,
    which the provirus table covers with coordinates.
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

    Accepts a geNomad run directory or `metawrap2 viruses` output (geNomad's output under
    ``genomad/<input>/`` plus ``plasmid_contigs.tsv``). Files are found by name anywhere in the
    tree, since geNomad nests them under ``<name>_summary/`` and ``<name>_find_proviruses/``.

    Returns empty containers when nothing is found. Callers treat that as "looked, found none",
    which differs from not running a caller.
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

#: RGI's column names. A missing Best_Hit_ARO column raises KeyError at parse time.
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

#: RGI's confidence tiers, strongest first: `Perfect` is an exact match to a curated reference,
#: `Strict` is above the curated bit-score cutoff, `Loose` below it (excluded by RGI unless asked).
#: The tier is recorded per call; on one bin the four non-Perfect calls had 27.8-46.1% identity.
RGI_TIERS: Tuple[str, ...] = ("Perfect", "Strict", "Loose")

#: RGI's `Nudged` value when it promoted a Loose hit to Strict because the gene looked partial.
NUDGED = "TRUE"


def parse_rgi(
    path: str,
    genome: str,
    contig_lengths: Optional[Mapping[str, int]] = None,
    database_version: str = "unknown",
) -> List[ArgCall]:
    """Read RGI's ``*.txt`` into the same calls :func:`parse_amrfinder` produces.

    Both backends share one dataclass, so context classification, evidence levels and clustering
    work the same for each.

    ``Cut_Off`` goes into ``method``, since Strict hits can be distant homologs (27.8-46.1%
    identity on one bin). ``Model_type`` goes into ``subtype``, since a protein variant model call
    is a point mutation, not an acquired gene.
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
                # RGI promoted a Loose hit because the gene looked partial; mark it in the method.
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
                    # The gene family; Best_Hit_ARO is already used as the symbol.
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

    RGI writes the bare number (``3002647``); argNorm, hAMRonization and CARD write
    ``ARO:3002647``. Prefixing here lets both backends join on this column.
    """
    text = str(raw or "").strip()
    if not text:
        return ""
    return text if text.upper().startswith("ARO:") else "ARO:%s" % text


def card_version(card_json: str) -> str:
    """CARD's own version string, from the data directory's ``card.json``.

    Read from the file, not inferred from a directory name. Returns "unknown" when the file cannot
    be read.
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
    symbol otherwise. Symbols alone undercount agreement: on five CAMI2 bins AMRFinderPlus calls
    ``cepA`` and RGI calls ``CepA-29``, the same beta-lactamase at 99.67% identity, and a symbol
    join reported one shared gene where two were shared.

    *aro_for* supplies ARO terms for AMRFinderPlus calls, from argNorm, keyed by (genome,
    lowercased symbol). RGI calls already carry ``aro``.
    """
    lookup = aro_for or {}
    by_tool: Dict[str, set] = {}
    for call in calls:
        symbol = call.symbol.lower()
        aro = call.aro or lookup.get((call.genome, symbol), "")
        by_tool.setdefault(call.tool, set()).add((call.genome, aro or symbol))
    shared: set = set.intersection(*by_tool.values()) if len(by_tool) > 1 else set()
    # Unique calls are reported by gene symbol, not by the ARO join key.
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
