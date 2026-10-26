"""Protein domains from Pfam, for the genes KEGG cannot place.

On CAMI2 bins, 41% of called genes get a curated KEGG ortholog. Many of the rest have a function
known at the domain level but not the pathway level. Pfam describes what kind of protein a gene
encodes; KEGG describes which reaction it catalyses.

**Cutoffs are Pfam's gathering thresholds (GA).** Each family has a curated GA, and a hit above it
is what Pfam defines as a match. A single e-value cutoff across 21,000 families would be too
strict for some and too loose for others. PyHMMER applies GA through
`bit_cutoffs="gathering"`; :data:`CUTOFF_MODE` records the choice.

**Results are domains, not proteins.** A protein can carry several domains, or the same domain
twice. :class:`DomainHit` keeps the coordinates and :func:`architecture` renders the ordered list.

InterProScan is not used: it is a Java wrapper around a dozen member databases, its bioconda
build is four years old, and its data is about 50 GB. Pfam-A is 399 MB and is searched with
PyHMMER, which this package already uses. The other member databases are not covered.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

#: How the bit-score cutoff is chosen: Pfam's per-family gathering threshold (HMMER `--cut_ga`,
#: PyHMMER `bit_cutoffs="gathering"`).
CUTOFF_MODE = "gathering"

#: Fallback when a profile has no gathering threshold. All Pfam-A profiles have one, so this is for
#: custom or truncated HMM files; each such hit is flagged (``curated_cutoff=False``).
FALLBACK_MAX_EVALUE = 1e-5

#: Pfam accessions look like PF00001.23 (family, then version). The version is kept in the output,
#: since families change between releases, and stripped for joins, so a version mismatch does not
#: drop rows.
_ACCESSION = re.compile(r"^(PF\d{5})(?:\.(\d+))?$")


@dataclass(frozen=True)
class DomainHit:
    """One Pfam domain found in one protein, with where it sits."""

    gene: str
    accession: str
    name: str
    #: Envelope coordinates on the protein, 1-based inclusive: the region HMMER assigns to the
    #: domain, wider than the alignment coordinates.
    start: int
    end: int
    bit_score: float
    evalue: float
    #: False when the profile had no gathering threshold and the e-value fallback was used.
    curated_cutoff: bool

    @property
    def family(self) -> str:
        """The accession without its version, for joining across releases."""
        match = _ACCESSION.match(self.accession)
        return match.group(1) if match else self.accession

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)


def _as_text(value) -> str:
    """Return pyhmmer names as str; some pyhmmer versions return bytes.

    kofam.py has the same helper. It is duplicated so that a pyhmmer change in one module cannot
    affect the other.
    """
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else str(value)


def _hits_from_hmmer(
    hmmer_hits: Iterable,
    gene_of: Dict[str, str],
    curated: Dict[str, bool],
) -> Iterator[DomainHit]:
    """Turn PyHMMER's TopHits into domain hits. Separated so the parsing can be tested alone."""
    for top in hmmer_hits:
        # `top.query`, not `top.query_name`: pyhmmer 0.11 replaced the flat attributes with the
        # query object, and the older name raises AttributeError.
        query = _as_text(getattr(top.query, "name", None))
        accession = _as_text(getattr(top.query, "accession", None)) or query
        for hit in top:
            gene = _as_text(hit.name)
            for domain in hit.domains.included:
                yield DomainHit(
                    gene=gene_of.get(gene, gene),
                    accession=accession,
                    name=query,
                    start=domain.env_from,
                    end=domain.env_to,
                    bit_score=domain.score,
                    evalue=domain.i_evalue,
                    curated_cutoff=curated.get(accession, True),
                )


def search(
    proteins: Sequence[Tuple[str, str]],
    hmm_path: str,
    threads: int = 1,
    max_evalue: float = FALLBACK_MAX_EVALUE,
) -> List[DomainHit]:
    """Search *proteins* against Pfam-A, honouring each family's gathering threshold.

    *proteins* is (gene id, sequence). Returns every included domain, so a protein with two
    domains gives two rows.

    The gathering thresholds are read from the HMM file. A profile without one falls back to an
    e-value cutoff, flagged in the row.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    # Typed loosely, as in kofam.py: pyhmmer's stubs type hmmsearch's sequence argument
    # invariantly in its alphabet, so an amino-acid block does not satisfy
    # Iterable[DigitalSequence[Alphabet]].
    sequences: Any = [
        pyhmmer.easel.TextSequence(name=gene.encode(), sequence=sequence).digitize(alphabet)
        for gene, sequence in proteins
    ]
    gene_of = {gene: gene for gene, _ in proteins}
    curated: Dict[str, bool] = {}
    hits: List[DomainHit] = []
    with pyhmmer.plan7.HMMFile(hmm_path) as handle:
        profiles = []
        for hmm in handle:
            accession = _as_text(hmm.accession or hmm.name)
            curated[accession] = hmm.cutoffs.gathering_available()
            profiles.append(hmm)
        # One pass over the profile set with all sequences, which lets the gathering thresholds
        # come from the profiles.
        results = pyhmmer.hmmsearch(
            profiles,
            sequences,
            cpus=threads,
            bit_cutoffs=CUTOFF_MODE,
        )
        hits.extend(_hits_from_hmmer(results, gene_of, curated))
    return hits


def architecture(hits: Sequence[DomainHit]) -> str:
    """One protein's domains in the order they appear, as `A + B + A`.

    Ordered and with repeats, since the order and count of domains distinguish proteins.
    """
    return " + ".join(hit.name for hit in sorted(hits, key=lambda hit: hit.start))


@dataclass(frozen=True)
class Coverage:
    """How much of a gene set Pfam placed, beside how much KEGG did."""

    genes: int
    genes_with_a_domain: int
    domains: int
    #: Genes with a Pfam domain and no KEGG ortholog.
    genes_rescued: Optional[int] = None

    @property
    def fraction_with_a_domain(self) -> Optional[float]:
        return self.genes_with_a_domain / self.genes if self.genes else None

    @property
    def mean_domains_per_annotated_gene(self) -> Optional[float]:
        if not self.genes_with_a_domain:
            return None
        return self.domains / self.genes_with_a_domain


def coverage(
    hits: Sequence[DomainHit], genes: Sequence[str], ko_genes: Optional[Sequence[str]] = None
) -> Coverage:
    """Count what Pfam placed, and how much of it KEGG had not.

    The increment over KEGG is the useful figure, because the two sets overlap heavily.
    """
    annotated = {hit.gene for hit in hits}
    rescued = None
    if ko_genes is not None:
        rescued = len(annotated - set(ko_genes))
    return Coverage(
        genes=len(set(genes)),
        genes_with_a_domain=len(annotated),
        domains=len(hits),
        genes_rescued=rescued,
    )


#: Columns of the per-domain table.
DOMAIN_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "pfam_accession",
    "pfam_family",
    "domain_name",
    "start",
    "end",
    "length",
    "bit_score",
    "evalue",
    "cutoff",
)


def domain_row(hit: DomainHit, genome: str) -> Dict[str, str]:
    return {
        "genome": genome,
        "gene": hit.gene,
        # Versioned accession for the release, bare family for joins across releases.
        "pfam_accession": hit.accession,
        "pfam_family": hit.family,
        "domain_name": hit.name,
        "start": str(hit.start),
        "end": str(hit.end),
        "length": str(hit.length),
        "bit_score": "%.1f" % hit.bit_score,
        "evalue": "%.2g" % hit.evalue,
        "cutoff": CUTOFF_MODE if hit.curated_cutoff else "evalue_fallback",
    }


#: Columns of the per-gene architecture table.
ARCHITECTURE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "domains",
    "architecture",
)


def architecture_rows(hits: Sequence[DomainHit], genome: str) -> List[Dict[str, str]]:
    """One row per gene, with its ordered domain architecture."""
    by_gene: Dict[str, List[DomainHit]] = {}
    for hit in hits:
        by_gene.setdefault(hit.gene, []).append(hit)
    return [
        {
            "genome": genome,
            "gene": gene,
            "domains": str(len(gene_hits)),
            "architecture": architecture(gene_hits),
        }
        for gene, gene_hits in sorted(by_gene.items())
    ]


#: Name of Pfam's release notes file, the only file read for the version.
NOTES_NAME = "relnotes.txt"

#: Size limit for the notes file. Release notes are a few KB; Pfam-A.hmm is 2.2 GB. Larger files are
#: not scanned, because the configured PFAM_DB is usually the HMM file and scanning it took minutes.
NOTES_MAX_BYTES = 1 << 20


def release_from_notes(path: str) -> str:
    """The Pfam release from its own `relnotes.txt`, or "unknown".

    Read from the notes, not the directory name, because Pfam renumbers families between
    releases.

    *path* may be the notes file, the directory holding it, or a sibling file in that directory.
    The last is the usual case: ``install-db pfam`` downloads relnotes.txt and records
    ``<dir>/Pfam-A.hmm`` in the config (``key_file="Pfam-A.hmm"``).
    """
    candidate = path
    if os.path.isdir(candidate):
        candidate = os.path.join(candidate, NOTES_NAME)
    elif os.path.basename(candidate) != NOTES_NAME:
        candidate = os.path.join(os.path.dirname(candidate), NOTES_NAME)
    if not os.path.isfile(candidate):
        return "unknown"
    try:
        if os.path.getsize(candidate) > NOTES_MAX_BYTES:
            return "unknown"
    except OSError:
        return "unknown"
    with open(candidate, errors="replace") as handle:
        for line in handle:
            match = re.search(r"RELEASE\s+([\d.]+)", line, re.IGNORECASE)
            if match:
                return match.group(1)
    return "unknown"


#: Printed with every Pfam result.
DOMAIN_CAVEAT = (
    "a Pfam domain gives the protein type, not the reaction it catalyses: a gene with "
    "a kinase domain and no KEGG ortholog is a kinase of unknown pathway. Hits are called at "
    "each family's own gathering threshold, "
    "Pfam's definition of a match"
)

#: Printed when Pfam is reported beside KEGG.
INCREMENT_CAVEAT = (
    "the useful number is how many genes Pfam annotated that KEGG did "
    "not; the two overlap heavily, so the total Pfam count overstates what the search adds"
)
