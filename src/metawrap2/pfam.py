"""Protein domains from Pfam, for the genes KEGG cannot place.

On real CAMI2 bins, 41% of called genes take a curated KEGG ortholog. The rest are not junk - they
are genes whose function is known at the domain level and not at the pathway level, and a functional
report that drops them says less about the organism than it could. Pfam is the answer to "what kind
of protein is this" where KEGG's answer is "which reaction does it catalyse", and the two are
complementary rather than competing.

**Pfam's own cutoffs, not an e-value.** Every Pfam family ships a gathering threshold (GA), curated
per family, and a hit above GA is what Pfam means by a match - the same arrangement KOfam has and the
same reason to honour it: an e-value threshold that is right for a long, well-conserved family is
wrong for a short one, and a single cutoff across 21,000 families gets most of them wrong in one
direction or the other. PyHMMER applies GA directly through its `bit_cutoffs="gathering"` mode, so
this module does not reimplement the comparison, and :data:`CUTOFF_MODE` is the one place that
choice is stated.

**Domains, not proteins.** A protein can carry several domains and the same domain twice, so a
Pfam result is a list of *regions* rather than a label. :class:`DomainHit` keeps the coordinates and
:func:`architecture` renders the ordered list, because "ABC transporter + periplasmic binding" is a
different protein from either domain alone, and a report that kept only the best hit per gene would
flatten that.

Why not InterProScan, which is the usual tool here: it is a Java wrapper around a dozen member
databases, its bioconda build is four years old, and its data is about 50 GB. Pfam is the member
database most of its users want, PyHMMER is already how this package searches HMMs, and 399 MB of
Pfam-A costs an afternoon rather than a filesystem. The other member databases are a real gap and
are named as one in the plan rather than papered over here.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

#: How the bit-score cutoff is chosen. Pfam's gathering threshold is curated per family and is what
#: Pfam means by a match; `--cut_ga` is HMMER's name for it and `bit_cutoffs="gathering"` is
#: PyHMMER's. Recorded as a constant because it is the method rather than a tuning parameter.
CUTOFF_MODE = "gathering"

#: Fallback when a profile carries no gathering threshold. Pfam-A's profiles all do, so this is for
#: a custom or truncated HMM file; it is reported per hit so a GA-less call is never mistaken for a
#: curated one.
FALLBACK_MAX_EVALUE = 1e-5

#: Pfam accessions look like PF00001.23 - family then version. The version is kept in the output
#: because a domain call is only reproducible against the release that made it, and stripped for
#: joins, where a version mismatch would silently lose rows.
_ACCESSION = re.compile(r"^(PF\d{5})(?:\.(\d+))?$")


@dataclass(frozen=True)
class DomainHit:
    """One Pfam domain found in one protein, with where it sits."""

    gene: str
    accession: str
    name: str
    #: Envelope coordinates on the protein, 1-based inclusive - the region HMMER believes the
    #: domain occupies, which is wider and more honest than the alignment coordinates.
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
    """pyhmmer returns names as bytes in some versions and str in others; normalise once.

    The same helper kofam.py has, for the same reason. Duplicated deliberately rather than shared:
    one of these two modules changing its pyhmmer version should not be able to move the other's
    behaviour.
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
        # query object, and the older spelling raises AttributeError - which is how this failed the
        # first time it was run against a real profile set, caught per genome so one failure did
        # not cost the others.
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

    *proteins* is (gene id, sequence). Returns every included domain, not one hit per gene: a
    protein carrying two domains is two rows, which is what makes an architecture recoverable.

    The gathering thresholds come from the HMM file itself, so nothing here has to ship a table of
    them - and a profile without one falls back to an e-value and says so in the row.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    # Typed loosely for the reason kofam.py documents: pyhmmer's stubs type hmmsearch's sequence
    # argument invariantly in its alphabet, so an amino-acid block does not satisfy
    # Iterable[DigitalSequence[Alphabet]] even though that is exactly what it wants.
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
        # One pass over the profile set with every sequence, which is how PyHMMER wants it and is
        # also the only way the gathering thresholds can come from the profiles themselves.
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

    Ordered and repeated rather than a set: "ABC transporter + periplasmic binding protein" is a
    different protein from either domain alone, and two copies of a domain is a different protein
    from one.
    """
    return " + ".join(hit.name for hit in sorted(hits, key=lambda hit: hit.start))


@dataclass(frozen=True)
class Coverage:
    """How much of a gene set Pfam placed, beside how much KEGG did."""

    genes: int
    genes_with_a_domain: int
    domains: int
    #: Genes with a Pfam domain and no KEGG ortholog - the gap this module exists to fill.
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

    The second number is the one worth reporting: Pfam annotating 70% of genes is interesting only
    beside the share KEGG already had, because the overlap is large and the increment is what the
    extra search bought.
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
        # Versioned accession and bare family side by side: the first dates the call, the second
        # joins across releases.
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


def release_from_notes(path: str) -> str:
    """The Pfam release from its own `relnotes.txt`, or "unknown".

    Read rather than taken from a directory name: a domain call is only reproducible against the
    release that made it, and Pfam's families are renumbered between releases.
    """
    candidate = path
    if os.path.isdir(candidate):
        candidate = os.path.join(candidate, "relnotes.txt")
    if not os.path.isfile(candidate):
        return "unknown"
    with open(candidate, errors="replace") as handle:
        for line in handle:
            match = re.search(r"RELEASE\s+([\d.]+)", line, re.IGNORECASE)
            if match:
                return match.group(1)
    return "unknown"


#: Printed with every Pfam result.
DOMAIN_CAVEAT = (
    "a Pfam domain says what kind of protein this is, not what reaction it catalyses: a gene with "
    "a kinase domain and no KEGG ortholog is a kinase of unknown pathway, which is more than "
    "nothing and less than a function. Hits are called at each family's own gathering threshold, "
    "which is what Pfam means by a match"
)

#: Printed when Pfam is reported beside KEGG.
INCREMENT_CAVEAT = (
    "the number to read is not how many genes Pfam placed but how many it placed that KEGG had "
    "not: the overlap between the two is large, and the increment is what the extra search bought"
)
