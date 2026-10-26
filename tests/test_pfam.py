"""Tests for Pfam domain handling: cutoffs, architecture and coverage reporting.

The search needs Pfam-A on disk, so these tests cover the code around it: hits to rows, domains to
architectures, and coverage beside KEGG's. A real run on one CAMI2 bin (in the release notes) gave
6,671 domains over 3,185 of 3,891 genes.
"""

from __future__ import annotations

from metawrap2 import pfam

# ── cutoffs ──────────────────────────────────────────────────────────────────────────────────


# ── accessions ───────────────────────────────────────────────────────────────────────────────


# ── domain architecture ──────────────────────────────────────────────────────────────────────


# ── coverage beside KEGG ─────────────────────────────────────────────────────────────────────


# ── release ──────────────────────────────────────────────────────────────────────────────────


# ── parsing pyhmmer results ──────────────────────────────────────────────────────────────────


class _Domain:
    def __init__(self, start, end, score):
        self.env_from, self.env_to, self.score, self.i_evalue = start, end, score, 1e-18


class _Domains:
    def __init__(self, domains):
        self.included = domains


class _Hit:
    def __init__(self, name, domains):
        self.name = name.encode()
        self.domains = _Domains(domains)


class _Query:
    def __init__(self, name, accession):
        self.name = name.encode()
        self.accession = accession.encode() if accession else None


class _TopHits:
    def __init__(self, name, accession, hits):
        self.query = _Query(name, accession)
        self._hits = hits

    def __iter__(self):
        return iter(self._hits)


def test_tophits_are_read_through_the_query_object():
    """pyhmmer 0.11 replaced `query_name`/`query_accession` with `query`.

    The old attributes raise AttributeError, which, caught per genome, gave zero domains.
    """
    results = [
        _TopHits("7tm_1", "PF00001.23", [_Hit("gene_1", [_Domain(10, 120, 55.0)])]),
        _TopHits("HTH", None, [_Hit("gene_2", [_Domain(5, 60, 30.0), _Domain(120, 175, 28.0)])]),
    ]
    hits = list(pfam._hits_from_hmmer(results, {}, {"PF00001.23": True}))
    assert len(hits) == 3
    assert hits[0].accession == "PF00001.23"
    assert hits[0].name == "7tm_1"
    assert hits[0].gene == "gene_1"
    # A profile without an accession uses its name.
    assert hits[1].accession == "HTH"
    assert [h.start for h in hits[1:]] == [5, 120]
