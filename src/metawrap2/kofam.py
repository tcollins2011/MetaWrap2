"""KEGG ortholog assignment for proteins, and the threshold policy that decides what is a hit.

Module completeness, the functional hierarchy and the per-sample distribution are all sums over
KO assignments made here. The design follows a comparison on 396 genomes with identical gene calls
(Kananen et al., *Bioinform Adv* 2025, 10.1093/bioadv/vbaf039):

=============================  ========================  ============================  =========================
tool                           calls confirmed by eggNOG  tool-specific calls confirmed  ... contradicted
=============================  ========================  ============================  =========================
KofamScan                      84%                        -                              -
anvi'o (adaptive thresholds)   79%                        37%                            7.5%
MicrobeAnnotator               74%                        0.3%                           **35%**
=============================  ========================  ============================  =========================

Design rules:

1. Use KOfam's curated per-profile thresholds against the score each one is calibrated for. KOfam
   gives a bit-score threshold per KO and its score type: ``full`` for 17,890 profiles, ``domain``
   for 8,697. Comparing against the wrong score moves calls in both directions; see
   :data:`ScoreType`.

2. 1,866 of 28,453 profiles (6.6%) have no threshold. M00851 "Carbapenem resistance" is a flat OR
   of 21 KOs, one of them thresholdless, so a single permissive call marks the module complete. By
   default these profiles are not assigned, and when enabled each call is flagged. Unflagged
   thresholdless hits are what MicrobeAnnotator's 35% contradicted rate consists of.

3. An optional relative fallback, labelled per call. anvi'o's adaptive heuristic finds 13% more
   genes per genome (p<2.2e-16), mostly in divergent lineages (38% of Nanosynbacteraceae cases vs
   4.9% of Enterobacteriaceae), where thresholds calibrated on cultured relatives are too strict.
   A gene with no curated-passing hit can take its best hit if it reaches a set fraction of the
   profile's threshold; the call carries ``evidence="relative"`` so it can be counted separately.

Nothing is downloaded here. The profiles (``profiles.tar.gz``, 1.5 GB) are a database configured
in ``metawrap2.toml`` and checked by ``doctor`` like the others (see :mod:`metawrap2.config`).
"""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from .io.seqio import smart_open

#: Where `metawrap2 install-db kofam` downloads the profiles and their thresholds from.
KOFAM_URL = "https://www.genome.jp/ftp/db/kofam/"

#: md5 of the threshold list this release's policy and tests were calibrated on (KOfam's
#: ko_list.gz of 2026-09-24). KOfam revises thresholds, which changes calls. A different list runs
#: with a warning, and every row records which list it used.
KO_LIST_MD5 = "e037ff1135d820c775da78e6d5f256c7"


def find_ko_list(kofam_db: str) -> Optional[str]:
    """KOfam's ``ko_list`` next to the profiles: in the KOfam directory, or one level above."""
    base = kofam_db if os.path.isdir(kofam_db) else os.path.dirname(kofam_db)
    for directory in (base, os.path.dirname(base)):
        for name in ("ko_list.gz", "ko_list"):
            path = os.path.join(directory, name)
            if os.path.isfile(path):
                return path
    return None


def ko_list_md5(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def ko_list_label(path: str) -> str:
    """Value for each row's ``ko_list_snapshot``: a prefix of the list's md5."""
    return "md5:%s" % ko_list_md5(path)[:8]


#: What KOfam writes in ``threshold`` when a profile has no curated cutoff.
NO_THRESHOLD = "-"


class ScoreType(Enum):
    """Which HMMER score a profile's threshold is calibrated against.

    ``full`` is the summed score over all domains of the alignment; ``domain`` is the best single
    domain's score. They differ a lot for multi-domain proteins, and KOfam sets it per profile.
    """

    FULL = "full"
    DOMAIN = "domain"


@dataclass(frozen=True)
class KoProfile:
    """One KOfam profile's metadata: what to compare, against what, and what the KO is."""

    knum: str
    #: None when KOfam publishes no cutoff for this profile.
    threshold: Optional[float]
    score_type: Optional[ScoreType]
    #: ``all`` or ``trim``; recorded because it describes how the profile was built.
    profile_type: str
    #: KOfam's F-measure for this profile at its threshold. A low value means the profile separates
    #: its family poorly even at the curated cutoff.
    f_measure: Optional[float]
    definition: str

    @property
    def has_threshold(self) -> bool:
        return self.threshold is not None


def load_ko_list(path: str) -> Dict[str, KoProfile]:
    """Parse KOfam's ``ko_list``. Keyed by K number."""
    profiles: Dict[str, KoProfile] = {}
    with smart_open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for line in fh:
            cut = line.rstrip("\n").split("\t")
            if len(cut) < len(header):
                continue
            knum = cut[index["knum"]]
            raw_threshold = cut[index["threshold"]].strip()
            raw_score = cut[index["score_type"]].strip()
            threshold = None if raw_threshold in ("", NO_THRESHOLD) else float(raw_threshold)
            score_type = None if raw_score in ("", NO_THRESHOLD) else ScoreType(raw_score)
            try:
                f_measure: Optional[float] = float(cut[index["F-measure"]])
            except (ValueError, KeyError):
                f_measure = None
            profiles[knum] = KoProfile(
                knum=knum,
                threshold=threshold,
                score_type=score_type,
                profile_type=cut[index["profile_type"]].strip(),
                f_measure=f_measure,
                definition=cut[index["definition"]].strip(),
            )
    return profiles


# ── hits and the policy that filters them ────────────────────────────────────────────────────


@dataclass(frozen=True)
class Hit:
    """One profile against one protein, with both full and best-domain scores."""

    gene_id: str
    knum: str
    full_score: float
    domain_score: float
    evalue: float

    def score_for(self, score_type: Optional[ScoreType]) -> float:
        """The score this profile's threshold is calibrated against.

        Defaults to the full score (which is >= the best domain score) when the profile has no
        score type. Such profiles also have no threshold, so this never decides a curated call.
        """
        if score_type is ScoreType.DOMAIN:
            return self.domain_score
        return self.full_score


class Evidence(Enum):
    """How a KO assignment was made. Carried on every call through the pipeline."""

    #: Passed the profile's curated threshold.
    CURATED = "curated"
    #: Best hit for a gene that nothing curated passed, reaching a stated fraction of the threshold.
    RELATIVE = "relative"
    #: A profile with no curated threshold, accepted on E-value alone.
    NO_THRESHOLD = "no_threshold"


@dataclass(frozen=True)
class ThresholdPolicy:
    """What counts as a KO assignment. Every field affects the result and is recorded.

    The defaults are strict (curated thresholds only, thresholdless profiles dropped), which
    reproduces KofamScan, the tool with the highest confirmation rate in the benchmark above.
    """

    #: Accept the best hit of an otherwise-unannotated gene when it reaches this fraction of the
    #: profile's curated threshold. None disables the fallback.
    relative_fraction: Optional[float] = None
    #: Accept hits to profiles that have no curated threshold when their E-value is at most this.
    #: None (the default) drops them; see the module docstring on M00851.
    thresholdless_max_evalue: Optional[float] = None
    #: E-value ceiling for every call. A hit that passes a bit-score threshold with a poor E-value
    #: is usually a short alignment to a long profile.
    max_evalue: float = 1e-5
    #: Keep at most this many KOs per gene, best first. A gene can match several profiles of one
    #: family, and keeping all inflates downstream counts; hierarchy.weights_from_gene_kos splits
    #: a gene's weight across the kept KOs.
    max_kos_per_gene: int = 1

    def describe(self) -> str:
        """One-line policy description recorded in each output row."""
        parts = ["curated"]
        if self.relative_fraction is not None:
            parts.append("relative>=%.2fx" % self.relative_fraction)
        if self.thresholdless_max_evalue is not None:
            parts.append("thresholdless<=%g" % self.thresholdless_max_evalue)
        parts.append("max_evalue<=%g" % self.max_evalue)
        parts.append("top%d" % self.max_kos_per_gene)
        return ",".join(parts)


@dataclass(frozen=True)
class KoCall:
    """One accepted KO assignment for one gene."""

    gene_id: str
    knum: str
    score: float
    threshold: Optional[float]
    evalue: float
    evidence: Evidence
    definition: str
    f_measure: Optional[float]


def assign_gene(
    hits: Sequence[Hit],
    profiles: Mapping[str, KoProfile],
    policy: ThresholdPolicy,
) -> List[KoCall]:
    """Accepted KO calls for one gene, best first.

    Branches are tried in order: curated hits, then the relative fallback, then thresholdless
    profiles. A later branch only runs when the earlier ones produce nothing, so a permissive
    option never replaces a curated call.
    """
    usable = [h for h in hits if h.evalue <= policy.max_evalue]
    if not usable:
        return []

    curated: List[KoCall] = []
    for hit in usable:
        profile = profiles.get(hit.knum)
        if profile is None or not profile.has_threshold:
            continue
        score = hit.score_for(profile.score_type)
        if score >= profile.threshold:  # type: ignore[operator]
            curated.append(
                KoCall(
                    gene_id=hit.gene_id,
                    knum=hit.knum,
                    score=score,
                    threshold=profile.threshold,
                    evalue=hit.evalue,
                    evidence=Evidence.CURATED,
                    definition=profile.definition,
                    f_measure=profile.f_measure,
                )
            )
    if curated:
        curated.sort(key=lambda c: (-c.score, c.evalue, c.knum))
        return curated[: policy.max_kos_per_gene]

    if policy.relative_fraction is not None:
        relative: List[KoCall] = []
        for hit in usable:
            profile = profiles.get(hit.knum)
            if profile is None or not profile.has_threshold:
                continue
            score = hit.score_for(profile.score_type)
            if score >= policy.relative_fraction * profile.threshold:  # type: ignore[operator]
                relative.append(
                    KoCall(
                        gene_id=hit.gene_id,
                        knum=hit.knum,
                        score=score,
                        threshold=profile.threshold,
                        evalue=hit.evalue,
                        evidence=Evidence.RELATIVE,
                        definition=profile.definition,
                        f_measure=profile.f_measure,
                    )
                )
        if relative:
            relative.sort(key=lambda c: (-c.score, c.evalue, c.knum))
            return relative[: policy.max_kos_per_gene]

    if policy.thresholdless_max_evalue is not None:
        loose: List[KoCall] = []
        for hit in usable:
            profile = profiles.get(hit.knum)
            if profile is None or profile.has_threshold:
                continue
            if hit.evalue <= policy.thresholdless_max_evalue:
                loose.append(
                    KoCall(
                        gene_id=hit.gene_id,
                        knum=hit.knum,
                        score=hit.full_score,
                        threshold=None,
                        evalue=hit.evalue,
                        evidence=Evidence.NO_THRESHOLD,
                        definition=profile.definition,
                        f_measure=profile.f_measure,
                    )
                )
        if loose:
            loose.sort(key=lambda c: (c.evalue, -c.score, c.knum))
            return loose[: policy.max_kos_per_gene]

    return []


def assign(
    hits: Iterable[Hit],
    profiles: Mapping[str, KoProfile],
    policy: Optional[ThresholdPolicy] = None,
) -> Iterator[KoCall]:
    """Group hits by gene and apply the policy to each gene's hits."""
    rules = policy if policy is not None else ThresholdPolicy()
    per_gene: Dict[str, List[Hit]] = {}
    for hit in hits:
        per_gene.setdefault(hit.gene_id, []).append(hit)
    for gene_hits in per_gene.values():
        yield from assign_gene(gene_hits, profiles, rules)


#: Columns of the per-gene KO table.
KO_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene_id",
    "ko",
    "score",
    "threshold",
    "evalue",
    "evidence",
    "f_measure",
    "definition",
    "threshold_policy",
    "ko_list_snapshot",
)


def ko_row(
    genome: str, call: KoCall, policy: ThresholdPolicy, ko_list: str = "unknown"
) -> Dict[str, object]:
    return {
        "genome": genome,
        "gene_id": call.gene_id,
        "ko": call.knum,
        "score": "%.2f" % call.score,
        "threshold": "" if call.threshold is None else "%.2f" % call.threshold,
        "evalue": "%.3g" % call.evalue,
        "evidence": call.evidence.value,
        "f_measure": "" if call.f_measure is None else "%.4f" % call.f_measure,
        "definition": call.definition,
        "threshold_policy": policy.describe(),
        "ko_list_snapshot": ko_list,
    }


# ── searching ────────────────────────────────────────────────────────────────────────────────


def _as_text(value) -> str:
    """Normalise a pyhmmer name, which is bytes in some versions and str in others."""
    if value is None:
        return ""
    return value.decode() if isinstance(value, bytes) else str(value)


def search(
    protein_fasta: str,
    hmm_files: Sequence[str],
    threads: int = 1,
    max_evalue: float = 1e-3,
) -> Iterator[Hit]:
    """Search proteins against HMM profiles with PyHMMER, yielding both scores per hit.

    PyHMMER is used instead of an ``hmmsearch`` subprocess (as Pyrodigal is in
    :mod:`metawrap2.genecall`) so the threshold policy is tested in code instead of parsed from
    ``--tblout``. It releases the GIL, so *threads* runs in parallel.

    *max_evalue* is a pre-filter to limit the hit list. It is looser than
    :attr:`ThresholdPolicy.max_evalue` so the policy sees every hit its cutoffs could accept.

    Each hit carries both the full-sequence and best-domain scores; see :class:`ScoreType`.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(protein_fasta, digital=True, alphabet=alphabet) as handle:
        # Typed Any because pyhmmer's stubs type hmmsearch's sequence argument invariantly in its
        # alphabet, so an amino-acid block does not satisfy Iterable[DigitalSequence[Alphabet]].
        sequences: Any = list(handle.read_block())

    # One search over all profiles. hmmsearch parallelises across queries, so looping over KOfam's
    # 27,930 single-profile files runs on one core. On 2,339 genes the per-file loop ran at 108%
    # CPU and had not finished after 10 minutes; one search takes about 6 minutes for 14,171 genes.
    def profiles() -> Iterator[Any]:
        # One file open at a time, to stay under the descriptor limit. Each HMM is fully parsed
        # when yielded, so it outlives its handle.
        for hmm_file in hmm_files:
            with pyhmmer.plan7.HMMFile(hmm_file) as handle:
                yield from handle

    queries: Any = profiles()
    for top_hits in pyhmmer.hmmsearch(queries, sequences, cpus=threads, E=max_evalue):
        knum = _as_text(top_hits.query.name)
        for hit in top_hits:
            if hit.evalue > max_evalue:
                continue
            best_domain = max((domain.score for domain in hit.domains), default=hit.score)
            yield Hit(
                gene_id=_as_text(hit.name),
                knum=knum,
                full_score=float(hit.score),
                domain_score=float(best_domain),
                evalue=float(hit.evalue),
            )
