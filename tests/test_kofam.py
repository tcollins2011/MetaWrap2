"""Tests for KO assignment and its threshold policy.

On 396 genomes with identical gene calls, eggNOG contradicted 35% of MicrobeAnnotator's
tool-specific KO calls and 7.5% of anvi'o's; the difference comes from what each accepts as a hit.
Checked: a curated threshold is compared against the score KOfam calibrated it on, a permissive call
never displaces a curated one, and every call records how it qualified.
"""

from __future__ import annotations

import pytest

from metawrap2 import kofam

#: Three rows in KOfam's ko_list layout: a domain-calibrated profile, a full-score one, and one
#: with no curated threshold (KOfam writes "-" in every threshold column).
KO_LIST = (
    "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq"
    "\tre/pos\tdefinition\n"
    "K00001\t320.40\tdomain\tall\t0.883405\t97\t92\t560\t288\t1.26\t0.589\talcohol dehydrogenase"
    " [EC:1.1.1.1]\n"
    "K00002\t457.17\tfull\tall\t0.447532\t3471\t3318\t7322\t528\t8.57\t0.590\talcohol"
    " dehydrogenase (NADP+) [EC:1.1.1.2]\n"
    "K00492\t-\t-\t-\t-\t1\t1\t396\t396\t0.39\t0.593\t1,3,7-trimethyluric acid 5-monooxygenase\n"
)


@pytest.fixture(scope="module")
def kofam_dir(tmp_path_factory):
    """A KOfam directory as `install-db kofam` lays it out: profiles/ and ko_list beside them."""
    directory = tmp_path_factory.mktemp("kofam")
    (directory / "profiles").mkdir()
    (directory / "ko_list").write_text(KO_LIST)
    return directory


def _profile(knum="K00001", threshold=100.0, score_type=kofam.ScoreType.FULL):
    return kofam.KoProfile(
        knum=knum,
        threshold=threshold,
        score_type=score_type,
        profile_type="all",
        f_measure=0.9,
        definition="test enzyme",
    )


def _hit(knum="K00001", full=200.0, domain=150.0, evalue=1e-50, gene="g1"):
    return kofam.Hit(gene_id=gene, knum=knum, full_score=full, domain_score=domain, evalue=evalue)


# ── the threshold list ───────────────────────────────────────────────────────────────────────


def test_the_thresholds_are_found_beside_the_profiles(kofam_dir):
    """Wherever KOFAM_DB points: the KOfam directory, or the profiles/ inside it."""
    expected = str(kofam_dir / "ko_list")
    assert kofam.find_ko_list(str(kofam_dir)) == expected
    assert kofam.find_ko_list(str(kofam_dir / "profiles")) == expected


# ── which score a threshold is compared against ──────────────────────────────────────────────


def test_a_domain_calibrated_threshold_is_compared_against_the_domain_score():
    """Domain-calibrated thresholds (8,697 profiles) are compared against the domain score.

    A multi-domain protein's full score exceeds its best domain's, so comparing it would accept
    hits KOfam rejects.
    """
    profile = _profile(threshold=180.0, score_type=kofam.ScoreType.DOMAIN)
    hit = _hit(full=200.0, domain=150.0)
    assert hit.score_for(profile.score_type) == 150.0
    assert kofam.assign_gene([hit], {profile.knum: profile}, kofam.ThresholdPolicy()) == []
    # A full-calibrated threshold is compared against the full score.
    profile = _profile(threshold=180.0, score_type=kofam.ScoreType.FULL)
    calls = kofam.assign_gene([hit], {profile.knum: profile}, kofam.ThresholdPolicy())
    assert [(c.knum, c.score) for c in calls] == [("K00001", 200.0)]


# ── the policy ───────────────────────────────────────────────────────────────────────────────


def test_a_curated_hit_is_never_displaced_by_a_permissive_one():
    """Permissive options add genes but do not replace curated calls.

    A gene with a curated hit keeps it even when another profile scores higher under the looser
    rule.
    """
    good = _profile("K00001", threshold=100.0)
    tempting = _profile("K00002", threshold=10_000.0)
    policy = kofam.ThresholdPolicy(relative_fraction=0.1)
    calls = kofam.assign_gene(
        [_hit("K00001", full=150.0), _hit("K00002", full=5000.0)],
        {"K00001": good, "K00002": tempting},
        policy,
    )
    assert [(c.knum, c.evidence) for c in calls] == [("K00001", kofam.Evidence.CURATED)]


# ── the search itself, where the threads go ──────────────────────────────────────────────────


def test_every_profile_file_is_searched_in_one_pass(tmp_path, monkeypatch):
    """All profile files are searched in one hmmsearch call, with the requested threads.

    KOfam has one file per KO, and hmmsearch parallelises across queries. One call per file (27,930
    single-query searches) used one core: on 2,339 genes it ran at 108% CPU and had not finished in
    ten minutes, against over 1000% batched.
    """
    import pyhmmer

    calls = []

    def fake_hmmsearch(queries, sequences, cpus=0, **kwargs):
        calls.append((list(queries), cpus))
        return iter(())

    faa = tmp_path / "p.faa"
    faa.write_text(">g1\nMKV\n")
    hmms = []
    for knum in ("K00001", "K00002", "K00003"):
        # A real profile file each, so HMMFile is exercised.
        path = tmp_path / ("%s.hmm" % knum)
        path.write_text(_MINIMAL_HMM % knum)
        hmms.append(str(path))

    monkeypatch.setattr(pyhmmer, "hmmsearch", fake_hmmsearch)
    assert list(kofam.search(str(faa), hmms, threads=8)) == []
    assert len(calls) == 1, "one search over every profile, not one search per file"
    queries, cpus = calls[0]
    assert len(queries) == 3
    assert cpus == 8


_MINIMAL_HMM = """HMMER3/f [3.3.2 | Nov 2020]
NAME  %s
LENG  2
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
STATS LOCAL MSV      -6.0  0.70
STATS LOCAL VITERBI  -6.5  0.70
STATS LOCAL FORWARD  -3.0  0.70
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.57544  1.78073  1.31293  1.75577  0.18968  0.00000        *
      1   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503      1 m - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.57544  1.78073  1.31293  1.75577  0.18968  0.00000        *
      2   2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503      2 m - - -
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503
          0.57544  1.78073  1.31293  1.75577  0.18968        *        *
//
"""
