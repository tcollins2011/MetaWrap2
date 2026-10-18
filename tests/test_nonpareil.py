"""Tests for metawrap2.nonpareil: the port of Nonpareil's model, and the module around it.

The first test is the one that matters. `metawrap2.nonpareil` reimplements a published numerical
method, and the only evidence that a reimplementation is right is agreement with the original on
data neither was tuned to. The expected values in `test_matches_r_package` came out of the R
package itself - see tests/fixtures/nonpareil/SOURCE.md for the invocation - not out of this code.
"""

from __future__ import annotations

import argparse
import os

import pytest

from metawrap2 import nonpareil as np_
from metawrap2.modules import nonpareil as module

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "nonpareil")
DATASETS = ("HumanGut", "IowaSoil", "LakeLanier")


def npo(name: str) -> str:
    return os.path.join(FIXTURES, name + ".npo")


# ─── the port against the original ───────────────────────────────────────────────────────


def test_matches_r_package():
    """Six significant figures against R's Nonpareil 3.5.5 on its own example data.

    Not "close enough": every quantity the pipeline reports is checked, including both gamma
    parameters, because C can agree while the fit behind it differs - and it is the fit that
    produces the diversity index and the projection.
    """
    expected = {
        "HumanGut": {
            "C": 0.806352,
            "kappa": 0.790210,
            "diversity": 17.260378,
            "LRstar": 2.62214e09,
            "r": 0.998624,
            "a": 54.2436,
            "b": 3.08473,
        },
        "IowaSoil": {
            "C": 0.557296,
            "kappa": 0.527510,
            "diversity": 22.795958,
            "LRstar": 6.62396e11,
            "r": 0.999673,
            "a": 88.9061,
            "b": 3.85622,
        },
        "LakeLanier": {
            "C": 0.683087,
            "kappa": 0.659060,
            "diversity": 19.213483,
            "LRstar": 2.05296e10,
            "r": 0.998621,
            "a": 62.8164,
            "b": 3.21734,
        },
    }
    for name in DATASETS:
        _check_against_r(name, expected[name])


def _check_against_r(name, expected):
    curve = np_.curve_from_npo(npo(name))
    assert curve.has_model
    assert curve.C == pytest.approx(expected["C"], rel=1e-5)
    assert curve.kappa == pytest.approx(expected["kappa"], rel=1e-5)
    assert curve.diversity == pytest.approx(expected["diversity"], rel=1e-5)
    assert curve.LRstar == pytest.approx(expected["LRstar"], rel=1e-4)
    assert curve.model_r == pytest.approx(expected["r"], rel=1e-5)
    assert curve.a == pytest.approx(expected["a"], rel=1e-4)
    assert curve.b == pytest.approx(expected["b"], rel=1e-4)


# ─── parsing ─────────────────────────────────────────────────────────────────────────────


def test_kmer_kernel_is_recognised_by_ksize(tmp_path):
    path = tmp_path / "k.npo"
    path.write_text(
        "# @version: 3.5.5\n# @ksize: 24\n# @L: 150\n# @AL: 148.5\n# @R: 1000\n"
        "# @logsampling: 0.7\n" + "\n".join("%d\t0.1\t0.01\t0\t0.1\t0.2" % x for x in (10, 20, 30))
    )
    curve = np_.read_npo(str(path))
    assert curve.kernel == "kmer"
    # The kmer kernel's overlap is fixed at 50 and its length comes from AL, not L: with the
    # alignment kernel's rules this would use L=150 and whatever overlap the header lacked.
    assert curve.overlap == 50.0
    # AL under the kmer kernel is k-mers per read (L - k + 1), not a length in bp, so the effort
    # axis and the projection are in k-mers. Reading it as a read length overstates the sequencing
    # already done - here by a factor of three.
    assert curve.read_length == pytest.approx(150.0)
    assert curve.effort_length == pytest.approx(148.5)
    assert curve.effort_unit == "kmers"
    assert curve.sequenced_bp == pytest.approx(150_000)
    assert curve.effort_at_full_depth == pytest.approx(148_500)


# ─── kernel selection, which has to be one choice per cohort ─────────────────────────────


def _lengths(mean, short=0):
    return np_.ReadLengths(
        scanned=1000,
        minimum=int(mean) - 10,
        maximum=int(mean) + 10,
        mean=mean,
        shorter_than_k=short,
    )


def test_auto_picks_alignment_when_any_sample_is_short():
    """One kernel for the cohort, chosen by its shortest-read member.

    The two kernels disagree by 33 points of coverage on the same reads, so a cohort where one
    sample used each would have coverage figures that cannot be compared - which is the only thing
    this module is for.
    """
    args = argparse.Namespace(kernel="auto", kmer_size=24)
    chosen = module._choose_kernel(args, {"long": _lengths(150.0), "short": _lengths(34.0, 81)})
    assert chosen == "alignment"
    # An explicit kernel is never overridden.
    for forced in ("kmer", "alignment"):
        args = argparse.Namespace(kernel=forced, kmer_size=24)
        assert module._choose_kernel(args, {"short": _lengths(30.0, 100)}) == forced


# ── composing with read_qc, which was impossible ─────────────────────────────────────────────


def test_colliding_basenames_fall_back_to_the_directory_name():
    """`read_qc` writes `<out>/<sample>/final_pure_reads_1.fastq`, so a cohort's read files all
    share one basename and the sample identity lives in the directory.

    That made the composition this project documents - `nonpareil -r
    READ_QC/*/final_pure_reads_1.fastq` - fail every time with a correct but unhelpful message
    about two files reducing to one name. The two modules could not be chained the way the docs
    chain them. The directory is what the user meant.
    """
    from metawrap2.modules import nonpareil as module

    paths = [
        "READ_QC/ERR011347/final_pure_reads_1.fastq",
        "READ_QC/ERR011348/final_pure_reads_1.fastq",
        "READ_QC/ERR011349/final_pure_reads_1.fastq",
    ]
    assert module.sample_names(paths) == ["ERR011347", "ERR011348", "ERR011349"]


def test_names_that_do_not_collide_are_left_alone():
    """A cohort whose files are named per sample keeps the names it chose."""
    from metawrap2.modules import nonpareil as module

    paths = ["reads/gut_early_1.fastq.gz", "reads/gut_late_1.fastq.gz"]
    assert module.sample_names(paths) == ["gut_early", "gut_late"]


def test_only_the_colliding_names_are_rewritten():
    """One collision in a cohort must not rename the files that were already unambiguous."""
    from metawrap2.modules import nonpareil as module

    paths = [
        "READ_QC/ERR011347/final_pure_reads_1.fastq",
        "READ_QC/ERR011348/final_pure_reads_1.fastq",
        "reads/sampC_1.fastq",
    ]
    assert module.sample_names(paths) == ["ERR011347", "ERR011348", "sampC"]


def test_a_collision_the_directory_cannot_break_is_still_a_collision():
    """Guessing a third time would be inventing a name, and two rows with one label are not data."""
    from metawrap2.modules import nonpareil as module

    paths = ["a/reads_1.fastq", "a/reads_1.fastq.gz"]
    assert module.sample_names(paths) == ["a", "a"]
