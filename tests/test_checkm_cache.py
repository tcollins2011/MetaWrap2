"""Making consolidation cheaper without changing a single reported number.

bin_refinement scores 2^n-1 hybrid bin sets and CheckM is the whole cost of it. Two measured
redundancies, both removable:

* 32% of the bins across the seven sets of a three-binner run are the *same bin* scored again - same
  contigs, so necessarily the same completeness.
* 0.3% of bins hold 5.7% of all the sequence, and they are the ones whose size already rules them out.
  metaBAT2's unbinnable remainder reaches 159 Mb; no prokaryotic genome is.

Both are only worth having if the answer is unchanged, which is what most of this file checks.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import checkm


def _write_bin(path, contigs, length=1000):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, "ACGT" * (length // 4)) for name in contigs)


# --- the plausibility filter -----------------------------------------------------------------


def test_a_normal_genome_sized_bin_is_scored():
    """3.8 Mb at a 50 kb N50 is the median passing bin; it must never be filtered."""
    assert checkm.implausible_reason(3_800_000, 50_000) is None


@pytest.mark.parametrize(
    "size,n50", [(990_000, 15_000), (6_660_000, 6_530_000), (1_380_000, 2_025)]
)
def test_the_most_extreme_bins_that_actually_pass_are_still_scored(size, n50):
    """The real limits of 4,002 passing bins: smallest 0.99 Mb, largest 6.66 Mb, worst N50 2,025 bp.

    The bounds are calibrated against those rather than guessed, so the extremes of the real
    distribution have to survive them or the filter is losing genomes.
    """
    assert checkm.implausible_reason(size, n50) is None


def test_a_bin_bigger_than_any_genome_is_not_scored():
    reason = checkm.implausible_reason(159_000_000, 50_000)
    assert reason and "larger than any prokaryotic genome" in reason


def test_a_bin_too_small_to_be_mostly_a_genome_is_not_scored():
    assert "too small" in (checkm.implausible_reason(90_000, 50_000) or "")


def test_a_fragmented_bin_of_plausible_size_is_still_scored():
    """There is no N50 floor, on purpose.

    One was tried at 2 kb. The worst N50 among bins that pass 70/10 is 2,025 bp, and raising the floor
    to 2,500 loses 36 real genomes, while the 2 kb version bought only 1.4 percentage points of sequence
    over the size bounds alone. Fragmentation is a property of the assembly, not evidence of a chimera.
    """
    assert checkm.implausible_reason(3_000_000, 900) is None
    assert not hasattr(checkm, "MIN_PLAUSIBLE_N50")


def test_the_bounds_are_loose_on_purpose():
    """20 Mb against a ~16 Mb largest known genome, and 200 kb against a 0.99 Mb smallest passer."""
    assert checkm.MAX_PLAUSIBLE_BIN_BP >= 16_000_000
    assert checkm.MIN_PLAUSIBLE_BIN_BP <= 900_000


# --- bin shape and identity ------------------------------------------------------------------


def test_bin_shape_measures_size_n50_and_gc(tmp_path):
    path = tmp_path / "b.fa"
    with open(path, "w") as fh:
        fh.write(">c1\n%s\n>c2\n%s\n" % ("G" * 1000, "A" * 500))
    total, n50, gc = checkm.bin_shape(str(path))
    assert total == 1500
    assert n50 == 1000, "the contig at which half the assembly is reached"
    assert gc == pytest.approx(1000 / 1500)


def test_the_same_contigs_in_a_different_order_are_the_same_bin(tmp_path):
    """A hybrid set can emit a bin's contigs in any order; it is still that bin."""
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    _write_bin(a, ["c1", "c2", "c3"])
    _write_bin(b, ["c3", "c1", "c2"])
    assert checkm.bin_identity(str(a)) == checkm.bin_identity(str(b))


def test_different_contigs_are_different_bins(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    _write_bin(a, ["c1", "c2"])
    _write_bin(b, ["c1", "c2", "c3"])
    assert checkm.bin_identity(str(a)) != checkm.bin_identity(str(b))


# --- the cache ------------------------------------------------------------------------------


def test_the_cache_survives_a_resume(tmp_path):
    path = str(tmp_path / "scores.json")
    first = checkm.BinScoreCache(path)
    first.rows["abc"] = ["bin_1", "95.0", "2.0", "0.5", "Bacteria", "50000", "3800000"]
    first.save()
    assert checkm.BinScoreCache(path).rows["abc"][1] == "95.0"


def test_a_corrupt_cache_is_ignored_rather_than_fatal(tmp_path):
    """The cache is an optimisation; a truncated one is a reason to re-measure, not to fail a run."""
    path = tmp_path / "scores.json"
    path.write_text("{not json")
    assert checkm.BinScoreCache(str(path)).rows == {}


def test_an_unwritable_cache_does_not_fail_the_run(tmp_path):
    cache = checkm.BinScoreCache(str(tmp_path / "nonexistent-dir" / "scores.json"))
    cache.rows["x"] = ["b", "1", "1", "1", "l", "1", "1"]
    cache.save()  # must not raise


def test_the_summary_accounts_for_every_bin():
    cache = checkm.BinScoreCache()
    cache.hits, cache.misses, cache.skipped = 3, 5, 2
    text = cache.summary()
    assert "10 bin(s)" in text and "5 scored" in text and "3 reused" in text and "2 skipped" in text


# --- the wrapper, end to end with a stubbed CheckM -------------------------------------------


def _fake_run_checkm(monkeypatch, scored):
    """Stub run_checkm so it writes a plausible .stats for whatever it is handed."""

    def fake(bins_dir, **kwargs):
        from metawrap2.io.seqio import fasta_files

        rows = []
        for path in fasta_files(bins_dir):
            name = os.path.splitext(os.path.basename(path))[0]
            scored.append(name)
            total, n50, gc = checkm.bin_shape(path)
            rows.append([name, "90.0", "3.0", "%.2f" % gc, "Bacteria", str(n50), str(total)])
        out = bins_dir + checkm.STATS_SUFFIX
        checkm.write_stats(rows, out)
        return out

    monkeypatch.setattr(checkm, "run_checkm", fake)
    return scored


def test_an_identical_bin_in_a_second_set_is_not_rescored(tmp_path, monkeypatch):
    scored = _fake_run_checkm(monkeypatch, [])
    cache = checkm.BinScoreCache()
    for name in ("setA", "setB"):
        d = tmp_path / name
        d.mkdir()
        _write_bin(d / "bin_1.fasta", ["c1", "c2"], length=2_000_000)
        _write_bin(d / "bin_2.fasta", ["c3"], length=2_000_000)
        checkm.run_checkm_cached(str(d), cache=cache, threads=1, mem_gb=8, env=None)

    assert cache.misses == 2, "two distinct bins scored once each"
    assert cache.hits == 2, "both reappear in the second set and are reused"
    assert len(scored) == 2, "CheckM saw each bin exactly once"


def test_the_stats_are_identical_whether_or_not_the_cache_was_used(tmp_path, monkeypatch):
    """The whole point: faster, not different."""
    _fake_run_checkm(monkeypatch, [])

    def build(root):
        d = tmp_path / root
        d.mkdir()
        _write_bin(d / "bin_1.fasta", ["c1", "c2"], length=2_000_000)
        _write_bin(d / "bin_2.fasta", ["c3"], length=3_000_000)
        return d

    cold = build("cold")
    checkm.run_checkm_cached(str(cold), cache=checkm.BinScoreCache(), threads=1, mem_gb=8, env=None)

    warm_cache = checkm.BinScoreCache()
    primed = build("primed")
    checkm.run_checkm_cached(str(primed), cache=warm_cache, threads=1, mem_gb=8, env=None)
    warm = build("warm")
    checkm.run_checkm_cached(str(warm), cache=warm_cache, threads=1, mem_gb=8, env=None)

    def rows(d):
        with open(str(d) + checkm.STATS_SUFFIX) as fh:
            return sorted(fh.read().splitlines())

    assert rows(cold) == rows(warm), "a cached run must report exactly what a cold one does"


def test_an_implausible_bin_is_recorded_as_failing_rather_than_omitted(tmp_path, monkeypatch):
    """Omitting it would leave downstream code guessing; 0/100 can never win a comparison."""
    scored = _fake_run_checkm(monkeypatch, [])
    d = tmp_path / "set"
    d.mkdir()
    _write_bin(d / "good.fasta", ["c1"], length=3_000_000)
    _write_bin(d / "huge.fasta", ["c2"], length=40_000_000)
    cache = checkm.BinScoreCache()
    stats = checkm.run_checkm_cached(str(d), cache=cache, threads=1, mem_gb=8, env=None)

    assert scored == ["good"], "the oversized bin never reached CheckM"
    assert cache.skipped == 1
    with open(stats) as fh:
        text = fh.read()
    assert "huge" in text and "skipped:" in text
    row = next(ln.split("\t") for ln in text.splitlines() if ln.startswith("huge"))
    assert float(row[1]) == 0.0 and float(row[2]) == 100.0


def test_the_filter_can_be_turned_off(tmp_path, monkeypatch):
    scored = _fake_run_checkm(monkeypatch, [])
    d = tmp_path / "set"
    d.mkdir()
    _write_bin(d / "huge.fasta", ["c2"], length=40_000_000)
    checkm.run_checkm_cached(
        str(d), cache=checkm.BinScoreCache(), prefilter=False, threads=1, mem_gb=8, env=None
    )
    assert scored == ["huge"]


def test_an_empty_bin_set_falls_through_to_the_plain_call(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(checkm, "run_checkm", lambda d, **k: calls.append(d) or (d + ".stats"))
    d = tmp_path / "empty"
    d.mkdir()
    checkm.run_checkm_cached(str(d), cache=checkm.BinScoreCache(), threads=1, mem_gb=8, env=None)
    assert calls == [str(d)], "nothing to optimise, so nothing clever should happen"
