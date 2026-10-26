"""CheckM caching and size pre-filtering in bin_refinement leave the reported numbers unchanged.

bin_refinement scores 2^n-1 hybrid bin sets, and CheckM is nearly all of the cost. Two measured
redundancies:

* 32% of the bins across the seven sets of a three-binner run are the same bin (same contigs)
  scored again.
* 0.3% of bins hold 5.7% of all sequence and are already too large to be a genome; metaBAT2's
  unbinnable remainder reaches 159 Mb.

Most of this file checks that the results are the same with both optimisations.
"""

from __future__ import annotations

import os

from metawrap2 import checkm


def _write_bin(path, contigs, length=1000):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, "ACGT" * (length // 4)) for name in contigs)


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


def test_the_stats_are_identical_whether_or_not_the_cache_was_used(tmp_path, monkeypatch):
    """Stats are identical with and without the cache."""
    scored = _fake_run_checkm(monkeypatch, [])

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
    # The warm run was served from the cache: CheckM saw each distinct bin once per cache.
    assert warm_cache.hits == 2 and warm_cache.misses == 2
    assert len(scored) == 4, "two for the cold run, two to prime the warm cache, none for warm"


def test_an_implausible_bin_is_recorded_as_failing_rather_than_omitted(tmp_path, monkeypatch):
    """An out-of-bounds bin is recorded as 0/100 rather than omitted, so it never wins."""
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

    # The bounds are calibrated against 4,002 passing bins (smallest 0.99 Mb, largest 6.66 Mb,
    # worst N50 2,025 bp), so those extremes must pass. There is no N50 floor, because
    # fragmentation does not indicate a chimera.
    for size, n50 in (
        (990_000, 15_000),
        (6_660_000, 6_530_000),
        (1_380_000, 2_025),
        (3_000_000, 900),
    ):
        assert checkm.implausible_reason(size, n50) is None, (size, n50)
    assert "too small" in (checkm.implausible_reason(90_000, 50_000) or "")
