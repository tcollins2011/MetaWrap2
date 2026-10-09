"""Thinning each bin's recruited reads to a target coverage before SPAdes sees them.

Read recruitment is per bin, so an abundant organism hands SPAdes hundreds of times the coverage
short-read assembly can use - and SPAdes time grows with reads, not with genome size, so a few
high-abundance outliers set the pace of the whole step.
"""

from __future__ import annotations

import os
import random

import pytest

from metawrap2.modules import reassemble_bins as rb


class _Args:
    def __init__(self, max_coverage=rb.DEFAULT_MAX_COVERAGE):
        self.max_coverage = max_coverage


def _pair(directory, name, pairs, read_len=100):
    """Write ``<name>_1.fastq``/``_2.fastq`` with *pairs* reads of *read_len* bases each."""
    paths = []
    for mate in ("1", "2"):
        path = os.path.join(directory, "%s_%s.fastq" % (name, mate))
        with open(path, "w") as fh:
            fh.writelines(
                "@r%d/%s\n%s\n+\n%s\n" % (i, mate, "A" * read_len, "I" * read_len)
                for i in range(pairs)
            )
        paths.append(path)
    return paths


def _bin(directory, name, length):
    with open(os.path.join(directory, "%s.fa" % name), "w") as fh:
        fh.write(">c1\n%s\n" % ("ACGT" * (length // 4)))


def _text(path):
    with open(path) as fh:
        return fh.read()


def _read_names(path):
    with open(path) as fh:
        return [ln.split("/")[0] for ln in fh if ln.startswith("@r")]


def _count_pairs(path):
    with open(path) as fh:
        return sum(1 for _ in fh) // 4


def test_the_default_is_fifty_times_and_zero_turns_it_off():
    assert rb.DEFAULT_MAX_COVERAGE == 50.0
    # Fixed, because a subsample that changed between runs would make reassembly irreproducible for
    # the same reason an unseeded binner is.
    assert rb.SUBSAMPLE_SEED == 1337


def test_bin_lengths_are_the_bins_own_assembled_size(tmp_path):
    """Against the bin, not an assumed genome size: 1.2 Mbp and 6 Mbp want different read counts."""
    bins = tmp_path / "bins"
    bins.mkdir()
    _bin(str(bins), "bin.1", 1000)
    _bin(str(bins), "bin.2", 4000)
    assert rb._bin_lengths(str(bins)) == {"bin.1": 1000, "bin.2": 4000}


def test_a_bin_already_below_the_target_is_not_touched(tmp_path):
    """The common case. The thinning exists for outliers, so it must cost nothing for everyone else."""
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    _bin(str(bins), "bin.1", 100000)
    r1, r2 = _pair(str(reads), "bin.1.strict", pairs=50)
    before = (_text(r1), _text(r2))

    rb._cap_recruited_reads(_Args(), str(reads), str(bins))
    assert (_text(r1), _text(r2)) == before


def test_an_abundant_bin_is_thinned_towards_the_target(tmp_path):
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    # 10 kb bin, 10,000 pairs of 100 bp: 2,000,000 bases over 10,000 bp, so 200x before thinning.
    _bin(str(bins), "bin.1", 10000)
    r1, _r2 = _pair(str(reads), "bin.1.strict", pairs=10000)

    rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
    kept = _count_pairs(r1)
    # 50x of a 10 kb bin is 500,000 bases, and each pair carries 200, so the target is 2,500 pairs.
    # Each pair is kept by its own coin flip, so the landing point is approximate by design - but it
    # has to be in the right neighbourhood rather than merely smaller.
    assert 2000 < kept < 3000, kept


def test_both_mates_keep_exactly_the_same_reads(tmp_path):
    """SPAdes is handed a paired library; a file where one mate survived is worse than a smaller one."""
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    _bin(str(bins), "bin.1", 10000)
    r1, r2 = _pair(str(reads), "bin.1.permissive", pairs=4000)

    rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
    assert _count_pairs(r1) == _count_pairs(r2)
    assert _read_names(r1) == _read_names(r2)


def test_max_coverage_zero_keeps_every_read(tmp_path):
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    _bin(str(bins), "bin.1", 1000)
    r1, _r2 = _pair(str(reads), "bin.1.strict", pairs=5000)

    rb._cap_recruited_reads(_Args(max_coverage=0), str(reads), str(bins))
    assert _count_pairs(r1) == 5000


def test_the_cut_off_is_not_part_of_the_bins_identity(tmp_path):
    """``bin.1.strict`` and ``bin.1.permissive`` are both ``bin.1``, measured against its length."""
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    _bin(str(bins), "bin.1", 10000)
    strict, _ = _pair(str(reads), "bin.1.strict", pairs=4000)
    permissive, _ = _pair(str(reads), "bin.1.permissive", pairs=4000)

    rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
    assert _count_pairs(strict) < 4000
    assert _count_pairs(permissive) < 4000


def test_a_read_set_with_no_matching_bin_is_left_alone(tmp_path):
    """No length to measure against means no decision to make - not a guess."""
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    _bin(str(bins), "bin.1", 10000)
    orphan, _ = _pair(str(reads), "bin.99.strict", pairs=4000)

    rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
    assert _count_pairs(orphan) == 4000


def test_thinning_is_reproducible(tmp_path):
    """Same seed, same reads kept: a re-run must not produce a different assembly input."""
    kept = []
    for attempt in ("a", "b"):
        bins, reads = tmp_path / ("bins_%s" % attempt), tmp_path / ("reads_%s" % attempt)
        bins.mkdir()
        reads.mkdir()
        _bin(str(bins), "bin.1", 10000)
        r1, _ = _pair(str(reads), "bin.1.strict", pairs=4000)
        rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
        kept.append(_read_names(r1))
    assert kept[0] == kept[1]


def test_the_flag_rejects_a_negative_target():
    with pytest.raises(SystemExit):
        rb._parse_args(
            ["-b", "bins", "-o", "out", "-1", "r1.fq", "-2", "r2.fq", "--max-coverage", "-5"]
        )


def test_thin_pair_keeps_everything_when_asked_to(tmp_path):
    """The fraction is a probability, so keep=1.0 has to be exact rather than nearly all."""
    r1, r2 = _pair(str(tmp_path), "x", pairs=100)
    kept = rb._thin_pair(r1, r2, 1.0, random.Random(0))
    assert kept == 100 and _count_pairs(r1) == 100 and _count_pairs(r2) == 100
