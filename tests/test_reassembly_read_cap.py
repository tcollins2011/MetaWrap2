"""Thinning each bin's recruited reads to a target coverage before SPAdes sees them.

Reads are recruited per bin, so an abundant organism can give SPAdes hundreds of times more
coverage than it needs. SPAdes time grows with read count, so a few abundant bins dominate the
step's runtime.
"""

from __future__ import annotations

import os

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


def test_both_mates_keep_exactly_the_same_reads(tmp_path):
    """Thinning keeps the same reads in both mate files, and leaves bins below the target alone."""
    bins, reads = tmp_path / "bins", tmp_path / "reads"
    bins.mkdir()
    reads.mkdir()
    # 10 kb bin, 10,000 pairs of 100 bp: 200x before thinning.
    _bin(str(bins), "bin.1", 10000)
    r1, r2 = _pair(str(reads), "bin.1.permissive", pairs=10000)
    # A bin already below the target (the common case) is left unchanged.
    _bin(str(bins), "bin.2", 100000)
    q1, q2 = _pair(str(reads), "bin.2.strict", pairs=50)
    before = (_text(q1), _text(q2))

    rb._cap_recruited_reads(_Args(max_coverage=50), str(reads), str(bins))
    assert _count_pairs(r1) == _count_pairs(r2)
    assert _read_names(r1) == _read_names(r2)
    # 50x of 10 kb is 500,000 bases at 200 per pair, so 2,500 pairs; each pair is its own coin flip.
    assert 2000 < _count_pairs(r1) < 3000, _count_pairs(r1)
    assert (_text(q1), _text(q2)) == before

    # Same seed, same reads kept on a re-run.
    again = tmp_path / "again"
    again.mkdir()
    a1, _ = _pair(str(again), "bin.1.permissive", pairs=10000)
    rb._cap_recruited_reads(_Args(max_coverage=50), str(again), str(bins))
    assert _read_names(a1) == _read_names(r1)
