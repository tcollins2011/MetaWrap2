"""The contig depth table, computed here instead of by jgi_summarize_bam_contig_depths.

These tests cover the arithmetic and the two defects found in comparison with jgi. Agreement
needs real BAMs and was verified separately over all 210,632 contigs of the CAMI2 GI assembly
(Pearson 0.999 on means, 0.995 on variances).
"""

from __future__ import annotations

import numpy as np
import pytest

from metawrap2 import depth


def _lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def test_a_contig_no_sample_covered_divides_by_its_length_not_by_one(tmp_path):
    """A contig with no reads in some samples is divided by its length, not by one.

    ``samtools depth -a`` emits nothing for a reference with no reads, so a divisor counted from
    emitted positions was 1, and a sample that did cover the contig reported the sum over the
    contig as its depth (1.6 million for a 4 Mbp contig).
    """
    widths = depth.trimmed_widths(np.array([4000000]), 75)
    total = np.array([[4000000.0]])  # a sample that does cover it, depth ~1
    total_sq = np.array([[4000000.0]])
    out = str(tmp_path / "d.txt")
    depth._write(out, ["big"], {"big": 4000000}, ["s1"], total, total_sq, widths)
    mean = float(_lines(out)[1].split("\t")[3])
    assert mean == pytest.approx(1.0, rel=1e-3)


def test_a_chatty_samtools_does_not_deadlock_the_depth_stream(tmp_path, monkeypatch):
    """More stderr than a pipe buffer holds, while stdout is still being read.

    Both stderr streams used to be read only after the depth stream was consumed. With one warning
    per read, the 64 KiB pipe buffer fills, samtools blocks on the write and the stdout read blocks
    forever. Uses a real subprocess, since the deadlock is in the pipes.
    """
    import sys

    noisy = tmp_path / "noisy_samtools.py"
    noisy.write_text(
        "import sys\n"
        "# 'view' writes only to stderr; 'depth' writes the rows and a flood of warnings.\n"
        "mode = 'depth' if 'depth' in sys.argv else 'view'\n"
        "sys.stderr.write(('[%s] a warning about a read\\n' % mode) * 4000)\n"
        "if mode == 'depth':\n"
        "    for pos in range(1, 101):\n"
        "        sys.stdout.write('c1\\t%d\\t3\\n' % pos)\n"
        "sys.stderr.flush()\n"
    )
    monkeypatch.setattr(depth, "_samtools", lambda env: [sys.executable, str(noisy)])

    total = np.zeros((1, 1))
    total_sq = np.zeros((1, 1))
    depth._accumulate_bam(
        str(tmp_path / "s.bam"),
        0,
        {"c1": 100},
        {"c1": 0},
        0,
        90.0,
        0,
        1,
        None,
        total,
        total_sq,
    )
    # It returned, and read every row past the point where the buffer would have filled.
    assert total[0, 0] == 300.0
