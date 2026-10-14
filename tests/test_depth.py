"""The contig depth table, computed here rather than by jgi_summarize_bam_contig_depths.

These cover the arithmetic and the two places it went wrong when checked against the tool it replaced.
The agreement itself is not testable without real BAMs - it was verified separately, over all 210,632
contigs of the CAMI2 GI assembly, at Pearson 0.999 on the means and 0.995 on the variances - so what is
pinned here is the behaviour that agreement depends on.
"""

from __future__ import annotations

import numpy as np
import pytest

from metawrap2 import depth


def _lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def test_a_contig_no_sample_covered_divides_by_its_length_not_by_one(tmp_path):
    """The bug the comparison against jgi caught, pinned.

    ``samtools depth -a`` says nothing at all about a reference with no reads, so a divisor counted
    from emitted positions was 1 for such a contig - and any sample that *did* cover it then reported
    the sum over the whole contig as its depth. One 4 Mbp contig came out at 1.6 million.
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

    Both stderrs used to be left unread until after the whole depth stream had been consumed. These
    are per-read tools over hundreds of millions of reads, so one warning per record fills the
    64 KiB pipe buffer; samtools then blocks on the write, stops producing rows, and the read of its
    stdout blocks forever. The run hung with no output and no timeout, which is why this uses a
    real subprocess rather than a fake: the deadlock is in the pipes, not in the logic.
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
    # It returned at all, and read every row rather than stopping where the buffer filled.
    assert total[0, 0] == 300.0
