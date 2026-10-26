"""Tests for the jgi-style per-contig depth table: three leading columns, then two per sample.

metaBAT2, MaxBin2 and the downstream abundance steps read this layout, which is not documented
upstream. Guards a past bug where every column after the leading ones was read as a sample, which
doubled the sample count and used the variances as depths.
"""

from __future__ import annotations

import pytest

from metawrap2 import depth

#: Two samples in the real layout: a mean and a `-var` per sample. The variances differ from the
#: means so a test cannot pass by reading the wrong column.
TWO_SAMPLES = (
    "contigName\tcontigLen\ttotalAvgDepth\tsample_0\tsample_0-var\tsample_1\tsample_1-var\n"
    "c1\t5000\t8.5\t3.5\t91.1\t5\t92.2\n"
    "c2\t2000\t1.25\t0.25\t93.3\t1\t94.4\n"
)


def _table(tmp_path, text, name="depth.txt"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_three_leading_columns_then_two_per_sample(tmp_path):
    """The leading columns are ``contigName``, ``contigLen``, ``totalAvgDepth``; samples start at 3.

    A two-sample table has seven columns.
    """
    assert depth.LEADING_COLUMNS == ("contigName", "contigLen", "totalAvgDepth")
    samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert len(frame.columns) == 7
    assert samples == ["sample_0", "sample_1"]
    assert list(frame.columns[:3]) == list(depth.LEADING_COLUMNS)


def test_a_variance_column_is_never_mistaken_for_a_sample(tmp_path):
    """Columns ending in ``-var`` are variances, not samples.

    As a consequence, a sample cannot be named with a `-var` suffix.
    """
    samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert len(samples) * 2 + 3 == len(frame.columns)
    assert not any(name.endswith("-var") for name in samples)
    # Here the variances are an order of magnitude above the means.
    assert frame["sample_0"][0] == pytest.approx(3.5)
    assert frame["sample_0-var"][0] == pytest.approx(91.1)


def test_totalavgdepth_is_the_sum_over_samples_and_not_a_mean(tmp_path):
    """``totalAvgDepth`` is the per-contig sum over samples (3.5 + 5.0 = 8.5), despite its name."""
    _samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert frame["totalAvgDepth"][0] == pytest.approx(frame["sample_0"][0] + frame["sample_1"][0])
    assert frame["totalAvgDepth"][0] == pytest.approx(8.5)
