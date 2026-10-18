"""The jgi-style per-contig depth table: three leading columns, then *two* per sample.

This is the layout metaBAT2, MaxBin2 and every downstream abundance step read, and its shape is
documented nowhere. The defect it caused: a reader that treated every column after the leading ones
as a sample **doubled the cohort**, inventing one extra sample per real one and giving each of them
a nonsense depth profile taken from the variances. Nothing errored - the table simply described
twice as many samples as were sequenced.
"""

from __future__ import annotations

import pytest

from metawrap2 import depth

#: Two samples in the real layout: a mean and a `-var` per sample, with the variances deliberately
#: unlike the means so a test cannot pass by reading either one.
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
    """``contigName``, ``contigLen``, ``totalAvgDepth`` - and the sample columns start at index 3.

    A two-sample table is seven columns wide. Counting all four trailing columns as samples is what
    doubled the cohort; counting from index 2 instead would make `totalAvgDepth` a sample whose
    depth is the sum of all the others.
    """
    assert depth.LEADING_COLUMNS == ("contigName", "contigLen", "totalAvgDepth")
    samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert len(frame.columns) == 7
    assert samples == ["sample_0", "sample_1"]
    assert list(frame.columns[:3]) == list(depth.LEADING_COLUMNS)


def test_a_variance_column_is_never_mistaken_for_a_sample(tmp_path):
    """The ``-var`` suffix is the only thing distinguishing the two, and it is a suffix not a prefix.

    A sample genuinely named `sample_0-var` cannot exist here, which is the trade; what matters is
    that a four-sample cohort never reports eight samples, each second one carrying variances in a
    column every downstream step reads as a depth.
    """
    samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert len(samples) * 2 + 3 == len(frame.columns)
    assert not any(name.endswith("-var") for name in samples)
    # The values behind the two kinds of column are not interchangeable: here the variances are an
    # order of magnitude above the means, so reading one as the other is not a small error.
    assert frame["sample_0"][0] == pytest.approx(3.5)
    assert frame["sample_0-var"][0] == pytest.approx(91.1)


def test_totalavgdepth_is_the_sum_over_samples_and_not_a_mean(tmp_path):
    """3.5 + 5.0 = 8.5. A per-contig *total*, despite "Avg" in the name.

    Treated as a mean it is half the depth on a two-sample cohort and a tenth on a ten-sample one,
    so a cohort-size-dependent error in any coverage derived from this column.
    """
    _samples, frame = depth.read_depth_table(_table(tmp_path, TWO_SAMPLES))
    assert frame["totalAvgDepth"][0] == pytest.approx(frame["sample_0"][0] + frame["sample_1"][0])
    assert frame["totalAvgDepth"][0] == pytest.approx(8.5)
