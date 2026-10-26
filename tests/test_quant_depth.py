"""Per-bin abundance by aggregating the depth table the pipeline already built.

The other three quantifiers align the reads again. Their scale also differs: on strain-rich data
CoverM's per-contig coverage is twice a jgi-style depth, because jgi discards reads below 97%
identity and CoverM does not. This quantifier reads the table the binner read, so it needs no
alignment and uses the binner's scale.
"""

from __future__ import annotations

import pytest

from metawrap2.modules import quant_bins

DEPTH = (
    "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\ts2\ts2-var\n"
    "c1\t1000\t12.0\t10.0\t1.0\t2.0\t0.5\n"
    "c2\t3000\t8.0\t6.0\t1.0\t2.0\t0.5\n"
    "c3\t500\t100.0\t100.0\t1.0\t0.0\t0.0\n"
)


class _Args:
    def __init__(self, bins, depth, output, quantifier="depth"):
        self.bins = bins
        self.depth = depth
        self.output = output
        self.quantifier = quantifier
        self.assembly = None


def _setup(tmp_path, bins):
    table = tmp_path / "contig_depth.txt"
    table.write_text(DEPTH)
    bin_dir = tmp_path / "bins"
    bin_dir.mkdir()
    for name, contigs in bins.items():
        with open(bin_dir / ("%s.fa" % name), "w") as fh:
            fh.writelines(">%s\n%s\n" % (c, "ACGT" * 10) for c in contigs)
    out = tmp_path / "out"
    out.mkdir()
    return _Args(str(bin_dir), str(table), str(out)), str(out / "t.tsv")


def _read(path):
    import pandas as pd

    return pd.read_csv(path, sep="\t").set_index("bin")


def test_abundance_is_length_weighted_not_a_plain_average(tmp_path):
    """A bin of one 1 kb contig at 10x and one 3 kb at 6x is 7x, not 8x.

    The length-weighted mean (total bases over total length) measures how much of the genome was
    in the sample. A plain average describes contigs, and differs most for a bin of one long
    contig plus short ones, which is typical.
    """
    args, table = _setup(tmp_path, {"binA": ["c1", "c2"]})
    quant_bins._quantify_depth(args, [], table, None, None)
    row = _read(table).loc["binA"]
    assert row["s1"] == pytest.approx((10 * 1000 + 6 * 3000) / 4000)  # 7.0
    assert row["s1"] != pytest.approx(8.0), "that would be the unweighted mean"
    assert row["length"] == 4000 and row["contigs"] == 2
    # Every sample column is carried through, the variance columns are not.
    cols = set(_read(table).columns)
    assert {"s1", "s2", "s1-median", "s2-median"} <= cols
    assert not any(c.endswith("-var") for c in cols) and "totalAvgDepth" not in cols
