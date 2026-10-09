"""Per-bin abundance by aggregating the depth table the pipeline already built.

The other three quantifiers align the reads a second time to answer a question the binning run already
answered. They also answer it on a different scale: CoverM's per-contig coverage is twice a jgi-style
depth on strain-rich data, because jgi discards reads below 97% identity and CoverM does not. This
quantifier reads the table the binner itself read, so it is both free and on the binner's own scale.
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

    The length-weighted mean is "total bases over total length", which is what "how much of this genome
    was in the sample" means. A plain average answers a question about contigs instead, and the gap
    widens exactly when a bin is one long contig plus a tail of short ones - which is most bins.
    """
    args, table = _setup(tmp_path, {"binA": ["c1", "c2"]})
    quant_bins._quantify_depth(args, [], table, None, None)
    row = _read(table).loc["binA"]
    assert row["s1"] == pytest.approx((10 * 1000 + 6 * 3000) / 4000)  # 7.0
    assert row["s1"] != pytest.approx(8.0), "that would be the unweighted mean"
    assert row["length"] == 4000 and row["contigs"] == 2


def test_the_median_is_reported_alongside_for_the_robust_reading(tmp_path):
    args, table = _setup(tmp_path, {"binA": ["c1", "c2", "c3"]})
    quant_bins._quantify_depth(args, [], table, None, None)
    row = _read(table).loc["binA"]
    # depths 10, 6, 100 -> median 10, while the length-weighted mean is pulled down by the 3 kb contig
    assert row["s1-median"] == pytest.approx(10.0)
    assert row["s1"] == pytest.approx((10 * 1000 + 6 * 3000 + 100 * 500) / 4500)


def test_every_sample_column_is_carried_through_and_variance_columns_are_not(tmp_path):
    args, table = _setup(tmp_path, {"binA": ["c1", "c2"]})
    quant_bins._quantify_depth(args, [], table, None, None)
    cols = set(_read(table).columns)
    assert {"s1", "s2", "s1-median", "s2-median"} <= cols
    assert not any(c.endswith("-var") for c in cols)
    assert "totalAvgDepth" not in cols


def test_a_bin_whose_contigs_are_absent_is_reported_not_silently_zero(tmp_path, monkeypatch):
    """Bins and table from different assemblies is a mistake to hear about, not to average over."""
    warned = []
    monkeypatch.setattr(quant_bins, "warning", lambda m, *a, **k: warned.append(m))
    args, table = _setup(tmp_path, {"binA": ["c1"], "ghost": ["nope1", "nope2"]})
    quant_bins._quantify_depth(args, [], table, None, None)
    frame = _read(table)
    assert "ghost" not in frame.index and "binA" in frame.index
    assert any("none of its" in m for m in warned)


def test_a_table_that_is_not_a_depth_table_is_refused(tmp_path, monkeypatch):
    raised = []
    monkeypatch.setattr(
        quant_bins,
        "error",
        lambda m, *a, **k: raised.append(m) or (_ for _ in ()).throw(SystemExit),
    )
    bogus = tmp_path / "x.tsv"
    bogus.write_text("contig\tcoverage\nc1\t5\n")
    bins = tmp_path / "bins"
    bins.mkdir()
    args = _Args(str(bins), str(bogus), str(tmp_path))
    with pytest.raises(SystemExit):
        quant_bins._quantify_depth(args, [], str(tmp_path / "t.tsv"), None, None)
    assert raised and "not a contig depth table" in raised[0]


def test_it_is_registered_and_needs_no_reads():
    assert quant_bins.QUANTIFIERS["depth"] is quant_bins._quantify_depth
    # reads are nargs="*" so the quantifier can run without them at all
    args = quant_bins._parse_args(
        ["-b", "bins", "-o", "out", "--quantifier", "depth", "--depth", "d.txt"]
    )
    assert args.reads == [] and args.depth == "d.txt"


def test_omitting_the_table_is_an_error_that_says_where_to_find_one(tmp_path, monkeypatch):
    raised = []
    monkeypatch.setattr(
        quant_bins,
        "error",
        lambda m, *a, **k: raised.append(m) or (_ for _ in ()).throw(SystemExit),
    )
    bins = tmp_path / "bins"
    bins.mkdir()
    args = _Args(str(bins), None, str(tmp_path))
    with pytest.raises(SystemExit):
        quant_bins._quantify_depth(args, [], str(tmp_path / "t.tsv"), None, None)
    assert raised and "work_files" in raised[0]
