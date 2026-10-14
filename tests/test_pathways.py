"""The two places `pathways` reads someone else's table, which is where it can silently produce nothing.

Module completeness and the hierarchy maths are tested in test_kegg.py and test_hierarchy.py. What is
left here is the plumbing: filtering KO calls by how they qualified, and matching an abundance
table's row labels against the KO table's genome names. Both fail quietly rather than loudly - the
wrong answer is an empty profile, not an exception - so both are pinned.
"""

from __future__ import annotations

from metawrap2.modules import pathways

# ── the abundance table ──────────────────────────────────────────────────────────────────────


def _write_abundance(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


def test_the_bin_extension_is_stripped_so_the_names_match_the_ko_table(tmp_path):
    """The KO table's genome names have no extension and quant_bins' row labels do.

    Requiring the user to reconcile them by hand would end in a silently empty community profile,
    which looks exactly like "this sample has no function".
    """
    for label in ("bin_001.fasta", "bin_001.fa", "bin_001.fna", "bin_001.fas", "bin_001"):
        path = _write_abundance(tmp_path / "a.tsv", "bins\ts1\n%s\t1.0\n" % label)
        assert "bin_001" in pathways._read_abundance(path)["s1"]


# ── the depth table, which is the route taken when no bins are passed ────────────────────────


def test_the_depth_tables_variance_columns_are_not_read_as_samples(tmp_path):
    """The layout is three leading columns and then *two* per sample: mean depth and its variance.

    Taking every column after the third as a sample doubles the cohort and gives each extra one a
    profile named `<sample>-var`, weighted by a variance - a plausible-looking table of nothing.
    """
    path = _write_abundance(
        tmp_path / "d.txt",
        "contigName\tcontigLen\ttotalAvgDepth\ts1\ts1-var\ts2\ts2-var\n"
        "contig_24\t158695\t20.0\t27.7\t35.6\t12.3\t4.1\n",
    )
    depth = pathways._read_contig_depth(path)
    assert set(depth) == {"s1", "s2"}
    assert depth["s1"]["contig_24"] == 27.7
    assert depth["s2"]["contig_24"] == 12.3
