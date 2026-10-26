"""`pathways` input handling for external tables.

Module completeness and hierarchy calculations are tested in test_kegg.py and test_hierarchy.py.
These tests cover filtering KO calls by how they qualified and matching abundance table row labels
to the KO table's genome names. A failure in either gives an empty profile, not an exception.
"""

from __future__ import annotations

from metawrap2.modules import pathways

# ── the abundance table ──────────────────────────────────────────────────────────────────────


def _write_abundance(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


def test_the_bin_extension_is_stripped_so_the_names_match_the_ko_table(tmp_path):
    """quant_bins' row labels have an extension that is stripped to match the KO table.

    Unmatched names would give an empty community profile.
    """
    for label in ("bin_001.fasta", "bin_001.fa", "bin_001.fna", "bin_001.fas", "bin_001"):
        path = _write_abundance(tmp_path / "a.tsv", "bins\ts1\n%s\t1.0\n" % label)
        assert "bin_001" in pathways._read_abundance(path)["s1"]


# ── the depth table (used when no bins are passed) ───────────────────────────────────────────


def test_the_depth_tables_variance_columns_are_not_read_as_samples(tmp_path):
    """The depth table has three leading columns, then two per sample: mean depth and variance.

    `<sample>-var` columns must not be read as extra samples.
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


# ── the hierarchy location ───────────────────────────────────────────────────────────────────


def test_the_hierarchy_can_be_given_on_the_command_line():
    """The KEGG BRITE hierarchy can be given with a flag as well as the config key.

    The hierarchy is downloaded, not shipped, and a revised one changes the percentages, so a run
    must be able to name its snapshot without editing metawrap2.toml.
    """
    args = pathways._parse_args(["-o", "out", "-k", "kos.tsv", "--brite-db", "/db/KEGG_BRITE"])
    assert args.brite_db == "/db/KEGG_BRITE"
    assert pathways._parse_args(["-o", "out", "-k", "kos.tsv"]).brite_db is None


def test_a_missing_hierarchy_names_every_way_of_supplying_one():
    """The error for a missing hierarchy lists all three ways to supply one."""
    import inspect

    source = inspect.getsource(pathways.main)
    assert "install-db kegg" in source
    assert "KEGG_BRITE" in source
    assert "--brite-db" in source
