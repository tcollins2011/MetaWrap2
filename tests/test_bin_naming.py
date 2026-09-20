"""Tests for the one bin-naming convention: bin_001.fasta, for every binner.

Each binner names its output its own way - metaBAT2 counts from 1, MaxBin2 writes zero-padded
``.001.fasta``, CONCOCT emits its own cluster ids - so without normalising, the same number means a
different thing in each folder. The properties worth pinning are that the convention is uniform,
that it counts from
1, and that padding makes filename order and bin order the same thing (which matters because
refined bins are numbered in sorted filename order).
"""

from __future__ import annotations

import os

import pytest

from metawrap2.constants import (
    BIN_EXTENSION,
    BIN_PREFIX,
    UNBINNED_NAME,
    bin_filename,
    bin_stem,
)
from metawrap2.modules import binning

# --- the naming helpers -------------------------------------------------------------------


def test_bins_are_named_from_one_and_zero_padded():
    assert bin_filename(1) == "bin_001.fasta"
    assert bin_filename(2) == "bin_002.fasta"
    assert bin_filename(42) == "bin_042.fasta"
    assert bin_filename(999) == "bin_999.fasta"


def test_the_stem_is_what_a_stats_row_holds():
    assert bin_stem(1) == "bin_001"
    assert bin_filename(1) == bin_stem(1) + BIN_EXTENSION


def test_padding_makes_filename_order_match_bin_order():
    """The reason for padding: unpadded, sorting puts bin.10 between bin.1 and bin.2."""
    names = [bin_filename(n) for n in range(1, 13)]
    assert sorted(names) == names
    unpadded = ["bin.%d.fa" % n for n in range(1, 13)]
    assert sorted(unpadded) != unpadded  # what it used to do


def test_numbers_past_the_padding_width_stay_unique():
    assert bin_filename(1000) == "bin_1000.fasta"
    assert bin_filename(1000) != bin_filename(100)


def test_the_extension_names_the_format_in_full():
    assert BIN_EXTENSION == ".fasta"
    assert UNBINNED_NAME == "unbinned.fasta"
    assert BIN_PREFIX == "bin_"


# --- renumbering --------------------------------------------------------------------------


def make_bins(directory, names):
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        (directory / name).write_text(">c\nACGT\n")
    return directory


def test_renumber_gives_every_bin_the_convention(tmp_path):
    bins = make_bins(tmp_path / "bins", ["bin.1.fa", "bin.2.fa", "bin.3.fa"])
    assert binning.renumber_bins(str(bins)) == 3
    assert sorted(os.listdir(bins)) == ["bin_001.fasta", "bin_002.fasta", "bin_003.fasta"]


def test_renumber_orders_numerically_not_lexicographically(tmp_path):
    """metaBAT2's bin.2 must stay ahead of its bin.10, so bin contents keep their order."""
    bins = tmp_path / "bins"
    make_bins(bins, [])
    for number in (1, 2, 10):
        (bins / ("bin.%d.fa" % number)).write_text(">contig_from_%d\nACGT\n" % number)
    binning.renumber_bins(str(bins))
    assert (bins / "bin_001.fasta").read_text().startswith(">contig_from_1")
    assert (bins / "bin_002.fasta").read_text().startswith(">contig_from_2")
    assert (bins / "bin_003.fasta").read_text().startswith(">contig_from_10")


def test_renumber_leaves_non_fasta_files_alone(tmp_path):
    bins = make_bins(tmp_path / "bins", ["bin.1.fa"])
    (bins / "checkm.log").write_text("...")
    (bins / "depths.tsv").write_text("...")
    assert binning.renumber_bins(str(bins)) == 1
    assert sorted(os.listdir(bins)) == ["bin_001.fasta", "checkm.log", "depths.tsv"]


def test_renumber_handles_names_that_already_match(tmp_path):
    """Idempotent: running it twice must not shuffle or lose anything."""
    bins = make_bins(tmp_path / "bins", ["bin_001.fasta", "bin_002.fasta"])
    (bins / "bin_001.fasta").write_text(">first\nACGT\n")
    (bins / "bin_002.fasta").write_text(">second\nACGT\n")
    binning.renumber_bins(str(bins))
    assert (bins / "bin_001.fasta").read_text().startswith(">first")
    assert (bins / "bin_002.fasta").read_text().startswith(">second")


def test_renumber_does_not_lose_a_bin_to_an_overlapping_name(tmp_path):
    """The reason renaming is staged: source and target namespaces overlap.

    A directory holding both ``bin.1.fa`` and ``bin_001.fasta`` would, renamed in place, have one
    silently overwrite the other - losing a genome.
    """
    bins = tmp_path / "bins"
    make_bins(bins, [])
    (bins / "bin.1.fa").write_text(">from_dot_one\nACGT\n")
    (bins / "bin_001.fasta").write_text(">from_underscore_one\nACGT\n")
    assert binning.renumber_bins(str(bins)) == 2
    contents = sorted(p.read_text() for p in bins.iterdir())
    assert len(contents) == 2
    assert ">from_dot_one\nACGT\n" in contents
    assert ">from_underscore_one\nACGT\n" in contents


def test_renumber_mixed_extensions(tmp_path):
    bins = make_bins(tmp_path / "bins", ["a.fa", "b.fasta", "c.fna", "d.fas"])
    assert binning.renumber_bins(str(bins)) == 4
    assert all(n.endswith(BIN_EXTENSION) for n in os.listdir(bins))


def test_renumber_an_empty_directory(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    assert binning.renumber_bins(str(bins)) == 0


# --- CONCOCT's cluster ids ----------------------------------------------------------------


def test_concoct_cluster_zero_becomes_bin_one(tmp_path):
    """CONCOCT numbers clusters from 0; bins are numbered from 1 like every other binner's."""
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>c2\nTTTT\n>c3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\nc2,1\nc3,5\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert (out / "bin_001.fasta").read_text() == ">c1\nACGT\n"
    assert (out / "bin_002.fasta").read_text() == ">c2\nTTTT\n"
    # Cluster ids are not necessarily contiguous, and are not renumbered here - the number still
    # identifies which CONCOCT cluster the bin came from.
    assert (out / "bin_006.fasta").read_text() == ">c3\nGGGG\n"


def test_concoct_unbinned_contigs_go_to_the_unbinned_file(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>c2\nTTTT\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert (out / UNBINNED_NAME).read_text() == ">c2\nTTTT\n"


# --- what the tools are told --------------------------------------------------------------


@pytest.mark.parametrize(
    "template",
    [
        "metawrap2.checkm.CHECKM2_PREDICT",
        "metawrap2.modules.reassemble_bins.CHECKM_QA_PLOT",
        "metawrap2.modules.classify_bins.GTDBTK",
    ],
)
def test_no_tool_is_told_the_old_extension(template):
    """Every -x flag must follow BIN_EXTENSION; a stale `-x fa` silently scores zero bins."""
    import importlib

    module_path, name = template.rsplit(".", 1)
    try:
        value = getattr(importlib.import_module(module_path), name)
    except AttributeError:
        pytest.skip("%s does not exist" % template)
    assert "-x fa " not in value and not value.endswith("-x fa")
    if "-x " in value:
        assert "-x fasta" in value or "{" in value.split("-x ", 1)[1][:8]
