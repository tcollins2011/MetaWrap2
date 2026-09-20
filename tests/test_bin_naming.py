"""Tests for the one bin-naming convention: bin_001.fasta, for every binner.

Each binner names its output its own way - metaBAT2 counts from 1, MaxBin2 writes zero-padded
``.001.fasta``, CONCOCT emits its own cluster ids - so without normalising, the same number means a
different thing in each folder. The properties worth pinning are that the convention is uniform,
that it counts from 1, and that padding makes filename order and bin order the same thing (which
matters because refined bins are numbered in sorted filename order).
"""

from __future__ import annotations

import pytest

from metawrap2.constants import bin_filename
from metawrap2.modules import binning

# --- the naming helpers -------------------------------------------------------------------


def test_bins_are_named_from_one_and_zero_padded():
    assert bin_filename(1) == "bin_001.fasta"
    assert bin_filename(2) == "bin_002.fasta"
    assert bin_filename(42) == "bin_042.fasta"
    assert bin_filename(999) == "bin_999.fasta"


# --- renumbering --------------------------------------------------------------------------


def make_bins(directory, names):
    directory.mkdir(parents=True, exist_ok=True)
    for name in names:
        (directory / name).write_text(">c\nACGT\n")
    return directory


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
