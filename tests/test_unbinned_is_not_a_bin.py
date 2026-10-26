"""Tests that ``unbinned.fasta`` in a bin set is never treated as a genome.

`binning` writes contigs not assigned to any CONCOCT cluster to ``concoct_bins/unbinned.fasta`` so
users can inspect them. Previously only the AMBER scorer skipped the file; `bin_refinement`,
`annotate_bins`, `classify_bins` and `reassemble_bins` would have treated it as a bin. On the CAMI2
GI run it was 69.8 MB and refinement's 20 MB size filter dropped it, but a smaller file would pass.

The fixtures use an unbinned file small enough to pass every size filter, so the tests check
exclusion by name.
"""

from __future__ import annotations

import os

from metawrap2.benchmark import assess
from metawrap2.constants import UNBINNED_NAME, is_unbinned_name
from metawrap2.io.seqio import fasta_files
from metawrap2.modules import bin_refinement
from metawrap2.modules._common import list_bins


def _bin_folder(tmp_path, unbinned_bytes=200_000):
    """A bin set of two bins plus an unbinned file small enough to pass every size filter."""
    d = tmp_path / "concoct_bins"
    d.mkdir()
    for n in (1, 2):
        (d / ("bin_%03d.fasta" % n)).write_text(">c%d\n%s\n" % (n, "ACGT" * 30_000))
    # 200 kb: between refinement's 50 kb minimum and 20 Mb maximum.
    (d / UNBINNED_NAME).write_text(">leftover\n%s\n" % ("ACGT" * (unbinned_bytes // 4)))
    return d


def test_refinement_does_not_copy_it_in_as_a_bin(tmp_path, monkeypatch):
    """_copy_bins, through which bin sets enter refinement, skips the unbinned file."""
    monkeypatch.setattr(bin_refinement, "comm", lambda *a, **k: None)
    src = _bin_folder(tmp_path)
    dest = tmp_path / "binsB"
    kept = bin_refinement._copy_bins(str(src), str(dest))
    assert kept == 2, "the catch-all is not one of the bins"
    assert sorted(os.listdir(dest)) == ["bin_001.fasta", "bin_002.fasta"]

    # Other modules list bins through these shared helpers.
    assert [os.path.basename(p) for p in fasta_files(str(src))] == [
        "bin_001.fasta",
        "bin_002.fasta",
    ]
    assert list_bins(str(src)) == ["bin_001.fasta", "bin_002.fasta"]
    # A gzipped unbinned file is matched; a bin whose name only contains "unbinned" is not.
    assert is_unbinned_name("unbinned.fasta.gz") is True
    assert is_unbinned_name("unbinned_remainder.fasta") is False

    # The AMBER conversion also excludes it.
    query = tmp_path / "q.binning"
    assess.bins_to_bioboxes(str(src), str(query))
    body = query.read_text()
    assert "leftover" not in body and "bin_001" in body and "bin_002" in body
