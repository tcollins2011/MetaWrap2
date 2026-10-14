"""A bin set may carry ``unbinned.fasta``, and nothing may treat it as a genome.

`binning` writes the contigs no CONCOCT cluster claimed to ``concoct_bins/unbinned.fasta`` on
purpose, so a user can see what was left over. Only the AMBER scorer knew that: `bin_refinement`
copied it in as a candidate genome, `annotate_bins` would have run prokka on it, `classify_bins`
would have asked taxator to place it, and `reassemble_bins` would have mapped reads to it and
reassembled it through SPAdes. On the CAMI2 GI run none of that happened, because CONCOCT's
leftovers came to 69.8 MB and refinement's size filter drops anything over 20 MB - the right
outcome, reached for a reason that has nothing to do with what the file is, and only while the
catch-all stays large.

So these tests are about the name, not the size: each one hands a module a bin folder whose unbinned
file is small enough to pass every other filter, and asserts it is still not treated as a bin.
"""

from __future__ import annotations

import os

from metawrap2.benchmark import assess
from metawrap2.constants import UNBINNED_NAME, is_unbinned_name
from metawrap2.io.seqio import fasta_files
from metawrap2.modules import bin_refinement
from metawrap2.modules._common import list_bins


def _bin_folder(tmp_path, unbinned_bytes=200_000):
    """A bin set of two real bins plus a catch-all small enough to pass every size filter."""
    d = tmp_path / "concoct_bins"
    d.mkdir()
    for n in (1, 2):
        (d / ("bin_%03d.fasta" % n)).write_text(">c%d\n%s\n" % (n, "ACGT" * 30_000))
    # 200 kb: above refinement's 50 kb floor and far below its 20 Mb ceiling, so if it is dropped it
    # is dropped for being the leftovers and not for its size.
    (d / UNBINNED_NAME).write_text(">leftover\n%s\n" % ("ACGT" * (unbinned_bytes // 4)))
    return d


def test_refinement_does_not_copy_it_in_as_a_bin(tmp_path, monkeypatch):
    """_copy_bins is the single point where a user's bin sets enter refinement."""
    monkeypatch.setattr(bin_refinement, "comm", lambda *a, **k: None)
    src = _bin_folder(tmp_path)
    dest = tmp_path / "binsB"
    kept = bin_refinement._copy_bins(str(src), str(dest))
    assert kept == 2, "the catch-all is not one of the bins"
    assert sorted(os.listdir(dest)) == ["bin_001.fasta", "bin_002.fasta"]

    # The shared listing is where the answer lives, and every other module reads bins through it.
    assert [os.path.basename(p) for p in fasta_files(str(src))] == [
        "bin_001.fasta",
        "bin_002.fasta",
    ]
    assert list_bins(str(src)) == ["bin_001.fasta", "bin_002.fasta"]
    # Gzipping a bin set happens after a binner wrote it, and a check against the full filename stops
    # holding there. A genome merely named for something containing the word is still a genome.
    assert is_unbinned_name("unbinned.fasta.gz") is True
    assert is_unbinned_name("unbinned_remainder.fasta") is False

    # The AMBER scorer would otherwise score the leftovers as a genome it was never offered.
    query = tmp_path / "q.binning"
    assess.bins_to_bioboxes(str(src), str(query))
    body = query.read_text()
    assert "leftover" not in body and "bin_001" in body and "bin_002" in body
