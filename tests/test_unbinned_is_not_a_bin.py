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

import pytest

from metawrap2.benchmark import assess
from metawrap2.constants import UNBINNED_NAME, is_unbinned_name
from metawrap2.io.seqio import fasta_files
from metawrap2.modules import bin_refinement
from metawrap2.modules._common import list_bins


@pytest.mark.parametrize(
    "name,expected",
    [
        ("unbinned.fasta", True),
        ("unbinned.fa", True),
        # Normalising the extension and gzipping a bin set are both things that happen to a folder
        # after a binner wrote it, and a check against the full filename stops holding at each.
        ("unbinned.fasta.gz", True),
        ("/some/where/unbinned.fasta", True),
        ("bin_001.fasta", False),
        # A genome may legitimately be named for something that merely contains the word.
        ("unbinned_remainder.fasta", False),
        ("not_unbinned.fasta", False),
        # Extensionless, which bin_label leaves alone rather than guessing at - so this matches. It
        # is the right way round: callers filter by extension before asking, so an extensionless
        # file never reaches here, and if one did it would be the same leftovers under another name.
        ("unbinned", True),
    ],
)
def test_is_unbinned_name(name, expected):
    assert is_unbinned_name(name) is expected


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


def test_fasta_files_skips_it(tmp_path):
    """The shared listing is where the answer lives, so it is the first thing to pin."""
    d = _bin_folder(tmp_path)
    assert [os.path.basename(p) for p in fasta_files(str(d))] == ["bin_001.fasta", "bin_002.fasta"]
    assert list_bins(str(d)) == ["bin_001.fasta", "bin_002.fasta"]


def test_refinement_does_not_copy_it_in_as_a_bin(tmp_path, monkeypatch):
    """_copy_bins is the single point where a user's bin sets enter refinement."""
    monkeypatch.setattr(bin_refinement, "comm", lambda *a, **k: None)
    src = _bin_folder(tmp_path)
    dest = tmp_path / "binsB"
    kept = bin_refinement._copy_bins(str(src), str(dest))
    assert kept == 2, "the catch-all is not one of the bins"
    assert sorted(os.listdir(dest)) == ["bin_001.fasta", "bin_002.fasta"]


def test_refinement_still_drops_oversized_bins(tmp_path, monkeypatch):
    """The size filter is not what excludes the catch-all, and it still has its own job.

    Written because the fix could have been made by tightening the size filter instead, which would
    have passed a test that only ever looked at the final count.
    """
    monkeypatch.setattr(bin_refinement, "comm", lambda *a, **k: None)
    src = tmp_path / "bins"
    src.mkdir()
    (src / "bin_001.fasta").write_text(">c1\n%s\n" % ("ACGT" * 30_000))
    (src / "bin_002.fasta").write_text(">c2\n%s\n" % ("ACGT" * 10))  # under 50 kb
    assert bin_refinement._copy_bins(str(src), str(tmp_path / "out")) == 1


def test_refinement_says_how_much_the_size_filter_took(tmp_path, monkeypatch):
    """Dropping half a bin set's sequence is worth a line, and the per-bin lines do not add up to one.

    The two ends are reported apart because they mean different things: over the ceiling is a
    catch-all that could never have been a MAG, under the floor is a fragment. Tracking the recall
    that refinement appears to lose back to this filter took a separate AMBER run per stage, which is
    what this line is here to save.
    """
    said: list[str] = []
    monkeypatch.setattr(bin_refinement, "comm", lambda m, *a, **k: said.append(m))
    monkeypatch.setattr(bin_refinement, "MAX_SIZE", 400_000)
    src = tmp_path / "bins"
    src.mkdir()
    (src / "keep.fasta").write_text(">c\n%s\n" % ("ACGT" * 30_000))  # ~120 kb
    (src / "catchall.fasta").write_text(">c\n%s\n" % ("ACGT" * 150_000))  # ~600 kb, over
    (src / "fragment.fasta").write_text(">c\nACGT\n")  # under the floor
    assert bin_refinement._copy_bins(str(src), str(tmp_path / "out")) == 1

    over = [m for m in said if "ceiling" in m]
    under = [m for m in said if "floor" in m]
    assert len(over) == 1 and "1 bin(s)" in over[0] and "0.6Mb" in over[0]
    assert "83." in over[0], "600kb of 720kb, so the share has to be the set's total: %r" % over[0]
    assert len(under) == 1 and "1 bin(s)" in under[0]


def test_nothing_dropped_means_nothing_said(tmp_path, monkeypatch):
    """A clean bin set must not grow two lines reporting 0%."""
    said: list[str] = []
    monkeypatch.setattr(bin_refinement, "comm", lambda m, *a, **k: said.append(m))
    src = tmp_path / "bins"
    src.mkdir()
    (src / "keep.fasta").write_text(">c\n%s\n" % ("ACGT" * 30_000))
    assert bin_refinement._copy_bins(str(src), str(tmp_path / "out")) == 1
    assert not [m for m in said if "ceiling" in m or "floor" in m]


def test_amber_query_excludes_it(tmp_path):
    """The scorer already skipped it; this pins that it still does after the shared listing."""
    d = _bin_folder(tmp_path)
    out = tmp_path / "q.binning"
    assess.bins_to_bioboxes(str(d), str(out))
    body = out.read_text()
    assert "leftover" not in body, "the leftovers would score as a genome AMBER was never offered"
    assert "bin_001" in body and "bin_002" in body


def test_a_folder_of_nothing_but_leftovers_is_not_a_bin_set(tmp_path):
    """Refusing is the right answer, rather than silently scoring one bin called "unbinned"."""
    d = tmp_path / "bins"
    d.mkdir()
    (d / UNBINNED_NAME).write_text(">leftover\nACGT\n")
    with pytest.raises(assess.AssessError):
        assess.bins_to_bioboxes(str(d), str(tmp_path / "q.binning"))
