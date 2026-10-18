"""The Pfam release, read from the four path shapes a caller can hand over.

`annotate --pfam` printed "Pfam release unknown" on every run it ever made, beside a database whose
relnotes.txt said `RELEASE 38.2` on line 2. The reason is a shape mismatch, not a parsing bug:
``install-db pfam`` downloads relnotes.txt for exactly this purpose and then records
``<dir>/Pfam-A.hmm`` in the config, so the configured value is a *file that is not the notes* -
``isdir`` is false, ``isfile`` is true, and the reader streamed 2.2 GB of profiles to the end
looking for a line that was in the file next door.

Three measurements from the report, which these tests now hold in place:

    config value (file)   -> unknown     (the only shape the tool uses)
    directory             -> 38.2
    relnotes.txt directly -> 38.2

The release is not cosmetic. Pfam renumbers and merges families between releases, so a domain
table that does not name its release cannot be reproduced or joined to another run's.
"""

from __future__ import annotations

import os

from metawrap2 import pfam

#: Pfam's own wording, from the top of the real file: the release is on the second line.
RELNOTES = """\
                                 Pfam release notes
                                     RELEASE 38.2

Pfam 38.2 contains 24,076 families.
"""


def _pfam_dir(tmp_path, with_notes: bool = True):
    """A Pfam install as `install-db pfam` leaves it: the HMMs, and the notes beside them."""
    directory = tmp_path / "pfam-38.2"
    directory.mkdir()
    (directory / "Pfam-A.hmm").write_text("HMMER3/f [3.4]\nNAME  1-cysPrx_C\n")
    if with_notes:
        (directory / "relnotes.txt").write_text(RELNOTES)
    return directory


def test_the_release_is_found_from_the_configured_hmm_file(tmp_path):
    """The shape the tool actually passes, and the one that returned "unknown".

    PFAM_DB is ``<dir>/Pfam-A.hmm`` because GUNC-style file-valued keys are what install-db
    writes; the notes are its sibling. Looking beside the file is what makes the provenance the
    installer went and downloaded reach the output table.
    """
    directory = _pfam_dir(tmp_path)
    assert pfam.release_from_notes(str(directory / "Pfam-A.hmm")) == "38.2"


def test_the_release_is_found_from_the_directory_and_from_the_notes_themselves(tmp_path):
    """The two shapes that already worked. They keep working - this is a widening, not a swap."""
    directory = _pfam_dir(tmp_path)
    assert pfam.release_from_notes(str(directory)) == "38.2"
    assert pfam.release_from_notes(str(directory / "relnotes.txt")) == "38.2"


def test_no_notes_beside_the_database_is_unknown_rather_than_a_guess(tmp_path):
    """A directory name is whatever the installer called it; an unread release is "unknown".

    Labelling a domain table with a release nobody read is worse than labelling it honestly,
    because the wrong release number is indistinguishable from the right one downstream.
    """
    directory = _pfam_dir(tmp_path, with_notes=False)
    assert pfam.release_from_notes(str(directory / "Pfam-A.hmm")) == "unknown"
    assert pfam.release_from_notes(str(directory)) == "unknown"
    assert pfam.release_from_notes(str(tmp_path / "not-installed-here")) == "unknown"


def test_the_hmm_file_is_never_read_looking_for_a_release(tmp_path, monkeypatch):
    """The cost of the old behaviour, pinned: 2.2 GB of profiles read to print "unknown".

    Checked by watching which files are opened rather than by timing, so this says what it means
    on a cold filesystem as well as a warm one.
    """
    directory = _pfam_dir(tmp_path)
    opened = []
    real_open = open

    def tracking_open(file, *args, **kwargs):
        opened.append(str(file))
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr("builtins.open", tracking_open)
    assert pfam.release_from_notes(str(directory / "Pfam-A.hmm")) == "38.2"
    assert opened == [str(directory / "relnotes.txt")], "opened %s" % opened


def test_a_file_too_large_to_be_release_notes_is_not_streamed(tmp_path):
    """The guard: `relnotes.txt` is a couple of kilobytes, so a huge one is not release notes.

    Without it, any caller that can get a large file into this function - a notes path pointed at
    the wrong file, a future key_file change - pays a full read to learn nothing. Sparse, so the
    test costs a size and not a gigabyte.
    """
    directory = tmp_path / "pfam-pretending"
    directory.mkdir()
    notes = directory / "relnotes.txt"
    with open(notes, "wb") as handle:
        handle.write(b"RELEASE 38.2\n")
        handle.truncate(pfam.NOTES_MAX_BYTES + 1)
    assert os.path.getsize(notes) > pfam.NOTES_MAX_BYTES
    assert pfam.release_from_notes(str(notes)) == "unknown"
