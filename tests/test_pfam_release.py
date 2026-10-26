"""Tests for reading the Pfam release from the path shapes a caller can pass.

`annotate --pfam` printed "Pfam release unknown" although relnotes.txt said `RELEASE 38.2` on line
2. ``install-db pfam`` downloads relnotes.txt and records ``<dir>/Pfam-A.hmm`` in the config, so
the configured value is the HMM file; the reader streamed 2.2 GB of profiles looking for the
release line.

Results before the fix:

    config value (file)   -> unknown     (the only shape the tool uses)
    directory             -> 38.2
    relnotes.txt directly -> 38.2

Pfam renumbers and merges families between releases, so a domain table needs its release to be
reproduced or joined to another run's.
"""

from __future__ import annotations

import os

from metawrap2 import pfam

#: Pfam's wording from the top of the real file; the release is on the second line.
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
    """The configured HMM file path, which used to return "unknown".

    install-db writes PFAM_DB as ``<dir>/Pfam-A.hmm`` (a file-valued key, like GUNC's); the
    notes are in the same directory.
    """
    directory = _pfam_dir(tmp_path)
    assert pfam.release_from_notes(str(directory / "Pfam-A.hmm")) == "38.2"


def test_the_release_is_found_from_the_directory_and_from_the_notes_themselves(tmp_path):
    """A directory and relnotes.txt itself still work."""
    directory = _pfam_dir(tmp_path)
    assert pfam.release_from_notes(str(directory)) == "38.2"
    assert pfam.release_from_notes(str(directory / "relnotes.txt")) == "38.2"


def test_no_notes_beside_the_database_is_unknown_rather_than_a_guess(tmp_path):
    """The release is not guessed from a directory name; without notes it is "unknown"."""
    directory = _pfam_dir(tmp_path, with_notes=False)
    assert pfam.release_from_notes(str(directory / "Pfam-A.hmm")) == "unknown"
    assert pfam.release_from_notes(str(directory)) == "unknown"
    assert pfam.release_from_notes(str(tmp_path / "not-installed-here")) == "unknown"


def test_the_hmm_file_is_never_read_looking_for_a_release(tmp_path, monkeypatch):
    """The HMM file is not read when looking for the release.

    Checked by recording which files are opened, not by timing.
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
    """A large file is not read as release notes (`relnotes.txt` is a few kilobytes).

    The test file is sparse, so it takes no real disk space.
    """
    directory = tmp_path / "pfam-pretending"
    directory.mkdir()
    notes = directory / "relnotes.txt"
    with open(notes, "wb") as handle:
        handle.write(b"RELEASE 38.2\n")
        handle.truncate(pfam.NOTES_MAX_BYTES + 1)
    assert os.path.getsize(notes) > pfam.NOTES_MAX_BYTES
    assert pfam.release_from_notes(str(notes)) == "unknown"
