"""An aligner index is reused only when it is complete, newer than its reference, and built by the
same command; otherwise it is rebuilt."""

import os
import time

from metawrap2 import aligners
from metawrap2.modules import _common

CMD = "minimap2 -x sr -d ref.fa.mmi ref.fa"


def _reference(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">c1\nACGT\n")
    return str(ref)


def _build(ref, aligner, command=CMD, files=None):
    for suffix in files if files is not None else aligners.INDEX_FILES[aligner]:
        with open(ref + suffix, "w") as fh:
            fh.write("index")
    with open(ref + aligners.INDEX_MARKER, "w") as fh:
        fh.write(command + "\n")


def _older(path):
    past = time.time() - 3600
    os.utime(path, (past, past))


def test_a_complete_current_index_is_reused(tmp_path):
    ref = _reference(tmp_path)
    _older(ref)
    _build(ref, "minimap2")
    assert _common._index_is_current(ref + ".mmi", ref, "minimap2", CMD)


def test_an_index_older_than_its_reference_is_rebuilt(tmp_path):
    ref = _reference(tmp_path)
    _build(ref, "minimap2")
    _older(ref + ".mmi")
    assert not _common._index_is_current(ref + ".mmi", ref, "minimap2", CMD)


def test_a_missing_or_empty_index_is_rebuilt(tmp_path):
    ref = _reference(tmp_path)
    assert not _common._index_is_current(ref + ".mmi", ref, "minimap2", CMD)
    _older(ref)
    _build(ref, "minimap2")
    with open(ref + ".mmi", "w"):
        pass
    assert not _common._index_is_current(ref + ".mmi", ref, "minimap2", CMD)


def test_bwa_needs_all_five_index_files(tmp_path):
    ref = _reference(tmp_path)
    _older(ref)
    _build(ref, "bwa", command="bwa index ref.fa", files=aligners.INDEX_FILES["bwa"][:4])
    assert not _common._index_is_current(ref, ref, "bwa", "bwa index ref.fa")
    _build(ref, "bwa", command="bwa index ref.fa")
    assert _common._index_is_current(ref, ref, "bwa", "bwa index ref.fa")


def test_an_index_without_a_marker_is_rebuilt(tmp_path):
    """No marker means the build did not finish, or another process is still writing it."""
    ref = _reference(tmp_path)
    _older(ref)
    _build(ref, "minimap2")
    os.remove(ref + aligners.INDEX_MARKER)
    assert not _common._index_is_current(ref + ".mmi", ref, "minimap2", CMD)


def test_an_index_built_with_another_preset_is_rebuilt(tmp_path):
    ref = _reference(tmp_path)
    _older(ref)
    _build(ref, "minimap2")
    assert not _common._index_is_current(
        ref + ".mmi", ref, "minimap2", "minimap2 -x map-ont -d ref.fa.mmi ref.fa"
    )


def test_every_aligner_lists_its_index_files():
    assert set(aligners.INDEX_FILES) == set(aligners.ALIGNERS)
    for aligner, suffixes in aligners.INDEX_FILES.items():
        assert suffixes and set(suffixes) <= set(aligners.INDEX_SUFFIXES), aligner


def test_the_build_path_writes_the_marker_and_the_second_call_reuses(tmp_path, monkeypatch):
    ref = _reference(tmp_path)
    _older(ref)
    built = []

    def fake_run(command, **_kw):
        built.append(command)
        with open(ref + ".mmi", "w") as fh:
            fh.write("index")

    monkeypatch.setattr(_common, "run", fake_run)
    _common._build_or_reuse_index(ref + ".mmi", ref, "minimap2", CMD, None)
    _common._build_or_reuse_index(ref + ".mmi", ref, "minimap2", CMD, None)
    assert built == [CMD]
    with open(ref + aligners.INDEX_MARKER) as fh:
        assert fh.read().strip() == CMD
