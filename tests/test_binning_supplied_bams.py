"""Tests for `--bams DIR`: binning from existing alignments.

SemiBin2, VAMB, COMEBin and CONCOCT read alignments themselves and refuse `--depth <table>`.
`--bams` lets them use BAMs from `coverage --keep-bams`, so alignment and binning can run on
different machines.

A BAM aligned against a concatenated catalogue is well-formed, and a binner that matches no contigs
produces no bins without an error, so supplied BAM headers are checked against the assembly first.
"""

import os

import pytest

from metawrap2.modules import binning


def _assembly(tmp_path, names=("k141_1", "k141_2")):
    path = tmp_path / "assembly.fa"
    path.write_text("".join(">%s\nACGTACGTACGT\n" % name for name in names))
    return str(path)


def _bam_dir(tmp_path, names=("sample_0.bam", "sample_1.bam"), index=False):
    bams = tmp_path / "bams"
    bams.mkdir()
    for name in names:
        (bams / name).write_text("bam")
        if index:
            (bams / (name + ".bai")).write_text("bai")
    return str(bams)


def _header(names):
    return "@HD\tVN:1.6\n" + "".join("@SQ\tSN:%s\tLN:12\n" % n for n in names)


def test_reference_names_come_off_the_header(monkeypatch):
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: _header(["k141_1", "k141_2"]))
    assert binning.bam_reference_names("x.bam", "e") == ["k141_1", "k141_2"]


def test_bams_from_the_right_assembly_pass(tmp_path, monkeypatch):
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: _header(["k141_1", "k141_2"]))
    overlap = binning.check_bam_references(["s.bam"], _assembly(tmp_path), "e")
    assert overlap == 2


def test_catalogue_prefixed_bams_are_refused_with_the_reason(tmp_path, monkeypatch):
    """Prefixed names, as written by `coverage -a water=asm --keep-bams`, raise an error."""
    monkeypatch.setattr(
        binning, "capture_in_env", lambda env, argv: _header(["water~k141_1", "water~k141_2"])
    )
    with pytest.raises(ValueError) as problem:
        binning.check_bam_references(["s.bam"], _assembly(tmp_path), "e")
    message = str(problem.value)
    assert "names none of" in message
    # The message shows names from both sides of the mismatch.
    assert "water~k141_1" in message and "k141_1" in message
    assert "reheader" in message


def test_a_partial_overlap_is_reported_not_refused(tmp_path, monkeypatch):
    """A partial overlap is reported but allowed: a catalogue slice can cover extra contigs."""
    said = []
    monkeypatch.setattr(binning, "comm", lambda message: said.append(message))
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: _header(["k141_1"]))
    assert binning.check_bam_references(["s.bam"], _assembly(tmp_path), "e") == 1
    assert said and "1 of 2 assembly contig(s)" in said[0]


def test_a_file_with_no_sq_lines_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr(binning, "capture_in_env", lambda env, argv: "@HD\tVN:1.6\n")
    with pytest.raises(ValueError) as problem:
        binning.check_bam_references(["s.bam"], _assembly(tmp_path), "e")
    assert "no @SQ" in str(problem.value)


def test_supplied_bams_are_linked_not_copied(tmp_path):
    """Supplied BAMs and their indexes are symlinked; they can be tens of GB per sample."""
    bam_dir = _bam_dir(tmp_path, index=True)
    work = tmp_path / "work_files"
    staged = binning.stage_supplied_bams(bam_dir, str(work))
    assert len(staged) == 2
    for link in staged:
        assert os.path.islink(link)
        assert os.path.realpath(link).startswith(bam_dir)
        # An existing index is linked too; samtools looks for it beside the BAM.
        assert os.path.islink(link + ".bai")


def test_staging_twice_is_harmless(tmp_path):
    """Staging a second time (as --resume does) returns the same links without raising."""
    bam_dir = _bam_dir(tmp_path)
    work = tmp_path / "work_files"
    first = binning.stage_supplied_bams(bam_dir, str(work))
    again = binning.stage_supplied_bams(bam_dir, str(work))
    assert first == again


def test_non_bam_files_in_the_directory_are_ignored(tmp_path):
    bam_dir = _bam_dir(tmp_path)
    open(os.path.join(bam_dir, "notes.txt"), "w").close()
    staged = binning.stage_supplied_bams(bam_dir, str(tmp_path / "work"))
    assert [os.path.basename(p) for p in staged] == ["sample_0.bam", "sample_1.bam"]
