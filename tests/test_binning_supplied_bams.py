"""`--bams DIR`: bin from alignments that already exist.

The gap a user hit: `--depth <table>` is refused for SemiBin2, VAMB, COMEBin and CONCOCT, because
those read the alignments themselves - and there was no way to hand them alignments either, even
though `coverage --keep-bams` exists to produce exactly that. Splitting a cohort across a CPU box
(align) and a GPU box (bin) was inexpressible, and had to be faked by pre-placing BAMs in
`work_files/` and seeding the checkpoints by hand.

The header check is the other half of the feature. A BAM aligned against a concatenated catalogue
is as well-formed as the right one, and a BAM binner that matches no contigs produces no bins
without erroring - so supplied BAMs are compared against the assembly before anything runs.
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
    """The silent failure made loud: `coverage -a water=asm --keep-bams` writes these."""
    monkeypatch.setattr(
        binning, "capture_in_env", lambda env, argv: _header(["water~k141_1", "water~k141_2"])
    )
    with pytest.raises(ValueError) as problem:
        binning.check_bam_references(["s.bam"], _assembly(tmp_path), "e")
    message = str(problem.value)
    assert "names none of" in message
    # Both sides of the mismatch, so the user can see which way round it is.
    assert "water~k141_1" in message and "k141_1" in message
    assert "reheader" in message


def test_a_partial_overlap_is_reported_not_refused(tmp_path, monkeypatch):
    """A catalogue slice legitimately covers more contigs than one assembly has."""
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
    """Tens of gigabytes per sample, and this run only reads them."""
    bam_dir = _bam_dir(tmp_path, index=True)
    work = tmp_path / "work_files"
    staged = binning.stage_supplied_bams(bam_dir, str(work))
    assert len(staged) == 2
    for link in staged:
        assert os.path.islink(link)
        assert os.path.realpath(link).startswith(bam_dir)
        # An existing index is linked with it; samtools finds it beside the BAM.
        assert os.path.islink(link + ".bai")


def test_staging_twice_is_harmless(tmp_path):
    """A --resume re-enters this path; a second link attempt must not raise."""
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
