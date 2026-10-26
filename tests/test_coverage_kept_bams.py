"""Tests that kept BAMs name contigs as the assembly does, or that a warning says they do not.

`coverage -a water=asm.fa --keep-bams` aligns against a catalogue whose contigs are
`water~<contig>` and splits the depth table back to plain names. Kept BAMs left prefixed matched
zero contigs in SemiBin2, VAMB, COMEBin and CONCOCT, which read the BAM; metaBAT2 reads the table.
"""

import os

from metawrap2.modules import coverage


def _header(prefix=""):
    return (
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:%sk141_1\tLN:4000\n"
        "@SQ\tSN:%sk141_2\tLN:1500\n"
        "@PG\tID:bwa\tCL:bwa mem catalogue.fa\n" % (prefix, prefix)
    )


def test_only_the_sq_names_are_rewritten():
    plain, renamed = coverage.plain_header(_header("water~"), "water", "~")
    assert renamed == 2
    assert "SN:k141_1" in plain and "SN:water~k141_1" not in plain
    # @PG records what ran and is kept as is.
    assert "CL:bwa mem catalogue.fa" in plain
    assert plain.startswith("@HD\tVN:1.6")


def test_a_header_that_is_already_plain_is_left_alone():
    plain, renamed = coverage.plain_header(_header(), "water", "~")
    assert renamed == 0
    assert plain.rstrip("\n") == _header().rstrip("\n")


def test_another_labels_prefix_is_not_stripped():
    """Only this assembly's prefix: a different label's contigs are different sequences."""
    plain, renamed = coverage.plain_header(_header("gut~"), "water", "~")
    assert renamed == 0
    assert "SN:gut~k141_1" in plain


def test_one_assembly_gets_its_bam_headers_rewritten(tmp_path, monkeypatch):
    work = tmp_path / "work_files"
    work.mkdir()
    (work / "sample_0.bam").write_text("bam")
    calls = []
    monkeypatch.setattr(coverage, "capture_in_env", lambda env, argv: _header("water~"))

    def fake_run(cmd, **kwargs):
        # fs_move is real, so create the reheadered file samtools would write.
        calls.append((cmd, kwargs))
        if "stdout_path" in kwargs:
            with open(kwargs["stdout_path"], "w") as handle:
                handle.write("bam")

    monkeypatch.setattr(coverage, "run", fake_run)

    coverage._deprefix_kept_bams([coverage.Assembly("water", "asm.fa")], str(work), "~", "e")

    cmd, kwargs = calls[0]
    assert cmd.startswith("samtools reheader ")
    # In-place reheader is CRAM-only: BAM is rewritten to a new file and then moved over.
    assert kwargs["stdout_path"].endswith(".plain.bam")
    assert os.path.isfile(str(work / "sample_0.bam"))
    assert not os.path.exists(str(work / "sample_0.bam.plain.sam"))
    # The kept BAM is indexed after the rewrite, since reheadering writes a new file.
    assert calls[1][0].startswith("samtools index")


def test_several_assemblies_are_warned_about_rather_than_rewritten(tmp_path, monkeypatch):
    """Stripping the prefix would merge `k141_1` of one assembly with `k141_1` of another."""
    work = tmp_path / "work_files"
    work.mkdir()
    (work / "sample_0.bam").write_text("bam")
    said = []
    ran = []
    monkeypatch.setattr(coverage, "warning", lambda message: said.append(message))
    monkeypatch.setattr(coverage, "run", lambda cmd, **kw: ran.append(cmd))

    coverage._deprefix_kept_bams(
        [coverage.Assembly("water", "a.fa"), coverage.Assembly("gut", "b.fa")], str(work), "~", "e"
    )

    assert not any(
        c.startswith("samtools reheader") for c in ran
    ), "a multi-assembly catalogue must not have its BAM headers rewritten"
    assert any(c.startswith("samtools index") for c in ran), "kept BAMs are indexed either way"
    assert said, "the mismatch must be reported rather than left silent"
    message = said[0]
    # Names the binners that are affected, and the ones that are not.
    for binner in ("SemiBin2", "VAMB", "COMEBin", "CONCOCT", "metaBAT2"):
        assert binner in message
    assert "samtools view -b" in message, "the remedy has to be in the message"


def test_no_bams_is_not_an_error(tmp_path):
    work = tmp_path / "work_files"
    work.mkdir()
    coverage._deprefix_kept_bams([coverage.Assembly("water", "a.fa")], str(work), "~", "e")
