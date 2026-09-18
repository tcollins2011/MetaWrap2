import os

from metawrap2.modules import blobology


def test_command_templates_format_cleanly():
    blastn = blobology.BLASTN.format(query="asm.fa", blastdb="/db", threads=8)
    assert "{" not in blastn and "-num_threads 8" in blastn and "/db/nt" in blastn
    bt2 = blobology.BOWTIE2.format(index="asm.fa", threads=4, r1="a_1.fq", r2="a_2.fq")
    assert "{" not in bt2 and "-p 4" in bt2 and "-1 a_1.fq -2 a_2.fq" in bt2
    plot = blobology.MAKEBLOBPLOT_BASE.format(table="t", prop=0.005, taxlevel="bin", base="Unbinned")
    assert "{" not in plot and plot.endswith("t 0.005 bin Unbinned")


def test_bams_are_sorted(tmp_path):
    for name in ("b.bowtie2.bam", "a.bowtie2.bam", "ignore.txt"):
        (tmp_path / name).write_text("")
    bams = blobology._bams(str(tmp_path))
    assert bams == "%s %s" % (tmp_path / "a.bowtie2.bam", tmp_path / "b.bowtie2.bam")


def test_annotate_with_bins_drops_when_no_match(tmp_path, monkeypatch):
    # Simulate add_bins output where nothing was binned: fewer than two non-"Unbinned" lines.
    blob = tmp_path / "s.blobplot"
    blob.write_text("seqid\tlen\n")
    binned = tmp_path / "s.binned.blobplot"

    def fake_run(cmd, **kw):
        with open(kw["log_path"], "w") as fh:
            fh.write("seqid\tbin\nc1\tUnbinned\nc2\tUnbinned\n")
        return 0

    monkeypatch.setattr(blobology, "run", fake_run)
    assert blobology._annotate_with_bins(str(blob), str(binned), str(tmp_path)) is False
    assert not os.path.exists(str(blob) + ".binned.tmp")


def test_annotate_with_bins_keeps_and_writes_binned(tmp_path, monkeypatch):
    blob = tmp_path / "s.blobplot"
    blob.write_text("old\n")
    binned = tmp_path / "s.binned.blobplot"

    def fake_run(cmd, **kw):
        with open(kw["log_path"], "w") as fh:
            fh.write("seqid\tbin\nc1\tbin.1\nc2\tbin.2\nc3\tUnbinned\n")
        return 0

    monkeypatch.setattr(blobology, "run", fake_run)
    assert blobology._annotate_with_bins(str(blob), str(binned), str(tmp_path)) is True
    # binned blobplot excludes the Unbinned contig
    text = binned.read_text()
    assert "c1" in text and "c2" in text and "Unbinned" not in text
