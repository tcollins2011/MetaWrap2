import bz2
import gzip

from metawrap2.io import seqio

FASTA = ">contig_1 some description\nACGTACGT\nACGT\n>contig_2\nTTTT\n"


def test_detect_and_read_plain(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(FASTA)
    assert seqio.detect_compression(str(p)) == "none"
    recs = list(seqio.iter_fasta(str(p)))
    assert recs == [("contig_1 some description", "ACGTACGTACGT"), ("contig_2", "TTTT")]


def test_detect_and_read_gzip(tmp_path):
    p = tmp_path / "x.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(FASTA)
    assert seqio.detect_compression(str(p)) == "gzip"
    recs = list(seqio.iter_fasta(str(p)))
    assert [h for h, _ in recs] == ["contig_1 some description", "contig_2"]


def test_detect_and_read_bzip2(tmp_path):
    p = tmp_path / "x.fa.bz2"
    with bz2.open(p, "wt") as fh:
        fh.write(FASTA)
    assert seqio.detect_compression(str(p)) == "bzip2"
    recs = list(seqio.iter_fasta(str(p)))
    assert recs[1] == ("contig_2", "TTTT")


def test_detection_uses_content_not_extension(tmp_path):
    # A gzip file mislabeled with a .fa extension is still read correctly.
    p = tmp_path / "mislabeled.fa"
    with gzip.open(p, "wt") as fh:
        fh.write(FASTA)
    assert seqio.detect_compression(str(p)) == "gzip"
    assert len(list(seqio.iter_fasta(str(p)))) == 2


def test_smart_open_write_roundtrip_gz(tmp_path):
    p = tmp_path / "out.txt.gz"
    with seqio.smart_open(str(p), "wt") as fh:
        fh.write("hello\n")
    with seqio.smart_open(str(p), "rt") as fh:
        assert fh.read() == "hello\n"
    assert seqio.detect_compression(str(p)) == "gzip"
