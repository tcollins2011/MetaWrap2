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






