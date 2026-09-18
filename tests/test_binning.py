import os

from metawrap2.modules import binning


def test_command_templates_format_cleanly():
    # The surfaced command templates must fill without leftover placeholders.
    cmd = binning.METABAT2.format(assembly="a.fa", depth="d.txt", out="out",
                                  metabat_len=1500, threads=4)
    assert "{" not in cmd and "metabat2 -i a.fa" in cmd and "-t 4" in cmd


def test_split_concoct_bins(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1 foo\nACGT\n>contig_2\nTTTT\n>contig_3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\ncontig_1,0\ncontig_2,0\ncontig_3,1\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert os.path.exists(out / "bin.0.fa")
    assert os.path.exists(out / "bin.1.fa")
    bin0 = (out / "bin.0.fa").read_text()
    assert ">contig_1" in bin0 and ">contig_2" in bin0
    assert ">contig_3" in (out / "bin.1.fa").read_text()


def test_split_concoct_unbinned(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>c2\nTTTT\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\n")  # c2 unbinned
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert os.path.exists(out / "bin.0.fa")
    assert ">c2" in (out / "unbinned.fa").read_text()
