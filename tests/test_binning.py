import os

from metawrap2.modules import binning


def test_command_templates_format_cleanly():
    # The surfaced command templates must fill without leftover placeholders.
    cmd = binning.METABAT2.format(
        assembly="a.fa", depth="d.txt", out="out", metabat_len=1500, threads=4
    )
    assert "{" not in cmd and "metabat2 -i a.fa" in cmd and "-t 4" in cmd


def test_split_concoct_bins(tmp_path):
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">contig_1 foo\nACGT\n>contig_2\nTTTT\n>contig_3\nGGGG\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\ncontig_1,0\ncontig_2,0\ncontig_3,1\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    # CONCOCT clusters 0 and 1 become bins 1 and 2: bin names count from 1 for every binner.
    assert os.path.exists(out / "bin_001.fasta")
    assert os.path.exists(out / "bin_002.fasta")
    first = (out / "bin_001.fasta").read_text()
    assert ">contig_1" in first and ">contig_2" in first
    assert ">contig_3" in (out / "bin_002.fasta").read_text()


