import argparse
import io

from metawrap2.modules import kraken2
from metawrap2.scripts import kraken_to_krona


def test_command_templates_format_cleanly():
    cmd = kraken2.KRAKEN2_PAIRED.format(db="/db", threads=4, out="s.krak2", r1="a", r2="b")
    assert "{" not in cmd and "--paired" in cmd and "--use-names" in cmd

    cmd = kraken2.KRAKEN2_SINGLE.format(db="/db", threads=4, out="s.krak2", seq="asm.fa")
    assert "{" not in cmd and "--paired" not in cmd

    cmd = kraken2.KTIMPORTTEXT.format(out="k.html", krona_files="a.krona b.krona")
    assert "{" not in cmd and "ktImportText -o k.html" in cmd


def test_no_preload_appends_memory_mapping():
    args = argparse.Namespace(no_preload=True)
    base = kraken2.KRAKEN2_SINGLE.format(db="/db", threads=1, out="o", seq="x.fa")
    assert kraken2._kraken2(base, args).endswith("--memory-mapping")
    args.no_preload = False
    assert not kraken2._kraken2(base, args).endswith("--memory-mapping")


def test_is_fasta_and_sample_names():
    assert kraken2._is_fasta("x.fa") and kraken2._is_fasta("x.fasta.gz")
    assert not kraken2._is_fasta("x_1.fastq")
    assert kraken2._fasta_sample("/p/assembly.fasta") == "assembly"
    assert kraken2._fasta_sample("/p/assembly.fa.gz") == "assembly"
    assert kraken2._fastq_sample("/p/sampleA_1.fastq") == "sampleA"
    assert kraken2._fastq_sample("/p/sampleA_1.fastq.gz") == "sampleA"


def test_subsample_pairs(tmp_path):
    def fastq(prefix, n):
        return "\n".join("@%s%d\n%s\n+\n%s" % (prefix, i, "ACGT", "IIII") for i in range(n)) + "\n"

    r1 = tmp_path / "r_1.fastq"
    r2 = tmp_path / "r_2.fastq"
    r1.write_text(fastq("f", 10))
    r2.write_text(fastq("r", 10))
    o1 = tmp_path / "o_1.fastq"
    o2 = tmp_path / "o_2.fastq"
    kraken2._subsample_pairs(str(r1), str(r2), 3, str(o1), str(o2))
    # 3 records * 4 lines each
    assert len([ln for ln in o1.read_text().splitlines() if ln]) == 12
    assert len([ln for ln in o2.read_text().splitlines() if ln]) == 12


def test_kraken_to_krona_weights(tmp_path):
    # a contig line (length*cov weight) and a read line (weight 1)
    kfile = tmp_path / "s.kraken2"
    kfile.write_text("NODE_1_length_100_cov_2.0\tBacteria;Firmicutes\n"
                     "read1\tBacteria;Firmicutes\n")
    out = io.StringIO()
    kraken_to_krona.to_krona(str(kfile), out)
    row = out.getvalue().strip().split("\t")
    # weight = 100*2.0 (contig) + 1 (read) = 201.0
    assert row[0] == "201.0"
    assert row[1] == "Bacteria" and row[2] == "Firmicutes"
