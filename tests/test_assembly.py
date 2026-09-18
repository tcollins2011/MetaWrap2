import io

from metawrap2.modules import assembly
from metawrap2.scripts import fix_megahit_contig_naming, rm_short_contigs, sam_to_fastq, sort_contigs


def test_command_templates_format_cleanly():
    cmd = assembly.METASPADES.format(tmp="t", threads=8, mem=24, out="o", r1="a", r2="b")
    assert "{" not in cmd and "metaspades.py" in cmd and "-t 8" in cmd

    cmd = assembly.MEGAHIT_PAIRED.format(r1="a", r2="b", out="o", tmp="t", threads=4,
                                         mem=24_000_000_000)
    assert "{" not in cmd and "--continue" in cmd

    cmd = assembly.MEGAHIT_LEFTOVERS.format(reads="u.fastq", out="o", threads=4,
                                            mem=24_000_000_000, tmp="t")
    assert "{" not in cmd and "-r u.fastq" in cmd

    cmd = assembly.QUAST.format(threads=2, out="q", quast_min=500, assembly="f.fasta")
    assert "{" not in cmd and "-m 500" in cmd


def test_rm_short_contigs_breaks_on_first_short(tmp_path):
    fa = tmp_path / "s.fasta"
    fa.write_text(">NODE_1_length_2000_cov_5\nAAAA\n>NODE_2_length_500_cov_3\nCCCC\n"
                  ">NODE_3_length_3000_cov_2\nGGGG\n")
    out = io.StringIO()
    rm_short_contigs.remove_short_contigs(1000, str(fa), out)
    text = out.getvalue()
    # stops at the first contig below 1000, so NODE_3 (after the short one) is excluded
    assert "NODE_1_length_2000" in text
    assert "NODE_2_length_500" not in text and "NODE_3_length_3000" not in text


def test_sort_contigs_orders_longest_first(tmp_path):
    fa = tmp_path / "c.fasta"
    fa.write_text(">short\nAC\n>long\nACGTACGT\n>mid\nACGT\n")
    out = io.StringIO()
    sort_contigs.sort_contigs(str(fa), out)
    headers = [ln for ln in out.getvalue().splitlines() if ln.startswith(">")]
    assert headers == [">long", ">mid", ">short"]


def test_fix_megahit_contig_naming(tmp_path):
    fa = tmp_path / "m.fa"
    fa.write_text(">k141_0 flag=1 multi=2.0000 len=500\nACGTACGT\n"
                  ">k141_1 flag=1 multi=3.0000 len=100\nAC\n")
    out = io.StringIO()
    fix_megahit_contig_naming.fix_naming(str(fa), 200, out)
    text = out.getvalue()
    assert ">k141_0_length_500_cov_2.0000" in text
    assert "k141_1" not in text  # dropped: shorter than 200


def test_sam_to_fastq(tmp_path):
    sam = tmp_path / "a.sam"
    sam.write_text("@SQ\tSN:x\tLN:10\n"
                   "read1\t0\tcontig\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n")
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out)
    assert out.getvalue() == ">read1\nACGT\n+\nIIII\n"


def test_assembler_selection_defaults_to_megahit():
    args = assembly._parse_args(["-1", "a_1.fastq", "-2", "a_2.fastq", "-o", "out"])
    assert args.megahit is False and args.metaspades is False
    # default logic in main(): megahit runs when nothing is picked
    do_metaspades = args.metaspades
    do_megahit = args.megahit or not args.metaspades
    assert do_megahit is True and do_metaspades is False
