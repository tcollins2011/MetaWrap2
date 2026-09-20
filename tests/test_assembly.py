import io

from metawrap2.modules import assembly
from metawrap2.scripts import (
    fix_megahit_contig_naming,
    rm_short_contigs,
    sam_to_fastq,
    sort_contigs,
)


def test_command_templates_format_cleanly():
    cmd = assembly.METASPADES.format(tmp="t", threads=8, mem=24, out="o", r1="a", r2="b")
    assert "{" not in cmd and "metaspades.py" in cmd and "-t 8" in cmd

    cmd = assembly.MEGAHIT_PAIRED.format(
        r1="a", r2="b", out="o", tmp="t", threads=4, mem=24_000_000_000
    )
    assert "{" not in cmd and "--continue" in cmd

    cmd = assembly.MEGAHIT_LEFTOVERS.format(
        reads="u.fastq", out="o", threads=4, mem=24_000_000_000, tmp="t"
    )
    assert "{" not in cmd and "-r u.fastq" in cmd

    cmd = assembly.QUAST.format(threads=2, out="q", quast_min=500, assembly="f.fasta")
    assert "{" not in cmd and "-m 500" in cmd


def test_rm_short_contigs_filters_by_actual_sequence_length(tmp_path):
    # Length comes from the sequence, not from parsing the header, so this works on any
    # assembler's naming and does not depend on the file being sorted longest-first.
    fa = tmp_path / "s.fasta"
    fa.write_text(
        ">NODE_1_length_2000_cov_5\n" + "A" * 2000 + "\n"
        ">NODE_2_length_500_cov_3\n" + "C" * 500 + "\n"
        ">NODE_3_length_3000_cov_2\n" + "G" * 3000 + "\n"
    )
    out = io.StringIO()
    rm_short_contigs.remove_short_contigs(1000, str(fa), out)
    text = out.getvalue()
    assert "NODE_1_length_2000" in text
    assert "NODE_3_length_3000" in text  # no longer lost to the early break
    assert "NODE_2_length_500" not in text


def test_rm_short_contigs_works_on_non_spades_headers(tmp_path):
    """The old header parse raised ValueError on anything but SPAdes naming."""
    fa = tmp_path / "s.fasta"
    fa.write_text(">contig_A\n" + "A" * 2000 + "\n>contig_B\n" + "C" * 100 + "\n")
    out = io.StringIO()
    rm_short_contigs.remove_short_contigs(1000, str(fa), out)
    text = out.getvalue()
    assert ">contig_A" in text and ">contig_B" not in text


def test_sort_contigs_orders_longest_first(tmp_path):
    fa = tmp_path / "c.fasta"
    fa.write_text(">short\nAC\n>long\nACGTACGT\n>mid\nACGT\n")
    out = io.StringIO()
    sort_contigs.sort_contigs(str(fa), out)
    headers = [ln for ln in out.getvalue().splitlines() if ln.startswith(">")]
    assert headers == [">long", ">mid", ">short"]


def test_fix_megahit_contig_naming_reads_fields_by_name(tmp_path):
    """MEGAHIT's header field order has changed between releases; parse by name."""
    fa = tmp_path / "m.fa"
    fa.write_text(">k141_0 multi=2.0000 len=500 flag=1\n" + "A" * 500 + "\n")
    out = io.StringIO()
    fix_megahit_contig_naming.fix_naming(str(fa), 200, out)
    assert ">k141_0_length_500_cov_2.0000" in out.getvalue()


def test_fix_megahit_contig_naming_without_len_field(tmp_path):
    fa = tmp_path / "m.fa"
    fa.write_text(">k141_7\n" + "A" * 300 + "\n")
    out = io.StringIO()
    fix_megahit_contig_naming.fix_naming(str(fa), 200, out)
    assert ">k141_7_length_300_cov_0" in out.getvalue()


def test_fix_megahit_contig_naming(tmp_path):
    fa = tmp_path / "m.fa"
    fa.write_text(
        ">k141_0 flag=1 multi=2.0000 len=500\nACGTACGT\n"
        ">k141_1 flag=1 multi=3.0000 len=100\nAC\n"
    )
    out = io.StringIO()
    fix_megahit_contig_naming.fix_naming(str(fa), 200, out)
    text = out.getvalue()
    assert ">k141_0_length_500_cov_2.0000" in text
    assert "k141_1" not in text  # dropped: shorter than 200


def test_sam_to_fastq(tmp_path):
    sam = tmp_path / "a.sam"
    sam.write_text(
        "@SQ\tSN:x\tLN:10\n" "read1\t0\tcontig\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\n"
    )
    out = io.StringIO()
    sam_to_fastq.sam_to_fastq(str(sam), out)
    # A valid FASTQ record: the original emitted ">" here, making the .fastq neither valid
    # FASTA nor valid FASTQ. MEGAHIT accepts both identically, so this is safe.
    assert out.getvalue() == "@read1\nACGT\n+\nIIII\n"


def test_assembler_selection_defaults_to_megahit():
    args = assembly._parse_args(["-1", "a_1.fastq", "-2", "a_2.fastq", "-o", "out"])
    assert args.megahit is False and args.metaspades is False
    # default logic in main(): megahit runs when nothing is picked
    do_metaspades = args.metaspades
    do_megahit = args.megahit or not args.metaspades
    assert do_megahit is True and do_metaspades is False
