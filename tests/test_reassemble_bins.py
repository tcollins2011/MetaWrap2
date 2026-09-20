from metawrap2.modules import reassemble_bins


def test_command_templates_format_cleanly():
    idx = reassemble_bins.BWA_INDEX.format(assembly="asm.fa")
    assert idx == "bwa index asm.fa"
    sp = reassemble_bins.SPADES.format(
        threads=4, mem=40, tmp="/t", contigs="c.fa", r1="a_1.fq", r2="a_2.fq", out="/o"
    )
    assert "{" not in sp and "--careful" in sp and "-t 4" in sp and "-o /o" in sp
    spn = reassemble_bins.SPADES_NANOPORE.format(
        threads=1, mem=8, tmp="/t", contigs="c.fa", r1="a_1.fq", r2="a_2.fq", nano="n.fq", out="/o"
    )
    assert "{" not in spn and "--nanopore n.fq" in spn
