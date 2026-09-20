from metawrap2.modules import blobology


def test_command_templates_format_cleanly():
    blastn = blobology.BLASTN.format(query="asm.fa", blastdb="/db/nt", threads=8)
    assert "{" not in blastn and "-num_threads 8" in blastn and "-db /db/nt" in blastn
    bt2 = blobology.BOWTIE2.format(index="asm.fa", threads=4, r1="a_1.fq", r2="a_2.fq")
    assert "{" not in bt2 and "-p 4" in bt2 and "-1 a_1.fq -2 a_2.fq" in bt2
    plot = blobology.MAKEBLOBPLOT_BASE.format(
        table="t", prop=0.005, taxlevel="bin", base="Unbinned"
    )
    assert "{" not in plot and plot.endswith("t 0.005 bin Unbinned")
