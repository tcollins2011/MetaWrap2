import io

from metawrap2.modules import read_qc
from metawrap2.scripts import select_human_reads, skip_human_reads


def test_command_templates_format_cleanly():
    cmd = read_qc.FASTQC.format(threads=4, outdir="out/pre", files="a_1.fastq a_2.fastq")
    assert "{" not in cmd and "fastqc -q -t 4" in cmd
    # FastQC takes any number of files, so single-end works through the same template
    cmd = read_qc.FASTQC.format(threads=1, outdir="out/pre", files="a.fastq")
    assert "{" not in cmd and cmd.endswith("a.fastq")

    cmd = read_qc.TRIM_GALORE.format(outdir="out", r1="a_1.fastq", r2="a_2.fastq", trim_threads=4)
    assert "{" not in cmd and "--paired" in cmd and "-j 4" in cmd

    cmd = read_qc.BMTAGGER.format(
        bitmask="db/hg38.bitmask",
        srprism="db/hg38.srprism",
        tmpdir="out/tmp",
        r1="a_1.fastq",
        r2="a_2.fastq",
        listfile="out/a.bmtagger.list",
    )
    assert "{" not in cmd and "-q1" in cmd

    # single-end variants
    cmd = read_qc.TRIM_GALORE_SINGLE.format(outdir="out", reads="a.fastq", trim_threads=2)
    assert "{" not in cmd and "--paired" not in cmd and "-j 2" in cmd
    cmd = read_qc.BMTAGGER_SINGLE.format(
        bitmask="b", srprism="s", tmpdir="t", reads="a.fastq", listfile="l"
    )
    assert "{" not in cmd and "-1 a.fastq" in cmd and "-2" not in cmd




def _fastq(names):
    lines = []
    for n in names:
        lines += ["@%s/1" % n, "ACGT", "+", "IIII"]
    return "\n".join(lines) + "\n"




