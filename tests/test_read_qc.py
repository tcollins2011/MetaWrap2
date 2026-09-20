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


def test_sample_name():
    assert read_qc._sample_name("/path/sampleA_1.fastq") == "sampleA"
    assert read_qc._sample_name("/path/foo_bar_1.fastq.gz") == "foo_bar"
    assert read_qc._sample_name("/path/x_1.fq") == "x"


def _fastq(names):
    lines = []
    for n in names:
        lines += ["@%s/1" % n, "ACGT", "+", "IIII"]
    return "\n".join(lines) + "\n"


def test_skip_human_reads_drops_host(tmp_path):
    fastq = tmp_path / "r.fastq"
    fastq.write_text(_fastq(["read1", "read2", "read3"]))
    listf = tmp_path / "host.list"
    listf.write_text("read2\n")
    out = io.StringIO()
    skip_human_reads.filter_reads(str(listf), str(fastq), out)
    text = out.getvalue()
    assert "@read1/1" in text and "@read3/1" in text and "@read2/1" not in text


def test_select_human_reads_keeps_only_host(tmp_path):
    fastq = tmp_path / "r.fastq"
    fastq.write_text(_fastq(["read1", "read2", "read3"]))
    listf = tmp_path / "host.list"
    listf.write_text("read2\n")
    out = io.StringIO()
    select_human_reads.select_reads(str(listf), str(fastq), out)
    text = out.getvalue()
    assert "@read2/1" in text and "@read1/1" not in text and "@read3/1" not in text
