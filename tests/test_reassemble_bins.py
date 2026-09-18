from metawrap2.modules import reassemble_bins


def test_command_templates_format_cleanly():
    idx = reassemble_bins.BWA_INDEX.format(assembly="asm.fa")
    assert idx == "bwa index asm.fa"
    sp = reassemble_bins.SPADES.format(threads=4, mem=40, tmp="/t", contigs="c.fa",
                                       r1="a_1.fq", r2="a_2.fq", out="/o")
    assert "{" not in sp and "--careful" in sp and "-t 4" in sp and "-o /o" in sp
    spn = reassemble_bins.SPADES_NANOPORE.format(threads=1, mem=8, tmp="/t", contigs="c.fa",
                                                 r1="a_1.fq", r2="a_2.fq", nano="n.fq", out="/o")
    assert "{" not in spn and "--nanopore n.fq" in spn


def test_bwa_mem_filter_pipeline_shape():
    cmd = reassemble_bins.BWA_MEM_FILTER.format(threads=2, assembly="asm.fa", r1="a_1.fq",
                                                r2="a_2.fq", bins="bins", reads_out="ro",
                                                strict=2, permissive=5)
    assert "{" not in cmd and "bwa mem" in cmd and "| python -m" in cmd
    assert cmd.endswith("bins ro 2 5")


def test_combine_bins_concatenates(tmp_path):
    orig = tmp_path / "original_bins"
    orig.mkdir()
    (orig / "bin.1.fa").write_text(">c1\nACGT\n")
    (orig / "bin.2.fa").write_text(">c2\nTTTT\n")
    assembly = tmp_path / "binned_assembly" / "assembly.fa"
    reassemble_bins._combine_bins(str(orig), str(assembly))
    text = assembly.read_text()
    assert ">c1" in text and ">c2" in text


def test_count_fa(tmp_path):
    d = tmp_path / "reassembled_bins"
    d.mkdir()
    (d / "bin.1.fa").write_text(">a\nAC\n")
    assert reassemble_bins._count_fa(str(d)) == 1
    assert reassemble_bins._count_fa(str(tmp_path / "missing")) == 0
