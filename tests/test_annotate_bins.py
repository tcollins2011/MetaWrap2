from metawrap2.modules import annotate_bins
from metawrap2.scripts import shorten_contig_names


def test_command_template_formats_cleanly():
    cmd = annotate_bins.PROKKA.format(threads=2, outdir="out/binA", prefix="binA", input="tmp.fa")
    assert "{" not in cmd and "--cpus 2" in cmd and "--outdir out/binA" in cmd
    assert "--prefix binA" in cmd and cmd.endswith("tmp.fa")


def test_shorten_contig_names_latches(tmp_path):
    fa = tmp_path / "bin.fa"
    fa.write_text(
        ">short_name\nACGT\n"
        ">NODE_1_length_5000_cov_10.5_extra_more\nGGGG\n"
        ">NODE_2_length_10_cov_1_x_y\nTTTT\n"
    )
    out = list(shorten_contig_names.shorten_lines(str(fa)))
    assert out == [
        ">short_name",
        "ACGT",
        ">NODE_1_length_5000",
        "GGGG",
        ">NODE_2_length_10",
        "TTTT",  # shorten latched on -> this short-ish header also trimmed
    ]


def test_grep_product(tmp_path):
    gff = tmp_path / "in.gff"
    gff.write_text(
        "##gff\nc1\tProdigal\tCDS\t1\t9\t.\t+\t0\tID=1;product=hypothetical\nc1\tfoo\tbar\n"
    )
    out = tmp_path / "out.gff"
    annotate_bins._grep_product(str(gff), str(out))
    text = out.read_text()
    assert "product=hypothetical" in text
    assert "foo\tbar" not in text
