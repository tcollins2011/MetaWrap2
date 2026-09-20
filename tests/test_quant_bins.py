import io
import math
import sys

from metawrap2.modules import quant_bins
from metawrap2.scripts import split_salmon_out_into_bins, summarize_salmon_files


def test_command_templates_format_cleanly():
    idx = quant_bins.SALMON_INDEX.format(threads=4, assembly="a.fa", index="idx")
    assert "{" not in idx and "salmon index -p 4" in idx and "-i idx" in idx
    q = quant_bins.SALMON_QUANT.format(index="idx", r1="a_1.fq", r2="a_2.fq", out="o", threads=4)
    assert "{" not in q and "--libType IU" in q and "--meta -p 4" in q
    # The heatmap helper is a metawrap2 module, so it must run in the *host* interpreter -
    # the quant_bins conda env holds salmon, not python+seaborn.
    heat = [a.format(table="t.tab", png="h.png") for a in quant_bins.MAKE_HEATMAP]
    assert heat == [sys.executable, "-m", "metawrap2.scripts.make_heatmap", "t.tab", "h.png"]


def test_summarize_writes_counts(tmp_path):
    align = tmp_path / "align"
    q = align / "sampleA.quant"
    q.mkdir(parents=True)
    (q / "quant.sf").write_text(
        "Name\tLength\tEffLength\tTPM\tNumReads\nc1\t100\t80\t5.0\t10\nc2\t200\t150\t3.0\t7\n"
    )
    names = summarize_salmon_files.summarize(str(align))
    assert names == ["sampleA.quant.counts"]
    counts = (align / "sampleA.quant.counts").read_text()
    assert counts.startswith("transcript\tcount\n")
    assert "c1\t5.0" in counts and "c2\t3.0" in counts


def test_median_matches_numpy_semantics():
    assert math.isnan(split_salmon_out_into_bins._median([]))
    assert split_salmon_out_into_bins._median([1.0, 2.0, 3.0]) == 2.0
    assert split_salmon_out_into_bins._median([1.0, 2.0, 3.0, 4.0]) == 2.5


def test_build_table_weighted_median(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">c1\n" + "A" * 2000 + "\n>c2\n" + "C" * 3000 + "\n")
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\n" + "A" * 2000 + "\n>c2\n" + "C" * 3000 + "\n")
    quant = tmp_path / "quant"
    quant.mkdir()
    (quant / "sampleA.quant.counts").write_text("transcript\tcount\nc1\t4.0\nc2\t4.0\n")

    out = io.StringIO()
    split_salmon_out_into_bins.build_table(str(quant), str(bins), str(assembly), out)
    text = out.getvalue()
    assert text.startswith("Genomic bins\tsampleA\n")
    assert "bin.1\t4.0" in text
