from metawrap2.modules import annotate_bins
from metawrap2.scripts import shorten_contig_names


def test_command_template_formats_cleanly():
    cmd = annotate_bins.PROKKA.format(threads=2, outdir="out/binA", prefix="binA", input="tmp.fa")
    assert "{" not in cmd and "--cpus 2" in cmd and "--outdir out/binA" in cmd
    assert "--prefix binA" in cmd and cmd.endswith("tmp.fa")




