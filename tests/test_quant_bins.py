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






