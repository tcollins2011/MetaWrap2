import argparse
import io

from metawrap2.modules import kraken2
from metawrap2.scripts import kraken_to_krona


def test_command_templates_format_cleanly():
    cmd = kraken2.KRAKEN2_PAIRED.format(db="/db", threads=4, out="s.krak2", r1="a", r2="b")
    assert "{" not in cmd and "--paired" in cmd and "--use-names" in cmd

    cmd = kraken2.KRAKEN2_SINGLE.format(db="/db", threads=4, out="s.krak2", seq="asm.fa")
    assert "{" not in cmd and "--paired" not in cmd

    cmd = kraken2.KTIMPORTTEXT.format(out="k.html", krona_files="a.krona b.krona")
    assert "{" not in cmd and "ktImportText -o k.html" in cmd








