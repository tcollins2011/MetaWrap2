"""Collapse a translated KRAKEN2 file into the tab-delimited weight/taxonomy KRONA format.

Ported from metaWRAP's ``kraken_to_krona.py``. Contigs (SPAdes-style headers carrying
``length`` and having more than five ``_``-split fields) are weighted by length*coverage;
everything else (reads) is weighted 1. The taxonomy lineage is expanded from ``;`` to tab
separators as ``ktImportText`` expects.
"""

from __future__ import annotations

import sys
from typing import IO


def to_krona(kraken2_file: str, out: IO) -> None:
    """Summarize *kraken2_file* into KRONA text rows written to *out*."""
    data = {}
    with open(kraken2_file) as fh:
        for line in fh:
            cut = line.strip().split("\t")
            if len(cut) < 2:
                continue
            name = cut[0]
            tax = "\t".join(cut[1].split(";"))
            if "length" in name and len(name.split("_")) > 5:
                weight = float(name.split("_")[3]) * float(name.split("_")[5])
            else:
                weight = 1
            data[tax] = data.get(tax, 0) + weight

    for tax in data:
        out.write(str(data[tax]) + "\t" + tax + "\n")


def main(argv) -> int:
    to_krona(argv[0], sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
