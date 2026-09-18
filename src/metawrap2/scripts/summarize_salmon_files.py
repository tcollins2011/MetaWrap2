#!/usr/bin/env python
"""Gather & convert Salmon output counts into edgeR-readable ``.counts`` files.

Walks a directory of Salmon output sub-directories (each holding a ``quant.sf``) and
writes one ``<dirname>.counts`` file per sample, with a ``transcript\\tcount`` header and
one row per transcript carrying its TPM value.

Ported faithfully from metaWRAP's ``summarize_salmon_files.py`` (originally
``gather-counts.py`` by C. Titus Brown, 11/2015 —
https://github.com/ngs-docs/2015-nov-adv-rna/). Attribution retained.
"""

from __future__ import annotations

import os
import sys
from typing import List


def process_quant_file(quant_path: str, out_path: str) -> None:
    """Convert one ``quant.sf`` into a ``.counts`` file (``transcript\\tcount``)."""
    sys.stderr.write("Loading counts from: %s\n" % quant_path)
    with open(out_path, "w") as outfp:
        outfp.write("transcript\tcount\n")
        with open(quant_path) as fh:
            for line in fh:
                if line.startswith("Name"):
                    continue
                name, length, eff_length, tpm, count = line.strip().split("\t")
                outfp.write("%s\t%s\n" % (name, float(tpm)))


def summarize(directory: str) -> List[str]:
    """Find every ``quant.sf`` under *directory* and write ``<dirname>.counts`` beside it.

    ``.counts`` files are written into *directory* (named after the sample sub-directory,
    e.g. ``sampleA.quant`` -> ``sampleA.quant.counts``). Returns the sorted list of names.
    """
    sys.stderr.write("Starting in: %s\n" % os.path.abspath(directory))
    quantlist = []
    for root, _dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith("quant.sf"):
                dirname = os.path.basename(root)
                outname = dirname + ".counts"
                process_quant_file(os.path.join(root, filename),
                                   os.path.join(directory, outname))
                quantlist.append(outname)
                break
    return sorted(quantlist)


def main(argv: List[str]) -> int:
    directory = argv[0] if argv else "."
    quantlist = summarize(directory)
    print(",\n".join('"%s"' % i for i in quantlist))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
