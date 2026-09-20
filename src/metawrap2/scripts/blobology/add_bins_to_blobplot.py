#!/usr/bin/env python
"""Annotate a blobplot table with which bin each contig ended up in.

Adds three columns to every row of a blobplot file, so the plots can be coloured by bin as
well as by taxonomy:

``bin``
    the bin this contig belongs to, or ``Unbinned``
``binned_yes_no``
    ``Binned`` / ``Unbinned`` - for a two-colour "what did we recover?" plot
``binned_phylum``
    the contig's phylum if it was binned, else ``Unbinned`` - so the phylum plot shows only
    the part of the community that made it into genomes

The annotated table goes to stdout; ``blobology`` redirects it back over the input file.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, TextIO

from metawrap2.constants import FASTA_EXTENSIONS
from metawrap2.io.seqio import contig_id, iter_fasta

UNBINNED = "Unbinned"
NEW_COLUMNS = ("bin", "binned_yes_no", "binned_phylum")


def load_contig_bins(bin_folder: str) -> Dict[str, str]:
    """Map each contig id to the name of the bin file that contains it.

    Keyed on :func:`contig_id` - the header up to the first whitespace - which is what the
    blobplot's ``seqid`` column holds. The original keyed on the whole header line, so a bin
    whose contigs carry tool-added annotations (metaBAT2 writes ``total_depth=..``) matched
    nothing and the blobplot came out entirely unbinned.
    """
    contig_bins: Dict[str, str] = {}
    for bin_file in sorted(os.listdir(bin_folder)):
        if not bin_file.endswith(FASTA_EXTENSIONS):
            continue
        for header, _seq in iter_fasta(os.path.join(bin_folder, bin_file)):
            contig_bins[contig_id(header)] = bin_file
    return contig_bins


def annotate(blobplot_path: str, contig_bins: Dict[str, str], out: TextIO) -> int:
    """Write *blobplot_path* to *out* with the bin columns appended. Returns rows annotated."""
    phylum_column: Optional[int] = None
    binned = 0
    with open(blobplot_path) as fh:
        for line in fh:
            if line == "\n":
                continue
            row = line.rstrip("\n")
            fields = row.split("\t")
            if fields[0] == "seqid":
                for i, field in enumerate(fields):
                    if field == "taxlevel_phylum":
                        phylum_column = i
                out.write(row + "\t" + "\t".join(NEW_COLUMNS) + "\n")
                continue

            bin_name = contig_bins.get(fields[0])
            if bin_name is None:
                out.write("%s\t%s\t%s\t%s\n" % (row, UNBINNED, UNBINNED, UNBINNED))
                continue
            phylum = fields[phylum_column] if phylum_column is not None else UNBINNED
            out.write("%s\t%s\tBinned\t%s\n" % (row, bin_name, phylum))
            binned += 1
    return binned


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        sys.stderr.write(
            "usage: python -m metawrap2.scripts.blobology."
            "add_bins_to_blobplot <blobplot> <bin folder>\n"
        )
        return 2
    blobplot, bin_folder = argv
    annotate(blobplot, load_contig_bins(bin_folder), sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
