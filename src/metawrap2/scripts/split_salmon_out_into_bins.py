#!/usr/bin/env python
"""Build a per-bin abundance table from Salmon per-sample count files.

Given a directory of Salmon ``.counts`` files (one per sample), a folder of bins, and the
metagenomic assembly, this computes each bin's abundance in each sample as the
length-weighted median contig coverage: every contig contributes its coverage once per
kilobase of length, and the bin/sample value is the median of that weighted list.

Ported faithfully from metaWRAP's ``split_salmon_out_into_bins.py`` (the weighted-median
algorithm and its edge behavior are preserved verbatim, including that a bin whose contigs
are all shorter than 1 kb yields an empty weighted list and therefore a NaN cell). ``.gz``
bins/assemblies are supported via the shared sequence reader.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, TextIO

from ..io.seqio import iter_fasta


def _median(values: List[float]) -> float:
    """Median matching numpy.median: mean of the two central values; NaN on empty input."""
    n = len(values)
    if n == 0:
        return float("nan")
    ordered = sorted(values)
    mid = n // 2
    if n % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def load_bins(bin_folder: str) -> Dict[str, str]:
    """Map each contig header to the bin filename that contains it."""
    bins = {}
    for bin_name in os.listdir(bin_folder):
        for header, _seq in iter_fasta(os.path.join(bin_folder, bin_name)):
            bins[header] = bin_name
    return bins


def load_contig_lengths(assembly: str) -> Dict[str, int]:
    """Map each contig (first whitespace-delimited id token) to its length."""
    lengths = {}
    for header, seq in iter_fasta(assembly):
        lengths[header.split()[0]] = len(seq)
    return lengths


def build_table(quant_dir: str, bin_folder: str, assembly: str, out: TextIO) -> None:
    """Write the ``Genomic bins`` abundance table to *out* (tab-separated)."""
    bins = load_bins(bin_folder)
    contig_lengths = load_contig_lengths(assembly)

    bin_abundances: Dict[str, dict] = {}
    for salmon_file in os.listdir(quant_dir):
        sample = ".".join(salmon_file.split(".")[:-2])
        ct = 0
        with open(os.path.join(quant_dir, salmon_file)) as fh:
            for line in fh:
                if "transcript" in line:
                    continue
                contig = line.split("\t")[0]
                if contig not in bins:
                    continue
                ct += 1
                bin_name = bins[contig]
                abun = float(line.strip().split("\t")[1])
                if bin_name not in bin_abundances:
                    bin_abundances[bin_name] = {
                        "total_len": 0, "total_cov": 0, "cov_list": [], "samples": {},
                    }
                length = contig_lengths[contig]
                weight = length // 1000
                bin_abundances[bin_name]["cov_list"].extend([abun] * weight)
                bin_abundances[bin_name]["total_len"] += length
                bin_abundances[bin_name]["total_cov"] += abun * length

        if ct == 0:
            sys.stderr.write(
                "\nNone of the contigs/scaffolds in the -a metagenomic assembly file were "
                "present in the bin files. Please make sure that the bins and total assembly "
                "have the exact same bins. One cause for this could be that you reassembled "
                "the bins, disrupting the contig naming. If you do not have the original total "
                "metagenomic assembly file, then you could not provide the -a option at all "
                "(but this is not ideal for abundance estimation).\n"
            )
            sys.exit(1)

        for bin_name in bin_abundances:
            bin_abundances[bin_name]["samples"][sample] = _median(
                bin_abundances[bin_name]["cov_list"])
            bin_abundances[bin_name]["total_len"] = 0
            bin_abundances[bin_name]["total_cov"] = 0
            bin_abundances[bin_name]["cov_list"] = []

    first = True
    for bin_name in bin_abundances:
        if first:
            out.write("Genomic bins")
            for sample in bin_abundances[bin_name]["samples"]:
                out.write("\t" + sample)
            out.write("\n")
            first = False
        out.write(".".join(bin_name.split(".")[:-1]))
        for sample in bin_abundances[bin_name]["samples"]:
            out.write("\t" + str(bin_abundances[bin_name]["samples"][sample]))
        out.write("\n")


def main(argv: List[str]) -> int:
    quant_dir, bin_folder, assembly = argv[0], argv[1], argv[2]
    build_table(quant_dir, bin_folder, assembly, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
