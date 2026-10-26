#!/usr/bin/env python
"""Build a per-bin abundance table from Salmon per-sample count files.

Given a directory of Salmon ``.counts`` files (one per sample), a folder of bins, and the
metagenomic assembly, this computes each bin's abundance in each sample as the
length-weighted median contig coverage: every contig contributes its coverage once per
kilobase of length, and the bin/sample value is the median of that weighted list.

The weighted median and per-contig weighting are unchanged from the original, so a bin whose
contigs are all shorter than 1 kb has an empty weighted list and a NaN cell. ``.gz`` bins and
assemblies are read through the shared sequence reader.

Changes from the original:

- the per-sample accumulator reset is fixed (see :func:`build_table`);
- row labels come from :func:`metawrap2.constants.bin_label`, matching the CoverM quantifier (a
  gzipped bin was labelled ``bin_001.fasta`` here and ``bin_001`` there);
- count files are read in sorted order, so sample columns have the same order on every run;
- only FASTA files in the bin folder are read, since it also holds tables and tool subdirectories.
"""

from __future__ import annotations

import os
import sys
from typing import Dict, List, TextIO

from ..constants import bin_label
from ..io.seqio import contig_id, fasta_files, iter_fasta


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
    """Map each contig id to the bin filename that contains it.

    Keyed on :func:`contig_id` to match load_contig_lengths, the assembly and salmon's transcript
    names. Keying on the full header missed every metaBAT2 bin.

    Only FASTA files are read (see :func:`fasta_files`). A bin folder from bin_refinement also
    holds tables and tool subdirectories, and a subdirectory raised ``IsADirectoryError``.
    """
    bins = {}
    for path in fasta_files(bin_folder):
        for header, _seq in iter_fasta(path):
            bins[contig_id(header)] = os.path.basename(path)
    return bins


def load_contig_lengths(assembly: str) -> Dict[str, int]:
    """Map each contig (first whitespace-delimited id token) to its length."""
    lengths = {}
    for header, seq in iter_fasta(assembly):
        lengths[contig_id(header)] = len(seq)
    return lengths


def build_table(quant_dir: str, bin_folder: str, assembly: str, out: TextIO) -> None:
    """Write the ``Genomic bins`` abundance table to *out* (tab-separated).

    Each cell is the length-weighted median contig coverage of that bin in that sample only. The
    weighted list is rebuilt per sample. The original reset the accumulators through a bin name
    left over from the inner loop, so only the bin of the last contig in the counts file was
    cleared. Every other bin kept accumulating (sample N's cell was the median over samples 1..N)
    and the cleared bin was NaN in every sample.
    """
    bins = load_bins(bin_folder)
    contig_lengths = load_contig_lengths(assembly)

    # bin -> sample -> abundance. Insertion-ordered by first appearance, as before.
    abundances: Dict[str, Dict[str, float]] = {}
    samples: List[str] = []
    # Sorted so the sample columns have the same order on every run and tables can be diffed.
    for salmon_file in sorted(os.listdir(quant_dir)):
        sample = ".".join(salmon_file.split(".")[:-2])
        samples.append(sample)
        # Rebuilt per sample so earlier samples do not leak into this one.
        weighted: Dict[str, List[float]] = {}
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
                abundances.setdefault(bin_name, {})
                # Each contig adds its coverage once per kilobase of length. Contigs under 1 kb
                # add nothing, so a bin made only of those is NaN, as in the original.
                weighted.setdefault(bin_name, []).extend([abun] * (contig_lengths[contig] // 1000))

        if ct == 0:
            sys.stderr.write(
                "\nNone of the contigs/scaffolds in the -a metagenomic assembly file were "
                "present in the bin files. Check that the bins were made from this assembly "
                "and use the same contig names. Reassembling the bins changes the contig "
                "names. If you do not have the original metagenomic assembly file, omit "
                "the -a option, though abundance estimates are less accurate without "
                "it.\n"
            )
            sys.exit(1)

        for bin_name, cov_list in weighted.items():
            abundances[bin_name][sample] = _median(cov_list)

    out.write("Genomic bins\t%s\n" % "\t".join(samples))
    for bin_name, by_sample in abundances.items():
        # bin_label, because stripping the last suffix leaves "bin_001.fasta" on a gzipped bin.
        # Columns follow the sample list, so a bin missing from one sample gets NaN there.
        row = [str(by_sample.get(sample, float("nan"))) for sample in samples]
        out.write("%s\t%s\n" % (bin_label(bin_name), "\t".join(row)))


def main(argv: List[str]) -> int:
    quant_dir, bin_folder, assembly = argv[0], argv[1], argv[2]
    build_table(quant_dir, bin_folder, assembly, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
