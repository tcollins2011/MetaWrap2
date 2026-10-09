#!/usr/bin/env python
"""Build a per-bin abundance table from Salmon per-sample count files.

Given a directory of Salmon ``.counts`` files (one per sample), a folder of bins, and the
metagenomic assembly, this computes each bin's abundance in each sample as the
length-weighted median contig coverage: every contig contributes its coverage once per
kilobase of length, and the bin/sample value is the median of that weighted list.

The weighted-median algorithm and its per-contig weighting are kept verbatim, including that a bin
whose contigs are all shorter than 1 kb yields an empty weighted list and therefore a NaN cell.
``.gz`` bins/assemblies are supported via the shared sequence reader.

What is *not* kept is the accumulator reset, which was wrong: see :func:`build_table`. Coverage
leaked from one sample into the next for all but one bin, and that one bin reported NaN throughout.

Three other departures, all about the table's shape rather than its numbers: row labels come from
:func:`metawrap2.constants.bin_label` so that they match the ones the CoverM quantifier writes (a
gzipped bin used to be labelled ``bin_001.fasta`` here and ``bin_001`` there); the count files are
read in sorted order so the sample columns land in the same order on every run; and only the FASTA
files in the bin folder are read, since a bin folder also holds tables and tool subdirectories.
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

    Keyed on :func:`contig_id` to match load_contig_lengths (and the assembly, and salmon's
    transcript names). Keying on the full header made every metaBAT2 bin invisible here.

    Only the FASTA files in the folder are read (see :func:`fasta_files`): a bin folder that has
    been through bin_refinement also holds tables and tool subdirectories, and a subdirectory here
    used to raise a bare ``IsADirectoryError``.
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

    Each cell is the length-weighted median contig coverage of that bin **in that sample alone**.
    The weighted list is built fresh per sample, which is the fix for a bug carried in from the
    original: the accumulators were reset through a bin name left over from the inner loop, so
    exactly one bin - whichever contig happened to be last in the counts file - had its list
    cleared, and cleared partway through the pass that reads it. Every other bin therefore kept
    accumulating, making sample N's cell the median over samples 1..N, and the one bin that was
    cleared reported NaN in every sample. On two bins and two disjoint samples the table came out
    as ``bin_1  1000000  500000`` / ``bin_2  nan  nan`` where the counts said each bin belonged
    entirely to one sample.
    """
    bins = load_bins(bin_folder)
    contig_lengths = load_contig_lengths(assembly)

    # bin -> sample -> abundance. Insertion-ordered by first appearance, as before.
    abundances: Dict[str, Dict[str, float]] = {}
    samples: List[str] = []
    # Sorted, because the sample columns are written in the order the count files are read: an
    # unsorted listdir means two runs over the same data can emit the same table with its columns
    # in a different order. Labels travel with their values either way, so this is reproducibility
    # rather than correctness - but a table that changes shape between identical runs is not one
    # anybody can diff.
    for salmon_file in sorted(os.listdir(quant_dir)):
        sample = ".".join(salmon_file.split(".")[:-2])
        samples.append(sample)
        # Per sample, and discarded at the end of it: a bin's coverage in this sample is not a
        # function of the samples read before it.
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
                # Every contig contributes its coverage once per kilobase of its length, so that a
                # 300 kb contig counts for more of the bin's median than a 2 kb one. Contigs under
                # 1 kb weigh nothing, and a bin made only of those still yields NaN - that part is
                # the original's behaviour and is kept.
                weighted.setdefault(bin_name, []).extend([abun] * (contig_lengths[contig] // 1000))

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

        for bin_name, cov_list in weighted.items():
            abundances[bin_name][sample] = _median(cov_list)

    out.write("Genomic bins\t%s\n" % "\t".join(samples))
    for bin_name, by_sample in abundances.items():
        # bin_label rather than a strip of the last dotted component: on a gzipped bin the latter
        # leaves "bin_001.fasta", which is not what the CoverM path calls the same bin.
        # Columns are written against the sample list, not against this bin's own keys, so a bin
        # missing from one sample's counts file shifts nothing - it gets a NaN in that column.
        row = [str(by_sample.get(sample, float("nan"))) for sample in samples]
        out.write("%s\t%s\n" % (bin_label(bin_name), "\t".join(row)))


def main(argv: List[str]) -> int:
    quant_dir, bin_folder, assembly = argv[0], argv[1], argv[2]
    build_table(quant_dir, bin_folder, assembly, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
