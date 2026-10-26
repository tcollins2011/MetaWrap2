"""Scoring a ``bin_abundance_table.tsv`` against a per-sample abundance truth.

The quantifier counterpart to :mod:`metawrap2.benchmark.assess`, which scores binning with AMBER.
A quantifier is scored on how much of each sample each bin accounted for. The truth comes from
CAMI2's ``reads/reads_mapping.tsv.gz``, which gives the source genome of every read per sample, and
so a read fraction per genome.

A bin holds part of a genome (often most of one plus a little of another), so its expected share of
a sample is

    read fraction of its majority genome  x  fraction of that genome's assembly the bin holds

Reads from parts of the genome outside the bin cannot be attributed to it. Without the second
factor a 60%-complete bin would be scored against the whole genome's abundance.

salmon reports length-weighted median contig coverage, CoverM a percentage of the sample, sylph a
k-mer-derived percentage. The metrics are chosen to be independent of scale:

* ``abundance_spearman`` and ``abundance_pearson``: rank and linear correlation across bins within
  a sample, averaged over samples. Spearman measures whether the bins are in the right order and is
  the main metric.
* ``abundance_l1``: L1 distance between estimate and truth after both are renormalised to sum to 1.
  Lower is better, so ``1 - l1/2`` is also reported as ``abundance_l1_similarity`` for the ranker,
  which assumes bigger is better.

The correlations are written out here to avoid a scipy dependency in the benchmark env.
"""

from __future__ import annotations

import gzip
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from ..constants import bin_label
from ..io.seqio import contig_id, iter_fasta


class AbundanceError(RuntimeError):
    """Raised when an abundance truth or estimate cannot be read or lined up with the other."""


#: Columns of a CAMI2 per-sample read mapping. The first is an anonymous read id ending in ``/1`` or
#: ``/2`` for the two mates of a pair; the second is the genome the read came from.
_READ_ID, _READ_GENOME = 0, 1


def read_fractions(mapping: str) -> Dict[str, float]:
    """The fraction of *mapping*'s reads coming from each genome: the per-sample abundance truth.

    Only forward mates are counted, so each pair counts once even when the truth file is missing
    some mates.
    """
    counts: Dict[str, int] = {}
    total = 0
    opener = gzip.open if mapping.endswith(".gz") else open
    with opener(mapping, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= _READ_GENOME or not fields[_READ_ID].endswith("/1"):
                continue
            genome = fields[_READ_GENOME].strip()
            if not genome:
                continue
            counts[genome] = counts.get(genome, 0) + 1
            total += 1
    if not total:
        raise AbundanceError(
            "No forward reads found in %s; cannot build the abundance truth" % mapping
        )
    return {genome: n / total for genome, n in counts.items()}


def bin_composition(
    bins_dir: str, contig_genome: Dict[str, str], contig_length: Dict[str, int]
) -> Dict[str, Dict[str, int]]:
    """For each bin, how many base pairs it holds of each genome, keyed by :func:`bin_label`.

    Lengths come from the truth's contig lengths, so bins and genomes are measured the same way.
    Contigs the truth does not know are not counted against any genome.
    """
    from ..io.seqio import fasta_files

    composition: Dict[str, Dict[str, int]] = {}
    for path in fasta_files(bins_dir):
        by_genome: Dict[str, int] = {}
        for header, seq in iter_fasta(path):
            contig = contig_id(header)
            genome = contig_genome.get(contig)
            if genome is None:
                continue
            by_genome[genome] = by_genome.get(genome, 0) + contig_length.get(contig, len(seq))
        composition[bin_label(os.path.basename(path))] = by_genome
    if not composition:
        raise AbundanceError("No bins found in %s" % bins_dir)
    return composition


def expected_fractions(
    composition: Dict[str, Dict[str, int]],
    genome_bp: Dict[str, int],
    genome_fraction: Dict[str, float],
) -> Dict[str, float]:
    """Expected fraction of one sample per bin: genome abundance x the bin's share of the genome.

    The bin's genome is the one it holds most base pairs of; minor genomes in the bin get no credit.
    A bin whose genome is absent from the sample expects 0.0, which is scored like any other value.
    """
    expected: Dict[str, float] = {}
    for label, by_genome in composition.items():
        if not by_genome:
            expected[label] = 0.0
            continue
        genome = max(by_genome, key=lambda g: by_genome[g])
        total_bp = genome_bp.get(genome, 0)
        recovered = by_genome[genome] / total_bp if total_bp else 0.0
        expected[label] = genome_fraction.get(genome, 0.0) * min(recovered, 1.0)
    return expected


def read_abundance_table(path: str) -> Tuple[List[str], Dict[str, List[float]]]:
    """Read a ``bin_abundance_table.tsv`` into its sample names and one row of values per bin.

    The format is the same for every quantifier. Non-numeric cells (salmon writes ``nan`` for a bin
    with nothing attributed to it) are read as zero so they do not turn the statistics into NaN.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if len(header) < 2:
            raise AbundanceError("%s has no sample columns" % path)
        samples = [s.strip() for s in header[1:]]
        rows: Dict[str, List[float]] = {}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0].strip():
                continue
            values = []
            for cell in fields[1 : len(samples) + 1]:
                try:
                    value = float(cell)
                except ValueError:
                    value = 0.0
                values.append(0.0 if math.isnan(value) else value)
            values += [0.0] * (len(samples) - len(values))
            rows[fields[0].strip()] = values
    if not rows:
        raise AbundanceError("%s has no bin rows" % path)
    return samples, rows


def build_truth(
    bins_dir: str,
    assembly: str,
    contig_mapping: str,
    read_mappings: Dict[str, str],
) -> Dict[str, Dict[str, float]]:
    """Assemble the whole truth: expected fraction per bin, per sample.

    Inputs: the bins being scored, the assembly they came from (for contig lengths), the dataset's
    contig-to-genome mapping, and one read-to-genome mapping per sample. All are files shipped with
    a CAMI2-style dataset; no quantifier output is used.
    """
    from . import assess

    lengths = assess.fasta_lengths(assembly)
    contig_genome, _missing = assess.contig_to_genome(contig_mapping, lengths)
    if not contig_genome:
        raise AbundanceError(
            "%s names no contig that is in %s; cannot build the abundance truth"
            % (contig_mapping, assembly)
        )
    genome_bp: Dict[str, int] = {}
    for contig, genome in contig_genome.items():
        genome_bp[genome] = genome_bp.get(genome, 0) + lengths[contig]
    composition = bin_composition(bins_dir, contig_genome, lengths)
    return {
        sample: expected_fractions(composition, genome_bp, read_fractions(path))
        for sample, path in sorted(read_mappings.items())
    }


def _ranks(values: Sequence[float]) -> List[float]:
    """Ranks with ties averaged, as Spearman's correlation requires."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def pearson(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    """Linear correlation, or None where it is undefined (fewer than two points, or no variance)."""
    n = len(xs)
    if n < 2 or n != len(ys):
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def spearman(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    """Spearman rank correlation, or None where it is undefined."""
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    return pearson(_ranks(xs), _ranks(ys))


def _l1(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    """L1 distance between two vectors renormalised to sum to 1, in [0, 2]."""
    sx, sy = sum(xs), sum(ys)
    if sx <= 0 or sy <= 0:
        return None
    return sum(abs(x / sx - y / sy) for x, y in zip(xs, ys))


def score(table: str, expected_by_sample: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    """Score an abundance table against a truth per sample, averaged over the samples.

    *expected_by_sample* maps a sample name, as spelled in the table header, to the expected
    fraction for every bin. Bins are matched by label; every quantifier's table names rows with
    :func:`metawrap2.constants.bin_label`.

    Samples absent from the truth are skipped. A sample where a metric is undefined (one genome, or
    all estimates identical) is skipped for that metric only.
    """
    samples, rows = read_abundance_table(table)
    unknown = [s for s in samples if s not in expected_by_sample]
    if len(unknown) == len(samples):
        raise AbundanceError(
            "None of %s's samples (%s) appear in the truth (%s)"
            % (table, ", ".join(samples), ", ".join(sorted(expected_by_sample)))
        )
    collected: Dict[str, List[float]] = {
        "spearman": [],
        "pearson": [],
        "l1": [],
        # Bins missing from the quantifier's table, averaged over samples. They are scored as zero
        # abundance and reported here so the omission is visible.
        "missing_bins": [],
    }
    for index, sample in enumerate(samples):
        expected = expected_by_sample.get(sample)
        if expected is None:
            continue
        # Score every bin in the truth; bins absent from the table count as zero. Scoring only the
        # intersection would reward omitting bins, and would compare quantifiers over different
        # bin sets (CoverM writes a row per bin, salmon only the rows it emitted).
        labels = sorted(expected)
        if len(labels) < 2:
            continue
        estimate = [rows[label][index] if label in rows else 0.0 for label in labels]
        truth = [expected[label] for label in labels]
        missing = [label for label in labels if label not in rows]
        if missing:
            collected["missing_bins"].append(float(len(missing)))
        for name, value in (
            ("spearman", spearman(estimate, truth)),
            ("pearson", pearson(estimate, truth)),
            ("l1", _l1(estimate, truth)),
        ):
            if value is not None:
                collected[name].append(value)
    metrics = {
        "abundance_%s" % name: sum(values) / len(values)
        for name, values in collected.items()
        if values
    }
    # The ranker treats bigger as better. L1 between two compositions is in [0, 2], so this maps it
    # onto [1, 0].
    if "abundance_l1" in metrics:
        metrics["abundance_l1_similarity"] = 1.0 - metrics["abundance_l1"] / 2.0
    if not metrics:
        raise AbundanceError(
            "Nothing to score: no sample in %s has two bins that are in the truth" % table
        )
    return metrics
