"""Scoring a ``bin_abundance_table.tsv`` against a per-sample abundance truth.

This is the counterpart to :mod:`metawrap2.benchmark.assess` for the quantifiers. Binning is scored
by asking which contigs went in which bin, and AMBER does that; a quantifier is scored by asking how
much of each sample each bin accounted for, which is a different question with a different truth
file. CAMI2 ships one: ``reads/reads_mapping.tsv.gz`` says which genome every read came from, per
sample, which is a read fraction per genome and therefore an abundance per genome.

Turning that into a truth per *bin* takes one more step, and it is the step worth reading carefully.
A bin is not a genome: it is part of one, and usually a large part of one plus a little of another.
So the sequence a bin should account for in a sample is

    the read fraction of the genome it mostly is  x  how much of that genome's assembly it holds

because reads from the parts of the genome that did not make it into the bin cannot be attributed to
the bin by any quantifier, however good. Skipping the second factor scores a 60%-complete bin as
though it should have reported the whole genome's abundance, which penalises the quantifier for the
binner's incompleteness.

**Three quantifiers, three units, and what that means for the metrics.** salmon reports a
length-weighted median contig coverage, CoverM a percentage of the sample, sylph a k-mer-derived
percentage. Only the last two are on the same scale as the truth at all, so the metric that decides
anything has to be one that does not care about scale:

* ``abundance_spearman`` and ``abundance_pearson`` are rank and linear correlation across bins,
  within a sample, averaged over samples. Both are invariant to the unit, so all three quantifiers
  are compared on the same footing. Spearman is the one to lean on - it asks whether the quantifier
  ordered the bins correctly, which is what an abundance table is usually read for.
* ``abundance_l1`` is the L1 distance between the estimate and the truth after both have been
  renormalised to sum to 1, which makes a coverage table and a percentage table comparable as
  *compositions*. Lower is better, so it is reported as ``1 - l1/2`` under
  ``abundance_l1_similarity`` as well, for the ranker, which assumes bigger is better.

No scipy: the correlations are short enough to write out, and the benchmark harness already avoids
pulling a numerical stack into an env that does not otherwise need one.
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

    Only forward mates are counted. Both mates of a pair come from the same genome, so counting both
    would give the same fractions - but a truth file that has lost some mates would then weight those
    pairs differently from complete ones, and counting pairs is what "how much of this sample is this
    organism" means.
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
            "No forward reads found in %s, so there is no abundance truth to score against"
            % mapping
        )
    return {genome: n / total for genome, n in counts.items()}


def bin_composition(
    bins_dir: str, contig_genome: Dict[str, str], contig_length: Dict[str, int]
) -> Dict[str, Dict[str, int]]:
    """For each bin, how many base pairs it holds of each genome, keyed by :func:`bin_label`.

    Lengths come from the truth's own contig lengths rather than from the bin file, so that a bin and
    the genome it came from are measured with one ruler; a contig the truth does not know is counted
    against no genome at all, which is the honest thing to do with sequence nobody can attribute.
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
    """What each bin should account for in one sample: its genome's abundance, times its share of it.

    The bin's genome is the one it holds most base pairs of. A bin holding a little of a second
    genome is not given credit for that second genome's abundance: a quantifier cannot be expected to
    separate the two, and a bin that is mostly genome A is read as genome A by everybody downstream.

    A bin whose genome is absent from this sample expects zero, which is a real expectation and not a
    missing value - a quantifier reporting a genome that is not there is wrong in a way that matters.
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

    Whichever quantifier wrote it: the whole point of the shared table is that this function does not
    have to know. A cell that is not a number - salmon writes ``nan`` for a bin it could attribute
    nothing to - is read as zero, because "no abundance detected" is what it means here and a NaN
    would poison every statistic computed from the column.
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

    The four inputs are the four things the truth is made of - the bins being scored, the assembly
    they were binned from (for contig lengths), the dataset's contig-to-genome mapping, and one
    read-to-genome mapping per sample. Each is a file a CAMI2-shaped dataset already ships; nothing
    here is computed from a quantifier's output, which is the point.
    """
    from . import assess

    lengths = assess.fasta_lengths(assembly)
    contig_genome, _missing = assess.contig_to_genome(contig_mapping, lengths)
    if not contig_genome:
        raise AbundanceError(
            "%s named no contig that %s contains, so there is no abundance truth to build"
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
    """Ranks with ties averaged, which is what makes the rank correlation Spearman's."""
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
    """Rank correlation: did the quantifier put the bins in the right order?"""
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

    *expected_by_sample* maps a sample name - as the table's own header spells it - to the expected
    fraction for every bin. Bins are matched by label, which is why every quantifier's table goes
    through :func:`metawrap2.constants.bin_label`: the row names agree by construction.

    A sample the truth says nothing about is skipped rather than scored as zero, because a missing
    truth file is not a quantifier's mistake. A sample where the metric is undefined - one genome, or
    every estimate identical - is skipped for that metric alone.
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
        # How many of the sample's bins the quantifier left out of its table, averaged over
        # samples. Reported rather than silently absorbed: it is scored as zero abundance, and a
        # reader comparing two quantifiers should see which one declined to answer.
        "missing_bins": [],
    }
    for index, sample in enumerate(samples):
        expected = expected_by_sample.get(sample)
        if expected is None:
            continue
        # Every bin the truth knows about, not the intersection with what the quantifier wrote.
        # A bin absent from the table is one the quantifier attributed nothing to, which is a zero
        # for the same reason salmon's own `nan` cells are read as zero - dropping it instead
        # scored the tool on a smaller, easier set of bins, and rewarded omitting the hard ones.
        # It also made two quantifiers incomparable: the CoverM path pre-seeds a row per bin while
        # the salmon path only writes rows salmon emitted, so salmon could be scored over 40 bins
        # and CoverM over 200, with nothing in the table saying so.
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
    # The ranker takes bigger as better, and L1 is the other way round. 2 is the largest an L1
    # between two compositions can be, so this maps [0, 2] onto [1, 0] without changing the order.
    if "abundance_l1" in metrics:
        metrics["abundance_l1_similarity"] = 1.0 - metrics["abundance_l1"] / 2.0
    if not metrics:
        raise AbundanceError(
            "Nothing could be scored: no sample in %s had two bins the truth also knows" % table
        )
    return metrics
