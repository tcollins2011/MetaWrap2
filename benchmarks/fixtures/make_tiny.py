#!/usr/bin/env python3
"""Regenerate the committed ``tiny`` fixture.

The fixture is committed rather than generated at test time so that the tests measure the harness
against fixed bytes. This script exists so those bytes are reproducible: the seed is fixed, so
re-running it produces the same files, and a results diff is never noise from a regenerated
fixture.

    python benchmarks/fixtures/make_tiny.py

Three synthetic genomes with deliberately different GC content, shredded into ~250 contigs, plus a
small paired read set sampled from them and a gold-standard contig-to-genome mapping in the same
two-column format CAMI uses. Far too small to say anything about a binner; large enough to prove
the harness computes and records what it claims to.
"""

from __future__ import annotations

import gzip
import os
import random

#: Fixed so the fixture is reproducible. Changing it changes every committed byte, which is a
#: bigger decision than it looks.
SEED = 20260929

#: Genome name -> GC fraction. Different GC is what gives a binner a signal to find at this scale;
#: coverage alone would not, with one sample.
GENOMES = {"genome_A": 0.30, "genome_B": 0.50, "genome_C": 0.68}

GENOME_LENGTH = 120_000
CONTIG_RANGE = (1000, 1800)
READ_PAIRS = 4000
READ_LENGTH = 150
LINE_WIDTH = 70


def main() -> None:
    root = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(root, "tiny")
    os.makedirs(os.path.join(out, "reads"), exist_ok=True)
    rng = random.Random(SEED)

    def sequence(length: int, gc: float) -> str:
        return "".join(
            rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(length)
        )

    contigs = []  # (contig id, genome, sequence)
    for name, gc in GENOMES.items():
        genome = sequence(GENOME_LENGTH, gc)
        position = index = 0
        while position < len(genome) - CONTIG_RANGE[1]:
            length = rng.randint(*CONTIG_RANGE)
            contigs.append(
                ("%s_contig_%03d" % (name, index), name, genome[position : position + length])
            )
            position += length
            index += 1

    with open(os.path.join(out, "assembly.fa"), "w") as fh:
        for contig_id, _genome, seq in contigs:
            fh.write(">%s\n" % contig_id)
            fh.writelines(
                seq[offset : offset + LINE_WIDTH] + "\n"
                for offset in range(0, len(seq), LINE_WIDTH)
            )

    # With the sequence lengths, because that is what AMBER weights its bp metrics by and it reads
    # them from the gold standard alone. Without the column it dies inside pandas; the harness can
    # add them back from the assembly, but a fixture may as well ship complete.
    with open(os.path.join(out, "gold_standard_binning.tsv"), "w") as fh:
        fh.write("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\t_LENGTH\n")
        fh.writelines(
            "%s\t%s\t%d\n" % (contig_id, genome, len(seq)) for contig_id, genome, seq in contigs
        )

    quality = "I" * READ_LENGTH
    r1_path = os.path.join(out, "reads", "sample_1_1.fq.gz")
    r2_path = os.path.join(out, "reads", "sample_1_2.fq.gz")
    with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
        for n in range(READ_PAIRS):
            _contig_id, _genome, seq = contigs[rng.randrange(len(contigs))]
            if len(seq) < 400:
                continue
            start = rng.randrange(0, len(seq) - 350)
            forward = seq[start : start + READ_LENGTH]
            reverse = seq[start + 200 : start + 200 + READ_LENGTH]
            reverse = reverse.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            r1.write("@read%05d/1\n%s\n+\n%s\n" % (n, forward, quality))
            r2.write("@read%05d/2\n%s\n+\n%s\n" % (n, reverse, quality))

    print("wrote %d contigs from %d genomes into %s" % (len(contigs), len(GENOMES), out))


if __name__ == "__main__":
    main()
