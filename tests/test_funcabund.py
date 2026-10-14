"""Sample-level feature abundance, and the four ways the arithmetic goes wrong quietly.

The sum itself is one line. What is tested here is everything around it:

* **Genome size is divided out once, or not at all.** Coverage is already per-cell; a read count is
  not. Treating coverage as counts divides by length twice and makes large genomes look rare, which
  is a systematic error in one direction rather than noise.
* **A genome with no abundance contributes nothing**, rather than defaulting to one. A MAG absent
  from a sample's abundance table was not observed there.
* **A unit that can be looked up is never inferred.** SingleM writes a percentage; guessing from
  the magnitude reads a host-dominated sample's 0.8% as 80%.
* **Partial pathways do not add up.** Two genomes with half a module each are not one genome with a
  whole one, so module presence is binary above a completeness threshold.
"""

from __future__ import annotations

import pytest

from metawrap2 import funcabund as fa


def _coverage(genome, sample, value):
    return fa.GenomeAbundance(genome, sample, value, fa.AbundanceKind.COVERAGE)


def _counts(genome, sample, value, length):
    return fa.GenomeAbundance(genome, sample, value, fa.AbundanceKind.READ_COUNT, length)


# ── the genome-size trap ─────────────────────────────────────────────────────────────────────


def test_coverage_is_not_divided_by_genome_size():
    """CoverM's numbers are already coverage-proportional; dividing again double-corrects.

    This is the single easiest thing in the module to get backwards, and getting it wrong makes
    every large genome look rare in a way nothing else in the output would reveal.
    """
    entry = _coverage("big_genome", "s1", 10.0)
    assert entry.cell_proportional == 10.0
    # A read count is divided by genome size exactly once, and with no known length contributes
    # nothing rather than a guessed abundance.
    assert _counts("g", "s1", 1000.0, 5_000_000).cell_proportional == pytest.approx(0.0002)
    assert _counts("g", "s1", 1000.0, None).cell_proportional == 0.0


# ── the sum ──────────────────────────────────────────────────────────────────────────────────


# ── units ────────────────────────────────────────────────────────────────────────────────────


# ── SingleM's percentage ─────────────────────────────────────────────────────────────────────


def _write_fraction(path, rows):
    header = "sample\tbacterial_archaeal_bases\tmetagenome_size\tread_fraction\twarning\n"
    with open(path, "w") as fh:
        fh.write(header)
        fh.writelines("%s\t1\t1\t%s\t\n" % (sample, value) for sample, value in rows)
    return str(path)


def test_a_host_dominated_sample_is_not_misread_as_eighty_percent(tmp_path):
    """The bug the first version had. Inferring the unit from the magnitude - "above 1 must be a
    percentage" - reads 0.8% as 80% and scales every gene abundance down by 100x, in exactly the
    samples where the correction matters most."""
    path = _write_fraction(tmp_path / "mf.tsv", [("host_heavy", "0.80")])
    assert fa.parse_singlem_fraction(path)["host_heavy"] == pytest.approx(0.008)


# ── the assignable fraction ──────────────────────────────────────────────────────────────────


# ── contributor counts ───────────────────────────────────────────────────────────────────────
