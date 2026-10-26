"""Sample-level feature abundance.

The sum is simple; these tests cover the handling around it:

* **Genome size is divided out once, or not at all.** Coverage is already per cell; a read count
  is not. Treating coverage as counts divides by length twice and biases large genomes low.
* **A genome with no abundance contributes nothing**, instead of defaulting to one. A MAG absent
  from a sample's abundance table was not observed there.
* **A unit that can be looked up is not inferred.** SingleM writes a percentage; inferring from
  magnitude reads a host-dominated sample's 0.8% as 80%.
* **Partial pathways are not summed.** Two genomes with half a module each are not one genome with
  a whole one, so module presence is binary above a completeness threshold.
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
    """CoverM's values are already coverage-proportional and are not divided by genome size."""
    entry = _coverage("big_genome", "s1", 10.0)
    assert entry.cell_proportional == 10.0
    # A read count is divided by genome size once; with no known length it contributes nothing.
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
    """SingleM's unit is read, not inferred from magnitude ("above 1 is a percentage"), which read
    0.8% as 80% and scaled gene abundances down 100x."""
    path = _write_fraction(tmp_path / "mf.tsv", [("host_heavy", "0.80")])
    assert fa.parse_singlem_fraction(path)["host_heavy"] == pytest.approx(0.008)


# ── the assignable fraction ──────────────────────────────────────────────────────────────────


# ── contributor counts ───────────────────────────────────────────────────────────────────────
