"""How much of an assembly reached a bin, in the unit that answers the question.

A bin set reports what it found and not what it missed. The number this adds is cheap and changes
what a reader concludes: across 92,187 metagenomes only ~10% of assembled bp reach a medium- or
high-quality MAG, and the unbinned remainder is biased rather than random.

Two things the tests are really about:

* **Base pairs, not contigs.** Binners discard short contigs first and short contigs are numerous,
  so the two numbers differ by a lot and only one of them answers "how much of the sample is this".
  Measured on a real CAMI2 refined bin set: 17.8% of contigs, 58.7% of base pairs - reporting
  contigs would understate recovery 3.3x.
* **Not-measured is not zero.** With no quality-filtered set to compare, the clean fraction is
  absent rather than 0.0%, because "I did not look" and "none of it is clean" read the same way at a
  glance and are opposite findings.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import recovery


def _fasta(path, records):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, seq) for name, seq in records)
    return str(path)


SUMMARY = os.path.join(os.path.dirname(__file__), "fixtures", "appraise", "appraise_summary.tsv")


@pytest.fixture
def assembly(tmp_path):
    # One long contig and four short ones: the shape that makes bp and contig counts disagree.
    return _fasta(
        tmp_path / "assembly.fa",
        [("long", "A" * 10_000)] + [("short%d" % i, "C" * 250) for i in range(4)],
    )


def test_the_fraction_is_of_base_pairs_not_contigs(assembly, tmp_path):
    """The long contig is 20% of the contigs and 90.9% of the sequence."""
    bins = tmp_path / "bins"
    bins.mkdir()
    _fasta(bins / "bin.1.fa", [("long", "A" * 10_000)])
    result = recovery.measure(assembly, str(bins))
    assert result.binned_contigs == 1
    assert result.assembly_contigs == 5
    assert result.binned_fraction == pytest.approx(10_000 / 11_000)
    assert result.binned_fraction > 0.9, "a contig count would have said 20%"
    # Lengths come from the assembly, not the bins, so a binner that trimmed a contig cannot change
    # the denominator; a contig in a bin that is not in the assembly is ignored; and a contig in
    # two bins is still one contig of the assembly.
    _fasta(bins / "bin.1.fa", [("long", "A" * 50), ("not_in_assembly", "T" * 99_999)])
    _fasta(bins / "bin.2.fa", [("long", "A" * 10_000)])
    result = recovery.measure(assembly, str(bins))
    assert result.assembly_bp == 11_000
    assert result.binned_bp == 10_000, "the assembly's length for 'long', not the bin's"
    assert result.binned_contigs == 1, "the foreign contig is not counted, nor 'long' twice"
    assert result.unbinned_bp == 1_000


# ── not measured versus zero ─────────────────────────────────────────────────────────────────


def test_the_clean_fraction_is_absent_when_no_clean_set_was_given(assembly, tmp_path):
    """ "I did not look" and "none of it is clean" read the same at a glance and are opposites."""
    bins = tmp_path / "bins"
    bins.mkdir()
    _fasta(bins / "bin.1.fa", [("long", "A" * 10_000)])
    result = recovery.measure(assembly, str(bins))
    assert result.clean_bp is None
    assert result.clean_fraction is None
    assert "quality thresholds" not in result.describe()
    row = recovery.recovery_row(result)
    assert row["clean_binned_fraction_bp"] == ""

    # The other side of the distinction: looked, found none.
    clean = tmp_path / "clean"
    clean.mkdir()
    result = recovery.measure(assembly, str(bins), str(clean))
    assert result.clean_bp == 0
    assert result.clean_fraction == 0.0
    assert recovery.recovery_row(result)["clean_binned_fraction_bp"] == "0.0000"

    # The same for SingleM's appraisal: `coverage_assembled` appears only when an assembly was
    # appraised too, so a missing column is None (no assembly given), not zero (accounted for
    # nothing). Its `total` and `average` rows are not samples; averaging them in double-counts.
    rows = recovery.parse_appraisal(SUMMARY)
    assert {row.sample for row in rows} == {"ERR011347_1"}
    bacteria = next(row for row in rows if row.domain == "d__Bacteria")
    assert bacteria.percent_binned == pytest.approx(11.56, abs=0.01)
    assert all(row.assembled is None for row in rows)
    assert all(recovery.appraisal_row(row)["coverage_assembled_not_binned"] == "" for row in rows)
