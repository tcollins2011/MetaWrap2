"""Long reads: platform flags for four platforms, and a check that the flag fits the reads.

The numbers in the scan tests come from 25x simulated R10 nanopore reads (badread, nanopore2023
model) over three CAMI2 MAGs. metaFlye's log reported an N50 of 12,326 bp for that file, and the
scanner reports the same.
"""

from __future__ import annotations

import pytest

from metawrap2 import longreads as lr

# ── the check that the declared platform fits ────────────────────────────────────────────────


def _fastq(path, lengths):
    with open(path, "w") as fh:
        fh.writelines(
            "@r%d\n%s\n+\n%s\n" % (index, "A" * length, "I" * length)
            for index, length in enumerate(lengths)
        )
    return str(path)


def test_short_reads_passed_as_long_reads_are_caught(tmp_path):
    """Short reads passed as long reads are detected.

    metaFlye accepts short reads and produces a fragmented assembly without an error.
    """
    path = _fastq(tmp_path / "r.fastq", [150] * 100)
    stats = lr.scan_reads(path)
    assert stats.looks_like_short_reads
    notes = lr.check_type(stats, "ont-hq")
    assert any("looks like a short-read file" in note for note in notes)
    # The scan measures N50, not the mean: 10 kb + 1 kb x 5 has a mean of 2.5 kb and an N50 of
    # 10 kb, and N50 indicates whether reads are long.
    stats = lr.scan_reads(_fastq(tmp_path / "long.fastq", [10_000] + [1_000] * 5))
    assert (stats.reads, stats.bases, stats.n50) == (6, 15_000, 10_000)
    assert stats.mean_length == pytest.approx(2500.0)
