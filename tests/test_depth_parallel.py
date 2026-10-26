"""The depth phase runs several pileups at once.

A user on a large machine measured `coverage` using every core while aligning, then processing BAMs
one at a time with `samtools view | samtools depth` at ~2 minutes each and ~25% utilisation (6h44m
for 36 water samples, 2h14m for 50 gut samples). `samtools depth` does not scale with threads, but
BAMs can run in parallel since each fills its own column of the running sums.
"""

import numpy as np

from metawrap2 import depth


def test_workers_divide_the_budget_rather_than_giving_one_pileup_everything():
    # 128 threads, 50 BAMs: capped by MAX_DEPTH_WORKERS, not by the thread count.
    assert depth.depth_workers(50, 128) == depth.MAX_DEPTH_WORKERS
    # 16 threads, 8 BAMs: four threads each.
    assert depth.depth_workers(8, 16) == 4
    # More BAMs than threads allow: one worker per DEPTH_THREADS_PER_BAM.
    assert depth.depth_workers(36, 12) == 3


def test_a_small_budget_or_one_bam_stays_serial():
    """A small budget or a single BAM runs serially."""
    assert depth.depth_workers(1, 128) == 1
    assert depth.depth_workers(36, 4) == 1
    assert depth.depth_workers(36, 1) == 1


def test_every_bam_lands_in_its_own_column(monkeypatch, tmp_path):
    """Each BAM writes only its own column of the output table."""
    seen = []

    def fake_accumulate(
        bam, column, lengths, index, edge, identity, mapq, threads, env, total, total_sq
    ):
        seen.append((bam, column, threads))
        total[:, column] += float(column + 1)
        total_sq[:, column] += float(column + 1)

    monkeypatch.setattr(depth, "_accumulate_bam", fake_accumulate)
    monkeypatch.setattr(depth, "contig_lengths", lambda bam, env=None: {"k141_1": 1000})
    out = tmp_path / "depth.txt"
    depth.write_depth_table(
        ["a.bam", "b.bam", "c.bam"], str(out), samples=["a", "b", "c"], threads=16
    )

    assert sorted(column for _bam, column, _t in seen) == [0, 1, 2]
    # Three BAMs share a 16-thread budget, five threads each. DEPTH_THREADS_PER_BAM sets how many
    # run at once, not how many threads each gets.
    assert {t for _bam, _c, t in seen} == {5}
    header = out.read_text().splitlines()[0].split("\t")
    assert header[:3] == ["contigName", "contigLen", "totalAvgDepth"]
    row = out.read_text().splitlines()[1].split("\t")
    # Column order follows the BAM list.
    assert float(row[3]) < float(row[5]) < float(row[7])


def test_a_failing_bam_stops_the_table_rather_than_writing_a_partial_one(monkeypatch, tmp_path):
    def explode(bam, column, *args, **kwargs):
        if column == 1:
            raise RuntimeError("samtools depth failed on %s" % bam)

    monkeypatch.setattr(depth, "_accumulate_bam", explode)
    monkeypatch.setattr(depth, "contig_lengths", lambda bam, env=None: {"k141_1": 1000})
    out = tmp_path / "depth.txt"
    try:
        depth.write_depth_table(["a.bam", "b.bam"], str(out), threads=16)
    except RuntimeError as problem:
        assert "b.bam" in str(problem)
    else:
        raise AssertionError("the failure has to reach the caller")
    assert not out.exists(), "a table measured from a failed pileup must not be written"


def test_the_numbers_do_not_depend_on_how_many_workers_ran(monkeypatch, tmp_path):
    """The table is the same for any number of workers."""
    monkeypatch.setattr(depth, "contig_lengths", lambda bam, env=None: {"k141_1": 1000})

    def fake_accumulate(
        bam, column, lengths, index, edge, identity, mapq, threads, env, total, total_sq
    ):
        total[:, column] += 100.0 * (column + 1)
        total_sq[:, column] += 20000.0 * (column + 1)

    monkeypatch.setattr(depth, "_accumulate_bam", fake_accumulate)
    bams = ["a.bam", "b.bam", "c.bam", "d.bam"]
    serial = tmp_path / "serial.txt"
    parallel = tmp_path / "parallel.txt"
    depth.write_depth_table(bams, str(serial), threads=1)
    depth.write_depth_table(bams, str(parallel), threads=64)
    assert np.array_equal(
        [l.split("\t") for l in serial.read_text().splitlines()],
        [l.split("\t") for l in parallel.read_text().splitlines()],
    )
