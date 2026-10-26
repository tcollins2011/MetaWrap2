"""Tests for the CheckM1 retry after a pplacer crash.

With seven pplacer threads, pplacer can fail with
`Uncaught exception: Dune__exe__Multiprocessing.Child_error(_)`, which gives no cause; CheckM then
raises `RuntimeError: CheckM did not produce .../storage/bin_stats_ext.tsv`. This happened in three
of six `bin_refinement` steps on CAMI2 skin samples. The same input succeeds with one pplacer
thread, so CheckM is retried on one thread and then against the reduced tree.
"""

import os

from metawrap2 import checkm


class _Call:
    def __init__(self, cmd, kwargs):
        self.cmd = cmd
        self.kwargs = kwargs

    @property
    def pplacer_threads(self):
        parts = self.cmd.split()
        return int(parts[parts.index("--pplacer_threads") + 1])


def _fake_run(calls, produce_on):
    """A `run` that writes CheckM's marker table only on the attempt numbered *produce_on*."""

    def run(cmd, **kwargs):
        calls.append(_Call(cmd, kwargs))
        if len(calls) == produce_on:
            out = cmd.split()[5]  # checkm lineage_wf -x fasta <bins> <out>
            storage = os.path.join(out, "storage")
            os.makedirs(storage, exist_ok=True)
            with open(os.path.join(storage, "bin_stats_ext.tsv"), "w") as handle:
                handle.write(
                    "bin.1\t{'marker lineage': 'k__Bacteria', 'Completeness': 90.0, "
                    "'Contamination': 1.0, 'GC': 0.42, 'N50 (contigs)': 21000, "
                    "'Genome size': 3100000}\n"
                )
        return 0

    return run


def _bins(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin.1.fasta").write_text(">c1\nACGT\n")
    return str(bins)


def test_the_reduced_tree_is_the_last_resort_and_only_that(tmp_path, monkeypatch):
    """The reduced tree is the third attempt, never the first.

    One 2,022-residue alignment crashed pplacer at one thread in three of three tries and placed
    against the reduced tree. The reduced tree can choose a different lineage-specific marker set,
    which changes completeness, so it is used only as the last resort.
    """
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=3))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    stats = checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")
    assert len(calls) == 3
    assert [("--reduced_tree" in c.cmd) for c in calls] == [False, False, True]
    assert [c.pplacer_threads for c in calls] == [7, 1, 1]
    assert os.path.isfile(stats)


def test_the_reduced_tree_from_the_start_has_no_third_attempt(tmp_path, monkeypatch):
    """With ``reduced_tree`` (bin_refinement's default under --checkm1) there is no fallback."""
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=0))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    try:
        checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e", reduced_tree=True)
    except RuntimeError:
        pass
    assert len(calls) == 2
    assert all("--reduced_tree" in c.cmd for c in calls)


def test_a_crashed_first_attempt_is_retried_on_one_pplacer_thread(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=2))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)

    stats = checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")

    assert len(calls) == 2, "the first attempt produced nothing and must be retried"
    # 300 GB allows seven pplacer threads; the retry uses one.
    assert calls[0].pplacer_threads == checkm.pplacer_threads(60, 300) == 7
    assert calls[1].pplacer_threads == 1
    # The rest of the command is unchanged.
    assert calls[0].cmd.replace("--pplacer_threads 7", "--pplacer_threads 1").split()[:6] == (
        calls[1].cmd.split()[:6]
    )
    assert os.path.isfile(stats)


def test_the_first_attempt_is_not_fatal_but_the_retry_is(tmp_path, monkeypatch):
    """The first attempt runs with check=False so it can be retried; the second with check=True."""
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=2))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")
    assert calls[0].kwargs["check"] is False
    assert calls[1].kwargs["check"] is True


def test_a_working_checkm_is_run_once(tmp_path, monkeypatch):
    """A successful first attempt runs once, at the computed pplacer thread count."""
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=1))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    checkm.run_checkm(_bins(tmp_path), threads=8, mem_gb=64, env="e")
    assert len(calls) == 1
    assert calls[0].pplacer_threads == 1  # 64 // 40


def test_both_attempts_failing_still_raises_naming_the_file(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=0))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    try:
        checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")
    except RuntimeError as problem:
        assert "bin_stats_ext.tsv" in str(problem)
    else:
        raise AssertionError("a CheckM that never produced its table must raise")
    # Three attempts: as configured, then one pplacer thread, then the reduced reference tree.
    assert len(calls) == 3
    assert "--reduced_tree" in calls[2].cmd and "--reduced_tree" not in calls[0].cmd


def test_the_retry_clears_what_the_crash_left_behind(tmp_path, monkeypatch):
    """The retry removes the crashed attempt's output, since pplacer would resume from it."""
    calls = []
    bins = _bins(tmp_path)
    leftover = bins + ".checkm"
    os.makedirs(os.path.join(leftover, "storage"), exist_ok=True)
    with open(os.path.join(leftover, "storage", "half_written.pkl"), "w") as handle:
        handle.write("x")
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=2))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    checkm.run_checkm(bins, threads=60, mem_gb=300, env="e")
    assert not os.path.exists(os.path.join(leftover, "storage", "half_written.pkl"))


def test_pplacer_threads_is_capped_at_eight():
    """pplacer gains little past 8 threads, and its placements do not depend on the count."""
    assert checkm.pplacer_threads(24, 700) == checkm.PPLACER_USEFUL_THREADS == 8
    assert checkm.pplacer_threads(64, 4000) == 8
    # Below the cap, memory bounds it first, then the thread count.
    assert checkm.pplacer_threads(60, 300) == 7
    assert checkm.pplacer_threads(8, 64) == 1
    assert checkm.pplacer_threads(2, 700) == 2
    assert checkm.pplacer_threads(1, 1) == 1


def test_the_reduced_tree_budgets_16_gb_per_pplacer_thread():
    assert checkm.pplacer_threads(60, 300, reduced_tree=True) == 8
    assert checkm.pplacer_threads(8, 64, reduced_tree=True) == 4
    assert checkm.pplacer_threads(8, 20, reduced_tree=True) == 1
