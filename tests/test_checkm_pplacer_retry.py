"""CheckM1's parallel tree placement can die without saying why, so it gets one retry.

What happened: on three CAMI2 skin samples, three of six `bin_refinement` steps ended with
`RuntimeError: CheckM did not produce .../storage/bin_stats_ext.tsv`, each after CheckM had logged
every marker gene aligned. Underneath, pplacer had crashed as
`Uncaught exception: Dune__exe__Multiprocessing.Child_error(_)` - an OCaml child failure that
carries no cause - at seven pplacer threads. The other three steps, over
the same kind of bin set on the same host, finished. CheckM is effectively the whole cost of
refinement, so an hour of it was being thrown away by a crash that a single thread survives.
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
    """One ordinary 2,022-residue alignment killed pplacer at one thread, three times of three,
    and placed fine against CheckM's reduced tree. So the reduced tree is tried - third, never
    first: it can choose a different lineage-specific marker set, which moves completeness.
    """
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=3))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    stats = checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")
    assert len(calls) == 3
    assert [("--reduced_tree" in c.cmd) for c in calls] == [False, False, True]
    assert [c.pplacer_threads for c in calls] == [7, 1, 1]
    assert os.path.isfile(stats)


def test_quick_runs_the_reduced_tree_from_the_start_and_does_not_fall_back_to_it(
    tmp_path, monkeypatch
):
    """`--quick` already is the reduced tree, so there is no third attempt to make."""
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=0))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    try:
        checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e", quick=True)
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
    # Seven threads is what 300 GB allows, and one is what the retry is for.
    assert calls[0].pplacer_threads == checkm.pplacer_threads(60, 300) == 7
    assert calls[1].pplacer_threads == 1
    # Everything else is the same command, so the retry cannot drift from the first attempt.
    assert calls[0].cmd.replace("--pplacer_threads 7", "--pplacer_threads 1").split()[:6] == (
        calls[1].cmd.split()[:6]
    )
    assert os.path.isfile(stats)


def test_the_first_attempt_is_not_fatal_but_the_retry_is(tmp_path, monkeypatch):
    """The first call must come back so the retry can happen; the second must not be swallowed."""
    calls = []
    monkeypatch.setattr(checkm, "run", _fake_run(calls, produce_on=2))
    monkeypatch.setattr(checkm, "checkm1_data_env", dict)
    checkm.run_checkm(_bins(tmp_path), threads=60, mem_gb=300, env="e")
    assert calls[0].kwargs["check"] is False
    assert calls[1].kwargs["check"] is True


def test_a_working_checkm_is_run_once(tmp_path, monkeypatch):
    """No retry, and no change to the common path: one invocation, at the computed thread count."""
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
    # Three attempts: as asked, then one pplacer thread, then CheckM's reduced reference tree.
    assert len(calls) == 3
    assert "--reduced_tree" in calls[2].cmd and "--reduced_tree" not in calls[0].cmd


def test_the_retry_clears_what_the_crash_left_behind(tmp_path, monkeypatch):
    """pplacer resumes from its tmpdir, so a crashed tree left in place reproduces the crash."""
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
