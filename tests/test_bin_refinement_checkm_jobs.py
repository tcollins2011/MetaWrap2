"""Scoring the hybrid bin sets in parallel (`--checkm-jobs`).

The sets are independent and share only the bin-score cache, which holds a lock. These tests pin
how many sets run at once, how bins are assigned to sets, and that a failure in one set stops the
run.
"""

import concurrent.futures
import json
import os
import threading

import pytest

from metawrap2 import checkm as _checkm
from metawrap2.modules import bin_refinement


def _args(*extra):
    return bin_refinement._parse_args(["-A", "a", "-B", "b", "-o", "o"] + list(extra))


def test_memory_bounds_concurrency_for_checkm1():
    """16 GB per CheckM1 call on the reduced tree (the default), 40 GB with --full-tree."""
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "120", "--checkm1"), 7) == 7
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "64", "--checkm1"), 7) == 4
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "8", "--checkm1"), 7) == 1
    full = ("--checkm1", "--full-tree")
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "120", *full), 7) == 3
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "700", *full), 7) == 7
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "64", *full), 7) == 1


def test_threads_also_bound_it():
    """At least one thread per call."""
    assert bin_refinement.checkm_jobs(_args("-t", "4", "-m", "700"), 7) == 4
    assert bin_refinement.checkm_jobs(_args("-t", "1", "-m", "700"), 7) == 1


def test_it_never_exceeds_the_number_of_sets():
    assert bin_refinement.checkm_jobs(_args("-t", "64", "-m", "700"), 3) == 3
    assert bin_refinement.checkm_jobs(_args("-t", "64", "-m", "700"), 1) == 1


def test_checkm2_is_budgeted_separately():
    """24 GB per CheckM2 call (the default); --full-tree does not change it."""
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "120"), 7) == 5
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "120", "--checkm2"), 7) == 5
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "120", "--full-tree"), 7) == 5


def test_the_flag_overrides_the_budget():
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "700", "--checkm-jobs", "1"), 7) == 1
    assert bin_refinement.checkm_jobs(_args("-t", "24", "-m", "64", "--checkm-jobs", "4"), 7) == 4
    assert bin_refinement.checkm_jobs(_args("--checkm-jobs", "9"), 2) == 2


def test_checkm2_is_the_default_and_checkm1_selects_checkm1():
    """CheckM2 scores by default; --checkm1 selects CheckM1, and --checkm2 is kept as a no-op."""
    assert _args().checkm2 is True
    assert _args("--checkm2").checkm2 is True
    assert _args("--checkm1").checkm2 is False
    assert bin_refinement.scorer_name(_args()) == "checkm2"
    assert bin_refinement.scorer_name(_args("--checkm1")) == "checkm1_reduced"
    assert bin_refinement.scorer_name(_args("--checkm1", "--full-tree")) == "checkm1"


def test_quick_is_hidden_and_changes_nothing(capsys):
    """--quick added the reduced tree before it became CheckM1's default; old scripts still run."""
    from metawrap2.modules._common import scorer_notes

    assert _args("--checkm1", "--quick").full_tree is False
    with pytest.raises(SystemExit):
        _args("-h")
    assert "--quick" not in capsys.readouterr().out
    scorer_notes(_args("--checkm1", "--quick"))
    captured = capsys.readouterr()
    assert "reduced tree by default" in captured.out + captured.err


def test_checkm1_and_checkm2_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        _args("--checkm1", "--checkm2")


def test_checkm2_ttable_is_opt_in(monkeypatch, tmp_path):
    assert _args().checkm2_ttable is None
    assert _args("--checkm2-ttable", "11").checkm2_ttable == 11
    with pytest.raises(SystemExit):
        _args("--checkm2-ttable", "5")

    seen = []
    monkeypatch.setattr(_checkm, "run", lambda cmd, **kw: seen.append(cmd))
    monkeypatch.setattr(_checkm, "summarize_checkm2", lambda *a, **kw: [])
    from metawrap2.command import runner

    monkeypatch.setattr(runner, "dry_run", True)
    _checkm.run_checkm2(str(tmp_path), threads=2, env="e", database="")
    _checkm.run_checkm2(str(tmp_path), threads=2, env="e", database="", ttable=4)
    assert "--ttable" not in seen[0]
    assert seen[1].endswith("--ttable 4")


def _bin_sets(root, layout):
    """Real bin folders, because the planning pass reads the bins to key and size them."""
    for name, bins in layout.items():
        directory = os.path.join(str(root), name)
        os.makedirs(directory, exist_ok=True)
        for bin_name, contigs in bins.items():
            with open(os.path.join(directory, bin_name + ".fasta"), "w") as handle:
                handle.writelines(">%s\n%s\n" % (contig, "ACGT" * 75_000) for contig in contigs)
    return list(layout)


def test_every_set_is_scored_once_and_the_resources_are_divided(monkeypatch, tmp_path):
    """Each set scores its own bins; a bin two sets share is scored by the first set only."""
    sets = _bin_sets(
        tmp_path,
        {
            "binsA": {"bin.1": ["c1", "c2"], "bin.2": ["c3"]},
            "binsB": {"bin.1": ["c4"], "bin.2": ["c3"]},  # binsB's bin.2 has binsA's bin.2 contigs
            "binsC": {"bin.1": ["c5"]},
        },
    )
    seen = {}
    lock = threading.Lock()

    def fake_score(directory, bins, cache, **kwargs):
        with lock:
            seen[os.path.basename(directory)] = (
                sorted(os.path.basename(b) for b in bins),
                kwargs.get("threads"),
                kwargs.get("mem_gb"),
            )
        for path in bins:  # record scores so the write pass finds them
            with cache.lock:
                cache.rows[_checkm.bin_identity(path)] = [
                    os.path.splitext(os.path.basename(path))[0],
                    "91.0",
                    "1.0",
                    "0.42",
                    "root",
                    "21000",
                    "3100000",
                ]

    monkeypatch.setattr(_checkm, "score_bins", fake_score)
    monkeypatch.setattr(bin_refinement, "_run_checkm_set", lambda *a, **kw: None)

    bin_refinement._score_bin_sets(
        str(tmp_path),
        sets,
        _args("-t", "12", "-m", "700", "--score-everything"),
        "env",
        None,
        _checkm.BinScoreCache(),
        jobs=3,
        per_threads=4,
        per_memory=233,
    )

    assert sorted(seen) == ["binsA", "binsB", "binsC"]
    assert seen["binsA"][0] == ["bin.1.fasta", "bin.2.fasta"]
    assert seen["binsB"][0] == ["bin.1.fasta"]
    assert {(t, m) for _b, t, m in seen.values()} == {(4, 233)}


def test_the_plan_does_not_depend_on_finishing_order(tmp_path):
    sets = _bin_sets(
        tmp_path,
        {
            "binsA": {"bin.1": ["shared"], "bin.2": ["a2"]},
            "binsB": {"bin.1": ["shared"], "bin.2": ["b2"]},
        },
    )
    dirs = [os.path.join(str(tmp_path), name) for name in sets]
    plan = _checkm.plan_scoring(dirs, _checkm.BinScoreCache(), prefilter=False)
    assert [os.path.basename(p) for p in plan[dirs[0]]] == ["bin.1.fasta", "bin.2.fasta"]
    assert [os.path.basename(p) for p in plan[dirs[1]]] == ["bin.2.fasta"]
    assert _checkm.plan_scoring(dirs, _checkm.BinScoreCache(), prefilter=False) == plan


def test_a_bin_already_in_the_cache_is_not_planned(tmp_path):
    sets = _bin_sets(tmp_path, {"binsA": {"bin.1": ["c1"], "bin.2": ["c2"]}})
    directory = os.path.join(str(tmp_path), sets[0])
    cache = _checkm.BinScoreCache()
    known = os.path.join(directory, "bin.1.fasta")
    cache.rows[_checkm.bin_identity(known)] = ["bin.1", "90.0", "1.0", "0.4", "root", "2", "3"]
    plan = _checkm.plan_scoring([directory], cache, prefilter=False)
    assert [os.path.basename(p) for p in plan[directory]] == ["bin.2.fasta"]


def test_one_job_scores_the_sets_in_turn(monkeypatch, tmp_path):
    """jobs=1 runs one call per set with the full thread and memory budget, and no planning pass."""
    seen = []
    monkeypatch.setattr(
        bin_refinement,
        "_run_checkm_set",
        lambda d, a, e, c=None, **kw: seen.append((os.path.basename(d), kw.get("threads"))),
    )
    monkeypatch.setattr(
        _checkm, "plan_scoring", lambda *a, **kw: pytest.fail("no planning with one job")
    )
    bin_refinement._score_bin_sets(
        str(tmp_path),
        ["binsA", "binsB"],
        _args("-t", "24", "-m", "120"),
        "env",
        None,
        _checkm.BinScoreCache(),
        jobs=1,
        per_threads=24,
        per_memory=120,
    )
    assert seen == [("binsA", None), ("binsB", None)]


def test_a_failed_set_stops_the_run(monkeypatch, tmp_path):
    """Consolidation compares the sets, so a set that failed to score must not be left out."""
    sets = _bin_sets(
        tmp_path,
        {"binsA": {"bin.1": ["c1"]}, "binsB": {"bin.1": ["c2"]}, "binsC": {"bin.1": ["c3"]}},
    )

    def fake_score(directory, bins, cache, **kwargs):
        if directory.endswith("binsB"):
            raise RuntimeError("CheckM did not produce bin_stats_ext.tsv")

    monkeypatch.setattr(_checkm, "score_bins", fake_score)
    monkeypatch.setattr(bin_refinement, "_run_checkm_set", lambda *a, **kw: None)

    with pytest.raises(RuntimeError, match="bin_stats_ext"):
        bin_refinement._score_bin_sets(
            str(tmp_path),
            sets,
            _args("-t", "8", "-m", "700", "--score-everything"),
            "env",
            None,
            _checkm.BinScoreCache(),
            jobs=3,
            per_threads=2,
            per_memory=233,
        )


def test_the_score_cache_survives_concurrent_use(tmp_path):
    """No rows or counts are lost when threads record scores at the same time."""
    cache = _checkm.BinScoreCache(str(tmp_path / "scores.json"))
    rows = ["bin_%03d" % i for i in range(400)]

    def record(name):
        with cache.lock:
            cache.rows[name] = [name, "90.0", "1.0", "0.42", "checkm2", "21000", "3100000"]
            cache.hits += 1
        cache.save()

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(record, rows))

    assert len(cache.rows) == 400 and cache.hits == 400
    with open(cache.path) as handle:
        assert len(json.load(handle)) == 400
    assert [f for f in os.listdir(tmp_path) if f.endswith(".tmp")] == []
