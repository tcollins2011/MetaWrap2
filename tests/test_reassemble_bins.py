import os

import pytest

from metawrap2.modules import reassemble_bins


def _finished_run(out, winners=("bin.1.strict", "bin.2.orig")):
    """Create a finished reassemble_bins layout: chosen bins promoted, candidates in work_files."""
    os.makedirs(os.path.join(out, "reassembled_bins"))
    for name in winners:
        with open(os.path.join(out, "reassembled_bins", name + ".fasta"), "w") as fh:
            fh.write(">old\nACGT\n")
    with open(os.path.join(out, "reassembled_bins" + ".stats.tsv"), "w") as fh:
        fh.write("bin\tcompleteness\n")
    for item in (
        "reassembled_bins.checkm",
        "reads_for_reassembly",
        "binned_assembly",
        "reassemblies",
    ):
        os.makedirs(os.path.join(out, "work_files", item))
        with open(os.path.join(out, "work_files", item, "leftover"), "w") as fh:
            fh.write("from the first run\n")


def test_the_shuffle_survives_a_destination_that_already_exists(tmp_path):
    """Moving intermediates into work_files works when a previous run's copies are there.

    ``os.replace`` refuses to replace a non-empty directory, and each ``work_files`` destination may
    hold the previous run's copy. A failure here would come after all reassemblies had finished.
    """
    out = str(tmp_path / "out")
    _finished_run(out)
    # this run's candidates and winners, where the step expects them
    for item in ("reassembled_bins", "reads_for_reassembly", "binned_assembly", "reassemblies"):
        os.makedirs(os.path.join(out, item), exist_ok=True)
        with open(os.path.join(out, item, "new"), "w") as fh:
            fh.write("from this run\n")
    best = os.path.join(out, "reassembled_best_bins")
    os.makedirs(best)
    with open(os.path.join(best, "bin.1.permissive.fasta"), "w") as fh:
        fh.write(">new\nACGTACGT\n")

    reassembled = os.path.join(out, "reassembled_bins")
    reassemble_bins._shuffle_intermediates(out, best, reassembled)

    assert sorted(os.listdir(reassembled)) == ["bin.1.permissive.fasta"]
    # and work_files holds this run's candidates, not the previous run's leftovers
    moved = os.path.join(out, "work_files", "reassemblies")
    assert sorted(os.listdir(moved)) == ["new"], os.listdir(moved)
    assert not os.path.exists(best)


def test_force_clears_a_finished_reassembly(tmp_path):
    """--force removes the previous run's results so they are not chosen again.

    reassembled_bins is added to, not replaced, so old promoted bins and .orig.fa copies would
    otherwise be scored with this run's candidates.
    """
    out = str(tmp_path / "out")
    _finished_run(out)
    # The run record, log and original bins are kept.
    os.makedirs(os.path.join(out, "original_bins"))
    for name in ("run_config.json", "metawrap2.log", "provenance.txt"):
        with open(os.path.join(out, name), "w") as fh:
            fh.write("{}\n")
    reassemble_bins._clean_previous_run(out)
    assert not os.path.exists(os.path.join(out, "reassembled_bins"))
    assert not os.path.exists(os.path.join(out, "work_files"))
    assert not os.path.exists(os.path.join(out, "reassembled_bins.stats.tsv"))
    kept = sorted(os.listdir(out))
    assert kept == ["metawrap2.log", "original_bins", "provenance.txt", "run_config.json"], kept


def _scored_stage(tmp_path, monkeypatch, *flags):
    """Run the CheckM stage as a dry run; return the scorer keywords and the commands issued."""
    args = reassemble_bins._parse_args(
        ["-b", "b", "-o", str(tmp_path), "-1", "r1", "-2", "r2"] + list(flags)
    )
    calls, commands = [], []
    monkeypatch.setattr(reassemble_bins, "dry_run", lambda: True)
    monkeypatch.setattr(reassemble_bins, "run", lambda cmd, **kw: commands.append(cmd) or 0)
    monkeypatch.setattr(
        reassemble_bins._checkm,
        "run_checkm",
        lambda bins_dir, **kw: calls.append(kw) or bins_dir + ".stats.tsv",
    )
    reassemble_bins._run_checkm_stage(args, str(tmp_path), "checkm2-env")
    return calls, commands


def test_reassembled_bins_are_scored_with_checkm2(tmp_path, monkeypatch):
    calls, commands = _scored_stage(tmp_path, monkeypatch)
    assert calls and all(kw["checkm2"] is True for kw in calls)
    assert calls[0]["checkm2_env"] == "checkm2-env"
    assert not any("bin_qa_plot" in cmd for cmd in commands)
    assert any("plot_reassembly" in cmd for cmd in commands)


def test_checkm2_is_still_accepted_and_changes_nothing(tmp_path, monkeypatch):
    calls, _ = _scored_stage(tmp_path, monkeypatch, "--checkm2")
    assert calls and all(kw["checkm2"] is True for kw in calls)


def test_checkm1_is_refused_with_one_sentence(capsys):
    reassemble = ["-b", "b", "-o", "o", "-1", "r1", "-2", "r2"]
    with pytest.raises(SystemExit):
        reassemble_bins._parse_args(reassemble + ["--checkm1"])
    err = capsys.readouterr().err
    assert "reassemble_bins scores bins with CheckM2 only" in err


def test_every_module_that_scores_bins_defaults_to_checkm2():
    """binning --run-checkm and bin_refinement take --checkm1; reassemble_bins is CheckM2 only."""
    from metawrap2.modules import bin_refinement, binning

    base = ["-a", "a.fa", "-o", "out", "--metabat2", "r_1.fastq"]
    assert binning._parse_args(base + ["--run-checkm"]).checkm2 is True
    assert binning._parse_args(base + ["--run-checkm", "--checkm1"]).checkm2 is False
    refine = ["-A", "a", "-B", "b", "-o", "o"]
    assert bin_refinement._parse_args(refine).checkm2 is True
    reassemble = ["-b", "b", "-o", "o", "-1", "r1", "-2", "r2"]
    assert reassemble_bins._parse_args(reassemble).checkm2 is True
