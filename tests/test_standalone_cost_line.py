"""What a module run on its own tells the user about itself.

A module invoked directly is a supported way to use MetaWrap2 - the docs say every module runs
on its own - but two things the ``run`` driver prints were missing from that path: where the
live log is, and what the run cost.
"""

from metawrap2 import command, usage
from metawrap2.modules import _common, binning


def _cost_lines(text: str) -> list:
    """Every line of *text* that is a cost summary, however it is worded."""
    return [line for line in text.splitlines() if "wall" in line and "cores busy" in line]


def _binning_dry_run(tmp_path, monkeypatch) -> None:
    """One standalone ``binning`` run that executes no tools. Leaves no global state set."""
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">contig_1\nACGTACGTACGT\n")
    reads = tmp_path / "reads"
    reads.mkdir()
    for mate in ("1", "2"):
        (reads / ("sample_1_%s.fq.gz" % mate)).write_bytes(b"")
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    monkeypatch.setattr(binning, "env_for", lambda module, settings: None)
    # --dry-run is a global flag cli.main strips, so the module parser never sees it; it has to
    # be set on the runner, and reset afterwards or it leaks into the next test.
    command.set_dry_run(True)
    try:
        assert (
            binning.main(
                [
                    "-a",
                    str(assembly),
                    "-o",
                    str(tmp_path / "out"),
                    "-t",
                    "4",
                    "--metabat2",
                    str(reads / "sample_1_1.fq.gz"),
                    str(reads / "sample_1_2.fq.gz"),
                ]
            )
            == 0
        )
    finally:
        command.set_dry_run(False)


def test_a_standalone_module_run_prints_one_cost_line(tmp_path, monkeypatch, capsys):
    """The run's cost is printed when it ends, exactly once.

    Regresses: a user running a module directly saw no cost line at all, so the one number that
    shows a 24-thread step using 1.0 core - the reason unused CPU budget gets noticed - was
    visible only to people who drove the module through ``run``.
    """
    _binning_dry_run(tmp_path, monkeypatch)
    lines = _cost_lines(capsys.readouterr().out)
    assert len(lines) == 1
    assert "peak memory" in lines[0]
    assert "binning" in lines[0]


def test_a_step_of_a_study_prints_no_cost_line_of_its_own(tmp_path, monkeypatch, capsys):
    """Under ``run`` the driver measures the step from outside and prints the only cost line.

    Regresses: two figures for the same step, one measured from inside the process and one from
    outside it, would never quite agree - and two cost lines that disagree are worse than the
    one that was missing.
    """
    monkeypatch.setenv(_common.DRIVEN_BY_RUN_ENV, "1")
    _binning_dry_run(tmp_path, monkeypatch)
    assert _cost_lines(capsys.readouterr().out) == []


def test_the_cost_line_comes_from_the_shared_summariser(tmp_path):
    """One rendering of the cost, shared with the driver, not a second one that can drift.

    Regresses: a hand-built "%.1f GB, %d cores" here would report the same run in different
    units from ``run``'s line for the same step.
    """
    rec = _common.RunRecorder(
        module="binning",
        version="test",
        parameters={},
        config_file=None,
        conda_env=None,
        databases={},
        inputs=[],
    )
    measured = _common.self_usage(rec)
    assert set(measured) == {"wall_seconds", "peak_rss_bytes", "user_seconds", "system_seconds"}
    # The whole process tree, as os.wait4 reports it to the driver: never less than this
    # process's own high-water mark.
    assert measured["peak_rss_bytes"] >= usage.self_peak_rss_bytes()
    assert _cost_lines(usage.summarise(measured))


def test_one_run_reports_its_cost_once_however_often_finish_run_is_called(tmp_path, capsys):
    """Two cost lines for one run would only raise the question of which number is real.

    Regresses: ``finish_run`` is not called exactly once per run on every path - pathways calls
    it before an early ``return 0`` that sits inside the ``try`` whose ``finally`` calls it
    again - so an unguarded line prints twice there, with two different figures.
    """
    rec = _common.RunRecorder(
        module="pathways",
        version="test",
        parameters={},
        config_file=None,
        conda_env=None,
        databases={},
        inputs=[],
    )
    out = tmp_path / "out"
    out.mkdir()
    command.set_dry_run(True)
    try:
        _common.finish_run(rec, str(out))
        _common.finish_run(rec, str(out))
    finally:
        command.set_dry_run(False)
    assert len(_cost_lines(capsys.readouterr().out)) == 1


def test_a_standalone_run_says_where_its_live_log_is(tmp_path, monkeypatch, capsys):
    """The run names its own log at the start, not only when something fails.

    Regresses: during a four-hour step nothing said ``metawrap2.log`` existed, so a user with a
    silent screen had to guess - or, as reported, go looking inside a tool's own working
    directory for a log the pipeline never mentioned.
    """
    _binning_dry_run(tmp_path, monkeypatch)
    out = capsys.readouterr().out
    assert str(tmp_path / "out" / "metawrap2.log") in out
    assert "run.stdout" in out
