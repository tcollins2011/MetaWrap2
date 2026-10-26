"""Tests for what a standalone module run prints: the live log location and the cost line.

Modules can be run directly as well as through the ``run`` driver.
"""

from metawrap2 import command, usage
from metawrap2.modules import _common, binning


def _cost_lines(text: str) -> list:
    """Return the lines of *text* that are cost summaries."""
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
    # cli.main strips the global --dry-run flag before the module parser sees it, so set it on
    # the runner and reset it afterwards so it does not leak into other tests.
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
    """A standalone run prints one cost line when it ends.

    Previously only runs through ``run`` printed it.
    """
    _binning_dry_run(tmp_path, monkeypatch)
    lines = _cost_lines(capsys.readouterr().out)
    assert len(lines) == 1
    assert "peak memory" in lines[0]
    assert "binning" in lines[0]


def test_a_step_of_a_study_prints_no_cost_line_of_its_own(tmp_path, monkeypatch, capsys):
    """Under ``run``, the module prints no cost line; the driver prints its own.

    Figures measured inside and outside the process would differ slightly.
    """
    monkeypatch.setenv(_common.DRIVEN_BY_RUN_ENV, "1")
    _binning_dry_run(tmp_path, monkeypatch)
    assert _cost_lines(capsys.readouterr().out) == []


def test_the_cost_line_comes_from_the_shared_summariser(tmp_path):
    """The cost line uses the same summariser as the driver, so units match."""
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
    # Covers the whole process tree, as os.wait4 reports to the driver, so it is at least this
    # process's own peak.
    assert measured["peak_rss_bytes"] >= usage.self_peak_rss_bytes()
    assert _cost_lines(usage.summarise(measured))


def test_one_run_reports_its_cost_once_however_often_finish_run_is_called(tmp_path, capsys):
    """The cost line prints once even when ``finish_run`` is called twice.

    pathways calls ``finish_run`` before an early ``return 0`` inside a ``try`` whose ``finally``
    calls it again.
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
    """A standalone run prints the path of its log at the start."""
    _binning_dry_run(tmp_path, monkeypatch)
    out = capsys.readouterr().out
    assert str(tmp_path / "out" / "metawrap2.log") in out
    assert "run.stdout" in out
