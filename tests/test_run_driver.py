"""Tests for the `metawrap2 run` pipeline driver.

``plan_steps`` does not touch the filesystem, so the wiring between steps (which directory feeds
which module, previously done by hand from the tutorial) can be asserted directly.
"""

from metawrap2.commands import run as driver


def _sheet(tmp_path, body):
    path = tmp_path / "samples.toml"
    path.write_text(body)
    return str(path)


PAIRED_SHEET = """
[settings]
coassemble = true
steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]

[modules]
binning = "--metabat2 --maxbin2 --concoct"
bin_refinement = "-c 50 -x 10"

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq.gz"
r2 = "raw/B_2.fastq.gz"
"""


# --- planning -----------------------------------------------------------------------------


def _plan(tmp_path, sheet_body=PAIRED_SHEET, **overrides):
    study = driver.load_sheet(_sheet(tmp_path, sheet_body))
    for key, value in overrides.items():
        setattr(study, key, value)
    layout = driver.Layout(str(tmp_path / "study"))
    return study, layout, driver.plan_steps(study, layout)


def _argv_for(steps, name):
    return next(s.argv for s in steps if s.name == name)


def test_plan_covers_every_requested_step(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    assert [s.name for s in steps] == [
        "read_qc:A",
        "read_qc:B",
        "assembly",
        "binning",
        "bin_refinement",
        "reassemble_bins",
        "quant_bins",
    ]
    argv = _argv_for(steps, "binning")
    assert layout.assembly in argv
    for sample in ("A", "B"):
        for path in layout.clean_pair(sample):
            assert path in argv


# --- resume -----------------------------------------------------------------------------


def _finished_step(tmp_path, state, name="read_qc:A"):
    """Record *name* as a completed step with one real input and one real output."""
    src = tmp_path / "in.fastq"
    src.write_text("@r\nACGT\n+\nIIII\n")
    out = tmp_path / "out.fastq"
    out.write_text("@r\nACGT\n+\nIIII\n")
    step = driver.Step(name=name, module="read_qc", argv=[], inputs=[str(src)], outputs=[str(out)])
    started = state.start(step, ["echo", "hi"])
    state.finish(step, started, "completed")
    return step, src, out


def test_resume_skips_a_step_whose_outputs_are_intact(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, src, out = _finished_step(tmp_path, state)
    assert state.skip_reason(step) is not None
    # Only a verified step is skipped; a changed input or missing output reruns it.
    src.write_text("@r\nACGTACGT\n+\nIIIIIIII\n")
    assert state.skip_reason(step) is None
    step, src, out = _finished_step(tmp_path, state, name="read_qc:B")
    out.unlink()
    assert state.skip_reason(step) is None


# --- the driver divides the resource budget between steps ---------------------------------


def _run_main(tmp_path, monkeypatch, extra_argv=(), limits_toml=""):
    """Run the driver far enough to see what it resolved, without executing any step."""
    config = tmp_path / "metawrap2.toml"
    config.write_text(limits_toml)
    sheet = _sheet(tmp_path, PAIRED_SHEET)
    seen = {}

    def fake_execute(
        steps, layout, study, threads, jobs, config_, flags, state, dry, memory_gb=None
    ):
        seen.update(threads=threads, jobs=jobs, memory_gb=memory_gb)
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)
    rc = driver.main([sheet, "-o", str(tmp_path / "out"), "--config", str(config), *extra_argv])
    assert rc == 0
    return seen


def test_the_driver_defaults_to_the_limits_budget(tmp_path, monkeypatch, capsys):
    # Each step gets a share of the budget. A module left to decide would take the whole
    # machine, and N concurrent steps would overcommit it N times.
    from metawrap2 import limits

    monkeypatch.setattr(limits, "machine_cores", lambda: 64)
    monkeypatch.setattr(limits, "machine_memory_bytes", lambda: 256 * 1024**3)
    monkeypatch.setattr(limits, "machine_gpus", lambda: 0)
    seen = _run_main(tmp_path, monkeypatch)
    assert seen["threads"] == 32  # 50% of 64
    assert seen["memory_gb"] == 204  # 80% of 256
    # An oversized request is capped to the budget, and the message names the setting.
    capsys.readouterr()
    seen = _run_main(tmp_path, monkeypatch, extra_argv=["-t", "512", "-m", "9999"])
    assert seen["threads"] == 32
    assert seen["memory_gb"] == 204
    out = capsys.readouterr().out
    assert "max_thread_fraction" in out and "max_memory_fraction" in out


# --- refinement inputs follow the binning outputs ----------------------------------------


def _refinement_bin_sets(argv):
    return [argv[i + 1] for i, token in enumerate(argv) if token == "--bins"]


def test_refinement_inputs_follow_the_selected_binners(tmp_path):
    """Refinement gets the bin sets of the selected binners and no others.

    selected_binners knew only the classical three, so ``--metabat2 --semibin2`` refined metaBAT2
    alone, and ``--semibin2`` alone fell back to the default pair and asked for two directories
    that were never created.

    Tested with all six, because zipping against -A/-B/-C dropped any beyond three.
    """
    for selection in (
        ("metabat2", "concoct"),
        ("metabat2", "semibin2"),
        ("metabat2", "maxbin2", "concoct", "semibin2", "vamb", "comebin"),
    ):
        body = PAIRED_SHEET.replace(
            'binning = "--metabat2 --maxbin2 --concoct"',
            "binning = %s" % '"%s"' % " ".join("--" + b for b in selection),
        )
        _study, layout, steps = _plan(tmp_path, body)
        argv = _argv_for(steps, "bin_refinement")
        assert _refinement_bin_sets(argv) == [layout.binner_bins(b) for b in selection]
        if "maxbin2" not in selection:
            assert "maxbin2" not in " ".join(argv)


# --- the default binner set --------------------------------------------------------------


def test_the_default_binners_are_the_pair_the_benchmark_chose():
    """The default is metaBAT2+SemiBin2, chosen on cami2_gi 2026-10-09 against metaWRAP 1.

    53 clean high-quality MAGs against metaWRAP 1's 48, for 43.8 CPU-hours against 52.5, the first
    configuration measured to beat the baseline on both. The previous default, metaBAT2+MaxBin2,
    gave 41 for 13.8 CPU-hours and was chosen on quality per cost.

    Asserted because adding a third binner back is a one-line edit that would roughly quadruple
    the cost of every `metawrap2 run`.

    The code and `benchmarks/configurations.toml` must name the same default, or the results table
    ranks challengers against a configuration that does not run (this happened with
    `refinement_completeness`).
    """
    from metawrap2.benchmark import load_configurations
    from metawrap2.modules import binning

    assert driver.DEFAULT_BINNERS == ("--metabat2", "--semibin2")
    assert driver.DEFAULT_BINNERS == tuple(
        "--" + b for b in binning.PRESETS[binning.DEFAULT_PRESET]
    )
    variants = load_configurations()["binner_set"].variants
    shipped = [v for v in variants if v.default]
    assert len(shipped) == 1, "exactly one binner_set variant may be the default"
    assert tuple(shipped[0].args) == driver.DEFAULT_BINNERS
