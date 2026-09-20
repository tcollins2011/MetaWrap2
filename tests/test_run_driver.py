"""Tests for the `metawrap2 run` pipeline driver.

``plan_steps`` is deliberately pure - it touches no filesystem - so the wiring between steps
(which directory feeds which module) can be asserted directly. That wiring is exactly the
plumbing users previously did by hand out of the tutorial, so it is worth pinning down.
"""

import os

import pytest

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


# --- sample sheet parsing -----------------------------------------------------------------


def test_loads_samples_and_settings(tmp_path):
    study = driver.load_sheet(_sheet(tmp_path, PAIRED_SHEET))
    assert [s.name for s in study.samples] == ["A", "B"]
    assert study.coassemble is True
    assert study.steps[0] == "read_qc"
    assert study.extra("binning") == ["--metabat2", "--maxbin2", "--concoct"]
    assert study.extra("bin_refinement") == ["-c", "50", "-x", "10"]
    assert study.extra("kraken2") == []


def test_sample_name_is_derived_when_omitted(tmp_path):
    body = (
        '[[samples]]\nr1 = "raw/ERR011347_R1_001.fastq.gz"\nr2 = "raw/ERR011347_R2_001.fastq.gz"\n'
    )
    study = driver.load_sheet(_sheet(tmp_path, body))
    assert study.samples[0].name == "ERR011347"


def test_single_end_and_interleaved_layouts(tmp_path):
    body = (
        '[[samples]]\nname = "S"\nr1 = "raw/S.fastq"\nlayout = "single"\n'
        '[[samples]]\nname = "I"\nr1 = "raw/I.fastq"\nlayout = "interleaved"\n'
    )
    study = driver.load_sheet(_sheet(tmp_path, body))
    assert study.samples[0].layout_flags() == ["--single-end"]
    assert study.samples[1].layout_flags() == ["--interleaved"]


@pytest.mark.parametrize(
    "body,message",
    [
        ("", "defines no samples"),
        ('[[samples]]\nname = "A"\n', "has no 'r1'"),
        ('[[samples]]\nname = "A"\nr1 = "a.fq"\nlayout = "paired"\n', "no 'r2'"),
        (
            '[[samples]]\nname = "A"\nr1 = "a.fq"\nlayout = "sideways"\n',
            "paired/single/interleaved",
        ),
        (
            (
                '[[samples]]\nname = "A"\nr1 = "a_1.fq"\nr2 = "a_2.fq"\n'
                '[[samples]]\nname = "A"\nr1 = "b_1.fq"\nr2 = "b_2.fq"\n'
            ),
            "more than once",
        ),
        (
            '[settings]\nsteps = ["nonsense"]\n[[samples]]\nname = "A"\nr1 = "a.fq"\n',
            "unknown step",
        ),
    ],
)
def test_malformed_sheets_are_rejected_with_a_useful_message(tmp_path, body, message):
    with pytest.raises((ValueError, TypeError), match=message):
        driver.load_sheet(_sheet(tmp_path, body))


def test_directory_of_reads_is_paired_automatically(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    for name in ("X_R1_001.fastq", "X_R2_001.fastq", "Y_1.fastq", "Y_2.fastq", "notes.txt"):
        (raw / name).write_text("@r\nA\n+\nI\n")
    study = driver.study_from_directory(str(raw))
    assert sorted(s.name for s in study.samples) == ["X", "Y"]


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
    _study, _layout, steps = _plan(tmp_path)
    assert [s.name for s in steps] == [
        "read_qc:A",
        "read_qc:B",
        "assembly",
        "binning",
        "bin_refinement",
        "reassemble_bins",
        "quant_bins",
    ]


def test_read_qc_steps_are_marked_parallel(tmp_path):
    _study, _layout, steps = _plan(tmp_path)
    qc = [s for s in steps if s.module == "read_qc"]
    assert {s.parallel_group for s in qc} == {"read_qc"}
    # whole-study steps are not parallelisable with each other
    assert next(s for s in steps if s.name == "binning").parallel_group is None


def test_assembly_consumes_the_concatenated_clean_reads(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    argv = _argv_for(steps, "assembly")
    assert os.path.join(layout.clean_reads, "ALL_READS_1.fastq") in argv
    assert os.path.join(layout.clean_reads, "ALL_READS_2.fastq") in argv


def test_binning_consumes_the_assembly_and_every_sample(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    argv = _argv_for(steps, "binning")
    assert layout.assembly in argv
    for sample in ("A", "B"):
        for path in layout.clean_pair(sample):
            assert path in argv


def test_refinement_is_handed_all_three_binner_outputs(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    argv = _argv_for(steps, "bin_refinement")
    for binner in ("metabat2", "maxbin2", "concoct"):
        assert layout.binner_bins(binner) in argv


def test_thresholds_from_the_sheet_name_the_refined_bin_directory(tmp_path):
    """bin_refinement names its output metawrap_<c>_<x>_bins, so -c/-x must be honoured."""
    _study, layout, steps = _plan(tmp_path)
    assert layout.refined_bins(50, 10) in _argv_for(steps, "reassemble_bins")


def test_default_thresholds_are_used_when_the_sheet_is_silent(tmp_path):
    body = PAIRED_SHEET.replace('bin_refinement = "-c 50 -x 10"\n', "")
    _study, layout, steps = _plan(tmp_path, body)
    assert layout.refined_bins(70, 10) in _argv_for(steps, "reassemble_bins")


def test_downstream_steps_use_the_reassembled_bins_when_reassembly_runs(tmp_path):
    body = PAIRED_SHEET.replace(
        'steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]',
        'steps = ["binning", "bin_refinement", "reassemble_bins", "quant_bins", '
        '"classify_bins", "annotate_bins", "blobology"]',
    )
    _study, layout, steps = _plan(tmp_path, body)
    for step in ("quant_bins", "classify_bins", "annotate_bins", "blobology"):
        assert layout.reassembled_bins in _argv_for(steps, step), step


def test_downstream_steps_fall_back_to_refined_bins_without_reassembly(tmp_path):
    body = PAIRED_SHEET.replace(
        'steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]',
        'steps = ["binning", "bin_refinement", "quant_bins"]',
    )
    _study, layout, steps = _plan(tmp_path, body)
    assert layout.refined_bins(50, 10) in _argv_for(steps, "quant_bins")


def test_without_read_qc_the_raw_inputs_are_used(tmp_path):
    body = PAIRED_SHEET.replace(
        'steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]',
        'steps = ["binning"]',
    )
    _study, _layout, steps = _plan(tmp_path, body)
    argv = _argv_for(steps, "binning")
    assert "raw/A_1.fastq.gz" in argv and "raw/B_2.fastq.gz" in argv


def test_per_sample_assembly_when_coassembly_is_off(tmp_path):
    _study, layout, steps = _plan(tmp_path, coassemble=False)
    names = [s.name for s in steps if s.module == "assembly"]
    assert names == ["assembly:A", "assembly:B"]
    argv = _argv_for(steps, "assembly:A")
    assert layout.sample_assembly_dir("A") in argv


def test_single_end_study_plans_no_second_read_file(tmp_path):
    body = (
        '[settings]\nsteps = ["read_qc", "assembly"]\n'
        '[[samples]]\nname = "S"\nr1 = "raw/S.fastq"\nlayout = "single"\n'
    )
    _study, _layout, steps = _plan(tmp_path, body)
    qc = _argv_for(steps, "read_qc:S")
    assert "--single-end" in qc and "-2" not in qc
    assert "-2" not in _argv_for(steps, "assembly")


# --- commands and state -------------------------------------------------------------------


def test_step_command_is_a_plain_module_invocation(tmp_path):
    _study, _layout, steps = _plan(tmp_path)
    cmd = steps[0].command(threads=8, config="/cfg.toml", global_flags=["--force"])
    assert cmd[1:3] == ["-m", "metawrap2.cli"]
    assert "--force" in cmd and "read_qc" in cmd
    assert cmd[cmd.index("-t") + 1] == "8"
    assert cmd[cmd.index("--config") + 1] == "/cfg.toml"


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
    step, _src, _out = _finished_step(tmp_path, state)
    assert state.skip_reason(step) is not None


def test_resume_reruns_a_step_whose_output_is_gone(tmp_path):
    """A marker file said "done"; the manifest checks the outputs are actually there."""
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, _src, out = _finished_step(tmp_path, state)
    out.unlink()
    assert state.skip_reason(step) is None
    assert "is gone" in state.rerun_reason(step)


def test_resume_reruns_a_step_whose_output_is_empty(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, _src, out = _finished_step(tmp_path, state)
    out.write_text("")
    assert state.skip_reason(step) is None
    assert "now empty" in state.rerun_reason(step)


def test_resume_reruns_a_step_whose_input_changed(tmp_path):
    """The stale-input case: skipping here would mix results from two different inputs."""
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, src, _out = _finished_step(tmp_path, state)
    src.write_text("@r\nTTTT\n+\nIIII\n")  # same length, different content
    assert state.skip_reason(step) is None
    assert "has changed" in state.rerun_reason(step)


def test_resume_reruns_an_interrupted_step(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    src = tmp_path / "in.fastq"
    src.write_text("@r\nACGT\n+\nIIII\n")
    step = driver.Step(name="binning", module="binning", argv=[], inputs=[str(src)], outputs=[])
    started = state.start(step, ["echo"])
    state.finish(step, started, "interrupted")
    assert state.skip_reason(step) is None
    assert "interrupted" in state.rerun_reason(step)


def test_state_is_ignored_without_resume(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, _src, _out = _finished_step(tmp_path, state, name="binning")
    assert driver.StepState(layout, resume=False).skip_reason(step) is None


def test_manifest_records_command_timing_and_files(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=True)
    step, src, out = _finished_step(tmp_path, state)
    record = state.manifest.step(step.name)
    assert record.status == "completed"
    assert record.data["command"] == ["echo", "hi"]
    assert "duration_seconds" in record.data and "started" in record.data
    assert [e["path"] for e in record.inputs] == [str(src)]
    assert [e["path"] for e in record.outputs] == [str(out)]
    # and it is on disk, at the documented location
    assert os.path.isfile(os.path.join(layout.state_dir, "manifest.json"))


def test_concatenation_is_skipped_when_the_output_is_current(tmp_path, capsys):
    study, layout, _steps = _plan(tmp_path)
    os.makedirs(layout.clean_reads)
    sources = []
    for sample in ("A", "B"):
        for path in layout.clean_pair(sample):
            with open(path, "w") as fh:
                fh.write("@r\nACGT\n+\nIIII\n")
            sources.append(path)
    driver.concatenate_for_coassembly(layout, study, dry_run=False)
    combined = os.path.join(layout.clean_reads, "ALL_READS_1.fastq")
    assert os.path.isfile(combined)
    first = os.path.getmtime(combined)

    driver.concatenate_for_coassembly(layout, study, dry_run=False)
    assert "already up to date" in capsys.readouterr().out
    assert os.path.getmtime(combined) == first


def test_example_sheet_is_valid(tmp_path):
    path = tmp_path / "example.toml"
    path.write_text(driver.EXAMPLE_SHEET)
    study = driver.load_sheet(str(path))
    assert len(study.samples) == 2
    assert study.extra("binning") == ["--metabat2", "--maxbin2", "--concoct"]


# --- refinement must only be given the bin sets that binning produced --------------------


def test_refinement_inputs_follow_the_selected_binners(tmp_path):
    """A study selecting two binners used to fail: -B pointed at a maxbin2 dir never created."""
    body = PAIRED_SHEET.replace(
        'binning = "--metabat2 --maxbin2 --concoct"', 'binning = "--metabat2 --concoct"'
    )
    _study, layout, steps = _plan(tmp_path, body)
    argv = _argv_for(steps, "bin_refinement")
    assert argv[argv.index("-A") + 1] == layout.binner_bins("metabat2")
    assert argv[argv.index("-B") + 1] == layout.binner_bins("concoct")
    assert "-C" not in argv
    assert "maxbin2" not in " ".join(argv)


def test_single_binner_study_passes_only_A(tmp_path):
    body = PAIRED_SHEET.replace(
        'binning = "--metabat2 --maxbin2 --concoct"', 'binning = "--concoct"'
    )
    _study, layout, steps = _plan(tmp_path, body)
    argv = _argv_for(steps, "bin_refinement")
    assert argv[argv.index("-A") + 1] == layout.binner_bins("concoct")
    assert "-B" not in argv and "-C" not in argv


def test_selected_binners_helper():
    assert driver.selected_binners(["--metabat2", "--concoct"]) == ["metabat2", "concoct"]
    assert driver.selected_binners(["--maxbin2"]) == ["maxbin2"]
    # nothing recognisable -> the defaults, in -A/-B/-C order
    assert driver.selected_binners(["-t", "8"]) == ["metabat2", "maxbin2", "concoct"]


def test_without_refinement_downstream_uses_the_first_selected_binner(tmp_path):
    body = PAIRED_SHEET.replace(
        'steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]',
        'steps = ["binning", "quant_bins"]',
    ).replace('binning = "--metabat2 --maxbin2 --concoct"', 'binning = "--concoct"')
    _study, layout, steps = _plan(tmp_path, body)
    assert layout.binner_bins("concoct") in _argv_for(steps, "quant_bins")
