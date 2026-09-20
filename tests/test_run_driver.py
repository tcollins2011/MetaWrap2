"""Tests for the `metawrap2 run` pipeline driver.

``plan_steps`` is deliberately pure - it touches no filesystem - so the wiring between steps
(which directory feeds which module) can be asserted directly. That wiring is exactly the
plumbing users previously did by hand out of the tutorial, so it is worth pinning down.
"""

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


def test_binning_consumes_the_assembly_and_every_sample(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    argv = _argv_for(steps, "binning")
    assert layout.assembly in argv
    for sample in ("A", "B"):
        for path in layout.clean_pair(sample):
            assert path in argv


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
