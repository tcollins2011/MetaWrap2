"""Tests for what a cohort run tells the user it is doing: the pinned-steps warning, the ``-j``
advice, and the order `metawrap2 status` reports a plan in.

`run` and `cohort` both read "a sample sheet" and the docs present them as siblings, so the same
file gets fed to both. A sheet that pins ``steps`` for a co-assembly has no ``coverage`` and no
``dereplicate`` in it - neither does anything when there is one assembly - and a cohort run from it
produced a pile of per-sample bins and no catalogue without a word either way, which is the
opposite of why anyone chooses ``cohort``. The mirror case, ``coassemble = true`` under ``cohort``,
has been overridden out loud since it was written; these pin the other direction and keep that one
honest. The last two cover the other half of a first cohort run: being told when the thread grant
went unused, and being shown the steps in the order they ran.
"""

from __future__ import annotations

from typing import List

from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.commands import status_cmd
from metawrap2.manifest import COMPLETED, Manifest

PINNED_SHEET = """
[settings]
steps = ["read_qc", "assembly", "binning", "bin_refinement", "reassemble_bins", "quant_bins"]

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"
"""

COHORT_SHEET = """
[settings]
coassemble = false

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"
"""

#: A sheet that asks for the cohort stages but says nothing about `coassemble` - so reading it
#: defaults to co-assembling, and the mode only becomes per-sample once `cohort` has had its say.
PINNED_COVERAGE_SHEET = """
[settings]
steps = ["read_qc", "assembly", "coverage", "binning", "bin_refinement", "dereplicate"]

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"
"""


def _flat(text: str) -> str:
    """The message banners wrap and pad, so compare the words and not the layout."""
    return " ".join(text.replace("*", " ").split())


def _planned(tmp_path, monkeypatch, capsys, body, command=cohort_cmd.main):
    """Run a mode as far as the plan, executing nothing. Returns the study it settled on and
    everything it said about it."""
    seen = {}

    def fake_execute(steps, layout, study, *args, **kwargs):
        seen["study"] = study
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)
    sheet = tmp_path / "samples.toml"
    sheet.write_text(body)
    assert command([str(sheet), "-o", str(tmp_path / "out")]) == 0
    return seen["study"], _flat(capsys.readouterr().out)


def _cohort_output(tmp_path, monkeypatch, capsys, body):
    """Run `metawrap2 cohort` as far as the plan, executing nothing, and return what it said."""
    return _planned(tmp_path, monkeypatch, capsys, body)[1]


def test_pinned_steps_warn_and_name_the_cohort_stages_they_drop(tmp_path, monkeypatch, capsys):
    """What regresses: a cohort silently producing no depth table and no catalogue."""
    out = _cohort_output(tmp_path, monkeypatch, capsys, PINNED_SHEET)
    assert "the sheet pins steps and omits coverage and dereplicate" in out
    assert "one depth table across every sample or a non-redundant MAG catalogue" in out
    # Actionable, not just a complaint: the way out, and what taking it would run.
    assert "Remove `steps`" in out
    assert " -> ".join(driver.COHORT_STEPS) in out


def test_dropping_only_dereplicate_names_only_the_catalogue(tmp_path, monkeypatch, capsys):
    """What regresses: one message for two different losses, so it fits neither case."""
    body = COHORT_SHEET.replace(
        "coassemble = false",
        'coassemble = false\nsteps = ["read_qc", "assembly", "coverage", "binning"]',
    )
    out = _cohort_output(tmp_path, monkeypatch, capsys, body)
    assert "omits dereplicate" in out
    assert "a non-redundant MAG catalogue" in out
    assert "depth table" not in out


def test_a_sheet_that_pins_no_steps_does_not_warn(tmp_path, monkeypatch, capsys):
    """What regresses: a warning on the default path, which is the path that is already right."""
    out = _cohort_output(tmp_path, monkeypatch, capsys, COHORT_SHEET)
    assert "pins steps" not in out


def test_a_sheet_asking_to_coassemble_is_still_overridden_out_loud(tmp_path, monkeypatch, capsys):
    """What regresses: the mirror-image override this warning was modelled on going quiet."""
    body = COHORT_SHEET.replace("coassemble = false", "coassemble = true")
    out = _cohort_output(tmp_path, monkeypatch, capsys, body)
    assert "the sheet sets coassemble = true" in out
    assert "pins steps" not in out


def test_a_pinned_cohort_stage_survives_a_cohort_run(tmp_path, monkeypatch, capsys):
    """What regresses: `coverage` dropped while the sheet was being read, blamed on co-assembly."""
    study, out = _planned(tmp_path, monkeypatch, capsys, PINNED_COVERAGE_SHEET)
    # The user asked for both stages and is running the mode they exist in, so they run.
    assert "coverage" in study.steps and "dereplicate" in study.steps
    assert "nothing to do" not in out
    assert "pins steps" not in out


def test_the_same_sheet_under_run_still_drops_them_and_says_why(tmp_path, monkeypatch, capsys):
    """What regresses: a co-assembly paying for a depth table of one assembly, or paying silently."""
    study, out = _planned(tmp_path, monkeypatch, capsys, PINNED_COVERAGE_SHEET, command=driver.main)
    assert "coverage" not in study.steps and "dereplicate" not in study.steps
    assert "Co-assembling, so coverage, dereplicate have nothing to do" in out
    assert "Set coassemble = false" in out


def test_a_near_serial_stage_says_what_j_would_have_used_the_budget(tmp_path, monkeypatch):
    """What regresses: three steps using 1.7 cores of many, one at a time, in silence."""
    layout = driver.Layout(str(tmp_path / "study"))
    study = driver.Study(
        samples=[driver.Sample(name=s, r1="raw/%s_1.fq" % s, r2="raw/%s_2.fq" % s) for s in "ABC"]
    )
    steps = [
        driver.Step(name="read_qc:%s" % s.name, module="read_qc", argv=[], parallel_group="read_qc")
        for s in study.samples
    ]
    state = driver.StepState(layout, resume=False)

    def fake_run_step(step, _layout, threads, _config, _flags, _dry, step_state=None, memory=None):
        """A step that finishes having kept 1.7 cores busy, as read_qc measurably does."""
        started = step_state.start(step, ["metawrap2", "read_qc"], grant={"threads": threads})
        step_state.finish(
            step,
            started,
            COMPLETED,
            usage={"wall_seconds": 100.0, "user_seconds": 170.0, "system_seconds": 0.0},
        )
        return step, COMPLETED, "log"

    said: List[str] = []
    monkeypatch.setattr(driver, "run_step", fake_run_step)
    # The message as emitted, not as the console lays it out: the banner wraps a long line, and
    # what matters here is the sentence and that there is exactly one of it.
    monkeypatch.setattr(driver, "comm", said.append)
    assert driver.execute(steps, layout, study, 64, 1, None, [], state, True) == 0
    assert [m for m in said if "granted cores" in m] == [
        (
            "read_qc used 1.7 of its 64 granted cores; with 3 sample(s) to run, -j 3 would "
            "spend the same CPU in less wall time."
        )
    ]


def test_a_step_that_uses_its_grant_is_not_advised_to_run_more_at_once(tmp_path):
    """What regresses: advising -j for a step that is already saturating the machine."""
    assert driver.parallelism_advice("assembly", 70.0, 64, 3, 1) is None
    # Nor when the suggestion is what is already happening, nor when the grant could not feed it.
    assert driver.parallelism_advice("read_qc", 1.7, 64, 3, 3) is None
    assert driver.parallelism_advice("read_qc", 1.0, 2, 8, 2) is None
    assert driver.cores_busy({"wall_seconds": 100.0}) is None


def test_status_lists_the_steps_in_execution_order(tmp_path, capsys):
    """What regresses: a table reading assembly, bin_refinement, binning, quant_bins, read_qc."""
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    manifest = Manifest(status_cmd.manifest_path(str(root)))
    # Recorded in no useful order, and saved with sorted keys - which is where the alphabetical
    # table came from - including one step no plan contains.
    for name in ("quant_bins", "pangenome", "bin_refinement", "read_qc:A", "binning", "assembly"):
        manifest.finish(name, manifest.start(name, ["metawrap2", name], []), COMPLETED)

    assert status_cmd.main([str(root)]) == 0
    out = capsys.readouterr().out
    planned = ("read_qc:A", "assembly", "binning", "bin_refinement", "quant_bins")
    positions = [out.index(name) for name in planned]
    assert positions == sorted(positions)
    # A step the plan does not contain still has to appear, after the ones it does.
    assert out.index("pangenome") > positions[-1]
