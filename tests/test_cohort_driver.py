"""Tests for the cohort path of the driver: per-sample assembly, one coverage pass, one catalogue.

The co-assembly path was the only one wired correctly. With ``coassemble = false`` the planner
produced per-sample assemblies and then pointed binning and quant_bins at
``ASSEMBLY/final_assembly.fasta`` - a file that mode never writes - and reassembly at
``CLEAN_READS/ALL_READS_*``, which is only concatenated when co-assembling. Every assertion here
about a per-sample path is pinning that down, not decorating it.
"""

import pytest

from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver

COHORT_SHEET = """
[settings]
coassemble = false

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq.gz"
r2 = "raw/B_2.fastq.gz"
"""


def _sheet(tmp_path, body=COHORT_SHEET):
    path = tmp_path / "samples.toml"
    path.write_text(body)
    return str(path)


def _plan(tmp_path, body=COHORT_SHEET):
    study = driver.load_sheet(_sheet(tmp_path, body))
    layout = driver.Layout(str(tmp_path / "study"))
    return study, layout, driver.plan_steps(study, layout)


def _argv_for(steps, name):
    return next(s.argv for s in steps if s.name == name)


def _names(steps):
    return [s.name for s in steps]


# --- which steps a cohort runs ------------------------------------------------------------


def test_a_cohort_defaults_to_the_cohort_steps(tmp_path):
    study, _layout, steps = _plan(tmp_path)
    assert study.steps == driver.COHORT_STEPS
    # The two stages that only exist in this mode are on by default; the expensive per-sample one
    # is not.
    assert "coverage" in study.steps and "dereplicate" in study.steps
    assert "reassemble_bins" not in study.steps
    assert _names(steps) == [
        "read_qc:A",
        "read_qc:B",
        "assembly:A",
        "assembly:B",
        "coverage",
        "binning:A",
        "binning:B",
        "bin_refinement:A",
        "bin_refinement:B",
        "dereplicate",
        "quant_bins",
    ]


def test_a_coassembly_still_defaults_to_the_coassembly_steps(tmp_path):
    study = driver.load_sheet(_sheet(tmp_path, COHORT_SHEET.replace("false", "true")))
    assert study.steps == driver.DEFAULT_STEPS
    assert "coverage" not in study.steps and "dereplicate" not in study.steps


def test_cohort_only_steps_are_dropped_from_a_coassembly_with_a_reason(tmp_path, capsys):
    body = COHORT_SHEET.replace(
        "coassemble = false",
        'coassemble = true\nsteps = ["assembly", "coverage", "binning", "dereplicate"]',
    )
    study = driver.load_sheet(_sheet(tmp_path, body))
    assert study.steps == ("assembly", "binning")
    out = capsys.readouterr().out
    assert "coverage, dereplicate" in out and "one assembly, one bin set" in out


def test_named_steps_are_not_replaced_by_the_mode_default(tmp_path):
    body = COHORT_SHEET.replace("coassemble = false", 'coassemble = false\nsteps = ["binning"]')
    study = driver.load_sheet(_sheet(tmp_path, body))
    assert study.steps == ("binning",)


# --- the wiring ---------------------------------------------------------------------------


def test_every_step_reads_a_path_the_cohort_actually_produces(tmp_path):
    study = driver.load_sheet(_sheet(tmp_path))
    study.steps = tuple(driver.ALL_STEPS)
    study.steps_chosen = True
    layout = driver.Layout(str(tmp_path / "study"))
    steps = driver.plan_steps(study, layout)
    # Nothing may reference the co-assembly's outputs: in this mode they are never written.
    forbidden = (
        layout.assembly,
        layout.binner_bins("metabat2"),
        layout.binner_bins("semibin2"),
        layout.refinement_dir + "/metawrap",
        layout.reassembled_bins,
        layout.clean_reads + "/ALL_READS",
    )
    for step in steps:
        for token in step.argv + list(step.inputs):
            assert not any(token.startswith(bad) for bad in forbidden), (step.name, token)


def test_coverage_takes_every_assembly_labelled_by_sample(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    argv = _argv_for(steps, "coverage")
    for sample in ("A", "B"):
        assert "%s=%s" % (sample, layout.sample_assembly(sample)) in argv
    # Every sample's reads, once: the whole point is N alignments against one catalogue.
    for sample in ("A", "B"):
        for path in layout.clean_pair(sample):
            assert argv.count(path) == 1


def test_binning_takes_its_own_assembly_and_its_slice_of_the_depth_table(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    for sample in ("A", "B"):
        argv = _argv_for(steps, "binning:%s" % sample)
        assert layout.sample_assembly(sample) in argv
        assert argv[argv.index("--depth") + 1] == layout.coverage_depth(sample)
        # No reads: the coverage stage already aligned them, which is the saving.
        assert not any(a.endswith(".fastq") for a in argv)
        assert layout.sample_binning_dir(sample) in argv


def test_binning_falls_back_to_its_own_reads_without_the_coverage_step(tmp_path):
    body = COHORT_SHEET.replace(
        "coassemble = false",
        'coassemble = false\nsteps = ["assembly", "binning"]',
    )
    _study, layout, steps = _plan(tmp_path, body)
    argv = _argv_for(steps, "binning:A")
    assert "--depth" not in argv
    assert layout.sample_assembly("A") in argv
    # read_qc is not in the steps, so the raw reads are what binning gets.
    assert "raw/A_1.fastq.gz" in argv


def test_refinement_consolidates_each_sample_on_its_own(tmp_path):
    _study, layout, steps = _plan(tmp_path)
    for sample in ("A", "B"):
        argv = _argv_for(steps, "bin_refinement:%s" % sample)
        assert layout.sample_refinement_dir(sample) in argv
        for binner in driver.selected_binners(list(driver.DEFAULT_BINNERS)):
            assert layout.sample_binner_bins(sample, binner) in argv
        # Never another sample's bins: consolidation across samples is dereplication's job.
        other = "B" if sample == "A" else "A"
        assert not any(("/%s/" % other) in a for a in argv)


def test_dereplicate_collects_every_sample_and_quantifies_the_catalogue(tmp_path):
    study, layout, steps = _plan(tmp_path)
    comp, cont = driver._quality_thresholds(study)
    argv = _argv_for(steps, "dereplicate")
    for sample in ("A", "B"):
        assert layout.sample_refined_bins(sample, comp, cont) in argv
    assert layout.catalogue_dir in argv

    quant = _argv_for(steps, "quant_bins")
    assert quant[quant.index("-b") + 1] == layout.catalogue_genomes
    # No -a: near-identical MAGs from different samples would divide one organism's reads between
    # them, and there is no single assembly to quantify against anyway.
    assert "-a" not in quant


def test_reassembly_runs_per_sample_against_that_sample_reads(tmp_path):
    body = COHORT_SHEET.replace(
        "coassemble = false",
        'coassemble = false\nsteps = ["assembly", "coverage", "binning", "bin_refinement",'
        ' "reassemble_bins", "dereplicate"]',
    )
    study, layout, steps = _plan(tmp_path, body)
    comp, cont = driver._quality_thresholds(study)
    argv = _argv_for(steps, "reassemble_bins:A")
    assert argv[argv.index("-b") + 1] == layout.sample_refined_bins("A", comp, cont)
    assert argv[argv.index("-1") + 1] == "raw/A_1.fastq.gz"
    # Dereplication then sees the reassembled bins, not the pre-reassembly ones.
    derep = _argv_for(steps, "dereplicate")
    assert layout.sample_reassembled_bins("A") in derep


def test_per_sample_stages_are_grouped_so_a_cohort_runs_a_stage_at_a_time(tmp_path):
    _study, _layout, steps = _plan(tmp_path)
    groups = {s.name: s.parallel_group for s in steps}
    assert groups["binning:A"] == groups["binning:B"] == "binning"
    assert groups["bin_refinement:A"] == "bin_refinement"
    # The two cohort-wide steps are barriers, not members of a group.
    assert groups["coverage"] is None and groups["dereplicate"] is None


# --- the alias ----------------------------------------------------------------------------


def test_cohort_forces_per_sample_assembly_over_the_sheet(tmp_path, monkeypatch, capsys):
    seen = {}

    def fake_execute(steps, layout, study, *args, **kwargs):
        seen["study"] = study
        seen["names"] = _names(steps)
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)
    sheet = _sheet(tmp_path, COHORT_SHEET.replace("false", "true"))
    assert cohort_cmd.main([sheet, "-o", str(tmp_path / "out")]) == 0
    assert seen["study"].coassemble is False
    assert "coverage" in seen["names"] and "dereplicate" in seen["names"]
    # A sheet that asked for the other mode is overridden, but not silently.
    assert "coassemble = true" in capsys.readouterr().out


def test_cohort_is_the_same_plan_as_run_no_coassembly(tmp_path, monkeypatch):
    plans = []

    def fake_execute(steps, layout, study, *args, **kwargs):
        plans.append([(s.name, tuple(s.argv)) for s in steps])
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)
    sheet = _sheet(tmp_path)
    out = str(tmp_path / "out")
    assert cohort_cmd.main([sheet, "-o", out]) == 0
    assert driver.main([sheet, "-o", out, "--no-coassembly"]) == 0
    assert plans[0] == plans[1]


def test_cohort_is_reachable_from_the_cli(tmp_path):
    from metawrap2 import cli

    assert "cohort" in cli.UTILITIES
    assert "cohort" in cli.help_message()


@pytest.mark.parametrize("skipped", ["coverage", "dereplicate"])
def test_skipping_a_cohort_step_leaves_a_runnable_plan(tmp_path, skipped):
    study = driver.load_sheet(_sheet(tmp_path))
    study.steps = tuple(s for s in study.steps if s != skipped)
    layout = driver.Layout(str(tmp_path / "study"))
    steps = driver.plan_steps(study, layout)
    assert skipped not in _names(steps)
    # Without dereplication there is no catalogue, so quant_bins must not ask for one.
    quant = _argv_for(steps, "quant_bins")
    if skipped == "dereplicate":
        assert layout.catalogue_genomes not in quant
