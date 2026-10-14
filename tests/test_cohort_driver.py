"""Tests for the cohort path of the driver: per-sample assembly, one coverage pass, one catalogue.

The co-assembly path was the only one wired correctly. With ``coassemble = false`` the planner
produced per-sample assemblies and then pointed binning and quant_bins at
``ASSEMBLY/final_assembly.fasta`` - a file that mode never writes - and reassembly at
``CLEAN_READS/ALL_READS_*``, which is only concatenated when co-assembling. Every assertion here
about a per-sample path is pinning that down, not decorating it.
"""

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


def test_a_cohort_defaults_to_the_cohort_steps(tmp_path, monkeypatch, capsys):
    study, layout, steps = _plan(tmp_path)
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

    # `metawrap2 cohort` is this plan, forced over a sheet that asked for co-assembly - not silently.
    seen = {}

    def fake_execute(steps, layout, study, *args, **kwargs):
        seen["study"] = study
        seen["names"] = _names(steps)
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)
    capsys.readouterr()
    sheet = _sheet(tmp_path, COHORT_SHEET.replace("false", "true"))
    assert cohort_cmd.main([sheet, "-o", str(tmp_path / "out")]) == 0
    assert seen["study"].coassemble is False
    assert seen["names"] == _names(steps)
    assert "coassemble = true" in capsys.readouterr().out


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
