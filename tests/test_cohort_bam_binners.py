"""Tests that `cohort` gives binning the alignments that BAM-reading binners need.

The default binners are metaBAT2 + SemiBin2. SemiBin2 reads alignments, and `binning` rejects
`--depth` alone for it, so a cohort plan with only the depth slice failed at the first binning
step, after read QC, assembly and coverage. `coverage --split-bams` now writes one BAM per assembly
per sample, and binning gets both the depth table and the BAMs.
"""

import os

from metawrap2.commands import run as driver
from metawrap2.modules import binning as _binning

SHEET = """
[settings]
coassemble = false

[[samples]]
name = "A"
r1 = "a_1.fastq"
r2 = "a_2.fastq"

[[samples]]
name = "B"
r1 = "b_1.fastq"
r2 = "b_2.fastq"
"""


def _plan(tmp_path, body=SHEET):
    path = tmp_path / "samples.toml"
    path.write_text(body)
    study = driver.load_sheet(str(path))
    layout = driver.Layout(str(tmp_path / "study"))
    return study, layout, driver.plan_steps(study, layout)


def _argv(steps, name):
    return next(s.argv for s in steps if s.name == name)


def test_the_default_cohort_gives_binning_both_the_table_and_the_alignments(tmp_path, capsys):
    _study, layout, steps = _plan(tmp_path)
    coverage = _argv(steps, "coverage")
    assert "--split-bams" in coverage

    binning = _argv(steps, "binning:A")
    assert "--depth" in binning and "--bams" in binning
    assert binning[binning.index("--bams") + 1] == layout.coverage_bams("A")
    # The BAMs are for this assembly, not the catalogue: one directory per assembly label.
    assert binning[binning.index("--bams") + 1].endswith(os.path.join("bams", "A"))


def test_a_table_only_binner_set_does_not_pay_for_the_split(tmp_path):
    """`--preset fast` (metaBAT2 + MaxBin2) reads only the table, so no BAMs are split."""
    body = SHEET.replace(
        "coassemble = false", 'coassemble = false\n\n[modules]\nbinning = "--preset fast"'
    )
    _study, _layout, steps = _plan(tmp_path, body)
    assert "--split-bams" not in _argv(steps, "coverage")
    binning = _argv(steps, "binning:A")
    assert "--depth" in binning and "--bams" not in binning


def test_every_bam_reading_binner_triggers_the_split(tmp_path):
    """Every binner in BAM_BINNERS (SemiBin2, CONCOCT, VAMB, COMEBin) triggers the split."""
    for binner in _binning.BAM_BINNERS:
        body = SHEET.replace(
            "coassemble = false",
            'coassemble = false\n\n[modules]\nbinning = "--metabat2 --%s"' % binner,
        )
        _study, _layout, steps = _plan(tmp_path, body)
        assert "--split-bams" in _argv(steps, "coverage"), binner
        assert "--bams" in _argv(steps, "binning:A"), binner


def test_binning_accepts_depth_and_bams_together():
    """`binning` accepts `--depth` with `--bams` for SemiBin2."""
    parser_opts = _binning._parse_args(
        ["-a", "asm.fa", "-o", "out", "--depth", "d.txt", "--bams", "bams", "--semibin2"]
    )
    assert parser_opts.depth and parser_opts.bams and parser_opts.semibin2
