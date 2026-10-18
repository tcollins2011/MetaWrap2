"""A default `cohort` planned a binning step that could not run, and paid for a study first.

With the cohort coverage stage on, every sample's binning got `--depth <slice>`. The default binner
pair is metaBAT2 + SemiBin2, and SemiBin2 reads the alignments themselves - `binning` refuses
`--depth` with it by name. So the shipped cohort path died at the first binning step, after read QC,
assembly and a cohort-wide alignment pass had all been paid for. A user hit it on 36 water samples.

The fix is the missing piece rather than a narrower default: `coverage --split-bams` writes one BAM
per assembly per sample, and the binning step is given both - the table for metaBAT2, the
alignments for SemiBin2.
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
    # The BAMs are this assembly's, not the catalogue's: one directory per assembly label.
    assert binning[binning.index("--bams") + 1].endswith(os.path.join("bams", "A"))


def test_a_table_only_binner_set_does_not_pay_for_the_split(tmp_path):
    """`--preset fast` is metaBAT2 + MaxBin2, both of which read the table. No BAMs needed."""
    body = SHEET.replace(
        "coassemble = false", 'coassemble = false\n\n[modules]\nbinning = "--preset fast"'
    )
    _study, _layout, steps = _plan(tmp_path, body)
    assert "--split-bams" not in _argv(steps, "coverage")
    binning = _argv(steps, "binning:A")
    assert "--depth" in binning and "--bams" not in binning


def test_every_bam_reading_binner_triggers_the_split(tmp_path):
    """Not just SemiBin2: CONCOCT, VAMB and COMEBin read alignments too."""
    for binner in _binning.BAM_BINNERS:
        body = SHEET.replace(
            "coassemble = false",
            'coassemble = false\n\n[modules]\nbinning = "--metabat2 --%s"' % binner,
        )
        _study, _layout, steps = _plan(tmp_path, body)
        assert "--split-bams" in _argv(steps, "coverage"), binner
        assert "--bams" in _argv(steps, "binning:A"), binner


def test_binning_accepts_depth_and_bams_together():
    """The refusal has to stay for `--depth` alone and go away when the alignments are there."""
    parser_opts = _binning._parse_args(
        ["-a", "asm.fa", "-o", "out", "--depth", "d.txt", "--bams", "bams", "--semibin2"]
    )
    assert parser_opts.depth and parser_opts.bams and parser_opts.semibin2
