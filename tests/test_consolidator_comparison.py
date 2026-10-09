"""A bin_refinement comparison has to hand every arm the assembly, not just the bins.

MetaWrap2's own consolidation works from the bin sets alone. DAS Tool reads the assembly contigs
itself and exits immediately without them, so the consolidator comparison ran its metawrap arm and
failed its DAS Tool arm in one second - which looks like DAS Tool is broken rather than like the
harness never passed `-a`.
"""

from __future__ import annotations

import pytest

from metawrap2.benchmark import runner as _runner


class _Comparison:
    def __init__(self, module="bin_refinement"):
        self.module = module
        self.name = "consolidator"
        self.common_args = ()


class _Variant:
    def __init__(self, args):
        self.args = args
        self.name = "dastool"


def _bin_sets(tmp_path):
    root = tmp_path / "bins_from"
    for name in ("metabat2_bins", "maxbin2_bins"):
        d = root / name
        d.mkdir(parents=True)
        (d / "bin_001.fasta").write_text(">c1\nACGT\n")
    return str(root)


class _Inputs:
    def __init__(self, bins_from, assembly):
        self.pairs = [("s", "r1.fastq", "r2.fastq")]
        self.assembly = assembly
        self.bins_from = bins_from

    @property
    def read_files(self):
        return ["r1.fastq", "r2.fastq"]


def test_the_assembly_is_passed_to_a_bin_refinement_arm(tmp_path):
    argv = _runner.build_argv(
        _Comparison(),
        _Variant(["--consolidator", "dastool"]),
        _Inputs(_bin_sets(tmp_path), "/data/gsa.fa"),
        out_dir=str(tmp_path / "out"),
        threads=8,
        memory_gb=100,
    )
    assert "-a" in argv and "/data/gsa.fa" in argv
    # And both bin sets, under the flag bin_refinement accepts repeatedly.
    assert argv.count("--bins") == 2


def test_no_assembly_still_builds_a_command_for_the_arms_that_do_not_need_one(tmp_path):
    """A dataset without a gold assembly should still answer the half of the question it can."""
    argv = _runner.build_argv(
        _Comparison(),
        _Variant(["--consolidator", "metawrap"]),
        _Inputs(_bin_sets(tmp_path), None),
        out_dir=str(tmp_path / "out"),
        threads=8,
        memory_gb=100,
    )
    assert "-a" not in argv
    assert "--bins" in argv


def test_resolve_inputs_finds_the_gold_assembly_for_bin_refinement(tmp_path, monkeypatch):
    """The command builder was already right; nothing was filling the field in."""
    gold = tmp_path / "gsa.fa"
    gold.write_text(">c1\nACGT\n")
    monkeypatch.setattr(_runner, "prepare_reads", lambda *a, **k: [("s", "1.fq", "2.fq")])
    monkeypatch.setattr(_runner, "gold_path", lambda *a, **k: str(gold))

    class _Dataset:
        name = "d"

    inputs = _runner.gather_inputs(
        _Dataset(),
        _Comparison(),
        data_dir=str(tmp_path),
        work_dir=str(tmp_path),
        bins_from=_bin_sets(tmp_path),
        say=lambda _m: None,
    )
    assert inputs.assembly == str(gold)


def test_a_bin_refinement_comparison_without_bins_is_still_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(_runner, "prepare_reads", lambda *a, **k: [("s", "1.fq", "2.fq")])
    monkeypatch.setattr(_runner, "gold_path", lambda *a, **k: None)

    class _Dataset:
        name = "d"

    with pytest.raises(_runner.RunnerError, match="needs the bins"):
        _runner.gather_inputs(
            _Dataset(),
            _Comparison(),
            data_dir=str(tmp_path),
            work_dir=str(tmp_path),
            bins_from=None,
            say=lambda _m: None,
        )
