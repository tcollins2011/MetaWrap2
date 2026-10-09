def test_a_bin_refinement_comparison_is_handed_the_assembly():
    """DAS Tool reads the assembly contigs itself, so without this its arm cannot start at all.

    MetaWrap2's own consolidation works from the bins alone, which is why this was missed: the default
    variant ran fine and the comparison looked merely half-broken rather than wrongly wired.
    """
    import os
    import tempfile

    from metawrap2.benchmark import runner as _runner
    from metawrap2.benchmark.registry import Comparison, Variant

    with tempfile.TemporaryDirectory() as root:
        for name in ("metabat2_bins", "maxbin2_bins"):
            os.makedirs(os.path.join(root, name))
        assembly = os.path.join(root, "asm.fa")
        with open(assembly, "w") as fh:
            fh.write(">c1\nACGT\n")
        comparison = Comparison(name="consolidator", module="bin_refinement")
        variant = Variant(name="dastool", args=("--consolidator", "dastool"))
        inputs = _runner.Inputs(pairs=[], assembly=assembly, bins_from=root)
        argv = _runner.build_argv(comparison, variant, inputs, os.path.join(root, "out"), 4, 8)
        assert "-a" in argv and argv[argv.index("-a") + 1] == assembly
        # Both input sets still reach it, under the flag bin_refinement accepts repeatedly.
        assert argv.count("--bins") == 2


# --- the depth quantifier arm ----------------------------------------------------------------


def test_depth_quantifier_arm_is_handed_the_table_binning_wrote(tmp_path, monkeypatch):
    """The arm exists to reuse binning's coverage; deriving a fresh table would measure the aligner."""
    from metawrap2.benchmark import runner as _runner
    from metawrap2.modules.binning import DEPTH_TABLE

    root = tmp_path / "binning_out"
    (root / "work_files").mkdir(parents=True)
    table = root / "work_files" / DEPTH_TABLE
    table.write_text("contigName\tcontigLen\ttotalAvgDepth\n")
    bins = root / "metawrap_70_10_bins"
    bins.mkdir()
    (bins / "bin_001.fasta").write_text(">c1\nACGT\n")

    assert _runner.depth_table_for(str(root)) == str(table)
    assert _runner.depth_table_for(str(bins)) == str(table), "also found from inside the run"


def test_a_missing_depth_table_says_which_run_should_have_written_it(tmp_path):
    import pytest

    from metawrap2.benchmark import runner as _runner

    bins = tmp_path / "just_bins"
    bins.mkdir()
    with pytest.raises(_runner.RunnerError) as excinfo:
        _runner.depth_table_for(str(bins))
    assert "metawrap2 binning" in str(excinfo.value)


def test_every_bin_set_reaches_refinement_and_a_leftover_answer_does_not(tmp_path):
    """The sets were zipped onto ("-A", "-B", "-C"), and zip truncates.

    A variant that binned with six binners was therefore refined from the three whose directory
    names sort first - metabat2_bins among those dropped - so the arm's numbers came from a
    different input set than the one it was named for, with nothing said. A consolidated set left
    behind by an earlier refinement into the same directory is also not an input: refining it
    together with the sets it came from is not a comparison anyone declared.
    """
    import os

    from metawrap2.benchmark import runner as _runner
    from metawrap2.benchmark.registry import Comparison, Variant

    binners = ("comebin", "concoct", "maxbin2", "metabat2", "semibin2", "vamb")
    for name in binners:
        (tmp_path / ("%s_bins" % name)).mkdir()
    (tmp_path / "metawrap_70_10_bins").mkdir()  # a previous refinement's own answer
    comparison = Comparison(name="consolidator", module="bin_refinement")
    variant = Variant(name="metawrap", args=())
    inputs = _runner.Inputs(pairs=[], assembly="", bins_from=str(tmp_path))
    argv = _runner.build_argv(comparison, variant, inputs, str(tmp_path / "out"), 4, 8)
    passed = {os.path.basename(argv[i + 1]) for i, tok in enumerate(argv) if tok == "--bins"}
    assert passed == {"%s_bins" % name for name in binners}
    assert "metawrap_70_10_bins" not in passed


def test_more_bin_sets_than_refinement_takes_is_an_error_not_a_silent_subset(tmp_path):
    """Seven sets is 127 CheckM runs and past bin_refinement's -A..-F; refusing beats guessing."""
    import pytest

    from metawrap2.benchmark import runner as _runner
    from metawrap2.benchmark.registry import Comparison, Variant

    for i in range(_runner.MAX_REFINEMENT_SETS + 1):
        (tmp_path / ("binner%d_bins" % i)).mkdir()
    comparison = Comparison(name="consolidator", module="bin_refinement")
    inputs = _runner.Inputs(pairs=[], assembly="", bins_from=str(tmp_path))
    with pytest.raises(_runner.RunnerError) as caught:
        _runner.build_argv(
            comparison, Variant(name="metawrap", args=()), inputs, str(tmp_path / "out"), 4, 8
        )
    assert "consolidates at most %d" % _runner.MAX_REFINEMENT_SETS in str(caught.value)
