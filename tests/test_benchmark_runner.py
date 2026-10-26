"""The runner's command for a bin_refinement comparison: every input set, and the assembly."""

import os

from metawrap2.benchmark import runner as _runner
from metawrap2.benchmark.registry import Comparison, Variant


def test_every_bin_set_reaches_refinement_and_a_leftover_answer_does_not(tmp_path):
    """All bin sets are passed with --bins, a previous consolidated set is not, and -a is passed.

    The sets used to be zipped onto ("-A", "-B", "-C"), so with six binners only the three that
    sorted first were refined. DAS Tool needs the assembly (``-a``) and exits immediately without
    it.
    """
    binners = ("comebin", "concoct", "maxbin2", "metabat2", "semibin2", "vamb")
    for name in binners:
        (tmp_path / ("%s_bins" % name)).mkdir()
    (tmp_path / "metawrap_70_10_bins").mkdir()  # a previous refinement's output
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n")
    comparison = Comparison(name="consolidator", module="bin_refinement")
    variant = Variant(name="dastool", args=("--consolidator", "dastool"))
    inputs = _runner.Inputs(pairs=[], assembly=str(assembly), bins_from=str(tmp_path))
    argv = _runner.build_argv(comparison, variant, inputs, str(tmp_path / "out"), 4, 8)
    passed = {os.path.basename(argv[i + 1]) for i, tok in enumerate(argv) if tok == "--bins"}
    assert passed == {"%s_bins" % name for name in binners}
    assert "metawrap_70_10_bins" not in passed
    assert argv[argv.index("-a") + 1] == str(assembly)
