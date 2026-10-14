"""The runner's command for a bin_refinement comparison: every input set, and the assembly."""

import os

from metawrap2.benchmark import runner as _runner
from metawrap2.benchmark.registry import Comparison, Variant


def test_every_bin_set_reaches_refinement_and_a_leftover_answer_does_not(tmp_path):
    """The sets were zipped onto ("-A", "-B", "-C"), and zip truncates.

    A variant that binned with six binners was therefore refined from the three whose directory
    names sort first - metabat2_bins among those dropped - so the arm's numbers came from a
    different input set than the one it was named for, with nothing said. A consolidated set left
    behind by an earlier refinement into the same directory is also not an input: refining it
    together with the sets it came from is not a comparison anyone declared.

    The arm is also handed the assembly: DAS Tool reads the assembly contigs itself and exits
    immediately without them, which looked like DAS Tool being broken rather than like the harness
    never passing ``-a``.
    """
    binners = ("comebin", "concoct", "maxbin2", "metabat2", "semibin2", "vamb")
    for name in binners:
        (tmp_path / ("%s_bins" % name)).mkdir()
    (tmp_path / "metawrap_70_10_bins").mkdir()  # a previous refinement's own answer
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
