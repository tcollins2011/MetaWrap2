"""plot_binning_results draws a figure when one bin set has no bins, instead of raising NameError."""

import os
import subprocess
import sys

from metawrap2.checkm import STATS_HEADER


def _stats(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(STATS_HEADER) + "\n")
        fh.writelines(
            "bin.%d\t%s\t%s\t0.5\troot\t10000\t2000000\n" % (i, comp, cont)
            for i, (comp, cont) in enumerate(rows)
        )
    return str(path)


def test_an_empty_bin_set_does_not_stop_the_plot(tmp_path):
    # The empty set comes first, so no x_pos is left over from an earlier set.
    empty = _stats(tmp_path / "binsA.stats.tsv", [])
    full = _stats(tmp_path / "binsB.stats.tsv", [(95, 1), (80, 3), (72, 8)])
    env = dict(os.environ, MPLBACKEND="Agg")
    result = subprocess.run(
        [sys.executable, "-m", "metawrap2.scripts.plot_binning_results", "70", "10", empty, full],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr[-2000:]
