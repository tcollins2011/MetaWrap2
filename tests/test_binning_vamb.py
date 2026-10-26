"""Tests for binning --vamb: the command line and the bins it produces.

Three VAMB defaults do not suit this module:

* ``--minfasta`` has no default, and without it VAMB writes a cluster TSV and no FASTA files.
* ``-o`` defaults to the separator ``C``, which splits bins on a letter found inside ordinary
  contig names.
* VAMB refuses an existing ``--outdir``, so ``--resume`` after an interrupted VAMB would fail on its
  leftovers.
"""

import os

from metawrap2.constants import bin_filename
from metawrap2.modules import binning


class _Args:
    def __init__(self, output, threads=4, min_len=1000, device=binning.DEFAULT_LEARNED_DEVICE):
        self.output = output
        self.threads = threads
        self.min_len = min_len
        self.device = device


def _capture(monkeypatch, work, names=()):
    """Stub run() to write *names* where VAMB would, and record the calls.

    The stub creates the bins because the module clears any existing --outdir before running VAMB.
    """
    calls = []
    produced = os.path.join(str(work), "vamb_out", binning.VAMB_BIN_DIR)

    def fake_run(cmd, **kw):
        calls.append((cmd, kw))
        os.makedirs(produced, exist_ok=True)
        for name in names:
            with open(os.path.join(produced, name), "w") as fh:
                fh.write(">c1\nACGT\n")

    monkeypatch.setattr(binning, "run", fake_run)
    return calls


def _work(tmp_path):
    work = tmp_path / "work_files"
    work.mkdir(parents=True)
    (work / "sample_0.bam").write_text("")
    return work


def test_the_command_asks_for_fastas_and_turns_binsplitting_off(tmp_path, monkeypatch):
    work = _work(tmp_path)
    # VAMB raises FileExistsError on an existing --outdir. The checkpoint has marked this stage
    # for a rerun, so the leftovers are cleared.
    (work / "vamb_out").mkdir()
    stale = work / "vamb_out" / "stale.log"
    stale.write_text("from the attempt that died")
    calls = _capture(monkeypatch, work, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")
    assert not stale.exists()
    assert os.listdir(tmp_path / "out" / "vamb_bins") == [bin_filename(1)]

    cmd = calls[0][0]
    # No default; without it VAMB writes only a cluster TSV.
    assert "--minfasta %d" % binning.VAMB_MIN_FASTA in cmd
    # VAMB's default separator is "C"; MetaWrap2's assemblies have no sample prefix to split on.
    assert "-o ''" in cmd
    assert "--seed %d" % binning.VAMB_SEED in cmd
    # VAMB reads the BAM directory directly instead of a depth table.
    assert "--bamdir %s" % work in cmd
    # No --cuda: VAMB stays on the CPU unless asked, which the cost model behind the defaults
    # assumes.
    assert "--cuda" not in cmd
