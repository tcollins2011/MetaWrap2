"""binning --vamb: the command it builds, and the bins it hands over.

Three of VAMB's defaults are wrong for this module in ways that do not look like failures, and each
has a test here. ``--minfasta`` has no default at all, so without it VAMB writes a cluster TSV and
zero FASTA files - a binner that exited 0 and found nothing. ``-o`` defaults to the separator ``C``,
which splits every bin apart on a letter that appears mid-identifier in ordinary contig names. And
VAMB refuses to run into an existing ``--outdir``, so a ``--resume`` after an interrupted VAMB would
fail on its own leftovers rather than on anything real.
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
    """Stub run() to write *names* where VAMB would, and record the calls it was asked to make.

    The bins have to be created *by* the stub rather than before it: VAMB creates its own --outdir,
    and the module clears any leftover one first, so anything staged in advance is correctly gone.
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
    # VAMB raises FileExistsError on its own --outdir, which would break every --resume; the
    # checkpoint has already decided this stage is to be redone, so the leftovers are cleared.
    (work / "vamb_out").mkdir()
    stale = work / "vamb_out" / "stale.log"
    stale.write_text("from the attempt that died")
    calls = _capture(monkeypatch, work, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")
    assert not stale.exists()
    assert os.listdir(tmp_path / "out" / "vamb_bins") == [bin_filename(1)]

    cmd = calls[0][0]
    # The flag with no default: without it there are no FASTA files to adopt, only a cluster TSV.
    assert "--minfasta %d" % binning.VAMB_MIN_FASTA in cmd
    # VAMB's separator default is "C", and MetaWrap2's assemblies carry no sample prefix to split.
    assert "-o ''" in cmd
    assert "--seed %d" % binning.VAMB_SEED in cmd
    # It reads the BAM directory itself rather than a depth table, so work/ is handed over whole.
    assert "--bamdir %s" % work in cmd
    # No --cuda: unlike SemiBin2, VAMB stays on the CPU unless asked, which is what the cost model
    # that picks MetaWrap2's defaults needs.
    assert "--cuda" not in cmd
