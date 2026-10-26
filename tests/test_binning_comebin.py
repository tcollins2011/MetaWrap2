"""Tests for binning --comebin: the command, the input filter, and the adopted bins.

COMEBin seeds nothing unless asked, and its device default is ``cuda`` (SemiBin2 uses ``auto``),
so it fails on a machine without a GPU. Both are set explicitly. ``run_comebin.sh`` does not
forward ``--contig_len``, so ``-l`` is applied by filtering the input.
"""

import os

from metawrap2.modules import binning


class _Args:
    def __init__(self, output, threads=4, min_len=1000, device=binning.DEFAULT_LEARNED_DEVICE):
        self.output = output
        self.threads = threads
        self.min_len = min_len
        self.device = device
        self.assembly = "asm.fa"


def _capture(monkeypatch, work, names=("1.fa",)):
    """Stub run() to write *names* where COMEBin publishes its winning bin set."""
    calls = []
    produced = os.path.join(str(work), "comebin_out", binning.COMEBIN_BIN_DIR)

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


def test_the_command_pins_a_seed_and_the_cpu(tmp_path, monkeypatch):
    work = _work(tmp_path)
    calls = _capture(monkeypatch, work)
    binning._run_comebin(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")

    cmd = calls[0][0]
    # COMEBin seeds only when -s is given.
    assert "-s %d" % binning.COMEBIN_SEED in cmd
    # Its device default is cuda, with no CPU fallback.
    assert "-d cpu" in cmd and binning.DEFAULT_LEARNED_DEVICE == "cpu"
    # The entry point is the shell script; there is no "comebin" executable.
    assert cmd.startswith("run_comebin.sh ")
    # It globs *.bam from a directory, like VAMB, so work_files is passed as a whole.
    assert "-p %s" % work in cmd
