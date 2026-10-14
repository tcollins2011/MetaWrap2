"""binning --comebin: the command it builds, the input it filters, and the bins it adopts.

COMEBin's two defaults that cannot be left alone are the opposite kind of problem from VAMB's. It
seeds nothing unless asked, so two runs over one assembly disagree; and its device default is plain
``cuda`` rather than SemiBin2's ``auto``, so on a machine without a GPU the default does not degrade
to the CPU, it fails. The third thing tested here is not a default at all but a missing flag:
``run_comebin.sh`` never forwards ``--contig_len``, so ``-l`` is honoured by filtering the input.
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
    # COMEBin forwards --seed only when -s is given, so without this nothing is seeded at all.
    assert "-s %d" % binning.COMEBIN_SEED in cmd
    # And its device default is cuda outright, not "a GPU if there is one".
    assert "-d cpu" in cmd and binning.DEFAULT_LEARNED_DEVICE == "cpu"
    # The entry point is the shell script; there is no "comebin" executable to call instead.
    assert cmd.startswith("run_comebin.sh ")
    # It globs *.bam out of a directory, like VAMB, so work_files goes over whole.
    assert "-p %s" % work in cmd
