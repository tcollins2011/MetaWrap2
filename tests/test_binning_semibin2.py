"""binning --semibin2: the command it builds, and the bins it hands over.

Two of these are guarding against failures that look like success, which is why they are asserted on
the command string rather than on an end-to-end run. SemiBin2 has a ``multi_easy_bin`` subcommand
that accepts exactly the same inputs as the one this module wants and interprets them differently -
it expects every contig name to carry a sample prefix, which MetaWrap2's assemblies do not - and it
gzips its bins by default, which leaves ``renumber_bins`` with nothing it recognises as a FASTA and
the bin set silently empty.
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


def _capture(monkeypatch):
    """Replace the module's run() with a recorder, returning the list it appends to."""
    calls = []

    def fake_run(cmd, **kw):
        calls.append((cmd, kw))

    monkeypatch.setattr(binning, "run", fake_run)
    return calls


def _work_with_bins(tmp_path, names):
    work = tmp_path / "work_files"
    produced = work / "semibin2_out" / binning.SEMIBIN2_BIN_DIR
    produced.mkdir(parents=True)
    for name in names:
        (produced / name).write_text(">c1\nACGT\n")
    (work / "sample_0.bam").write_text("")
    return work


def test_the_command_is_single_easy_bin_and_does_not_gzip_its_bins(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = _work_with_bins(tmp_path, ["bin.1.fa", "bin.2.fa", "bin.10.fa"])
    out = tmp_path / "out"
    binning._run_semibin2(_Args(str(out)), "e", str(work), "asm.fa")

    cmd = calls[0][0]
    assert "single_easy_bin" in cmd and "multi_easy_bin" not in cmd
    # Without this the bins arrive as .fa.gz and the move below finds no FASTA at all.
    assert "--compression none" in cmd
    # Fixed seed: SemiBin2 trains a net, so two runs over one input otherwise disagree about which
    # contigs form a bin, and a benchmark measures that noise as part of the difference.
    assert "--random-seed %d" % binning.SEMIBIN2_SEED in cmd
    # And the same reproducibility argument, plus the cost model the defaults are chosen on:
    # SemiBin2 would otherwise train on whatever GPU torch finds, which is neither a CPU-hour nor a
    # device the next machine necessarily has.
    assert "--engine cpu" in cmd and binning.DEFAULT_LEARNED_DEVICE == "cpu"

    # The bins are moved out and renamed like every other binner's: bin.10 after bin.2, and
    # nothing left behind in SemiBin2's own tree to be scored twice.
    assert sorted(os.listdir(out / "semibin2_bins")) == [bin_filename(i) for i in (1, 2, 3)]
    assert os.listdir(work / "semibin2_out" / binning.SEMIBIN2_BIN_DIR) == []
