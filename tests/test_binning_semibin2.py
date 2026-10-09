"""binning --semibin2: the command it builds, and the bins it hands over.

Two of these are guarding against failures that look like success, which is why they are asserted on
the command string rather than on an end-to-end run. SemiBin2 has a ``multi_easy_bin`` subcommand
that accepts exactly the same inputs as the one this module wants and interprets them differently -
it expects every contig name to carry a sample prefix, which MetaWrap2's assemblies do not - and it
gzips its bins by default, which leaves ``renumber_bins`` with nothing it recognises as a FASTA and
the bin set silently empty.
"""

import os

import pytest

from metawrap2.config import OPTIONAL_ENV_TOOLS, OPTIONAL_ENVS
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
    work = _work_with_bins(tmp_path, ["bin.1.fa"])
    binning._run_semibin2(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")

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


def test_it_runs_in_its_own_env_rather_than_the_binning_one(tmp_path, monkeypatch):
    """pytorch is why SemiBin2 is not in bins.yaml, so it cannot be run from that env."""
    calls = _capture(monkeypatch)
    work = _work_with_bins(tmp_path, ["bin.1.fa"])
    binning._run_semibin2(
        _Args(str(tmp_path / "out")), OPTIONAL_ENVS["semibin2"], str(work), "asm.fa"
    )
    assert calls[0][1]["env"] == "metawrap2-learned" != binning.CONDA_ENV
    # And doctor has to probe the executable's real, capitalised name inside that env - plus the
    # other two binners sharing it, since a half-built env should not report as healthy.
    assert OPTIONAL_ENV_TOOLS["semibin2"] == ["SemiBin2", "vamb", "run_comebin.sh"]


def test_the_bins_are_moved_out_and_renamed_like_every_other_binners(tmp_path, monkeypatch):
    """SemiBin2's ``bin.10.fa`` has to end up after its ``bin.2.fa``, not before it."""
    _capture(monkeypatch)
    work = _work_with_bins(tmp_path, ["bin.1.fa", "bin.2.fa", "bin.10.fa"])
    out = tmp_path / "out"
    binning._run_semibin2(_Args(str(out)), "e", str(work), "asm.fa")

    bins = out / "semibin2_bins"
    assert sorted(os.listdir(bins)) == [bin_filename(i) for i in (1, 2, 3)]
    # Nothing left behind in SemiBin2's own output tree to be scored twice.
    assert os.listdir(work / "semibin2_out" / binning.SEMIBIN2_BIN_DIR) == []


def test_producing_no_bins_is_an_error_rather_than_an_empty_bin_set(tmp_path, monkeypatch):
    """SemiBin2 exits 0 having written nothing when no cluster clears its 200 kb floor.

    Left alone that becomes an empty ``semibin2_bins``, which bin_refinement consolidates happily
    and which looks in the output exactly like a sample with no binnable genomes in it.
    """
    _capture(monkeypatch)
    work = _work_with_bins(tmp_path, [])
    with pytest.raises(SystemExit):
        binning._run_semibin2(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")


def test_a_missing_output_directory_says_so_instead_of_raising_oserror(tmp_path, monkeypatch):
    _capture(monkeypatch)
    work = tmp_path / "work_files"
    work.mkdir()
    with pytest.raises(SystemExit):
        binning._run_semibin2(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")


def test_semibin2_alone_satisfies_the_pick_a_binner_check(tmp_path, monkeypatch):
    """``--semibin2`` on its own is a complete request; the guard listed three binners.

    Under a dry run every tool is only printed, so this reaches the SemiBin2 command and stops -
    which is what makes it a test of the guard rather than of SemiBin2.
    """
    from metawrap2.command import runner

    monkeypatch.setattr(runner, "dry_run", True)
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\n%s\n" % ("ACGT" * 500))
    reads = []
    for mate in ("1", "2"):
        read = tmp_path / ("r_%s.fastq" % mate)
        read.write_text("@r1\nACGT\n+\n!!!!\n")
        reads.append(str(read))
    rc = binning.main(
        ["--semibin2", "-a", str(assembly), "-o", str(tmp_path / "out"), "-t", "2"] + reads
    )
    assert rc == 0
