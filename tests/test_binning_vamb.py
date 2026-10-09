"""binning --vamb: the command it builds, and the bins it hands over.

Three of VAMB's defaults are wrong for this module in ways that do not look like failures, and each
has a test here. ``--minfasta`` has no default at all, so without it VAMB writes a cluster TSV and
zero FASTA files - a binner that exited 0 and found nothing. ``-o`` defaults to the separator ``C``,
which splits every bin apart on a letter that appears mid-identifier in ordinary contig names. And
VAMB refuses to run into an existing ``--outdir``, so a ``--resume`` after an interrupted VAMB would
fail on its own leftovers rather than on anything real.
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
    calls = _capture(monkeypatch, work, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")

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


def test_the_minimum_contig_length_is_floored_at_vambs_own_limit(tmp_path, monkeypatch):
    """VAMB rejects anything under 250, so a user's -l 100 has to be raised rather than passed."""
    work = _work(tmp_path)
    calls = _capture(monkeypatch, work, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out"), min_len=100), "e", str(work), "asm.fa")
    assert "-m %d" % binning.VAMB_MIN_CONTIG_LEN in calls[0][0]

    work2 = _work(tmp_path / "second")
    calls2 = _capture(monkeypatch, work2, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out2"), min_len=2500), "e", str(work2), "asm.fa")
    assert "-m 2500" in calls2[0][0]


def test_an_existing_output_directory_is_cleared_rather_than_tripped_over(tmp_path, monkeypatch):
    """VAMB raises FileExistsError on its own --outdir, which would break every --resume.

    The directory is only cleared because the checkpoint has already decided this stage is to be
    redone, so nothing in there is a result anyone wanted.
    """
    work = _work(tmp_path)
    (work / "vamb_out").mkdir()
    stale = work / "vamb_out" / "stale.log"
    stale.write_text("from the attempt that died")
    _capture(monkeypatch, work, ["cluster_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")
    assert not stale.exists()
    assert os.listdir(tmp_path / "out" / "vamb_bins") == [bin_filename(1)]


def test_fna_bins_are_moved_out_and_renamed_like_every_other_binners(tmp_path, monkeypatch):
    """VAMB writes ``.fna``, which has to be recognised and renumbered like metaBAT2's ``.fa``."""
    work = _work(tmp_path)
    _capture(monkeypatch, work, ["c_1.fna", "c_2.fna", "c_10.fna"])
    out = tmp_path / "out"
    binning._run_vamb(_Args(str(out)), "e", str(work), "asm.fa")

    assert sorted(os.listdir(out / "vamb_bins")) == [bin_filename(i) for i in (1, 2, 3)]
    # Moved, not copied: a bin left in both trees is scored twice by anything walking the output.
    assert os.listdir(work / "vamb_out" / binning.VAMB_BIN_DIR) == []


def test_it_runs_in_the_learned_binner_env_not_the_binning_one(tmp_path, monkeypatch):
    """pytorch is why VAMB is not in bins.yaml - and why it shares an env with the other two.

    The three learned binners were once three envs. They solve together, so that was three copies of
    the same ~1.5 GB graph; what has to stay true is only that they are not in the *binning* env,
    which the three default binners share and which never touches pytorch.
    """
    work = _work(tmp_path)
    calls = _capture(monkeypatch, work, ["c_1.fna"])
    binning._run_vamb(_Args(str(tmp_path / "out")), OPTIONAL_ENVS["vamb"], str(work), "asm.fa")
    assert calls[0][1]["env"] == "metawrap2-learned" != binning.CONDA_ENV
    assert OPTIONAL_ENVS["vamb"] == OPTIONAL_ENVS["semibin2"] == OPTIONAL_ENVS["comebin"]
    # doctor probes all three executables, so a half-built env reports as broken rather than fine.
    assert OPTIONAL_ENV_TOOLS["vamb"] == ["SemiBin2", "vamb", "run_comebin.sh"]


def test_producing_no_bins_is_an_error_rather_than_an_empty_bin_set(tmp_path, monkeypatch):
    """Writing the TSV and no FASTAs is VAMB's likeliest silent failure - see the docstring."""
    work = _work(tmp_path)
    calls = []

    def tsv_only(cmd, **kw):
        calls.append(cmd)
        os.makedirs(os.path.join(str(work), "vamb_out"))
        with open(os.path.join(str(work), "vamb_out", "vae_clusters.tsv"), "w") as fh:
            fh.write("cluster\tcontig\n")

    monkeypatch.setattr(binning, "run", tsv_only)
    with pytest.raises(SystemExit):
        binning._run_vamb(_Args(str(tmp_path / "out")), "e", str(work), "asm.fa")


def test_vamb_alone_satisfies_the_pick_a_binner_check(tmp_path, monkeypatch):
    """``--vamb`` on its own is a complete request, and the guard has to list it."""
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
        ["--vamb", "-a", str(assembly), "-o", str(tmp_path / "out"), "-t", "2"] + reads
    )
    assert rc == 0
