"""binning --comebin: the command it builds, the input it filters, and the bins it adopts.

COMEBin's two defaults that cannot be left alone are the opposite kind of problem from VAMB's. It
seeds nothing unless asked, so two runs over one assembly disagree; and its device default is plain
``cuda`` rather than SemiBin2's ``auto``, so on a machine without a GPU the default does not degrade
to the CPU, it fails. The third thing tested here is not a default at all but a missing flag:
``run_comebin.sh`` never forwards ``--contig_len``, so ``-l`` is honoured by filtering the input.
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


def test_a_higher_min_len_filters_the_assembly_because_comebin_has_no_flag_for_it(
    tmp_path, monkeypatch
):
    """``-l 2000`` has to reach COMEBin somehow, and the only lever is the fasta it is handed."""
    work = _work(tmp_path)
    assembly = work / "assembly.fa"
    assembly.write_text(">short\n%s\n>long\n%s\n" % ("A" * 500, "C" * 3000))
    calls = _capture(monkeypatch, work)
    binning._run_comebin(_Args(str(tmp_path / "out"), min_len=2000), "e", str(work), str(assembly))

    filtered = work / "assembly_comebin.fa"
    assert "-a %s" % filtered in calls[0][0]
    assert filtered.read_text() == ">long\n%s\n" % ("C" * 3000)


def test_the_shared_assembly_copy_is_used_when_min_len_is_at_comebins_own_floor(
    tmp_path, monkeypatch
):
    """Nothing is written in the common case - the filter exists only to keep a promise -l made."""
    work = _work(tmp_path)
    calls = _capture(monkeypatch, work)
    binning._run_comebin(
        _Args(str(tmp_path / "out"), min_len=binning.COMEBIN_MIN_LEN), "e", str(work), "shared.fa"
    )
    assert "-a shared.fa" in calls[0][0]
    assert not (work / "assembly_comebin.fa").exists()


def test_filtering_away_every_contig_is_an_error_rather_than_an_empty_fasta(tmp_path, monkeypatch):
    """An empty fasta would reach COMEBin and fail inside it, several stages from the cause."""
    work = _work(tmp_path)
    assembly = work / "assembly.fa"
    assembly.write_text(">short\n%s\n" % ("A" * 500))
    _capture(monkeypatch, work)
    with pytest.raises(SystemExit):
        binning._run_comebin(
            _Args(str(tmp_path / "out"), min_len=5000), "e", str(work), str(assembly)
        )


def test_only_the_published_bin_set_is_adopted_not_the_candidates(tmp_path, monkeypatch):
    """COMEBin leaves several clusterings under cluster_res and picks between them; one wins."""
    work = _work(tmp_path)
    _capture(monkeypatch, work, names=("1.fa", "2.fa", "10.fa"))
    candidates = work / "comebin_out" / "comebin_res" / "cluster_res"
    candidates.mkdir(parents=True)
    (candidates / "weight_seed_kmeans_k_30_result.tsv_bins").mkdir()
    out = tmp_path / "out"
    binning._run_comebin(_Args(str(out)), "e", str(work), "asm.fa")

    assert sorted(os.listdir(out / "comebin_bins")) == [bin_filename(i) for i in (1, 2, 3)]
    assert os.listdir(work / "comebin_out" / binning.COMEBIN_BIN_DIR) == []
    # The candidate directory is left exactly where COMEBin put it, not adopted alongside.
    assert (candidates / "weight_seed_kmeans_k_30_result.tsv_bins").is_dir()


def test_it_runs_in_the_one_env_the_three_learned_binners_share(tmp_path, monkeypatch):
    """One env for all three, and not the binning env. See the VAMB test for why they merged."""
    work = _work(tmp_path)
    calls = _capture(monkeypatch, work)
    binning._run_comebin(_Args(str(tmp_path / "out")), OPTIONAL_ENVS["comebin"], str(work), "a.fa")
    assert calls[0][1]["env"] == "metawrap2-learned" != binning.CONDA_ENV
    assert len({OPTIONAL_ENVS[k] for k in ("semibin2", "vamb", "comebin")}) == 1
    assert OPTIONAL_ENV_TOOLS["comebin"] == ["SemiBin2", "vamb", "run_comebin.sh"]


def test_comebin_alone_satisfies_the_pick_a_binner_check(tmp_path, monkeypatch):
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
        ["--comebin", "-a", str(assembly), "-o", str(tmp_path / "out"), "-t", "2"] + reads
    )
    assert rc == 0
