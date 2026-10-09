import os
import shlex

from metawrap2.modules import read_qc


def test_command_templates_format_cleanly():
    cmd = read_qc.FASTQC.format(threads=4, outdir="out/pre", files="a_1.fastq a_2.fastq")
    assert "{" not in cmd and "fastqc -q -t 4" in cmd
    # FastQC takes any number of files, so single-end works through the same template
    cmd = read_qc.FASTQC.format(threads=1, outdir="out/pre", files="a.fastq")
    assert "{" not in cmd and cmd.endswith("a.fastq")

    cmd = read_qc.TRIM_GALORE.format(outdir="out", r1="a_1.fastq", r2="a_2.fastq", trim_threads=4)
    assert "{" not in cmd and "--paired" in cmd and "-j 4" in cmd

    cmd = read_qc.BMTAGGER.format(
        bitmask="db/hg38.bitmask",
        srprism="db/hg38.srprism",
        tmpdir="out/tmp",
        r1="a_1.fastq",
        r2="a_2.fastq",
        listfile="out/a.bmtagger.list",
    )
    assert "{" not in cmd and "-q1" in cmd

    # single-end variants
    cmd = read_qc.TRIM_GALORE_SINGLE.format(outdir="out", reads="a.fastq", trim_threads=2)
    assert "{" not in cmd and "--paired" not in cmd and "-j 2" in cmd
    cmd = read_qc.BMTAGGER_SINGLE.format(
        bitmask="b", srprism="s", tmpdir="t", reads="a.fastq", listfile="l"
    )
    assert "{" not in cmd and "-1 a.fastq" in cmd and "-2" not in cmd


def _fastq(names):
    lines = []
    for n in names:
        lines += ["@%s/1" % n, "ACGT", "+", "IIII"]
    return "\n".join(lines) + "\n"


def test_fastp_templates_format_cleanly():
    cmd = read_qc.FASTP.format(
        r1="a_1.fastq",
        r2="a_2.fastq",
        o1="out/trimmed_1.fastq",
        o2="out/trimmed_2.fastq",
        trim_threads=16,
        json="out/fastp.json",
        html="out/fastp.html",
    )
    assert "{" not in cmd and "-w 16" in cmd and "-I a_2.fastq" in cmd
    cmd = read_qc.FASTP_SINGLE.format(
        reads="a.fastq",
        o1="out/trimmed_1.fastq",
        trim_threads=2,
        json="out/fastp.json",
        html="out/fastp.html",
    )
    assert "{" not in cmd and "-I " not in cmd and "-O " not in cmd


def test_fastp_is_the_default_the_benchmark_chose():
    """The default moved on evidence, not on expectation.

    `benchmark run cami2_gi trimmer` over the 10 CAMI2 GI samples put fastp ahead on both halves of
    the ratio at once - marginally more reads, bases and Q30 bases kept, for fewer CPU-hours - and
    that is what unseated trim-galore. The rows are in benchmarks/results/cami2_gi.tsv, so this
    assertion is checkable rather than a preference: re-measure and it either holds or it does not.
    """
    assert read_qc.DEFAULT_TRIMMER == "fastp"
    args = read_qc._parse_args(["-1", "a_1.fastq", "-2", "a_2.fastq", "-o", "out"])
    assert args.trimmer == "fastp"
    # The loser stays reachable: it is lighter on memory, which is a real reason to ask for it.
    other = read_qc._parse_args(
        ["-1", "a_1.fastq", "-2", "a_2.fastq", "-o", "out", "--trimmer", "trim-galore"]
    )
    assert other.trimmer == "trim-galore"


def _write_read(path):
    with open(path, "w") as fh:
        fh.write("@r\nACGT\n+\nIIII\n")


def _fake_trimmer(tmp_path, monkeypatch, calls):
    """Stand in for the real trimmer: record the command and write what it promises to write."""

    def fake_run(cmd, env=None, tool=None, **kwargs):
        calls.append(cmd)
        parts = shlex.split(cmd)
        if tool == "fastp":
            for flag in ("-o", "-O"):
                if flag in parts:
                    _write_read(parts[parts.index(flag) + 1])
        else:
            outdir = parts[parts.index("-o") + 1]
            for i, path in enumerate(parts[-2:]):
                stem = read_qc._trim_stem(path)
                _write_read(os.path.join(outdir, "%s_val_%d.fq" % (stem, i + 1)))

    monkeypatch.setattr(read_qc, "run", fake_run)
    monkeypatch.setattr(read_qc, "env_for", lambda module, settings: None)


def test_fastp_writes_straight_to_the_final_names(tmp_path, monkeypatch):
    """fastp takes output paths, so its result needs no rename - and none is attempted."""
    from metawrap2 import limits

    monkeypatch.setattr(limits, "machine_cores", lambda: 64)  # -t 32 must not be clamped first
    calls = []
    _fake_trimmer(tmp_path, monkeypatch, calls)
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    r2.write_text("@r\nTGCA\n+\nIIII\n")
    out = tmp_path / "qc"

    rc = read_qc.main(
        [
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(out),
            "-t",
            "32",
            "--trimmer",
            "fastp",
            "--skip-bmtagger",
            "--skip-pre-qc-report",
            "--skip-post-qc-report",
        ]
    )
    assert rc == 0
    (cmd,) = calls
    assert cmd.startswith("fastp ")
    # Capped at fastp's own useful ceiling rather than handed all 32 threads.
    assert "-w %d" % read_qc.FASTP_MAX_THREADS in cmd
    assert str(out / "trimmed_1.fastq") in cmd
    assert (out / "final_pure_reads_1.fastq").is_file()
    assert (out / "final_pure_reads_2.fastq").is_file()


def test_trim_galore_output_is_renamed_from_its_own_convention(tmp_path, monkeypatch):
    """The non-default branch still works: _val_1.fq is moved onto trimmed_1.fastq.

    Asked for explicitly, since fastp holds the default slot now.
    """
    calls = []
    _fake_trimmer(tmp_path, monkeypatch, calls)
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    r2.write_text("@r\nTGCA\n+\nIIII\n")
    out = tmp_path / "qc"

    rc = read_qc.main(
        [
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(out),
            "-t",
            "32",
            "--skip-bmtagger",
            "--skip-pre-qc-report",
            "--skip-post-qc-report",
            "--trimmer",
            "trim-galore",
        ]
    )
    assert rc == 0
    (cmd,) = calls
    assert cmd.startswith("trim_galore ")
    # Still capped at 4: cutadapt underneath stops scaling there whatever -t says.
    assert "-j %d" % read_qc.TRIM_GALORE_MAX_THREADS in cmd
    assert (out / "final_pure_reads_1.fastq").is_file()
