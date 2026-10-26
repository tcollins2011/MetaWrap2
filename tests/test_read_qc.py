import os
import shlex

from metawrap2.modules import read_qc


def _write_read(path):
    with open(path, "w") as fh:
        fh.write("@r\nACGT\n+\nIIII\n")


def _fake_trimmer(tmp_path, monkeypatch, calls):
    """Fake trimmer: record the command and write its expected outputs."""

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
    """fastp writes to the given output paths, so no rename is attempted."""
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
    # Capped at fastp's useful thread ceiling instead of all 32 threads.
    assert "-w %d" % read_qc.FASTP_MAX_THREADS in cmd
    assert str(out / "trimmed_1.fastq") in cmd
    assert (out / "final_pure_reads_1.fastq").is_file()
    assert (out / "final_pure_reads_2.fastq").is_file()


def test_trim_galore_output_is_renamed_from_its_own_convention(tmp_path, monkeypatch):
    """The non-default branch still works: _val_1.fq is moved onto trimmed_1.fastq.

    Asked for explicitly, since fastp holds the default slot now.
    """
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
    # Capped at 4, where the underlying cutadapt stops scaling.
    assert "-j %d" % read_qc.TRIM_GALORE_MAX_THREADS in cmd
    assert (out / "final_pure_reads_1.fastq").is_file()
