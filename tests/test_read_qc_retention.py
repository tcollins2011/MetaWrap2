"""read_qc has to say how much of the library it threw away, and FastQC has to start.

Three defects a first-time user hit on the MetaHIT tutorial cohort (ERR011347/8/9, 2010-era 75 bp
reads), one test per behaviour:

* fastp's defaults discarded 58.5% of the reads whole and the module said nothing, so the run's
  nonpareil estimate asked for 20x more sequencing where the same sample trimmed with trim-galore
  asked for 4.06x. Retention is now reported always, and warned about when it is low.
* FastQC 0.13.0 dropped the short ``-o``, which killed the first step of the pipeline on every
  fresh env solve.
* FastQC was handed every thread on the host and answered with a warning about it.
"""

import json
import os
import shlex

from metawrap2.modules import read_qc

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "read_qc")


def _write_read(path):
    with open(path, "w") as fh:
        fh.write("@r\nACGT\n+\nIIII\n")


def _run_read_qc(tmp_path, monkeypatch, extra_argv=(), report=None, fastqc=False):
    """Run read_qc end to end with the trimmer and FastQC faked out.

    Returns (commands, said, warned): the commands that would have run, the ``comm`` lines and the
    ``warning`` lines. *report* is a fixture filename (or raw text) dropped at ``fastp.json`` in
    place of the report the real fastp would have written.

    The "<dir> already exists" warning is dropped from *warned*: the disk-space check creates the
    output directory before the module looks to see whether it was already there, so every run
    emits it and no assertion here is about it.
    """
    calls, said, warned = [], [], []

    def fake_run(cmd, env=None, tool=None, **kwargs):
        calls.append(cmd)
        parts = shlex.split(cmd)
        if tool == "fastp":
            for flag in ("-o", "-O"):
                if flag in parts:
                    _write_read(parts[parts.index(flag) + 1])
            if report is not None:
                text = report
                if os.path.isfile(os.path.join(FIXTURES, report)):
                    with open(os.path.join(FIXTURES, report)) as fh:
                        text = fh.read()
                with open(parts[parts.index("-j") + 1], "w") as fh:
                    fh.write(text)
        elif tool == "trim_galore":
            outdir = parts[parts.index("-o") + 1]
            for i, path in enumerate(parts[-2:]):
                _write_read(
                    os.path.join(outdir, "%s_val_%d.fq" % (read_qc._trim_stem(path), i + 1))
                )
            if report is not None:
                with open(os.path.join(FIXTURES, report)) as src:
                    text = src.read()
                with open(os.path.join(outdir, "%s_trimming_report.txt" % "r_1.fastq"), "w") as fh:
                    fh.write(text)

    monkeypatch.setattr(read_qc, "run", fake_run)
    monkeypatch.setattr(read_qc, "env_for", lambda module, settings: None)
    monkeypatch.setattr(read_qc, "comm", said.append)
    monkeypatch.setattr(read_qc, "warning", warned.append)

    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    _write_read(str(r1))
    _write_read(str(r2))
    out = tmp_path / "qc"
    argv = ["-1", str(r1), "-2", str(r2), "-o", str(out), "-t", "32", "--skip-bmtagger"]
    if not fastqc:
        argv += ["--skip-pre-qc-report", "--skip-post-qc-report"]
    assert read_qc.main(argv + list(extra_argv)) == 0
    return calls, said, [m for m in warned if not m.endswith("already exists.")]


def test_retention_line_is_read_out_of_fastps_own_json(tmp_path, monkeypatch):
    """Regresses the silent loss: the pairs-in/pairs-out line and its reason breakdown.

    The counts are fastp's, halved from reads to pairs, and the reason is the one fastp gave -
    58.5% filtered on quality with zero reads shortened.
    """
    _, said, _ = _run_read_qc(tmp_path, monkeypatch, report="fastp_low_retention.json")
    (line,) = [m for m in said if m.startswith("Read retention:")]
    assert line == (
        "Read retention: 12,621,104 pairs in -> 5,243,141 out (41.5%): "
        "58.5% removed by quality filtering"
    )


def test_low_retention_warns_and_names_the_trimmer_that_trims(tmp_path, monkeypatch):
    """Regresses a loss that is reported but not acted on: the warning must say what to do.

    It has to carry the measured fraction (not just "low"), name --trimmer trim-galore, and say
    that the difference shows up on older or lower-quality data.
    """
    _, _, warned = _run_read_qc(tmp_path, monkeypatch, report="fastp_low_retention.json")
    (message,) = warned
    assert "41.5%" in message and "5,243,141 of 12,621,104" in message
    assert "--trimmer trim-galore" in message
    assert "--fastp-quality-trim" in message
    assert "lower-quality data" in message


def test_high_retention_does_not_warn(tmp_path, monkeypatch):
    """Regresses a warning that cries wolf: a 99.1% run is reported, never warned about.

    The threshold has to leave normal libraries alone or it stops meaning anything.
    """
    _, said, warned = _run_read_qc(tmp_path, monkeypatch, report="fastp_high_retention.json")
    assert warned == []
    (line,) = [m for m in said if m.startswith("Read retention:")]
    assert "10,000,000 pairs in -> 9,910,000 out (99.1%)" in line


def test_a_missing_fastp_json_is_a_quiet_skip(tmp_path, monkeypatch):
    """Regresses a QC step failed by its own report: no fastp.json means no line, not an error."""
    _, said, warned = _run_read_qc(tmp_path, monkeypatch, report=None)
    assert [m for m in said if m.startswith("Read retention:")] == []
    assert warned == []
    assert (tmp_path / "qc" / "final_pure_reads_1.fastq").is_file()


def test_a_truncated_fastp_json_is_a_quiet_skip(tmp_path, monkeypatch):
    """Regresses the same failure from the other direction: a half-written report is not fatal.

    fastp killed mid-run leaves exactly this behind, and the reads it already wrote are fine.
    """
    _, said, warned = _run_read_qc(
        tmp_path, monkeypatch, report='{"summary": {"before_filtering": {"total_'
    )
    assert [m for m in said if m.startswith("Read retention:")] == []
    assert warned == []
    assert (tmp_path / "qc" / "final_pure_reads_1.fastq").is_file()


def test_trim_galore_retention_comes_from_its_trimming_report(tmp_path, monkeypatch):
    """Regresses a retention line that only works on the default path.

    Trim Galore has no JSON, so the pair-level totals are read out of the text report - which means
    --no_report_file must stay off the command.
    """
    calls, said, warned = _run_read_qc(
        tmp_path,
        monkeypatch,
        extra_argv=["--trimmer", "trim-galore"],
        report="trim_galore_paired_report.txt",
    )
    assert "--no_report_file" not in calls[0]
    (line,) = [m for m in said if m.startswith("Read retention:")]
    # 12,621,104 pairs analysed, 584,125 pairs dropped for length: 95.4% kept, so no warning.
    assert line == (
        "Read retention: 12,621,104 pairs in -> 12,036,979 out (95.4%): "
        "4.6% removed by length filtering"
    )
    assert warned == []


def test_fastp_quality_trim_puts_the_cut_flags_on_the_command(tmp_path, monkeypatch):
    """Regresses the opt-in that lets fastp trim instead of discard."""
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, extra_argv=["--fastp-quality-trim"])
    (cmd,) = calls
    assert cmd.endswith(" --cut_front --cut_tail --cut_mean_quality 20")


def test_without_the_flag_the_fastp_command_is_unchanged(tmp_path, monkeypatch):
    """Regresses a default moved by accident: the shipped command must stay byte-identical.

    Changing which reads fastp keeps takes a measured benchmark arm, so the flags only ever arrive
    with the opt-in - nothing may leak into the default invocation.
    """
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch)
    out = tmp_path / "qc"
    assert calls == [
        "fastp -i %s -I %s -o %s -O %s -w %d -j %s -h %s"
        % (
            tmp_path / "s_1.fastq",
            tmp_path / "s_2.fastq",
            out / "trimmed_1.fastq",
            out / "trimmed_2.fastq",
            read_qc.FASTP_MAX_THREADS,
            out / "fastp.json",
            out / "fastp.html",
        )
    ]


def test_fastqc_is_told_outdir_not_dash_o(tmp_path, monkeypatch):
    """Regresses the dead first step: FastQC 0.13.0 dropped -o and read_qc died on it.

    --outdir is accepted by 0.12.x as well, so one spelling has to serve both.
    """
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, fastqc=True)
    fastqc_calls = [c for c in calls if c.startswith("fastqc ")]
    assert len(fastqc_calls) == 2  # pre-QC and post-QC
    for cmd in fastqc_calls:
        parts = shlex.split(cmd)
        assert "--outdir" in parts
        assert "-o" not in parts


def test_fastqc_threads_are_capped(tmp_path, monkeypatch):
    """Regresses `[WARNING] You can't use more than 4 threads per file` on every read_qc run.

    -t 32 is what the user asked for; 4 is all FastQC will use.
    """
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, fastqc=True)
    for cmd in [c for c in calls if c.startswith("fastqc ")]:
        assert "-t %d" % read_qc.FASTQC_MAX_THREADS in cmd


def test_single_end_retention_counts_reads_not_pairs(tmp_path, monkeypatch):
    """Regresses a halved count on single-end data, where fastp's reads *are* the units."""
    with open(os.path.join(FIXTURES, "fastp_low_retention.json")) as fh:
        doc = json.load(fh)
    path = tmp_path / "fastp.json"
    path.write_text(json.dumps(doc))
    retention = read_qc._fastp_retention(str(path), paired=False)
    assert (retention.before, retention.after) == (25242208, 10486282)
    assert retention.reasons == (("quality filtering", 14755926),)
