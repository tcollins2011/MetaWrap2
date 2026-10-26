"""read_qc retention reporting and FastQC invocation.

Covers three issues found on the MetaHIT tutorial cohort (ERR011347/8/9, 2010-era 75 bp reads):

* fastp's defaults discarded 58.5% of the reads whole without a message; nonpareil then estimated
  20x more sequencing, against 4.06x for the same sample trimmed with trim-galore. Retention is now
  always reported, with a warning when it is low.
* FastQC 0.13.0 dropped the short ``-o``, so the first step failed on a fresh env solve.
* FastQC was given every available thread and warned about it.
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
    place of the report fastp would have written.

    The "<dir> already exists" warning is dropped from *warned*: the disk-space check creates the
    output directory first, so every run emits it.
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
    """The retention line gives pairs in and out and the reason breakdown from fastp's JSON.

    Counts are halved from reads to pairs; the reason is 58.5% quality filtering, none shortened.
    """
    _, said, _ = _run_read_qc(tmp_path, monkeypatch, report="fastp_low_retention.json")
    (line,) = [m for m in said if m.startswith("Read retention:")]
    assert line == (
        "Read retention: 12,621,104 pairs in -> 5,243,141 out (41.5%): "
        "58.5% removed by quality filtering"
    )


def test_low_retention_warns_and_names_the_trimmer_that_trims(tmp_path, monkeypatch):
    """The low-retention warning gives the fraction and suggests the alternatives.

    It names --trimmer trim-galore and --fastp-quality-trim and mentions lower-quality data.
    """
    _, _, warned = _run_read_qc(tmp_path, monkeypatch, report="fastp_low_retention.json")
    (message,) = warned
    assert "41.5%" in message and "5,243,141 of 12,621,104" in message
    assert "--trimmer trim-galore" in message
    assert "--fastp-quality-trim" in message
    assert "lower-quality data" in message


def test_high_retention_does_not_warn(tmp_path, monkeypatch):
    """A 99.1% retention run is reported without a warning."""
    _, said, warned = _run_read_qc(tmp_path, monkeypatch, report="fastp_high_retention.json")
    assert warned == []
    (line,) = [m for m in said if m.startswith("Read retention:")]
    assert "10,000,000 pairs in -> 9,910,000 out (99.1%)" in line


def test_a_missing_fastp_json_is_a_quiet_skip(tmp_path, monkeypatch):
    """Without fastp.json there is no retention line and no error."""
    _, said, warned = _run_read_qc(tmp_path, monkeypatch, report=None)
    assert [m for m in said if m.startswith("Read retention:")] == []
    assert warned == []
    assert (tmp_path / "qc" / "final_pure_reads_1.fastq").is_file()


def test_a_truncated_fastp_json_is_a_quiet_skip(tmp_path, monkeypatch):
    """A truncated fastp.json (as left by an interrupted fastp) is skipped without an error."""
    _, said, warned = _run_read_qc(
        tmp_path, monkeypatch, report='{"summary": {"before_filtering": {"total_'
    )
    assert [m for m in said if m.startswith("Read retention:")] == []
    assert warned == []
    assert (tmp_path / "qc" / "final_pure_reads_1.fastq").is_file()


def test_trim_galore_retention_comes_from_its_trimming_report(tmp_path, monkeypatch):
    """With trim-galore, retention is read from its text trimming report.

    So --no_report_file must not be on the command.
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
    """--fastp-quality-trim appends the cut flags to the fastp command."""
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, extra_argv=["--fastp-quality-trim"])
    (cmd,) = calls
    assert cmd.endswith(" --cut_front --cut_tail --cut_mean_quality 20")


def test_without_the_flag_the_fastp_command_is_unchanged(tmp_path, monkeypatch):
    """Without --fastp-quality-trim the default fastp command is unchanged."""
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
    """FastQC is called with --outdir, not -o, which 0.13.0 dropped (0.12.x accepts both)."""
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, fastqc=True)
    fastqc_calls = [c for c in calls if c.startswith("fastqc ")]
    assert len(fastqc_calls) == 2  # pre-QC and post-QC
    for cmd in fastqc_calls:
        parts = shlex.split(cmd)
        assert "--outdir" in parts
        assert "-o" not in parts


def test_fastqc_threads_are_capped(tmp_path, monkeypatch):
    """FastQC's -t is capped at FASTQC_MAX_THREADS when read_qc gets -t 32."""
    calls, _, _ = _run_read_qc(tmp_path, monkeypatch, fastqc=True)
    for cmd in [c for c in calls if c.startswith("fastqc ")]:
        assert "-t %d" % read_qc.FASTQC_MAX_THREADS in cmd


def test_single_end_retention_counts_reads_not_pairs(tmp_path, monkeypatch):
    """Single-end retention counts reads without halving."""
    with open(os.path.join(FIXTURES, "fastp_low_retention.json")) as fh:
        doc = json.load(fh)
    path = tmp_path / "fastp.json"
    path.write_text(json.dumps(doc))
    retention = read_qc._fastp_retention(str(path), paired=False)
    assert (retention.before, retention.after) == (25242208, 10486282)
    assert retention.reasons == (("quality filtering", 14755926),)
