"""The contig depth table, computed here rather than by jgi_summarize_bam_contig_depths.

These cover the arithmetic and the two places it went wrong when checked against the tool it replaced.
The agreement itself is not testable without real BAMs - it was verified separately, over all 210,632
contigs of the CAMI2 GI assembly, at Pearson 0.999 on the means and 0.995 on the variances - so what is
pinned here is the behaviour that agreement depends on.
"""

from __future__ import annotations

import numpy as np
import pytest

from metawrap2 import depth


def _lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def test_the_divisor_is_the_trimmed_length_not_the_whole_contig():
    """75 bp comes off each end, because coverage falls off where reads cannot hang over the edge."""
    widths = depth.trimmed_widths(np.array([1000, 10000]), 75)
    assert list(widths) == [850, 9850]


@pytest.mark.parametrize("length", [1, 10, 150])
def test_a_contig_too_short_to_trim_is_measured_whole(length):
    """Trimming 150 bp off a 150 bp contig would leave nothing to average over.

    Measuring none of it is not an option - the binner was asked to bin it - so the whole contig is
    used. Dropping it would silently lose a contig the user asked for.
    """
    assert depth.trimmed_widths(np.array([length]), 75)[0] == length


def test_a_contig_just_over_the_trim_threshold_keeps_one_base():
    """151 bp is the first length with a middle, and it has exactly one base of it."""
    assert depth.trimmed_widths(np.array([151]), 75)[0] == 1


def test_the_identity_filter_is_expressed_over_nm_and_query_length():
    """97% identity, as a multiplication so a zero-length query cannot divide by zero.

    Over ``[NM]`` rather than over a ratio of its own so that a BAM with no NM tags fails at the
    filter instead of quietly counting every read.
    """
    assert depth._identity_expression(97.0) == ["-e", "[NM] <= (qlen * 0.03)"]
    assert depth._identity_expression(100.0) == ["-e", "[NM] <= (qlen * 0)"]
    # 0 turns the filter off rather than demanding perfection of every read.
    assert depth._identity_expression(0) == []


def test_the_excluded_flags_drop_secondary_alignments_but_not_multimappers():
    """A repeat shared between two genomes is real coverage of both; the binner decides, not this."""
    unmapped, secondary, qcfail, duplicate = 0x4, 0x100, 0x200, 0x400
    assert depth.EXCLUDED_FLAGS == unmapped | secondary | qcfail | duplicate
    # MAPQ 0 is kept: raising it would throw away the primary record of a multi-mapping read.
    assert depth.MIN_MAP_QUAL == 0


def test_contig_lengths_refuses_a_header_with_no_sq_lines(monkeypatch):
    """A BAM naming no contigs is not something to average over - it is a broken input.

    This is the error that caught the real bug: a stub BAM produced no header, and the alternative to
    raising here was a table of contigs that came from somewhere else.
    """

    class _Done:
        stdout = "@HD\tVN:1.6\n"

    monkeypatch.setattr(depth.subprocess, "run", lambda *a, **k: _Done())
    with pytest.raises(ValueError, match="no @SQ lines"):
        depth.contig_lengths("x.bam")


def test_contig_lengths_reads_names_and_lengths_in_header_order(monkeypatch):
    class _Done:
        stdout = "@HD\tVN:1.6\n@SQ\tSN:c2\tLN:200\n@SQ\tSN:c1\tLN:100\n@PG\tID:bwa\n"

    monkeypatch.setattr(depth.subprocess, "run", lambda *a, **k: _Done())
    # Order matters: it decides the row order of the table, and two tables built from the same BAMs
    # have to be comparable row by row.
    assert list(depth.contig_lengths("x.bam").items()) == [("c2", 200), ("c1", 100)]


def test_sample_names_default_to_the_bam_filenames():
    assert depth._sample_name("/a/b/sample_10.bam") == "sample_10"


def test_a_mismatched_sample_name_count_is_refused(tmp_path):
    with pytest.raises(ValueError, match="correspond one to one"):
        depth.write_depth_table(["a.bam", "b.bam"], str(tmp_path / "d.txt"), samples=["only_one"])


def test_no_bams_is_an_error_rather_than_an_empty_table(tmp_path):
    with pytest.raises(ValueError, match="nothing to measure"):
        depth.write_depth_table([], str(tmp_path / "d.txt"))


def test_the_table_layout_and_arithmetic(tmp_path):
    """Mean, variance and totalAvgDepth from known per-base depths, in the documented column order."""
    order = ["c1", "c2"]
    lengths = {"c1": 1000, "c2": 400}
    names = ["s1", "s2"]
    widths = depth.trimmed_widths(np.array([1000, 400]), 75)  # 850, 250

    # c1: depth 2 everywhere in s1 (variance 0), and in s2 a total of 850 with sum of squares 1700,
    # which is mean 1 and variance 1. c2: all zeros, which has to come out as zeros and not NaN.
    total = np.array([[2.0 * 850, 850.0], [0.0, 0.0]])
    total_sq = np.array([[4.0 * 850, 1700.0], [0.0, 0.0]])

    out = str(tmp_path / "d.txt")
    depth._write(out, order, lengths, names, total, total_sq, widths)
    rows = [ln.split("\t") for ln in _lines(out)]

    assert rows[0] == ["contigName", "contigLen", "totalAvgDepth", "s1", "s1-var", "s2", "s2-var"]
    # contigLen is the real length, not the trimmed one that was averaged over.
    assert rows[1][:3] == ["c1", "1000", "3"]  # totalAvgDepth is the sum of the per-sample means
    assert rows[1][3:] == ["2", "0", "1", "1"]
    assert rows[2] == ["c2", "400", "0", "0", "0", "0", "0"]


def test_a_variance_that_rounds_negative_is_written_as_zero(tmp_path):
    """E[x^2] - E[x]^2 goes slightly negative on a constant column, and metaBAT2 would read it."""
    widths = np.array([1000.0])
    total = np.array([[1000.0]])
    # Nudged below the exact square so the subtraction underflows, as it does on real data.
    total_sq = np.array([[1000.0 * (1.0 - 1e-16)]])
    out = str(tmp_path / "d.txt")
    depth._write(out, ["c1"], {"c1": 1075}, ["s1"], total, total_sq, widths)
    variance = float(_lines(out)[1].split("\t")[4])
    assert variance == 0.0


def test_a_contig_no_sample_covered_divides_by_its_length_not_by_one(tmp_path):
    """The bug the comparison against jgi caught, pinned.

    ``samtools depth -a`` says nothing at all about a reference with no reads, so a divisor counted
    from emitted positions was 1 for such a contig - and any sample that *did* cover it then reported
    the sum over the whole contig as its depth. One 4 Mbp contig came out at 1.6 million.
    """
    widths = depth.trimmed_widths(np.array([4000000]), 75)
    total = np.array([[4000000.0]])  # a sample that does cover it, depth ~1
    total_sq = np.array([[4000000.0]])
    out = str(tmp_path / "d.txt")
    depth._write(out, ["big"], {"big": 4000000}, ["s1"], total, total_sq, widths)
    mean = float(_lines(out)[1].split("\t")[3])
    assert mean == pytest.approx(1.0, rel=1e-3)


# --- the identity floor, which is the one place this departs from jgi's defaults ---------------


def test_the_default_identity_floor_is_90_not_jgis_97():
    """Measured, not preferred: 97 costs 26% of the clean bins on CAMI2 GI.

    One set of alignments, six depth tables from it, metaBAT2 on each, AMBER and GUNC on the results.
    At a 90% floor: 88 bins against 74, higher purity, and 6.8% chimeric against 12.2% - so 82 clean
    bins against 65. The objection to loosening it is that a relative's reads inflate a contig's depth
    and the binner merges two strains; both AMBER purity and GUNC say that is not what happens here.
    """
    assert depth.PERCENT_IDENTITY == 90.0
    assert depth.JGI_PERCENT_IDENTITY == 97.0, "kept named, for reproducing a jgi table"


def test_a_looser_floor_keeps_strictly_more_reads():
    """The filter is monotonic in the threshold, which is what makes the sweep a sweep."""
    kept = {}
    for floor in (0, 90, 95, 97, 99):
        expr = depth._identity_expression(float(floor))
        kept[floor] = expr
    assert kept[0] == [], "no filter at all, rather than a filter that passes everything"
    # Every other floor is an expression whose allowance shrinks as the floor rises.
    allowances = [float(kept[f][1].split("* ")[1].rstrip(")")) for f in (90, 95, 97, 99)]
    assert allowances == sorted(allowances, reverse=True)
    assert allowances[0] == 0.1, "a 90% floor allows 10% of the query length in mismatches"


def test_the_floor_is_an_argument_so_a_jgi_table_is_still_reachable():
    """A changed default must not become a changed capability."""
    import inspect

    sig = inspect.signature(depth.write_depth_table)
    assert sig.parameters["percent_identity"].default == depth.PERCENT_IDENTITY


def test_a_chatty_samtools_does_not_deadlock_the_depth_stream(tmp_path, monkeypatch):
    """More stderr than a pipe buffer holds, while stdout is still being read.

    Both stderrs used to be left unread until after the whole depth stream had been consumed. These
    are per-read tools over hundreds of millions of reads, so one warning per record fills the
    64 KiB pipe buffer; samtools then blocks on the write, stops producing rows, and the read of its
    stdout blocks forever. The run hung with no output and no timeout, which is why this uses a
    real subprocess rather than a fake: the deadlock is in the pipes, not in the logic.
    """
    import sys

    noisy = tmp_path / "noisy_samtools.py"
    noisy.write_text(
        "import sys\n"
        "# 'view' writes only to stderr; 'depth' writes the rows and a flood of warnings.\n"
        "mode = 'depth' if 'depth' in sys.argv else 'view'\n"
        "sys.stderr.write(('[%s] a warning about a read\\n' % mode) * 4000)\n"
        "if mode == 'depth':\n"
        "    for pos in range(1, 101):\n"
        "        sys.stdout.write('c1\\t%d\\t3\\n' % pos)\n"
        "sys.stderr.flush()\n"
    )
    monkeypatch.setattr(depth, "_samtools", lambda env: [sys.executable, str(noisy)])

    total = np.zeros((1, 1))
    total_sq = np.zeros((1, 1))
    depth._accumulate_bam(
        str(tmp_path / "s.bam"),
        0,
        {"c1": 100},
        {"c1": 0},
        0,
        90.0,
        0,
        1,
        None,
        total,
        total_sq,
    )
    # It returned at all, and read every row rather than stopping where the buffer filled.
    assert total[0, 0] == 300.0
