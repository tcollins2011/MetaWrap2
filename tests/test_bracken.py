"""Bracken's table, and the four ways an abundance estimate here stops being one.

Fixtures are real Bracken 3.0.1 output over hand-written inputs; see
``tests/fixtures/bracken/SOURCE.md`` for exactly which file came from where. What is tested:

* **The column layout**, which Bracken's documentation describes in prose without naming, and the
  fact that ``fraction_total_reads`` is a fraction in 0-1 rather than a percentage. Established from
  the one line of ``est_abundance.py`` that writes the header and the one that writes the value;
  three other tools' tables have already cost this project a percentage-for-fraction confusion.
* **A rank mismatch**, because ``--resume`` over a directory whose tables were written at a
  different ``--bracken-rank`` is how a genus table gets read as a species table.
* **A read-length mismatch**, which Bracken reports and then exits 0 on, so it has to be caught
  before the call rather than after it.
* **Redistributing nothing**, which is what a database built from a different Kraken2 database looks
  like, and which is invisible in every column of the output except ``added_reads``.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import bracken as br

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "bracken")


def fixture(name: str) -> str:
    return os.path.join(FIXTURES, name)


# ── the real column layout ───────────────────────────────────────────────────────────────────


def test_column_layout_is_brackens_own():
    """The seven columns, in order, as est_abundance.py writes them."""
    assert br.BRACKEN_COLUMNS == (
        "name",
        "taxonomy_id",
        "taxonomy_lvl",
        "kraken_assigned_reads",
        "added_reads",
        "new_est_reads",
        "fraction_total_reads",
    )

    rows = br.read_table(fixture("sample1_species.bracken"))
    assert len(rows) == 3
    pneumoniae = next(row for row in rows if row.taxonomy_id == 1313)
    assert pneumoniae.name == "Streptococcus pneumoniae"
    assert pneumoniae.taxonomy_lvl == "S"
    assert pneumoniae.kraken_assigned_reads == 150
    assert pneumoniae.added_reads == 151
    assert pneumoniae.new_est_reads == 301
    # added_reads is what Bracken wrote, and it agrees with the arithmetic it claims.
    assert pneumoniae.added_reads == pneumoniae.new_est_reads - pneumoniae.kraken_assigned_reads


def test_the_correction_is_large_which_is_why_the_module_exists():
    """150 of 1,000 reads at its own node, 301 estimated: the species count was off by 2x.

    Kept as an assertion rather than a comment because it is the entire justification for the
    feature. If a change ever makes the parser read `kraken_assigned_reads` into the abundance
    column, every number stays plausible and only this ratio moves.
    """
    rows = br.read_table(fixture("sample1_species.bracken"))
    pneumoniae = next(row for row in rows if row.taxonomy_id == 1313)
    assigned_share = pneumoniae.kraken_assigned_reads / 1000.0
    assert assigned_share == pytest.approx(0.15)
    assert pneumoniae.fraction_total_reads == pytest.approx(0.31684)
    assert pneumoniae.added_share > 0.5


def test_a_taxon_below_the_threshold_is_absent_not_zero():
    """S. mutans has a row in sample 1 and none in sample 2, so the two tables differ in length.

    A consumer that assumes every sample lists the same taxa in the same order mis-aligns them.
    """
    sample1 = {row.name for row in br.read_table(fixture("sample1_species.bracken"))}
    sample2 = {row.name for row in br.read_table(fixture("sample2_species.bracken"))}
    assert "Streptococcus mutans" in sample1
    assert "Streptococcus mutans" not in sample2


# ── fraction, not percentage ──────────────────────────────────────────────────────────────────


def test_fraction_total_reads_is_a_fraction_in_0_to_1():
    """Verified in est_abundance.py: it writes new_est_reads / sum_all_reads at five decimals.

    So the column sums to 1 across one sample's table, not to 100. The 0.99894 here is the rounding
    loss from those five decimals, which is why the parser's tolerance is a band rather than a
    comparison against exactly 1.
    """
    rows = br.read_table(fixture("sample1_species.bracken"))
    total = sum(row.fraction_total_reads for row in rows)
    assert total == pytest.approx(0.99894)
    assert all(0.0 <= row.fraction_total_reads <= 1.0 for row in rows)
    # And the unit travels with the data rather than living in documentation.
    assert br.FRACTION_UNIT == "fraction_total_reads"


def test_a_percentage_column_is_refused_rather_than_scaled():
    """The same table with the last column x100. Nothing may guess a unit it can look up.

    Guessing is what makes this class of bug expensive: a percentage silently treated as a fraction
    makes every abundance 100x too large, and the magnitudes of a typical row make it look fine.
    """
    with pytest.raises(br.BrackenError) as caught:
        br.read_table(fixture("percent_not_fraction.bracken"))
    message = str(caught.value)
    assert "above 1" in message
    assert "fraction in 0-1" in message


def test_the_minus_w_report_is_not_mistaken_for_the_table():
    """Bracken's two outputs look alike and are not. The header is checked, not assumed."""
    with pytest.raises(br.BrackenError) as caught:
        br.read_table(fixture("sample1_bracken_outreport.kreport"))
    assert "-w report" in str(caught.value)


# ── a rank the table was not written at ───────────────────────────────────────────────────────


def test_rank_mismatch_raises_and_says_both_ranks():
    """A genus table read as a species table. The message has to name what it found and wanted."""
    rows = br.read_table(fixture("sample1_genus.bracken"))
    assert [row.taxonomy_lvl for row in rows] == ["G", "G"]

    br.check_table_rank(rows, "G")  # the matching case is silent

    with pytest.raises(br.BrackenError) as caught:
        br.check_table_rank(rows, "S", "sample1_genus.bracken")
    message = str(caught.value)
    assert "estimated at rank G" in message
    assert "species" in message
    assert "sample1_genus.bracken" in message


def test_an_unknown_rank_is_rejected_with_the_list_of_real_ones():
    with pytest.raises(br.BrackenError) as caught:
        br.check_rank("K")
    message = str(caught.value)
    assert "D,P,C,O,F,G,S" in message
    assert "S=species" in message
    # Case and whitespace are normalised rather than rejected; the letter is what matters.
    assert br.check_rank(" s ") == "S"


# ── a read length the database was not built for ──────────────────────────────────────────────


def _database(tmp_path, *read_lengths, kmer_length=None):
    """A stand-in Bracken database: the kmer_distrib file names are all the guard reads."""
    db = tmp_path / "kraken2_db"
    db.mkdir()
    for length in read_lengths:
        (db / (br.KMER_DISTRIB % length)).write_text("mapped_taxid\tgenome_taxids\n")
    if kmer_length is not None:
        (db / "inspect.txt").write_text(
            "# Database options: nucleotide db, k = %d, l = 31\n" % kmer_length
        )
    return str(db)


def test_read_length_mismatch_is_reported_with_what_the_database_has(tmp_path):
    """The lengths the standard Kraken2 database actually ships, and a request for one it does not.

    Reported rather than raised: the caller is a module that wants to say which sample it happened
    on. Bracken prints the same complaint itself and then exits 0, so this must run before the call.
    """
    db = _database(tmp_path, 50, 75, 100, 150, 200, 250, 300)
    assert br.available_read_lengths(db) == [50, 75, 100, 150, 200, 250, 300]
    assert br.check_read_length(db, 150) is None

    problem = br.check_read_length(db, 125)
    assert problem is not None
    assert "no k-mer distribution for 125 bp" in problem
    assert "50 bp, 75 bp, 100 bp, 150 bp" in problem
    assert "database125mers.kmer_distrib" in problem
    assert "bracken-build" in problem


def test_a_database_with_no_distributions_is_a_different_failure(tmp_path):
    """Never built, versus built for other lengths. One is a missing step, the other a wrong flag."""
    db = _database(tmp_path)
    problem = br.check_read_length(db, 100)
    assert problem is not None
    assert "holds no Bracken k-mer distribution" in problem

    missing = br.check_read_length(str(tmp_path / "nope"), 100)
    assert missing is not None
    assert "does not exist" in missing


def test_a_read_shorter_than_the_databases_k_is_rejected_first(tmp_path):
    """A 25 bp read against k=35 yields no complete k-mer, so it cannot be classified at all."""
    db = _database(tmp_path, 25, 100, kmer_length=35)
    assert br.database_kmer_length(db) == 35
    problem = br.check_read_length(db, 25)
    assert problem is not None
    assert "shorter than this Kraken2 database's k of 35" in problem
    assert br.check_read_length(db, 100) is None


# ── Bracken redistributing nothing ────────────────────────────────────────────────────────────


def test_a_mismatched_database_redistributes_nothing_and_is_caught():
    """added_reads is 0 on every row while the fractions still sum to exactly 1.

    This is the quiet failure the guard exists for: the table is complete, the fractions are
    well-formed, and the numbers are just Kraken2's assigned counts renormalised.
    """
    rows = br.read_table(fixture("sample1_mismatched_db.bracken"))
    assert all(row.added_reads == 0 for row in rows)
    assert sum(row.fraction_total_reads for row in rows) == pytest.approx(1.0)

    movable = br.reads_above_rank(fixture("sample1.kreport"), "S")
    assert movable == 450  # 100 at Firmicutes, 200 at Streptococcus, 50 + 100 on the other branch

    summary = br.redistribution(rows, "S", movable)
    assert summary.added == 0
    assert summary.suspicious
    problem = br.check_redistribution(summary)
    assert problem is not None
    assert "Almost nothing was redistributed" in problem
    assert "built from a different Kraken2 database" in problem


def test_a_matched_database_is_not_flagged():
    rows = br.read_table(fixture("sample1_species.bracken"))
    movable = br.reads_above_rank(fixture("sample1.kreport"), "S")
    summary = br.redistribution(rows, "S", movable)
    assert summary.added == 449  # the 450 movable reads, less one lost to integer truncation
    assert summary.moved_share > 0.9
    assert not summary.suspicious
    assert br.check_redistribution(summary) is None
    assert "Bracken added 449" in summary.describe()


def test_nothing_to_move_is_not_a_failure():
    """A sample classified cleanly at the rank has nothing above it, and a 0% share is correct.

    Without this the guard would fire on exactly the samples Bracken was least needed for.
    """
    rows = br.read_table(fixture("sample1_mismatched_db.bracken"))
    summary = br.redistribution(rows, "S", movable=0)
    assert not summary.suspicious
    assert br.check_redistribution(summary) is None


def test_unclassified_reads_are_not_counted_as_movable():
    """The report's 50 unclassified reads are absent from the 450. Bracken cannot place them."""
    assert br.reads_above_rank(fixture("sample1.kreport"), "S") == 450
    # At genus, the genus nodes' own reads are no longer above the rank, so only the phyla remain.
    assert br.reads_above_rank(fixture("sample1.kreport"), "G") == 150


# ── the long-format abundance table ───────────────────────────────────────────────────────────


def test_long_table_has_the_same_shape_as_the_rest_of_the_pipeline(tmp_path):
    """sample, feature, value, unit - the four columns `catalogue` writes and `diff_abundance` reads.

    Checked as a header equality rather than a description, because the whole value of matching is
    that no conversion step stands between this module and the statistics one.
    """
    assert br.LONG_COLUMNS == ("sample", "feature", "value", "unit")

    per_sample = {
        "sample1": br.read_table(fixture("sample1_species.bracken")),
        "sample2": br.read_table(fixture("sample2_species.bracken")),
    }
    out = tmp_path / "bracken_abundance.tsv"
    written = br.write_long_table(str(out), per_sample)
    assert written == 5  # three taxa in sample1, two in sample2

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(br.LONG_COLUMNS)
    assert all(len(line.split("\t")) == 4 for line in lines)

    rows = [line.split("\t") for line in lines[1:]]
    assert {row[0] for row in rows} == {"sample1", "sample2"}
    assert {row[3] for row in rows} == {"fraction_total_reads"}
    # The feature is the taxon name, which is what diff_abundance reports back to the user.
    assert ("sample1", "Escherichia coli") in {(row[0], row[1]) for row in rows}
    # Values are the fractions, unscaled, and still fractions on the way out.
    assert all(0.0 < float(row[2]) <= 1.0 for row in rows)


def test_long_table_is_byte_stable_across_runs(tmp_path):
    """A table that reorders itself cannot be diffed, which is how a changed number gets missed."""
    per_sample = {"s1": br.read_table(fixture("sample1_species.bracken"))}
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    br.write_long_table(str(first), per_sample)
    br.write_long_table(str(second), dict(per_sample))
    assert first.read_bytes() == second.read_bytes()
    # Descending by estimate, so the most abundant taxon is the first row a reader sees.
    names = [line.split("\t")[1] for line in first.read_text().splitlines()[1:]]
    assert names == ["Escherichia coli", "Streptococcus pneumoniae", "Streptococcus mutans"]


def test_the_read_count_unit_is_available_and_labelled():
    """Counts rather than fractions, for a caller that wants them - but never silently.

    The unit is on every row either way; that is the point of having two named units instead of a
    boolean, and of the writer refusing a unit it does not produce.
    """
    rows = br.read_table(fixture("sample1_species.bracken"))
    counts = br.long_rows("s1", rows, br.READS_UNIT)
    assert [row[3] for row in counts] == ["est_reads"] * 3
    assert {row[2] for row in counts} == {"301", "248", "400"}

    with pytest.raises(br.BrackenError) as caught:
        br.long_rows("s1", rows, "percent")
    assert "not one this writer produces" in str(caught.value)


# ── the module's command and flags ────────────────────────────────────────────────────────────


def test_bracken_command_template_formats_cleanly():
    from metawrap2.modules import kraken2

    cmd = kraken2.BRACKEN.format(
        db="/db",
        report="s.kreport",
        out="s.bracken",
        outreport="s_bracken_species.kreport",
        read_length=150,
        rank="S",
        threshold=10,
    )
    assert "{" not in cmd
    assert cmd.startswith("bracken -d /db")
    # -i takes the aggregated report, never the per-read --output file; Bracken exits 1 on the latter.
    assert "-i s.kreport" in cmd
    assert "-r 150" in cmd and "-l S" in cmd and "-t 10" in cmd

    # The report flag is appended to kraken2 only under --bracken, so a run without the flag is
    # byte-identical to one from before the flag existed.
    assert kraken2.KRAKEN2_REPORT.format(report="s.kreport") == "--report s.kreport"
    assert "--report" not in kraken2.KRAKEN2_PAIRED
    assert "--report" not in kraken2.KRAKEN2_SINGLE


def test_help_renders_with_its_percent_signs_escaped():
    """A bare `10%` in an argparse help string makes every --help in the package fail.

    Rendering the help is the only test that catches it, because argparse interpolates help strings
    lazily - the module imports and runs fine right up until someone asks for --help.
    """
    import contextlib
    import io

    from metawrap2.modules import kraken2

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer), pytest.raises(SystemExit):
        kraken2._parse_args(["--help"])
    help_text = buffer.getvalue()

    assert "--bracken-rank" in help_text
    assert "{D,P,C,O,F,G,S}" in help_text
    # The escaped %% renders as one literal %, and no %% survives into the output.
    assert "10% of" in help_text
    assert "%%" not in help_text
    # argparse's own %(default)s interpolation still works alongside the escaping.
    assert "default: 100" in help_text
