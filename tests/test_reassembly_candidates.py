"""The third candidate: reassemble_bins choosing the *original* bin over both reassemblies.

The module exists to pick the best of three versions of every bin - the strict reassembly, the
permissive reassembly, and the bin it started from. The original was copied in for scoring under a
hardcoded ``.orig.fa`` while CheckM was invoked ``-x fasta``, so it was never scored, never a
candidate, and never chosen: on a 26-bin run, 4 bins (15%) were replaced with a reassembly that
scored *worse* than the original, and the module announced "0 bins were not improved by any
reassembly" - the defect restating itself as a result.

So the fixture these tests are built on (``tests/fixtures/reassembly``, see its SOURCE.md) is one
where the original wins, which is the case that failed silently. The downstream crash it caused -
an empty series in ``plot_reassembly``, which is where a user actually saw this - is pinned at the
bottom of the file.
"""

import os
import subprocess
import sys

import pytest

from metawrap2.constants import BIN_EXTENSION
from metawrap2.modules import reassemble_bins
from metawrap2.scripts import choose_best_bin as cbb

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "reassembly")
CANDIDATES = os.path.join(FIXTURES, "candidates.stats.tsv")
WINNERS = os.path.join(FIXTURES, "reassembled_bins.stats.tsv")
HEADER_ONLY = os.path.join(FIXTURES, "original_bins.stats.tsv")

#: The module's own defaults, which is what the fixture's margins are drawn around.
MIN_COMPLETENESS, MAX_CONTAMINATION = 70, 10

#: What the fixture's 18 candidate rows come to. ``bin_006`` is in the stats file but fails the
#: thresholds in all three versions, so nothing is chosen for it.
EXPECTED_WINNERS = [
    "bin_001.strict",
    "bin_002.orig",
    "bin_003.permissive",
    "bin_004.orig",
    "bin_005.orig",
]


def _versions(path=CANDIDATES):
    with open(path) as fh:
        return cbb.parse_stats(fh)


def _chosen():
    return cbb.best_bin_names(CANDIDATES, MIN_COMPLETENESS, MAX_CONTAMINATION)


# --- the original as a candidate --------------------------------------------------------


def test_the_original_bin_wins_when_the_reassemblies_score_worse():
    """The test the reporter asked for: "a fixture where the original bin is the best of the three".

    ``bin_002``'s permissive reassembly is *more* complete than the original (93.0% against 91.5%)
    and has twice its N50 (52 kb against 26 kb), and is still the wrong answer, because it is 3.4
    points more contaminated and contamination is penalised five times completeness: 585.5 against
    570.0. Choosing it is the measurable cost of the defect - a bin handed on as improved when the
    input was better.
    """
    assert _chosen() == EXPECTED_WINNERS

    chosen = cbb.choose_best(_versions(), MIN_COMPLETENESS, MAX_CONTAMINATION)
    winner = chosen["bin_002"]
    assert winner.style == "orig"
    for style in ("strict", "permissive"):
        loser = next(v for v in _versions() if v.bin_name == "bin_002" and v.style == style)
        assert winner.beats(loser), style
        # and each loser has the longer contigs, so a selection that ranked on N50 - or one that
        # could not see the original at all - would pick one of them
        assert loser.n50 > winner.n50, style


def test_a_bin_whose_reassemblies_both_fail_the_thresholds_keeps_its_original():
    """``bin_005``: strict drops to 68.0% complete and permissive rises to 14.0% contaminated.

    Neither is a candidate at all, so with the original excluded from scoring this bin had no
    qualifying version and dropped out of the output entirely - a bin lost, not merely degraded.
    """
    chosen = _chosen()
    assert "bin_005.orig" in chosen
    assert [n for n in chosen if n.startswith("bin_005.")] == ["bin_005.orig"]


def test_a_bin_with_no_qualifying_version_is_not_chosen():
    """``bin_006`` is under 70% complete in all three versions, so it is in none of the counts."""
    assert not [n for n in _chosen() if n.startswith("bin_006")]


# --- what gets staged for scoring -------------------------------------------------------


def _bin_set(tmp_path, count=3):
    """An ``original_bins``/``reassembled_bins`` pair in the layout the CheckM stage expects."""
    original = tmp_path / "original_bins"
    reassembled = tmp_path / "reassembled_bins"
    original.mkdir()
    reassembled.mkdir()
    for number in range(1, count + 1):
        stem = "bin_%03d" % number
        (original / (stem + BIN_EXTENSION)).write_text(">c1\nACGTACGTAC\n")
        for style in ("strict", "permissive"):
            (reassembled / ("%s.%s%s" % (stem, style, BIN_EXTENSION))).write_text(">c1\nACGTA\n")
    return str(original), str(reassembled)


def test_every_candidate_staged_for_scoring_carries_the_bin_extension(tmp_path):
    """One CheckM run scores the whole directory with ``-x fasta``, so a ``.fa`` copy is invisible.

    This is the root cause. The staged originals used to be ``bin_001.orig.fa``, which left the
    reporter's directory holding 26 ``.fa`` files alongside 52 ``.fasta`` ones, and only the 52 in
    the stats table that the choice was made from.
    """
    original, reassembled = _bin_set(tmp_path)
    assert reassemble_bins._stage_original_candidates(original, reassembled) == 3

    staged = sorted(os.listdir(reassembled))
    assert len(staged) == 9  # 3 bins x (strict, permissive, orig)
    assert [n for n in staged if ".orig" in n] == [
        "bin_001.orig.fasta",
        "bin_002.orig.fasta",
        "bin_003.orig.fasta",
    ]
    # Compared against the constant the CheckM command is built from rather than against a spelled
    # out ".fasta", so renaming the extension again cannot reintroduce this.
    assert all(name.endswith(BIN_EXTENSION) for name in staged), staged


def test_the_name_choose_best_bin_returns_resolves_to_a_staged_file(tmp_path):
    """The two halves have to agree: the winner is copied by ``name + BIN_EXTENSION``.

    A staged name the selection cannot produce, or a selected name that is not on disk, both end as
    the same FileNotFoundError at the promote step - after every reassembly is done.
    """
    original, reassembled = _bin_set(tmp_path, count=5)
    reassemble_bins._stage_original_candidates(original, reassembled)
    for name in EXPECTED_WINNERS:
        assert os.path.isfile(os.path.join(reassembled, name + BIN_EXTENSION)), name


# --- the closing summary ----------------------------------------------------------------


def test_the_summary_counts_match_the_selection():
    """Counted off the chosen names' style suffixes, and cross-checked against choose_best_bin."""
    counts = reassemble_bins._style_counts(EXPECTED_WINNERS)
    assert (counts["orig"], counts["strict"], counts["permissive"]) == (3, 1, 1)
    # the selection function's own tally over the same stats file, as (orig, strict, permissive)
    assert cbb.summarize(CANDIDATES, MIN_COMPLETENESS, MAX_CONTAMINATION) == (3, 1, 1)
    assert sum(counts.values()) == len(EXPECTED_WINNERS)


def test_a_bin_named_for_a_style_is_not_counted_as_that_style():
    """The counts were substring tests over filenames, so one bin could be counted three times."""
    counts = reassemble_bins._style_counts(["orig_strict_permissive_site.orig"])
    assert (counts["orig"], counts["strict"], counts["permissive"]) == (1, 0, 0)


def test_the_summary_reports_the_bins_that_kept_their_original():
    """It used to say "0 bins were not improved by any reassembly", which was the bug speaking."""
    message = reassemble_bins._summary_message(EXPECTED_WINNERS)
    assert "Of 5 bins" in message
    assert "1 were improved by the 'strict' reassembly" in message
    assert "1 by the 'permissive' reassembly" in message
    assert "3 kept their original bin" in message
    assert "not improved by any reassembly" not in message


def test_the_summary_total_is_the_bins_chosen_not_the_bins_offered():
    """``bin_006`` qualified in no version, so a total of 6 would claim a bin that is not there."""
    assert "Of 5 bins" in reassemble_bins._summary_message(_chosen())


# --- the empty series that crashed the run ----------------------------------------------


def _run_plot(out, *stats_paths):
    """Run plot_reassembly the way the module does: host interpreter, ``-m``, thresholds first.

    A subprocess because the script is module-level code driven by ``sys.argv`` - importing it
    would plot - and because that is exactly how reassemble_bins invokes it.
    """
    src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
    env = dict(os.environ, PYTHONPATH=src, MPLBACKEND="agg")
    return subprocess.run(
        [sys.executable, "-m", "metawrap2.scripts.plot_reassembly", str(out)]
        + [str(MIN_COMPLETENESS), str(MAX_CONTAMINATION)]
        + [str(p) for p in stats_paths],
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,  # the exit status is the thing under test
    )


@pytest.mark.parametrize(
    "stats",
    [
        pytest.param((WINNERS, HEADER_ONLY), id="one_empty_series"),
        pytest.param((HEADER_ONLY,), id="every_series_empty"),
    ],
)
def test_an_empty_series_does_not_crash_any_of_the_three_plots(tmp_path, stats):
    """The visible symptom: IndexError in all three panels, 62m53s into an otherwise finished run.

    ``original_bins.stats.tsv`` is header-only whenever no ``.orig`` candidate was scored, and each
    panel wrote its label at a fraction of the series' length - ``series[x_pos - 1]``, which is
    ``series[-1]`` on an empty list, for the N50 panel, and ``series[x_pos]`` for the completion and
    contamination ones. A header-only table is still legitimate once the extension is fixed (every
    original can fail ``-c``/``-x``), so the guard has to hold rather than the file being prevented.
    """
    done = _run_plot(tmp_path, *stats)
    assert done.returncode == 0, done.stderr
    assert "IndexError" not in done.stderr, done.stderr
    # skipped with a note, once per panel, naming the series that had nothing in it
    for panel in ("N50", "completion", "contamination"):
        assert "left out of the %s panel" % panel in done.stdout, done.stdout
    assert done.stdout.count("No bin in original_bins") == 3, done.stdout
    for name in ("reassembly_results.png", "reassembly_results.eps"):
        assert (tmp_path / name).is_file(), name


def test_a_series_of_three_bins_is_labelled_without_wrapping(tmp_path):
    """``len // 4 - 1`` is ``-1`` for a series of fewer than four bins.

    Not fatal - it wrapped to the end of the series - but the N50 panel then drew the label at x=0
    carrying the *lowest* N50 in the set, which is the same off-by-one as the crash. Three bins is
    the smallest realistic bin set, and the only way to reach the label code is to run the script.
    """
    with open(WINNERS) as fh:
        lines = fh.readlines()
    short = tmp_path / "reassembled_bins.stats.tsv"
    short.write_text("".join(lines[:4]))  # header + 3 bins

    done = _run_plot(tmp_path, short)
    assert done.returncode == 0, done.stderr
    assert "left out of the" not in done.stdout, done.stdout  # it has points; it is not skipped
    assert (tmp_path / "reassembly_results.png").is_file()
