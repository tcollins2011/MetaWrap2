"""Tests for reassemble_bins choosing the original bin when it beats both reassemblies.

reassemble_bins picks the best of three versions of each bin: strict, permissive and original. The
original was once staged as ``.orig.fa`` while CheckM ran with ``-x fasta``, so it was never scored.
On a 26-bin run, 4 bins were replaced by worse reassemblies and the summary reported "0 bins were
not improved by any reassembly".

The fixture (``tests/fixtures/reassembly``, see its SOURCE.md) includes bins where the original
wins. The tests at the end cover the resulting empty-series crash in ``plot_reassembly``.
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

#: The module's default thresholds, which the fixture is built around.
MIN_COMPLETENESS, MAX_CONTAMINATION = 70, 10

#: Expected winners from the fixture's 18 candidate rows. ``bin_006`` fails the thresholds in all
#: three versions, so nothing is chosen for it.
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
    """The original is chosen when it scores higher than both reassemblies.

    ``bin_002``'s permissive reassembly is more complete than the original (93.0% against 91.5%)
    with twice the N50 (52 kb against 26 kb), but is 3.4 points more contaminated. With
    contamination weighted five times, the original scores 585.5 against 570.0.
    """
    assert _chosen() == EXPECTED_WINNERS

    chosen = cbb.choose_best(_versions(), MIN_COMPLETENESS, MAX_CONTAMINATION)
    winner = chosen["bin_002"]
    assert winner.style == "orig"
    for style in ("strict", "permissive"):
        loser = next(v for v in _versions() if v.bin_name == "bin_002" and v.style == style)
        assert winner.beats(loser), style
        # each loser has the higher N50, so ranking on N50 or ignoring the original would pick it
        assert loser.n50 > winner.n50, style


def test_a_bin_whose_reassemblies_both_fail_the_thresholds_keeps_its_original():
    """``bin_005``: strict is 68.0% complete and permissive 14.0% contaminated, so both fail.

    Without the original as a candidate this bin would be dropped from the output.
    """
    chosen = _chosen()
    assert "bin_005.orig" in chosen
    assert [n for n in chosen if n.startswith("bin_005.")] == ["bin_005.orig"]


def test_a_bin_with_no_qualifying_version_is_not_chosen():
    """``bin_006`` is under 70% complete in all three versions, so it is in none of the counts."""
    assert not [n for n in _chosen() if n.startswith("bin_006")]


# --- what gets staged for scoring -------------------------------------------------------


def _bin_set(tmp_path, count=3):
    """Create ``original_bins`` and ``reassembled_bins`` in the layout the CheckM stage expects."""
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
    """Staged originals use BIN_EXTENSION, since CheckM scores the directory with ``-x fasta``.

    Originals staged as ``bin_001.orig.fa`` were not scored.
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
    # Compared against the constant the CheckM command uses, not a literal ".fasta".
    assert all(name.endswith(BIN_EXTENSION) for name in staged), staged


def test_the_name_choose_best_bin_returns_resolves_to_a_staged_file(tmp_path):
    """Every selected name, plus BIN_EXTENSION, is a staged file.

    Otherwise the promote step fails with FileNotFoundError after all reassemblies have finished.
    """
    original, reassembled = _bin_set(tmp_path, count=5)
    reassemble_bins._stage_original_candidates(original, reassembled)
    for name in EXPECTED_WINNERS:
        assert os.path.isfile(os.path.join(reassembled, name + BIN_EXTENSION)), name


# --- the closing summary ----------------------------------------------------------------


def test_the_summary_counts_match_the_selection():
    """Style counts from the chosen names match choose_best_bin's own tally."""
    counts = reassemble_bins._style_counts(EXPECTED_WINNERS)
    assert (counts["orig"], counts["strict"], counts["permissive"]) == (3, 1, 1)
    # the selection function's own tally over the same stats file, as (orig, strict, permissive)
    assert cbb.summarize(CANDIDATES, MIN_COMPLETENESS, MAX_CONTAMINATION) == (3, 1, 1)
    assert sum(counts.values()) == len(EXPECTED_WINNERS)


def test_a_bin_named_for_a_style_is_not_counted_as_that_style():
    """Styles are read from the suffix only, not by substring search in the name."""
    counts = reassemble_bins._style_counts(["orig_strict_permissive_site.orig"])
    assert (counts["orig"], counts["strict"], counts["permissive"]) == (1, 0, 0)


def test_the_summary_reports_the_bins_that_kept_their_original():
    """The summary states how many bins kept their original."""
    message = reassemble_bins._summary_message(EXPECTED_WINNERS)
    assert "Of 5 bins" in message
    assert "1 were improved by the 'strict' reassembly" in message
    assert "1 by the 'permissive' reassembly" in message
    assert "3 kept their original bin" in message
    assert "not improved by any reassembly" not in message


def test_the_summary_total_is_the_bins_chosen_not_the_bins_offered():
    """The total counts chosen bins (5), not every bin in the stats file (6)."""
    assert "Of 5 bins" in reassemble_bins._summary_message(_chosen())


# --- plotting empty series --------------------------------------------------------------


def _run_plot(out, *stats_paths):
    """Run plot_reassembly the way the module does: host interpreter, ``-m``, thresholds first.

    Run as a subprocess because the script runs at import time from ``sys.argv``, and
    reassemble_bins invokes it the same way.
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
    """An empty series is skipped with a note in each panel instead of raising IndexError.

    ``original_bins.stats.tsv`` is header-only when no ``.orig`` candidate passes ``-c``/``-x``.
    Each panel places its label by indexing into the series (``series[x_pos - 1]`` for N50,
    ``series[x_pos]`` for completion and contamination), which fails on an empty list.
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
    """A three-bin series is plotted and labelled.

    ``len // 4 - 1`` is ``-1`` for fewer than four bins, which wrapped the N50 label position to
    the end of the series.
    """
    with open(WINNERS) as fh:
        lines = fh.readlines()
    short = tmp_path / "reassembled_bins.stats.tsv"
    short.write_text("".join(lines[:4]))  # header + 3 bins

    done = _run_plot(tmp_path, short)
    assert done.returncode == 0, done.stderr
    assert "left out of the" not in done.stdout, done.stdout  # not skipped
    assert (tmp_path / "reassembly_results.png").is_file()
