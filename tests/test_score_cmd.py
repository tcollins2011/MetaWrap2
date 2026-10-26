"""Tests for ``metawrap2 score``.

The tests check that the table matches the one ``bin_refinement`` writes, column for column, and
that every bin set in one run is scored by the same tool in the same env. ``metawrap2.checkm.run``
and ``metawrap2.gunc.run`` are replaced with stand-ins that record the command and conda env and
write the report the real tool would write.
"""

from __future__ import annotations

import json
import os
import shlex
from typing import Dict, List, Tuple

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2 import config as config_mod
from metawrap2 import gunc as gunc_mod
from metawrap2.commands import score_cmd
from metawrap2.constants import BIN_EXTENSION
from metawrap2.modules import _common, bin_refinement

#: Stand-in scores per genome: completeness, contamination, GUNC clade separation score, GUNC
#: pass. One genome each: clean near-complete, chimeric near-complete, between 70/10 and 90/5, and
#: below both.
SCORES: Dict[str, Tuple[float, float, float, bool]] = {
    "bin_001": (97.5, 1.2, 0.10, True),
    "bin_002": (95.0, 2.0, 0.80, False),
    "bin_003": (78.0, 6.0, 0.05, True),
    "bin_004": (31.0, 1.0, 0.02, True),
}

#: Sequence per genome, above :data:`metawrap2.checkm.MIN_PLAUSIBLE_BIN_BP` (200 kb) so the
#: default size pre-filter passes it to the stand-in scorer.
GENOME_BP = "ACGT" * 51_000

#: GUNC's report name, with the database version in the middle as in the real one (the finder
#: matches on prefix and suffix).
GUNC_REPORT = gunc_mod.REPORT_PREFIX + "progenomes_2.1" + gunc_mod.REPORT_SUFFIX


# ── fixtures on disk ─────────────────────────────────────────────────────────────────────────


def _bin_set(path, extension=BIN_EXTENSION, genomes=tuple(SCORES), tag="c"):
    """A directory of genomes named for :data:`SCORES`, under *extension*.

    *tag* goes into the contig names, so two bin sets can hold genomes with the same filenames and
    different contents (and different CheckM cache keys).
    """
    path.mkdir(parents=True, exist_ok=True)
    for name in genomes:
        (path / (name + extension)).write_text(">%s_%s\n%s\n" % (tag, name, GENOME_BP))
    return str(path)


def _config(tmp_path, use_conda_envs=True):
    """A config with a GUNC database that exists, since GUNC runs unless told not to."""
    db = tmp_path / "gunc_db.dmnd"
    db.write_text("not a real diamond database\n")
    path = tmp_path / "metawrap2.toml"
    path.write_text(
        "[settings]\nuse_conda_envs = %s\n\n[databases]\nGUNC_DB = %r\n"
        % ("true" if use_conda_envs else "false", str(db))
    )
    return str(path)


# ── stand-ins for CheckM1, CheckM2 and GUNC ──────────────────────────────────────────────────


def _staged(bins_dir: str) -> List[str]:
    """The genome names a scorer was pointed at: the stems of the links `_stage` left."""
    return sorted(
        os.path.splitext(name)[0] for name in os.listdir(bins_dir) if name.endswith(BIN_EXTENSION)
    )


def _write_checkm2_report(out: str, names: List[str]) -> None:
    os.makedirs(out, exist_ok=True)
    columns = (
        "Name",
        "Completeness",
        "Contamination",
        "Genome_Size",
        "GC_Content",
        "Contig_N50",
    )
    with open(os.path.join(out, "quality_report.tsv"), "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for name in names:
            completeness, contamination, _css, _passed = SCORES[name]
            fh.write(
                "%s\t%s\t%s\t%d\t0.51\t%d\n"
                % (name, completeness, contamination, len(GENOME_BP), 4000)
            )


def _write_checkm1_report(out: str, names: List[str]) -> None:
    storage = os.path.join(out, "storage")
    os.makedirs(storage, exist_ok=True)
    with open(os.path.join(storage, "bin_stats_ext.tsv"), "w") as fh:
        for name in names:
            completeness, contamination, _css, _passed = SCORES[name]
            fh.write(
                "%s\t%r\n"
                % (
                    name,
                    {
                        "marker lineage": "k__Bacteria",
                        "Completeness": completeness,
                        "Contamination": contamination,
                        "GC": 0.51,
                        "N50 (contigs)": 4000,
                        "Genome size": len(GENOME_BP),
                    },
                )
            )


def _write_gunc_report(out: str, names: List[str]) -> None:
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, GUNC_REPORT), "w") as fh:
        fh.write("genome\ttaxonomic_level\tclade_separation_score\tpass.GUNC\n")
        for name in names:
            _completeness, _contamination, css, passed = SCORES[name]
            fh.write("%s\tgenus\t%s\t%s\n" % (name, css, passed))


class Calls(list):
    """Every scorer invocation, as ``(command, tool, conda env)``."""

    def tools(self) -> List[str]:
        return [tool for _cmd, tool, _env in self]

    def env_of(self, tool: str):
        return next(env for _cmd, name, env in self if name == tool)

    def command_of(self, tool: str) -> str:
        return next(cmd for cmd, name, _env in self if name == tool)


def _fake_scorers(monkeypatch) -> Calls:
    """Replace the ``run`` that CheckM and GUNC call, and let their own parsers do the rest.

    Only ``run`` is replaced, so the commands, conda envs, report finding and all three
    annotation passes are the real ones.
    """
    calls = Calls()

    def record(cmd, tool, env):
        calls.append((cmd if isinstance(cmd, str) else " ".join(cmd), tool, env))
        return shlex.split(cmd) if isinstance(cmd, str) else list(cmd)

    def fake_checkm(cmd, env=None, tool=None, **kwargs):
        parts = record(cmd, tool, env)
        if tool == "checkm2":
            bins = parts[parts.index("--input") + 1]
            _write_checkm2_report(parts[parts.index("--output-directory") + 1], _staged(bins))
        elif tool == "checkm":
            # `checkm lineage_wf -x <ext> <bins> <out> ...`
            at = parts.index("lineage_wf")
            _write_checkm1_report(parts[at + 4], _staged(parts[at + 3]))
        return 0

    def fake_gunc(cmd, env=None, tool=None, **kwargs):
        parts = record(cmd, tool, env)
        bins = parts[parts.index("--input_dir") + 1]
        _write_gunc_report(parts[parts.index("--out_dir") + 1], _staged(bins))
        return 0

    monkeypatch.setattr(checkm_mod, "run", fake_checkm)
    monkeypatch.setattr(gunc_mod, "run", fake_gunc)
    # The env names are asserted on, so the env manager and env existence checks are stubbed.
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)
    monkeypatch.setattr(config_mod, "conda_env_exists", lambda env: True)
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: True)
    return calls


def _score(tmp_path, monkeypatch, *bin_dirs, extra=(), out="SCORED", use_conda_envs=True):
    """One ``metawrap2 score`` run. Returns the output directory and the recorded calls."""
    calls = _fake_scorers(monkeypatch)
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    output = str(tmp_path / out)
    argv = ["-o", output, "-t", "1", "-m", "8", "--config", _config(tmp_path, use_conda_envs)]
    argv += list(extra) + list(bin_dirs)
    assert score_cmd.main(argv) == 0
    return output, calls


def _header(path: str) -> List[str]:
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")


def _rows(path: str) -> List[List[str]]:
    with open(path) as fh:
        return [line.rstrip("\n").split("\t") for line in fh.read().splitlines()[1:] if line]


# ── the schema ───────────────────────────────────────────────────────────────────────────────


def test_the_table_has_the_columns_bin_refinement_writes(tmp_path, monkeypatch):
    """The header equals the one a ``bin_refinement`` run writes for the same bins."""
    bins = _bin_set(tmp_path / "metawrap_70_10_bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)

    # bin_refinement with one bin set runs only CheckM -> GUNC -> MIMAG.
    monkeypatch.setattr(bin_refinement, "run", lambda *a, **kw: 0)  # the ranking plot
    refined = str(tmp_path / "REFINED")
    assert (
        bin_refinement.main(
            [
                "-o",
                refined,
                "-A",
                bins,
                "-t",
                "1",
                "-m",
                "8",
                "--checkm2",
                "--config",
                _config(tmp_path),
            ]
        )
        == 0
    )

    assert _header(os.path.join(scored, "metawrap_70_10_bins.stats.tsv")) == _header(
        os.path.join(refined, "metawrap_70_10_bins.stats.tsv")
    )


def test_the_schema_is_composed_from_the_modules_that_write_it(tmp_path, monkeypatch):
    """``stats_columns`` matches the header the CheckM, GUNC and MIMAG passes write."""
    bins = _bin_set(tmp_path / "bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    header = _header(os.path.join(scored, "bins.stats.tsv"))
    assert header == score_cmd.stats_columns()
    assert header[: len(checkm_mod.STATS_HEADER)] == list(checkm_mod.STATS_HEADER)
    assert header[-len(score_cmd._mimag.STATS_COLUMNS) :] == list(score_cmd._mimag.STATS_COLUMNS)

    # The columns gunc.annotate_stats appends to a table.
    table = tmp_path / "plain.stats.tsv"
    table.write_text("\t".join(checkm_mod.STATS_HEADER) + "\nbin_001\t90\t1\t0.5\tBacteria\t1\t1\n")
    gunc_mod.annotate_stats(str(table), {})
    appended = _header(str(table))[len(checkm_mod.STATS_HEADER) :]
    assert tuple(appended) == score_cmd.GUNC_STATS_COLUMNS


def test_every_genome_gets_a_row_with_its_scores_and_its_mimag_tag(tmp_path, monkeypatch):
    """Each genome's row holds its CheckM scores, GUNC result and MIMAG tag."""
    bins = _bin_set(tmp_path / "bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    header = _header(os.path.join(scored, "bins.stats.tsv"))
    rows = {row[0]: dict(zip(header, row)) for row in _rows(os.path.join(scored, "bins.stats.tsv"))}
    assert set(rows) == set(SCORES)
    assert rows["bin_001"]["completeness"] == "97.5"
    assert rows["bin_001"]["gunc_pass"] == "yes"
    assert rows["bin_002"]["gunc_pass"] == "no"
    # No rRNA evidence, so a genome over 90/5 is tagged high quality, unverified.
    assert rows["bin_001"]["mimag_tag"] == score_cmd._mimag.MimagTag.HIGH_UNVERIFIED.value
    assert rows["bin_001"]["mimag_rrna"] == score_cmd._mimag.UNKNOWN


# ── which scorer, and where it runs ──────────────────────────────────────────────────────────


def test_checkm2_scores_by_default_in_its_own_env(tmp_path, monkeypatch):
    """CheckM2 is the default scorer and runs in the CheckM2 env; GUNC in the module env."""
    bins = _bin_set(tmp_path / "bins")
    _scored, calls = _score(tmp_path, monkeypatch, bins)
    assert "checkm2" in calls.tools() and "checkm" not in calls.tools()
    assert calls.env_of("checkm2") == config_mod.OPTIONAL_ENVS["checkm2"]
    assert calls.env_of("gunc") == config_mod.MODULE_ENVS[score_cmd.CHECKM_MODULE]


def test_checkm1_selects_checkm1(tmp_path, monkeypatch):
    """``--checkm1`` runs ``checkm lineage_wf`` in bin_refinement's env and no CheckM2."""
    bins = _bin_set(tmp_path / "bins")
    scored, calls = _score(tmp_path, monkeypatch, bins, extra=["--checkm1"])
    assert "checkm" in calls.tools() and "checkm2" not in calls.tools()
    assert calls.env_of("checkm") == config_mod.MODULE_ENVS[score_cmd.CHECKM_MODULE]
    # CheckM1 reports a marker lineage; CheckM2 writes "checkm2" in that column.
    lineage = _header(os.path.join(scored, "bins.stats.tsv")).index("lineage")
    assert {row[lineage] for row in _rows(os.path.join(scored, "bins.stats.tsv"))} == {"Bacteria"}


def test_the_two_estimators_cannot_share_a_score_cache(tmp_path, monkeypatch):
    """The score cache file is named for the estimator that filled it."""
    bins = _bin_set(tmp_path / "bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    assert os.path.isfile(os.path.join(scored, score_cmd.SCORE_CACHE % "checkm2"))
    assert not os.path.exists(os.path.join(scored, score_cmd.SCORE_CACHE % "checkm1"))


# ── GUNC ─────────────────────────────────────────────────────────────────────────────────────


def test_gunc_runs_by_default_over_the_staged_genomes(tmp_path, monkeypatch):
    bins = _bin_set(tmp_path / "bins")
    scored, calls = _score(tmp_path, monkeypatch, bins)
    command = calls.command_of("gunc")
    assert "--input_dir %s" % os.path.join(scored, "bins") in command
    # GUNC names genomes with basename.split(file_suffix)[0], so the suffix needs its dot.
    assert "--file_suffix %s" % BIN_EXTENSION in command
    assert gunc_mod.find_report(os.path.join(scored, "bins.gunc"))


def test_no_gunc_omits_gunc_and_its_columns(tmp_path, monkeypatch):
    """``--no-gunc`` runs no GUNC and leaves the ``gunc_*`` columns out of the table."""
    bins = _bin_set(tmp_path / "bins")
    scored, calls = _score(tmp_path, monkeypatch, bins, extra=["--no-gunc"])
    assert "gunc" not in calls.tools()
    header = _header(os.path.join(scored, "bins.stats.tsv"))
    assert header == score_cmd.stats_columns(with_gunc=False)
    assert not [name for name in header if name.startswith("gunc")]


def test_a_missing_gunc_database_stops_the_run_before_checkm(tmp_path, monkeypatch):
    """An unconfigured GUNC database stops the run before any scorer runs."""
    bins = _bin_set(tmp_path / "bins")
    calls = _fake_scorers(monkeypatch)
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    config = tmp_path / "no_db.toml"
    config.write_text("[settings]\nuse_conda_envs = true\n")
    with pytest.raises(SystemExit):
        score_cmd.main(["-o", str(tmp_path / "out"), "-t", "1", "--config", str(config), bins])
    assert calls.tools() == []


# ── the extension ────────────────────────────────────────────────────────────────────────────


def test_a_non_default_extension_is_honoured(tmp_path, monkeypatch):
    """``-x .fa`` scores every genome in a metaWRAP 1 bin set."""
    bins = _bin_set(tmp_path / "metawrap1_bins", extension=".fa")
    scored, calls = _score(tmp_path, monkeypatch, bins, extra=["-x", ".fa"])
    rows = _rows(os.path.join(scored, "metawrap1_bins.stats.tsv"))
    assert sorted(row[0] for row in rows) == sorted(SCORES)
    # The genomes are staged under the default extension, which the scorers are given.
    assert "-x %s" % BIN_EXTENSION.lstrip(".") in calls.command_of("checkm2")
    assert sorted(os.listdir(os.path.join(scored, "metawrap1_bins"))) == sorted(
        name + BIN_EXTENSION for name in SCORES
    )


def test_the_default_extension_finds_nothing_in_a_metawrap1_bin_set(tmp_path, monkeypatch):
    """Without ``-x``, a ``.fa`` bin set is scored as empty and nothing runs."""
    bins = _bin_set(tmp_path / "metawrap1_bins", extension=".fa")
    scored, calls = _score(tmp_path, monkeypatch, bins)
    assert calls.tools() == []
    assert _rows(os.path.join(scored, "metawrap1_bins.stats.tsv")) == []


def test_an_extension_without_its_dot_is_the_same_answer():
    """Extensions with and without the leading dot normalise to the dotted form."""
    assert score_cmd.normalise_extension("fa") == ".fa"
    assert score_cmd.normalise_extension(".fa") == ".fa"
    assert score_cmd.normalise_extension(" fasta ") == ".fasta"
    assert score_cmd.normalise_extension("") == BIN_EXTENSION


# ── empty, odd and colliding inputs ──────────────────────────────────────────────────────────


def test_an_empty_directory_is_reported_not_fatal(tmp_path, monkeypatch, capsys):
    """An empty bin set gets a header-only table and a warning; the other sets are scored."""
    empty = tmp_path / "empty_bins"
    empty.mkdir()
    full = _bin_set(tmp_path / "full_bins")
    scored, _calls = _score(tmp_path, monkeypatch, str(empty), full)

    table = os.path.join(scored, "empty_bins.stats.tsv")
    assert _header(table) == score_cmd.stats_columns()
    assert _rows(table) == []
    assert len(_rows(os.path.join(scored, "full_bins.stats.tsv"))) == len(SCORES)
    printed = capsys.readouterr().out
    assert "nothing to score" in printed
    # The warning mentions -x.
    assert "-x" in printed


def test_the_unbinned_catch_all_is_not_scored_as_a_genome(tmp_path, monkeypatch):
    """``unbinned.fasta`` gets no row in the table."""
    bins = _bin_set(tmp_path / "bins")
    (tmp_path / "bins" / "unbinned.fasta").write_text(">leftovers\n%s\n" % GENOME_BP)
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    assert sorted(row[0] for row in _rows(os.path.join(scored, "bins.stats.tsv"))) == sorted(SCORES)


def test_input_directories_that_share_a_basename_stay_distinct(tmp_path):
    """Inputs with the same basename get distinct names from their parent or a counter."""
    names = score_cmd.set_names(
        [
            "/runs/metawrap1/metawrap_70_10_bins",
            "/runs/metawrap2/metawrap_70_10_bins",
            "/runs/published/bins",
        ]
    )
    assert names == [
        "metawrap1.metawrap_70_10_bins",
        "metawrap2.metawrap_70_10_bins",
        "bins",
    ]
    assert len(set(names)) == 3
    # Same parent: a counter separates them. A trailing slash does not change the parent.
    assert score_cmd.set_names(["/runs/a/bins", "/runs/a/bins/"]) == ["a.bins", "a.bins.2"]


def test_a_missing_directory_is_reported_before_anything_runs(tmp_path, monkeypatch):
    calls = _fake_scorers(monkeypatch)
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    with pytest.raises(SystemExit):
        score_cmd.main(
            [
                "-o",
                str(tmp_path / "out"),
                "--config",
                _config(tmp_path),
                str(tmp_path / "not_there"),
            ]
        )
    assert calls.tools() == []


# ── the comparison table ─────────────────────────────────────────────────────────────────────


def test_several_directories_produce_one_comparison_table(tmp_path, monkeypatch, capsys):
    """Two bin sets give one comparison table, in command-line order, from one scorer."""
    one = _bin_set(tmp_path / "metawrap1_bins", tag="one")
    two = _bin_set(tmp_path / "metawrap2_bins", tag="two")
    _scored, calls = _score(tmp_path, monkeypatch, one, two)
    printed = capsys.readouterr().out

    table = printed[printed.index("Scored with CheckM2") :]
    assert "metawrap1_bins" in table and "metawrap2_bins" in table
    assert table.index("metawrap1_bins") < table.index("metawrap2_bins")  # command-line order
    # Both scored by CheckM2 in the same env.
    assert calls.tools().count("checkm2") == 2
    assert len({env for _cmd, tool, env in calls if tool == "checkm2"}) == 1
    assert score_cmd.COMPARISON_NOTE.split("\n")[0][:40] in printed


def test_one_directory_prints_no_comparison_table(tmp_path, monkeypatch, capsys):
    """A single bin set prints no comparison table."""
    bins = _bin_set(tmp_path / "bins")
    _score(tmp_path, monkeypatch, bins)
    assert "Scored with" not in capsys.readouterr().out


def test_clean_near_complete_is_counted_per_genome(tmp_path, monkeypatch):
    """Clean near-complete is counted per genome: 1 here, where a pass fraction gives 1.5.

    2 genomes clear 90/5 and 3 of 4 pass GUNC; the chimeric genome is one of the near-complete two.
    """
    bins = _bin_set(tmp_path / "bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    result = score_cmd.summarise_table(
        "bins",
        bins,
        os.path.join(scored, "bins.stats.tsv"),
        completeness=70.0,
        contamination=10.0,
        with_gunc=True,
    )
    assert result.genomes == 4
    assert result.passing == 3
    assert result.tiers[score_cmd.NEAR_COMPLETE_TIER] == 2
    assert (result.gunc_pass, result.gunc_scored) == (3, 4)
    assert result.clean_near_complete == 1


def test_the_thresholds_the_counts_use_are_the_ones_asked_for(tmp_path, monkeypatch, capsys):
    """``-c``/``--contamination`` set the passing count and its heading."""
    bins = _bin_set(tmp_path / "bins")
    _score(tmp_path, monkeypatch, bins, extra=["-c", "90", "--contamination", "5"])
    printed = capsys.readouterr().out
    assert ">=90/<=5" in printed
    passing = [line for line in printed.splitlines() if ">=90/<=5" in line]
    assert passing and passing[0].split()[-1] == "2"


# ── the run's own record of itself ───────────────────────────────────────────────────────────


def test_the_run_records_its_provenance_and_its_cost(tmp_path, monkeypatch, capsys):
    """The run writes a log, a ``run_config.json`` provenance record and a cost line."""
    bins = _bin_set(tmp_path / "bins")
    scored, _calls = _score(tmp_path, monkeypatch, bins)
    assert os.path.isfile(os.path.join(scored, "metawrap2.log"))
    with open(os.path.join(scored, "run_config.json")) as fh:
        record = json.load(fh)
    assert record["module"] == "score"
    # The scorer, env and database are recorded.
    assert record["parameters"]["checkm1"] is False
    assert record["conda_env"] == config_mod.MODULE_ENVS[score_cmd.CHECKM_MODULE]
    assert record["databases"]["GUNC_DB"].endswith("gunc_db.dmnd")
    assert record["inputs"] == [bins]
    cost = [
        line for line in capsys.readouterr().out.splitlines() if "wall" in line and "score" in line
    ]
    assert cost
