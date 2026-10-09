"""GUNC: report parsing, the stats annotation, and the wiring into bin_refinement."""

import os

import pytest

from metawrap2 import gunc
from metawrap2.command import runner
from metawrap2.modules import bin_refinement

REPORT_HEADER = (
    "genome\tn_genes_called\tn_genes_mapped\ttaxonomic_level\t"
    "clade_separation_score\tpass.GUNC\n"
)


def _report(tmp_path, rows, name="GUNC.progenomes_2.1.maxCSS_level.tsv"):
    out = tmp_path / "gunc"
    out.mkdir(exist_ok=True)
    (out / name).write_text(REPORT_HEADER + "".join(rows))
    return str(out)


def test_command_template_formats_cleanly():
    cmd = gunc.GUNC_RUN.format(bins="b", ext=".fasta", db="db.dmnd", out="o", threads=4)
    assert "{" not in cmd
    # The dot matters: GUNC derives the genome name with basename.split(file_suffix)[0], so
    # "fasta" would name every bin "bin_001." and match nothing in the stats table.
    assert "--file_suffix .fasta" in cmd


def test_threshold_matches_guncs_own():
    assert gunc.CSS_THRESHOLD == 0.45


def test_report_is_found_whatever_the_database_version_is(tmp_path):
    out = _report(tmp_path, [], name="GUNC.some_future_db_9.9.maxCSS_level.tsv")
    assert gunc.find_report(out).endswith("GUNC.some_future_db_9.9.maxCSS_level.tsv")
    assert gunc.find_report(str(tmp_path / "nothing-here")) is None


def test_parse_report_prefers_guncs_own_verdict_over_the_threshold(tmp_path):
    # Row two is below the CSS threshold but GUNC still failed it - too few of its genes were
    # represented in the reference for the score to mean anything. Recomputing from CSS alone
    # would call it clean, which is the opposite of what GUNC said.
    out = _report(
        tmp_path,
        [
            "bin_001\t900\t850\tkingdom\t0.9100\tFalse\n",
            "bin_002\t900\t20\tphylum\t0.0100\tFalse\n",
            "bin_003\t900\t880\tphylum\t0.0200\tTrue\n",
        ],
    )
    verdicts = gunc.parse_report(gunc.find_report(out))
    assert [v.passed for v in (verdicts["bin_001"], verdicts["bin_002"], verdicts["bin_003"])] == [
        False,
        False,
        True,
    ]
    assert verdicts["bin_001"].css == pytest.approx(0.91)
    assert verdicts["bin_001"].level == "kingdom"


def test_parse_report_falls_back_to_the_threshold_without_a_pass_column(tmp_path):
    out = tmp_path / "gunc"
    out.mkdir()
    (out / "GUNC.progenomes_2.1.maxCSS_level.tsv").write_text(
        "genome\ttaxonomic_level\tclade_separation_score\n"
        "clean\tphylum\t0.4000\n"
        "chimeric\tphylum\t0.5000\n"
    )
    verdicts = gunc.parse_report(str(out / "GUNC.progenomes_2.1.maxCSS_level.tsv"))
    assert verdicts["clean"].passed is True
    assert verdicts["chimeric"].passed is False


def test_parse_report_skips_a_bin_it_could_not_score(tmp_path):
    out = _report(
        tmp_path, ["bin_001\t0\t0\tNA\tNA\tFalse\n", "bin_002\t9\t9\tphylum\t0.1\tTrue\n"]
    )
    verdicts = gunc.parse_report(gunc.find_report(out))
    # Absent, not present-and-passing: an unscored bin must read as "not assessed".
    assert "bin_001" not in verdicts
    assert "bin_002" in verdicts


def _stats(tmp_path):
    path = tmp_path / "metawrap_70_10_bins.stats.tsv"
    path.write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin_001\t95.0\t1.0\t0.5\tBacteria\t10000\t3000000\n"
        "bin_002\t90.0\t2.0\t0.4\tBacteria\t20000\t4000000\n"
        "bin_003\t80.0\t3.0\t0.6\tBacteria\t30000\t5000000\n"
    )
    return path


def test_annotate_stats_adds_two_columns_and_leaves_unscored_bins_blank(tmp_path):
    stats = _stats(tmp_path)
    verdicts = {
        "bin_001": gunc.Verdict(css=0.01, passed=True, level="phylum"),
        "bin_002": gunc.Verdict(css=0.91, passed=False, level="kingdom"),
    }
    assert gunc.annotate_stats(str(stats), verdicts) == 2

    rows = [line.split("\t") for line in stats.read_text().splitlines()]
    assert rows[0][-2:] == ["gunc_css", "gunc_pass"]
    assert rows[1][-2:] == ["0.0100", "yes"]
    assert rows[2][-2:] == ["0.9100", "no"]
    # bin_003 was never scored, so it is blank rather than passed or failed.
    assert rows[3][-2:] == ["", ""]


def test_annotate_stats_is_idempotent(tmp_path):
    stats = _stats(tmp_path)
    verdicts = {"bin_001": gunc.Verdict(css=0.01, passed=True, level="phylum")}
    gunc.annotate_stats(str(stats), verdicts)
    first = stats.read_text()
    # A --resume that re-runs the step must overwrite the columns, not append a second pair and
    # leave a table nothing can parse.
    gunc.annotate_stats(str(stats), verdicts)
    assert stats.read_text() == first
    assert stats.read_text().splitlines()[0].count("gunc_css") == 1


def test_annotate_stats_can_change_a_verdict_on_re_run(tmp_path):
    stats = _stats(tmp_path)
    gunc.annotate_stats(str(stats), {"bin_001": gunc.Verdict(css=0.9, passed=False, level="k")})
    gunc.annotate_stats(str(stats), {"bin_001": gunc.Verdict(css=0.1, passed=True, level="p")})
    assert stats.read_text().splitlines()[1].endswith("\t0.1000\tyes")


def test_annotate_stats_tolerates_a_missing_table(tmp_path):
    assert gunc.annotate_stats(str(tmp_path / "not-there.stats.tsv"), {}) == 0


def test_run_gunc_returns_nothing_under_dry_run(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "dry_run", True)
    assert (
        gunc.run_gunc(
            str(tmp_path / "bins"), threads=2, db="db.dmnd", env=None, out_dir=str(tmp_path / "o")
        )
        is None
    )
    # Nothing was created, because nothing ran.
    assert not os.path.isdir(str(tmp_path / "o"))


def test_bin_refinement_has_a_way_to_opt_out():
    args = bin_refinement._parse_args(["-o", "out", "-A", "a", "-B", "b"])
    assert args.run_gunc is True, "GUNC is on by default; --skip-gunc opts out"
    assert (
        bin_refinement._parse_args(["-o", "o", "-A", "a", "-B", "b", "--skip-gunc"]).run_gunc
        is False
    )


def test_gunc_ships_in_bin_refinements_own_env_because_it_runs_by_default():
    from metawrap2 import config

    # GUNC used to have an env of its own. It runs unless --skip-gunc is passed, so that made a
    # default step depend on a second install command - and it solves with CheckM1 anyway.
    assert config.MODULE_ENVS["bin_refinement"] == "metawrap2-binqc"
    assert "gunc" in config.MODULE_TOOLS["bin_refinement"]
    assert "gunc" not in config.OPTIONAL_ENVS
    assert "GUNC_DB" in config.DB_KEYS


def test_install_db_knows_gunc_and_points_the_key_at_the_file():
    from metawrap2.commands import install_db

    db = install_db.DATABASES["gunc"]
    assert db.key == "GUNC_DB"
    # GUNC takes --db_file, so the config value must name the .dmnd, not the directory holding it.
    assert db.key_file and db.key_file.endswith(".dmnd")
    assert "gunc" in install_db.INSTALLERS


@pytest.mark.conda_envs
def test_gunc_scoring_resolves_its_env_through_module_envs():
    """The benchmark scorer must not look GUNC up by an opt-in key that no longer exists.

    GUNC had an env of its own, and `assess.gunc_settings` indexed OPTIONAL_ENVS by "gunc". Folding
    GUNC into bin_refinement's env turned that into a KeyError - and one that only fired on a variant
    whose bins had not already been assessed, so scoring crashed partway through a benchmark chain
    after the expensive part had succeeded.
    """
    from metawrap2 import config
    from metawrap2.benchmark import assess

    env, _db = assess.gunc_settings()
    assert env == config.MODULE_ENVS["bin_refinement"] == "metawrap2-binqc"


def test_gunc_scoring_honours_use_conda_envs_false(monkeypatch):
    from metawrap2.benchmark import assess

    class _Settings:
        use_conda_envs = False

        def db(self, _key):
            return "/db/gunc.dmnd"

    monkeypatch.setattr(assess._config, "load_settings", lambda: _Settings())
    assert assess.gunc_settings() == (None, "/db/gunc.dmnd")
