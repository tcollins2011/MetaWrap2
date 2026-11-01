"""`check --for "<module command>"` and `install-db --for` / `--for-sheet`.

Both read which databases a command needs from MODULE_DB_NEEDS and the flags given. No conda env,
tool, network or database is used: env lookups and installers are stubbed, and databases are
temporary directories.
"""

from __future__ import annotations

import pytest

from metawrap2 import command, dbcheck
from metawrap2.commands import check, install_db


def _config(tmp_path, **databases):
    path = tmp_path / "metawrap2.toml"
    lines = ["[databases]"] + ['%s = "%s"' % (k, v) for k, v in databases.items()]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# -- reading the command line ---------------------------------------------------------------------


def test_the_program_name_and_global_flags_are_dropped():
    module, argv = dbcheck.split_command(
        "metawrap2 --dry-run --max-threads 4 --max-gpus=1 annotate --kofam -b 'my bins' -o OUT"
    )
    assert module == "annotate"
    assert argv == ["--kofam", "-b", "my bins", "-o", "OUT"]


def test_an_unknown_module_is_named_with_a_suggestion():
    with pytest.raises(ValueError, match="'annotat' is not a module. Did you mean 'annotate'"):
        dbcheck.split_command("annotat -b BINS")
    with pytest.raises(ValueError, match="names no module"):
        dbcheck.split_command("metawrap2 --dry-run")
    with pytest.raises(ValueError, match="cannot split"):
        dbcheck.split_command("annotate -b 'unclosed")


# -- check --for ----------------------------------------------------------------------------------


@pytest.fixture
def envs(monkeypatch):
    """Every env exists unless its name is added to the returned set."""
    absent: set = set()
    monkeypatch.setattr(check, "conda_env_exists", lambda env: env not in absent)
    return absent


def test_check_for_lists_only_what_the_command_needs(envs, tmp_path, capsys):
    kofam = tmp_path / "kofam"
    kofam.mkdir()
    rc = check.main(
        [
            "--config",
            _config(tmp_path, KOFAM_DB=kofam),
            "--for",
            "annotate --kofam --cazy -b BINS -o OUT",
        ]
    )
    out = capsys.readouterr().out
    assert rc == 1
    rows = [line.split()[:2] for line in out.splitlines() if line.startswith("  ")]
    assert rows == [["OK", "KOFAM_DB"], ["MISSING", "DBCAN_DB"], ["OK", "metawrap2-annotate"]]
    assert "-> metawrap2 install-db dbcan" in out
    assert "PFAM_DB" not in out and "BAKTA_DB" not in out
    assert out.rstrip().splitlines()[-1] == "ready: no (1 missing)"


def test_check_for_is_ready_when_everything_is_present(envs, tmp_path, capsys):
    genomad = tmp_path / "genomad"
    genomad.mkdir()
    config = _config(tmp_path, GENOMAD_DB=genomad)
    assert check.main(["--config", config, "--for", "viruses --skip-checkv -a A -o O"]) == 0
    assert capsys.readouterr().out.rstrip().splitlines()[-1] == "ready: yes"


def test_check_for_counts_a_missing_env(envs, tmp_path, capsys):
    envs.add("metawrap2-gtdbtk")
    gtdb = tmp_path / "gtdb"
    (gtdb / "taxonomy").mkdir(parents=True)
    config = _config(tmp_path, GTDBTK_DATA_PATH=gtdb)
    assert check.main(["--config", config, "--for", "classify_bins --gtdbtk -b B -o O"]) == 1
    out = capsys.readouterr().out
    assert "-> metawrap2 install-env gtdbtk" in out
    assert "BLASTDB" not in out
    assert out.rstrip().splitlines()[-1] == "ready: no (1 missing)"


def test_checkm2_env_follows_the_scorer_flags():
    assert "metawrap2-checkm2" in check.envs_for_command("bin_refinement", [])
    assert "metawrap2-checkm2" not in check.envs_for_command("bin_refinement", ["--checkm1"])


def test_an_optional_database_does_not_count(envs, tmp_path, capsys):
    gunc = tmp_path / "gunc.dmnd"
    gunc.write_text("x")
    config = _config(tmp_path, GUNC_DB=gunc)
    assert check.main(["--config", config, "--for", "bin_refinement -A a -o o"]) == 0
    out = capsys.readouterr().out
    assert "NOT SET  CHECKM2DB" in out and "ready: yes" in out


def test_check_for_rejects_an_unknown_module(envs, capsys):
    assert check.main(["--for", "binnig -a A"]) == 1
    assert "Did you mean 'binning'" in capsys.readouterr().out


# -- install-db --for / --for-sheet ---------------------------------------------------------------


@pytest.fixture
def installing(monkeypatch):
    """Record the databases install-db would install, without installing anything."""
    installed: list = []

    def fake(name, target, small, jobs):
        installed.append(name)
        return True, []

    monkeypatch.setattr(install_db, "_install_capturing", fake)
    monkeypatch.setattr(install_db.shutil, "which", lambda program: "/usr/bin/" + program)
    monkeypatch.setattr(install_db._interrupt, "install_handlers", lambda: None)
    return installed


def test_install_for_a_module_installs_what_it_needs(installing, tmp_path, capsys):
    config = _config(tmp_path)
    args = ["--config", config, "--dir", str(tmp_path / "dbs"), "--no-config-update"]
    assert install_db.main(args + ["--for", "viruses", "--for", "read_qc"]) == 0
    assert sorted(installing) == ["checkv", "genomad", "host"]
    out = capsys.readouterr().out
    assert "Databases needed by viruses:" in out
    assert "total" in out


def test_install_for_skips_a_database_already_set_up(installing, tmp_path, capsys):
    genomad = tmp_path / "genomad"
    genomad.mkdir()
    config = _config(tmp_path, GENOMAD_DB=genomad)
    args = ["--config", config, "--dir", str(tmp_path / "dbs"), "--no-config-update"]
    assert install_db.main(args + ["--for", "viruses --skip-checkv"]) == 0
    assert installing == []
    out = capsys.readouterr().out
    assert "GENOMAD_DB: already set up at %s" % genomad in out
    assert "Nothing to install" in out


def test_install_for_a_module_whose_databases_are_opt_in(installing, tmp_path, capsys):
    args = ["--config", _config(tmp_path), "--dir", str(tmp_path / "dbs")]
    assert install_db.main(args + ["--for", "annotate"]) == 0
    assert installing == []
    assert 'as in --for "annotate --bakta"' in capsys.readouterr().out


def test_install_for_an_unknown_module_fails(installing, tmp_path, capsys):
    assert install_db.main(["--config", _config(tmp_path), "--for", "virus"]) == 1
    assert "Did you mean 'viruses'" in capsys.readouterr().out


def _sheet(tmp_path, modules: str) -> str:
    path = tmp_path / "samples.toml"
    path.write_text(
        "[settings]\ncoassemble = true\n\n[modules]\n%s\n\n"
        '[[samples]]\nname = "s1"\nr1 = "a_1.fq"\nr2 = "a_2.fq"\n' % modules
    )
    return str(path)


def test_install_for_sheet_uses_the_planned_steps_and_sheet_flags(
    installing, tmp_path, capsys, monkeypatch
):
    monkeypatch.setattr(command.runner, "dry_run", True)
    sheet = _sheet(tmp_path, 'bin_refinement = "--skip-gunc"\nread_qc = "--skip-bmtagger"')
    args = ["--config", _config(tmp_path), "--dir", str(tmp_path / "dbs")]
    assert install_db.main(args + ["--for-sheet", sheet]) == 0
    out = capsys.readouterr().out
    # CheckM2 scores the bins by default; GUNC and the host index are skipped by sheet flags.
    plan = out.split("[dry-run]")[0]
    assert "checkm2" in plan
    assert "gunc" not in plan and "host" not in plan
    assert installing == []


def test_install_for_sheet_reports_an_unreadable_sheet(installing, tmp_path, capsys):
    args = ["--config", _config(tmp_path), "--for-sheet", str(tmp_path / "missing.toml")]
    assert install_db.main(args) == 1
    assert "Cannot read" in capsys.readouterr().out


def test_check_for_binning_follows_the_preset():
    """With no binner named, the preset decides the envs: default and deep include SemiBin2."""
    from metawrap2.commands import check

    learned = "metawrap2-learned"
    assert learned in check.envs_for_command("binning", ["-a", "asm.fa", "-o", "out"])
    assert learned in check.envs_for_command("binning", ["--preset", "deep", "-o", "out"])
    assert learned in check.envs_for_command("binning", ["--preset=gpu"])
    assert learned not in check.envs_for_command("binning", ["--preset", "very-fast"])
    assert learned not in check.envs_for_command("binning", ["--metabat2", "--maxbin2"])
