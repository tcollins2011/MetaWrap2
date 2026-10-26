"""Tests for `install-db` fixes: CARD extraction, incremental config writes, failure reasons,
the size table, and the tool-managed databases.

Nothing here uses the network. Tool-managed installers are tested by replacing ``_mamba_run``.
"""

import io
import json
import os
import shutil
import sys
import tarfile

import pytest

from metawrap2.commands import install_db
from metawrap2.config import DB_KEYS, MODULE_ENVS, OPTIONAL_ENVS

#: Databases downloaded by their tool's own installer, with the config key each module reads.
TOOL_MANAGED = (
    ("amrfinder", "AMRFINDER_DB"),
    ("genomad", "GENOMAD_DB"),
    ("checkv", "CHECKV_DB"),
    ("bakta", "BAKTA_DB"),
    ("gtdbtk", "GTDBTK_DATA_PATH"),
)


@pytest.fixture
def mamba_present(monkeypatch):
    """Pretend mamba is on PATH, so results do not depend on the host.

    `main` refuses up front when a tool-managed database is requested and mamba is missing.
    """
    monkeypatch.setattr(install_db.shutil, "which", lambda program: "/usr/bin/" + program)


@pytest.fixture
def envs_present(monkeypatch):
    """Pretend the conda env that owns each downloader exists.

    `_run_downloader` checks the env before shelling out, so without this the results would
    depend on which envs are built on the host.
    """
    monkeypatch.setattr(install_db, "env_prefix", lambda env: "/fake/envs/" + env)


def test_card_installs_from_a_tarball_whose_root_member_is_dot(tmp_path, monkeypatch, capsys):
    """CARD's tarball has a `.` root member, which the path escape check must accept.

    Rejecting it made `install-db card` always fail with "refusing to extract '.'".
    """
    archive = tmp_path / "card_data.tar.bz2"
    payload = json.dumps({"_version": "4.0.2"}).encode()
    with tarfile.open(archive, "w:bz2") as tf:
        root = tarfile.TarInfo(".")
        root.type = tarfile.DIRTYPE
        tf.addfile(root)
        for name in ("./card.json", "./nucleotide_fasta_protein_homolog_model.fasta"):
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            tf.addfile(member, io.BytesIO(payload))

    def fetch(url, dest):
        shutil.copy(archive, dest)
        return True

    monkeypatch.setattr(install_db, "_fetch", fetch)
    target = tmp_path / "CARD"

    assert install_db._install_card(str(target), False) is True
    assert (target / "card.json").is_file()
    out = capsys.readouterr().out
    assert "refusing to extract" not in out
    assert "CARD version 4.0.2" in out


def test_an_escaping_member_is_still_refused(tmp_path):
    """`_extract` still refuses `..` and absolute members; archives are untrusted input."""
    archive = tmp_path / "db.tar"
    with tarfile.open(archive, "w") as tf:
        for name in ("../escaped.dmp", "/absolute.dmp", "./names.dmp"):
            member = tarfile.TarInfo(name)
            member.size = 2
            tf.addfile(member, io.BytesIO(b"ok"))
    into = tmp_path / "db"

    assert install_db._extract(str(archive), str(into)) is True
    assert (into / "names.dmp").is_file()
    assert not (tmp_path / "escaped.dmp").exists()
    assert not (into / "absolute.dmp").exists()


def test_config_keys_are_written_as_each_database_lands(tmp_path, monkeypatch, capsys):
    """Each database's config key is written when it finishes, so an interrupted run keeps them.

    The config used to be written once after all installs, so interrupting a run lost the keys
    for databases already on disk.
    """
    config = tmp_path / "config.toml"
    databases = tmp_path / "dbs"

    def completes(target, small, jobs):
        os.makedirs(target, exist_ok=True)
        return True

    def interrupted(target, small, jobs):
        raise KeyboardInterrupt

    monkeypatch.setitem(install_db.INSTALLERS, "kegg", completes)
    monkeypatch.setitem(install_db.INSTALLERS, "card", interrupted)

    with pytest.raises(KeyboardInterrupt):
        install_db.main(
            ["kegg", "card", "--dir", str(databases), "--config", str(config), "-t", "1"]
        )

    capsys.readouterr()
    assert 'KEGG_BRITE = "%s"' % (databases / "KEGG") in config.read_text()


def test_a_key_written_as_it_landed_is_not_written_twice(tmp_path, monkeypatch, capsys):
    """The end-of-run summary must not duplicate the incremental writes into [databases]."""
    config = tmp_path / "config.toml"
    databases = tmp_path / "dbs"

    def completes(target, small, jobs):
        os.makedirs(target, exist_ok=True)
        return True

    monkeypatch.setitem(install_db.INSTALLERS, "kegg", completes)
    assert install_db.main(["kegg", "--dir", str(databases), "--config", str(config)]) == 0

    capsys.readouterr()
    text = config.read_text()
    assert text.count("KEGG_BRITE") == 1


def test_a_failed_database_quotes_the_output_that_explains_it(tmp_path, monkeypatch, capsys):
    """A failed database is reported with the tool's error output and a note that a retry skips
    completed databases.

    Subprocess errors used to go to a shared stdout mixed with other databases' progress output.
    """
    message = "run_dbcan: ERROR 503 from the pinned S3 release, giving up"

    def failing(target, small, jobs):
        install_db._run(
            [sys.executable, "-c", "import sys; sys.stderr.write(%r); sys.exit(1)" % (message)]
        )
        return False

    monkeypatch.setitem(install_db.INSTALLERS, "dbcan", failing)
    monkeypatch.setattr(install_db.shutil, "which", lambda program: "/usr/bin/" + program)

    code = install_db.main(
        ["dbcan", "--dir", str(tmp_path / "dbs"), "--config", str(tmp_path / "config.toml")]
    )
    out = capsys.readouterr().out

    assert code == 1
    assert "  FAILED: dbcan" in out
    assert "    " + message in out  # indented under the line it explains
    assert "1 database(s) failed: dbcan" in out
    assert "completed databases are detected and skipped" in out


def test_every_database_reports_a_download_and_an_on_disk_size(capsys):
    """Every `--list` row has a download size and an on-disk size, and both columns total.

    The list used to have one mixed "full size" column; GUNC showed "6.7 GB download" for about
    13 GB on disk.
    """
    for name, db in install_db.DATABASES.items():
        download = install_db._gigabytes(db.download)
        on_disk = install_db._gigabytes(db.on_disk)
        assert download > 0, name
        assert on_disk > 0, name
        # Unpacking and indexing do not shrink a database.
        assert on_disk >= download, name
        assert install_db._gigabytes(db.small_size or db.download) > 0, name

    install_db._print_list()
    out = capsys.readouterr().out

    assert "download" in out and "on disk" in out
    assert "6.7 GB" in out and "~13 GB" in out  # gunc, both numbers
    for name in install_db.DATABASES:
        assert name in out
    totals = [line for line in out.splitlines() if line.strip().startswith("total")]
    assert len(totals) == 1
    assert len(totals[0].split()) == 1 + 3 * 2  # "total" plus three "~N GB" figures


def test_singlem_says_its_download_is_slow_rather_than_stuck(capsys):
    """The list notes SingleM's slow download (2.2 GB at 150-300 KB/s is 2-4 hours)."""
    install_db._print_list()
    out = capsys.readouterr().out
    assert "150-300 KB/s" in out


def test_the_five_advertised_databases_are_installable_by_name():
    """The tool-managed databases are registered with install-db under the modules' config keys.

    They were listed in the README but `install-db` did not know them, so `amr`, `viruses`,
    `annotate --bakta` and `classify_bins --gtdbtk` had no documented way to get their data.
    """
    for name, key in TOOL_MANAGED:
        assert name in install_db.DATABASES, name
        assert name in install_db.INSTALLERS, name
        assert install_db.DATABASES[name].key == key
        # The key must be the one the modules read.
        assert key in DB_KEYS, key
        assert install_db.DATABASES[name].used_by


def test_each_tool_managed_installer_runs_its_downloader_in_its_own_env(
    tmp_path, monkeypatch, envs_present
):
    """Each tool-managed installer runs the tool's downloader in the env that owns the tool.

    Env names come from config.py so that an env rename cannot leave a stale hardcoded name.
    """
    expected = {
        "amrfinder": (
            MODULE_ENVS["amr"],
            ["amrfinder_update", "--database"],
            "2024-07-22.1",
        ),
        "genomad": (
            MODULE_ENVS["viruses"],
            ["genomad", "download-database"],
            "genomad_db",
        ),
        "checkv": (
            MODULE_ENVS["viruses"],
            ["checkv", "download_database"],
            "checkv-db-v1.5",
        ),
        "bakta": (
            MODULE_ENVS["annotate"],
            ["bakta_db", "download", "--output"],
            "db",
        ),
        "gtdbtk": (
            OPTIONAL_ENVS["gtdbtk"],
            ["download-db.sh"],
            "taxonomy",
        ),
    }
    assert set(expected) == {name for name, _ in TOOL_MANAGED}

    for name, (env, argv, creates) in expected.items():
        target = tmp_path / name.upper()
        calls = []

        def downloader(used_env, cmd, _calls=calls, _target=target, _creates=creates):
            _calls.append((used_env, cmd))
            # Mimic the real downloader's output; its presence makes a second install a no-op.
            os.makedirs(os.path.join(str(_target), _creates), exist_ok=True)
            with open(os.path.join(str(_target), _creates, "data"), "w") as handle:
                handle.write("x")
            return 0

        monkeypatch.setattr(install_db, "_mamba_run", downloader)
        assert install_db.INSTALLERS[name](str(target), False, 1) is True, name
        assert len(calls) == 1, name
        assert calls[0][0] == env, name
        assert calls[0][1][: len(argv)] == argv, name
        assert str(target) in calls[0][1], name

        # A second run does not call the downloader.
        calls.clear()
        assert install_db.INSTALLERS[name](str(target), False, 1) is True, name
        assert calls == [], name


def test_a_tool_managed_install_configures_the_path_the_module_reads(
    tmp_path, monkeypatch, capsys, mamba_present, envs_present
):
    """GENOMAD_DB is set to the inner genomad_db/ directory that geNomad writes.

    `genomad end-to-end` takes that directory itself and fails if given its parent.
    """
    config = tmp_path / "config.toml"
    databases = tmp_path / "dbs"

    def downloader(env, cmd):
        inner = os.path.join(str(databases), "GENOMAD", "genomad_db")
        os.makedirs(inner, exist_ok=True)
        with open(os.path.join(inner, "genomad_db.dbtype"), "w") as handle:
            handle.write("x")
        return 0

    monkeypatch.setattr(install_db, "_mamba_run", downloader)
    code = install_db.main(["genomad", "--dir", str(databases), "--config", str(config)])

    capsys.readouterr()
    assert code == 0
    assert 'GENOMAD_DB = "%s"' % (databases / "GENOMAD" / "genomad_db") in config.read_text()


def test_bakta_small_installs_the_light_database_and_configures_that(
    tmp_path, monkeypatch, capsys, mamba_present, envs_present
):
    """`bakta --small` installs Bakta's 1.3 GB light database and configures db-light/.

    The full and light databases are written to different directories (db/ and db-light/).
    """
    config = tmp_path / "config.toml"
    databases = tmp_path / "dbs"
    calls = []

    def downloader(env, cmd):
        calls.append(cmd)
        inner = os.path.join(str(databases), "BAKTA", "db-light")
        os.makedirs(inner, exist_ok=True)
        with open(os.path.join(inner, "bakta.db"), "w") as handle:
            handle.write("x")
        return 0

    monkeypatch.setattr(install_db, "_mamba_run", downloader)
    code = install_db.main(["bakta", "--small", "--dir", str(databases), "--config", str(config)])

    capsys.readouterr()
    assert code == 0
    assert calls[0][-2:] == ["--type", "light"]
    assert 'BAKTA_DB = "%s"' % (databases / "BAKTA" / "db-light") in config.read_text()


# --- the global --dry-run flag -----------------------------------------------------------


@pytest.fixture
def planning(monkeypatch):
    """Put this command under the global ``--dry-run``, the way ``cli`` does.

    Set on the runner instance with monkeypatch, not through ``command.set_dry_run``, so it is
    undone when the test ends. A leaked dry-run flag would make later tests pass without running.
    """
    from metawrap2 import command

    monkeypatch.setattr(command.runner, "dry_run", True)


@pytest.fixture
def no_subprocess(monkeypatch):
    """Make any attempt to run a child process fail the test."""

    def boom(*args, **kwargs):
        raise AssertionError("a dry run must not run anything, but ran %r" % (args,))

    monkeypatch.setattr(install_db, "_run", boom)
    monkeypatch.setattr(install_db.subprocess, "call", boom)
    monkeypatch.setattr(install_db.subprocess, "Popen", boom)


def test_a_dry_run_runs_nothing_and_writes_nothing(
    tmp_path, monkeypatch, capsys, planning, no_subprocess
):
    """`metawrap2 --dry-run install-db` runs no subprocess and writes no files.

    `_run` used to ignore command.runner, so a dry run of `gtdbtk` downloaded ~110 GB.
    """
    config = tmp_path / "config.toml"
    databases = tmp_path / "dbs"
    names = [name for name, _ in TOOL_MANAGED] + ["card"]

    code = install_db.main(names + ["--dir", str(databases), "--config", str(config)])
    out = capsys.readouterr().out

    assert code == 0
    assert not config.exists()  # no config write
    assert not databases.exists()  # not even an empty install directory
    for name in names:
        assert "=== %s" % name in out


def test_a_dry_run_names_the_command_the_destination_and_the_key(
    tmp_path, capsys, planning, no_subprocess
):
    """The dry-run plan shows the command, destination, sizes and config key for each database."""
    databases = tmp_path / "dbs"
    config = tmp_path / "config.toml"

    install_db.main(["gtdbtk", "card", "--dir", str(databases), "--config", str(config)])
    out = capsys.readouterr().out

    # The downloader command, in its env, for a tool-managed database.
    assert "download-db.sh %s" % (databases / "GTDBTK") in out
    assert OPTIONAL_ENVS["gtdbtk"] in out
    # The URL, for one fetched directly.
    assert install_db.CARD_URL in out
    # Sizes from the two-column table, and the key each would set.
    assert "~110 GB download, ~110 GB on disk" in out
    assert "would set GTDBTK_DATA_PATH = %s in %s" % (databases / "GTDBTK", config) in out
    assert "would set CARD_DB = %s in %s" % (databases / "CARD" / "card.json", config) in out
    # A total and the same headroom sentence the real run prints.
    assert "to install: would download" in out
    assert "of headroom)" in out


def test_a_dry_run_reports_a_database_that_is_already_present(
    tmp_path, capsys, planning, no_subprocess
):
    """A dry run checks for an existing install and reports that it would skip it."""
    databases = tmp_path / "dbs"
    (databases / "KRAKEN2").mkdir(parents=True)
    (databases / "KRAKEN2" / "hash.k2d").write_text("x")

    install_db.main(["kraken2", "--dir", str(databases), "--config", str(tmp_path / "c.toml")])
    out = capsys.readouterr().out

    assert "already present, would skip" in out
    assert "genome-idx" not in out  # the Kraken2 URL is not in the plan
    assert "Nothing to do" in out


def test_a_dry_run_plans_even_when_the_conda_env_is_missing(
    tmp_path, monkeypatch, capsys, planning, no_subprocess
):
    """A dry run still prints a plan when mamba or the env is missing."""
    monkeypatch.setattr(install_db.shutil, "which", lambda program: None)

    code = install_db.main(["gtdbtk", "--dir", str(tmp_path / "dbs"), "--config", "/dev/null"])
    out = capsys.readouterr().out

    assert code == 0
    assert "mamba is not on PATH" in out
    assert "install-env --all" in out
    assert "would run `download-db.sh" in out


def test_a_missing_env_names_install_env_rather_than_mambas_prefix(tmp_path, monkeypatch, capsys):
    """A missing env is reported with the `install-env` command, before the downloader runs.

    `quickstart` builds only the default envs, so opt-in envs are missing on a fresh machine.
    mamba's own message ("The given prefix does not exist") does not name the fix.
    """
    monkeypatch.setattr(install_db, "env_prefix", lambda env: None)
    ran = []
    monkeypatch.setattr(install_db, "_mamba_run", lambda env, cmd: ran.append(cmd) or 0)

    assert install_db._run_downloader("checkm2", str(tmp_path), False) == 1
    out = capsys.readouterr().out
    assert ran == [], "the downloader must not be started in an env that does not exist"
    assert "metawrap2 install-env checkm2" in out
    assert "metawrap2-checkm2" in out


def test_a_present_env_still_runs_the_downloader(tmp_path, monkeypatch):
    """The env check lets the downloader run when the env exists."""
    monkeypatch.setattr(install_db, "env_prefix", lambda env: "/fake/envs/" + env)
    ran = []
    monkeypatch.setattr(install_db, "_mamba_run", lambda env, cmd: ran.append((env, cmd)) or 0)
    assert install_db._run_downloader("checkm2", str(tmp_path), False) == 0
    assert ran and ran[0][0] == OPTIONAL_ENVS["checkm2"]


def test_every_database_has_a_plan():
    """Every registered database has a dry-run plan with a URL or conda env step."""
    for name in install_db.DATABASES:
        steps = install_db._plan_steps(name, "/nonexistent/%s" % name.upper(), False, False)
        assert steps, name
        assert steps != ["install %s" % name], name
        assert any("http" in step or "conda env" in step for step in steps), name
