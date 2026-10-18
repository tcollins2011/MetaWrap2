"""`install-db` defects found in one real install session, each with the symptom it produced.

Every test here stands for a database install that a user could not complete, or completed and
then could not use: CARD refusing its own tarball, an interrupted run losing the config for what
it had already downloaded, a failure that named no reason, a size table that could not be added
up, and five advertised databases the command did not know about.

Nothing here touches the network. The tool-managed installers are exercised by replacing
``_mamba_run``, which is the only thing between them and a 110 GB download.
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

#: The five that were advertised in the README's table and in usage.md but were not installable,
#: with the config key each module actually reads.
TOOL_MANAGED = (
    ("amrfinder", "AMRFINDER_DB"),
    ("genomad", "GENOMAD_DB"),
    ("checkv", "CHECKV_DB"),
    ("bakta", "BAKTA_DB"),
    ("gtdbtk", "GTDBTK_DATA_PATH"),
)


@pytest.fixture
def mamba_present(monkeypatch):
    """Pretend mamba is on PATH, so a test of install-db does not depend on this machine's.

    `main` refuses up front when a requested database is downloaded by a tool from a conda env and
    mamba is missing, which is right for a user and would otherwise make these tests pass or fail
    according to whether the developer happens to have conda installed.
    """
    monkeypatch.setattr(install_db.shutil, "which", lambda program: "/usr/bin/" + program)


@pytest.fixture
def envs_present(monkeypatch):
    """Pretend the conda env that owns each downloader exists.

    `_run_downloader` checks the env before shelling out, so a user on a fresh machine is told to
    run `install-env` rather than reading mamba's "The given prefix does not exist". Without this
    fixture these tests would pass on a developer's machine - where the envs happen to be built -
    and fail in CI, which is the worse of the two outcomes.
    """
    monkeypatch.setattr(install_db, "env_prefix", lambda env: "/fake/envs/" + env)


def test_card_installs_from_a_tarball_whose_root_member_is_dot(tmp_path, monkeypatch, capsys):
    """CARD's tarball names its root as `.`, and rejecting that made `amr --rgi` uninstallable.

    The guard this exercises is the escape check, which must keep refusing `..` and absolute
    members; what it must not do is refuse the one member that cannot escape. If this breaks,
    `install-db card` fails 100% of the time with "refusing to extract '.'".
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
    """The guard CARD now shares has to keep doing its job: these archives are untrusted input."""
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
    """An interrupted run must keep the config for everything that already finished.

    The reporter stopped a 9-database run after 4 had completed - 17 GB on disk - and the config
    named none of them, because the single write happened after the whole pool. If this breaks,
    every interrupted `--all` leaves databases that are installed and unusable.
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
    """A failure has to carry its reason, and say that a retry only redoes what is missing.

    The reporter's whole diagnostic for a failed database was "FAILED: dbcan": the subprocess's
    error went to a shared stdout interleaved with two other databases' progress meters. Re-running
    the identical command by hand worked, so the information that it was transient - and that a
    retry skips the 4 databases that had succeeded - was exactly what got lost.
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
    """`--list` is authoritative, so every row needs both numbers and a column that totals.

    One "full size" column used to mix a bare size, a download size, a download-and-unpacked pair
    and "~20 GB + indexing"; only 4 of 13 rows said how much disk the install needed, and the
    two-number rows overflowed into the `--small` column. GUNC was the expensive case: "6.7 GB
    download" for about 13 GB on disk, in a database that runs by default in bin_refinement.
    """
    for name, db in install_db.DATABASES.items():
        download = install_db._gigabytes(db.download)
        on_disk = install_db._gigabytes(db.on_disk)
        assert download > 0, name
        assert on_disk > 0, name
        # Unpacking and indexing never shrink a database, so a row the wrong way round is a typo.
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
    """2.2 GB at 150-300 KB/s is 2-4 hours; without the note it reads as a hung install."""
    install_db._print_list()
    out = capsys.readouterr().out
    assert "150-300 KB/s" in out


def test_the_five_advertised_databases_are_installable_by_name():
    """The README's table had 18 rows and install-db knew 13 of them.

    `install-db gtdbtk` answered "Unknown database(s): gtdbtk", which left `amr`, `viruses`,
    `annotate --bakta` and `classify_bins --gtdbtk` with no documented way to get their data.
    """
    for name, key in TOOL_MANAGED:
        assert name in install_db.DATABASES, name
        assert name in install_db.INSTALLERS, name
        assert install_db.DATABASES[name].key == key
        # The key has to be the one the modules already read, not a new spelling of it.
        assert key in DB_KEYS, key
        assert install_db.DATABASES[name].used_by


def test_each_tool_managed_installer_runs_its_downloader_in_its_own_env(
    tmp_path, monkeypatch, envs_present
):
    """Each of the five runs the tool's own downloader, in the env that owns that tool.

    The env comes from config.py rather than being spelled out here for the reason the CheckM
    installer records: an env rename left a hardcoded name pointing at one that no longer existed,
    and a `mamba run -n` against a missing env fails in a way that is easy to miss.
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
            # What the real downloader leaves behind, which is what makes the install idempotent.
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

        # And a second run asks the tool for nothing at all.
        calls.clear()
        assert install_db.INSTALLERS[name](str(target), False, 1) is True, name
        assert calls == [], name


def test_a_tool_managed_install_configures_the_path_the_module_reads(
    tmp_path, monkeypatch, capsys, mamba_present, envs_present
):
    """geNomad writes genomad_db/ inside the destination, and that inner path is what runs.

    `genomad end-to-end` takes the database directory itself, so GENOMAD_DB pointing at its parent
    is a config that looks installed and fails at the first viral call.
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
    """`--small` for bakta is a real 1.3 GB Bakta configuration, and it lands in db-light/.

    The two types live in different directories under one output directory, so a `--small` install
    configured as `db/` points at a database that is not there.
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

    Set on the runner instance rather than through ``command.set_dry_run`` so that it is undone
    when the test ends: a leaked dry-run flag makes every later test pass without running
    anything, which is the one failure mode a test suite cannot report on itself.
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
    """`metawrap2 --dry-run install-db` must touch nothing at all.

    `_run` went straight to subprocess and nothing in this command consulted command.runner, so
    the global flag was silently ignored here: the first user to type it with `gtdbtk` pulled
    ~110 GB from the one flag whose purpose is to let them decide first.
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
    """The plan has to be specific enough to decide on: what runs, where, how big, which key."""
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
    """The marker checks still run under a dry run, or the plan lies about a 90 GB download."""
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
    """Planning an install is what you do *before* building the envs, so this cannot be fatal."""
    monkeypatch.setattr(install_db.shutil, "which", lambda program: None)

    code = install_db.main(["gtdbtk", "--dir", str(tmp_path / "dbs"), "--config", "/dev/null"])
    out = capsys.readouterr().out

    assert code == 0
    assert "mamba is not on PATH" in out
    assert "install-env --all" in out
    assert "would run `download-db.sh" in out


def test_a_missing_env_names_install_env_rather_than_mambas_prefix(tmp_path, monkeypatch, capsys):
    """`quickstart` builds the default envs, so the opt-in ones are missing on a fresh machine.

    What that looked like: `critical libmamba The given prefix does not exist:
    ".../envs/metawrap2-checkm2"` followed by `FAILED: checkm2`. True, and useless - it names the
    prefix, not the command that creates it. The check also has to happen *before* the downloader
    runs, or the user waits for a shell-out that cannot work.
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
    """The check must not become a second failure mode of its own."""
    monkeypatch.setattr(install_db, "env_prefix", lambda env: "/fake/envs/" + env)
    ran = []
    monkeypatch.setattr(install_db, "_mamba_run", lambda env, cmd: ran.append((env, cmd)) or 0)
    assert install_db._run_downloader("checkm2", str(tmp_path), False) == 0
    assert ran and ran[0][0] == OPTIONAL_ENVS["checkm2"]


def test_every_database_has_a_plan():
    """A registered database with no plan would print a placeholder step instead of its URL."""
    for name in install_db.DATABASES:
        steps = install_db._plan_steps(name, "/nonexistent/%s" % name.upper(), False, False)
        assert steps, name
        assert steps != ["install %s" % name], name
        assert any("http" in step or "conda env" in step for step in steps), name
