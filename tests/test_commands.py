def test_install_env_under_dry_run_prints_the_build_and_runs_nothing(monkeypatch, capsys):
    """`--dry-run install-env` used to solve and build every environment it named.

    The global flag promises the commands are printed rather than run, and every other command
    keeps that promise. This one went straight to `mamba env create`, so checking what an install
    would do cost the install - four "dry runs" built 5 GB of environments.
    """
    from metawrap2 import command
    from metawrap2.commands import install_env

    monkeypatch.setattr(command.runner, "dry_run", True)
    monkeypatch.setattr(install_env, "ensure_mamba", lambda: True)
    monkeypatch.setattr(install_env, "conda_env_exists", lambda env: False)

    def must_not_run(*_a, **_k):
        raise AssertionError("a dry run started a real environment build")

    monkeypatch.setattr(install_env, "_create_one", must_not_run)
    monkeypatch.setattr(install_env.subprocess, "call", must_not_run)
    assert install_env.main(["checkm2", "--no-test"]) == 0
    out = capsys.readouterr().out
    assert "[dry-run] mamba env create -y -f" in out and "checkm2.yaml" in out


def test_the_help_table_keeps_every_name_apart_from_its_description():
    """A fixed 16-character column ran `quantify_function` straight into its description."""
    from metawrap2 import cli

    for line in cli.help_message().splitlines():
        if not line.startswith("\t") or line.lstrip().startswith("-"):
            continue
        name = line.strip().split()[0]
        if name in cli.MODULES or name in cli.UTILITIES:
            assert line.strip()[len(name) : len(name) + 2] == "  ", line
