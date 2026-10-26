def test_install_env_under_dry_run_prints_the_build_and_runs_nothing(monkeypatch, capsys):
    """`--dry-run install-env` prints the environment builds and runs nothing.

    It used to call `mamba env create` directly; four dry runs built 5 GB of environments.
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
    """Long command names such as `quantify_function` stay separate from their descriptions."""
    from metawrap2 import cli

    for line in cli.help_message().splitlines():
        if not line.startswith("\t") or line.lstrip().startswith("-"):
            continue
        name = line.strip().split()[0]
        if name in cli.MODULES or name in cli.UTILITIES:
            assert line.strip()[len(name) : len(name) + 2] == "  ", line
