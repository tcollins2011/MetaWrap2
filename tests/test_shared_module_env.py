"""Two modules, one conda environment.

`coverage` runs the aligner, samtools and jgi_summarize_bam_contig_depths that `binning` installs, so
it maps to metawrap2-bins rather than to an env of its own: a second solve holding the same four
tools is a second thing to keep current for no gain.

That only works if install-env resolves a module's yaml from the *environment* it maps to rather than
from the module's name. It used to use the module name, which would have made `install-env --all`
print "No env file for module 'coverage'" and return 1 - the one command the installation docs tell
people to run. These tests pin both halves: the resolution, and that a shared env is built once.
"""

from __future__ import annotations

import os

from metawrap2.commands import install_env as ie
from metawrap2.config import MODULE_ENVS, OPTIONAL_ENVS


def test_a_modules_yaml_is_named_after_its_env_not_after_the_module():
    assert os.path.basename(ie._env_yaml("metawrap2-bins")) == "bins.yaml"
    assert os.path.basename(ie._env_yaml(ie._env_for("coverage"))) == "bins.yaml"
    # Opt-in envs are named after their yaml and reached by that name, not through MODULE_ENVS.
    assert os.path.basename(ie._env_yaml(ie._env_for("learned"))) == ("learned.yaml")


def test_every_installable_name_has_a_yaml_so_install_env_all_cannot_fail():
    """`install-env --all` iterates MODULE_ENVS, and a key with no yaml makes it return 1."""
    missing = [
        module
        for module in sorted(MODULE_ENVS)
        if not os.path.isfile(ie._env_yaml(ie._env_for(module)))
    ]
    assert not missing, "no env yaml for %s" % ", ".join(missing)
    missing_optional = [
        key for key, env in sorted(OPTIONAL_ENVS.items()) if not os.path.isfile(ie._env_yaml(env))
    ]
    assert not missing_optional, "no env yaml for %s" % ", ".join(missing_optional)


def test_a_shared_env_is_created_once_not_once_per_module(monkeypatch, capsys):
    """Two mamba solves writing the same prefix at the same time, otherwise."""
    built = []
    monkeypatch.setattr(ie, "ensure_mamba", lambda: True)
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: False)
    monkeypatch.setattr(
        ie,
        "_create_one",
        lambda module, env, yaml, force, log_dir: built.append(env) or (module, True, ""),
    )

    rc = ie.main(["binning", "coverage", "--no-test"])

    assert rc == 0
    assert built == ["metawrap2-bins"]
    assert "shares it" in capsys.readouterr().out


def test_coverage_alone_builds_the_binning_env(monkeypatch):
    built = []
    monkeypatch.setattr(ie, "ensure_mamba", lambda: True)
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: False)
    monkeypatch.setattr(
        ie,
        "_create_one",
        lambda module, env, yaml, force, log_dir: built.append((env, os.path.basename(yaml)))
        or (module, True, ""),
    )

    assert ie.main(["coverage", "--no-test"]) == 0
    assert built == [("metawrap2-bins", "bins.yaml")]


def test_doctor_still_reports_the_opt_in_envs_after_they_were_consolidated(monkeypatch):
    """The opt-in envs are owned by a stated module, not by a dash in their name.

    `doctor` used to find the owning module with `env.split("-")[1]`, which worked only while every
    opt-in env was called `metawrap2-<module>-<tool>`. Consolidating the three learned binners into
    `metawrap2-learned` left that with no second dash, so the lookup produced an empty owner, every
    opt-in env was filtered out, and the whole section printed nothing - without an error.
    """
    from metawrap2.commands import doctor
    from metawrap2.config import OPTIONAL_ENV_OWNERS, OPTIONAL_ENVS

    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rows = doctor._check_optional_envs(["binning", "bin_refinement", "classify_bins"])
    flags = {row["flag"] for row in rows}
    assert {"semibin2", "vamb", "comebin", "checkm2", "gtdbtk"} <= flags
    # AMBER has no module: it is reached by `benchmark score`, so it is not a module's install step.
    assert "amber" not in flags
    # Every flag resolves to a real env, and every owner to a module doctor knows about.
    assert set(OPTIONAL_ENV_OWNERS) == set(OPTIONAL_ENVS)
