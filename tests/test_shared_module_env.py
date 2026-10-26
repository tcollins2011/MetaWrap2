"""Two modules, one conda environment.

`coverage` uses the aligner, samtools and jgi_summarize_bam_contig_depths that `binning` installs,
so it maps to metawrap2-bins instead of its own env with the same four tools.

install-env must therefore resolve a module's yaml from the env it maps to, not the module name;
by module name, `install-env --all` would print "No env file for module 'coverage'" and return 1.
This file tests that doctor still finds the opt-in envs after consolidation.
"""

from __future__ import annotations

from metawrap2.config import OPTIONAL_ENVS


def test_doctor_still_reports_the_opt_in_envs_after_they_were_consolidated(monkeypatch):
    """Each opt-in env has an explicit owning module instead of one parsed from its name.

    `doctor` used `env.split("-")[1]`, which assumed `metawrap2-<module>-<tool>`. After the three
    learned binners moved to `metawrap2-learned`, the owner came out empty, every opt-in env was
    filtered out, and the section printed nothing.
    """
    from metawrap2.commands import doctor
    from metawrap2.config import OPTIONAL_ENV_OWNERS

    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rows = doctor._check_optional_envs(["binning", "bin_refinement", "classify_bins"])
    flags = {row["flag"] for row in rows}
    assert {"semibin2", "vamb", "comebin", "checkm2", "gtdbtk"} <= flags
    # AMBER has no module; it is used by `benchmark score`.
    assert "amber" not in flags
    # Every flag resolves to a real env, and every owner to a module doctor knows about.
    assert set(OPTIONAL_ENV_OWNERS) == set(OPTIONAL_ENVS)


def test_doctor_counts_a_missing_checkm2_env_as_a_failure(monkeypatch):
    """CheckM2 is the default bin scorer, so its env is not reported as opt-in."""
    from metawrap2.commands import doctor

    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rows = {row["flag"]: row for row in doctor._check_optional_envs(["reassemble_bins"])}
    assert rows["checkm2"]["default"] is True
    assert "install-env checkm2" in rows["checkm2"]["detail"]
    rows = {row["flag"]: row for row in doctor._check_optional_envs(["classify_bins"])}
    assert rows["gtdbtk"]["default"] is False
    assert "checkm2" not in rows
