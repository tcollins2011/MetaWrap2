"""Two modules, one conda environment.

`coverage` runs the aligner, samtools and jgi_summarize_bam_contig_depths that `binning` installs, so
it maps to metawrap2-bins rather than to an env of its own: a second solve holding the same four
tools is a second thing to keep current for no gain.

That only works if install-env resolves a module's yaml from the *environment* it maps to rather than
from the module's name. It used to use the module name, which would have made `install-env --all`
print "No env file for module 'coverage'" and return 1 - the one command the installation docs tell
people to run. What is pinned here is that doctor still finds the opt-in envs after consolidation.
"""

from __future__ import annotations

from metawrap2.config import OPTIONAL_ENVS


def test_doctor_still_reports_the_opt_in_envs_after_they_were_consolidated(monkeypatch):
    """The opt-in envs are owned by a stated module, not by a dash in their name.

    `doctor` used to find the owning module with `env.split("-")[1]`, which worked only while every
    opt-in env was called `metawrap2-<module>-<tool>`. Consolidating the three learned binners into
    `metawrap2-learned` left that with no second dash, so the lookup produced an empty owner, every
    opt-in env was filtered out, and the whole section printed nothing - without an error.
    """
    from metawrap2.commands import doctor
    from metawrap2.config import OPTIONAL_ENV_OWNERS

    monkeypatch.setattr(doctor, "conda_env_exists", lambda env: False)
    rows = doctor._check_optional_envs(["binning", "bin_refinement", "classify_bins"])
    flags = {row["flag"] for row in rows}
    assert {"semibin2", "vamb", "comebin", "checkm2", "gtdbtk"} <= flags
    # AMBER has no module: it is reached by `benchmark score`, so it is not a module's install step.
    assert "amber" not in flags
    # Every flag resolves to a real env, and every owner to a module doctor knows about.
    assert set(OPTIONAL_ENV_OWNERS) == set(OPTIONAL_ENVS)
