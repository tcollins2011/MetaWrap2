"""External benchmark arms: another pipeline run in place of MetaWrap2 in a comparison.

Covers the argv built from the runner's inputs, the load-time checks on a declaration, cost
measurement, publishing of the arm's bins for scoring, and the run record.
"""

import json
import os

import pytest

from metawrap2.benchmark import assess, registry
from metawrap2.benchmark import runner as runner_mod
from metawrap2.benchmark.registry import Comparison, External, RegistryError, Variant

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks")

#: One step, so a test can read the resolved argv without a shell. The metaWRAP 1 example at the
#: end of configurations.toml is two steps inside one ``bash -c``.
MW1_COMMAND = "metawrap binning -o {output}/out -t {threads} -m {memory} -a {assembly} {reads}"


def _tiny():
    return registry.get_dataset("tiny", os.path.join(FIXTURES, "datasets.toml"))


def _external_variant(command=MW1_COMMAND, **kw):
    options = {
        "command": command,
        "bins": "out/final_bins",
        "env": "metawrap-1.3.0",
        "version": "1.3.0",
    }
    options.update(kw)
    return Variant(name="mw1", external=External(**options))


def _mixed_comparison():
    """A binning comparison with one native arm and one external arm."""
    return Comparison(
        name="binner_set",
        module="binning",
        metrics=("bins_high_quality",),
        variants=(
            Variant(name="metabat2", args=("--metabat2",)),
            _external_variant(),
        ),
    )


def _inputs(tmp_path):
    return runner_mod.Inputs(
        pairs=[("sample_1", str(tmp_path / "s_1.fastq"), str(tmp_path / "s_2.fastq"))],
        assembly=str(tmp_path / "gold.fa"),
    )


def test_an_external_arms_command_is_filled_from_the_runners_own_inputs(tmp_path):
    """An external arm's command is filled with the runner's inputs, budget and output dir."""
    comparison = _mixed_comparison()
    inputs = _inputs(tmp_path)
    out_dir = str(tmp_path / "work" / "binner_set" / "mw1")

    argv = runner_mod.build_argv(comparison, comparison.variant("mw1"), inputs, out_dir, 6, 48)

    assert argv[:2] == ["metawrap", "binning"]
    # The same thread and memory grant as the other arms.
    assert argv[argv.index("-t") + 1] == "6"
    assert argv[argv.index("-m") + 1] == "48"
    # The inputs are the runner's, not the declaration's.
    assert argv[argv.index("-a") + 1] == str(tmp_path / "gold.fa")
    assert argv[-2:] == [str(tmp_path / "s_1.fastq"), str(tmp_path / "s_2.fastq")]
    # {reads} is several files, so substitution happens before the argv is split.
    assert all("{" not in token for token in argv)
    # Its own output goes one level down, away from the published bins.
    assert argv[argv.index("-o") + 1] == os.path.join(out_dir, "external", "out")


def test_a_native_arm_builds_exactly_the_command_it_did_before(tmp_path):
    """A native arm's full argv is unchanged by the external-arm branch in build_argv."""
    comparison = _mixed_comparison()
    inputs = _inputs(tmp_path)

    argv = runner_mod.build_argv(comparison, comparison.variant("metabat2"), inputs, "out", 6, 48)

    assert argv[1:] == [
        "-m",
        "metawrap2.cli",
        "binning",
        "--metabat2",
        "-o",
        "out",
        "-t",
        "6",
        "-a",
        str(tmp_path / "gold.fa"),
        str(tmp_path / "s_1.fastq"),
        str(tmp_path / "s_2.fastq"),
    ]


@pytest.mark.conda_envs
def test_an_external_arm_runs_inside_its_declared_env(tmp_path, monkeypatch):
    """An external arm runs through ``mamba run -n <env>``, the wrapper MetaWrap2 tools use."""
    from metawrap2 import command as command_mod

    monkeypatch.setattr(command_mod.shutil, "which", lambda name: "/opt/mamba/bin/" + name)
    comparison = _mixed_comparison()

    argv = runner_mod.build_argv(
        comparison, comparison.variant("mw1"), _inputs(tmp_path), "out", 6, 48
    )

    assert argv[0] == "mamba" and argv[1] == "run"
    assert argv[argv.index("-n") + 1] == "metawrap-1.3.0"
    assert argv[argv.index("-n") + 2 :][:2] == ["metawrap", "binning"]


# ── refusals at load time ─────────────────────────────────────────────────────────────────────

#: A comparison with one native arm and one external arm, as TOML. The external arm's table is
#: formatted in, so each case below changes one field.
_CONFIG = """
schema = 1

[binner_set]
module = "binning"
metrics = ["bins_high_quality"]

[[binner_set.variants]]
name = "metabat2"
args = ["--metabat2"]
default = true

[[binner_set.variants]]
name = "mw1"
%s

[binner_set.variants.external]
%s
"""

_VALID_EXTERNAL = (
    'command = "metawrap binning -o {output} -t {threads} -a {assembly} {reads}"\n'
    'bins = "out/final_bins"\n'
    'env = "metawrap-1.3.0"\n'
    'version = "1.3.0"\n'
)


def _load(tmp_path, external=_VALID_EXTERNAL, variant_keys=""):
    path = tmp_path / "configurations.toml"
    path.write_text(_CONFIG % (variant_keys, external))
    return registry.load_configurations(str(path))


def _refusal(tmp_path, external=_VALID_EXTERNAL, variant_keys=""):
    with pytest.raises(RegistryError) as caught:
        _load(tmp_path, external, variant_keys)
    message = str(caught.value)
    # Every message names the arm.
    assert "'mw1'" in message
    return message


def test_a_valid_external_declaration_loads(tmp_path):
    """A valid external declaration loads; positive control for the refusals below."""
    comparison = _load(tmp_path)["binner_set"]
    arm = comparison.variant("mw1")
    assert arm.is_external
    assert arm.external.env == "metawrap-1.3.0"
    assert arm.external.version == "1.3.0"
    assert arm.external.placeholders == ("output", "threads", "assembly", "reads")
    # The native arm is unaffected.
    assert comparison.variant("metabat2").is_external is False
    assert comparison.incumbent.name == "metabat2"


def test_an_external_arm_with_nowhere_to_read_bins_from_is_refused(tmp_path):
    """A declaration without ``bins`` is refused, since the scorer would find no bins."""
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace('bins = "out/final_bins"\n', ""))
    assert "needs 'bins'" in message


def test_an_external_arm_without_an_env_is_refused(tmp_path):
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace('env = "metawrap-1.3.0"\n', ""))
    assert "needs 'env'" in message


def test_an_external_arm_without_a_command_is_refused(tmp_path):
    message = _refusal(
        tmp_path,
        "\n".join(line for line in _VALID_EXTERNAL.splitlines() if "command" not in line) + "\n",
    )
    assert "needs a 'command'" in message


def test_a_placeholder_the_runner_cannot_fill_is_refused_with_the_list_that_can(tmp_path):
    """An unknown placeholder is refused, and the message lists the known ones."""
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace("{assembly}", "{contigs}"))
    assert "{contigs}" in message and "{assembly}" in message


def test_an_arm_that_does_not_reference_its_modules_inputs_is_refused(tmp_path):
    """A binning command that omits ``{assembly}`` and ``{reads}`` is refused."""
    message = _refusal(
        tmp_path,
        _VALID_EXTERNAL.replace(
            "-a {assembly} {reads}", "-a /scratch/mine/gold.fa /scratch/mine/*.fq"
        ),
    )
    assert "{assembly}, {reads}" in message


def test_an_unversioned_external_arm_is_refused(tmp_path):
    """A declaration with neither ``version`` nor ``version_command`` is refused."""
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace('version = "1.3.0"\n', ""))
    assert "'version' or 'version_command'" in message


def test_a_command_that_returns_before_its_work_is_refused(tmp_path):
    """Submitted or backgrounded commands (sbatch, nohup, trailing ``&``) are refused at load."""
    for word, replacement in (
        ("sbatch", "sbatch metawrap binning"),
        ("nohup", "nohup metawrap binning"),
    ):
        message = _refusal(tmp_path, _VALID_EXTERNAL.replace("metawrap binning", replacement))
        assert word in message and "measure the submission" in message
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace("{reads}", "{reads} &"))
    assert "measure the submission" in message


def test_bins_must_not_escape_the_arms_output_directory(tmp_path):
    message = _refusal(tmp_path, _VALID_EXTERNAL.replace('"out/final_bins"', '"../../elsewhere"'))
    assert "must be relative to {output}" in message


def test_bins_named_like_a_refined_set_at_the_top_level_is_refused(tmp_path):
    """A top-level ``metawrap_<c>_<x>_bins`` path is refused: its bins would be scored twice."""
    message = _refusal(
        tmp_path, _VALID_EXTERNAL.replace('"out/final_bins"', '"metawrap_70_10_bins"')
    )
    assert "scored twice" in message


def test_an_arm_cannot_be_both_metawrap2_and_foreign(tmp_path):
    message = _refusal(tmp_path, variant_keys='args = ["--metabat2"]')
    assert "never both" in message


def test_an_external_arm_cannot_be_the_incumbent(tmp_path):
    """An external arm with ``default = true`` is refused."""
    message = _refusal(tmp_path, variant_keys="default = true")
    assert "cannot be 'default = true'" in message


def test_the_shipped_example_declaration_loads_when_uncommented(tmp_path):
    """The commented-out metaWRAP 1.3.0 example in configurations.toml loads once uncommented.

    The declaration starts at the example's `[[binner_set.variants]]` line.
    """
    with open(os.path.join(FIXTURES, "configurations.toml")) as fh:
        lines = fh.read().splitlines()
    start = [i for i, line in enumerate(lines) if line.startswith("# [[binner_set.variants]]")]
    assert len(start) == 1, "the example declaration moved or was duplicated"
    uncommented = "\n".join(
        line[2:] if line.startswith("# ") else ("" if line.strip() == "#" else line)
        for line in lines[start[0] :]
    )
    path = tmp_path / "configurations.toml"
    path.write_text(
        'schema = 1\n\n[binner_set]\nmodule = "binning"\nmetrics = ["bins_high_quality"]\n\n'
        '[[binner_set.variants]]\nname = "metabat2"\nargs = ["--metabat2"]\ndefault = true\n\n'
        + uncommented
    )

    arm = registry.load_configurations(str(path))["binner_set"].variant("mw1")
    assert arm.is_external
    assert arm.external.env == "metawrap-1.3.0"
    assert arm.external.version == "1.3.0"
    assert arm.external.bins == "refine/metawrap_70_10_bins"
    # One `bash -c`, so the two metaWRAP steps are one measured process tree.
    assert arm.external.command.startswith("bash -c")
    assert arm.external.placeholders == ("output", "threads", "assembly", "reads")


def test_the_shipped_registry_has_no_live_external_arm():
    """The shipped registry declares no live external arm.

    metaWRAP 1 needs an env this repository does not build.
    """
    comparisons = registry.load_configurations(os.path.join(FIXTURES, "configurations.toml"))
    live = [
        "%s.%s" % (name, variant.name)
        for name, comparison in comparisons.items()
        for variant in comparison.variants
        if variant.is_external
    ]
    assert not live, "declared external arm(s) in the shipped registry: %s" % ", ".join(live)


# ── measured, published, recorded ─────────────────────────────────────────────────────────────


def _fake_run(work, bins=("bin.1.fa", "bin.2.fa")):
    """A run_measured that writes an external arm's bins where its declaration says they land."""
    calls = []

    def fake(argv, **_kw):
        calls.append(argv)
        if "metawrap2.cli" not in argv:  # the external arm, not the native one
            target = work / "binner_set" / "mw1" / "external" / "out" / "final_bins"
            target.mkdir(parents=True, exist_ok=True)
            for name in bins:
                (target / name).write_text(">contig_1\nACGTACGTAC\n")
        return 0, {
            "wall_seconds": 100.0,
            "user_seconds": 3400.0,
            "system_seconds": 200.0,
            "peak_rss_bytes": 4096,
        }

    return fake, calls


def test_publishing_rebuilds_its_own_links_and_refuses_to_delete_anything_else(tmp_path):
    """Re-publishing replaces the arm's own links, and refuses a directory holding real files."""
    arm = tmp_path / "mw1"
    source = arm / "external" / "out" / "final_bins"
    source.mkdir(parents=True)
    (source / "bin.1.fa").write_text(">c\nACGT\n")
    variant = _external_variant()

    first = runner_mod.publish_external_bins(variant, str(arm))
    assert runner_mod.publish_external_bins(variant, str(arm)) == first
    assert sorted(os.listdir(first)) == ["bin.1.fasta"]

    (tmp_path / "mw1" / "mw1_bins" / "someones_data.fasta").write_text(">c\nACGT\n")
    with pytest.raises(runner_mod.RunnerError) as caught:
        runner_mod.publish_external_bins(variant, str(arm))
    assert "someones_data.fasta" in str(caught.value)
    assert (tmp_path / "mw1" / "mw1_bins" / "someones_data.fasta").is_file()


def test_an_external_arm_is_measured_and_recorded_like_a_native_one(tmp_path, monkeypatch):
    """An external arm gets the same cost rows as a native arm, and the run record has its version.

    cpu_hours comes from user+system time, not wall time.
    """
    work = tmp_path / "work"
    fake, calls = _fake_run(work)
    monkeypatch.setattr(runner_mod._usage, "run_measured", fake)

    rows = runner_mod.run_comparison(
        _tiny(),
        _mixed_comparison(),
        os.path.join(FIXTURES, "fixtures", "tiny"),
        str(work),
        say=lambda _m: None,
    )

    assert len(calls) == 2
    by_arm = {}
    for row in rows:
        by_arm.setdefault(row.variant, {})[row.metric] = row.value
    assert set(by_arm) == {"metabat2", "mw1"}
    # Same cost fields and arithmetic: 3600 CPU-seconds is one CPU-hour.
    assert by_arm["mw1"] == by_arm["metabat2"]
    assert by_arm["mw1"] == {"wall_seconds": 100.0, "cpu_hours": 1.0, "peak_rss_bytes": 4096.0}

    with open(work / "binner_set.run.json") as fh:
        record = json.load(fh)
    arm = record["variants"]["mw1"]
    assert arm["returncode"] == 0
    assert arm["external"] is True
    assert arm["external_version"] == "1.3.0"
    assert arm["external_env"] == "metawrap-1.3.0"
    assert arm["external_command_template"] == MW1_COMMAND
    # A native arm has no "external" key, so it reads the same as an older record.
    assert "external" not in record["variants"]["metabat2"]
    # Both arms ran on the dataset's one sample, recorded once for the comparison.
    assert record["samples"] == ["sample_1"]


def test_a_foreign_bin_set_reaches_the_scoring_path_unchanged(tmp_path, monkeypatch):
    """An external arm's bins are published as ``<arm>_bins`` with ``.fasta`` links.

    ``bin_dirs`` then finds exactly one bin set for the arm.
    """
    work = tmp_path / "work"
    fake, _calls = _fake_run(work)
    monkeypatch.setattr(runner_mod._usage, "run_measured", fake)

    runner_mod.run_comparison(
        _tiny(),
        _mixed_comparison(),
        os.path.join(FIXTURES, "fixtures", "tiny"),
        str(work),
        variants=["mw1"],
        say=lambda _m: None,
    )

    arm_dir = work / "binner_set" / "mw1"
    assert assess.bin_dirs(str(arm_dir)) == [str(arm_dir / "mw1_bins")]
    published = sorted(os.listdir(arm_dir / "mw1_bins"))
    assert published == ["bin.1.fasta", "bin.2.fasta"]
    link = arm_dir / "mw1_bins" / "bin.1.fasta"
    # Relative, so the links still resolve after the work directory is moved.
    assert not os.path.isabs(os.readlink(link))
    assert link.read_text().startswith(">contig_1")
    # The pipeline's own output is unchanged.
    assert (arm_dir / "external" / "out" / "final_bins" / "bin.1.fa").is_file()


def test_an_external_arm_that_wrote_no_bins_is_a_failed_arm_not_an_empty_row(tmp_path, monkeypatch):
    """An arm that exits 0 with no bins is recorded as failed and gives no rows.

    The other arms' rows and the run record are kept.
    """
    work = tmp_path / "work"
    fake, _calls = _fake_run(work, bins=())
    monkeypatch.setattr(runner_mod._usage, "run_measured", fake)

    rows = runner_mod.run_comparison(
        _tiny(),
        _mixed_comparison(),
        os.path.join(FIXTURES, "fixtures", "tiny"),
        str(work),
        say=lambda _m: None,
    )

    assert {row.variant for row in rows} == {"metabat2"}
    with open(work / "binner_set.run.json") as fh:
        record = json.load(fh)
    assert record["variants"]["mw1"]["returncode"] == 1
    assert "final_bins" in record["variants"]["mw1"]["error"]


def test_a_version_that_cannot_be_read_stops_the_comparison_before_it_starts(tmp_path, monkeypatch):
    """``version_command`` is run in the env before any arm starts; no answer stops the run."""
    answers = {"metawrap --version": "metaWRAP v=1.3.0\n"}
    monkeypatch.setattr(
        runner_mod, "_capture_in_env", lambda _env, argv: answers.get(" ".join(argv), "")
    )

    probed = _external_variant(version="", version_command="metawrap --version")
    assert runner_mod.external_version(probed) == "metaWRAP v=1.3.0"
    # A dry run does not probe, since the env may not be installed yet.
    assert runner_mod.external_version(probed, probe=False) == ""

    silent = _external_variant(version="", version_command="metawrap --nope")
    with pytest.raises(runner_mod.RunnerError) as caught:
        runner_mod.external_version(silent)
    assert "'mw1'" in str(caught.value) and "not reproducible" in str(caught.value)

    def explode(*_a, **_k):
        raise AssertionError("started an arm whose version is unknown")

    monkeypatch.setattr(runner_mod._usage, "run_measured", explode)
    comparison = Comparison(
        name="binner_set",
        module="binning",
        variants=(Variant(name="metabat2", args=("--metabat2",)), silent),
    )
    with pytest.raises(runner_mod.RunnerError):
        runner_mod.run_comparison(
            _tiny(),
            comparison,
            os.path.join(FIXTURES, "fixtures", "tiny"),
            str(tmp_path / "work"),
            say=lambda _m: None,
        )
