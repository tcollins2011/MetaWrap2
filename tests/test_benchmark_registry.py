"""The registry and the fetcher: what the benchmark says it will download, and what it does.

The registry's whole value is that it is trusted - a fetch reads it and commits an afternoon and a
hundred gigabytes to what it says. So these tests are mostly about the loader refusing bad
declarations rather than accepting good ones, and about the fetcher refusing to start rather than
finishing.
"""

import io
import os
import tarfile

import pytest

from metawrap2 import limits
from metawrap2.benchmark import fetch, registry

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARKS = os.path.join(HERE, "..", "benchmarks")


def write(tmp_path, text, name="datasets.toml"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


GOOD = """
schema = 1
[demo]
title = "Demo"
base_url = "https://example.invalid/data/"
samples = [1]
[[demo.files]]
sample = 1
path = "sample_1.tar.gz"
md5 = "0123456789abcdef0123456789abcdef"
bytes = 1024
[[demo.gold_standard]]
role = "binning"
path = "truth.tsv"
md5 = "fedcba9876543210fedcba9876543210"
bytes = 16
"""


# ── the shipped registry ──────────────────────────────────────────────────────────────────────


def test_the_shipped_registry_loads():
    datasets = registry.load_datasets(os.path.join(BENCHMARKS, "datasets.toml"))
    assert "cami2_gi" in datasets
    assert "tiny" in datasets


def test_every_remote_file_in_the_shipped_registry_is_verifiable():
    """A remote file with no checksum is a file nobody can tell is intact.

    This is the guard on the registry itself rather than on the loader: the loader would reject a
    missing md5, but only for a dataset someone tried to load. This asserts it for all of them.
    """
    for dataset in registry.load_datasets(os.path.join(BENCHMARKS, "datasets.toml")).values():
        if dataset.is_local:
            continue
        for remote in list(dataset.files) + list(dataset.gold_standard):
            assert remote.md5, "%s/%s has no md5" % (dataset.name, remote.path)
            assert remote.size and remote.size > 0, "%s/%s has no size" % (
                dataset.name,
                remote.path,
            )


def test_the_cami2_download_total_is_the_sum_of_its_files():
    dataset = registry.get_dataset("cami2_gi", os.path.join(BENCHMARKS, "datasets.toml"))
    expected = sum(f.size for f in list(dataset.files) + list(dataset.gold_standard))
    assert dataset.download_bytes == expected
    # The ten gastrointestinal samples, not the twenty in the published directory.
    assert len(dataset.files) == 10
    assert list(dataset.samples) == [0, 1, 2, 3, 4, 5, 9, 10, 11, 12]


def test_the_committed_fixture_is_actually_committed():
    """The fixture exists so the harness is testable without a 95 GB download.

    Which means the files the registry promises have to be in the repository, not generated on
    demand: a fixture regenerated per run would make every results diff noise.
    """
    dataset = registry.get_dataset("tiny", os.path.join(BENCHMARKS, "datasets.toml"))
    for remote in list(dataset.files) + list(dataset.gold_standard):
        path = dataset.local_path_for(remote)
        assert os.path.isfile(path), "%s is missing" % path
        assert os.path.getsize(path) > 0


def test_the_shipped_comparisons_each_have_exactly_one_incumbent():
    """A comparison with no incumbent cannot show a challenger what it has to beat."""
    comparisons = registry.load_configurations(os.path.join(BENCHMARKS, "configurations.toml"))
    assert comparisons
    for comparison in comparisons.values():
        defaults = [v for v in comparison.variants if v.default]
        assert len(defaults) == 1, "%s has %d defaults" % (comparison.name, len(defaults))
        assert comparison.metrics, "%s declares no metrics to rank on" % comparison.name


def test_every_shipped_comparison_names_a_real_module():
    from metawrap2.cli import MODULES

    comparisons = registry.load_configurations(os.path.join(BENCHMARKS, "configurations.toml"))
    for comparison in comparisons.values():
        assert comparison.module in MODULES, "%s names %r" % (comparison.name, comparison.module)


# ── what the loader refuses ───────────────────────────────────────────────────────────────────


def test_a_remote_file_without_a_checksum_is_refused(tmp_path):
    text = GOOD.replace('md5 = "0123456789abcdef0123456789abcdef"\n', "")
    with pytest.raises(registry.RegistryError, match="no md5"):
        registry.load_datasets(write(tmp_path, text))


def test_a_remote_file_without_a_size_is_refused(tmp_path):
    """The size is what the disk pre-flight compares against; without it there is no pre-flight."""
    text = GOOD.replace("bytes = 1024\n", "")
    with pytest.raises(registry.RegistryError, match="no 'bytes'"):
        registry.load_datasets(write(tmp_path, text))


def test_an_escaping_path_is_refused(tmp_path):
    text = GOOD.replace('path = "sample_1.tar.gz"', 'path = "../../etc/passwd"')
    with pytest.raises(registry.RegistryError, match="must not contain"):
        registry.load_datasets(write(tmp_path, text))


def test_an_absolute_path_is_refused(tmp_path):
    text = GOOD.replace('path = "sample_1.tar.gz"', 'path = "/etc/passwd"')
    with pytest.raises(registry.RegistryError, match="must be relative"):
        registry.load_datasets(write(tmp_path, text))


def test_an_unknown_gold_standard_role_is_refused(tmp_path):
    """A typo in a role means scoring against the wrong file, silently."""
    text = GOOD.replace('role = "binning"', 'role = "bining"')
    with pytest.raises(registry.RegistryError, match="unknown gold-standard role"):
        registry.load_datasets(write(tmp_path, text))


def test_a_missing_schema_number_is_refused(tmp_path):
    with pytest.raises(registry.RegistryError, match="schema = 1"):
        registry.load_datasets(write(tmp_path, GOOD.replace("schema = 1", "")))


def test_a_sample_listed_but_never_declared_is_refused(tmp_path):
    """Otherwise `samples = [0..9]` with nine files quietly benchmarks nine samples."""
    text = GOOD.replace("samples = [1]", "samples = [1, 2]")
    with pytest.raises(registry.RegistryError, match="no file declares them"):
        registry.load_datasets(write(tmp_path, text))


def test_a_file_referencing_an_undeclared_sample_is_refused(tmp_path):
    text = GOOD.replace("sample = 1\npath", "sample = 7\npath")
    with pytest.raises(registry.RegistryError, match="does not list"):
        registry.load_datasets(write(tmp_path, text))


def test_a_dataset_that_is_both_local_and_remote_is_refused(tmp_path):
    text = GOOD.replace('title = "Demo"', 'title = "Demo"\nlocal = "somewhere"')
    with pytest.raises(registry.RegistryError, match="not both"):
        registry.load_datasets(write(tmp_path, text))


def test_an_unknown_reads_layout_is_refused(tmp_path):
    text = GOOD.replace('title = "Demo"', 'title = "Demo"\nreads_layout = "pared"')
    with pytest.raises(registry.RegistryError, match="unknown reads_layout"):
        registry.load_datasets(write(tmp_path, text))


COMPARISONS = """
schema = 1
[choice]
module = "binning"
[[choice.variants]]
name = "a"
default = true
[[choice.variants]]
name = "b"
"""


def test_a_comparison_with_one_variant_is_refused(tmp_path):
    text = COMPARISONS.replace('[[choice.variants]]\nname = "b"\n', "")
    with pytest.raises(registry.RegistryError, match="not comparing anything"):
        registry.load_configurations(write(tmp_path, text, "c.toml"))


def test_two_incumbents_are_refused(tmp_path):
    text = COMPARISONS.replace('name = "b"', 'name = "b"\ndefault = true')
    with pytest.raises(registry.RegistryError, match="Only one option can be the incumbent"):
        registry.load_configurations(write(tmp_path, text, "c.toml"))


def test_a_duplicated_variant_name_is_refused(tmp_path):
    """Two variants with one name overwrite each other's rows in the results TSV."""
    text = COMPARISONS.replace('name = "b"', 'name = "a"')
    with pytest.raises(registry.RegistryError, match="twice"):
        registry.load_configurations(write(tmp_path, text, "c.toml"))


def test_an_unknown_variant_is_named_with_the_known_ones(tmp_path):
    comparison = registry.load_configurations(write(tmp_path, COMPARISONS, "c.toml"))["choice"]
    with pytest.raises(registry.RegistryError, match="Known: a, b"):
        comparison.variant("nope")


# ── verification and planning ─────────────────────────────────────────────────────────────────


def _remote(path="f.bin", md5=None, size=None):
    return registry.Remote(path=path, md5=md5, size=size)


def test_verify_catches_a_truncated_download_by_size_alone(tmp_path):
    """The cheap half of verification, and the failure that actually happens."""
    path = tmp_path / "f.bin"
    path.write_bytes(b"short")
    reason = fetch.verify(str(path), _remote(size=1000, md5="x" * 32))
    assert "wrong size" in reason
    assert "1000 expected" in reason


def test_verify_catches_a_corrupt_download_by_checksum(tmp_path):
    path = tmp_path / "f.bin"
    path.write_bytes(b"hello")
    reason = fetch.verify(str(path), _remote(size=5, md5="0" * 32))
    assert "md5 mismatch" in reason


def test_verify_accepts_an_intact_file(tmp_path):
    path = tmp_path / "f.bin"
    path.write_bytes(b"hello")
    md5 = fetch.md5_of(str(path))
    assert fetch.verify(str(path), _remote(size=5, md5=md5)) == fetch.OK


def test_verify_reports_a_missing_file_as_missing_not_as_corrupt(tmp_path):
    assert fetch.verify(str(tmp_path / "nope"), _remote(size=5, md5="x" * 32)) == "missing"


def test_a_plan_distinguishes_never_downloaded_from_downloaded_and_broken(tmp_path):
    """Different problems calling for different responses: fetch it, versus look at the host."""
    dataset = registry.load_datasets(write(tmp_path, GOOD))["demo"]
    dest = tmp_path / "dest"
    (dest / "archives").mkdir(parents=True)
    (dest / "archives" / "sample_1.tar.gz").write_bytes(b"truncated")
    plan = fetch.plan(dataset, str(dest))
    reasons = {remote.path: reason for remote, reason in plan.needed}
    assert "wrong size" in reasons["sample_1.tar.gz"]
    assert reasons["truth.tsv"] == "missing"


def test_a_plan_over_the_disk_ceiling_is_refused(tmp_path, monkeypatch):
    dataset = registry.load_datasets(write(tmp_path, GOOD))["demo"]
    plan = fetch.plan(dataset, str(tmp_path))
    monkeypatch.setattr(
        fetch._limits,
        "disk_state",
        lambda path, budget: limits.DiskState(
            path=str(tmp_path), total_bytes=1000, used_bytes=900, ceiling_bytes=950
        ),
    )
    ok, message = fetch.disk_verdict(plan, limits.resolve())
    assert not ok
    assert "of headroom" in message


def test_a_plan_inside_the_disk_ceiling_is_allowed(tmp_path, monkeypatch):
    dataset = registry.load_datasets(write(tmp_path, GOOD))["demo"]
    plan = fetch.plan(dataset, str(tmp_path))
    monkeypatch.setattr(
        fetch._limits,
        "disk_state",
        lambda path, budget: limits.DiskState(
            path=str(tmp_path), total_bytes=10**9, used_bytes=0, ceiling_bytes=10**9
        ),
    )
    ok, _message = fetch.disk_verdict(plan, limits.resolve())
    assert ok


def test_nothing_is_downloaded_when_the_disk_check_fails(tmp_path, monkeypatch):
    """The point of the pre-flight: a refusal has to happen before the first byte."""
    dataset = registry.load_datasets(write(tmp_path, GOOD))["demo"]

    def explode(*_a, **_k):
        raise AssertionError("downloaded despite the disk check failing")

    monkeypatch.setattr(fetch, "_download", explode)
    monkeypatch.setattr(
        fetch._limits,
        "disk_state",
        lambda path, budget: limits.DiskState(
            path=str(tmp_path), total_bytes=1000, used_bytes=999, ceiling_bytes=999
        ),
    )
    said = []
    rc = fetch.fetch(dataset, str(tmp_path / "dest"), say=said.append)
    assert rc == 1
    assert any("disk ceiling" in line for line in said)


def test_a_dry_run_downloads_nothing(tmp_path, monkeypatch):
    dataset = registry.load_datasets(write(tmp_path, GOOD))["demo"]

    def explode(*_a, **_k):
        raise AssertionError("downloaded during a dry run")

    monkeypatch.setattr(fetch, "_download", explode)
    said = []
    rc = fetch.fetch(dataset, str(tmp_path / "dest"), dry_run=True, say=said.append)
    assert rc == 0
    assert any("--dry-run" in line for line in said)


# ── extraction ────────────────────────────────────────────────────────────────────────────────


def _archive(tmp_path, entries, links=()):
    path = tmp_path / "a.tar.gz"
    with tarfile.open(path, "w:gz") as tf:
        for name, body in entries:
            member_path = tmp_path / "staging" / name
            member_path.parent.mkdir(parents=True, exist_ok=True)
            member_path.write_bytes(body)
            tf.add(str(member_path), arcname=name)
        for name, target in links:
            info = tarfile.TarInfo(name)
            info.type = tarfile.SYMTYPE
            info.linkname = target
            tf.addfile(info)
    return str(path)


def test_extraction_keeps_only_the_members_the_registry_asks_for(tmp_path):
    """~95% of a CAMI2 archive is BAMs no binning benchmark reads."""
    archive = _archive(
        tmp_path,
        [
            ("run_sample_0/reads/anonymous_reads.fq.gz", b"reads"),
            ("run_sample_0/bam/OTU_97.1.bam", b"x" * 5000),
            ("run_sample_0/contigs/anonymous_gsa.fasta.gz", b"contigs"),
        ],
    )
    into = tmp_path / "out"
    members, written = fetch.extract(archive, str(into), keep=["*/reads/*", "*/contigs/*"])
    assert members == 2
    assert written == len(b"reads") + len(b"contigs")
    assert os.path.isfile(into / "run_sample_0" / "reads" / "anonymous_reads.fq.gz")
    assert not os.path.exists(into / "run_sample_0" / "bam" / "OTU_97.1.bam")


def test_extraction_keeps_everything_when_no_keep_globs_are_declared(tmp_path):
    archive = _archive(tmp_path, [("a/one.txt", b"1"), ("a/two.txt", b"22")])
    members, written = fetch.extract(archive, str(tmp_path / "out"))
    assert members == 2
    assert written == 3


def test_extraction_refuses_a_member_that_escapes_the_destination(tmp_path):
    """These are archives from the internet, and extractall writes wherever a member says to."""
    # Built by hand rather than with tarfile.add, which normalises an escaping name away before it
    # ever reaches the archive - and so would test nothing.
    evil = tmp_path / "evil.tar"
    with tarfile.open(evil, "w") as tf:
        for name in ("../escaped.txt", "/absolute.txt", "a/keep.txt"):
            info = tarfile.TarInfo(name)
            info.size = 2
            tf.addfile(info, io.BytesIO(b"no"))
    into = tmp_path / "out"
    members, _written = fetch.extract(str(evil), str(into))
    assert members == 1
    assert os.path.isfile(into / "a" / "keep.txt")
    assert not os.path.exists(tmp_path / "escaped.txt")


def test_extraction_skips_symlinks(tmp_path):
    """The other half of how an archive escapes its destination."""
    archive = _archive(tmp_path, [("a/keep.txt", b"ok")], links=[("a/link", "/etc/passwd")])
    into = tmp_path / "out"
    fetch.extract(archive, str(into))
    assert os.path.isfile(into / "a" / "keep.txt")
    assert not os.path.lexists(into / "a" / "link")
