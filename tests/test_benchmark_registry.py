"""The registry and the fetcher: what the benchmark says it will download, and what it does.

The registry's whole value is that it is trusted - a fetch reads it and commits an afternoon and a
hundred gigabytes to what it says. So these tests are about the shipped declarations being
verifiable and the extractor writing only what it was asked to, only where it was asked to.
"""

import io
import os
import tarfile

from metawrap2.benchmark import fetch, registry

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARKS = os.path.join(HERE, "..", "benchmarks")


def test_every_remote_file_in_the_shipped_registry_is_verifiable():
    """A remote file with no checksum is a file nobody can tell is intact.

    This is the guard on the registry itself rather than on the loader: the loader would reject a
    missing md5, but only for a dataset someone tried to load. This asserts it for all of them, and
    that every shipped comparison has exactly one incumbent - a comparison with none cannot show a
    challenger what it has to beat - and names a module that exists.
    """
    from metawrap2.cli import MODULES

    datasets = registry.load_datasets(os.path.join(BENCHMARKS, "datasets.toml"))
    assert "cami2_gi" in datasets and "tiny" in datasets
    for dataset in datasets.values():
        if dataset.is_local:
            continue
        for remote in list(dataset.files) + list(dataset.gold_standard):
            assert remote.md5, "%s/%s has no md5" % (dataset.name, remote.path)
            assert remote.size and remote.size > 0, "%s/%s has no size" % (
                dataset.name,
                remote.path,
            )

    comparisons = registry.load_configurations(os.path.join(BENCHMARKS, "configurations.toml"))
    assert comparisons
    for comparison in comparisons.values():
        defaults = [v for v in comparison.variants if v.default]
        assert len(defaults) == 1, "%s has %d defaults" % (comparison.name, len(defaults))
        assert comparison.metrics, "%s declares no metrics to rank on" % comparison.name
        assert comparison.module in MODULES, "%s names %r" % (comparison.name, comparison.module)


def test_extraction_keeps_only_the_members_the_registry_asks_for(tmp_path):
    """~95% of a CAMI2 archive is BAMs no binning benchmark reads.

    And these are archives from the internet, where extractall writes wherever a member says to, so a
    member that escapes the destination is refused. The escaping archive is built by hand rather than
    with tarfile.add, which normalises an escaping name away before it ever reaches the archive.
    """
    staging = tmp_path / "staging"
    archive = tmp_path / "a.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        for name, body in (
            ("run_sample_0/reads/anonymous_reads.fq.gz", b"reads"),
            ("run_sample_0/bam/OTU_97.1.bam", b"x" * 5000),
            ("run_sample_0/contigs/anonymous_gsa.fasta.gz", b"contigs"),
        ):
            member_path = staging / name
            member_path.parent.mkdir(parents=True, exist_ok=True)
            member_path.write_bytes(body)
            tf.add(str(member_path), arcname=name)
    into = tmp_path / "out"
    members, written = fetch.extract(str(archive), str(into), keep=["*/reads/*", "*/contigs/*"])
    assert members == 2
    assert written == len(b"reads") + len(b"contigs")
    assert os.path.isfile(into / "run_sample_0" / "reads" / "anonymous_reads.fq.gz")
    assert not os.path.exists(into / "run_sample_0" / "bam" / "OTU_97.1.bam")

    evil = tmp_path / "evil.tar"
    with tarfile.open(evil, "w") as tf:
        for name in ("../escaped.txt", "/absolute.txt", "a/keep.txt"):
            info = tarfile.TarInfo(name)
            info.size = 2
            tf.addfile(info, io.BytesIO(b"no"))
    into = tmp_path / "out2"
    members, _written = fetch.extract(str(evil), str(into))
    assert members == 1
    assert os.path.isfile(into / "a" / "keep.txt")
    assert not os.path.exists(tmp_path / "escaped.txt")
