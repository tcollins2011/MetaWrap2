"""Tests for the benchmark registry and fetcher.

A fetch can take hours and a hundred gigabytes, so the shipped declarations must be verifiable, and
the extractor must write only the requested members, only inside the destination.
"""

import io
import os
import tarfile

from metawrap2.benchmark import fetch, registry

HERE = os.path.dirname(os.path.abspath(__file__))
BENCHMARKS = os.path.join(HERE, "..", "benchmarks")


def test_every_remote_file_in_the_shipped_registry_is_verifiable():
    """Every shipped remote file has an md5 and size, and every comparison is well formed.

    The loader rejects a missing md5 only for the dataset being loaded; this checks all of them.
    Each comparison must have exactly one default variant, declare metrics, and name an existing
    module.
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
    """Extraction keeps only matching members and refuses members that escape the destination.

    About 95% of a CAMI2 archive is BAMs that no binning benchmark reads. The escaping archive is
    built with addfile, because tarfile.add normalises an escaping name away.
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
