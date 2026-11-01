"""Datasets of plain (non-archive) read files, such as `scn_bioreactor`."""

from __future__ import annotations

import os

from metawrap2.benchmark import fetch, registry, runner


def _scn() -> registry.Dataset:
    return registry.load_datasets()["scn_bioreactor"]


def test_scn_bioreactor_is_registered_with_pinned_checksums():
    dataset = _scn()
    assert dataset.reads_layout == "paired"
    assert len(dataset.files) == 6
    assert all(remote.md5 and remote.size for remote in dataset.files)
    assert not any(remote.is_archive for remote in dataset.files)
    assert "trimmer" in registry.load_configurations()


def test_plain_data_files_land_where_run_reads_them(tmp_path):
    dataset = _scn()
    landed = fetch._landing_path(dataset, dataset.files[0], str(tmp_path))
    assert landed == os.path.join(str(tmp_path), "extracted", "SRR3901694_1.fastq.gz")
    cami = registry.load_datasets()["cami2_gi"]
    assert os.path.dirname(fetch._landing_path(cami, cami.files[0], "D")) == "D/archives"
    assert os.path.dirname(fetch._landing_path(cami, cami.gold_standard[0], "D")) == "D/gold"


def test_plain_data_files_count_as_present_once_downloaded(tmp_path):
    dataset = _scn()
    for remote in dataset.files:
        path = fetch._landing_path(dataset, remote, str(tmp_path))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "wb") as fh:
            fh.truncate(remote.size)
    assert fetch.plan(dataset, str(tmp_path)).complete


def test_sample_labels_come_from_file_names_and_match_samples(tmp_path):
    """The labels `benchmark run` gives the reads are the accessions `samples` lists."""
    dataset = _scn()
    for remote in dataset.files:
        path = fetch._landing_path(dataset, remote, str(tmp_path))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").close()
    pairs = runner.prepare_reads(
        dataset, str(tmp_path), str(tmp_path / "work"), say=lambda _m: None
    )
    assert sorted(pair[0] for pair in pairs) == sorted(dataset.samples)
