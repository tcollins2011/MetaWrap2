"""Regression tests for bugs found reviewing the initial MetaWrap2 import.

Each test here pins down behaviour that was previously wrong, so the fix cannot silently
regress. They are grouped by the bug they cover, with the original symptom noted.

This is a curated subset: one strong test per distinct past bug / invariant, rather than the
full set of near-duplicate permutations that all re-assert the same code path.
"""

import os

from metawrap2.modules import _common, binning, read_qc

# --- read_qc must never move the user's input files -------------------------------------


def test_finalize_copies_rather_than_moves_raw_input(tmp_path, monkeypatch):
    """--skip-trimming --skip-bmtagger used to os.replace() the user's raw FASTQ away."""
    # Every stage is skipped, so no tool runs - but env_for checks the env manager is findable
    # before the module does any work, and mamba is not on PATH in the shell pytest runs in.
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)

    raw1 = tmp_path / "sample_1.fastq"
    raw2 = tmp_path / "sample_2.fastq"
    raw1.write_text("@r1\nACGT\n+\nIIII\n")
    raw2.write_text("@r1\nTGCA\n+\nIIII\n")
    out = tmp_path / "qc"

    rc = read_qc.main(
        [
            "-1",
            str(raw1),
            "-2",
            str(raw2),
            "-o",
            str(out),
            "--skip-trimming",
            "--skip-bmtagger",
            "--skip-pre-qc-report",
            "--skip-post-qc-report",
        ]
    )
    assert rc == 0
    # The originals must still be there, untouched.
    assert raw1.is_file() and raw1.read_text() == "@r1\nACGT\n+\nIIII\n"
    assert raw2.is_file()
    assert (out / "final_pure_reads_1.fastq").is_file()


# --- CONCOCT bin splitting must not truncate contig names at the first dot --------------


def test_concoct_split_keeps_dotted_contig_names(tmp_path):
    """metaSPAdes names contain a dot (cov_3.5); splitting on the first one merged bins."""
    assembly = tmp_path / "asm.fa"
    assembly.write_text(
        ">NODE_1_length_5000_cov_3.5\nAAAA\n"
        ">NODE_2_length_4000_cov_3.9\nCCCC\n"
        ">NODE_3_length_3000_cov_7.1\nGGGG\n"
    )
    clustering = tmp_path / "clust.csv"
    clustering.write_text(
        "contig_id,cluster_id\n" "NODE_1_length_5000_cov_3.5,0\n" "NODE_2_length_4000_cov_3.9,1\n"
    )
    out = tmp_path / "concoct_bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))

    assert (out / "bin_001.fasta").read_text() == ">NODE_1_length_5000_cov_3.5\nAAAA\n"
    assert (out / "bin_002.fasta").read_text() == ">NODE_2_length_4000_cov_3.9\nCCCC\n"
    # The contig absent from the clustering is unbinned, not silently dropped.
    assert ">NODE_3_length_3000_cov_7.1" in (out / "unbinned.fasta").read_text()


# --- install-db: config keys and archives off the internet ------------------------------


def test_install_db_untar_refuses_a_member_that_escapes_the_destination(tmp_path):
    """These are archives off the public internet, and extractall writes where a member says to."""
    import io
    import tarfile

    from metawrap2.commands import install_db

    archive = tmp_path / "db.tar"
    # Built by hand: tarfile.add normalises an escaping name away before it reaches the archive.
    with tarfile.open(archive, "w") as tf:
        for name in ("../escaped.dmp", "/absolute.dmp", "names.dmp"):
            info = tarfile.TarInfo(name)
            info.size = 2
            tf.addfile(info, io.BytesIO(b"ok"))
    into = tmp_path / "db"
    install_db._untar(str(archive), str(into))
    assert os.path.isfile(into / "names.dmp")
    assert not os.path.exists(tmp_path / "escaped.dmp")
    assert not os.path.exists(into / "absolute.dmp")
    assert not os.path.exists(archive)  # unpacked archives are still removed

    # Writing the installed path back updates the existing key without destroying the user's file.
    cfg = tmp_path / "config.toml"
    cfg.write_text(
        "# my notes\n[databases]\n# where kraken lives\n"
        'KRAKEN2_DB = "/old/path"\nTAXDUMP = "/tax"\n'
    )
    install_db._set_config_keys(str(cfg), {"KRAKEN2_DB": "/new/path"})
    text = cfg.read_text()
    assert 'KRAKEN2_DB = "/new/path"' in text
    assert "/old/path" not in text
    assert 'TAXDUMP = "/tax"' in text  # other keys survive
    assert "# my notes" in text and "# where kraken lives" in text


# --- dereplication must survive a bin with no stats row ---------------------------------


def test_dereplicate_tolerates_a_bin_missing_from_stats(tmp_path, capsys):
    """This used to raise KeyError and abort the whole refinement."""
    from metawrap2 import refinement

    bins = tmp_path / "binsM"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">shared\nAAAA\n>only1\nCCCC\n")
    (bins / "bin.2.fa").write_text(">shared\nAAAA\n>only2\nGGGG\n")  # no stats row
    stats = tmp_path / "binsM.stats"
    stats.write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin.1\t90\t1\t0.5\tBacteria\t10000\t3000000\n"
    )

    out = tmp_path / "binsO"
    refinement.dereplicate(str(stats), str(bins), str(out), mode="best")

    # The scored bin keeps the contested contig; the unscored bin keeps its own only.
    assert ">shared" in (out / "bin.1.fa").read_text()
    assert ">shared" not in (out / "bin.2.fa").read_text()
    assert ">only2" in (out / "bin.2.fa").read_text()
    assert "no row in" in capsys.readouterr().err


# --- metaBAT2 writes non-genome files into its output prefix -----------------------------


def test_metabat2_non_bin_outputs_are_moved_out_of_the_bin_set(tmp_path, monkeypatch):
    """lowDepth/tooShort/unbinned/summaries are not genomes; CheckM used to score them."""
    from metawrap2.modules import binning as binning_mod

    out = tmp_path / "out"
    bins_dir = out / "metabat2_bins"
    bins_dir.mkdir(parents=True)
    for name in (
        "bin.1.fa",
        "bin.2.fa",
        "bin.unbinned.fa",
        "bin.lowDepth.fa",
        "bin.tooShort.fa",
        "bin.BinInfo.txt",
        "bin.BinMembers.txt",
    ):
        (bins_dir / name).write_text(">c\nACGT\n")

    calls = []
    monkeypatch.setattr(binning_mod, "run", lambda *a, **k: calls.append(a))
    monkeypatch.setattr(binning_mod, "dry_run", lambda: False)

    class Args:
        output = str(out)
        threads = 1
        min_len = 1000

    binning_mod._run_metabat2(
        Args(), None, str(tmp_path / "work"), str(tmp_path / "asm.fa"), "depth.txt"
    )

    remaining = sorted(os.listdir(bins_dir))
    # metaBAT2's own bin.1.fa/bin.2.fa are renumbered to the canonical naming once the
    # non-bins are out of the way, so a rejected-contigs file can never take a bin number.
    assert remaining == ["bin_001.fasta", "bin_002.fasta"], remaining
    moved = sorted(os.listdir(out / "metabat2_extra"))
    assert moved == [
        "bin.BinInfo.txt",
        "bin.BinMembers.txt",
        "bin.lowDepth.fa",
        "bin.tooShort.fa",
        "bin.unbinned.fa",
    ], moved
