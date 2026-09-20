"""Regression tests for bugs found reviewing the initial MetaWrap2 import.

Each test here pins down behaviour that was previously wrong, so the fix cannot silently
regress. They are grouped by the bug they cover, with the original symptom noted.

This is a curated subset: one strong test per distinct past bug / invariant, rather than the
full set of near-duplicate permutations that all re-assert the same code path.
"""

import os
import random
import sys
from typing import ClassVar

import pytest

from metawrap2 import checkm
from metawrap2.modules import (
    binning,
    kraken2,
    quant_bins,
    read_qc,
)

# --- MetaWrap2's own helpers must run in the host interpreter ---------------------------
# The per-module conda envs contain the module's external tools only - no python, no
# metawrap2, no matplotlib. Invoking "python -m metawrap2..." inside them always failed.




# --- read_qc must never move the user's input files -------------------------------------


def test_finalize_copies_rather_than_moves_raw_input(tmp_path, monkeypatch):
    """--skip-trimming --skip-bmtagger used to os.replace() the user's raw FASTQ away."""
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


# --- kraken2 --subsample must be streaming, not load the whole library ------------------


def test_subsample_is_bounded_and_keeps_pairs_together(tmp_path):
    r1 = tmp_path / "a_1.fastq"
    r2 = tmp_path / "a_2.fastq"
    r1.write_text("".join("@r%d\nACGT\n+\nIIII\n" % i for i in range(100)))
    r2.write_text("".join("@r%d\nTGCA\n+\nIIII\n" % i for i in range(100)))
    out1, out2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"

    random.seed(0)
    kraken2._subsample_pairs(str(r1), str(r2), 10, str(out1), str(out2))
    names1 = [ln[1:].strip() for ln in out1.read_text().splitlines()[::4]]
    names2 = [ln[1:].strip() for ln in out2.read_text().splitlines()[::4]]
    assert len(names1) == 10
    assert names1 == names2  # mates stayed paired
    assert len(set(names1)) == 10  # no record sampled twice


# --- CheckM2 produces the same .stats contract as CheckM1 -------------------------------


def test_checkm2_summary_matches_stats_contract(tmp_path):
    report = tmp_path / "quality_report.tsv"
    report.write_text(
        "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tContig_N50\t"
        "Genome_Size\tGC_Content\n"
        "bin.1\t95.42\t1.23\tneural\t45000\t3800000\t0.512\n"
        "bin.2\t61.00\t8.90\tneural\t12000\t2100000\t0.401\n"
    )
    rows = checkm.summarize_checkm2(str(report))
    assert [r[0] for r in rows] == ["bin.1", "bin.2"]
    # Same column order/arity as the CheckM1 path, so every consumer works unchanged.
    assert len(rows[0]) == len(checkm.STATS_HEADER)
    assert rows[0][1:3] == ["95.42", "1.23"]
    assert rows[0][4] == "checkm2"
    assert rows[0][5:7] == ["45000", "3800000"]

    out = tmp_path / "bins.stats"
    checkm.write_stats(rows, str(out))
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == checkm.STATS_HEADER
    assert lines[1].startswith("bin.1")  # sorted by completeness, descending


# --- the BLAST database must not be assumed to be called "nt" with 2-digit volumes -------


def _settings(**dbs):
    from metawrap2.config import Settings

    return Settings(databases=dbs)


@pytest.mark.parametrize(
    "filename",
    [
        "nt.nsq",  # single volume, no numbering
        "nt.00.nsq",  # legacy 2-digit volume
        "nt.000.nsq",  # BLAST v5 3-digit volume
        "nt.nal",  # alias file standing in for volumes
    ],
)
def test_blastdb_is_found_for_every_real_layout(tmp_path, filename):
    from metawrap2 import blastdb

    (tmp_path / filename).write_text("")
    ok, problem = blastdb.find(_settings(BLASTDB=str(tmp_path)))
    assert ok, problem


# --- install-db writes config keys without destroying the user's file --------------------


def test_install_db_updates_existing_key_and_keeps_comments(tmp_path):
    from metawrap2.commands import install_db

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


# --- --dry-run must print a module's WHOLE pipeline, not just the first command ----------
# Previously the orchestration consumed outputs that a dry run never produced, so every
# module died with a raw traceback after one or two commands. A representative spread of
# modules (a simple one, the widest fan-out, and a bin-set consumer) pins the plumbing.

DRY_RUN_CASES = [
    ("read_qc", ["-1", "{r1}", "-2", "{r2}", "-o", "{out}"], 3),
    (
        "binning",
        ["-a", "{asm}", "-o", "{out}", "--metabat2", "--maxbin2", "--concoct", "{r1}", "{r2}"],
        8,
    ),
    ("bin_refinement", ["-o", "{out}", "-A", "{bins}", "-B", "{bins2}"], 3),
]


@pytest.fixture
def dry_inputs(tmp_path, monkeypatch):
    """Minimal inputs plus a config that satisfies every module's database checks."""
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">NODE_1_length_5000_cov_3.5\n" + "ACGT" * 1250 + "\n")
    bins2 = tmp_path / "bins2"
    bins2.mkdir()
    (bins2 / "bin.1.fa").write_text(">NODE_1_length_5000_cov_3.5\n" + "ACGT" * 1250 + "\n")
    asm = tmp_path / "asm.fa"
    asm.write_text(">NODE_1_length_5000_cov_3.5\n" + "ACGT" * 1250 + "\n")
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@r1\nACGT\n+\nIIII\n")
    r2.write_text("@r1\nTGCA\n+\nIIII\n")

    # Stand-in databases: the modules only check that these paths/files exist.
    dbs = tmp_path / "dbs"
    (dbs / "KRAKEN2").mkdir(parents=True)
    (dbs / "BLAST").mkdir(parents=True)
    (dbs / "BLAST" / "nt.nsq").write_text("")
    (dbs / "TAX").mkdir(parents=True)
    for name in ("names.dmp", "nodes.dmp", "citations.dmp"):
        (dbs / "TAX" / name).write_text("")
    (dbs / "HOST").mkdir(parents=True)
    (dbs / "HOST" / "hg38.bitmask").write_text("")

    cfg = tmp_path / "config.toml"
    cfg.write_text(
        "[settings]\nthreads = 4\nuse_conda_envs = true\n\n[databases]\n"
        'KRAKEN2_DB = "%s"\nBMTAGGER_DB = "%s"\nBLASTDB = "%s"\nTAXDUMP = "%s"\n'
        % (dbs / "KRAKEN2", dbs / "HOST", dbs / "BLAST", dbs / "TAX")
    )

    # Every module env "exists" so env_for does not abort; nothing is executed anyway.
    from metawrap2.modules import _common

    monkeypatch.setattr(_common, "conda_env_exists", lambda env: True)

    return {
        "bins": str(bins),
        "bins2": str(bins2),
        "asm": str(asm),
        "r1": str(r1),
        "r2": str(r2),
        "out": str(tmp_path / "out"),
        "config": str(cfg),
    }


@pytest.mark.parametrize("module_name,argv_template,min_commands", DRY_RUN_CASES)
def test_dry_run_reaches_the_end_of_every_module(
    module_name, argv_template, min_commands, dry_inputs, capsys, monkeypatch
):
    import importlib

    from metawrap2 import command as command_mod

    recorded = []

    class Recorder:
        def record_command(self, cmd):
            recorded.append(cmd)

    module = importlib.import_module("metawrap2.modules.%s" % module_name)
    monkeypatch.setattr(command_mod.runner, "dry_run", True)
    monkeypatch.setattr(command_mod.runner, "force", True)

    argv = [a.format(**dry_inputs) for a in argv_template]
    argv += ["--config", dry_inputs["config"], "-o", dry_inputs["out"] + "-" + module_name]
    # de-duplicate the -o we may have added twice
    if argv.count("-o") > 1:
        first = argv.index("-o")
        del argv[first : first + 2]

    try:
        rc = module.main(argv)
    finally:
        command_mod.set_recorder(None)

    out = capsys.readouterr()
    assert "Traceback" not in out.out + out.err
    assert rc == 0, out.err[-2000:]
    printed = [ln for ln in out.err.splitlines() if ln.startswith("[dry-run]")]
    assert len(printed) >= min_commands, "%s only got to %d command(s) under --dry-run:\n%s" % (
        module_name,
        len(printed),
        "\n".join(printed),
    )


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

    binning_mod._run_metabat2(Args(), None, str(tmp_path / "work"), str(tmp_path / "asm.fa"))

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


# --- log noise filtering: command data output must never be filtered ---------------------




# --- relative -o must work even for tools that run with cwd=<output dir> -----------------


