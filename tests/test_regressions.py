"""Regression tests for bugs found reviewing the initial MetaWrap2 import.

Each test here pins down behaviour that was previously wrong, so the fix cannot silently
regress. They are grouped by the bug they cover, with the original symptom noted.
"""

import os
import random
import sys
from typing import ClassVar

import pytest

from metawrap2 import checkm
from metawrap2.modules import (
    annotate_bins,
    bin_refinement,
    binning,
    blobology,
    kraken2,
    quant_bins,
    read_qc,
    reassemble_bins,
)

# --- MetaWrap2's own helpers must run in the host interpreter ---------------------------
# The per-module conda envs contain the module's external tools only - no python, no
# metawrap2, no matplotlib. Invoking "python -m metawrap2..." inside them always failed.

HELPER_TEMPLATES = [
    bin_refinement.BINNING_REFINER,
    bin_refinement.BINNING_REFINER_3,
    bin_refinement.PLOT_BINNING,
    reassemble_bins.RM_SHORT_CONTIGS,
    reassemble_bins.CHOOSE_BEST_BIN,
    reassemble_bins.PLOT_REASSEMBLY,
    reassemble_bins.BWA_MEM_FILTER,
    reassemble_bins.MINIMAP2_FILTER,
    blobology.FASTAQUAL_SELECT,
    blobology.GC_COV_ANNOTATE,
    blobology.ADD_BINS,
    blobology.MAKEBLOBPLOT,
    blobology.MAKEBLOBPLOT_BASE,
]


@pytest.mark.parametrize("template", HELPER_TEMPLATES)
def test_helpers_use_host_interpreter_by_absolute_path(template):
    assert "metawrap2." in template
    # Never a bare "python": inside `mamba run -n env` that resolves to the env's python
    # (or nothing at all), which cannot import metawrap2.
    assert " python -m " not in template
    assert not template.startswith("python ")
    assert sys.executable in template


def test_make_heatmap_argv_is_host_interpreter():
    assert quant_bins.MAKE_HEATMAP[0] == sys.executable
    assert quant_bins.MAKE_HEATMAP[1:3] == ["-m", "metawrap2.scripts.make_heatmap"]


# --- read_qc must never move the user's input files -------------------------------------


def test_is_ours_only_claims_paths_inside_the_output_dir(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    assert read_qc._is_ours(str(out / "trimmed_1.fastq"), str(out))
    assert not read_qc._is_ours(str(tmp_path / "raw_1.fastq"), str(out))
    # A sibling directory sharing a name prefix must not be mistaken for being inside.
    assert not read_qc._is_ours(str(tmp_path / "out2" / "x.fastq"), str(out))


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


def test_concoct_split_is_idempotent(tmp_path):
    """Appending per contig meant a re-run duplicated every sequence."""
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nAAAA\n")
    clustering = tmp_path / "clust.csv"
    clustering.write_text("contig_id,cluster_id\nc1,0\n")
    out = tmp_path / "bins"
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    first = (out / "bin_001.fasta").read_text()
    binning._split_concoct_bins(str(clustering), str(assembly), str(out))
    assert (out / "bin_001.fasta").read_text() == first


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


def test_subsample_of_more_than_available_keeps_everything(tmp_path):
    r1 = tmp_path / "a_1.fastq"
    r2 = tmp_path / "a_2.fastq"
    r1.write_text("@r0\nACGT\n+\nIIII\n")
    r2.write_text("@r0\nTGCA\n+\nIIII\n")
    out1, out2 = tmp_path / "s_1.fastq", tmp_path / "s_2.fastq"
    kraken2._subsample_pairs(str(r1), str(r2), 50, str(out1), str(out2))
    assert out1.read_text().count("@r0") == 1


def test_truncated_fastq_is_reported(tmp_path):
    bad = tmp_path / "bad.fastq"
    bad.write_text("@r0\nACGT\n+\n")  # missing the quality line
    with pytest.raises(ValueError, match="truncated"):
        list(kraken2._iter_records(str(bad)))


# --- annotate_bins must only annotate fasta files ---------------------------------------


def test_annotate_bins_rejects_a_folder_with_no_fasta(tmp_path, capsys):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bins.stats").write_text("bin\tcompleteness\n")  # not a genome
    (bins / "bins.contigs").write_text("c1\tbin.1\n")
    with pytest.raises(SystemExit):
        annotate_bins.main(["-b", str(bins), "-o", str(tmp_path / "out")])


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


def test_checkm2_rejects_unexpected_columns(tmp_path):
    report = tmp_path / "quality_report.tsv"
    report.write_text("Genome\tComplete\n bin.1\t95\n")
    with pytest.raises(RuntimeError, match="Unexpected CheckM2"):
        checkm.summarize_checkm2(str(report))


# --- optional tool paths are actually reachable from the CLI ----------------------------


@pytest.mark.parametrize(
    "module,flag",
    [
        (bin_refinement, "--checkm2"),
        (annotate_bins, "--bakta"),
    ],
)
def test_optional_flags_are_accepted(module, flag):
    """These were documented in envs/*.yaml and docs but never implemented."""
    parser_argv = {
        bin_refinement: ["-o", "o", "-A", "a", flag],
        annotate_bins: ["-o", "o", "-b", "b", flag],
    }[module]
    args = module._parse_args(parser_argv)
    assert getattr(args, flag.lstrip("-").replace("-", "_")) is True


def test_classify_bins_accepts_gtdbtk():
    from metawrap2.modules import classify_bins

    args = classify_bins._parse_args(["-b", "bins", "-o", "out", "--gtdbtk"])
    assert args.gtdbtk is True


# --- envs/ discovery must not depend on the installation layout ------------------------


def test_envs_dir_is_found_and_has_every_module_yaml():
    from metawrap2.commands import install_env
    from metawrap2.config import MODULE_ENVS, OPTIONAL_ENVS

    d = install_env.envs_dir()
    assert os.path.isdir(d), d
    for module in MODULE_ENVS:
        assert os.path.isfile(os.path.join(d, "%s.yaml" % module)), module
    # The opt-in envs ship alongside them.
    for name in OPTIONAL_ENVS.values():
        yaml = os.path.join(d, name[len("metawrap2-") :] + ".yaml")
        assert os.path.isfile(yaml), yaml


def test_envs_dir_honours_override(monkeypatch, tmp_path):
    from metawrap2.commands import install_env

    monkeypatch.setenv("METAWRAP2_ENVS_DIR", str(tmp_path))
    assert install_env.envs_dir() == str(tmp_path)


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


def test_blastdb_name_is_configurable(tmp_path):
    from metawrap2 import blastdb

    (tmp_path / "core_nt.000.nsq").write_text("")
    s = _settings(BLASTDB=str(tmp_path))
    ok, _ = blastdb.find(s)
    assert not ok  # no database called "nt" here
    s = _settings(BLASTDB=str(tmp_path), BLASTDB_NAME="core_nt")
    ok, problem = blastdb.find(s)
    assert ok, problem
    assert blastdb.db_path(s) == os.path.join(str(tmp_path), "core_nt")


def test_blastdb_missing_gives_an_actionable_message(tmp_path):
    from metawrap2 import blastdb

    ok, problem = blastdb.find(_settings(BLASTDB=str(tmp_path)))
    assert not ok
    assert "BLASTDB_NAME" in problem and str(tmp_path) in problem
    ok, problem = blastdb.find(_settings())
    assert not ok and "BLASTDB is not set" in problem


# --- install-db writes config keys without destroying the user's file --------------------


def test_install_db_adds_databases_section(tmp_path):
    from metawrap2.commands import install_db

    cfg = tmp_path / "config.toml"
    cfg.write_text("[settings]\nthreads = 8\n")
    install_db._set_config_keys(str(cfg), {"CHECKM_DB": "/db/checkm"})
    text = cfg.read_text()
    assert "[settings]" in text and "threads = 8" in text  # untouched
    assert "[databases]" in text
    assert 'CHECKM_DB = "/db/checkm"' in text


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


def test_install_db_inserts_before_a_following_section(tmp_path):
    from metawrap2.commands import install_db

    cfg = tmp_path / "config.toml"
    cfg.write_text('[databases]\nTAXDUMP = "/tax"\n\n[modules.binning]\nMIN_CONTIG_LEN = 2000\n')
    install_db._set_config_keys(str(cfg), {"BLASTDB": "/blast"})
    lines = [ln.strip() for ln in cfg.read_text().splitlines() if ln.strip()]
    assert lines.index('BLASTDB = "/blast"') < lines.index("[modules.binning]")
    assert "MIN_CONTIG_LEN = 2000" in cfg.read_text()


def test_install_db_creates_a_missing_config(tmp_path):
    from metawrap2.commands import install_db

    cfg = tmp_path / "sub" / "config.toml"
    install_db._set_config_keys(str(cfg), {"BMTAGGER_DB": "/host"})
    assert 'BMTAGGER_DB = "/host"' in cfg.read_text()


def test_every_database_has_an_installer():
    from metawrap2.commands import install_db

    assert set(install_db.DATABASES) == set(install_db.INSTALLERS)
    # Every database that writes a config key must use a key `check` knows about.
    from metawrap2.config import DB_KEYS

    for name, db in install_db.DATABASES.items():
        if db.key:
            assert db.key in DB_KEYS, (name, db.key)


# --- progress bars must never corrupt data or logs --------------------------------------


def test_progress_bar_passes_items_through_unchanged():
    from metawrap2 import progress

    items = list(range(5))
    assert list(progress.bar(iter(items), desc="x")) == items


def test_progress_is_off_when_stderr_is_not_a_tty(monkeypatch):
    """Bars must not fill run.stderr with thousands of redraw lines."""
    from metawrap2 import progress

    monkeypatch.delenv("METAWRAP2_PROGRESS", raising=False)

    class NotATty:
        def isatty(self):
            return False

    monkeypatch.setattr(sys, "stderr", NotATty())
    assert progress.enabled() is False


def test_progress_env_override(monkeypatch):
    from metawrap2 import progress

    monkeypatch.setenv("METAWRAP2_PROGRESS", "0")
    assert progress.enabled() is False
    monkeypatch.setenv("METAWRAP2_PROGRESS", "1")
    # True only if tqdm is importable; either way it must not raise.
    assert isinstance(progress.enabled(), bool)


def test_progress_bar_writes_to_stderr_not_stdout(monkeypatch, capsys):
    """A bar on stdout would land inside a FASTA/SAM that a module redirects there."""
    from metawrap2 import progress

    monkeypatch.setenv("METAWRAP2_PROGRESS", "1")
    list(progress.bar(iter(range(3)), desc="test-bar"))
    captured = capsys.readouterr()
    assert "test-bar" not in captured.out


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


def test_dereplicate_remove_mode_drops_shared_contigs(tmp_path):
    from metawrap2 import refinement

    bins = tmp_path / "binsM"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">shared\nAAAA\n>only1\nCCCC\n")
    (bins / "bin.2.fa").write_text(">shared\nAAAA\n>only2\nGGGG\n")
    stats = tmp_path / "binsM.stats"
    stats.write_text(
        "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
        "bin.1\t90\t1\t0.5\tBacteria\t10000\t3000000\n"
        "bin.2\t80\t2\t0.4\tBacteria\t9000\t2500000\n"
    )
    out = tmp_path / "binsO"
    refinement.dereplicate(str(stats), str(bins), str(out), mode="remove")
    assert ">shared" not in (out / "bin.1.fa").read_text()
    assert ">shared" not in (out / "bin.2.fa").read_text()


def test_read_qc_removes_decompressed_input_copies(tmp_path):
    """Gzipped inputs were decompressed into the output dir and left there (GBs per sample)."""
    import gzip

    raw1 = tmp_path / "s_1.fastq.gz"
    raw2 = tmp_path / "s_2.fastq.gz"
    for path, seq in ((raw1, "ACGT"), (raw2, "TGCA")):
        with gzip.open(path, "wt") as fh:
            fh.write("@r1\n%s\n+\nIIII\n" % seq)
    out = tmp_path / "qc"

    assert (
        read_qc.main(
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
        == 0
    )
    assert (out / "final_pure_reads_1.fastq").is_file()
    # the decompressed intermediates are gone...
    assert not (out / "s_1.fastq").exists()
    assert not (out / "s_2.fastq").exists()
    # ...and the user's gzipped originals are untouched
    assert raw1.is_file() and raw2.is_file()


# --- --dry-run must print a module's WHOLE pipeline, not just the first command ----------
# Previously the orchestration consumed outputs that a dry run never produced, so every
# module died with a raw traceback (or a misleading "something went wrong with X") after one
# or two commands.

DRY_RUN_CASES = [
    ("read_qc", ["-1", "{r1}", "-2", "{r2}", "-o", "{out}"], 3),
    ("assembly", ["-1", "{r1}", "-2", "{r2}", "-o", "{out}", "--metaspades", "--megahit"], 4),
    ("kraken2", ["-o", "{out}", "{asm}", "{r1}"], 3),
    (
        "binning",
        ["-a", "{asm}", "-o", "{out}", "--metabat2", "--maxbin2", "--concoct", "{r1}", "{r2}"],
        8,
    ),
    ("bin_refinement", ["-o", "{out}", "-A", "{bins}", "-B", "{bins2}"], 3),
    ("reassemble_bins", ["-b", "{bins}", "-o", "{out}", "-1", "{r1}", "-2", "{r2}"], 6),
    ("quant_bins", ["-b", "{bins}", "-a", "{asm}", "-o", "{out}", "{r1}", "{r2}"], 2),
    ("classify_bins", ["-b", "{bins}", "-o", "{out}"], 4),
    ("annotate_bins", ["-b", "{bins}", "-o", "{out}"], 1),
    ("blobology", ["-a", "{asm}", "-o", "{out}", "{r1}", "{r2}"], 4),
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


def test_dry_run_creates_no_output_data(dry_inputs, monkeypatch, capsys):
    """A dry run may write its provenance record, but must not fabricate result data."""
    import os as _os

    from metawrap2 import command as command_mod
    from metawrap2.modules import binning as binning_mod

    monkeypatch.setattr(command_mod.runner, "dry_run", True)
    out = dry_inputs["out"] + "-nodata"
    binning_mod.main(
        [
            "-a",
            dry_inputs["asm"],
            "-o",
            out,
            "--metabat2",
            "--config",
            dry_inputs["config"],
            dry_inputs["r1"],
            dry_inputs["r2"],
        ]
    )

    # Recording the *plan* is the point of --dry-run, so these are expected.
    provenance = {
        "run_commands.txt",
        "run_config.json",
        "provenance.txt",
        "run.stdout",
        "run.stderr",
        "metawrap2.log",
    }
    produced = []
    for root, _dirs, files in _os.walk(out):
        for f in files:
            if f in provenance or root.endswith("steps"):
                continue
            produced.append(_os.path.relpath(_os.path.join(root, f), out))
    assert produced == [], produced

    # And the commands it would run are recorded for reproducibility.
    with open(_os.path.join(out, "run_commands.txt")) as fh:
        assert "metabat2" in fh.read()


# --- NCBI renamed "superkingdom" to "domain" --------------------------------------------
# A current taxdump contains zero "superkingdom" nodes, so asking for that rank silently
# produced "Not annotated" for every contig - losing the Bacteria/Archaea/Eukaryota split.


def test_rank_synonyms_cover_superkingdom_and_domain():
    from metawrap2.scripts.blobology import gc_cov_annotate as gca

    assert gca.normalize_ranks(["superkingdom"]) == {"superkingdom", "domain"}
    assert gca.normalize_ranks(["domain"]) == {"domain", "superkingdom"}
    assert gca.normalize_ranks(["genus"]) == {"genus"}
    # a taxdump "domain" node satisfies a request for "superkingdom", reported under the
    # name the caller asked for so the output columns do not change
    assert gca.canonical_rank("domain", ["superkingdom"]) == "superkingdom"
    assert gca.canonical_rank("superkingdom", ["domain"]) == "domain"
    assert gca.canonical_rank("genus", ["superkingdom"]) is None


def test_gc_cov_annotate_resolves_domain_as_superkingdom(tmp_path):
    """End-to-end over a miniature taxdump that uses the modern 'domain' rank."""
    from metawrap2.scripts.blobology import gc_cov_annotate as gca

    taxdump = tmp_path / "tax"
    taxdump.mkdir()
    # 561 (genus Escherichia) -> 1224 (phylum) -> 2 (domain Bacteria) -> 1 (root)
    (taxdump / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tdomain\t|\n"
        "1224\t|\t2\t|\tphylum\t|\n"
        "561\t|\t1224\t|\tgenus\t|\n"
    )
    (taxdump / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "1224\t|\tPseudomonadota\t|\t\t|\tscientific name\t|\n"
        "561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n"
    )

    assembly = tmp_path / "asm.fa"
    assembly.write_text(">NODE_1_length_8_cov_3.5\nACGTACGT\n")
    hits = tmp_path / "hits.tab"
    hits.write_text("NODE_1_length_8_cov_3.5\t561\n")
    out = tmp_path / "blobplot"

    argv = [
        "gc_cov_annotate",
        "--blasttaxid",
        str(hits),
        "--assembly",
        str(assembly),
        "--out",
        str(out),
        "--taxdump",
        str(taxdump),
        "--taxlist",
        "genus",
        "phylum",
        "superkingdom",
    ]
    monkey = sys.argv
    sys.argv = argv
    try:
        gca.main()
    finally:
        sys.argv = monkey

    header, row = out.read_text().splitlines()[:2]
    cols = dict(zip(header.split("\t"), row.split("\t")))
    # the column is still named after what the caller asked for...
    assert "taxlevel_superkingdom" in cols
    # ...and it is populated from the taxdump's "domain" node
    assert cols["taxlevel_superkingdom"] == "Bacteria"
    assert cols["taxlevel_genus"] == "Escherichia"
    assert cols["taxlevel_phylum"] == "Pseudomonadota"


def test_prune_blast_hits_keeps_domain_ranked_ids(tmp_path):
    from metawrap2.scripts import prune_blast_hits

    assert "domain" in prune_blast_hits.INCLUDE
    assert "superkingdom" in prune_blast_hits.INCLUDE

    nodes = tmp_path / "nodes.dmp"
    nodes.write_text("2\t|\t1\t|\tdomain\t|\n999\t|\t1\t|\tno rank\t|\n")
    raw = tmp_path / "raw.tab"
    # col 6 (index 5) is the staxids field
    raw.write_text(
        "q1\t1\t10\t100\ts1\t999;2\t1\t10\t50\t1e-9\t9\t10\n"
        "q2\t1\t10\t100\ts2\t999\t1\t10\t50\t1e-9\t9\t10\n"
    )
    kept = list(prune_blast_hits.prune(str(nodes), str(raw)))
    # q1's domain-ranked id survives; q2 has no ranked id at all and is dropped
    assert len(kept) == 1
    assert kept[0].split("\t")[5] == "2"


# --- prebuilt Kraken2 databases keep names.dmp at the top level, not under taxonomy/ -----


def test_kraken2_taxdump_found_in_classic_layout(tmp_path):
    from metawrap2.scripts import kraken2_translate as kt

    db = tmp_path / "db"
    (db / "taxonomy").mkdir(parents=True)
    (db / "taxonomy" / "names.dmp").write_text("")
    assert kt.find_taxdump_file(str(db), "names.dmp") == str(db / "taxonomy" / "names.dmp")


def test_kraken2_taxdump_found_in_prebuilt_layout(tmp_path):
    """This is the layout of every database downloaded from genome-idx."""
    from metawrap2.scripts import kraken2_translate as kt

    db = tmp_path / "db"
    db.mkdir()
    (db / "names.dmp").write_text("")
    assert kt.find_taxdump_file(str(db), "names.dmp") == str(db / "names.dmp")


def test_kraken2_taxdump_falls_back_to_configured_taxdump(tmp_path):
    from metawrap2.scripts import kraken2_translate as kt

    db = tmp_path / "db"
    db.mkdir()
    tax = tmp_path / "NCBI_tax"
    tax.mkdir()
    (tax / "nodes.dmp").write_text("")
    assert kt.find_taxdump_file(str(db), "nodes.dmp", (str(tax),)) == str(tax / "nodes.dmp")


def test_kraken2_taxdump_missing_names_every_location(tmp_path):
    from metawrap2.scripts import kraken2_translate as kt

    db = tmp_path / "db"
    db.mkdir()
    with pytest.raises(FileNotFoundError) as exc:
        kt.find_taxdump_file(str(db), "names.dmp")
    message = str(exc.value)
    assert "taxonomy" in message and "install-db taxdump" in message


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


# --- an empty consolidated bin set must explain itself, not fail inside CheckM -----------


def test_require_bins_message_is_actionable(tmp_path):
    from metawrap2 import checkm as checkm_mod

    empty = tmp_path / "binsO"
    empty.mkdir()
    assert checkm_mod.count_bins(str(empty)) == 0
    with pytest.raises(RuntimeError, match="No bins to analyse"):
        checkm_mod.require_bins(str(empty), "extra context here")

    (empty / "bin_001.fasta").write_text(">c\nACGT\n")
    assert checkm_mod.count_bins(str(empty)) == 1
    checkm_mod.require_bins(str(empty))  # no longer raises


def test_bin_refinement_explains_thresholds_when_nothing_qualifies(tmp_path, monkeypatch, capsys):
    """Previously this died with CheckM's misleading 'Check the extension (-x)' error."""
    from metawrap2.modules import bin_refinement as br

    out = tmp_path / "out"
    (out / "binsO").mkdir(parents=True)  # consolidation produced nothing
    (out / "binsA").mkdir(parents=True)
    (out / "binsA" / "bin.1.fa").write_text(">c\nACGT\n")

    # error() raises SystemExit; make sure we get there with a useful message, and that we
    # never reach CheckM.
    monkeypatch.setattr(
        br._checkm, "run_checkm", lambda *a, **k: pytest.fail("CheckM should not be reached")
    )

    with pytest.raises(SystemExit):
        br.error(
            "Consolidation produced no bins: not one bin from any bin set met your "
            "quality thresholds of >50%% completeness and <10%% contamination."
        )
    assert "quality thresholds" in capsys.readouterr().out


# --- commands that run with cwd= must not be given paths relative to the *old* cwd -------


def test_bin_refinement_plot_paths_are_absolute(tmp_path, monkeypatch):
    """With a relative -o this produced "cd out && ... out/binsA.stats" (= out/out/...)."""
    from metawrap2.modules import bin_refinement as br

    out = tmp_path / "REFINE"
    out.mkdir()
    (out / "binsA.stats").write_text("bin\tcompleteness\tcontamination\n")
    (out / "binsO.stats").write_text("bin\tcompleteness\tcontamination\n")

    recorded = {}

    def fake_run(cmd, **kw):
        recorded["cmd"] = cmd
        recorded["cwd"] = kw.get("cwd")
        return 0

    monkeypatch.setattr(br, "run", fake_run)
    monkeypatch.chdir(tmp_path)

    stats_files = " ".join(
        sorted(
            os.path.abspath(os.path.join("REFINE", f))
            for f in os.listdir("REFINE")
            if f.endswith(".stats")
        )
    )
    br.run(
        br.PLOT_BINNING.format(comp=50, cont=10, stats=stats_files),
        env=None,
        tool="plot_binning_results",
        cwd="REFINE",
    )

    # every .stats argument must be absolute so cwd= cannot change what it resolves to
    args = [a for a in recorded["cmd"].split() if a.endswith(".stats")]
    assert args, recorded["cmd"]
    for a in args:
        assert os.path.isabs(a), a
        assert os.path.isfile(a), a


# --- log noise filtering -----------------------------------------------------------------
# Real measurement: INITIAL_BINNING/run.stderr was 1.8 MB for a 3-sample run, 23,770 of its
# 24,770 lines being two metaBAT2 warnings repeated ~11,900 times each.


def test_signature_groups_lines_differing_only_in_numbers():
    from metawrap2.logfilter import signature

    a = "WARNING: calculated a huge mean=1.2e+03. correctedLen=456 contigDepth=789"
    b = "WARNING: calculated a huge mean=9.9e+02. correctedLen=12 contigDepth=3"
    assert signature(a) == signature(b)
    assert signature("different message entirely") != signature(a)


def test_repeated_messages_are_capped_and_counted():
    from metawrap2.logfilter import LogFilter

    f = LogFilter(max_repeats=3)
    written = [f.feed("huge mean=%d contigDepth=%d\n" % (i, i * 2)) for i in range(100)]
    kept = [w for w in written if w is not None]
    # 3 verbatim + 1 "further occurrences are being counted" notice
    assert len(kept) == 4
    distinct, suppressed = f.stats()
    assert distinct == 1
    assert suppressed == 96
    summary = f.summary()
    assert len(summary) == 1
    assert "suppressed 96 further occurrences" in summary[0]


def test_distinct_messages_are_never_suppressed():
    from metawrap2.logfilter import LogFilter

    f = LogFilter(max_repeats=2)
    # Bodies must differ in their *non-numeric* text: signature() masks numbers, so
    # "step 1 of 10" and "step 2 of 10" are deliberately the same message.
    words = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa"]
    kept = [f.feed("processing %s\n" % w) for w in words]
    assert all(k is not None for k in kept)
    assert f.stats() == (10, 0)


def test_progress_bar_redraws_collapse_to_final_state():
    from metawrap2.logfilter import LogFilter, collapse_redraws

    assert collapse_redraws("10%\r50%\r100%") == "100%"
    assert collapse_redraws("no carriage returns") == "no carriage returns"
    assert collapse_redraws("\r\r") == ""

    f = LogFilter()
    out = f.feed("progress: 10%\rprogress: 55%\rprogress: 100%\n")
    assert out == "progress: 100%\n"


def test_blank_lines_pass_through():
    from metawrap2.logfilter import LogFilter

    f = LogFilter(max_repeats=1)
    assert f.feed("\n") == "\n"
    assert f.feed("\n") == "\n"  # blank lines are not counted as repeats


def test_filter_can_be_disabled():
    from metawrap2.logfilter import LogFilter

    f = LogFilter(max_repeats=1, enabled=False)
    assert all(f.feed("same line\n") == "same line\n" for _ in range(50))
    assert f.summary() == []


def test_command_data_output_is_never_filtered(tmp_path, monkeypatch):
    """stdout redirected to a data file (bwa mem > x.sam) must be byte-exact."""
    from metawrap2.command import CommandRunner

    data = tmp_path / "out.txt"
    runner = CommandRunner()
    # 100 identical lines: as a *log* these would be collapsed; as data they must all survive.
    script = "for i in $(seq 100); do echo 'identical data line'; done"
    rc = runner.run(["bash", "-c", script], stdout_path=str(data))
    assert rc == 0
    assert data.read_text().count("identical data line") == 100


def test_command_run_log_is_filtered(tmp_path):
    from metawrap2.command import CommandRunner

    log = tmp_path / "run.stderr"
    runner = CommandRunner(run_stderr_path=str(log))
    script = 'for i in $(seq 200); do echo "noise number $i" >&2; done'
    assert runner.run(["bash", "-c", script]) == 0
    text = log.read_text()
    assert text.count("noise number") < 200  # collapsed
    assert "suppressed" in text  # and said so


def test_verbose_logs_keeps_everything(tmp_path):
    from metawrap2.command import CommandRunner

    log = tmp_path / "run.stderr"
    runner = CommandRunner(run_stderr_path=str(log), verbose_logs=True)
    script = 'for i in $(seq 200); do echo "noise number $i" >&2; done'
    assert runner.run(["bash", "-c", script]) == 0
    # The log also opens with a "$ <command>" header that quotes the script, so count the
    # emitted lines rather than substring occurrences.
    emitted = [ln for ln in log.read_text().splitlines() if ln.startswith("noise number ")]
    assert len(emitted) == 200
    assert "suppressed" not in log.read_text()


def test_error_tail_is_unfiltered(tmp_path):
    """The ToolError message must show what the tool actually last printed."""
    from metawrap2.command import CommandRunner, ToolError

    runner = CommandRunner()
    script = (
        'for i in $(seq 100); do echo "repeated noise $i" >&2; done; '
        "echo 'THE REAL ERROR' >&2; exit 3"
    )
    with pytest.raises(ToolError) as exc:
        runner.run(["bash", "-c", script], tool="thing")
    assert "THE REAL ERROR" in str(exc.value)


# --- relative -o must work even for tools that run with cwd=<output dir> -----------------


def test_absolutize_paths_makes_every_path_argument_absolute(tmp_path, monkeypatch):
    from metawrap2.modules._common import absolutize_paths

    monkeypatch.chdir(tmp_path)

    class Args:
        output = "OUT"
        assembly = "asm.fa"
        bins = "bins"
        bins_a = "binsA"
        bins_b = None
        reads_1 = "s_1.fastq"
        reads_2 = "s_2.fastq"
        reads: ClassVar = ["a_1.fastq", "a_2.fastq"]
        seqs: ClassVar = ["x.fa"]
        threads = 4  # untouched non-path argument

    args = Args()
    absolutize_paths(args)
    assert os.path.isabs(args.output) and args.output.endswith("OUT")
    assert os.path.isabs(args.assembly)
    assert os.path.isabs(args.bins_a)
    assert args.bins_b is None  # unset stays unset
    assert all(os.path.isabs(r) for r in args.reads)
    assert all(os.path.isabs(s) for s in args.seqs)
    assert args.threads == 4


def test_absolutize_paths_expands_user(monkeypatch):
    from metawrap2.modules._common import absolutize_paths

    class Args:
        output = "~/somewhere/out"

    args = Args()
    absolutize_paths(args)
    assert "~" not in args.output
    assert os.path.isabs(args.output)


@pytest.mark.parametrize(
    "module_name,argv_template",
    [
        ("classify_bins", ["-b", "{bins}", "-o", "REL_OUT"]),
        ("bin_refinement", ["-A", "{bins}", "-o", "REL_OUT"]),
    ],
)
def test_relative_output_produces_absolute_command_paths(
    module_name, argv_template, dry_inputs, monkeypatch, capsys
):
    """A relative -o used to yield "cd REL_OUT && ... REL_OUT/file" (= REL_OUT/REL_OUT/file)."""
    import importlib

    from metawrap2 import command as command_mod

    monkeypatch.chdir(os.path.dirname(dry_inputs["bins"]))
    monkeypatch.setattr(command_mod.runner, "dry_run", True)
    monkeypatch.setattr(command_mod.runner, "force", True)

    module = importlib.import_module("metawrap2.modules.%s" % module_name)
    argv = [a.format(**dry_inputs) for a in argv_template] + ["--config", dry_inputs["config"]]
    try:
        module.main(argv)
    finally:
        command_mod.set_recorder(None)

    printed = [ln for ln in capsys.readouterr().err.splitlines() if ln.startswith("[dry-run]")]
    assert printed
    # no command may contain the doubled relative path
    for line in printed:
        assert "REL_OUT/REL_OUT" not in line, line


# --- run logs must describe the run that produced the data beside them -------------------


def test_fresh_run_truncates_the_previous_run_log(tmp_path, monkeypatch):
    """Without --resume, run.stdout/run.stderr used to accumulate across runs."""
    from metawrap2 import command as command_mod
    from metawrap2.config import Settings
    from metawrap2.modules._common import start_run

    out = tmp_path / "out"
    out.mkdir()
    (out / "run.stderr").write_text("OUTPUT FROM A PREVIOUS RUN\n")
    (out / "run.stdout").write_text("OUTPUT FROM A PREVIOUS RUN\n")

    class Args:
        output = str(out)
        config = None

    monkeypatch.setattr(command_mod.runner, "resume", False)
    monkeypatch.setattr(command_mod.runner, "force", True)
    try:
        start_run("binning", Args(), None, Settings(), inputs=[])
    finally:
        command_mod.set_recorder(None)
        command_mod.set_run_logs(None, None)

    assert (out / "run.stderr").read_text() == ""
    assert (out / "run.stdout").read_text() == ""


def test_resumed_run_keeps_the_existing_log(tmp_path, monkeypatch):
    from metawrap2 import command as command_mod
    from metawrap2.config import Settings
    from metawrap2.modules._common import start_run

    out = tmp_path / "out"
    out.mkdir()
    (out / "run.stderr").write_text("EARLIER PART OF THIS RUN\n")

    class Args:
        output = str(out)
        config = None

    monkeypatch.setattr(command_mod.runner, "resume", True)
    try:
        start_run("binning", Args(), None, Settings(), inputs=[])
    finally:
        command_mod.set_recorder(None)
        command_mod.set_run_logs(None, None)

    assert "EARLIER PART OF THIS RUN" in (out / "run.stderr").read_text()


# --- the config template must not drift from the shipped example -------------------------


def test_example_config_matches_the_template():
    """`config init` and metawrap2.toml.example were two hand-maintained copies that diverged."""
    from metawrap2.commands.config_cmd import DEFAULT_CONFIG

    root = os.path.dirname(os.path.dirname(os.path.abspath(__import__("metawrap2").__file__)))
    example = os.path.join(os.path.dirname(root), "metawrap2.toml.example")
    if not os.path.isfile(example):
        pytest.skip("not running from a source checkout")
    with open(example) as fh:
        assert (
            fh.read() == DEFAULT_CONFIG
        ), "metawrap2.toml.example has drifted from config_cmd.DEFAULT_CONFIG; regenerate it"


def test_template_is_valid_toml_and_covers_every_database_key():
    try:
        import tomllib as toml
    except ModuleNotFoundError:
        import tomli as toml

    from metawrap2.commands.config_cmd import DEFAULT_CONFIG
    from metawrap2.config import DB_KEYS

    data = toml.loads(DEFAULT_CONFIG)
    assert data["settings"]["threads"] >= 1
    # every key the checker knows about is at least mentioned (set or commented) in the template
    for key in DB_KEYS:
        assert key in DEFAULT_CONFIG, key


# --- fail fast on empty and truncated inputs ---------------------------------------------
# An empty FASTQ used to produce an empty assembly, then zero bins, then a CheckM error about
# file extensions - hours and three modules away from the actual cause.


def test_require_file_rejects_missing_empty_and_directories(tmp_path):
    from metawrap2 import validate

    assert "does not exist" in validate.require_file(str(tmp_path / "nope.fastq"), "reads")
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    assert "is empty" in validate.require_file(str(empty), "reads")
    assert "is a directory" in validate.require_file(str(tmp_path), "reads")
    assert validate.require_file(__file__, "this test") is None
    assert "no reads was given" in validate.require_file("", "reads")


def test_require_nonempty_dir(tmp_path):
    from metawrap2 import validate

    empty = tmp_path / "bins"
    empty.mkdir()
    assert "is empty" in validate.require_nonempty_dir(str(empty), "bin folder")
    (empty / "notes.txt").write_text("hello")
    assert "no .fa/.fasta files" in validate.require_nonempty_dir(
        str(empty), "bin folder", (".fa", ".fasta")
    )
    (empty / "bin_001.fasta").write_text(">c\nACGT\n")
    assert validate.require_nonempty_dir(str(empty), "bin folder", (".fa", ".fasta")) is None


def test_check_fastq_accepts_a_complete_file(tmp_path):
    from metawrap2 import validate

    good = tmp_path / "good.fastq"
    good.write_text("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n")
    assert validate.check_fastq(str(good)) is None


@pytest.mark.parametrize(
    "content,expected",
    [
        ("@r1\nACGT\n+\n", "ends mid-record"),  # missing quality line
        ("@r1\nACGT\n+\nIIII", "does not end with a newline"),  # truncated mid-line
        (">r1\nACGT\n", "does not look like FASTQ"),  # actually FASTA
    ],
)
def test_check_fastq_rejects_truncated_and_wrong_format(tmp_path, content, expected):
    from metawrap2 import validate

    bad = tmp_path / "bad.fastq"
    bad.write_text(content)
    problem = validate.check_fastq(str(bad))
    assert problem and expected in problem


def test_check_fastq_handles_gzip(tmp_path):
    import gzip

    from metawrap2 import validate

    good = tmp_path / "good.fastq.gz"
    with gzip.open(good, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")
    assert validate.check_fastq(str(good)) is None


def test_check_fastq_reports_a_corrupt_gzip_stream(tmp_path):
    from metawrap2 import validate

    broken = tmp_path / "broken.fastq.gz"
    broken.write_bytes(b"\x1f\x8b\x08\x00garbage-not-a-real-gzip-stream")
    problem = validate.check_fastq(str(broken))
    assert problem and ("truncated or corrupt" in problem or "could not be read" in problem)


def test_check_fasta(tmp_path):
    from metawrap2 import validate

    good = tmp_path / "a.fa"
    good.write_text(">c1\nACGT\n")
    assert validate.check_fasta(str(good)) is None

    headerless = tmp_path / "b.fa"
    headerless.write_text("ACGT\n")
    assert "does not look like FASTA" in validate.check_fasta(str(headerless))

    dangling = tmp_path / "c.fa"
    dangling.write_text(">c1\nACGT\n>c2\n")
    assert "ends with a header line and no sequence" in validate.check_fasta(str(dangling))


def test_require_inputs_reports_every_problem_at_once(tmp_path):
    from metawrap2 import validate

    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    with pytest.raises(validate.InputProblem) as exc:
        validate.require_inputs(
            [
                (validate.check_fastq, str(tmp_path / "missing.fastq"), "forward reads"),
                (validate.check_fastq, str(empty), "reverse reads"),
            ],
            context="read_qc",
        )
    message = str(exc.value)
    assert "read_qc" in message
    assert "2 input problem(s)" in message
    assert "missing.fastq" in message and "empty.fastq" in message


def test_module_refuses_to_start_on_an_empty_read_file(tmp_path, capsys):
    """End to end through the module: the error names the file and never runs a tool."""
    from metawrap2.modules import read_qc as read_qc_mod

    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@r1\nACGT\n+\nIIII\n")
    r2.write_text("")  # the mate is empty
    with pytest.raises(SystemExit):
        read_qc_mod.main(
            ["-1", str(r1), "-2", str(r2), "-o", str(tmp_path / "out"), "--skip-bmtagger"]
        )
    out = capsys.readouterr().out
    assert "input problem" in out and "s_2.fastq" in out


def test_validation_is_skipped_under_dry_run(monkeypatch, tmp_path):
    """--dry-run prints a plan; it must not demand that the data already be valid."""
    from metawrap2 import command as command_mod
    from metawrap2.modules import _common

    monkeypatch.setattr(command_mod.runner, "dry_run", True)
    # a path that does not exist at all - would fail validation if it were run
    _common.validate_inputs("test", [(_common.check_fastq, str(tmp_path / "nope.fastq"), "reads")])
