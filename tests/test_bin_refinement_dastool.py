"""The --consolidator dastool route through bin_refinement.

DAS_Tool itself is not installed here (it is an opt-in env), so it is replaced by a stand-in that
writes the bins a real run would write. That is the right thing to fake: what is under test is
MetaWrap2's half - which tables DAS_Tool is handed, which assembly, where its bins are picked up
from, and that everything downstream of ``binsO`` behaves as it does for the default consolidator.
Whether DAS_Tool's own scoring is any good is a question for the benchmark, not for a unit test.

CheckM is faked the same way the golden tests fake it, and for the same reason.
"""

import hashlib
import os

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2.constants import BIN_EXTENSION, CONTIGS_SUFFIX, STATS_SUFFIX
from metawrap2.io.seqio import contig_id, iter_fasta
from metawrap2.modules import _common, bin_refinement

CONTIGS = {"c%02d" % i: ("ACGT" * (50000 + i * 2500)) for i in range(1, 13)}

BIN_SETS = {
    "A": {"bin.1": ["c01", "c02", "c03", "c04"], "bin.2": ["c05", "c06", "c07", "c08"]},
    "B": {"bin.1": ["c01", "c02", "c03"], "bin.2": ["c05", "c06", "c07", "c08", "c09"]},
}


def _fake_run_checkm(bins_dir, **kwargs):
    """A .stats file for *bins_dir*, with scores derived from each bin's contig set."""
    rows = []
    for name in sorted(os.listdir(bins_dir)):
        if not name.endswith(BIN_EXTENSION):
            continue
        path = os.path.join(bins_dir, name)
        contigs = sorted(contig_id(h) for h, _ in iter_fasta(path))
        digest = hashlib.sha256("|".join(contigs).encode()).digest()
        rows.append(
            [
                os.path.splitext(name)[0],
                str(round(80.0 + (digest[0] / 255.0) * 20.0, 2))[:5],  # 80..100, always "good"
                str(round((digest[1] / 255.0) * 5.0, 3))[:5],
                "0.500",
                "Bacteria",
                str(1000 + digest[2] * 37),
                str(sum(len(s) for _, s in iter_fasta(path))),
            ]
        )
    stats_path = bins_dir + checkm_mod.STATS_SUFFIX
    checkm_mod.write_stats(rows, stats_path)
    return stats_path


@pytest.fixture
def fake_tools(monkeypatch):
    """Fake CheckM and DAS_Tool. Returns the list of commands that reached :func:`run`.

    The DAS_Tool stand-in reads the tables it was given and writes one bin per table group, so a
    test can tell apart "MetaWrap2 handed it the wrong tables" from "MetaWrap2 looked in the wrong
    place for the output" - both of which would otherwise look like an empty consolidation.
    """
    monkeypatch.setattr(bin_refinement._checkm, "run_checkm", _fake_run_checkm)
    monkeypatch.setattr(checkm_mod, "run_checkm", _fake_run_checkm)
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)
    # The opt-in env is not built in CI, and nothing here enters it: DAS_Tool is faked below.
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: True)
    commands = []
    real_run = bin_refinement.run

    def fake_run(cmd, **kwargs):
        commands.append(cmd)
        if not cmd.startswith("DAS_Tool "):
            return real_run(cmd, **kwargs)
        parts = cmd.split()
        tables = parts[parts.index("-i") + 1].split(",")
        assembly = parts[parts.index("-c") + 1]
        prefix = parts[parts.index("-o") + 1]
        seqs = {contig_id(h): s for h, s in iter_fasta(assembly)}
        # The first table wins every contig it claims, which is a plausible-shaped answer: one
        # non-overlapping bin set drawn from the assembly's own sequences.
        chosen = {}
        for table in tables:
            with open(table) as fh:
                for line in fh:
                    contig, bin_name = line.rstrip("\n").split("\t")
                    chosen.setdefault(contig, bin_name)
        bins = {}
        for contig, bin_name in chosen.items():
            bins.setdefault(bin_name, []).append(contig)
        produced = prefix + "_DASTool_bins"
        os.makedirs(produced, exist_ok=True)
        for bin_name, contigs in bins.items():
            with open(os.path.join(produced, bin_name + ".fa"), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, seqs[c]) for c in sorted(contigs))
        return 0

    monkeypatch.setattr(bin_refinement, "run", fake_run)
    return commands


def _inputs(tmp_path, letters="AB"):
    """Write the bin sets and the assembly they came from; returns (dirs, assembly path)."""
    dirs = []
    for letter in letters:
        directory = tmp_path / ("binner_%s" % letter)
        directory.mkdir()
        for bin_name, contigs in BIN_SETS[letter].items():
            with open(directory / ("%s.fa" % bin_name), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, CONTIGS[c]) for c in contigs)
        dirs.append(str(directory))
    assembly = tmp_path / "assembly.fa"
    with open(assembly, "w") as fh:
        fh.writelines(">%s\n%s\n" % (c, s) for c, s in sorted(CONTIGS.items()))
    return dirs, str(assembly)


def _run(tmp_path, extra=(), letters="AB"):
    dirs, assembly = _inputs(tmp_path, letters)
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10", "--skip-gunc"]
    for flag, path in zip(("-A", "-B", "-C"), dirs):
        argv += [flag, path]
    argv += ["--consolidator", "dastool", "-a", assembly, *extra]
    assert bin_refinement.main(argv) == 0
    return str(out)


# --- the command ---------------------------------------------------------------------------


def test_the_command_names_an_engine_that_is_actually_in_the_env():
    """DAS_Tool's default --search_engine is usearch, which is not free software and not installed.

    Left at the default the command does not fail at solve time or import time - it fails partway
    into a run, which on real data is after prodigal has read every bin.
    """
    cmd = bin_refinement.DAS_TOOL.format(
        tables="a.tsv,b.tsv",
        labels="binsA,binsB",
        assembly="asm.fa",
        prefix="p",
        threads=4,
        score=0.5,
    )
    assert "{" not in cmd
    assert "--search_engine diamond" in cmd and "usearch" not in cmd
    assert "-i a.tsv,b.tsv" in cmd and "-l binsA,binsB" in cmd
    # Without --write_bins DAS_Tool reports which bins it would keep and writes no FASTA at all.
    assert "--write_bins" in cmd


# --- the contig table ----------------------------------------------------------------------


def test_the_contig_table_is_contig_to_bin_with_no_header(tmp_path):
    """DAS_Tool reads this with no header and expects exactly two columns.

    A third column or a header row is not rejected: the header is read as a contig that is not in
    the assembly, and the run ends having consolidated nothing.
    """
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / ("bin.1" + BIN_EXTENSION)).write_text(">c1 len=10\nAAAA\n>c2\nCCCC\n")
    (bins / ("bin.2" + BIN_EXTENSION)).write_text(">c3\nGGGG\n")
    (bins / "notes.txt").write_text("ignored\n")
    # A stats sidecar sits next to every bin set, and it is not a bin.
    (bins / ("bin.3" + STATS_SUFFIX)).write_text("bin\tcompleteness\n")
    table = tmp_path / "t.tsv"
    assert bin_refinement.write_contig_table(str(bins), str(table)) == 3
    assert table.read_text() == "c1\tbin.1\nc2\tbin.1\nc3\tbin.2\n"


def test_the_same_table_is_written_beside_the_final_bins(tmp_path, fake_tools):
    """One writer for both consumers, so a change for DAS_Tool cannot silently break Anvio."""
    out = _run(tmp_path)
    final = [d for d in os.listdir(out) if d.startswith("metawrap_") and d.endswith("_bins")]
    assert len(final) == 1
    table = os.path.join(out, final[0] + CONTIGS_SUFFIX)
    with open(table) as fh:
        rows = [line.split("\t") for line in fh.read().splitlines()]
    assert rows and all(len(r) == 2 for r in rows)


# --- the assembly --------------------------------------------------------------------------


def test_an_assembly_whose_headers_already_match_is_passed_through_untouched(tmp_path):
    """metaSPAdes and MEGAHIT names normalise to themselves, so the common case copies nothing."""
    assembly = tmp_path / "a.fa"
    assembly.write_text(">NODE_1_length_5_cov_2 extra\nAAAAA\n>k141_7 flag=1\nCCCCC\n")
    work = tmp_path / "work"
    assert bin_refinement.assembly_for_bins(str(assembly), str(work)) == str(assembly)
    assert not work.exists()


def test_an_assembly_is_renamed_when_the_bins_were(tmp_path):
    """'=' inside a contig id is what _fix_contig_naming rewrites, so the assembly must match.

    Unmatched names are not an error anywhere: DAS_Tool simply finds none of the bins' contigs in
    the assembly and keeps nothing.
    """
    assembly = tmp_path / "a.fa"
    assembly.write_text(">ctg=1\nAAAAA\n")
    work = tmp_path / "work"
    result = bin_refinement.assembly_for_bins(str(assembly), str(work))
    assert result != str(assembly)
    assert [contig_id(h) for h, _ in iter_fasta(result)] == ["ctg_1"]


# --- routing -------------------------------------------------------------------------------


def test_the_consolidated_set_lands_in_binsO_and_runs_the_normal_downstream_half(
    tmp_path, fake_tools
):
    """What makes the arms comparable: only the consolidation step differs.

    binsO is re-scored with CheckM, cut on -c/-x, plotted and laid out by the same code the
    default consolidator reaches, so a comparison between the two measures consolidation and not
    two different definitions of a final bin set.
    """
    out = _run(tmp_path)
    # Working bin sets are tidied into work_files at the end, binsO along with the rest.
    work = os.path.join(out, "work_files")
    # Numbered the way every other MetaWrap2 bin set is, not left under DAS_Tool's own names:
    # .stats.tsv rows are keyed on these stems, and the arms have to differ only in contents.
    assert sorted(os.listdir(os.path.join(work, "binsO"))) == ["bin_001.fasta", "bin_002.fasta"]
    with open(os.path.join(work, "dastool_bin_names.tsv")) as fh:
        names = fh.read().splitlines()
    assert names[0] == "bin\tdastool_name"
    assert [row.split("\t")[0] for row in names[1:]] == ["bin_001", "bin_002"]
    assert os.path.isfile(os.path.join(work, "binsO" + STATS_SUFFIX))
    final = os.path.join(out, "metawrap_50_10_bins")
    assert sorted(os.listdir(final)) == ["bin_001.fasta", "bin_002.fasta"]
    # binsM is this module's own intermediate; the external route never builds one.
    assert not os.path.isdir(os.path.join(work, "binsM"))


def test_checkm_runs_once_and_not_on_every_input_set(tmp_path, fake_tools, monkeypatch):
    """CheckM on the input sets exists to feed MetaWrap2's own choice between bin versions.

    DAS_Tool makes that choice from its own scoring, so running CheckM anyway would add the
    module's largest cost to the arm and change nothing - and would make the cost comparison
    between the two consolidators meaningless.
    """
    scored = []

    def counting_checkm(bins_dir, **kwargs):
        scored.append(os.path.basename(bins_dir))
        return _fake_run_checkm(bins_dir, **kwargs)

    monkeypatch.setattr(bin_refinement._checkm, "run_checkm", counting_checkm)
    monkeypatch.setattr(checkm_mod, "run_checkm", counting_checkm)
    _run(tmp_path)
    assert scored == ["binsO"]


def test_no_bin_set_is_hybridised(tmp_path, fake_tools):
    """Binning_refiner builds every combination up front; DAS_Tool consumes the sets as given.

    Skipping it is most of why this route is cheap, and it is also a correctness point: DAS_Tool
    scoring a hybrid of two binners alongside the binners themselves is not what its scoring is
    for, and would feed it the same contigs several times over.
    """
    out = _run(tmp_path)
    assert not [cmd for cmd in fake_tools if "binning_refiner" in cmd]
    assert "binsAB" not in os.listdir(out)


def test_the_tables_handed_over_are_the_input_sets_in_label_order(tmp_path, fake_tools):
    """-i and -l are positional with respect to each other: a mismatch mislabels every bin."""
    out = _run(tmp_path)
    cmd = next(c for c in fake_tools if c.startswith("DAS_Tool "))
    parts = cmd.split()
    tables = parts[parts.index("-i") + 1].split(",")
    labels = parts[parts.index("-l") + 1].split(",")
    assert labels == ["binsA", "binsB"]
    assert [os.path.basename(t) for t in tables] == ["binsA.dastool.tsv", "binsB.dastool.tsv"]
    assert all(t.startswith(os.path.join(out, "work_files")) for t in tables)


# --- preconditions -------------------------------------------------------------------------


def test_a_missing_assembly_stops_the_run_naming_the_flag_to_pass(tmp_path, fake_tools, capsys):
    """-a is only needed by this route, so it cannot be a required argument; it has to be checked.

    The message has to name -a. Left to the generic input validation the run still stops in the
    first second, but it stops saying an assembly path is empty, which for someone who never
    passed one reads as a bug rather than a missing argument.
    """
    dirs, _assembly = _inputs(tmp_path)
    out = tmp_path / "refined"
    with pytest.raises(SystemExit):
        bin_refinement.main(
            ["-o", str(out), "-t", "1", "-A", dirs[0], "-B", dirs[1], "--consolidator", "dastool"]
        )
    assert "-a/--assembly" in capsys.readouterr().out
    assert not os.path.isdir(os.path.join(str(out), "binsA"))


def test_skip_consolidation_with_a_consolidator_is_reported_as_a_contradiction(
    tmp_path, fake_tools
):
    """Both flags were typed on purpose, and either precedence silently does the other thing."""
    with pytest.raises(SystemExit):
        _run(tmp_path, extra=["--skip-consolidation"])


def test_an_assembly_with_an_empty_record_is_filtered_before_a_consolidator_sees_it(tmp_path):
    """DAS_Tool's assembly reader stops on a header with no sequence, and blames the contigs.

    R's ``fread`` gets one field where it expects two, stops part way through the ``.seqlength`` table,
    and DAS_Tool then reports every contig past that point as "not found in assembly" - 36,854 of them
    in the run that found this, in a message naming none of the real problem. CAMI2's gold-standard
    assembly ships two such records; `binning` already drops them from its own copy, so a consolidator
    reading the original was getting what binning had deliberately not passed on.
    """
    from metawrap2.modules import bin_refinement

    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>_from__to_-42_total_0\n\n>c2\nTTTT\n")
    work = tmp_path / "work"

    out = bin_refinement.assembly_for_bins(str(assembly), str(work))
    assert out != str(assembly), "a copy is needed when a record has to be dropped"
    with open(out) as fh:
        text = fh.read()
    assert "_from__to_" not in text
    assert ">c1\nACGT\n" in text and ">c2\nTTTT\n" in text


def test_a_clean_assembly_is_passed_through_without_a_copy(tmp_path):
    """Assemblies are gigabytes; a copy nothing needed would be the module's most expensive no-op."""
    from metawrap2.modules import bin_refinement

    assembly = tmp_path / "asm.fa"
    assembly.write_text(">k141_1 len=4\nACGT\n>k141_2\nTTTT\n")
    assert bin_refinement.assembly_for_bins(str(assembly), str(tmp_path / "work")) == str(assembly)
