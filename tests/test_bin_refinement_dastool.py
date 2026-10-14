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
from metawrap2.constants import BIN_EXTENSION, STATS_SUFFIX
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

    (cmd,) = [c for c in fake_tools if c.startswith("DAS_Tool ")]
    # DAS_Tool's default --search_engine is usearch, which is not installed, and it fails partway
    # into a run; without --write_bins it reports which bins it would keep and writes no FASTA.
    assert "--search_engine diamond" in cmd and "usearch" not in cmd and "--write_bins" in cmd
    # -i and -l are positional with respect to each other: a mismatch mislabels every bin.
    parts = cmd.split()
    tables = parts[parts.index("-i") + 1].split(",")
    assert parts[parts.index("-l") + 1].split(",") == ["binsA", "binsB"]
    assert [os.path.basename(t) for t in tables] == ["binsA.dastool.tsv", "binsB.dastool.tsv"]
    # The contig tables carry no header and exactly two columns: a header row would be read as a
    # contig, and the run would end having consolidated nothing.
    with open(os.path.join(work, "binsA.dastool.tsv")) as fh:
        rows = [line.split("\t") for line in fh.read().splitlines()]
    assert rows and all(len(r) == 2 for r in rows) and rows[0][0] == "c01"


def test_an_assembly_with_an_empty_record_is_filtered_before_a_consolidator_sees_it(tmp_path):
    """DAS_Tool's assembly reader stops on a header with no sequence, and blames the contigs.

    R's ``fread`` gets one field where it expects two, stops part way through the ``.seqlength`` table,
    and DAS_Tool then reports every contig past that point as "not found in assembly" - 36,854 of them
    in the run that found this, in a message naming none of the real problem. CAMI2's gold-standard
    assembly ships two such records; `binning` already drops them from its own copy, so a consolidator
    reading the original was getting what binning had deliberately not passed on.
    """
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\nACGT\n>_from__to_-42_total_0\n\n>c2\nTTTT\n")
    work = tmp_path / "work"

    out = bin_refinement.assembly_for_bins(str(assembly), str(work))
    assert out != str(assembly), "a copy is needed when a record has to be dropped"
    with open(out) as fh:
        text = fh.read()
    assert "_from__to_" not in text
    assert ">c1\nACGT\n" in text and ">c2\nTTTT\n" in text
