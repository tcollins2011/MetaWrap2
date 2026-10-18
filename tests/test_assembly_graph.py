"""``assembly --keep-graph``: the one output a graph-aware binner needs, and its disk cost.

PlasMAAG, GraphMB and Hi-C placement all read the assembly graph rather than the contigs, and
until this flag the graphs sat inside the assemblers' own directories under three different names -
with none at all for MEGAHIT, which writes no graph and only leaves the per-k contigs a FASTG can be
built from. These tests pin the three things that makes usable: the flag is opt-in, the kept file
lands at a path a downstream tool can be told rather than have to glob for, and an assembler that
left no graph is reported instead of ending the run.

Note on what is *not* here: nothing in this module deletes an assembly graph, and these tests pin
that too (``test_a_run_without_keep_graph_keeps_and_deletes_nothing``). metaSPAdes' and metaFlye's
graphs have always survived in their own directories, as has MEGAHIT's ``intermediate_contigs``;
what was missing was a named path and, for MEGAHIT, a graph at all.
"""

from __future__ import annotations

import os
import shlex
from typing import List, Optional, Tuple

from metawrap2.modules import assembly

#: A MEGAHIT contig long enough to survive the module's 1 kb length filter, in MEGAHIT's own header
#: layout - the formatting step parses ``len=``/``multi=`` out of it.
MEGAHIT_CONTIG = ">k59_0 flag=1 multi=5.0000 len=2000\n%s\n" % ("ACGT" * 500)

#: The ks a real auto-k MEGAHIT run leaves behind. 59 is the largest, and is where final.contigs.fa
#: came from; a measured reference assembly stopped at exactly this k under a k_max of 141.
MEGAHIT_KS = (21, 39, 59)


def _reads(tmp_path) -> Tuple[str, str]:
    r1 = tmp_path / "s_1.fastq"
    r2 = tmp_path / "s_2.fastq"
    r1.write_text("@r\nACGT\n+\nIIII\n")
    r2.write_text("@r\nTGCA\n+\nIIII\n")
    return str(r1), str(r2)


def _fake_megahit_run(
    monkeypatch, tmp_path, intermediates: bool = True
) -> List[Tuple[str, Optional[str]]]:
    """Stand in for MEGAHIT, its toolkit and QUAST: record each call, write what it promises.

    The history file is redirected because ``main`` records a real run when it finishes, and a test
    has no business appending to the developer's ``~/.metawrap2/history.jsonl``.
    """
    calls: List[Tuple[str, Optional[str]]] = []

    def fake_run(cmd, env=None, tool=None, stdout_path=None, **kwargs):
        calls.append((cmd if isinstance(cmd, str) else " ".join(cmd), tool))
        parts = shlex.split(cmd) if isinstance(cmd, str) else list(cmd)
        if tool == "megahit":
            out = parts[parts.index("-o") + 1]
            os.makedirs(out, exist_ok=True)
            with open(os.path.join(out, "final.contigs.fa"), "w") as fh:
                fh.write(MEGAHIT_CONTIG)
            if intermediates:
                per_k = os.path.join(out, assembly.MEGAHIT_INTERMEDIATES)
                os.makedirs(per_k, exist_ok=True)
                for k in MEGAHIT_KS:
                    with open(os.path.join(per_k, "k%d.contigs.fa" % k), "w") as fh:
                        fh.write(MEGAHIT_CONTIG)
        elif tool == "megahit_toolkit contig2fastg":
            with open(stdout_path, "w") as fh:
                fh.write("FASTG\n")
        elif tool == "quast":
            out = parts[parts.index("-o") + 1]
            os.makedirs(out, exist_ok=True)
            with open(os.path.join(out, "report.html"), "w") as fh:
                fh.write("<html></html>\n")
        return 0

    monkeypatch.setattr(assembly, "run", fake_run)
    monkeypatch.setattr(assembly, "env_for", lambda module, settings: None)
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    return calls


def _assemble(
    tmp_path, monkeypatch, *extra, intermediates: bool = True
) -> Tuple[str, List[Tuple[str, Optional[str]]]]:
    r1, r2 = _reads(tmp_path)
    out = str(tmp_path / "asm")
    calls = _fake_megahit_run(monkeypatch, tmp_path, intermediates=intermediates)
    assert assembly.main(["-1", r1, "-2", r2, "-o", out, "-t", "2", "-m", "8", *extra]) == 0
    return out, calls


def test_keep_graph_is_off_by_default():
    """Opt-in, or every run pays for a graph most runs never read.

    Regresses: a default-on flag would add a FASTG conversion and (measured on an oral
    sample) 235 MB of disk to every single-sample MEGAHIT assembly in a cohort.
    """
    args = assembly._parse_args(["-1", "a_1.fq", "-2", "a_2.fq", "-o", "out"])
    assert args.keep_graph is False
    assert assembly._parse_args(["-L", "l.fq", "-o", "out", "--keep-graph"]).keep_graph is True


def test_a_run_without_keep_graph_keeps_and_deletes_nothing(tmp_path, monkeypatch):
    """Without the flag the assembler's output directory is exactly what it always was.

    Regresses both directions of the same change: no graph may appear at the kept path, and the
    material the graph would be built from - MEGAHIT's ``intermediate_contigs``, which this module
    has never deleted - must still be there untouched.
    """
    out, calls = _assemble(tmp_path, monkeypatch)
    assert not [name for name in os.listdir(out) if name.startswith("assembly_graph")]
    assert [tool for _cmd, tool in calls if tool == "megahit_toolkit contig2fastg"] == []
    intermediates = os.path.join(out, "megahit", assembly.MEGAHIT_INTERMEDIATES)
    assert sorted(os.listdir(intermediates)) == sorted("k%d.contigs.fa" % k for k in MEGAHIT_KS)


def test_keep_graph_names_the_kept_megahit_graph_and_its_size(tmp_path, monkeypatch, capsys):
    """The kept graph is at a fixed path and that path is in the run's own output.

    Regresses the whole point of the flag: a downstream tool has to be able to be *told* where the
    graph is, and a user has to be able to see what keeping it cost, without globbing an assembler's
    directory or knowing which assembler ran.
    """
    out, calls = _assemble(tmp_path, monkeypatch, "--keep-graph")
    graph = os.path.join(out, assembly.KEPT_GRAPHS["megahit"])
    assert os.path.getsize(graph)
    # Built from the largest k MEGAHIT actually wrote, not from options.json's k_max.
    cmd = next(cmd for cmd, tool in calls if tool == "megahit_toolkit contig2fastg")
    assert cmd == assembly.MEGAHIT_CONTIG2FASTG.format(
        k=59, contigs=os.path.join(out, "megahit", assembly.MEGAHIT_INTERMEDIATES, "k59.contigs.fa")
    )
    printed = capsys.readouterr().out
    assert assembly.KEPT_GRAPHS["megahit"] in printed
    assert "6.0 B" in printed  # the size of what was kept, not just that something was


def test_an_existing_graph_is_kept_without_a_second_copy_of_it(tmp_path, capsys):
    """metaSPAdes' and metaFlye's graphs are linked, not copied, so keeping one costs no disk.

    Regresses a copy creeping back in: the metaSPAdes GFA measured 121 MB against a 54 MB assembly,
    and a flag whose cost is "the graph twice" is a flag nobody can afford to use.
    """
    src = tmp_path / "metaspades" / assembly.SPADES_GRAPH
    src.parent.mkdir()
    src.write_text("H\tVN:Z:1.0\n")
    dest = tmp_path / assembly.KEPT_GRAPHS["metaspades"]

    assert assembly._keep_existing_graph(str(src), str(dest), "metaSPAdes") == [str(dest)]
    assert dest.stat().st_ino == src.stat().st_ino
    assert dest.stat().st_nlink == 2
    assert str(dest) in capsys.readouterr().out


def test_the_largest_k_megahit_wrote_is_the_one_converted(tmp_path):
    """MEGAHIT's auto-k stops early, so the k to convert is the largest one on disk.

    Regresses reading k_max out of options.json: the measured TEST leftover assembly ran k21..k59
    under a k_max of 141, so converting 141 would have read a file that does not exist. An empty
    k file is the same case - MEGAHIT touches one per k it starts and does not finish.
    """
    intermediates = tmp_path / assembly.MEGAHIT_INTERMEDIATES
    intermediates.mkdir()
    for k in MEGAHIT_KS:
        (intermediates / ("k%d.contigs.fa" % k)).write_text(MEGAHIT_CONTIG)
    (intermediates / "k141.contigs.fa").write_text("")
    (intermediates / "k59.addi.fa").write_text(MEGAHIT_CONTIG)

    assert assembly._megahit_graph_k(str(intermediates)) == 59
    assert assembly._megahit_graph_k(str(tmp_path / "not_there")) is None


def test_an_assembler_that_left_no_graph_is_reported_not_fatal(tmp_path, monkeypatch, capsys):
    """A missing graph warns and the run finishes; it does not throw the assembly away.

    Regresses turning an optional extra into a run-ending failure. By the time this runs the
    assembly is finished and correct, and hours of it are on disk - the graph is a convenience that
    nothing downstream has consumed yet, so "there is no graph" is news, not a failure.
    """
    # A MEGAHIT that produced contigs and left no per-k intermediates: nothing to convert.
    out, calls = _assemble(tmp_path, monkeypatch, "--keep-graph", intermediates=False)
    assert not os.path.exists(os.path.join(out, assembly.KEPT_GRAPHS["megahit"]))
    printed = capsys.readouterr().out
    assert "no k*.contigs.fa" in printed
    assert os.path.getsize(os.path.join(out, "final_assembly.fasta"))
    assert [tool for _cmd, tool in calls if tool == "quast"]  # the run carried on to the end
