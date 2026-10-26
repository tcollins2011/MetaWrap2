"""``assembly --keep-graph``: keeping the assembly graph for graph-aware binners, and its cost.

PlasMAAG, GraphMB and Hi-C placement read the assembly graph. metaSPAdes and metaFlye write graphs
under different names in their own directories; MEGAHIT writes none, only the per-k contigs a FASTG
can be built from. Tested: the flag is opt-in, the graph is kept at a fixed path, and a missing
graph is reported without ending the run.

The module never deletes an assembly graph or MEGAHIT's ``intermediate_contigs``
(``test_a_run_without_keep_graph_keeps_and_deletes_nothing``).
"""

from __future__ import annotations

import os
import shlex
from typing import List, Optional, Tuple

from metawrap2.modules import assembly

#: A MEGAHIT contig long enough to pass the module's 1 kb length filter, in MEGAHIT's header
#: layout (the formatting step parses ``len=``/``multi=``).
MEGAHIT_CONTIG = ">k59_0 flag=1 multi=5.0000 len=2000\n%s\n" % ("ACGT" * 500)

#: The ks a real auto-k MEGAHIT run leaves. 59 is the largest and produced final.contigs.fa; a
#: reference assembly stopped at this k with a k_max of 141.
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
    """Fake MEGAHIT, its toolkit and QUAST: record each call and write its expected outputs.

    The history file is redirected so the test does not append to ``~/.metawrap2/history.jsonl``.
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
    """--keep-graph is off by default.

    When on, it adds a FASTG conversion and (on an oral sample) 235 MB of disk per MEGAHIT
    assembly.
    """
    args = assembly._parse_args(["-1", "a_1.fq", "-2", "a_2.fq", "-o", "out"])
    assert args.keep_graph is False
    assert assembly._parse_args(["-L", "l.fq", "-o", "out", "--keep-graph"]).keep_graph is True


def test_a_run_without_keep_graph_keeps_and_deletes_nothing(tmp_path, monkeypatch):
    """Without the flag, no graph is written and MEGAHIT's ``intermediate_contigs`` are kept."""
    out, calls = _assemble(tmp_path, monkeypatch)
    assert not [name for name in os.listdir(out) if name.startswith("assembly_graph")]
    assert [tool for _cmd, tool in calls if tool == "megahit_toolkit contig2fastg"] == []
    intermediates = os.path.join(out, "megahit", assembly.MEGAHIT_INTERMEDIATES)
    assert sorted(os.listdir(intermediates)) == sorted("k%d.contigs.fa" % k for k in MEGAHIT_KS)


def test_keep_graph_names_the_kept_megahit_graph_and_its_size(tmp_path, monkeypatch, capsys):
    """The kept graph is at a fixed path, and the path and its size are printed."""
    out, calls = _assemble(tmp_path, monkeypatch, "--keep-graph")
    graph = os.path.join(out, assembly.KEPT_GRAPHS["megahit"])
    assert os.path.getsize(graph)
    # Built from the largest k MEGAHIT wrote, not options.json's k_max.
    cmd = next(cmd for cmd, tool in calls if tool == "megahit_toolkit contig2fastg")
    assert cmd == assembly.MEGAHIT_CONTIG2FASTG.format(
        k=59, contigs=os.path.join(out, "megahit", assembly.MEGAHIT_INTERMEDIATES, "k59.contigs.fa")
    )
    printed = capsys.readouterr().out
    assert assembly.KEPT_GRAPHS["megahit"] in printed
    assert "6.0 B" in printed  # the size of the kept graph


def test_an_existing_graph_is_kept_without_a_second_copy_of_it(tmp_path, capsys):
    """metaSPAdes' and metaFlye's graphs are hard-linked, not copied.

    A metaSPAdes GFA measured 121 MB for a 54 MB assembly.
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
    """The k converted is the largest non-empty one on disk, not k_max from options.json.

    MEGAHIT's auto-k can stop early (k21..k59 with a k_max of 141), and it leaves an empty file
    for a k it started but did not finish.
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
    """A missing graph gives a warning and the run finishes with its assembly."""
    # MEGAHIT produced contigs but no per-k intermediates, so there is nothing to convert.
    out, calls = _assemble(tmp_path, monkeypatch, "--keep-graph", intermediates=False)
    assert not os.path.exists(os.path.join(out, assembly.KEPT_GRAPHS["megahit"]))
    printed = capsys.readouterr().out
    assert "no k*.contigs.fa" in printed
    assert os.path.getsize(os.path.join(out, "final_assembly.fasta"))
    assert [tool for _cmd, tool in calls if tool == "quast"]  # the run continued to the end
