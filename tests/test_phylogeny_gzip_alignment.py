"""A gzipped alignment is decompressed before the tree builder sees it.

GTDB-Tk writes ``*.user_msa.fasta.gz``. FastTree and IQ-TREE read only plain text, and were given
the .gz, so ``phylogeny -d <GTDB-Tk folder>`` never built a tree.
"""

from __future__ import annotations

import gzip
import shlex

import pytest

from metawrap2 import command
from metawrap2.modules import phylogeny as module

ROWS = [
    ("g1", "MKVLIWGAGHDEQFPL"),
    ("g2", "MKVLIWGAGHDENFPL"),
    ("g3", "MRVLIYGAGHDEQFPI"),
    ("g4", "MRVLIYGSGHDEQFPI"),
]


@pytest.fixture(autouse=True)
def _history(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))


def _fake_fasttree(monkeypatch, seen):
    def fake_run(cmd, env=None, tool=None, stdout_path=None, **kwargs):
        alignment = shlex.split(cmd)[-1]
        with open(alignment) as fh:  # text mode, as the tool reads it
            seen.append((alignment, fh.read()))
        with open(stdout_path, "w") as fh:
            fh.write("((g1:0.1,g2:0.1)0.95:0.2,g3:0.1,g4:0.1);\n")
        return 0

    monkeypatch.setattr(module, "run", fake_run)
    monkeypatch.setattr(module, "env_for", lambda name, settings: None)


def test_a_gzipped_alignment_is_given_to_fasttree_decompressed(tmp_path, monkeypatch):
    align_dir = tmp_path / "align"
    align_dir.mkdir()
    gz = align_dir / "gtdbtk.bac120.user_msa.fasta.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write("".join(">%s\n%s\n" % row for row in ROWS))
    seen = []
    _fake_fasttree(monkeypatch, seen)
    out = tmp_path / "out"
    command.runner.configure(dry_run=False)
    assert module.main(["-a", str(gz), "-o", str(out), "--fast", "-t", "2"]) == 0
    ((path, text),) = seen
    assert not path.endswith(".gz")
    assert path.startswith(str(out / module.ALIGNMENTS_DIR))
    assert text.startswith(">g1\n")
    assert (out / module.SUMMARY_TABLE).is_file()


def test_a_plain_alignment_is_used_where_it_is(tmp_path, monkeypatch):
    plain = tmp_path / "core.faa"
    plain.write_text("".join(">%s\n%s\n" % row for row in ROWS))
    seen = []
    _fake_fasttree(monkeypatch, seen)
    out = tmp_path / "out"
    assert module.main(["-a", str(plain), "-o", str(out), "--fast", "-t", "2"]) == 0
    ((path, _text),) = seen
    assert path == str(plain)
    assert not (out / module.ALIGNMENTS_DIR).exists()
