"""FastTree data type: the module passed `-nt` for protein alignments.

`fasttree -nt` declares nucleotide input. The module's inputs are protein (its pangenome core-gene
alignments and GTDB-Tk's MSA). With `-nt` FastTree treats every residue other than A, C, G or T
as ambiguous and returns a tree from almost no data, without an error. A user found this by
reading the command.
"""

from metawrap2 import phylogeny
from metawrap2.modules import phylogeny as module


def _write(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text("".join(">%s\n%s\n" % (taxon, seq) for taxon, seq in rows))
    return str(path)


def test_a_protein_alignment_is_called_protein(tmp_path):
    rows = [("g1", "MKVLIWGAGH-DEQ"), ("g2", "MKVLIWGAGH-DEN"), ("g3", "MRVLIYGAGH-DEQ")]
    assert phylogeny.read_alignment(_write(tmp_path, "core.faa", rows)).kind == "protein"


def test_a_nucleotide_alignment_is_called_nucleotide(tmp_path):
    rows = [("g1", "ACGTACGT--ACGT"), ("g2", "ACGTACGTTTACGT"), ("g3", "ACGTNCGT--ACGA")]
    assert phylogeny.read_alignment(_write(tmp_path, "core.fna", rows)).kind == "nucleotide"


def test_dual_use_letters_do_not_flip_a_real_protein_alignment(tmp_path):
    """Fourteen amino-acid letters are also IUPAC base codes, so detection uses a fraction.

    L, I, E, Q, F and P make up about a third of a real protein and are not base codes. A sequence
    made only of dual-use letters is called nucleotide.
    """
    rows = [
        ("g1", "MKVLIWGAGHDEQFPLIEKSTNRYACG"),
        ("g2", "MKVLIWGAGHDENFPLIEKSTNRYACG"),
    ]
    assert phylogeny.read_alignment(_write(tmp_path, "dual.faa", rows)).kind == "protein"
    # A sequence of only dual-use letters reads as nucleotide.
    assert phylogeny.alignment_kind({"A": 5, "C": 5, "G": 5, "T": 5, "W": 2}) == "nucleotide"


def test_an_empty_alignment_defaults_to_protein():
    """An empty alignment defaults to protein, which every upstream aligner writes."""
    assert phylogeny.alignment_kind({}) == "protein"


def test_fasttree_gets_lg_on_protein_and_gtr_nt_on_nucleotides():
    """Protein uses LG and nucleotide uses GTR."""
    assert "-nt" in module.FASTTREE_NT and "-gtr" in module.FASTTREE_NT
    assert "-nt" not in module.FASTTREE_AA and "-lg" in module.FASTTREE_AA


def test_the_summary_records_the_model_that_was_used(tmp_path, monkeypatch, capsys):
    """The summary records the substitution model used."""
    rows = [("g1", "MKVLIWGAGHDEQ"), ("g2", "MKVLIWGAGHDEN"), ("g3", "MRVLIYGAGHDEQ")]
    alignment = phylogeny.read_alignment(_write(tmp_path, "core.faa", rows))
    assert alignment.kind == "protein"
    template = module.FASTTREE_NT if alignment.kind == "nucleotide" else module.FASTTREE_AA
    assert template is module.FASTTREE_AA
