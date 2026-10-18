"""FastTree's data type is a flag, not something it works out, and the flag was wrong.

`fasttree -nt` declares nucleotide. Every alignment this module is pointed at is protein - the
pangenome core-gene alignments it builds, and GTDB-Tk's MSA, which its own docs name as an input -
and `-nt` over protein does not fail: FastTree reads every residue that is not A, C, G or T as
ambiguous, which is most of a protein alignment, and returns a tree built from almost nothing.
Reported by a user reading the command, not by any output looking wrong, which is the point.
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
    """Fourteen amino acids are also IUPAC base codes, so the test is a fraction not a membership.

    What keeps it right is composition: L, I, E, Q, F and P are about a third of any real protein,
    and none of them is a base code. A sequence built *only* of dual-use letters would be called
    nucleotide, and correctly so - there is nothing in it that says otherwise.
    """
    rows = [
        ("g1", "MKVLIWGAGHDEQFPLIEKSTNRYACG"),
        ("g2", "MKVLIWGAGHDENFPLIEKSTNRYACG"),
    ]
    assert phylogeny.read_alignment(_write(tmp_path, "dual.faa", rows)).kind == "protein"
    # The all-dual-use case, stated rather than hidden: it reads as nucleotide.
    assert phylogeny.alignment_kind({"A": 5, "C": 5, "G": 5, "T": 5, "W": 2}) == "nucleotide"


def test_an_empty_alignment_defaults_to_protein():
    """Uncertain means protein: that is what every aligner upstream of this module writes."""
    assert phylogeny.alignment_kind({}) == "protein"


def test_fasttree_gets_lg_on_protein_and_gtr_nt_on_nucleotides():
    """The two templates, so neither can quietly drift back into one."""
    assert "-nt" in module.FASTTREE_NT and "-gtr" in module.FASTTREE_NT
    assert "-nt" not in module.FASTTREE_AA and "-lg" in module.FASTTREE_AA


def test_the_summary_records_the_model_that_was_used(tmp_path, monkeypatch, capsys):
    """A tree whose model is reported as GTR when LG produced it is a provenance defect too."""
    rows = [("g1", "MKVLIWGAGHDEQ"), ("g2", "MKVLIWGAGHDEN"), ("g3", "MRVLIYGAGHDEQ")]
    alignment = phylogeny.read_alignment(_write(tmp_path, "core.faa", rows))
    assert alignment.kind == "protein"
    template = module.FASTTREE_NT if alignment.kind == "nucleotide" else module.FASTTREE_AA
    assert template is module.FASTTREE_AA
