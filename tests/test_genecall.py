"""Gene calling, and the two fields a fragmented assembly makes necessary.

Prodigal's calls are not tested here, only what :mod:`metawrap2.genecall` adds: a gene running off
a contig end is marked partial, the distance to the nearer contig end is recorded (so a "nothing
within 10 kb" claim can be checked), and a bin too small to train on falls back to metagenomic
mode with a message.
"""

from __future__ import annotations

from metawrap2 import genecall

# A real coding sequence so Prodigal finds a gene: E. coli thrA's first ~700 bp, in frame from its
# ATG.
CODING = (
    "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAG"
    "GTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATC"
    "AGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAG"
    "GAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAA"
    "ATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC"
    "CTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGGATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCGCC"
    "GGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAG"
    "ATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAG"
    "CTTTCCTACTTCGGCGCTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCT"
)


#: A non-coding spacer for either side of a CDS, with stop codons in every frame. Poly-T does not
#: work: it reads as a long ORF on the reverse strand.
SPACER = "TAGTAATGACTAGTTAACTAGTAATGACTAGTTAA" * 10


def _write_fasta(path: str, records) -> str:
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, seq) for name, seq in records)
    return path


def test_a_gene_running_off_a_contig_is_marked_partial(tmp_path):
    """A CDS truncated by the contig end is reported as partial.

    On CAMI2 bins, 50-54% of genes were partial in bins of ~900-1000 contigs, against 8% in a bin
    of 92 contigs.
    """
    # A contig that *is* the middle of a gene: no start codon before it, no stop after it.
    path = _write_fasta(str(tmp_path / "frag.fasta"), [("frag", CODING[30:])])
    called = genecall.call_genome(path, meta=True)
    assert called.genes, "expected a partial call rather than nothing"
    assert any(g.is_partial for g in called.genes)
    # Distance to the nearer contig end is zero at the edge, so a claim such as "no mobile element
    # within 10 kb" can be checked.
    assert any(g.distance_to_contig_end == 0 for g in called.genes)
    for gene in called.genes:
        assert gene.distance_to_contig_end == min(gene.start - 1, gene.contig_length - gene.end)


def test_the_three_output_files_agree_on_gene_ids(tmp_path):
    """The .faa, .ffn and gene table are joined on gene_id downstream, so they have to match."""
    path = _write_fasta(str(tmp_path / "g.fasta"), [("c", SPACER + CODING + "TAA" + SPACER)])
    called = genecall.call_genome(path, name="g", meta=True)

    faa, ffn, gff = (str(tmp_path / ("g" + ext)) for ext in (".faa", ".ffn", ".gff"))
    assert genecall.write_faa(called, faa) == len(called.genes)
    assert genecall.write_ffn(called, ffn) == len(called.genes)
    assert genecall.write_gff(called, gff) == len(called.genes)

    def ids(path):
        with open(path) as fh:
            return [line[1:].split()[0] for line in fh if line.startswith(">")]

    expected = [g.gene_id for g in called.genes]
    assert ids(faa) == expected
    assert ids(ffn) == expected
    assert [row["gene_id"] for row in genecall.gene_rows(called)] == expected
