"""Gene calling over a bin set or a whole assembly, done once and shared.

All functional analyses downstream (genome content, KEGG module completeness, resistance genes and
their neighbourhoods, abundance) start from the same gene list. ``annotate_bins`` used to call
genes separately per analysis, which duplicated compute and produced gene sets that could not be
joined when two calls differed.

This is a library, not a module. Genes are called with Pyrodigal (Cython bindings to Prodigal
2.6.3, so the calls are Prodigal's) and returned as data with two fields needed for fragmented
assemblies:

* ``partial``: Prodigal's flag for a gene whose start or stop codon was not found because the
  contig ended. The gene is real; its boundary is unknown.
* ``distance_to_contig_end``: needed to interpret a negative result. "No mobile element within
  10 kb of this resistance gene" says nothing if the gene is 2 kb from its contig end. Computed
  here once for every consumer that reports a neighbourhood.

Single-genome or metagenomic mode is chosen per input. Prodigal trains its coding model on the
input, which works best on one organism. Each bin trains its own model; a whole assembly uses
Prodigal's pre-trained metagenomic models. Below Prodigal's 20 kb training floor
(:data:`pyrodigal.MIN_SINGLE_GENOME`) a bin falls back to meta mode, and the mode is recorded.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import pyrodigal

from .io.seqio import contig_id, iter_fasta, smart_open

#: Translation table used to train a single-genome coding model: 11 (bacterial, archaeal, plant
#: plastid), Prodigal's default. A model trained under a different table finds different genes.
#:
#: Proteins are not always translated with this table. Some of Prodigal's pre-trained metagenomic
#: models use table 4 (TGA codes for tryptophan), and translating those genes with table 11 would
#: truncate them at the first TGA. Each gene is translated with the table it was called under, which
#: is recorded per gene.
TRAINING_TABLE = 11

#: Minimum sequence for training a bin's own coding model instead of using Prodigal's metagenomic
#: models. Prodigal's floor; below it a trained model is worse.
MIN_TRAIN_BP = pyrodigal.MIN_SINGLE_GENOME

#: Written into the GFF header so a file says which caller and which mode produced it.
SOURCE = "Pyrodigal_%s" % pyrodigal.PRODIGAL_VERSION


@dataclass(frozen=True)
class Gene:
    """One predicted gene, with the provenance a fragmented assembly needs.

    Coordinates are 1-based and inclusive, as in GFF and Prodigal output.
    """

    gene_id: str
    contig: str
    start: int
    end: int
    strand: int
    #: Prodigal's two-character flag: "00" complete, "10" missing a start, "01" missing a stop,
    #: "11" running off both ends of its contig.
    partial: str
    #: Bases between this gene and the nearer end of its contig. 0 for a gene that reaches the end.
    distance_to_contig_end: int
    contig_length: int
    #: The genetic code this gene was called and translated under. Usually 11; some of Prodigal's
    #: metagenomic models use table 4.
    translation_table: int
    protein: str
    nucleotide: str

    @property
    def is_partial(self) -> bool:
        return self.partial != "00"

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1


def _partial_flag(gene: pyrodigal.Gene) -> str:
    """Prodigal's partial flag, in Prodigal's own order: begin then end, in genome orientation."""
    return "%d%d" % (int(gene.partial_begin), int(gene.partial_end))


def _genes_of_contig(finder: pyrodigal.GeneFinder, name: str, sequence: str) -> Iterator[Gene]:
    contig_length = len(sequence)
    for index, found in enumerate(finder.find_genes(sequence), start=1):
        # Distance to the nearer end, which limits how much neighbourhood can be seen. Both ends
        # can be recovered from start/end and contig_length.
        distance = min(found.begin - 1, contig_length - found.end)
        yield Gene(
            gene_id="%s_%d" % (name, index),
            contig=name,
            start=found.begin,
            end=found.end,
            strand=found.strand,
            partial=_partial_flag(found),
            distance_to_contig_end=distance,
            contig_length=contig_length,
            translation_table=found.translation_table,
            # No translation_table= argument: translate under the code the gene was called with.
            protein=found.translate(),
            nucleotide=found.sequence(),
        )


def _trained_finder(sequences: Sequence[str]) -> Tuple[pyrodigal.GeneFinder, str]:
    """A finder for one organism, trained on it, or the metagenomic one if there is too little.

    Returns the finder and its mode. Meta mode gives systematically different gene boundaries, so
    the mode is recorded for comparing gene counts across bins.
    """
    total = sum(len(s) for s in sequences)
    if total < MIN_TRAIN_BP:
        return pyrodigal.GeneFinder(meta=True), "meta"
    finder = pyrodigal.GeneFinder(meta=False)
    # pyrodigal's overloads accept the literal; the tests assert TRAINING_TABLE == 11.
    finder.train(*sequences, translation_table=11)
    return finder, "single"


def _covered_bases(genes: Sequence[Gene]) -> int:
    """Bases inside at least one gene, counting each base once.

    Computed per contig, since coordinates restart on each contig.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for gene in genes:
        by_contig.setdefault(gene.contig, []).append((gene.start, gene.end))
    covered = 0
    for spans in by_contig.values():
        spans.sort()
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                covered += current_end - current_start + 1
                current_start, current_end = start, end
        covered += current_end - current_start + 1
    return covered


@dataclass(frozen=True)
class GenomeGenes:
    """The genes of one genome (or one assembly), with how they were called."""

    name: str
    genes: List[Gene]
    #: "single" if a coding model was trained on this sequence, "meta" if Prodigal's pre-trained
    #: metagenomic bins were used instead.
    mode: str
    n_contigs: int
    total_bp: int

    @property
    def n_partial(self) -> int:
        return sum(1 for g in self.genes if g.is_partial)

    @property
    def coding_density(self) -> float:
        """Fraction of the genome inside a predicted gene, as a sanity check.

        Bacterial genomes are about 0.85-0.90. Much lower values usually mean a viral contig with
        overlapping reading frames or a bin of very short fragments.

        Computed over the union of gene intervals. Prodigal allows up to ``max_overlap`` bases of
        overlap between neighbours, so summing lengths could exceed 1.0.
        """
        if not self.total_bp:
            return 0.0
        return _covered_bases(self.genes) / float(self.total_bp)


def call_genome(path: str, name: Optional[str] = None, meta: bool = False) -> GenomeGenes:
    """Call genes in one FASTA.

    *meta* forces metagenomic mode; otherwise a coding model is trained on the file when it has
    enough sequence. Pass ``meta=True`` for a whole assembly, which contains many organisms.
    """
    records = [(contig_id(header), seq) for header, seq in iter_fasta(path)]
    sequences = [seq for _name, seq in records]
    if meta:
        finder, mode = pyrodigal.GeneFinder(meta=True), "meta"
    else:
        finder, mode = _trained_finder(sequences)

    genes: List[Gene] = []
    for contig, sequence in records:
        genes.extend(_genes_of_contig(finder, contig, sequence))
    return GenomeGenes(
        name=name if name is not None else os.path.splitext(os.path.basename(path))[0],
        genes=genes,
        mode=mode,
        n_contigs=len(records),
        total_bp=sum(len(s) for s in sequences),
    )


# ── writing ──────────────────────────────────────────────────────────────────────────────────


def write_faa(called: GenomeGenes, path: str) -> int:
    """Write translated genes, with headers carrying the fields needed for downstream joins.

    ``>gene_id # start # end # strand # partial=00;dist_to_end=1234`` keeps Prodigal's
    ``#``-delimited format, which eggNOG-mapper, PyHMMER and other consumers treat as a
    description, and adds the partial flag and distance. Parsers that read only the first token are
    unaffected.
    """
    with smart_open(path, "wt") as out:
        for gene in called.genes:
            out.write(
                ">%s # %d # %d # %d # partial=%s;dist_to_end=%d\n"
                % (
                    gene.gene_id,
                    gene.start,
                    gene.end,
                    gene.strand,
                    gene.partial,
                    gene.distance_to_contig_end,
                )
            )
            out.write(gene.protein + "\n")
    return len(called.genes)


def write_ffn(called: GenomeGenes, path: str) -> int:
    """Untranslated genes, with the same headers as the .faa so the two files join on gene_id."""
    with smart_open(path, "wt") as out:
        for gene in called.genes:
            out.write(
                ">%s # %d # %d # %d # partial=%s;dist_to_end=%d\n"
                % (
                    gene.gene_id,
                    gene.start,
                    gene.end,
                    gene.strand,
                    gene.partial,
                    gene.distance_to_contig_end,
                )
            )
            out.write(gene.nucleotide + "\n")
    return len(called.genes)


def write_gff(called: GenomeGenes, path: str) -> int:
    """GFF3 of the calls, with partial flag and distance as attributes."""
    with smart_open(path, "wt") as out:
        out.write("##gff-version 3\n")
        out.write("# source: %s (%s mode)\n" % (SOURCE, called.mode))
        for gene in called.genes:
            out.write(
                "\t".join(
                    (
                        gene.contig,
                        SOURCE,
                        "CDS",
                        str(gene.start),
                        str(gene.end),
                        ".",
                        "+" if gene.strand > 0 else "-",
                        "0",
                        "ID=%s;partial=%s;dist_to_contig_end=%d"
                        % (gene.gene_id, gene.partial, gene.distance_to_contig_end),
                    )
                )
                + "\n"
            )
    return len(called.genes)


#: Columns of the per-gene table, in order. One row per gene per genome; later joins (KOs, Pfams,
#: ARG calls) can be many-to-one per gene.
GENE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "length_bp",
    "partial",
    "distance_to_contig_end",
    "contig_length",
    "translation_table",
    "calling_mode",
)


def gene_rows(called: GenomeGenes) -> Iterator[Dict[str, object]]:
    for gene in called.genes:
        yield {
            "genome": called.name,
            "gene_id": gene.gene_id,
            "contig": gene.contig,
            "start": gene.start,
            "end": gene.end,
            "strand": "+" if gene.strand > 0 else "-",
            "length_bp": gene.length_bp,
            "partial": gene.partial,
            "distance_to_contig_end": gene.distance_to_contig_end,
            "contig_length": gene.contig_length,
            "translation_table": gene.translation_table,
            "calling_mode": called.mode,
        }
