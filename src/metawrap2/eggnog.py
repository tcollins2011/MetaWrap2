"""Orthologous groups from eggNOG, for the COG category, GO terms and EC numbers nothing else here
produces.

The gap this fills was measured in 2.1.2 rather than assumed. Of 3,891 genes in one bin, **1,578
carry a Pfam domain and no KO, while only 13 carry a KO and no Pfam**: Pfam nearly subsumes KEGG's
gene set, so running both buys 13 genes and a second opinion. The remaining value of KEGG is not
coverage at all - it is that its call is a *reaction*, which a domain name is not. What is actually
missing is the annotation that is orthogonal to both, and Pfam produces none of it: no GO term, no
COG functional category, no EC number.

eggNOG gives all three from one search, because it annotates by *orthology* rather than by profile
match. A gene is placed in an orthologous group, and the group's members carry GO, EC, COG and KEGG
annotations transferred from characterised proteins. That is a different inference from "this
protein contains a kinase domain", and it fails differently, which is the point of having it beside
the other two rather than instead of them.

**eggNOG's KEGG orthologs are not KOfam's, and are not merged with them.** eggNOG transfers a KO
from an orthologous group; KOfam calls one from a per-profile curated bit-score threshold. They
disagree, and which is right varies by gene. :func:`compare_kos` reports the three counts - called
by both, by eggNOG alone, by KOfam alone - and this module never prefers one. Picking a winner
silently would destroy the only measurement that justifies running both.

This file does not download anything. The eggNOG data (~50 GB: ``eggnog.db`` and
``eggnog_proteins.dmnd``) is a database like every other one here - named on the command line or in
``metawrap2.toml``, never fetched behind a user's back.

Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2:
functional annotation, orthology assignments, and domain prediction at the metagenomic scale.
*Mol Biol Evol* 38(12):5825-5829 (2021). doi:10.1093/molbev/msab293
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

#: eggNOG-mapper writes a literal ``-`` into every field it has no value for, and its row writer
#: substitutes ``-`` for anything empty or whitespace-only as well - so *every* column can be this,
#: not only the annotation ones. Parsed naively it becomes a GO term called "-", an EC number "-"
#: and a COG category "-", each of which counts as an annotation in any downstream tally. The same
#: defect was found and fixed in dbCAN's output, where ``-|-`` was being reported as the EC number
#: of genes that have none; :func:`_split` and :func:`_scalar` are the only way a value gets out of
#: a row here, so the sentinel cannot leak past them.
EMPTY = "-"

#: The columns ``emapper.annotations`` carries, in order. Taken from eggNOG-mapper 2.1.13's own
#: writer (``eggnogmapper/annotation/output.py``: ``HIT_HEADER + ANNOTATIONS_HEADER``) rather than
#: from documentation, because the column set is what a parser is pinned to and the docs lag it.
ANNOTATION_COLUMNS: Tuple[str, ...] = (
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)

#: ``--md5`` appends one more column. Accepted so a user who passed it is not told their file is
#: malformed, but never required.
MD5_COLUMN = "md5"

#: emapper glues a ``#`` onto the first column name, so the header line reads ``#query<TAB>...``.
#: It also writes ``## `` provenance lines above the header (date, version, argv) and ``## N
#: queries scanned`` footer lines below the last row. Those footers are tab-free single fields and
#: would parse as a data row with a gene called "## 3 queries scanned" if comment lines were not
#: skipped at both ends of the file.
COMMENT = "#"

#: Fields emapper joins with commas, from a set, so they arrive already sorted. Split and
#: deduplicated here so a count of GO terms counts terms rather than strings, and the order is left
#: as emapper wrote it so a re-run over the same file writes a byte-identical table.
MULTI_VALUED: Tuple[str, ...] = (
    "eggNOG_OGs",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "BRITE",
    "PFAMs",
)

#: ``KEGG_ko`` arrives namespace-tagged: ``ko:K00059``, not ``K00059``. Every KO table in this
#: package - ``annotate --kofam``'s ``kos.tsv``, ``pathways``' BRITE hierarchy - is keyed on the
#: bare K number, so an unstripped prefix makes the join that is the entire reason for running both
#: annotations return zero rows while looking like a legitimate disagreement.
KO_PREFIX = "ko:"

#: Contrast with ``KO_PREFIX``, because the two look alike and only one is strippable. ``GO:`` is
#: part of a GO accession - ``GO:0003674`` *is* the term's name - and ``KEGG_Pathway``'s ``ko00010``
#: carries no colon at all. So exactly one prefix is removed, and it is removed from exactly one
#: column.
GO_PREFIX = "GO:"


def _scalar(value: Optional[str]) -> str:
    """One field as text, with eggNOG's absent-value sentinel resolved to empty."""
    text = (value or "").strip()
    return "" if text == EMPTY else text


def _split(value: Optional[str]) -> Tuple[str, ...]:
    """One comma-joined field as its members, sentinels and empties dropped, order preserved.

    Order is emapper's own (it joins a sorted set), not re-sorted here: a parser that imposes its
    own order makes its output depend on its own version rather than on the input, and byte-identity
    across re-runs is what lets a table be diffed.
    """
    text = _scalar(value)
    if not text:
        return ()
    members = [piece.strip() for piece in text.split(",")]
    return tuple(dict.fromkeys(piece for piece in members if piece and piece != EMPTY))


def _number(value: Optional[str]) -> Optional[float]:
    """An e-value or bit score, or None.

    None rather than 0.0 for an unparseable field: emapper's ``-m no_search`` and cache modes can
    leave the hit columns as ``-``, and a 0.0 e-value there would read as a perfect hit.
    """
    text = _scalar(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


#: COG functional categories for one gene arrive *concatenated*, not comma-separated: a gene in
#: three categories is the single token ``EGP``. Splitting this column on commas yields one
#: "category" called EGP, which belongs to no category and silently drops two real ones, so it is
#: deliberately absent from MULTI_VALUED and split per character instead.
COG_CATEGORY_WIDTH = 1


@dataclass(frozen=True)
class EggnogRow:
    """One gene's eggNOG annotation, as one line of ``emapper.annotations`` describes it."""

    genome: str
    gene: str
    #: The best-matching eggNOG protein, ``taxid.protein``. Kept because every annotation in this
    #: row was transferred from this hit's orthologous groups, so it is the row's evidence.
    seed_ortholog: str
    evalue: Optional[float]
    score: Optional[float]
    #: Orthologous groups, finest first, each ``OG@taxid|taxon``.
    og_names: Tuple[str, ...]
    #: The taxonomic level the annotations were transferred at. A phylum-level transfer and a
    #: genus-level one are not equally strong evidence, which is why it travels with the row.
    max_annot_lvl: str
    #: Concatenated one-letter categories, as written - see COG_CATEGORY_WIDTH.
    cog_category: str
    description: str
    preferred_name: str
    #: GO terms, ``GO:`` prefix intact because it is part of the accession.
    gos: Tuple[str, ...]
    ecs: Tuple[str, ...]
    #: Bare K numbers, ``ko:`` already stripped - see KO_PREFIX.
    kegg_kos: Tuple[str, ...]
    kegg_pathways: Tuple[str, ...]
    brites: Tuple[str, ...]
    pfams: Tuple[str, ...]

    @property
    def cog_categories(self) -> Tuple[str, ...]:
        """``"EGP"`` -> ``("E", "G", "P")``; ``""`` -> ``()``.

        A gene with no COG category yields nothing rather than a category named "-", which is what
        the sentinel becomes if this column is counted as written.
        """
        return tuple(
            self.cog_category[i : i + COG_CATEGORY_WIDTH]
            for i in range(0, len(self.cog_category), COG_CATEGORY_WIDTH)
        )

    @property
    def has_orthogonal_annotation(self) -> bool:
        """Whether this row carries anything KOfam and Pfam cannot produce.

        The three things this module exists for. A row with a seed ortholog, a description and
        nothing else is a hit without an annotation, and counting it as one overstates what eggNOG
        added over the searches that already ran.
        """
        return bool(self.gos or self.ecs or self.cog_category)


class AnnotationFormatError(ValueError):
    """``emapper.annotations`` did not have the columns this parser is pinned to.

    Its own class so a caller can tell a schema change from an IO error, and raised rather than
    skipped: a silent zero-row table from a renamed column looks exactly like a genome with no
    annotations, and the two want opposite responses.
    """


def _check_header(line: str, path: str) -> Tuple[str, ...]:
    """Validate emapper's header line and return the column names it declares.

    Names what was found rather than failing later on a KeyError, because the failure this catches
    is an eggNOG-mapper version whose columns differ from 2.1.13's - and the useful message there is
    which columns moved, not which row could not be read.
    """
    found = tuple(line.lstrip(COMMENT).rstrip("\n").split("\t"))
    expected = ANNOTATION_COLUMNS
    if found == expected or found == expected + (MD5_COLUMN,):
        return found
    missing = [name for name in expected if name not in found]
    extra = [name for name in found if name not in expected + (MD5_COLUMN,)]
    raise AnnotationFormatError(
        "%s does not look like eggNOG-mapper %d-column output: found %d column(s)%s%s"
        % (
            path,
            len(expected),
            len(found),
            "; missing %s" % ", ".join(missing) if missing else "",
            "; unexpected %s" % ", ".join(extra) if extra else "",
        )
    )


def parse_annotations(path: str, genome: str) -> List[EggnogRow]:
    """Read one ``emapper.annotations`` file.

    Only genes emapper placed appear in it, so its length is the annotated-gene count and not the
    gene count - the share of a genome eggNOG annotated has to be computed against the gene table,
    the way the Pfam and CAZyme passes do.
    """
    rows: List[EggnogRow] = []
    columns: Optional[Tuple[str, ...]] = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(COMMENT):
                # The header is the only comment line that is not a comment. Everything else
                # starting with "#" is emapper's own provenance block or its footer counts.
                if columns is None and line.startswith(COMMENT + ANNOTATION_COLUMNS[0]):
                    columns = _check_header(line, path)
                continue
            if columns is None:
                raise AnnotationFormatError(
                    "%s has data rows before any %s%s header line"
                    % (path, COMMENT, ANNOTATION_COLUMNS[0])
                )
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(ANNOTATION_COLUMNS):
                continue
            cell = dict(zip(columns, fields))
            rows.append(
                EggnogRow(
                    genome=genome,
                    gene=_scalar(cell.get("query")),
                    seed_ortholog=_scalar(cell.get("seed_ortholog")),
                    evalue=_number(cell.get("evalue")),
                    score=_number(cell.get("score")),
                    og_names=_split(cell.get("eggNOG_OGs")),
                    max_annot_lvl=_scalar(cell.get("max_annot_lvl")),
                    cog_category=_scalar(cell.get("COG_category")),
                    description=_scalar(cell.get("Description")),
                    preferred_name=_scalar(cell.get("Preferred_name")),
                    gos=_split(cell.get("GOs")),
                    ecs=_split(cell.get("EC")),
                    kegg_kos=strip_ko_prefixes(_split(cell.get("KEGG_ko"))),
                    kegg_pathways=_split(cell.get("KEGG_Pathway")),
                    brites=_split(cell.get("BRITE")),
                    pfams=_split(cell.get("PFAMs")),
                )
            )
    if columns is None:
        raise AnnotationFormatError(
            "%s has no %s%s header line; it is not an eggNOG-mapper annotations file"
            % (path, COMMENT, ANNOTATION_COLUMNS[0])
        )
    return rows


def strip_ko_prefixes(kos: Iterable[str]) -> Tuple[str, ...]:
    """``("ko:K00059",)`` -> ``("K00059",)``, leaving anything already bare alone.

    Separate and public because this is the join that the whole feature is for: an unstripped
    ``ko:`` makes every comparison with ``kos.tsv`` come back empty, and an empty intersection is
    indistinguishable from total disagreement between the two annotators.
    """
    out: List[str] = []
    for ko in kos:
        bare = ko.removeprefix(KO_PREFIX)
        if bare and bare != EMPTY and bare not in out:
            out.append(bare)
    return tuple(out)


# ── this pipeline's own tables ────────────────────────────────────────────────────────────────

#: Columns of the per-gene table. The orthogonal annotations lead, because they are the reason this
#: step ran: a reader scanning for what eggNOG added over KOfam and Pfam should not have to pass the
#: hit statistics to reach them.
GENE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "cog_category",
    "go_terms",
    "ec_numbers",
    "kegg_ko",
    "description",
    "preferred_name",
    "orthologous_groups",
    "max_annot_lvl",
    "seed_ortholog",
    "evalue",
    "score",
    "kegg_pathway",
    "brite",
    "pfams",
    "eggnog_db_version",
)


def gene_row(row: EggnogRow, db_version: str = "unknown") -> Dict[str, str]:
    """One row of the per-gene table.

    ``;`` separated rather than ``,``, matching the CAZyme and Pfam tables: a comma inside a
    tab-separated field is legal and so is a comma-joined list, and having one separator across this
    module's outputs means a downstream split is the same split everywhere.
    """
    return {
        "genome": row.genome,
        "gene": row.gene,
        "cog_category": row.cog_category,
        "go_terms": ";".join(row.gos),
        "ec_numbers": ";".join(row.ecs),
        # Bare K numbers, so this column joins kos.tsv directly - see KO_PREFIX.
        "kegg_ko": ";".join(row.kegg_kos),
        "description": row.description,
        "preferred_name": row.preferred_name,
        "orthologous_groups": ";".join(row.og_names),
        "max_annot_lvl": row.max_annot_lvl,
        "seed_ortholog": row.seed_ortholog,
        "evalue": "" if row.evalue is None else "%.2g" % row.evalue,
        "score": "" if row.score is None else "%.1f" % row.score,
        "kegg_pathway": ";".join(row.kegg_pathways),
        "brite": ";".join(row.brites),
        "pfams": ";".join(row.pfams),
        # Stamped per row for the reason kofam.py stamps its threshold list: an orthology transfer
        # is only reproducible against the eggNOG release that held the group, and eggNOG renumbers
        # groups between releases.
        "eggnog_db_version": db_version,
    }


#: The namespaces the long table splits out, in the order its rows appear. Fixed rather than derived
#: from a dict so the table's row order is a property of this constant and not of iteration order.
TERM_NAMESPACES: Tuple[str, ...] = ("COG", "GO", "EC", "KEGG_ko")

#: Columns of the long per-term table: one row per (genome, gene, namespace, term). Long rather than
#: wide because a gene has any number of GO terms and exactly one row per term is what a group-by
#: can count - the wide column in GENE_COLUMNS is for reading, this is for arithmetic.
TERM_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "namespace",
    "term",
    "eggnog_db_version",
)


def _terms(row: EggnogRow, namespace: str) -> Tuple[str, ...]:
    if namespace == "COG":
        return row.cog_categories
    if namespace == "GO":
        return row.gos
    if namespace == "EC":
        return row.ecs
    return row.kegg_kos


def term_rows(rows: Sequence[EggnogRow], db_version: str = "unknown") -> List[Dict[str, str]]:
    """The long table: one row per term a gene carries, no row for a gene that carries none.

    Deterministic in input order and in TERM_NAMESPACES order, with no re-sorting, so two runs over
    the same annotations file write byte-identical output and a diff between two runs is a real
    difference rather than a reordering.
    """
    out: List[Dict[str, str]] = []
    for row in rows:
        for namespace in TERM_NAMESPACES:
            for term in _terms(row, namespace):
                out.append(
                    {
                        "genome": row.genome,
                        "gene": row.gene,
                        "namespace": namespace,
                        "term": term,
                        "eggnog_db_version": db_version,
                    }
                )
    return out


# ── what eggNOG added, and where it disagrees ─────────────────────────────────────────────────


@dataclass(frozen=True)
class Coverage:
    """How much of a gene set eggNOG annotated, and how much of that was new.

    ``genes_with_orthogonal_annotation`` is the number this feature is judged on. Pfam already
    covers nearly every gene KEGG covers - 1,578 Pfam-only genes against 13 KEGG-only ones in the
    bin measured in 2.1.2 - so "eggNOG annotated 78% of genes" says almost nothing. What it added is
    the GO, EC and COG that neither of the other two emits at all.
    """

    genes: int
    genes_annotated: int
    genes_with_orthogonal_annotation: int
    go_terms: int
    ec_numbers: int
    cog_categories: int

    @property
    def fraction_annotated(self) -> Optional[float]:
        return self.genes_annotated / self.genes if self.genes else None


def coverage(rows: Sequence[EggnogRow], genes: Sequence[str]) -> Coverage:
    """Count what eggNOG placed over a known gene set.

    *genes* is every gene called, not every gene in the annotations file: emapper writes no line for
    a gene it could not place, so the denominator has to come from the gene table.
    """
    gos: Set[str] = set()
    ecs: Set[str] = set()
    cogs: Set[str] = set()
    for row in rows:
        gos.update(row.gos)
        ecs.update(row.ecs)
        cogs.update(row.cog_categories)
    return Coverage(
        genes=len(set(genes)),
        genes_annotated=len({row.gene for row in rows}),
        genes_with_orthogonal_annotation=len(
            {row.gene for row in rows if row.has_orthogonal_annotation}
        ),
        go_terms=len(gos),
        ec_numbers=len(ecs),
        cog_categories=len(cogs),
    )


@dataclass(frozen=True)
class KoOverlap:
    """Where eggNOG's KEGG orthologs and KOfam's agree, per gene.

    Three counts and no verdict. eggNOG transfers a KO from an orthologous group; KOfam calls one
    from that profile's curated bit-score threshold. Those are different inferences with different
    failure modes, and merging them - taking either, or preferring one - throws away the only thing
    the pair measures. Both tables stay on disk; this says how far apart they are.
    """

    #: Genes both annotators gave at least one KO to, whether or not it is the same KO.
    both: int
    eggnog_only: int
    kofam_only: int
    #: Of `both`, the genes where the two KO sets actually intersect. The gap between this and
    #: `both` is the genes where they agree that there is a KO and disagree about which one, which
    #: is the subset worth looking at by hand.
    agreeing: int

    @property
    def total(self) -> int:
        return self.both + self.eggnog_only + self.kofam_only

    def describe(self) -> str:
        """One line for the console, phrased so neither annotator looks authoritative."""
        return (
            "KEGG orthologs: %d gene(s) called by both eggNOG and KOfam (%d on the same KO), "
            "%d by eggNOG only, %d by KOfam only; both tables are kept and neither is preferred"
            % (self.both, self.agreeing, self.eggnog_only, self.kofam_only)
        )


def compare_kos(rows: Sequence[EggnogRow], kofam_kos: Mapping[str, Iterable[str]]) -> KoOverlap:
    """Compare this genome's eggNOG KOs against KOfam's, per gene.

    *kofam_kos* is gene -> the K numbers ``--kofam`` assigned it. Both sides must be bare K numbers
    for the intersection to mean anything, which is what :func:`strip_ko_prefixes` is for: eggNOG
    writes ``ko:K00059`` and ``kos.tsv`` holds ``K00059``, so an unstripped prefix would report
    perfect disagreement on genes that agree exactly.
    """
    mine: Dict[str, Set[str]] = {}
    for row in rows:
        if row.kegg_kos:
            mine.setdefault(row.gene, set()).update(row.kegg_kos)
    theirs = {
        gene: set(strip_ko_prefixes(kos))
        for gene, kos in kofam_kos.items()
        if set(strip_ko_prefixes(kos))
    }
    shared = set(mine) & set(theirs)
    return KoOverlap(
        both=len(shared),
        eggnog_only=len(set(mine) - set(theirs)),
        kofam_only=len(set(theirs) - set(mine)),
        agreeing=sum(1 for gene in shared if mine[gene] & theirs[gene]),
    )


#: Printed with every eggNOG result.
TRANSFER_CAVEAT = (
    "an eggNOG annotation is transferred from an orthologous group, not matched to a profile: the "
    "GO, EC and COG terms on a gene are what characterised members of its group carry, and "
    "max_annot_lvl says how far up the tree that transfer came from - a phylum-level transfer and a "
    "genus-level one are not equally strong evidence"
)

#: Printed when eggNOG is reported beside KOfam and Pfam.
ORTHOGONAL_CAVEAT = (
    "the number to read is not how many genes eggNOG annotated but how many it gave a GO term, EC "
    "number or COG category to: Pfam already covers nearly every gene KEGG covers - 1578 genes had "
    "a Pfam domain and no KO against 13 the other way round on the bin this was measured on - so "
    "coverage is not what the third search buys, orthogonal annotation is"
)
