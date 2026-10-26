"""Parsing eggNOG-mapper's annotations table, and the `annotate --eggnog` step.

Parser cases, each of which would otherwise give a plausible table instead of an error:

* `-` is the absent-value sentinel in every column and must not become a GO term, EC number or COG
  category (the same defect as dbCAN's `-|-` EC numbers);
* `## ` comment lines bracket the rows, and the footer `## N queries scanned` has no tab and would
  parse as a gene;
* `COG_category` is concatenated, not comma-joined: `EGP` is three categories;
* `KEGG_ko` is `ko:`-prefixed and `kos.tsv` is not, so the KOfam join needs the prefix stripped;
* a renamed column must raise instead of giving a zero-row table;
* the long table must be byte-identical across runs.

Module cases: the `--eggnog` block runs after KOfam, because the overlap summary reads `kos.tsv`
from disk.

The fixtures are modelled on eggNOG-mapper 2.1.13's writer (see tests/fixtures/eggnog/SOURCE.md).
``run`` is stubbed, so eggNOG-mapper is never launched.
"""

from __future__ import annotations

import csv
import io
import os
import shutil
from typing import List

import pytest

from metawrap2 import eggnog
from metawrap2.modules import annotate

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "eggnog")
ANNOTATIONS = os.path.join(FIXTURES, "bin_006.emapper.annotations")
WITH_MD5 = os.path.join(FIXTURES, "bin_md5.emapper.annotations")
RENAMED = os.path.join(FIXTURES, "bin_renamed.emapper.annotations")


def _rows():
    return eggnog.parse_annotations(ANNOTATIONS, "bin_006")


def _by_gene():
    return {row.gene: row for row in _rows()}


# ── parsing ──────────────────────────────────────────────────────────────────────────────────


def test_the_annotations_file_is_parsed():
    rows = _rows()
    assert len(rows) == 6
    assert all(row.genome == "bin_006" for row in rows)
    gene = _by_gene()["k141_100_length_5000_cov_8.0000_1"]
    assert gene.seed_ortholog == "272563.CD630_02740"
    assert gene.evalue == pytest.approx(3.19e-243)
    assert gene.score == pytest.approx(676.0)
    assert gene.max_annot_lvl == "1239|Firmicutes"
    assert gene.preferred_name == "aldB"
    assert gene.description == "Aldehyde dehydrogenase family protein"


def test_comment_lines_at_both_ends_are_skipped():
    """The four `## ` lines above the header and the three below the last row are skipped.

    The footer `## 6 queries scanned` has no tab, so it would otherwise parse as a gene.
    """
    genes = [row.gene for row in _rows()]
    assert not any(gene.startswith("#") for gene in genes)
    assert not any("queries scanned" in gene for gene in genes)
    # The header line also starts with "#" and must be read.
    assert genes[0] == "k141_100_length_5000_cov_8.0000_1"


# ── the `-` sentinel ─────────────────────────────────────────────────────────────────────────


def test_the_dash_sentinel_is_never_an_annotation():
    """A gene with `-` in every annotation column parses as having no annotations."""
    gene = _by_gene()["k141_100_length_5000_cov_8.0000_2"]
    assert gene.gos == ()
    assert gene.ecs == ()
    assert gene.kegg_kos == ()
    assert gene.pfams == ()
    assert gene.cog_category == ""
    assert gene.cog_categories == ()
    assert gene.description == ""
    assert gene.preferred_name == ""
    assert not gene.has_orthogonal_annotation
    # It contributes no row to the long table.
    assert not [row for row in eggnog.term_rows(_rows()) if row["gene"] == gene.gene]


def test_no_table_this_module_writes_contains_the_sentinel():
    """`-` does not appear as a whole field in any written table."""
    rows = _rows()
    cells = []
    for row in rows:
        cells.extend(eggnog.gene_row(row, "5.0.2").values())
    for term in eggnog.term_rows(rows, "5.0.2"):
        cells.extend(term.values())
    assert eggnog.EMPTY not in cells
    # Nor inside a joined list (as `-;peptidoglycan` did in dbCAN output).
    assert not any(eggnog.EMPTY in str(cell).split(";") for cell in cells)


# ── multi-valued fields ──────────────────────────────────────────────────────────────────────


def test_multi_valued_fields_are_split_and_deduplicated():
    gene = _by_gene()["k141_200_length_3000_cov_5.0000_1"]
    assert gene.ecs == ("1.1.1.1", "1.1.1.2")
    assert gene.kegg_kos == ("K00001", "K00002")
    assert gene.gos == ("GO:0008270", "GO:0016491")
    assert gene.pfams == ("ADH_N", "ADH_zinc_N")
    # Duplicates collapse; emapper's order is kept.
    assert eggnog._split("b,a,b,,a") == ("b", "a")
    assert eggnog._split("-") == ()
    assert eggnog._split("a,-,b") == ("a", "b")


def test_cog_categories_are_concatenated_not_comma_separated():
    """`EGP` is split into three categories.

    COG_category is the one multi-valued column emapper does not comma-join, so it is not in
    MULTI_VALUED.
    """
    gene = _by_gene()["k141_100_length_5000_cov_8.0000_1"]
    assert gene.cog_category == "EGP"
    assert gene.cog_categories == ("E", "G", "P")
    assert "COG_category" not in eggnog.MULTI_VALUED
    single = _by_gene()["k141_200_length_3000_cov_5.0000_2"]
    assert single.cog_categories == ("S",)
    # A COG category alone counts as orthogonal annotation (neither KOfam nor Pfam emits it).
    assert single.has_orthogonal_annotation


# ── the `ko:` prefix ─────────────────────────────────────────────────────────────────────────


def test_the_ko_prefix_is_stripped_so_the_join_against_kos_tsv_works():
    """eggNOG's `ko:K00059` is stripped to `K00059`, as `annotate --kofam` writes it.

    Unstripped, the two KO sets would never intersect.
    """
    gene = _by_gene()["k141_100_length_5000_cov_8.0000_1"]
    assert gene.kegg_kos == ("K00059",)
    assert not any(ko.startswith(eggnog.KO_PREFIX) for row in _rows() for ko in row.kegg_kos)
    assert eggnog.strip_ko_prefixes(["ko:K00059", "K00060"]) == ("K00059", "K00060")
    assert eggnog.strip_ko_prefixes(["-", ""]) == ()


def test_the_go_prefix_is_kept_because_it_is_part_of_the_accession():
    """The `GO:` prefix is kept, since it is part of the accession."""
    gene = _by_gene()["k141_100_length_5000_cov_8.0000_1"]
    assert gene.gos == ("GO:0003674", "GO:0008152", "GO:0016620")
    assert all(go.startswith(eggnog.GO_PREFIX) for go in gene.gos)
    # KEGG_Pathway's ko-prefix has no colon and is not touched either.
    assert gene.kegg_pathways == ("ko00010", "ko01100", "map00010")


# ── unexpected headers ───────────────────────────────────────────────────────────────────────


def test_a_renamed_column_is_an_error_that_names_what_it_found():
    with pytest.raises(eggnog.AnnotationFormatError) as caught:
        eggnog.parse_annotations(RENAMED, "bin_006")
    message = str(caught.value)
    assert "COG_category" in message and "GOs" in message
    assert "COG_cat" in message and "GO_terms" in message


def test_a_file_with_no_header_is_an_error_not_a_zero_row_table(tmp_path):
    """A headerless or comment-only file raises instead of giving a zero-row table.

    A zero-row table would look the same as a genome in which eggNOG placed nothing.
    """
    path = tmp_path / "headerless.emapper.annotations"
    path.write_text("## emapper-2.1.13\ngene_1\t272563.CD630_02740\t1e-9\t10.0\n")
    with pytest.raises(eggnog.AnnotationFormatError):
        eggnog.parse_annotations(str(path), "bin_006")

    empty = tmp_path / "comments_only.emapper.annotations"
    empty.write_text("## emapper-2.1.13\n## 0 queries scanned\n")
    with pytest.raises(eggnog.AnnotationFormatError):
        eggnog.parse_annotations(str(empty), "bin_006")


def test_the_md5_column_is_accepted():
    """The extra column that `--md5` appends is accepted."""
    rows = eggnog.parse_annotations(WITH_MD5, "bin_md5")
    assert len(rows) == 2
    assert rows[0].kegg_kos == ("K00059",)


# ── overlap with KOfam ───────────────────────────────────────────────────────────────────────


#: KOfam assignments as _kofam_kos() reads them from kos.tsv. Covers every case: one gene where the
#: two agree on the KO, one where they call different KOs, two eggNOG-only and two KOfam-only.
KOFAM = {
    "k141_100_length_5000_cov_8.0000_1": ["K00059"],
    "k141_200_length_3000_cov_5.0000_1": ["K09999"],
    "k141_400_length_2000_cov_4.0000_1": ["K00928"],
    "k141_400_length_2000_cov_4.0000_2": ["K01626"],
    "k141_500_length_1000_cov_2.0000_1": [],
}


def test_ko_overlap_counts_both_sides_and_prefers_neither():
    overlap = eggnog.compare_kos(_rows(), KOFAM)
    # eggNOG calls a KO on 4 genes, KOfam on 4 (the fifth entry is empty and is not a call).
    assert (overlap.both, overlap.eggnog_only, overlap.kofam_only) == (2, 2, 2)
    # Of the 2 both called, 1 is the same KO and 1 is a disagreement.
    assert overlap.agreeing == 1
    assert overlap.total == 6
    sentence = overlap.describe()
    assert "2 gene(s) called by both eggNOG and KOfam (1 on the same KO)" in sentence
    assert "2 by eggNOG only, 2 by KOfam only" in sentence
    # The wording does not prefer either tool.
    assert "neither is preferred" in sentence


def test_a_gene_with_no_ko_on_either_side_is_in_no_cell_of_the_overlap():
    """An empty KO list is not counted as a KOfam call."""
    overlap = eggnog.compare_kos(_rows(), KOFAM)
    assert overlap.total == 6
    assert "k141_500_length_1000_cov_2.0000_1" not in KOFAM or overlap.kofam_only == 2


def test_the_overlap_survives_a_prefix_on_either_side():
    """The `ko:` prefix is stripped on both sides of the comparison."""
    bare = eggnog.compare_kos(_rows(), {"k141_100_length_5000_cov_8.0000_1": ["K00059"]})
    prefixed = eggnog.compare_kos(_rows(), {"k141_100_length_5000_cov_8.0000_1": ["ko:K00059"]})
    assert bare == prefixed
    assert prefixed.agreeing == 1


# ── coverage ─────────────────────────────────────────────────────────────────────────────────


def test_coverage_separates_annotated_from_orthogonally_annotated():
    """Coverage counts placed genes and, separately, genes with GO, EC or COG terms.

    The second count is the relevant one: Pfam already covers nearly every gene KEGG covers (1,578
    Pfam-only against 13 KEGG-only genes on the bin measured in 2.1.2). Two of the six placed genes
    here have no GO, EC or COG term.
    """
    genes = [row.gene for row in _rows()] + ["unplaced_1", "unplaced_2"]
    found = eggnog.coverage(_rows(), genes)
    assert found.genes == 8
    assert found.genes_annotated == 6
    assert found.genes_with_orthogonal_annotation == 4
    assert found.fraction_annotated == pytest.approx(0.75)
    assert (found.go_terms, found.ec_numbers, found.cog_categories) == (7, 4, 5)


# ── the tables ───────────────────────────────────────────────────────────────────────────────


def _write(columns, rows):
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def test_the_long_table_is_byte_identical_across_runs():
    """Two parses of one file write the same bytes.

    Neither the split fields nor the long table are re-sorted; emapper's order is kept.
    """
    first = _write(eggnog.TERM_COLUMNS, eggnog.term_rows(_rows(), "5.0.2"))
    second = _write(eggnog.TERM_COLUMNS, eggnog.term_rows(_rows(), "5.0.2"))
    assert first == second
    genes = _write(eggnog.GENE_COLUMNS, [eggnog.gene_row(row, "5.0.2") for row in _rows()])
    assert genes == _write(eggnog.GENE_COLUMNS, [eggnog.gene_row(row, "5.0.2") for row in _rows()])


def test_the_long_table_has_one_row_per_term():
    rows = eggnog.term_rows(_rows(), "5.0.2")
    # 8 terms on the first gene, 0 on the unannotated one, 7, 1, 1 and 5 on the rest.
    assert len(rows) == 22
    counts = {namespace: 0 for namespace in eggnog.TERM_NAMESPACES}
    for row in rows:
        counts[row["namespace"]] += 1
    assert counts == {"COG": 6, "GO": 7, "EC": 4, "KEGG_ko": 5}
    assert all(row["eggnog_db_version"] == "5.0.2" for row in rows)


def test_the_gene_table_carries_the_bare_ko_so_it_joins_kos_tsv():
    row = eggnog.gene_row(_by_gene()["k141_100_length_5000_cov_8.0000_1"], "5.0.2")
    assert row["kegg_ko"] == "K00059"
    assert row["cog_category"] == "EGP"
    assert row["go_terms"] == "GO:0003674;GO:0008152;GO:0016620"
    assert row["ec_numbers"] == "1.2.1.3"
    assert row["eggnog_db_version"] == "5.0.2"
    assert set(row) == set(eggnog.GENE_COLUMNS)


def test_the_column_set_is_pinned_to_the_tools_own_writer():
    """ANNOTATION_COLUMNS matches emapper's 21 columns, in order."""
    assert len(eggnog.ANNOTATION_COLUMNS) == 21
    assert eggnog.ANNOTATION_COLUMNS[0] == "query"
    assert eggnog.ANNOTATION_COLUMNS[-1] == "PFAMs"
    with open(ANNOTATIONS) as handle:
        header = next(line for line in handle if line.startswith("#query"))
    assert tuple(header.lstrip("#").rstrip("\n").split("\t")) == eggnog.ANNOTATION_COLUMNS


# ── the annotate --eggnog step ───────────────────────────────────────────────────────────────


GENOME = "bin_006"

#: Two genes emapper did not place. emapper writes no line for them, so the coverage share must use
#: the protein file as its denominator.
UNPLACED = ("k141_900_length_800_cov_1.0000_1", "k141_900_length_800_cov_1.0000_2")


def _said(capsys):
    """Everything printed, with banner borders and line wrapping removed.

    ``comm`` wraps messages in a fixed-width banner, so substring checks on the raw capture would
    depend on where lines break.
    """
    return " ".join(capsys.readouterr().out.replace("-", " ").replace("*", " ").split())


def _setup(tmp_path, kofam_table=None):
    """An output dir with one genome's proteins in it, and optionally a kos.tsv beside them."""
    out = tmp_path / "out"
    genes_dir = out / annotate.GENES_DIR
    genes_dir.mkdir(parents=True)
    # The six genes the fixture annotates, plus two it does not place.
    genes = [row.gene for row in _rows()] + list(UNPLACED)
    with open(genes_dir / (GENOME + ".faa"), "w") as handle:
        handle.writelines(">%s\nMKRISTTITTTITITTGNGAG\n" % gene for gene in genes)
    (tmp_path / "metawrap2.toml").write_text("[settings]\nuse_conda_envs = false\n")
    db = tmp_path / "eggnog_data"
    db.mkdir()
    for name in annotate.EGGNOG_DB_FILES:
        (db / name).touch()
    if kofam_table is not None:
        with open(out / annotate.KOS_TABLE, "w") as handle:
            handle.write("genome\tgene_id\tko\n")
            handle.writelines(kofam_table)
    return out, genes_dir, db


def _namespace(tmp_path, out, db):
    import argparse

    return argparse.Namespace(
        output=str(out),
        config=str(tmp_path / "metawrap2.toml"),
        threads=2,
        eggnog_db=str(db),
        bins="bins",
        assembly=None,
    )


@pytest.fixture
def fake_emapper(monkeypatch) -> List[str]:
    """Record the command, and drop the fixture annotations where emapper would have written it."""
    calls: List[str] = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        outdir = cmd.split("--output_dir ", 1)[1].split(" ", 1)[0]
        shutil.copy(ANNOTATIONS, os.path.join(outdir, "%s.emapper.annotations" % GENOME))
        return 0

    monkeypatch.setattr(annotate, "run", fake_run)
    return calls


def test_the_emapper_command_is_the_one_the_tool_documents(tmp_path, fake_emapper):
    """The emapper command has the flags required by emapper 2.1.13's argument parser."""
    out, genes_dir, db = _setup(tmp_path)
    annotate._run_eggnog(_namespace(tmp_path, out, db), str(genes_dir), [(GENOME, "x.fa", False)])
    assert len(fake_emapper) == 1
    cmd = fake_emapper[0]
    # Our called proteins; with --itype genome emapper would run Prodigal and call its own genes.
    assert "--itype proteins" in cmd
    assert "-m diamond" in cmd
    # Without --override emapper refuses to start when its outputs exist, as on every re-run.
    assert "--override" in cmd
    assert "--data_dir %s" % db in cmd
    assert "--cpu 2" in cmd
    # These are type=existing_dir in emapper's parser, so they must exist beforehand.
    for flag in ("--output_dir", "--temp_dir"):
        directory = cmd.split(flag + " ", 1)[1].split(" ", 1)[0]
        assert os.path.isdir(directory), "%s must exist before emapper is launched" % flag


def test_both_tables_are_written_with_the_declared_columns(tmp_path, fake_emapper):
    out, genes_dir, db = _setup(tmp_path)
    annotate._run_eggnog(_namespace(tmp_path, out, db), str(genes_dir), [(GENOME, "x.fa", False)])
    with open(out / annotate.EGGNOG_TABLE) as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(genes[0]) == eggnog.GENE_COLUMNS
    assert len(genes) == 6
    assert {row["genome"] for row in genes} == {GENOME}
    with open(out / annotate.EGGNOG_TERM_TABLE) as handle:
        terms = list(csv.DictReader(handle, delimiter="\t"))
    assert tuple(terms[0]) == eggnog.TERM_COLUMNS
    assert len(terms) == 22
    # The sentinel is not written to either table.
    assert eggnog.EMPTY not in [value for row in genes + terms for value in row.values()]


def test_the_overlap_is_reported_only_when_kofam_ran_first(tmp_path, fake_emapper, capsys):
    """The KOfam overlap is reported only when kos.tsv exists.

    Without kos.tsv no overlap line is printed. (The Pfam block once ran before KOfam and its
    increment was always zero.)
    """
    out, genes_dir, db = _setup(tmp_path)
    args = _namespace(tmp_path, out, db)
    annotate._run_eggnog(args, str(genes_dir), [(GENOME, "x.fa", False)], annotate._kofam_kos(args))
    without = _said(capsys)
    assert "KEGG orthologs:" not in without, "no KO table means no overlap claim, not a zero one"
    # The coverage line is still reported, against the protein file's denominator.
    assert "6 of 8 gene(s) (75.0%) placed in an orthologous group" in without
    assert "4 of them carry a GO term, EC number or COG category" in without

    out2, genes_dir2, db2 = _setup(
        tmp_path / "second",
        kofam_table=[
            "%s\tk141_100_length_5000_cov_8.0000_1\tK00059\n" % GENOME,
            "%s\tk141_200_length_3000_cov_5.0000_1\tK09999\n" % GENOME,
            "%s\tk141_400_length_2000_cov_4.0000_1\tK00928\n" % GENOME,
            "%s\tk141_400_length_2000_cov_4.0000_2\tK01626\n" % GENOME,
        ],
    )
    args2 = _namespace(tmp_path / "second", out2, db2)
    found = annotate._kofam_kos(args2)
    assert found[GENOME]["k141_100_length_5000_cov_8.0000_1"] == ["K00059"]
    annotate._run_eggnog(args2, str(genes_dir2), [(GENOME, "x.fa", False)], found)
    with_kofam = _said(capsys)
    assert "KEGG orthologs: 2 gene(s) called by both eggNOG and KOfam (1 on the same KO)" in (
        with_kofam
    )
    assert "2 by eggNOG only, 2 by KOfam only" in with_kofam
    assert "both tables are kept and neither is preferred" in with_kofam


def test_a_genome_emapper_fails_on_is_warned_about_not_fatal(tmp_path, monkeypatch, capsys):
    """A genome emapper fails on is warned about and the others are still written."""
    out, genes_dir, db = _setup(tmp_path)
    (genes_dir / "bin_007.faa").write_text(">bin_007_1\nMKRISTT\n")
    attempts: List[str] = []

    def fail_the_second(cmd, **kwargs):
        attempts.append(cmd)
        if len(attempts) == 1:
            outdir = cmd.split("--output_dir ", 1)[1].split(" ", 1)[0]
            shutil.copy(ANNOTATIONS, os.path.join(outdir, "%s.emapper.annotations" % GENOME))
            return 0
        return 1

    monkeypatch.setattr(annotate, "run", fail_the_second)
    annotate._run_eggnog(
        _namespace(tmp_path, out, db),
        str(genes_dir),
        [(GENOME, "x.fa", False), ("bin_007", "y.fa", False)],
    )
    printed = _said(capsys)
    assert "eggNOG mapper failed on bin_007" in printed
    with open(out / annotate.EGGNOG_TABLE) as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 6


def test_every_genome_failing_fails_the_step(tmp_path, monkeypatch, capsys):
    """When every genome fails, the step exits with status 1."""
    out, genes_dir, db = _setup(tmp_path)
    monkeypatch.setattr(annotate, "run", lambda cmd, **kwargs: 1)
    with pytest.raises(SystemExit) as stopped:
        annotate._run_eggnog(
            _namespace(tmp_path, out, db), str(genes_dir), [(GENOME, "x.fa", False)]
        )
    assert stopped.value.code == 1
    assert "no eggNOG mapper output was produced" in _said(capsys)
