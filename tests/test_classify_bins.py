from metawrap2.modules import classify_bins
from metawrap2.scripts import classify_bins as classify_helper
from metawrap2.scripts import prune_blast_hits


def test_command_templates_format_cleanly():
    b = classify_bins.BLASTN.format(threads=8, blastdb="/db", query="q.fa")
    assert "{" not in b and "-num_threads 8" in b and "-db /db/nt" in b and "-query q.fa" in b
    t = classify_bins.TAXATOR.format(taxdump="/tax", taxator_t=0.3, taxator_e=0.01,
                                     mapping="m.tax", input="in.tab", output="out.gff3")
    assert "{" not in t and "TAXATORTK_TAXONOMY_NCBI=/tax" in t and "< in.tab > out.gff3" in t
    bn = classify_bins.BINNER.format(predictions="p.gff3", genus_cutoff=0.6, output="b.txt")
    assert "{" not in bn and "genus:0.6" in bn
    tk = classify_bins.TAXKNIFE.format(binned="b.txt", output="c.tab")
    assert "{" not in tk and "grep -v 'Could not'" in tk


def _raw_line(qseqid, staxids):
    cols = [qseqid, "1", "100", "500", "s1", staxids, "10", "110", "200", "0.001", "95", "100"]
    return "\t".join(cols)


def test_prune_keeps_first_ranked_id(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    nodes.write_text("562\t|\t561\t|\tspecies\t|\txx\n999\t|\t1\t|\tno rank\t|\txx\n")
    raw = tmp_path / "raw.tab"
    raw.write_text(_raw_line("c1", "999;562") + "\n" + _raw_line("c2", "999") + "\n")

    lines = list(prune_blast_hits.prune(str(nodes), str(raw)))
    # c1 keeps first ranked id (562, since 999 is unranked); c2 has no ranked id -> dropped
    assert len(lines) == 1
    cols = lines[0].split("\t")
    assert len(cols) == 12 and cols[0] == "c1" and cols[5] == "562"


def test_cut_and_mapping_columns(tmp_path):
    pruned = tmp_path / "pruned.tab"
    pruned.write_text(_raw_line("c1", "562") + "\n")
    tab = tmp_path / "out.tab"
    mapping = tmp_path / "map.tax"
    classify_bins._write_pruned_columns(str(pruned), str(tab))
    classify_bins._write_mapping(str(pruned), str(mapping))
    assert tab.read_text().strip().split("\t") == \
        ["c1", "1", "100", "500", "s1", "10", "110", "200", "0.001", "95", "100"]
    assert mapping.read_text().strip() == "s1\t562"


def test_consensus_taxonomy(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">c1\nAAAA\n>c2\nCCCC\n")
    tax = tmp_path / "contig_taxonomy.tab"
    tax.write_text("c1\tBacteria;Firmicutes\nc2\tBacteria;Firmicutes\n")
    res = list(classify_helper.consensus_taxonomy(str(tax), str(bins)))
    assert res == [("bin.1.fa", "Bacteria;Firmicutes")]
