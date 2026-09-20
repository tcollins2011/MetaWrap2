from metawrap2.modules import classify_bins


def test_command_templates_format_cleanly():
    # blastdb is now the full -db argument (directory + database name), so a database that
    # is not called "nt" works too.
    b = classify_bins.BLASTN.format(threads=8, blastdb="/db/nt", query="q.fa")
    assert "{" not in b and "-num_threads 8" in b and "-db /db/nt" in b and "-query q.fa" in b
    b = classify_bins.BLASTN.format(threads=1, blastdb="/db/core_nt", query="q.fa")
    assert "-db /db/core_nt" in b
    t = classify_bins.TAXATOR.format(
        taxdump="/tax",
        taxator_t=0.3,
        taxator_e=0.01,
        mapping="m.tax",
        input="in.tab",
        output="out.gff3",
    )
    assert "{" not in t and "TAXATORTK_TAXONOMY_NCBI=/tax" in t and "< in.tab > out.gff3" in t
    # binner and taxknife are taxator-tk tools too, so they also need the taxonomy env var -
    # without it binner aborts with "Specify the folder containing the NCBI taxonomy dump
    # files as TAXATORTK_TAXONOMY_NCBI environment variable".
    bn = classify_bins.BINNER.format(
        taxdump="/tax", predictions="p.gff3", genus_cutoff=0.6, output="b.txt"
    )
    assert "{" not in bn and "genus:0.6" in bn
    assert "TAXATORTK_TAXONOMY_NCBI=/tax" in bn
    tk = classify_bins.TAXKNIFE.format(taxdump="/tax", binned="b.txt", output="c.tab")
    assert "{" not in tk and "grep -v 'Could not'" in tk
    assert "TAXATORTK_TAXONOMY_NCBI=/tax" in tk


def _raw_line(qseqid, staxids):
    cols = [qseqid, "1", "100", "500", "s1", staxids, "10", "110", "200", "0.001", "95", "100"]
    return "\t".join(cols)
