# GECCO output

GECCO 0.11.0 output, unedited, over the genes `metawrap2 annotate` called from two CAMI2 bins.

`clusters.tsv` is bin_003: four clusters. Two start within 400 bp of a contig end. Two are
classified `Unknown`, with ripp at p=0.46 and saccharide at p=0.34 as their best classes. The file
tests truncated clusters and unclassified clusters.

`clusters_empty.tsv` is bin_006, which has no clusters. With `--force-tsv`, GECCO writes a
placeholder row instead of a header alone: every field is blank except `type` (`Unknown`) and
`proteins`/`domains` (`""`). Parsed literally, this is one unclassified cluster of zero genes on no
contig, and it was counted as a BGC before the parser recognised it. The header is also shorter
because the per-class probability columns are absent when nothing is classified, so the parser must
read columns by name in both files.
