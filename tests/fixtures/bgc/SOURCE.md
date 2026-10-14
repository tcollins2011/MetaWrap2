# GECCO output, as it really looks

Real GECCO 0.11.0 output over the genes `metawrap2 annotate` called from two CAMI2 bins, unedited.

`clusters.tsv` is bin_003: four clusters, two of them starting within 400 bp of a contig end, two
classified `Unknown` with ripp at p=0.46 and saccharide at p=0.34 as their best classes. It is the
fixture for the cases that matter - a truncated cluster, and an unclassified cluster that is not an
unreal one.

`clusters_empty.tsv` is bin_006, which has no clusters, and is the more important fixture. With
`--force-tsv` GECCO writes a **placeholder row** rather than a header alone: every field blank
except `type`, which says `Unknown`, and `proteins`/`domains`, which say `""`. Read literally that
is one unclassified cluster of zero genes on no contig, and it was counted as a BGC until the
parser learned to recognise it. The header is also shorter here - the per-class probability columns
are absent when there is nothing to classify - so a positional parser would mis-read both files.
