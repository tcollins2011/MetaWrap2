# MetaWrap2 dependencies

MetaWrap2's core is lightweight: it only needs `python>=3.10` and `biopython` (installed by
`pip install -e .`). The heavy third-party tools each module drives (SPAdes, CheckM, salmon, KRAKEN2, ...)
live in their own per-module conda environment, defined by the YAML files in [`src/metawrap2/envs/`](../src/metawrap2/envs/). Splitting
them up means each tool's dependencies resolve independently, so adding or updating one tool no longer has
to satisfy every other tool at the same time.

You do not install these tools by hand. Create the environment(s) for the module(s) you plan to run with:
```bash
metawrap2 install-env <module>     # e.g. metawrap2 install-env binning
metawrap2 install-env --all        # create every environment
```
Each environment is named `metawrap2-<module>` and holds only that module's tools. You only need the
environments for the modules you use. Verify a module's tools (and databases) are present before running
it with `metawrap2 check <module>`.

## Per-module environments and their tools

| Module | Environment | Tools |
|---|---|---|
| read_qc | `metawrap2-read_qc` | trim-galore, fastqc, bmtagger, bowtie2 |
| assembly | `metawrap2-assembly` | spades (metaSPAdes), megahit, quast, bwa |
| kraken2 | `metawrap2-kraken2` | kraken2, krona (KronaTools) |
| binning | `metawrap2-binning` | metabat2, maxbin2, concoct, bwa, bowtie2, samtools |
| bin_refinement | `metawrap2-bin_refinement` | checkm-genome (CheckM1), pplacer |
| reassemble_bins | `metawrap2-reassemble_bins` | spades, bwa, samtools, minimap2 |
| quant_bins | `metawrap2-quant_bins` | salmon, bwa, samtools |
| blobology | `metawrap2-blobology` | blast (MEGABLAST), bowtie2, samtools, biopython, pandas, matplotlib |
| classify_bins | `metawrap2-classify_bins` | blast (MEGABLAST), taxator-tk |
| annotate_bins | `metawrap2-annotate_bins` | prokka |

## Optional (opt-in) environments
Some modules ship newer engines in separate environments that are off by default so existing runs
reproduce the same results:

| Environment | Tool | Module |
|---|---|---|
| `metawrap2-bin_refinement-checkm2` | CheckM2 (avoids the CheckM1/pplacer memory cost) | bin_refinement |
| `metawrap2-classify_bins-gtdbtk` | GTDB-Tk | classify_bins |
| `metawrap2-annotate_bins-bakta` | Bakta | annotate_bins |

See [`src/metawrap2/envs/README.md`](../src/metawrap2/envs/README.md) for the authoritative, up-to-date list of environments and the
exact packages each one pins.
