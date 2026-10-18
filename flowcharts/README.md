# Flowcharts

Diagrams of MetaWrap2's modules and workflows, drawn in draw.io.

| Diagram | What it shows |
|---|---|
| [overview](png/overview.png) | Every module and how they connect (the README figure) |
| [end_to_end](png/end_to_end.png) | A whole study: reads, co-assembly or per-sample recovery, then genomes, function and beyond |
| [reads_and_assembly](png/reads_and_assembly.png) | READ_QC, ASSEMBLY, KRAKEN2 and NONPAREIL |
| [binning](png/binning.png) | COVERAGE and BINNING, from reads to one bin set per binner |
| [bin_refinement](png/bin_refinement.png) | BIN_REFINEMENT, from up to six bin sets to one refined set |
| [reassemble_bins](png/reassemble_bins.png) | REASSEMBLE_BINS |
| [genome_catalogue](png/genome_catalogue.png) | DEREPLICATE, QUANT_BINS, CLASSIFY_BINS, APPRAISE and BLOBOLOGY |
| [function](png/function.png) | ANNOTATE, PATHWAYS, QUANTIFY_FUNCTION, DIFF_ABUNDANCE, CATALOGUE and AMR |
| [beyond_bins](png/beyond_bins.png) | VIRUSES, STRAINS, PANGENOME and PHYLOGENY |

## Layout

- `diagrams/` holds one Python file per diagram, which places the boxes and arrows.
- `drawio/` holds the generated `.drawio` files. Each opens and edits in draw.io.
- `png/` holds the PNGs, exported by draw.io itself.
- `drawio.py` is the library the diagrams use. It sizes boxes from their text, attaches arrows to
  the boxes they connect, checks that no two boxes overlap, and exports.

## Rebuilding

```bash
brew install --cask drawio        # once: the exporter is the draw.io desktop app
python flowcharts/build.py        # every diagram
python flowcharts/build.py binning
```

Building the overview also copies it to `docs/figures/metawrap2_overview.png`, which the README and
the docs site use. The `.drawio` files are generated, so to change a diagram, edit its Python file and rebuild.
