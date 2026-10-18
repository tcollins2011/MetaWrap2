# Plan: functional annotation, pathways, AMR, and abundance

Plan to replace `annotate_bins`, a Prokka/Bakta wrapper that writes GFF, `.faa` and `.ffn`, with a
subsystem covering:

1. What genes does this genome carry? (structural and functional annotation)
2. What can it do? (orthology, then pathway/module completeness)
3. What resistance and virulence does it carry, and in what genetic context?
4. How much of each function is in the sample? (gene and pathway abundance, per genome and
   community-wide, with and without bins)

Each section starts with its status; §11 has the state of every phase. "No status recorded" means this
plan does not say whether the item is done. Numbers quoted for other tools describe those tools and are
not targets for MetaWrap2.

---

## 1. Current state and existing tools

Context.

`annotate_bins` runs one tool per bin and reorganises its output. It has no identifiers joinable
across genomes, no pathways, AMR or abundance, and requires bins.

Other pipelines: SqueezeMeta writes KEGG/COG/Pfam abundance tables for reads and bins, with R
bindings. DRAM condenses KEGG/Pfam/dbCAN/MEROPS into a metabolism summary. VEBA scores BGC novelty.
MGnify has a mobilome pipeline. nf-core/funcscan screens contigs for AMR/AMP/BGC/CAZy. MetaWrap2
currently provides gene calling.

Permissively licensed implementations to read first:

- mettannotator (Apache-2.0, EMBL-EBI, validated on MAGs), the closest analogue.
- MGnify's `assembly-analysis-pipeline`: the bin-free path, with exact tool and database versions.
- nf-core/funcscan 4.0.0: database handling and version pinning for the AMR and BGC arms.

---

## 2. Licensing

Context and two rules (no status recorded for the rules).

For non-profit academic use all candidate tools are available: CARD and RGI, geNomad, VFDB, KEGG's
KOfam profiles and module definitions, MetaCyc through HUMAnN.

- Do not call `rest.kegg.jp` from a pipeline run. It is rate-limited to 3 requests/s. Module
  definitions come from a vendored, dated snapshot.
- Pin and record every database version; a KOfam or CARD snapshot changes results as much as a
  tool version. `archives/` on the KEGG FTP has dated KOfam snapshots.

Commercial redistribution would need a review of CARD/RGI (prohibited for commercial organisations
without McMaster's permission), geNomad (LBNL academic-use licence), VFDB and integron_finder's data
files (CC BY-NC), and MetaCyc (paid).

Planned tool chain: Pyrodigal → Bakta → PyHMMER on KOfam/Pfam/NCBIfam → eggNOG-mapper →
kegg-pathways-completeness → AMRFinderPlus + RGI + hAMRonization → dbCAN → CoverM/SingleM.

---

## 3. Architecture

Status: done; all five modules exist (§11, phases 1–7). No status recorded for the `doctor` database
probes.

`annotate_bins` is replaced by five modules. None of them requires bins.

```
                  ┌─ annotate (genes + function) ──┬── per-genome tables
assembly ─┬─ bins ─┤                                │
          │        ├─ pathways (modules/completeness)┤
          │        ├─ amr (ARG/VF + context)        ├── join on protein id
          │        └─ bgc/cazy (opt-in)             │
          └─ contigs (bin-free) ── gene catalogue ──┴── abundance ── sample × feature matrices
```

| Module | Input | Output |
|---|---|---|
| `annotate` | bins or a whole assembly | genes, products, orthologs, domains |
| `pathways` | `annotate` output | module/pathway completeness per genome |
| `amr` | bins or contigs | ARGs, VFs, and their genetic context |
| `quantify_function` | `annotate` + reads | gene and pathway abundance, per genome and per sample |
| `catalogue` | contigs from any number of samples | non-redundant gene catalogue for the bin-free path |

`annotate_bins` stays as a deprecated alias for `annotate`, with a warning, because published
scripts call it.

Gene prediction and homology searches run once over all contigs, and downstream modules score
groupings of the cached features (as Binette does for refinement). `annotate_bins` re-runs Prokka per
bin and re-predicts genes on shared contigs. With the shared front end, the second of the binned and
bin-free paths adds little compute.

### 3.1 CLI surface

Each module follows existing conventions: a `COMMANDS` block at the top of the file with the exact
command each tool runs, a `CONSTANTS` block for tunable numbers, `-o/--output`, `-t/--threads`,
`--dry-run`, `--resume`, and one conda env per workflow.

As planned (the shipped flags differ; see `mw2 <module> --help`):

```
mw2 annotate            -o OUT (-b BINS | -a ASSEMBLY) [--bakta-db {light,full}]
                        [--eggnog] [--kofam] [--interproscan]
mw2 pathways            -o OUT --features features.tsv [--definitions SNAPSHOT]
                        [--min-completeness 0.0]
mw2 amr                 -o OUT (-b BINS | -a ASSEMBLY) [--resfinder] [--rgi]
                        [--mobility-window 10000] [--vf]  [--plasmids]
mw2 quantify_function   -o OUT --features features.tsv --reads ... [--genomes BINS]
                        [--units cpg,per_gb,relab]
mw2 catalogue           -o OUT -a ASSEMBLY [ASSEMBLY ...] [--identity 0.95] [--coverage 0.90]
```

Two new envs: `metawrap2-annotate` (Pyrodigal, PyHMMER, Bakta, eggNOG-mapper) and `metawrap2-amr`
(AMRFinderPlus, hAMRonization, ISEScan, Platon). `pathways` and `quantify_function` are pure Python plus
CoverM/SingleM and use existing envs.

`doctor` gets a probe per new database (KOfam snapshot, Pfam, eggNOG 5, AMRFinderPlus DB,
module-definition snapshot): installed / missing / broken, with version and date. The express
installer does not download them; `doctor --fix` prompts.

---

## 4. `annotate`

### 4.1 Gene calling

Status: Pyrodigal with per-gene `partial` and `distance_to_contig_end` is done (`metawrap2.genecall`).
No status recorded for a second caller.

Pyrodigal (MIT, maintained Cython bindings to Prodigal, which is no longer developed) runs
in-process and releases the GIL.

- MAG contigs are short, so genes running off a contig end are called partial or missed. `partial` and
  `distance_to_contig_end` are recorded per gene and carried into every downstream result (module
  steps, ARGs, BGC borders).
- A second caller would help with fragmented genes; MGnify's assembly pipeline merges two.
  FragGeneScanRs is the candidate. The merge rule must not count a gene twice.

### 4.2 Functional annotation

Status: open. The Bakta pin is not done (`bakta >=1.9`). No status recorded for the light-database
default.

Bakta (GPL-3, actively developed) stays the per-genome annotator. Its tRNAscan-SE, Aragorn and
Infernal/Rfam output gives the rRNA and tRNA counts MIMAG requires (bin QC, companion plan).

Bakta's full database is ~70 GB and `--db-light` is ~2 GB. The difference is mostly UniRef cluster
coverage, which improves product names but not the identifiers that modules join on. Default to light,
make full a flag, and record which was used, since product names differ between the two.

`envs/annotate.yaml` asks for `bakta >=1.9`, which resolves to 1.12.x, so a later rerun can give
different products with no record of the version. Pin it and record the version, as `provenance` does
for databases. The same applies to `gtdbtk >=2.4.0` (companion plan §1.5).

### 4.3 Orthology and domains

Status: KO assignment done (`metawrap2.kofam`; needs the 1.5 GB profile set to run on real data). Pfam
done (`annotate --pfam`). No status recorded for eggNOG-mapper. InterProScan is inconsistent: "opt-in"
and "decided against" in this section, "deferred rather than rejected" in §11.

Most annotation error comes from this step. On 396 genomes with identical gene calls (Kananen et al.,
*Bioinform Adv* 2025, `10.1093/bioadv/vbaf039`):

| | calls confirmed by eggNOG | tool-specific calls confirmed | tool-specific calls contradicted |
|---|---|---|---|
| KofamScan | 84% | - | - |
| anvi'o (adaptive thresholds) | 79% | 37% | 7.5% |
| MicrobeAnnotator | 74% | 0.3% | 35% |

anvi'o's adaptive bit-score heuristic annotates 13% more genes per genome (p<2.2e-16), mostly in
under-represented taxa (38% of Nanosynbacteraceae cases against 4.9% of Enterobacteriaceae). Without
the heuristic, anvi'o matches KofamScan.

KOs are assigned with PyHMMER over KOfam profiles, with an anvi'o-style adaptive threshold where the
curated cutoff is missing or fails. MicrobeAnnotator's permissive path is excluded (35% of its unique
calls contradicted). 1,866 of 28,453 KOfam profiles (6.6%) have no bit-score threshold and are a
known source of false positives; drop them or give them a conservative threshold, and record which.

PyHMMER (MIT) replaces `hmmsearch` as a subprocess and keeps the threshold policy in code. There is
no published hmmsearch-vs-Pfam timing at ~1M proteins, so this step's compute must be measured before
the module ships.

eggNOG-mapper is the second opinion and the source of COG/GO/EC/CAZy. Pin 2.1.15 with eggNOG 5.
v3.0.0-beta6 is a beta and uses eggNOG 7, which is incompatible with eggNOG 5 data. Its
orthology-transferred KOs ignore KOfam's per-KO cutoffs and reflect the KEGG version eggNOG was built
on, so eggNOG is used to check KO calls and not to assign them.

InterProScan 6: opt-in only, because of cost (~14 GB on disk, too slow for every MAG in a cohort).

Decided against InterProScan; Pfam built instead. `annotate --pfam` searches Pfam-A with PyHMMER at each
family's gathering threshold, similar to KOfam's per-profile cutoffs. It is 399 MB and runs in-process.
On one CAMI2 bin, KEGG placed 1,620 of 3,891 genes (41.6%) and Pfam placed 3,185 (81.9%); 1,578 of the
Pfam-placed genes have no KEGG ortholog.

Not covered: InterProScan's other members (TIGRFAM, SMART, PANTHER, SUPERFAMILY), signal peptides,
GO terms and transmembrane topology. The bioconda InterProScan build is four years old (5.54_87.0,
2021).

---

## 5. `pathways`

### 5.1 Module completeness in existing tools

Context.

| Tool | Semantics | Bias |
|---|---|---|
| anvi'o | unrolls every alternative path; complexes score as the fraction of *essential* subunits; `--` steps always 0; module = max over paths | optimistic; the only one that also offers a stepwise score |
| EBI kegg-pathways-completeness | definition → weighted DiGraph, traversal minimising weight ratio; `--` substituted at weight 0 | a `--` step does not cap the module, unlike anvi'o |
| MicrobeAnnotator | three separate algorithms; complexes get proportional partial credit | most permissive on complexes |
| DRAM | a node is kept if any of its KOs is present; only 13 hardcoded modules | most permissive overall |
| METABOLIC | substitutes 1/0 and `eval`s the boolean; present if ≥0.75 of steps | module table dated 2023-06-12 |

On the same KO set, anvi'o reports 11–12% more modules at ≥80% completeness than MicrobeAnnotator or
KofamScan, so MetaWrap2 counts are only comparable with published counts made with the same algorithm.

Issues in the 573-module EBI snapshot:

- 47 modules (8.2%) have a single top-level step and 164 (28.6%) have two or fewer, so their
  completeness is close to binary and prone to false positives. M00851 "Carbapenem resistance" is a flat
  OR of 21 KOs: any one hit makes it 100% complete, and one of the 21 has no threshold. In the benchmark
  above, MicrobeAnnotator called it complete in all but two bacterial families.
- 14 modules contain a `--` step (a reaction with no KO) and cannot reach 100% under anvi'o.
- 45 KOs appear in five or more modules (K00024 in ten), which inflates unrelated modules.
- 8 modules reference other modules and need a second resolution pass.
- Definitions change over time: 573 in the pinned EBI snapshot against 582 live.

### 5.2 Design

Status: done (`metawrap2.kegg`, `modules/pathways.py`).

Completeness is computed in MetaWrap2 from vendored definitions by enumerating every path: 3 of 573
modules have more than 1,000 alternative routes, the maximum is 33,600 for M00001, and all 573 score
in 0.1 s.

A `+`-joined protein complex is one step whose subunits are all essential. M00142, NADH dehydrogenase,
has seven subunits; read as seven steps, a genome with three would score 3/7 instead of 0. 31 of the
573 modules are a single complex. Counting top-level spaces in the definition string gives the
published "47 single-step modules" through two errors that cancel (a complex read as a sequence, a
parenthesised group read as one step). The parser handles both and also gets 47.

Both scores are always written, `completeness_pathwise` (max over alternative paths) and
`completeness_stepwise` (fraction of top-level steps satisfied), plus the path chosen. The two can
differ by as much as the biological differences being compared.

| Flag | Meaning |
|---|---|
| `n_top_level_steps` | lets a reader discount the 8.2% of modules that are effectively binary |
| `is_signature_module` | AMR and signature modules do not indicate metabolic capability and are not reported as such |
| `has_no_ko_step` | the module cannot reach 100% under pathwise semantics |
| `proportion_unique_enzymes` | anvi'o's mitigation for promiscuous KOs; no other tool has one |
| `nt_ko_calls` | how many supporting KOs had no curated threshold |
| `definition_snapshot` | the vendored module-table date |

Vendor the EBI `modules_table.tsv` and `graphs.pkl` (Apache-2.0, dated), pin the version, write the
snapshot date on every row, and never mix snapshots within a run.

The output is labelled metabolic potential, because completeness does not predict phenotype. Across 56
strains in gnotobiotic honey bees, completeness "did not fully explain colonization patterns", and a
complete methylmalonyl-CoA pathway was never found in Cloacimonadota because alternative enzymes
substitute for the canonical ones.

### 5.3 Module completeness against MAG completeness

Status: done for one genome; two more would strengthen it. No status recorded for the
`benchmarks/configurations.toml` arm.

No paper quantifies module completeness error as a function of MAG completeness. A 70%-complete MAG is
missing 30% of its genes, so its module completeness is biased downward by an unknown amount. With
CAMI2's gold assembly and truth table, a known genome can be degraded to 95, 90, 80 and 70%
completeness and module completeness recomputed. The resulting curve gives each module call a
completeness-dependent confidence interval. Declare it as an arm in `benchmarks/configurations.toml`.

Degradation has to be done within one genome. Across five CAMI2 bins, the number of fully complete
modules tracks genome size more than CheckM2 completeness: a 100%-complete 2.0 Mbp bin has 23, a
97.3%-complete 4.7 Mbp bin has 44, and the two most fragmented bins (82% and 74% complete) have 9 each.

Measured on CAMI2 `bin_003` (4.66 Mbp, 58 complete modules intact), dropping whole contigs to reach
each target, three independent draws:

| MAG completeness | complete modules | spread | loss |
|---|---|---|---|
| 100% | 58.0 | - | - |
| 95% | 52.7 | 49–55 | 9% |
| 90% | 48.0 | 47–50 | 17% |
| 80% | 37.0 | 31–41 | 36% |
| 70% | 24.0 | 20–28 | 59% |
| 50% | 15.0 | 13–18 | 74% |

- Module loss is about twice the gene loss, because a module fails if any one of its genes is on a
  dropped contig.
- At the 70/10 threshold, 59% of complete modules are lost, so CheckM2 completeness goes in every
  row of `modules.tsv`.
- The mean pathwise score falls 25% across the range while the count of complete modules falls 74%.
  The count is reported, since it is the number usually cited.
- At 80% completeness the three draws gave 31, 39 and 41, a wider range than the gap between the 80%
  and 90% levels.

### 5.4 Out of scope

Status: decided against (both items).

Genome-scale metabolic models are not part of v1. On 70%-complete MAGs, gapfilling cannot tell a
missing gene from a real auxotrophy, and results change with the assumed medium. The best of the
three tools tested (gapseq: best false-negative rate, evidence-restricted gapfilling) still
over-predicted formate as a carbon source 29 times, and no paper quantifies draft-GEM accuracy against
MAG completeness. If offered later: require ≥90% completeness, record the medium, and never show
gapfilled reactions next to sequence-supported KOs without a provenance column.

HUMAnN reports community MetaCyc pathway abundance from reads against a pangenome catalogue, not
per-MAG capability. Its stratified-by-taxon view is covered in §8 with MMseqs2 taxonomy or
Woltka. HUMAnN 4 is alpha with no paper, and its main branch still ships MetaCyc v24 definitions, six
releases behind. Run it alongside MetaWrap2 for a read-based view; do not wrap it.

---

## 6. `amr`

### 6.1 Backends

Status: AMRFinderPlus, RGI (`--rgi`), hAMRonization and argNorm done. No status recorded for
ResFinder/PointFinder.

- AMRFinderPlus is the default: NCBI, quarterly database releases, curated per-gene thresholds, point
  mutations, and a `--plus` mode that writes 1,034 VIRULENCE entries in the same schema.
- RGI/CARD is the second backend. Every call maps to an ARO term, so results are comparable with the
  resistome literature. `--include_loose` is exposed and off by default; most of RGI's false positives
  are loose hits.
- ResFinder/PointFinder is a third, for acquired genes and resistance mutations.
- hAMRonization merges the three into one table. Its 35-field spec has no virulence fields, so VF
  results go to a separate table.
- argNorm maps every backend's calls onto ARO terms, so a gene found by two backends is counted once.

An audit of CARD's curated cutoffs found over-prediction rare (5 of 3,900 models above 5% false
positives) and under-prediction common (739 of 3,900 above 5% false negatives), concentrated in
β-lactamases, the largest ARG class. RGI's main error is missed calls.

### 6.2 Evidence levels (Cao et al.)

Status: done (Cao evidence levels, §11 phase 5).

Cao et al., *J Bacteriol* 2026, `10.1128/jb.00090-26` (CC BY), used as the schema:

| Level | Evidence | Permitted language |
|---|---|---|
| L1 Presence | ARG detected, with controls | "resistance determinants detected in the sample" |
| L2 Context | ARG linked to a host and genetic neighbourhood | "resistance potential" |
| L3 Activity | expression or protein evidence | "activity under the conditions tested" |
| L4 Phenotype | isolate AST, functional metagenomics, genetics | "confers resistance" |

MetaWrap2 reaches L1 and L2. Every ARG row has an `evidence_level`, and the report does not use L4
language for an L1 call. Required per-row provenance (Cao's L1 validity requirements): database name,
version and date; tool version; identity and coverage thresholds; the multi-mapping rule; the
normalisation denominator.

### 6.3 Contig edges

Status: the contig-edge rule is done. No status recorded for a BGC `contig_edge` flag.

Assemblies break at ARGs and mobile elements. Direct read mapping found 85 ARGs where individual
assemblers found 20–44, only 12% of ARG genomic contexts were assembled correctly, and MEGAHIT's
ARG-containing contigs averaged 284 bp, so a 300 bp length filter would discard 52% of them.

1. Annotate ARGs on contigs before any length filter.
2. Write `distance_to_nearest_contig_end` on every ARG call, and classify anything closer than the
   mobility window as `context_unknown` instead of `no_MGE`. "No mobile element within 10 kb" cannot be
   claimed when the ARG is 2 kb from a contig end. antiSMASH's `contig_edge` flag does the same, and
   BGC calls here should have it too.

### 6.4 Mobility and context

Status: geNomad mobility done (`amr --mobile-elements`, `--mobility-window` default 10 kb). No status
recorded for PlasmidFinder, ISEScan or mobileOG-db.

| Need | Tool | Note |
|---|---|---|
| plasmid and virus contigs | **geNomad** | best available caller, and one tool for both; shares its database with the `viruses` module in the companion plan |
| replicon typing | PlasmidFinder | a negative result means little on short contigs |
| insertion sequences | ISEScan | ships its own ~621 pHMMs |
| MGE families | mobileOG-db | its authors' rule: require multiple co-located hits, and "hits solely to RRR modules are not necessarily indicative of an MGE" |

Platon and MOB-suite are weaker alternatives if geNomad's database is unavailable. Plasmid calls are
reported as annotation, not reconstruction.

Published distance thresholds for "ARG on a mobile element" range from 1 kb (IS adjacency) through
10 kb flanks to same-contig co-location, and one paper calls its own 1 kb choice "somehow arbitrary".
The window is a parameter, default 10 kb, recorded in provenance, and the distance itself is reported.

### 6.5 Virulence

Status: AMRFinderPlus `--plus` is the virulence screen. VFDB decided against pending a measured arm.

VF results use the same search, context and quantification code with a different reference set. They
go to their own table, are limited to L1–L2 language ("virulence determinant detected"), and are never
aggregated into an AMR risk score. AMRFinderPlus `--plus` VIRULENCE comes with the AMR run; VFDB
setA/setB would add breadth.

VFDB distributes FASTA without per-gene thresholds, so a VFDB screen needs a chosen identity and
coverage floor, and the calls depend on it. AMRFinderPlus' virulence entries have NCBI's
per-reference cutoffs, so a row can be labelled `L1_presence`. VFDB would be added as a measured arm:
screen with both, count disagreements, and report the effect of the floor. Until then, VFDB's
coverage is a known gap. For comparison, the two resistance backends (`amr --rgi`) disagree on 6 of 8
genes on one bin.

VFDB headers use `WP_` protein accessions, the same namespace AMRFinderPlus writes, so hits on the same
protein join directly.

### 6.6 Prior art

Context.

MetaRanker (*AEM* 2026, `10.1128/aem.02422-25`, CC BY) combines read-level ARG/MGE/VF abundance with
contig co-occurrence into a risk index, validated on 353 metagenomes and more than 50% faster than
MetaCompare 2.0. Review it before implementing a risk index.

---

## 7. `quantify_function`

Status: done (`metawrap2.funcabund`), including the SingleM denominator and the assignable fraction.

### 7.1 Calculation

For gene family *g* in sample *s*, summing over genomes *k*:

```
A(g,s)  =  Σ_k  n(g,k) · abundance(k,s)
```

`n(g,k)` is the copy number of *g* in genome *k*; `abundance(k,s)` is that genome's abundance.

Do not divide a coverage-based abundance by genome size; CoverM's `relative_abundance` and
`trimmed_mean` are already proportional to cell abundance. Divide only read-count values (raw counts,
salmon `NumReads`). For MAGs, incompleteness cancels in `aligned_bases / MAG_length`; it biases the
gene-content term `n(g,k)`, not the abundance term. Put this in a code comment.

### 7.2 Units

| Unit | When | Caveat |
|---|---|---|
| copies per genome equivalent | comparable across samples | denominator from SingleM `prokaryotic_fraction`, not MicrobeCensus, which is abandoned (last release 2016) |
| length-normalised depth per Gb of bacterial/archaeal bases | robust to host or eukaryotic DNA | the denominator must exclude host, which SingleM provides |
| relative abundance / CPM | within-sample composition | compositional and sum-constrained; needs CoDA methods, never a raw t-test |

Length-normalise for any comparison between genes. Write at least two units plus the full provenance
tuple, and never divide one unit by another.

### 7.3 The functional hierarchy: gene → chemical process → pathway

Status: done (`metawrap2.hierarchy`, vendored BRITE projection, 1.1 MB gzipped), including Krona plots.

KEGG BRITE `ko00001` has four levels; the middle two are "chemical process" and "pathway":

| level | example | count |
|---|---|---|
| A | `Metabolism` | 8 |
| B (chemical process) | `Carbohydrate metabolism` | 57 |
| C (pathway) | `Glycolysis / Gluconeogenesis` | 570 |
| D (gene) | `K00844` hexokinase | 28,545 |

Two weightings:

- Per genome, every gene weighted equally (e.g. the fraction of the MAG's annotated genes in
  carbohydrate metabolism).
- Per sample, every genome weighted by abundance: the fraction of the community's functional capacity,
  which is not the mean of per-genome percentages. The abundance must be coverage-proportional and not
  divided by genome size again (§7.1).

Implementation:

- 12,025 of 28,545 KOs belong to more than one pathway, up to 128. Each gene's weight is split across
  its pathways, so every level sums to the input total and no percentage exceeds 100.
- Two of BRITE's eight top-level branches are not pathways. `09180 Brite Hierarchies` is the largest
  branch by a factor of two (22,768 of 64,956 rows) and holds groupings such as "Protein families:
  genetic information processing"; `09190 Not Included in Pathway or Brite` holds the remainder. Both
  are kept and flagged, because dropping them inflates the other percentages. Every row has two
  percentages: of all branches, and of the pathway branches only.

Krona: `ktImportText` is already a dependency (`kraken2` uses it); its input is a weight plus one
column per level. Two plots per unit: the full hierarchy, and pathway branches only, since the BRITE
branch otherwise dominates.

### 7.4 Assignable fraction

Status: assignable fraction done. Not recorded here for `singlem appraise` or Nonpareil; the companion
plan lists `singlem appraise` as open and Nonpareil as done (`metawrap2 nonpareil`).

Report the fraction of gene and ARG read abundance that can be assigned to a MAG or a classified
contig. In published work 29–41% of ARG abundance cannot be placed, and across 92,187 metagenomes only
about 10% of assembled base pairs end up in medium- or high-quality MAGs. The unbinned fraction is not
a random sample of the community.

Outputs: a sample-level read-based profile, a MAG-level per-genome profile, and the assignable
fraction. Never divide MAG counts into read abundance.

`singlem appraise` gives the per-cell version (fraction of the community the bins represent, and
lineages with no genome), shown next to a MAG-only functional profile.

Nonpareil (Artistic-2.0, small) in `read_qc` or `coverage` measures diversity coverage, a cheap
check before cross-sample comparison. Most published comparisons reportedly ignore it.

---

## 8. Bin-free path

Status: done, including quantification (`catalogue --quantify`, salmon). No status recorded for the
MMseqs2 licence-metadata question or per-gene taxonomy (step 5).

A user with an assembly and no bins gets the same functional outputs, as in SqueezeMeta and MGnify's
assembly pipeline.

1. Gene calling over all contigs: the shared pass from §3.
2. A gene catalogue across samples via MMseqs2 `linclust` (resolve its MIT-vs-GPLv3 metadata conflict
   first) at 95% identity and 90% coverage, keeping the longest representative. Write the catalogue and
   a member→representative map.
3. Annotate representatives once, not per sample.
4. Quantify by mapping reads to the catalogue (CoverM, or salmon/kallisto on the nucleotide catalogue),
   giving a gene × sample matrix in the units of §7.2.
5. Optionally assign taxonomy per gene with MMseqs2 taxonomy or Woltka, which gives HUMAnN's
   stratified-by-taxon view without HUMAnN's licence problem.

Binned and bin-free outputs share identifiers, units and provenance schema, so the fraction of
community function in the MAGs can be computed (§7.4).

---

## 9. Output contract

Status: no status recorded. `quantify_function --biom` exists, and the companion plan lists BIOM as
done.

Tidy long TSVs plus a joinable wide matrix, and BIOM 2.1 for feature tables (no other pipeline writes a
standard format).

```
FUNCTIONAL_ANNOTATION/
  genes/<genome>.{gff,faa,ffn}          per-genome structural annotation
  features.tsv                          gene × (KO, EC, COG, GO, Pfam, CAZy) long table
  modules.tsv                           genome × module, both completeness scores + the §5.2 flags
  amr.tsv                               hAMRonization schema + evidence_level + context columns
  virulence.tsv                         separate namespace, L1-L2 only
  bgc.tsv / cazy.tsv                    opt-in
  abundance/gene_{cpg,per_gb,relab}.tsv sample × gene family, one file per unit
  abundance/pathway_*.tsv               sample × module
  hierarchy/genome_hierarchy.tsv        genome × (level A, level B, pathway), both percentages
  hierarchy/sample_hierarchy.tsv        the same, abundance-weighted per sample
  hierarchy/krona/<unit>.html           interactive wheel per genome and per sample
  abundance/assignable_fraction.tsv     the §7.4 metric
  catalogue/                            bin-free path: representatives + member map
  provenance.json                       the tuple below, for every table
```

Provenance tuple for functional and pathway calls: tool, tool version, KOfam snapshot date,
module-definition snapshot date, threshold policy, completeness definition. For AMR: tool version,
database name with version and date, identity and coverage thresholds, multi-mapping rule,
normalisation unit and denominator, evidence level, decontamination state.

Published comparisons: ten pipelines on the same data differ 45-fold in functional output, two
resistome tools agree at 16% Jaccard, and threshold policy alone changes module counts by 11–12%.

---

## 10. Benchmarking

Status: item 1 done for one genome (§5.3). No status recorded for items 2 and 3.

There is no ground truth for what a genome can do: no CAMI functional-annotation track (absent from
CAMI I, II, the portal and CAMI III), no benchmark of eggNOG-mapper, KofamScan, InterProScan and DRAM on
the same genomes, and no published head-to-head of read-based functional profilers. What can be
measured:

1. Module completeness against MAG completeness (§5.3); not published elsewhere.
2. Annotation-layer concordance: KOfam with adaptive thresholds, eggNOG, and both, on CAMI2 genomes with
   known gene content. Reproduces a published result, then extends it to MAGs.
3. ARG recall from reads, assembly and MAGs, with the assignable-fraction metric. The only existing
   comparison is a preprint on read-based tools, whose best mean F1 under realistic uneven coverage was
   0.122.

No runtime or RAM figures are published for antiSMASH per genome, run_dbcan, GECCO, METABOLIC,
microTrait or AMRFinderPlus, so any figure quoted for them is an extrapolation, including the one in
§4.3.

---

## 11. Order of work

Each phase is usable before the next one exists.

| Phase | Deliverable | State |
|---|---|---|
| 1 | `annotate`: shared gene-calling pass over bins or an assembly, Pyrodigal, per-gene `partial` and `distance_to_contig_end`. Deprecate `annotate_bins` | **done** - `metawrap2.genecall`, `modules/annotate.py` |
| 2 | KO assignment: PyHMMER search, KOfam curated thresholds with the right score per profile, relative fallback, every call carrying its evidence | **done** - `metawrap2.kofam`; needs the 1.5 GB profile set to run on real data |
| 3 | `pathways`: vendored definitions, pathwise and stepwise completeness, all per-module flags | **done** - `metawrap2.kegg`, `modules/pathways.py` |
| 4 | the functional hierarchy and Krona plots, per genome and abundance-weighted per sample | **done** - `metawrap2.hierarchy` |
| 5 | `amr`: AMRFinderPlus, Cao evidence levels, the contig-edge rule | **done**: hAMRonization, argNorm, geNomad mobility (`--mobile-elements`) and RGI as a second backend (`--rgi`), joined on ARO terms - 2 of 8 genes shared on one bin |
| 6 | `quantify_function`: per-genome and per-sample abundance, SingleM denominator, assignable fraction | **done** - `metawrap2.funcabund` |
| 7 | bin-free path: `catalogue` plus annotation of representatives | **done**, including quantification (`catalogue --quantify`, salmon, gene abundance in the long format the rest of the pipeline reads) |
| 8 | breadth: dbCAN, antiSMASH, InterProScan, VFDB, geNomad mobility | dbCAN **done** (`--cazy`); BGCs **done** (`--bgc`, GECCO); geNomad mobility **done** (`amr --mobile-elements`); InterProScan **decided against** and Pfam built instead (`--pfam`, +1,578 genes KEGG missed) - see §4; VFDB **decided against** pending a measured arm, because it ships no per-gene thresholds - see §6.5 |
| 9 | the module-completeness-against-MAG-completeness benchmark | **done for one genome** - see §5.3; two more would strengthen it |

Phases 1–6 cover genome-resolved annotation, 7 the bin-free path, 8 breadth, 9 the new measurement.

Phase 8 notes:

- BGCs: GECCO is the default: a conditional random field over Pyrodigal genes and PyHMMER domains,
  no multi-gigabyte database, minutes per genome. antiSMASH takes tens of minutes per genome and a
  10 GB install, and is better for studying a specific cluster; GECCO is practical across hundreds of
  MAGs.
- InterProScan: described here as "deferred rather than rejected" (the table and §4.3 say "decided
  against"). A ~200 GB install here (§4.3 says ~14 GB on disk) that overlaps heavily with KOfam and
  dbCAN for a prokaryotic MAG; the additional Pfam domain architecture is available with PyHMMER
  against Pfam. Revisit as `--pfam`.
- VFDB: AMRFinderPlus `--plus` already writes virulence entries in the same schema, to their own
  table. A second virulence database would need its own normalisation and a measurement.

Development uses a three-bin fixture (`devset/`, 6 MB, 6.2 Mbp: one bin of 92 contigs and two of ~900,
so fragmented bins are included). Gene calling over all three takes under a second.
