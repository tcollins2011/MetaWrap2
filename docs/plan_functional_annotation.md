# Plan: functional annotation, pathways, AMR, and abundance

A plan to replace `annotate_bins`, a Prokka/Bakta wrapper that emits GFF, `.faa` and `.ffn`, with a
subsystem that answers four questions:

1. What genes does this genome carry? (structural and functional annotation)
2. What can it do? (orthology, then pathway/module completeness)
3. What resistance and virulence does it carry, and in what genetic context?
4. How much of each function is in the sample? (gene and pathway abundance, per genome and
   community-wide, with and without bins)

Each section starts with its status; §11 has the state of every phase. "No status recorded" means this
plan does not say whether the item is done. Numbers attributed to other tools are evidence about those
tools, not targets for MetaWrap2.

---

## 1. What is there now, and what the field has

Context.

`annotate_bins` runs one tool per bin and reorganises its output: no identifiers joinable across
genomes, no pathway inference, no AMR, no abundance, and nothing for a user who did not bin.

Peers: SqueezeMeta emits KEGG/COG/Pfam abundance tables for reads and bins with native R bindings. DRAM
distils KEGG/Pfam/dbCAN/MEROPS into a metabolism summary. VEBA scores BGC novelty. MGnify ships a
mobilome pipeline. nf-core/funcscan screens AMR/AMP/BGC/CAZy on contigs. MetaWrap2 ships gene calling.

Permissively licensed implementations to read first:

- mettannotator (Apache-2.0, EMBL-EBI, validated on MAGs), the closest analogue.
- MGnify's `assembly-analysis-pipeline`: the bin-free path, with exact tool and database versions.
- nf-core/funcscan 4.0.0: database handling and version pinning for the AMR and BGC arms.

---

## 2. Licensing

Context and two rules (no status recorded for the rules).

For non-profit academic use every tool is available: CARD and RGI, geNomad, VFDB, KEGG's KOfam profiles
and module definitions, MetaCyc through HUMAnN. Tools are chosen on evidence.

- Do not call `rest.kegg.jp` from a pipeline run. It is rate-limited to 3 requests/s. Module
  definitions come from a vendored, dated snapshot.
- Pin and record every database version; a KOfam or CARD snapshot moves results as much as a tool
  version. `archives/` on the KEGG FTP has dated KOfam snapshots.

Commercial redistribution would need to revisit CARD/RGI (prohibited for commercial organisations
without McMaster's permission), geNomad (LBNL academic-use licence), VFDB and integron_finder's data
files (CC BY-NC), and MetaCyc (paid).

The spine: Pyrodigal → Bakta → PyHMMER on KOfam/Pfam/NCBIfam → eggNOG-mapper →
kegg-pathways-completeness → AMRFinderPlus + RGI + hAMRonization → dbCAN → CoverM/SingleM.

---

## 3. Architecture

Status: done; all five modules exist (§11, phases 1–7). No status recorded for the `doctor` database
probes.

One module becomes five, none requiring that binning happened.

```
                  ┌─ annotate (genes + function) ──┬── per-genome tables
assembly ─┬─ bins ─┤                                │
          │        ├─ pathways (modules/completeness)┤
          │        ├─ amr (ARG/VF + context)        ├── join on protein id
          │        └─ bgc/cazy (opt-in)             │
          └─ contigs (bin-free) ── gene catalogue ──┴── abundance ── sample × feature matrices
```

| Module | Input | Answers |
|---|---|---|
| `annotate` | bins or a whole assembly | genes, products, orthologs, domains |
| `pathways` | `annotate` output | module/pathway completeness per genome |
| `amr` | bins or contigs | ARGs, VFs, and their genetic context |
| `quantify_function` | `annotate` + reads | gene and pathway abundance, per genome and per sample |
| `catalogue` | contigs from any number of samples | non-redundant gene catalogue for the bin-free path |

`annotate_bins` stays as a deprecated alias that forwards to `annotate` with a warning, because it is in
published scripts.

Gene prediction and the expensive homology searches run once over all contigs, and every downstream
consumer scores groupings of the cached features (Binette does the same for refinement). `annotate_bins`
re-runs Prokka per bin, re-predicting genes for shared contigs. With a shared front end, the bin-free
path is nearly free once the binned path has run, and the reverse.

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

`doctor` gains a probe per new database (KOfam snapshot, Pfam, eggNOG 5, AMRFinderPlus DB,
module-definition snapshot), reporting installed / missing / broken with version and date. The express
installer must not auto-download any of them; `doctor --fix` prompts.

---

## 4. `annotate`

### 4.1 Gene calling

Status: Pyrodigal with per-gene `partial` and `distance_to_contig_end` is done (`metawrap2.genecall`).
No status recorded for a second caller.

Pyrodigal (MIT, Cython bindings to Prodigal, maintained) is importable, releases the GIL, and removes a
subprocess. Prodigal is effectively frozen.

- MAG contigs are short, so genes running off a contig end are called partial or missed. `partial` and
  `distance_to_contig_end` are recorded per gene and propagated, because every downstream claim (a
  module step, an ARG, a BGC border) inherits that uncertainty.
- Fragmented genes are better served by a second caller; MGnify's assembly pipeline merges two.
  FragGeneScanRs is the candidate. The merge rule must not count a gene twice.

### 4.2 Functional annotation

Status: open. The Bakta pin is not done (`bakta >=1.9`). No status recorded for the light-database
default.

Bakta (GPL-3, actively developed) stays the per-genome annotator. It runs tRNAscan-SE, Aragorn and
Infernal/Rfam internally, so the rRNA and tRNA counts MIMAG requires are a parse of existing output; that
belongs to bin QC in the companion plan.

Bakta's full database is ~70 GB against `--db-light`'s ~2 GB; the extra is mostly UniRef cluster coverage,
which improves product names, not the identifiers modules join on. Default to light, make full a flag,
and record which ran: genomes annotated against different Bakta databases are not comparable on product
names.

`envs/annotate.yaml` asks for `bakta >=1.9`, which resolves to 1.12.x, so a rerun months later can give
different products with no record of why. Pin it exactly and record the version, as `provenance` does
for databases. The same applies to `gtdbtk >=2.4.0` (companion plan §1.5).

### 4.3 Orthology and domains

Status: KO assignment done (`metawrap2.kofam`; needs the 1.5 GB profile set to run on real data). Pfam
done (`annotate --pfam`). No status recorded for eggNOG-mapper. InterProScan is inconsistent: "opt-in"
and "decided against" in this section, "deferred rather than rejected" in §11.

Most of the error budget is here. On 396 genomes with identical gene calls (Kananen et al.,
*Bioinform Adv* 2025, `10.1093/bioadv/vbaf039`):

| | calls confirmed by eggNOG | tool-specific calls confirmed | tool-specific calls contradicted |
|---|---|---|---|
| KofamScan | 84% | — | — |
| anvi'o (adaptive thresholds) | 79% | 37% | 7.5% |
| MicrobeAnnotator | 74% | 0.3% | 35% |

anvi'o's adaptive bit-score heuristic annotates 13% more genes per genome (p<2.2e-16), mostly in
under-represented taxa: 38% of Nanosynbacteraceae cases against 4.9% of Enterobacteriaceae. Disabling it
reproduces KofamScan exactly.

KOs are assigned with PyHMMER over KOfam profiles, with an anvi'o-style adaptive threshold where the
curated cutoff is missing or fails. MicrobeAnnotator's permissive path is excluded (35% of its unique
calls contradicted). 1,866 of 28,453 KOfam profiles (6.6%) have no bit-score threshold; drop them or give
them a conservative substitute, and record the choice, since they are a documented source of false
positives.

PyHMMER (MIT) runs in-process instead of `hmmsearch` as a subprocess, and keeps the threshold policy in
code. No published hmmsearch-vs-Pfam timing at ~1M proteins exists, so compute for this step must be
measured before the module ships.

eggNOG-mapper is the second opinion and the source of COG/GO/EC/CAZy. Pin 2.1.15 with eggNOG 5:
v3.0.0-beta6 is on eggNOG 7, incompatible with eggNOG 5 data, and a beta. Its orthology-transferred KOs
ignore KOfam's per-KO cutoffs and reflect the KEGG version eggNOG was built on, so eggNOG validates and
does not assign.

InterProScan 6: opt-in only, on cost (~14 GB on disk, too slow for every MAG in a cohort).

Decided against InterProScan; Pfam built instead. `annotate --pfam` searches Pfam-A with PyHMMER at each
family's gathering threshold, like KOfam's per-profile cutoffs. It is 399 MB and in-process. On one CAMI2
bin, KEGG placed 1,620 of 3,891 genes (41.6%), Pfam placed 3,185 (81.9%), and 1,578 of those have no
KEGG ortholog.

Lost: InterProScan's other members (TIGRFAM, SMART, PANTHER, SUPERFAMILY) and the signal-peptide and
transmembrane predictors. This pipeline does not produce GO terms or transmembrane topology. The bioconda
build is also four years old (5.54_87.0, 2021).

---

## 5. `pathways`

### 5.1 Five tools, five definitions of "complete"

Context.

| Tool | Semantics | Bias |
|---|---|---|
| anvi'o | unrolls every alternative path; complexes score as the fraction of *essential* subunits; `--` steps always 0; module = max over paths | optimistic by construction; the only one that also offers a stepwise score |
| EBI kegg-pathways-completeness | definition → weighted DiGraph, traversal minimising weight ratio; `--` substituted at weight 0 | a `--` step does not cap the module, unlike anvi'o |
| MicrobeAnnotator | three separate algorithms; complexes get proportional partial credit | most permissive on complexes |
| DRAM | a node survives if any of its KOs is present; only 13 hardcoded modules | most permissive overall |
| METABOLIC | substitutes 1/0 and `eval`s the boolean; present if ≥0.75 of steps | module table dated 2023-06-12 |

On the same KO set, anvi'o reports 11–12% more modules at ≥80% completeness than MicrobeAnnotator or
KofamScan, so comparing MetaWrap2's number with a paper's compares algorithms.

Structural traps in the 573-module EBI snapshot:

- 47 modules (8.2%) have a single top-level step and 164 (28.6%) have two or fewer, so completeness is
  effectively binary and prone to false positives. M00851 "Carbapenem resistance" is a flat OR of 21
  KOs: any one hit makes it 100% complete, and one of the 21 has no threshold. In the benchmark above
  MicrobeAnnotator called it complete in every bacterial family but two.
- 14 modules contain a `--` step (a reaction with no KO) and cannot reach 100% under anvi'o.
- 45 KOs appear in five or more modules (K00024 in ten), inflating unrelated modules.
- 8 modules reference other modules and need a second resolution pass.
- Definitions drift: 573 in the pinned EBI snapshot against 582 live.

### 5.2 Design

Status: done (`metawrap2.kegg`, `modules/pathways.py`).

Completeness is computed in MetaWrap2 from vendored definitions. Brute-force path unrolling is
tractable: 3 of 573 modules exceed 1,000 alternative routes, the maximum is 33,600 for M00001, and all
573 score in 0.1 s.

A `+`-joined protein complex is one step whose subunits are all essential. M00142, NADH dehydrogenase,
has seven subunits; read as seven steps, a genome with three scores 3/7 when the answer is 0. 31 of the
573 modules are a single complex. Counting top-level spaces in the definition string reaches the
published "47 single-step modules" through two errors that cancel (a complex read as a sequence, a
parenthesised group read as one step). The parser distinguishes both and reproduces 47 correctly.

Both scores are always emitted, `completeness_pathwise` (max over alternative paths) and
`completeness_stepwise` (fraction of top-level steps satisfied), plus the path chosen. They can disagree
by as much as the biological signal.

| Flag | Meaning |
|---|---|
| `n_top_level_steps` | lets a reader discount the 8.2% of modules that are effectively binary |
| `is_signature_module` | AMR and signature modules are not metabolic capability and must not be reported as such |
| `has_no_ko_step` | the module cannot reach 100% under pathwise semantics |
| `proportion_unique_enzymes` | anvi'o's mitigation for promiscuous KOs, and the only existing defence |
| `nt_ko_calls` | how many supporting KOs had no curated threshold |
| `definition_snapshot` | the vendored module-table date |

Vendor the EBI `modules_table.tsv` and `graphs.pkl` (Apache-2.0, dated), pin the version, emit the
snapshot date on every row, and never mix snapshots within a run.

Call the output metabolic potential. Completeness does not predict phenotype: across 56 strains in
gnotobiotic honey bees, completeness "did not fully explain colonization patterns", and a complete
methylmalonyl-CoA pathway was never found in Cloacimonadota because alternative enzymes substitute for
the canonical ones.

### 5.3 Module completeness against MAG completeness

Status: done for one genome; two more would strengthen it. No status recorded for the
`benchmarks/configurations.toml` arm.

No paper quantifies module completeness error as a function of MAG completeness. A 70%-complete MAG is
missing 30% of its genes, so every module it reports is biased downward by an unknown amount. CAMI2's
gold assembly and truth table make this measurable: degrade a known genome to 95, 90, 80 and 70%
completeness, recompute module completeness, and report the curve, which would give every module call a
completeness-aware confidence interval. Declare it as an arm in `benchmarks/configurations.toml`.

Degradation has to be within one genome. On five CAMI2 bins the count of fully complete modules tracks
genome size, not CheckM2 completeness: a 100%-complete 2.0 Mbp bin has 23, a 97.3%-complete 4.7 Mbp bin
has 44, and the two most fragmented bins (82% and 74% complete) have 9 each.

Measured with `bench_work/completeness_degradation.py` on CAMI2 `bin_003` (4.66 Mbp, 58 complete modules
intact), dropping whole contigs to each target, three independent draws:

| MAG completeness | complete modules | spread | loss |
|---|---|---|---|
| 100% | 58.0 | — | — |
| 95% | 52.7 | 49–55 | 9% |
| 90% | 48.0 | 47–50 | 17% |
| 80% | 37.0 | 31–41 | 36% |
| 70% | 24.0 | 20–28 | **59%** |
| 50% | 15.0 | 13–18 | 74% |

- The loss is amplified about two-fold: a module fails if any one of its genes is on a dropped contig.
- At the 70/10 threshold, 59% of complete modules are gone. A bin that just passes reports two fifths of
  what the intact genome reports, so CheckM2 completeness belongs in every row of `modules.tsv`.
- The mean pathwise score falls 25% across the range while the count of complete modules falls 74%.
  Report the count; it is what people quote.
- At 80% completeness three draws gave 31, 39 and 41, wider than the gap between the 80% and 90% levels.
  Which contigs are lost matters nearly as much as how many.

Write-up and caveats: `bench_work/degradation/FINDINGS.md`. One genome so far.

### 5.4 What not to build

Status: decided against (both items).

Genome-scale metabolic models are not a v1 layer. On 70%-complete MAGs, gapfilling makes a missing gene
indistinguishable from a real auxotrophy, and results move with the assumed medium. The best of the three
tools (gapseq: best false-negative rate, evidence-restricted gapfilling) still over-predicted formate as
a carbon source 29 times, and no paper quantifies draft-GEM accuracy against MAG completeness. If offered:
gate on ≥90% completeness, record the medium, and never show gapfilled reactions beside sequence-supported
KOs without a provenance column.

HUMAnN reports community MetaCyc pathway abundance from reads against a pangenome catalogue; it does not
say what any MAG can do. Its stratified-by-taxon view is matched in §8 with MMseqs2 taxonomy or Woltka.
HUMAnN 4 is alpha with no paper, and master still ships MetaCyc v24 definitions, six releases behind. Run
it alongside MetaWrap2 for a read-based view; do not wrap it.

---

## 6. `amr`

### 6.1 Backends

Status: AMRFinderPlus, RGI (`--rgi`), hAMRonization and argNorm done. No status recorded for
ResFinder/PointFinder.

- AMRFinderPlus is the default: NCBI, quarterly database releases, curated per-gene thresholds, point
  mutations, and a `--plus` mode that emits 1,034 VIRULENCE entries in the same schema.
- RGI/CARD is the second backend. Every call lands on an ARO term, which makes results comparable to the
  resistome literature. `--include_loose` is exposed and off by default; loose hits hold its false
  positives.
- ResFinder/PointFinder is a third, for acquired genes and resistance mutations.
- hAMRonization normalises the three into one table. Its 35-field spec has no virulence concept, so VF
  results go to their own table.
- argNorm maps every backend's call onto ARO terms, so a gene found by two backends is one gene.

An audit of CARD's curated cutoffs found over-prediction rare (5 of 3,900 models above 5% false
positives) and under-prediction common (739 of 3,900 above 5% false negatives), concentrated in
β-lactamases, the largest ARG class. RGI's dominant error is the missed call.

### 6.2 Output contract: the Cao et al. evidence ladder

Status: done (Cao evidence levels, §11 phase 5).

Cao et al., *J Bacteriol* 2026, `10.1128/jb.00090-26` (CC BY), used as the schema:

| Level | Evidence | Permitted language |
|---|---|---|
| L1 Presence | ARG detected, with controls | "resistance determinants detected in the sample" |
| L2 Context | ARG linked to a host and genetic neighbourhood | "resistance potential" |
| L3 Activity | expression or protein evidence | "activity under the conditions tested" |
| L4 Phenotype | isolate AST, functional metagenomics, genetics | "confers resistance" |

MetaWrap2 reaches L1 and L2. Every ARG row carries `evidence_level`, and the report must not use L4
language for an L1 call. Required per-row provenance (Cao's L1 validity requirements): database name,
version and date; tool version; identity and coverage thresholds; the multi-mapping rule; the
normalisation denominator.

### 6.3 Contig edges

Status: the contig-edge rule is done. No status recorded for a BGC `contig_edge` flag.

Assemblies break at ARGs and mobile elements. Direct read mapping found 85 ARGs where individual
assemblers found 20–44; only 12% of ARG genomic contexts were captured correctly; MEGAHIT's
ARG-containing contigs averaged 284 bp, so a 300 bp length filter would discard 52% of them.

1. Annotate ARGs on contigs before any length filter.
2. Emit `distance_to_nearest_contig_end` on every ARG call, and classify anything closer than the
   mobility window as `context_unknown`, not `no_MGE`: "no mobile element within 10 kb" is unsupported
   when the ARG sits 2 kb from a contig end. antiSMASH's `contig_edge` flag does the same, and BGC calls
   here should too.

### 6.4 Mobility and context

Status: geNomad mobility done (`amr --mobile-elements`, `--mobility-window` default 10 kb). No status
recorded for PlasmidFinder, ISEScan or mobileOG-db.

| Need | Tool | Note |
|---|---|---|
| plasmid and virus contigs | **geNomad** | the best caller available, and one tool for both; shares its database with the `viruses` module in the companion plan |
| replicon typing | PlasmidFinder | absence proves nothing on short contigs |
| insertion sequences | ISEScan | ships its own ~621 pHMMs |
| MGE families | mobileOG-db | its authors' own rule: require multiple co-located hits, and "hits solely to RRR modules are not necessarily indicative of an MGE" |

Platon and MOB-suite are weaker alternatives if geNomad's database is unavailable. Plasmid calls are
reported as annotation, not reconstruction.

Published distance thresholds for "ARG on a mobile element" run from 1 kb (IS adjacency) through 10 kb
flanks to same-contig co-location, and one paper calls its own 1 kb choice "somehow arbitrary". The
window is a parameter, default 10 kb, recorded in provenance, and the distance itself is reported.

### 6.5 Virulence

Status: AMRFinderPlus `--plus` is the virulence screen. VFDB decided against pending a measured arm.

VF results use the same search, context and quantification engine with a different reference set. They
go to their own table, capped at L1–L2 language ("virulence determinant detected"), and are never
aggregated into an AMR risk score. AMRFinderPlus `--plus` VIRULENCE comes with the AMR run; VFDB
setA/setB would add breadth.

VFDB distributes FASTA without per-gene thresholds, so a VFDB screen means choosing an identity and
coverage floor, and that choice decides the call. AMRFinderPlus' virulence entries carry NCBI's
per-reference cutoffs, which is what lets a row say `L1_presence`. To add VFDB, run it as a measured arm:
screen both, count disagreements, report what a floor is worth. Until then its breadth is a stated gap.
The two resistance backends (`amr --rgi`) disagree on 6 of 8 genes on one bin, which shows both the value
of a second reference set and the need to measure it.

VFDB headers carry `WP_` protein accessions, the namespace AMRFinderPlus emits, so hits on the same
protein join directly.

### 6.6 Prior art

Context.

MetaRanker (*AEM* 2026, `10.1128/aem.02422-25`, CC BY) combines read-level ARG/MGE/VF abundance with
contig co-occurrence into a risk index, validated on 353 metagenomes and more than 50% faster than
MetaCompare 2.0. Read it before building.

---

## 7. `quantify_function`

Status: done (`metawrap2.funcabund`), including the SingleM denominator and the assignable fraction.

Per-genome gene and pathway abundance, extended to the whole sample.

### 7.1 The arithmetic

For gene family *g* in sample *s*, summing over genomes *k*:

```
A(g,s)  =  Σ_k  n(g,k) · abundance(k,s)
```

`n(g,k)` is the copy number of *g* in genome *k*; `abundance(k,s)` is that genome's abundance.

Do not divide a coverage-based abundance by genome size. CoverM's `relative_abundance` and `trimmed_mean`
are already proportional to cell abundance. Divide only read-count numbers (raw counts, salmon
`NumReads`). For MAGs, incompleteness cancels in `aligned_bases / MAG_length`; it biases the
gene-content term `n(g,k)`, not the abundance term. Put this in a code comment.

### 7.2 Units

| Unit | When | Caveat |
|---|---|---|
| copies per genome equivalent | the cross-sample comparable unit | denominator from SingleM `prokaryotic_fraction`, not MicrobeCensus, which is abandoned (last release 2016) |
| length-normalised depth per Gb of bacterial/archaeal bases | robust to host or eukaryotic DNA | the denominator must exclude host, which is what SingleM gives |
| relative abundance / CPM | within-sample composition | compositional and sum-constrained; needs CoDA methods, never a raw t-test |

Length-normalise for any comparison between genes. Emit at least two units plus the full provenance
tuple, and never divide one unit into another.

### 7.3 The functional hierarchy: gene → chemical process → pathway

Status: done (`metawrap2.hierarchy`, vendored BRITE projection, 1.1 MB gzipped), including Krona plots.

KEGG BRITE `ko00001` has four levels; the middle two are "chemical process" and "pathway":

| level | example | count |
|---|---|---|
| A | `Metabolism` | 8 |
| B — the chemical process | `Carbohydrate metabolism` | 57 |
| C — the pathway | `Glycolysis / Gluconeogenesis` | 570 |
| D — the gene | `K00844` hexokinase | 28,545 |

Two weightings:

- Per genome, every gene weighted equally: the fraction of this MAG's annotated genes in, say,
  carbohydrate metabolism. Describes the organism.
- Per sample, every genome weighted by abundance: the fraction of the community's functional capacity.
  This is not the mean of per-genome percentages. The abundance must be coverage-proportional and not
  divided by genome size again (§7.1).

Implementation:

- 12,025 of 28,545 KOs belong to more than one pathway, up to 128. Each gene's weight is split across its
  pathways, so every level sums to the input total and percentages do not exceed 100.
- Two of BRITE's eight top-level branches are not biology. `09180 Brite Hierarchies` is the largest
  branch by a factor of two (22,768 of 64,956 rows) and holds groupings such as "Protein families:
  genetic information processing"; `09190 Not Included in Pathway or Brite` is the leftovers. Both stay
  and are flagged, since dropping them inflates everything else. Every row carries two percentages: of
  everything, and of the pathway branches only.

Krona: `ktImportText` is already a dependency (`kraken2` uses it), and its input is a weight plus one
column per level. Two plots per unit: the full hierarchy, and pathway branches only, since the BRITE
branch otherwise dominates the wheel.

### 7.4 Report what fraction is assignable

Status: assignable fraction done. Not recorded here for `singlem appraise` or Nonpareil; the companion
plan lists `singlem appraise` as open and Nonpareil as done (`metawrap2 nonpareil`).

Report the fraction of gene and ARG read-abundance assignable to a MAG or a classified contig. In
published work 29–41% of ARG abundance cannot be placed, and across 92,187 metagenomes only about 10% of
assembled base pairs end up in medium- or high-quality MAGs. The unbinned remainder is biased, not
random.

Outputs: a sample-level read-based profile, a MAG-level per-genome profile, and the assignable fraction.
Never divide MAG counts into read abundance.

`singlem appraise` gives the per-cell version: what fraction of the community the bins represent, and
which lineages have no genome. It belongs beside a MAG-only functional profile.

Nonpareil (Artistic-2.0, small) in `read_qc` or `coverage` measures diversity coverage, which confounds
cross-sample comparison and is reportedly ignored by most published comparisons. It answers "can I
compare these samples?" cheaply.

---

## 8. The bin-free path

Status: done, including quantification (`catalogue --quantify`, salmon). No status recorded for the
MMseqs2 licence-metadata question or per-gene taxonomy (step 5).

A user with an assembly and no bins gets the same functional answers, as in SqueezeMeta and MGnify's
assembly pipeline.

1. Gene calling over all contigs: the shared pass from §3.
2. A gene catalogue across samples via MMseqs2 `linclust` (resolve its MIT-vs-GPLv3 metadata conflict
   first) at 95% identity and 90% coverage, keeping the longest representative. Emit the catalogue and a
   member→representative map.
3. Annotate representatives once, not per sample.
4. Quantify by mapping reads to the catalogue (CoverM, or salmon/kallisto on the nucleotide catalogue),
   giving a gene × sample matrix in the units of §7.2.
5. Optionally attribute taxonomy per gene with MMseqs2 taxonomy or Woltka, which gives the
   stratified-by-taxon view HUMAnN provides without HUMAnN's licence problem.

Binned and bin-free answers share identifiers, units and provenance schema, so "what fraction of the
community's function did my MAGs capture" is answerable (§7.4).

---

## 9. Output contract

Status: no status recorded. `quantify_function --biom` exists, and the companion plan lists BIOM as
done.

Tidy long TSVs plus a joinable wide matrix, and BIOM 2.1 for feature tables (no peer emits a standard
object).

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

Ten pipelines on the same data disagree 45-fold on functional output, two resistome tools agree at 16%
Jaccard, and threshold policy alone moves module counts by 11–12%.

---

## 10. Benchmarking

Status: item 1 done for one genome (§5.3). No status recorded for items 2 and 3.

There is no ground truth for "what can this genome do": no CAMI functional-annotation track (absent from
CAMI I, II, the portal and CAMI III), no benchmark of eggNOG-mapper, KofamScan, InterProScan and DRAM on
the same genomes, and no published head-to-head of read-based functional profilers. Measurable:

1. Module completeness against MAG completeness (§5.3). Nothing comparable exists.
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

Ordered so each stage is usable before the next exists.

| Phase | Deliverable | State |
|---|---|---|
| 1 | `annotate`: shared gene-calling pass over bins or an assembly, Pyrodigal, per-gene `partial` and `distance_to_contig_end`. Deprecate `annotate_bins` | **done** — `metawrap2.genecall`, `modules/annotate.py` |
| 2 | KO assignment: PyHMMER search, KOfam curated thresholds with the right score per profile, relative fallback, every call carrying its evidence | **done** — `metawrap2.kofam`; needs the 1.5 GB profile set to run on real data |
| 3 | `pathways`: vendored definitions, pathwise *and* stepwise completeness, all per-module flags | **done** — `metawrap2.kegg`, `modules/pathways.py` |
| 4 | the functional hierarchy and Krona plots, per genome and abundance-weighted per sample | **done** — `metawrap2.hierarchy` |
| 5 | `amr`: AMRFinderPlus, Cao evidence levels, the contig-edge rule | **done**: hAMRonization, argNorm, geNomad mobility (`--mobile-elements`) and RGI as a second backend (`--rgi`), joined on ARO terms - 2 of 8 genes shared on one bin |
| 6 | `quantify_function`: per-genome and per-sample abundance, SingleM denominator, assignable fraction | **done** — `metawrap2.funcabund` |
| 7 | bin-free path: `catalogue` plus annotation of representatives | **done**, including quantification (`catalogue --quantify`, salmon, gene abundance in the long format the rest of the pipeline reads) |
| 8 | breadth: dbCAN, antiSMASH, InterProScan, VFDB, geNomad mobility | dbCAN **done** (`--cazy`); BGCs **done** (`--bgc`, GECCO); geNomad mobility **done** (`amr --mobile-elements`); InterProScan **decided against** and Pfam built instead (`--pfam`, +1,578 genes KEGG missed) - see §4; VFDB **decided against** pending a measured arm, because it ships no per-gene thresholds - see §6.5 |
| 9 | the module-completeness-against-MAG-completeness benchmark | **done for one genome** — see §5.3; two more would strengthen it |

Phases 1–6 cover genome-resolved annotation, 7 the bin-free path, 8 breadth, 9 the new measurement.

Phase 8 notes:

- BGCs: GECCO is the default. It is a conditional random field over Pyrodigal genes and PyHMMER domains,
  needs no multi-gigabyte database, and runs in minutes per genome; antiSMASH takes tens of minutes and a
  10 GB install. antiSMASH is the deeper tool when a cluster is the result; GECCO makes "how many BGCs
  are in these 200 MAGs" answerable.
- InterProScan: described here as "deferred rather than rejected" (the table and §4.3 say "decided
  against"). A ~200 GB install here (§4.3 says ~14 GB on disk) that overlaps heavily with KOfam and
  dbCAN for a prokaryotic MAG; the additional Pfam domain architecture is reachable with PyHMMER against
  Pfam. Revisit as `--pfam`.
- VFDB: AMRFinderPlus `--plus` already emits virulence entries in the same schema, written to their own
  table. A second virulence database would need its own normalisation, with no measurement behind it.

Development uses a three-bin fixture (`devset/`, 6 MB, 6.2 Mbp: one bin of 92 contigs and two of ~900,
so the fragmentation cases are present). Gene calling over all three takes under a second.
