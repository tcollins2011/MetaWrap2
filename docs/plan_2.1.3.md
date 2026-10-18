# 2.1.3 plan

This is the plan for 2.1.3. For what shipped, see the [2.1.3 release notes](releases/2.1.3.md).

2.1.2 added thirteen modules, measured on one body site of one simulated dataset. 2.1.3 adds no
modules; it benchmarks more body sites and adds a PDF report.

Decisions already made:

- The release is 2.1.3, not 2.2.0.
- Long reads and HiFi are deferred (§7). Few users need them and the work is large.
- MIMAG tagging resumes (§3).
- The main result is more CAMI2 body sites, scored against metaWRAP 1 (§1).
- The report is a PDF instead of a set of TSVs (§2).
- Bracken runs after Kraken2 (§3b), because Kraken2 read counts are not abundances.

## 1. More CAMI2 body sites `[moderate]`

All defaults were measured on three CAMI2 gastrointestinal samples and may be gut-specific.

**Available data** (checked). The deposit our gastrointestinal data comes from (`frl:6425518`,
CC BY 4.0) has 49 short-read samples across five body sites, with per-contig gold standards, in the
same layout and under the same licence:

| Directory in the deposit | Samples | Body sites |
|---|---:|---|
| `gastrooral/` | 20 | gastrointestinal (the 10 in `datasets.toml`) and oral (the other 10) |
| `airskinurogenital/` | 29 | airways, skin, urogenital |
| `gastrooral_pbsim/`, `airskinurogenital_pbsim/` | 49 | the long-read counterparts of both, for later |

Archives are 9.7–10.3 GB each (from HEAD requests; recorded in `datasets.toml`). Extraction keeps
only reads, contigs and truth tables, about 5.7 GB per sample.

**Method (unchanged).** The existing comparison used gastrointestinal samples 0, 1 and 10:
per-sample assembly, metaWRAP 1.3.0 and MetaWrap2's `fast`, `default` and `slow` presets, all scored
by the same CheckM2 ∩ GUNC join. Each new body site gets three samples, the same four arms and the
same scoring.

**Deliverables.**

- `cami2_oral` and `cami2_skin` registered in `benchmarks/datasets.toml`, with md5s pinned from the
  publisher's manifest.
- All arms run through `metawrap2 benchmark run`, so CPU-hours are measured from the child process.
- `mw1_vs_mw2.png` and `headline_mw1_vs_mw2.png` get one panel per body site; a new figure shows
  clean HQ MAGs and CPU-hours per body site.
- The README comparison table gets one row per body site, whatever the margin.

**Cost.** The metaWRAP 1 arm took 132.0 CPU-hours on three gastrointestinal samples; the three
MetaWrap2 presets together take about the same. Budget ~220 CPU-hours per body site, one site at a
time (only one benchmark runs at once). Disk: ~30 GB of archives, deleted after extraction, ~17 GB
extracted, plus benchmark output. Oral first; skin only if oral results call for it.

**Stopping rule.** If oral matches gastrointestinal within noise on every arm, skip the third site
(220 CPU-hours) and report that the defaults hold across two body sites.

*Outcome:* oral gave a larger margin than gastrointestinal (+13% against +10%), and `--preset slow`
recovered more than the gastrointestinal preset table predicts, so skin was run as well. Skin gave
the largest margin (+22%) and the largest cost increase over metaWRAP 1 (47%). All three sites are in
[benchmarks.md §2.1](benchmarks.md) and [releases/2.1.3.md §11](releases/2.1.3.md).

## 2. PDF report `[moderate]`

A finished run writes tables but no single summary. Three parts, in order:

**2a. Per-MAG summary table `[low]`.** One row per genome: quality tier, taxonomy, size, N50,
coverage per sample, rRNA and tRNA counts, and any module flags. Every column already exists in
some TSV, so this is a join. The report is built around it.

**2b. Static wheels of the functional hierarchy `[low]`: done.** `ktImportText` writes only
interactive HTML, so there was no static image of the functional hierarchy for documents.
Rasterising Krona would need a headless browser. Instead, matplotlib draws a nested sunburst over
levels A/B/C from the same hierarchy table (deterministic, no new dependency). Krona's HTML is still
written. One PNG per genome and per sample weighting.

**2c. The PDF `[moderate]`.** One document per run: provenance (versions, databases, commands), read
QC, assembly statistics, the per-MAG table, quality distribution, abundance across samples, the
functional wheels from 2b, and AMR, virus and strain summaries for modules that ran. It also lists
the modules that were not run, so a missing section is not mistaken for a negative result. Built
with matplotlib's `PdfPages`, already a dependency (no HTML toolchain, browser or weasyprint), so it
builds on a cluster with no network.

## 3. MIMAG quality tags `[low]`: done

High, medium or low per Bowers et al. (*Nat Biotechnol* 2017), in `.stats.tsv` and the per-MAG
table. Besides completeness and contamination, MIMAG uses rRNA and tRNA counts, which are parsed
from the tRNAscan-SE, Aragorn and Infernal output of `annotate --bakta`.

MIMAG high quality requires the 5S, 16S and 23S rRNAs plus at least 18 tRNAs. A genome annotated
without Bakta has no rRNA evidence, so it is tagged `high_quality_unverified_rrna` with
`mimag_rrna=unknown` instead of medium.

What shipped:

- Three columns appended to `.stats.tsv`: `mimag_tag`, `mimag_rrna`, `mimag_trnas`.
- A one-line tier summary from `bin_refinement`.
- Five tags: `high_quality`, `high_quality_unverified_rrna`, `medium_quality`, `low_quality`, and
  `not_mimag` for bins above 10% contamination, which fall outside every MIMAG tier.
- Evidence is read from Bakta's TSV or GFF, with Infernal `--tblout` and tRNAscan-SE output as
  fallback when `annotate` wrote them separately.
- tRNAs are counted as distinct amino acids, not as features.
- The thresholds are Bowers' inequalities (>90, <5, >=50). They differ from the existing
  `quality_counts` tiers at 90.00 and 5.00; tests cover both.

## 3b. Bracken after Kraken2 `[low]`: done

Requested by users. Kraken2 reports reads assigned at each node, and reads matching conserved
sequence are assigned at genus or family, so a species at 10% of a community can have only a small
fraction of reads at its own node. Bracken (Lu et al., *PeerJ Comput Sci* 2017) redistributes those
reads to a chosen rank using the k-mer distribution of the Kraken2 database.

Shipped as `kraken2 --bracken`, with `--bracken-rank` (species by default) and
`--bracken-read-length`. Bracken needs Kraken2's `--report`, which the `kraken2` module did not
write before, so `--bracken` now adds it. Kraken2's own output is byte-identical without
`--bracken`.

Two conditions are checked at run time: the Bracken database must be built from the same Kraken2
database, and its read length must match the reads. Otherwise Bracken produces a table that
redistributes almost nothing, so the run reports the observed values next to the database's.
Output is the usual long-format abundance table, which `diff_abundance` reads directly.

## 4. Parser regression tests `[low]`

2.1.2 reduced the suite from ~1,360 tests to 181. This plan adds about forty back, all for parsers
of tables written by other tools. Parser bugs from 2.1.2's defect list:

- GUNC writes `pass.GUNC` as `1.0`/`0.0`, which scored a 166-bin set as 0% non-chimeric.
- inStrain stores identities as fractions where its paper says percentages.
- SingleM's read fraction is already a percentage.
- The depth table interleaves `-var` columns.

Each gave a wrong number without an error, and each had been covered by a fixture-based unit test
that was removed. An upstream release can bring any of them back.

Scope: one test per external table format, using fixtures already under `tests/fixtures/`. No new
frameworks, no machine-dependent assertions, no database.

## 5. eggNOG-mapper `[moderate]`, optional

In 2.1.2, of 3,891 genes in one bin, 1,578 had a Pfam domain and no KO, and 13 had a KO and no
Pfam. Pfam covers nearly all of KEGG's gene set but gives no GO terms, COGs or EC numbers.
eggNOG-mapper adds all three in one pass over `annotate`'s proteins, behind `--eggnog`, with its own
environment and `install-db` entry.

Lowest priority; moves to the next release if §1 and §2 fill this one.

## 6. Small items

| Item | Size | Why |
|---|---|---|
| GTDB-Tk as the `classify_bins` default | `[low]` | taxator-tk is deprecated with no upstream release; needs one measured arm |
| `assembly --keep-graph` | `[low]` | done. metaSPAdes' and metaFlye's graphs are hard-linked to a named path (121 MB and 2.3 MB measured, no extra disk) and MEGAHIT's FASTG is built from its intermediates (235 MB, 2 s). Needed for plasmid-aware and graph-aware binning. No graph was being deleted before; MEGAHIT never wrote one |
| Record the catalogue strategy in the cohort summary | `[low]` | done. `strategy`, `samples` and `dereplicated` in the study manifest, since co-assembly and per-sample catalogues are not comparable |
| Lockfiles for the eight new environments | `[low]` | amr, viruses, strains, comparative, singlem, rgi, tracs, benchmark-amber have none, so `--from-lock` falls back to solving from yaml without a warning |
| MultiQC across samples | `[low]` | one QC page per cohort instead of one per sample |
| MicrobeCensus genome equivalents | `[low]` | the denominator `quantify_function` needs for cross-sample comparison without MAGs |

## 7. Deferred

**Long reads and HiFi.** metaMDBG, polishing with racon and medaka, and hybrid assembly. The CAMI2
deposit has PacBio-simulated versions of all 49 samples (`*_pbsim/`), so a later long-read arm needs
only machine time. `--keep-graph` (§6) stays in 2.1.3 because plasmid-aware binning also needs it on
short-read assemblies.

**Per-sample versus co-assembly at twenty samples.** Not yet measured. It needs the same machine
as §1, which has priority. Next release.

**Workflow engines, cloud execution, metabolic models, protein and DNA language models, Hi-C
binning, HUMAnN 3, MOB-suite, Strainberry.** Out of scope. For metabolic models: 2.1.2 measured
that a genome at 70% completeness has lost 59% of its complete KEGG modules.

## Status

| Thread | Size | State |
|---|---|---|
| 1. Oral body site registered, run, figures extended | `[moderate]` | **done** - 51 clean HQ against metaWRAP 1's 45 (+13%); the gastrointestinal arms were re-scored with the same per-bin join (53 against 48) so all four numbers come from one method |
| 2a. Per-MAG summary table | `[low]` | **done** - `report --mag-table`, 15 columns joined from four inputs |
| 2b. Static hierarchy wheels | `[low]` | **done** - `metawrap2.sunburst`, one PNG per Krona wheel |
| 2c. PDF report | `[moderate]` | **done** - `metawrap2 report`, 23 pages on the reference run |
| 3. MIMAG tags | `[low]` | **done** - three columns, five tags, unverified rRNA tagged separately |
| 3b. Bracken after Kraken2 | `[low]` | **done** - `kraken2 --bracken`, four guards, no extra download |
| 4. Parser regression tests | `[low]` | **done** - 45 tests over 17 external formats; found a CAZyme EC defect |
| 5. eggNOG-mapper | `[moderate]` | **done** - `annotate --eggnog` in its own environment; 12,634 of 14,171 genes on five bins get a GO term, EC number or COG category |
| 6. Small items | `[low]` | `--keep-graph` and the catalogue-strategy record **done**; GTDB-Tk default still needs an arm; MultiQC, MicrobeCensus and the eight missing lockfiles not started |
| 7. Fresh-install fixes | `[moderate]` | **done** - twelve findings, four of them wrong-answer bugs; not in the original plan, and most of the release |
| 1b. Skin body site, and presets on oral | `[moderate]` | **done** - skin 50 clean HQ against metaWRAP 1's 41 (+22%); on oral `--preset slow` 59 against `default`'s 51 |

Thread 7 was unplanned. A tester installed from scratch on a bare machine. A conda environment
from another installation was reported as present but could not run, so `quickstart` built nothing
and `doctor` reported every module broken. Five listed databases could not be installed, CARD never
installed, and an interrupted `install-db` lost its config.

## Risks

**MetaWrap2 does worse on a body site.** The result goes in the same table.

**Disk.** A body site needs ~17 GB extracted plus benchmark output, against a 60% usage ceiling. A
recent cleanup freed 94.5 GB, enough for oral and skin but not for all five sites plus output.

**The two datasets disagree.** Then the documentation states which defaults depend on body site.

**PDF scope creep.** The report only assembles existing outputs and does no new computation.
