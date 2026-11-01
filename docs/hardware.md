# Hardware

8 cores and 64 GB of RAM is a reasonable minimum; requirements scale with the data. The figures
below come from the runs' usage records. Steps are reaped with `os.wait4`, which reports their peak
memory, and a step that uses a GPU has its peak VRAM sampled while it runs.

## Memory

Peak RAM per step on the three CAMI2 skin samples (`--preset default`):

| Step | Peak RAM | Driven by |
|:---|---:|:---|
| `read_qc` | 1.3 GB | streaming, so read volume barely matters |
| `assembly` (MEGAHIT) | 4.1-5.5 GB | the k-mer graph |
| `binning` (metaBAT2 + SemiBin2) | 18.4-21.7 GB | depth tables and the learned binner's model |
| `bin_refinement` (CheckM2 + GUNC) | 34.0-34.1 GB | the two reference databases |

These settings affect memory more than the sample does:

- **Assembler.** On three co-assembled CAMI2 gastrointestinal samples metaSPAdes peaked at 140 GB
  and MEGAHIT at 15 GB; metaSPAdes gave a 38% longer N50. `assembly -m` sets the ceiling.
- **Co-assembly.** It holds every sample's reads at once, so memory grows with the study.
  Per-sample assembly with `cohort` avoids this.
- **GTDB-Tk.** `classify_bins --gtdbtk` needs about 140 GB for `classify_wf`, by GTDB-Tk's own
  documentation.

## GPUs

GPUs are optional. Only the three learned binners can use one (SemiBin2, VAMB and COMEBin, in the
`learned` environment). They default to the CPU (`--device cpu`) so that results are the same on any
machine, except under `--preset deep`, which trains SemiBin2 and COMEBin with `--device cuda`.
`--device cuda` or `--device auto` turns the GPU on for any binner, and an explicit `--device`
overrides the preset.

Peak VRAM on a single device, sampled every 5 s on NVIDIA GPUs:

| Binner | Peak VRAM | Measured on | Grows with the assembly |
|:---|---:|:---|:---|
| VAMB | 0.56 GB | MetaHIT, 16,527 contigs | not measured |
| SemiBin2 | 0.41-2.75 GB | CAMI2 skin, 16,936 and 63,376 contigs | yes, faster than linearly |
| COMEBin | 5.71 GB | the same two skin samples | no, set by batch size |

An 8 GB card runs all three on this data. For assemblies larger than these, allow more than
2.75 GB for SemiBin2. COMEBin, which `--preset deep` adds, needs the most VRAM and gains the most
from a GPU. On CPU it ran for over 280 CPU-hours on three CAMI2 samples without finishing; metaBAT2
and MaxBin2 together took 7.1 CPU-hours.

The per-step cost line and `.metawrap2/run_steps/manifest.json` of each run give peak memory, CPU
time and, for a GPU step, `peak_gpu_memory_bytes`.

## Refinement

`bin_refinement` scores its hybrid bin sets with CheckM2 several at a time when memory allows
(`--checkm-jobs`): about 24 GB per CheckM2 call, or with `--checkm1` 16 GB per call on CheckM1's
reduced tree and 40 GB with `--full-tree`. Timings are in [Benchmarks](benchmarks.md#26-refinement-speed).
