# Usage

```bash
metawrap2 -h                       # list modules and utilities
metawrap2 <module> -h              # per-module options
```

Each module is run separately and writes to a single `-o` output directory. Example:

```bash
metawrap2 binning -a assembly.fa -o binning_out -t 24 --metabat2 --maxbin2 --concoct \
    sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq
```

Gzipped inputs (`.gz`) are accepted throughout. Every module is also runnable directly, e.g.
`python -m metawrap2.modules.binning ...`.

## Global options (before the module name)

| Option | Effect |
|---|---|
| `--dry-run` | print/record the commands without running them |
| `--force`   | overwrite an existing MetaWrap2 output directory |
| `--resume`  | reuse an existing output directory, re-running anything it cannot verify |
| `--verbose-logs` | don't collapse repeated tool messages in the captured logs |
| `--skip-space-check` | start even if the disk-space estimate says there is no room |

## Tinkering with the commands

Open a module in `src/metawrap2/modules/<module>.py`: the `COMMANDS` block holds the tool
command templates and the `CONSTANTS` block the tunable numbers. Edit them to change what
runs, or override any surfaced value in `~/.metawrap2/config.toml` without touching the source.

## Utilities

| Command | Purpose |
|---|---|
| `metawrap2 run <sheet\|dir> -o study/` | run the whole pipeline for a study from one sample sheet |
| `metawrap2 check <module>` | verify tools/databases/conda env for a module |
| `metawrap2 install-env <module\|--all> [-t N]` | create (and verify) conda env(s), N at a time |
| `metawrap2 install-db <name\|--all> [-t N] [--small]` | download, index, and configure databases |
| `metawrap2 doctor [module]` | probe envs, software and databases; print a status map (alias: `test`) |
| `metawrap2 doctor --fix` | recreate exactly the broken environments, leaving healthy ones alone |
| `metawrap2 status <study>` | what has run, and what `--resume` would do next |
| `metawrap2 history` | every run on this machine, and when |
| `metawrap2 calibrate <study>...` | check the built-in disk/memory estimates against real runs |
| `metawrap2 config init\|show` | scaffold or inspect the config |
| `metawrap2 completion bash\|zsh` | print a shell-completion script |

### Checking on a run

`metawrap2 status` reads the study's run manifest, so it can say what a resume would do without
doing it - the `--resume` column comes from the same check `--resume` itself makes:

```
metawrap2 status study/
metawrap2 status study/ --step binning     # inputs, outputs, command, and measured cost
metawrap2 status study/ --audit            # is every recorded output still there and non-empty?
```

The same audit runs automatically at the end of every run, so a step that reported success while
producing nothing is reported then rather than silently mis-skipped later.

### Comparing two runs

```
metawrap2 history                    # the 20 most recent runs
metawrap2 history --failed           # what didn't finish
metawrap2 history --diff RUN RUN     # version, env, databases and command line, side by side
```

Run ids come from `metawrap2 history --json`; a unique prefix is enough.

## Bin file names

Bins are named `bin_001.fasta`, for every binner, counted from 1 and zero-padded. The number is a
label within one binner's output, not an identity shared across binners - `metabat2_bins/bin_001`
and `maxbin2_bins/bin_001` are unrelated. If you pass bins to CheckM or GTDB-Tk yourself, the
extension flag is `-x fasta`.
