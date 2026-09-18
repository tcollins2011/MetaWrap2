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
| `--resume`  | reuse an existing output directory, skipping finished steps |

## Tinkering with the commands

Open a module in `src/metawrap2/modules/<module>.py`: the `COMMANDS` block holds the tool
command templates and the `CONSTANTS` block the tunable numbers. Edit them to change what
runs, or override any surfaced value in `~/.metawrap2/config.toml` without touching the source.

## Utilities

| Command | Purpose |
|---|---|
| `metawrap2 check <module>` | verify tools/databases/conda env for a module |
| `metawrap2 install-env <module\|--all>` | create (and verify) conda env(s) |
| `metawrap2 test [module]` | probe envs + software, run unit tests, print a status map |
| `metawrap2 config init\|show` | scaffold or inspect the config |
| `metawrap2 completion bash\|zsh` | print a shell-completion script |
