# Installation

MetaWrap2 has a light Python core plus one small conda environment per module.

```bash
git clone https://github.com/ursky/MetaWrap2.git
cd MetaWrap2
pip install -e .
```

Create the environment(s) for the module(s) you'll use. Environments are built and driven
with **mamba** (much faster than conda); `install-env` bootstraps mamba if it isn't present
and verifies each environment after building it:

```bash
metawrap2 install-env binning      # one module
metawrap2 install-env --all        # everything
```

Each environment is named `metawrap2-<module>` and is built from `envs/<module>.yaml`.

Check that a module is ready before running:

```bash
metawrap2 test binning             # env + software probe + unit tests, with a status map
```

## Configuration

```bash
metawrap2 config init              # writes ~/.metawrap2/config.toml
metawrap2 config show              # prints the resolved settings
```

Set database paths and settings in `metawrap2.toml`. To manage tools yourself instead of per
-module conda envs, set `use_conda_envs = false` and MetaWrap2 will call tools off your PATH.
