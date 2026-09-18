# Downloading and configuring MetaWrap2 databases

Databases are only needed for the modules that use them. Configure their locations in your MetaWrap2
config file (copy `metawrap2.toml.example` to `~/.metawrap2/config.toml`, or pass one with `--config`).
There is no more `config-metawrap` file - all database paths live in the `[databases]` section of
`metawrap2.toml`.

| Database | Size | Config key | Used by module(s) |
|---|---|---|---|
| CheckM DB | 1.4 GB | (set via `checkm data setRoot`) | binning (`--run-checkm`), bin_refinement, reassemble_bins |
| KRAKEN2 standard DB | ~125 GB | `KRAKEN2_DB` | kraken2 |
| NCBI_nt | ~71 GB | `BLASTDB` | blobology, classify_bins |
| NCBI_tax | ~283 MB | `TAXDUMP` | blobology, classify_bins |
| Indexed host genome (e.g. hg38) | ~20 GB | `BMTAGGER_DB` | read_qc (host removal) |

After downloading a database, set the matching key in `~/.metawrap2/config.toml`, for example:
```toml
[databases]
KRAKEN2_DB = "/path/to/my/kraken2_database"
BMTAGGER_DB = "/path/to/my/bmtagger_database"
BLASTDB = "/path/to/my/NCBI_nt"
TAXDUMP = "/path/to/my/NCBI_tax"
```
Confirm a module can find its databases with `metawrap2 check <module>`.

## CheckM database
Used by `bin_refinement`, `reassemble_bins`, and `binning --run-checkm`.
```bash
mkdir MY_CHECKM_FOLDER
cd MY_CHECKM_FOLDER
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
cd ../

# Tell CheckM where to find this data before running anything:
checkm data setRoot /path/to/your/dir/MY_CHECKM_FOLDER
```
CheckM stores this path itself, so there is no CheckM key in `metawrap2.toml`. (Run `checkm data setRoot`
from inside the `metawrap2-bin_refinement` environment so the path is set for the CheckM that the modules
use.)

## KRAKEN2 standard database
Used by the `kraken2` module. Compared to the old KRAKEN1 database, the KRAKEN2 database is considerably
more compact, so the download and indexing are much faster and less taxing on the system. You will need an
estimated 120 GB of RAM and ~128 GB of space. This is only needed if you intend to run the `kraken2` module.
```bash
kraken2-build --standard --threads 24 --db MY_KRAKEN2_DB
```
Then set `KRAKEN2_DB` in `~/.metawrap2/config.toml`:
```toml
KRAKEN2_DB = "/path/to/my/database/MY_KRAKEN2_DB"
```

## NCBI_nt BLAST database
Used by the `blobology` and `classify_bins` modules.
```bash
mkdir NCBI_nt
cd NCBI_nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
for a in nt.*.tar.gz; do tar xzf $a; done
```
Then set `BLASTDB` in `~/.metawrap2/config.toml`:
```toml
BLASTDB = "/your/location/of/database/NCBI_nt"
```

## NCBI taxonomy
Used by the `blobology` and `classify_bins` modules.
```bash
mkdir NCBI_tax
cd NCBI_tax
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```
Then set `TAXDUMP` in `~/.metawrap2/config.toml`:
```toml
TAXDUMP = "/your/location/of/database/NCBI_tax"
```

## Host genome index for bmtagger
Used by the `read_qc` module to remove host (e.g. human) reads. See the official bmtagger manual for
detailed instructions: https://www.hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf

First, download and merge the human genome hg38:
```bash
mkdir BMTAGGER_INDEX
cd BMTAGGER_INDEX
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
gunzip *fa.gz
cat *fa > hg38.fa
rm chr*.fa
```
Now index the human genome. Note that the index file names must be exactly as specified for MetaWrap2 to
recognize them! Indexing takes considerable memory and time (here we pass 100 GB of RAM as the `-M`
parameter):
```bash
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
```
MetaWrap2 looks for `hg38.bitmask` and `hg38.srprism` - make sure they are named exactly like this. Then
set `BMTAGGER_DB` in `~/.metawrap2/config.toml`:
```toml
BMTAGGER_DB = "/path/to/your/index/BMTAGGER_INDEX"
```

### Instructions for non-human hosts
For non-human hosts, the protocol for building and using the bmtagger index is the same. In the end you
need the `.srprism` and `.bitmask` index files in the directory you point `BMTAGGER_DB` at. When you run
`metawrap2 read_qc`, use the `-x/--host` option to specify the prefix of your host. For example, if your
`BMTAGGER_INDEX` directory has files `pig.srprism` and `pig.bitmask`, use `-x pig`. The default host
prefix is `hg38`.
