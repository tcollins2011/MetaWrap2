"""`metawrap2 install-db [name ...] [--all]` - download and configure MetaWrap2's databases.

Database setup used to be the roughest part of getting MetaWrap2 working: five separate
manual downloads copy-pasted out of the docs, each with its own unpacking and indexing step,
followed by hand-editing ``metawrap2.toml`` with the right key names, plus a
``checkm data setRoot`` that is easy to forget and silently breaks bin_refinement later.

This command does all of that: it fetches each database, unpacks it, runs whatever indexing
step it needs, writes the matching key into the config file, and verifies the result. Every
step is skipped if it already looks complete, so re-running it is cheap and safe.

    metawrap2 install-db --list            # what is available, how big, who needs it
    metawrap2 install-db checkm taxdump    # just these two
    metawrap2 install-db --all             # everything (large: see --list for sizes)

Sizes are the real thing, so ``--all`` is a many-hour, several-hundred-GB commitment. The
per-database ``--small`` variants download a capped-but-real stand-in instead, which is what
you want for a smoke test of the pipeline rather than a production analysis.
"""

from __future__ import annotations

import concurrent.futures
import glob
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import threading
from typing import Callable, Dict, List, NamedTuple, Optional

from ..cliformat import make_parser
from ..config import MODULE_ENVS, OPTIONAL_ENVS, find_config

# Databases are independent of each other and mostly network-bound, so fetching several at
# once is close to a linear speedup. Modest default so a slow link is not swamped.
DEFAULT_PARALLEL_DBS = 3
# Serialises the progress prints from concurrent installers so lines do not interleave.
_PRINT_LOCK = threading.Lock()


def _say(*lines: str) -> None:
    with _PRINT_LOCK:
        for line in lines:
            print(line)
        sys.stdout.flush()


class Database(NamedTuple):
    key: Optional[str]  # metawrap2.toml [databases] key, None if the tool stores its own path
    used_by: str
    size: str
    small_size: Optional[str]
    describe: str
    #: Written into the config instead of the install directory, when the tool wants a file rather
    #: than a directory - GUNC takes ``--db_file``, so GUNC_DB has to name the .dmnd itself.
    key_file: Optional[str] = None
    #: Like ``key_file``, for a path whose name is not known until it is downloaded. SingleM's
    #: metapackage directory carries its own version (``S6.5.0.GTDB_r232.metapackage_...``), and
    #: that version is part of the answer - a per-cell fraction is only as current as the reference
    #: behind it - so the versioned directory is what goes into the config rather than a `latest`
    #: symlink that would hide it.
    key_glob: Optional[str] = None


DATABASES: Dict[str, Database] = {
    "checkm": Database(
        key="CHECKM_DB",
        used_by="bin_refinement, reassemble_bins, binning --run-checkm",
        size="1.4 GB",
        small_size=None,
        describe="CheckM1 marker-gene reference data",
    ),
    "taxdump": Database(
        key="TAXDUMP",
        used_by="blobology, classify_bins",
        size="550 MB",
        small_size=None,
        describe="NCBI taxonomy dump (names.dmp / nodes.dmp)",
    ),
    "kraken2": Database(
        key="KRAKEN2_DB",
        used_by="kraken2",
        size="~90 GB",
        small_size="8 GB",
        describe="Kraken2 standard database (--small: the capped 8 GB standard build)",
    ),
    "blast": Database(
        key="BLASTDB",
        used_by="blobology, classify_bins",
        size="~300 GB",
        small_size="~2 GB",
        describe="NCBI nt (--small: complete representative prokaryote genomes)",
    ),
    "checkm2": Database(
        key="CHECKM2DB",
        used_by="bin_refinement --checkm2, binning --run-checkm --checkm2",
        size="1.7 GB download",
        small_size=None,
        describe="CheckM2's DIAMOND reference (uniref100.KO, ~3 GB unpacked; not CheckM1's data)",
        key_file=os.path.join("CheckM2_database", "uniref100.KO.1.dmnd"),
    ),
    "gunc": Database(
        key="GUNC_DB",
        used_by="bin_refinement (unless --skip-gunc)",
        size="6.7 GB download",
        small_size=None,
        describe="GUNC proGenomes 2.1 reference (diamond .dmnd; unpacks larger than it downloads)",
        key_file="gunc_db_progenomes2.1.dmnd",
    ),
    "dbcan": Database(
        key="DBCAN_DB",
        used_by="annotate --cazy",
        size="7.4 GB",
        small_size=None,
        describe="dbCAN CAZyme HMMs, CAZy DIAMOND database and substrate mapping",
    ),
    "pfam": Database(
        key="PFAM_DB",
        used_by="annotate --pfam",
        size="399 MB download, 1.6 GB unpacked",
        small_size=None,
        describe="Pfam-A HMMs with their per-family gathering thresholds, plus the release notes",
        key_file="Pfam-A.hmm",
    ),
    "card": Database(
        key="CARD_DB",
        used_by="amr --rgi",
        size="4 MB download",
        small_size=None,
        describe="CARD's resistance ontology and reference sequences (card.json plus its FASTAs)",
        key_file="card.json",
    ),
    "singlem": Database(
        key="SINGLEM_METAPACKAGE",
        used_by="appraise",
        size="2.2 GB download, 9.7 GB unpacked",
        small_size=None,
        describe="SingleM GTDB marker metapackage (single-copy marker genes + their taxonomy)",
        key_glob="*.smpkg.zb",
    ),
    "kofam": Database(
        key="KOFAM_DB",
        used_by="annotate --kofam",
        size="1.5 GB download",
        small_size=None,
        describe="KOfam HMM profiles and their per-profile thresholds (ko_list)",
    ),
    "kegg": Database(
        key="KEGG_BRITE",
        used_by="pathways (functional hierarchy and Krona plots)",
        size="4 MB",
        small_size=None,
        describe="KEGG BRITE ko00001: KO -> pathway -> category, projected to one table",
    ),
    "host": Database(
        key="BMTAGGER_DB",
        used_by="read_qc host removal",
        size="~20 GB + indexing",
        small_size="~9 GB",
        describe="bmtagger index of the host genome (--small: human chr21 only)",
    ),
}


def _print_list() -> None:
    print("\nMetaWrap2 databases\n")
    print("  %-10s %-20s %-10s %s" % ("name", "full size", "--small", "used by"))
    print("  " + "-" * 84)
    for name, db in DATABASES.items():
        print("  %-10s %-20s %-10s %s" % (name, db.size, db.small_size or "-", db.used_by))
    print()
    for name, db in DATABASES.items():
        print("  %-10s %s" % (name, db.describe))
    print("\n  Install into a directory of your choice with --dir (default ~/metawrap2_dbs).")
    print("  Each database's path is written to your config file automatically.\n")


# --- shared helpers ---------------------------------------------------------------------


def _run(cmd: List[str], **kw) -> int:
    print("  $ %s" % " ".join(cmd))
    return subprocess.call(cmd, **kw)


def _fetch(url: str, dest: str) -> bool:
    """Download *url* to *dest* with curl, resuming a partial file."""
    if os.path.isfile(dest) and os.path.getsize(dest):
        print("  already downloaded: %s" % os.path.basename(dest))
        return True
    tmp = dest + ".part"
    rc = _run(["curl", "-fL", "--retry", "3", "-C", "-", "-o", tmp, url])
    if rc != 0:
        print("  download failed: %s" % url)
        return False
    os.replace(tmp, dest)
    return True


def _untar(archive: str, into: str) -> None:
    """Unpack *archive* into *into*, refusing members that would land outside it.

    These archives come off the public internet, and a bare `extractall` writes wherever a member's
    name says to - including `/etc/...` and `../../..`. Python 3.12's `data` filter rejects exactly
    that (and strips the ownership and mode bits a database tarball has no business setting); the
    explicit check below is the same guard for older versions, so the refusal does not depend on
    which interpreter happens to be running.
    """
    print("  unpacking %s" % os.path.basename(archive))
    with tarfile.open(archive) as tf:
        safe = []
        for member in tf.getmembers():
            if _escapes(member):
                print("  skipping %s: it would be written outside %s" % (member.name, into))
                continue
            if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                # A device node or fifo in a database tarball is not something to write out.
                print("  skipping %s: not a file, directory or link" % member.name)
                continue
            safe.append(member)
        if hasattr(tarfile, "data_filter"):
            tf.extractall(into, members=safe, filter="data")
        else:
            tf.extractall(into, members=safe)
    os.remove(archive)


def _gunzip(archive: str, dest: str) -> bool:
    """Decompress *archive* to *dest*, removing the archive. Streamed, not read into memory."""
    print("  unpacking %s" % os.path.basename(archive))
    tmp = dest + ".part"
    with gzip.open(archive, "rb") as src, open(tmp, "wb") as out:
        shutil.copyfileobj(src, out, length=8 << 20)
    os.replace(tmp, dest)
    os.remove(archive)
    return os.path.isfile(dest) and os.path.getsize(dest) > 0


def _escapes(member: tarfile.TarInfo) -> bool:
    """Would extracting *member* write outside the destination directory?

    Checked on the link target as well as the name: a symlink to `/etc` is how the next member in
    the archive gets written there.
    """
    for name in (member.name, member.linkname or ""):
        if name.startswith(("/", "\\")):
            return True
        if ".." in name.replace("\\", "/").split("/"):
            return True
    return False


def _mamba_run(env: str, cmd: List[str]) -> int:
    # Flags come from _run_flags: --no-capture-output is conda-only and mamba 2.x chokes on it.
    from ..command import _run_flags

    return _run(["mamba", "run"] + _run_flags("mamba") + ["-n", env] + cmd)


# --- per-database installers -------------------------------------------------------------


def _install_checkm(target: str, small: bool, jobs: int = 8) -> bool:
    if os.path.isfile(os.path.join(target, "taxon_marker_sets.tsv")):
        print("  CheckM data already present.")
    else:
        os.makedirs(target, exist_ok=True)
        archive = os.path.join(target, "checkm_data.tar.gz")
        if not _fetch(
            "https://data.ace.uq.edu.au/public/CheckM_databases/" "checkm_data_2015_01_16.tar.gz",
            archive,
        ):
            return False
        _untar(archive, target)

    # CheckM stores this path inside its own env, so this is belt and braces rather than the
    # mechanism: every CheckM1 call also gets CHECKM_DATA_PATH from the configured CHECKM_DB (see
    # checkm.checkm1_data_env), which is what makes the setting survive the env being rebuilt. Doing it
    # anyway means `checkm` run by hand from the env also works.
    #
    # The env name comes from MODULE_ENVS rather than being written out: it used to be spelled
    # "metawrap2-bin_refinement" here, and consolidating the envs left this line pointing at one that
    # no longer exists - which fails quietly, since the message below is only a NOTE.
    checkm_env = MODULE_ENVS["bin_refinement"]
    print("  pointing CheckM at the data (checkm data setRoot)")
    if _mamba_run(checkm_env, ["checkm", "data", "setRoot", target]) != 0:
        print(
            "  NOTE: could not run `checkm data setRoot` automatically. CHECKM_DATA_PATH is passed\n"
            "        per call so bin_refinement will work regardless, but to use `checkm` directly:\n"
            "        mamba run -n %s checkm data setRoot %s" % (checkm_env, target)
        )
    return os.path.isfile(os.path.join(target, "taxon_marker_sets.tsv"))


def _install_taxdump(target: str, small: bool, jobs: int = 8) -> bool:
    if os.path.isfile(os.path.join(target, "names.dmp")):
        print("  taxdump already present.")
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "taxdump.tar.gz")
    if not _fetch("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", archive):
        return False
    _untar(archive, target)
    return os.path.isfile(os.path.join(target, "names.dmp"))


def _install_kraken2(target: str, small: bool, jobs: int = 8) -> bool:
    if os.path.isfile(os.path.join(target, "hash.k2d")):
        print("  Kraken2 database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if small:
        # A real standard build, capped to 8 GB - right for testing, not for production.
        archive = os.path.join(target, "k2_standard_08gb.tar.gz")
        if not _fetch(
            "https://genome-idx.s3.amazonaws.com/kraken/" "k2_standard_08gb_20240904.tar.gz",
            archive,
        ):
            return False
        _untar(archive, target)
    else:
        archive = os.path.join(target, "k2_standard.tar.gz")
        if not _fetch(
            "https://genome-idx.s3.amazonaws.com/kraken/" "k2_standard_20240904.tar.gz", archive
        ):
            return False
        _untar(archive, target)
    return os.path.isfile(os.path.join(target, "hash.k2d"))


def _install_blast(target: str, small: bool, jobs: int = 8) -> bool:
    from .. import blastdb
    from ..config import Settings

    os.makedirs(target, exist_ok=True)
    # Every other installer short-circuits when its database is already there; without the
    # same check here, re-running install-db re-ran makeblastdb over the whole database.
    present, _ = blastdb.find(Settings(databases={"BLASTDB": target}))
    if present:
        print("  BLAST database already present.")
        return True
    if small:
        print("  building a representative-genome BLAST database (stand-in for nt)")
        script = os.path.join(target, "build_small_nt.sh")
        with open(script, "w") as fh:
            fh.write(_SMALL_BLAST_SCRIPT)
        os.chmod(script, 0o755)
        # Thousands of small genome downloads: run them concurrently or this dominates.
        if _run(["bash", script, str(jobs)], cwd=target) != 0:
            return False
    else:
        print("  downloading NCBI nt with update_blastdb.pl (this is very large)")
        if (
            _mamba_run(MODULE_ENVS["classify_bins"], ["update_blastdb.pl", "--decompress", "nt"])
            != 0
        ):
            return False
    return True


_SMALL_BLAST_SCRIPT = r"""#!/bin/bash
# Build a small but real BLAST nucleotide database from complete RefSeq representative
# prokaryote genomes, named "nt" so it matches MetaWrap2's default BLASTDB_NAME.
set -euo pipefail
JOBS="${1:-8}"
curl -fsSL -o assembly_summary.txt \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F'\t' '$12=="Complete Genome" && $11=="latest" && ($5=="reference genome" || $5=="representative genome")' \
  assembly_summary.txt | awk -F'\t' '{print $1"\t"$6"\t"$20}' > download_list.tsv
mkdir -p genomes

# Fetch the genomes JOBS at a time - serially this is the slowest part of the whole build.
fetch_one() {
  url="$1"; base=$(basename "$url"); f="genomes/${base}_genomic.fna.gz"
  [ -s "$f" ] && return 0
  curl -fsSL --retry 2 -o "$f.part" "${url}/${base}_genomic.fna.gz" && mv "$f.part" "$f" \
    || { rm -f "$f.part"; return 0; }
}
export -f fetch_one
cut -f3 download_list.tsv | xargs -P "$JOBS" -I{} bash -c 'fetch_one "$@"' _ {}

# Build the seqid -> taxid map after the downloads, so a retry never appends twice.
: > taxid_map.txt
while IFS=$'\t' read -r acc taxid url; do
  base=$(basename "$url"); f="genomes/${base}_genomic.fna.gz"
  [ -s "$f" ] || continue
  zcat "$f" | grep '^>' | sed 's/^>//' | awk -v t="$taxid" '{print $1"\t"t}' >> taxid_map.txt
done < download_list.tsv
zcat genomes/*.fna.gz > all_genomes.fna
makeblastdb -in all_genomes.fna -dbtype nucl -out nt -parse_seqids -taxid_map taxid_map.txt \
  -title "RefSeq representative prokaryote genomes"
rm -f all_genomes.fna
"""


def _install_checkm2(target: str, small: bool, jobs: int = 8) -> bool:
    """CheckM2's DIAMOND reference, fetched by CheckM2 itself.

    ``checkm2 database --download`` is the only supported way to get this: the file lives on Zenodo
    under a versioned record that CheckM2 pins, and a hand-picked URL goes stale the next time
    CheckM2 is upgraded. So this installer runs the tool rather than fetching a URL of its own.

    It is a different database from ``checkm``, which is the one thing about CheckM2 that reliably
    catches people out - both are "the CheckM database", they share no files, and having one
    installed does nothing for the other. ``--small`` is ignored: completeness estimates come from a
    model trained against this exact reference, so a trimmed copy would not be a faster CheckM2 but
    a differently-calibrated one.
    """
    dmnd = os.path.join(target, DATABASES["checkm2"].key_file or "")
    if os.path.isfile(dmnd) and os.path.getsize(dmnd):
        print("  CheckM2 database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if (
        _mamba_run(
            OPTIONAL_ENVS["checkm2"],
            ["checkm2", "database", "--download", "--path", target],
        )
        != 0
    ):
        return False
    return os.path.isfile(dmnd) and os.path.getsize(dmnd) > 0


def _install_gunc(target: str, small: bool, jobs: int = 8) -> bool:
    """The GUNC proGenomes 2.1 diamond database.

    One file, so there is no indexing step and nothing to configure inside a conda env - which makes
    this the simplest installer here. ``--small`` is ignored because GUNC's verdict is only as good
    as the breadth of what it compares against: a capped reference would report bins as
    under-represented rather than as clean, which is a worse answer than a slower one.
    """
    dmnd = os.path.join(target, DATABASES["gunc"].key_file or "")
    if os.path.isfile(dmnd) and os.path.getsize(dmnd):
        print("  GUNC database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    archive = dmnd + ".gz"
    if not _fetch("https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz", archive):
        return False
    return _gunzip(archive, dmnd)


def _install_dbcan(target: str, small: bool, jobs: int = 8) -> bool:
    """dbCAN's databases, fetched by dbCAN itself.

    ``--aws_s3`` rather than the default HTTP path, and that is the whole decision here: dbCAN's
    own help calls the HTTP route "db_current (moving snapshot)" while the S3 route is a pinned
    release. A CAZyme family call is only reproducible against the release that produced it, so a
    database that changes under a stable name is worse than a large one.

    ``--small`` is ignored: the HMM sets are the method, and the three-method agreement this module
    reports needs all three present.
    """
    marker = os.path.join(target, "dbCAN.hmm")
    if os.path.isfile(marker) and os.path.getsize(marker):
        print("  dbCAN database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if (
        _mamba_run(
            MODULE_ENVS["annotate"],
            ["run_dbcan", "database", "--db_dir", target, "--aws_s3"],
        )
        != 0
    ):
        return False
    return os.path.isfile(marker) and os.path.getsize(marker) > 0


def _install_pfam(target: str, small: bool, jobs: int = 8) -> bool:
    """Pfam-A, from EBI's current release.

    The HMMs carry their own gathering thresholds, which is why no separate cutoff table is needed -
    and why this is one file rather than the half-dozen InterProScan's data would be.

    `relnotes.txt` is fetched with them and is not decoration: Pfam renumbers and merges families
    between releases, so a domain call means "PF00001 as of release 37.0" and the release is what
    makes it reproducible. `--small` is ignored; the profile set is the method.
    """
    marker = os.path.join(target, "Pfam-A.hmm")
    if os.path.isfile(marker) and os.path.getsize(marker):
        print("  Pfam already present (release %s)." % _pfam_release(target))
        return True
    os.makedirs(target, exist_ok=True)
    base = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
    for name in ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz"):
        if not _fetch("%s/%s" % (base, name), os.path.join(target, name)):
            return False
    _fetch("%s/relnotes.txt" % base, os.path.join(target, "relnotes.txt"))
    for name in ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz"):
        if _run(["gunzip", "-f", os.path.join(target, name)]) != 0:
            return False
    if not (os.path.isfile(marker) and os.path.getsize(marker)):
        return False
    print("  Pfam release %s" % _pfam_release(target))
    return True


def _pfam_release(target: str) -> str:
    from ..pfam import release_from_notes

    return release_from_notes(target)


def _install_card(target: str, small: bool, jobs: int = 8) -> bool:
    """CARD's data, from CARD's own `latest/data` endpoint.

    Fetched by URL rather than by a tool, because RGI has no downloader: `rgi load` takes a
    card.json that is already on disk. The tarball is 4 MB, which is why this one is listed with no
    `--small` variant - there is nothing to trim.

    The version inside card.json is what dates every call made against it, so it is printed here and
    recorded per run in the AMR table. CARD releases every few months and the URL always points at
    the newest one, so two installs months apart are different databases under one path: rename the
    directory after the version if that matters to you.
    """
    marker = os.path.join(target, "card.json")
    if os.path.isfile(marker) and os.path.getsize(marker):
        print("  CARD database already present (version %s)." % _card_version(marker))
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "card_data.tar.bz2")
    if not _fetch("https://card.mcmaster.ca/latest/data", archive):
        return False
    import tarfile

    try:
        with tarfile.open(archive, "r:bz2") as tar:
            # No member path escapes the target: a tarball is untrusted input, and `extractall`
            # without this will happily write to ../ if a member says so.
            for member in tar.getmembers():
                destination = os.path.realpath(os.path.join(target, member.name))
                if not destination.startswith(os.path.realpath(target) + os.sep):
                    print("  refusing to extract %r: it escapes %s" % (member.name, target))
                    return False
            tar.extractall(target)
    except (OSError, tarfile.TarError) as problem:
        print("  could not unpack CARD: %s" % problem)
        return False
    if not (os.path.isfile(marker) and os.path.getsize(marker)):
        return False
    print("  CARD version %s" % _card_version(marker))
    return True


def _card_version(card_json: str) -> str:
    from ..amr import card_version

    return card_version(card_json)


def _install_singlem(target: str, small: bool, jobs: int = 8) -> bool:
    """The SingleM metapackage, fetched by SingleM itself.

    ``singlem data`` is the supported route for the same reason ``checkm2 database --download`` is:
    the package lives in a Zenodo record whose DOI SingleM pins per release, it is distributed as a
    ZenodoBackpack with checksums SingleM verifies after unpacking, and a hand-picked URL goes stale
    at the next upgrade.

    ``--small`` is ignored. The metapackage is the marker set and its taxonomy; a trimmed copy would
    not be a faster appraisal but a different reference, and the whole figure is relative to that
    reference.
    """
    existing = glob.glob(os.path.join(target, "*.smpkg.zb"))
    if existing:
        print("  SingleM metapackage already present: %s" % os.path.basename(existing[0]))
        return True
    os.makedirs(target, exist_ok=True)
    if (
        _mamba_run(
            MODULE_ENVS["appraise"],
            ["singlem", "data", "--output-directory", target],
        )
        != 0
    ):
        return False
    return bool(glob.glob(os.path.join(target, "*.smpkg.zb")))


def _install_host(target: str, small: bool, jobs: int = 8) -> bool:
    fasta = os.path.join(target, "hg38.fa")
    os.makedirs(target, exist_ok=True)
    if not os.path.isfile(fasta):
        if small:
            # Real human sequence, one chromosome - enough to exercise host removal.
            gz = os.path.join(target, "chr21.fa.gz")
            if not _fetch(
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/" "chr21.fa.gz", gz
            ):
                return False
            print("  decompressing to hg38.fa")
            import gzip

            with gzip.open(gz, "rb") as fin, open(fasta, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.remove(gz)
        else:
            gz = os.path.join(target, "hg38.fa.gz")
            if not _fetch("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz", gz):
                return False
            import gzip

            with gzip.open(gz, "rb") as fin, open(fasta, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.remove(gz)

    bitmask = os.path.join(target, "hg38.bitmask")
    if not os.path.isfile(bitmask):
        print("  building the bmtagger bitmask (bmtool)")
        if (
            _mamba_run(
                MODULE_ENVS["read_qc"],
                ["bmtool", "-d", fasta, "-o", bitmask, "-A", "0", "-w", "18"],
            )
            != 0
        ):
            return False
    srprism = os.path.join(target, "hg38.srprism")
    if not os.path.isfile(srprism + ".idx"):
        print("  building the srprism index (this is the slow step)")
        if (
            _mamba_run(
                MODULE_ENVS["read_qc"],
                ["srprism", "mkindex", "-i", fasta, "-o", srprism, "-M", "60000"],
            )
            != 0
        ):
            return False
    return os.path.isfile(bitmask) and os.path.isfile(srprism + ".idx")


def _install_kofam(target: str, small: bool, jobs: int = 8) -> bool:
    """KOfam: the HMM profiles and ko_list, their per-profile thresholds, from KEGG's FTP.

    The two belong together: a profile's threshold is what makes a hit a KO call, so they are
    downloaded into one directory and `annotate` reads the thresholds from beside the profiles.
    `--small` is ignored; the profile set is the method.
    """
    from ..kofam import KO_LIST_MD5, KOFAM_URL, ko_list_md5

    ko_list = os.path.join(target, "ko_list.gz")
    if os.path.isfile(ko_list) and os.path.isdir(os.path.join(target, "profiles")):
        print("  KOfam already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if not _fetch(KOFAM_URL + "ko_list.gz", ko_list):
        return False
    if ko_list_md5(ko_list) != KO_LIST_MD5:
        print("  note: this ko_list differs from the one this release was validated on.")
    archive = os.path.join(target, "profiles.tar.gz")
    if not _fetch(KOFAM_URL + "profiles.tar.gz", archive):
        return False
    _untar(archive, target)
    return os.path.isdir(os.path.join(target, "profiles"))


def _install_kegg(target: str, small: bool, jobs: int = 8) -> bool:
    """KEGG BRITE ko00001, projected to the per-(KO, pathway) table `pathways` reads.

    The download is recorded in SNAPSHOT.json (date and md5), and every hierarchy row names it:
    BRITE is revised, so a percentage is only comparable with one over the same hierarchy.
    """
    from ..hierarchy import BRITE_URL, HIERARCHY_FILE, project_htext, write_snapshot

    if os.path.isfile(os.path.join(target, HIERARCHY_FILE)):
        print("  KEGG hierarchy already present.")
        return True
    os.makedirs(target, exist_ok=True)
    htext = os.path.join(target, "ko00001.keg")
    if not _fetch(BRITE_URL, htext):
        return False
    rows = project_htext(htext, os.path.join(target, HIERARCHY_FILE))
    if rows == 0:
        print("  the download did not look like BRITE ko00001: no KO rows found")
        return False
    write_snapshot(target, htext, BRITE_URL, rows)
    os.remove(htext)
    print("  %d KO-pathway rows" % rows)
    return True


INSTALLERS: Dict[str, Callable[..., bool]] = {
    "checkm": _install_checkm,
    "taxdump": _install_taxdump,
    "kraken2": _install_kraken2,
    "blast": _install_blast,
    "checkm2": _install_checkm2,
    "gunc": _install_gunc,
    "card": _install_card,
    "pfam": _install_pfam,
    "dbcan": _install_dbcan,
    "singlem": _install_singlem,
    "kofam": _install_kofam,
    "kegg": _install_kegg,
    "host": _install_host,
}


# --- writing the config ------------------------------------------------------------------


def _config_path(explicit: Optional[str]) -> str:
    return explicit or find_config() or os.path.expanduser("~/.metawrap2/config.toml")


def _set_config_keys(path: str, values: Dict[str, str]) -> None:
    """Set ``key = "value"`` under [databases] in *path*, creating the file/section if needed.

    Deliberately a line edit rather than a parse-and-rewrite: TOML comments in a user's
    config are worth preserving, and the standard library has no TOML writer.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lines = []
    if os.path.isfile(path):
        with open(path) as fh:
            lines = fh.read().splitlines()

    if not any(ln.strip() == "[databases]" for ln in lines):
        if lines and lines[-1].strip():
            lines.append("")
        lines.append("[databases]")

    start = next(i for i, ln in enumerate(lines) if ln.strip() == "[databases]")
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].strip().startswith("["):
            end = i
            break

    for key, value in values.items():
        entry = '%s = "%s"' % (key, value)
        for i in range(start + 1, end):
            stripped = lines[i].lstrip()
            if stripped.split("=")[0].strip() == key and not stripped.startswith("#"):
                lines[i] = entry
                break
        else:
            lines.insert(end, entry)
            end += 1

    with open(path, "w") as fh:
        fh.write("\n".join(lines).rstrip("\n") + "\n")


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 install-db",
        description="Download, index, and configure the databases MetaWrap2 modules need.",
    )
    ap.add_argument("names", nargs="*", help="databases to install (see --list)")
    ap.add_argument("--all", action="store_true", help="install every database")
    ap.add_argument("--list", action="store_true", help="show available databases and sizes")
    ap.add_argument(
        "--dir",
        default=os.path.expanduser("~/metawrap2_dbs"),
        help="parent directory to install into (default ~/metawrap2_dbs)",
    )
    ap.add_argument(
        "--small",
        action="store_true",
        help="where a capped-but-real variant exists, install that instead "
        "(much faster; good for testing, not for production analysis)",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_PARALLEL_DBS,
        help="how many databases to fetch concurrently (default %d). Also sets "
        "the parallelism of the many small downloads inside a single "
        "database build." % DEFAULT_PARALLEL_DBS,
    )
    ap.add_argument("--config", help="config file to update (default: the one in effect)")
    ap.add_argument(
        "--no-config-update",
        action="store_true",
        help="don't write the resulting paths into the config file",
    )
    args = ap.parse_args(argv)

    if args.list or (not args.names and not args.all):
        _print_list()
        return 0 if args.list else 1

    names = list(DATABASES) if args.all else args.names
    unknown = [n for n in names if n not in DATABASES]
    if unknown:
        print("Unknown database(s): %s\nKnown: %s" % (", ".join(unknown), ", ".join(DATABASES)))
        return 1

    if not shutil.which("mamba") and any(
        n in ("checkm", "checkm2", "blast", "host") for n in names
    ):
        print(
            "mamba is not on PATH, and checkm/checkm2/blast/host need tools from the module envs.\n"
            "Install the envs first: metawrap2 install-env --all"
        )
        return 1

    installed: Dict[str, str] = {}
    failed: List[str] = []

    def install(name: str):
        db = DATABASES[name]
        target = os.path.join(args.dir, name.upper())
        size = db.small_size if (args.small and db.small_size) else db.size
        note = (
            "  (no smaller variant for this one; installing the full database)"
            if args.small and not db.small_size
            else ""
        )
        _say("=== %s (%s) -> %s" % (name, size, target), *([note] if note else []))
        try:
            ok = INSTALLERS[name](target, args.small, args.threads)
        except Exception as exc:  # noqa: BLE001 - one bad database must not stop the others
            _say("  ERROR installing %s: %s" % (name, exc))
            ok = False
        value = os.path.join(target, db.key_file) if db.key_file else target
        if db.key_glob:
            # The path is not known until the download names it; take the newest match so an
            # upgraded metapackage beside an old one points at the new one.
            matches = sorted(glob.glob(os.path.join(target, db.key_glob)))
            if matches:
                value = matches[-1]
                if len(matches) > 1:
                    _say(
                        "  %s holds %d versions of this database; configuring %s"
                        % (target, len(matches), os.path.basename(value))
                    )
        return name, value, db.key, ok

    workers = max(1, min(args.threads, len(names)))
    print(
        "\nInstalling %d database(s) into %s with %d worker(s).\n" % (len(names), args.dir, workers)
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for name, value, key, ok in pool.map(install, names):
            if ok:
                _say("  OK: %s installed" % name)
                if key:
                    installed[key] = value
            else:
                _say("  FAILED: %s" % name)
                failed.append(name)

    if installed and not args.no_config_update:
        path = _config_path(args.config)
        _set_config_keys(path, installed)
        print("\nUpdated %s:" % path)
        for key, value in installed.items():
            print("  %s = %s" % (key, value))

    print()
    if failed:
        print("%d database(s) failed: %s" % (len(failed), ", ".join(failed)))
        return 1
    print("All requested databases are installed. Verify with: metawrap2 check\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
