"""MetaWrap2 classify_bins module: assign a taxonomy to each genomic bin.

All contigs from the bins are aligned to the NCBI nt database with MEGABLAST; hits are
pruned to ranked tax ids, then taxator-tk (megan-lca) assigns per-contig taxonomy, which is
consolidated per contig with binner/taxknife. Finally each bin's consensus taxonomy is
called by building a length-weighted taxonomy tree and following the majority path.

This file is meant to be read and edited: the commands each tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block. Change a flag there and it
takes effect — no need to follow the orchestration logic underneath.
"""

from __future__ import annotations

import argparse
import glob
import os
from typing import List

from .. import blastdb as _blastdb
from ..config import OPTIONAL_ENVS, load_settings
from ..constants import BIN_EXTENSION, FASTA_EXTENSIONS, TSV_SUFFIX
from ..io.seqio import contig_id, iter_fasta
from ..scripts import classify_bins as _classify_helper
from ..scripts import prune_blast_hits
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
)

CONDA_ENV = "metawrap2-classify_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BLASTN = (
    "blastn -task megablast -num_threads {threads} -db {blastdb} "
    "-outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length' "
    "-query {query}"
)
# taxator/binner/taxknife read stdin and are chained through pipes, so they run via bash -c.
# All three are taxator-tk tools and all three need TAXATORTK_TAXONOMY_NCBI - only taxator
# had it set, so `binner` failed with "Specify the folder containing the NCBI taxonomy dump
# files as TAXATORTK_TAXONOMY_NCBI environment variable" after the slow megablast finished.
_TAXONOMY_ENV = "TAXATORTK_TAXONOMY_NCBI={taxdump}"
TAXATOR = (
    _TAXONOMY_ENV
    + " taxator -a megan-lca -t {taxator_t} -e {taxator_e} -g {mapping} < {input} > {output}"
)
BINNER = (
    "sort -k1,1 {predictions} | "
    + _TAXONOMY_ENV
    + " binner -n classification -i genus:{genus_cutoff} > {output}"
)
TAXKNIFE = (
    _TAXONOMY_ENV
    + " taxknife -f 2 --mode annotate -s path < {binned} | grep -v 'Could not' | cut -f1,2 > {output}"
)
# --gtdbtk path: one call replaces the whole megablast -> taxator -> binner -> taxknife chain.
GTDBTK = (
    "gtdbtk classify_wf --genome_dir {bins} --out_dir {out} -x "
    + BIN_EXTENSION.lstrip(".")
    + " --cpus {threads} --skip_ani_screen"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
TAXATOR_LCA_THRESHOLD = 0.3  # taxator -t (megan-lca support threshold)
TAXATOR_MAX_EVALUE = 0.01  # taxator -e
BINNER_GENUS_CUTOFF = 0.6  # binner -i genus:<cutoff>
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 classify_bins",
        usage="metawrap2 classify_bins [options] -b bin_folder -o output_dir",
        description="Assign a consensus taxonomy to each bin with MEGABLAST + taxator-tk.",
    )
    p.add_argument(
        "-b", "--bins", required=True, help="folder with the bins to be classified (fasta)"
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "--gtdbtk",
        action="store_true",
        help="classify with GTDB-Tk against GTDB instead of MEGABLAST+taxator-tk "
        "(needs the metawrap2-classify_bins-gtdbtk env and GTDBTK_DATA_PATH)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _concat_all_contigs(bin_folder: str, all_contigs: str) -> None:
    """Merge every bin's contigs into one fasta (supports .gz bins)."""
    if os.path.isfile(all_contigs):
        os.remove(all_contigs)
    with open(all_contigs, "w") as out:
        for f in sorted(os.listdir(bin_folder)):
            out.writelines(
                ">%s\n%s\n" % (contig_id(header), seq)
                for header, seq in iter_fasta(os.path.join(bin_folder, f))
            )


def _write_pruned_columns(pruned: str, tab_out: str) -> None:
    """Drop the staxids column (keep cols 1-5,7-12) — the input taxator consumes."""
    with open(pruned) as fh, open(tab_out, "w") as out:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            out.write("\t".join(cols[i] for i in (0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11)) + "\n")


def _write_mapping(pruned: str, mapping_out: str) -> None:
    """Write the sseqid -> taxid mapping file (cols 5,6)."""
    with open(pruned) as fh, open(mapping_out, "w") as out:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            out.write("%s\t%s\n" % (cols[4], cols[5]))


def _classify_with_gtdbtk(args, settings, ckpt) -> None:
    """Classify bins with GTDB-Tk, writing the same bin_taxonomy.tab contract.

    GTDB-Tk already does the marker placement and taxonomy call, so none of the
    megablast/taxator-tk chain is needed. Its ``classify_wf`` writes per-domain summary
    files; we merge them into the single two-column bin_taxonomy.tab the rest of MetaWrap2
    (and users' downstream scripts) expect from this module.
    """
    env = OPTIONAL_ENVS["gtdbtk"] if settings.use_conda_envs else None
    ensure_dir(args.output)
    gtdb_out = os.path.join(args.output, "gtdbtk_out")

    announcement("CLASSIFYING BINS WITH GTDB-TK")
    if ckpt.todo("gtdbtk"):
        comm("running gtdbtk classify_wf on %s" % args.bins)
        run(
            GTDBTK.format(bins=args.bins, out=gtdb_out, threads=args.threads),
            env=env,
            tool="gtdbtk",
            hint="GTDB-Tk needs its reference data; set GTDBTK_DATA_PATH in the "
            "metawrap2-classify_bins-gtdbtk env (gtdbtk download-db.sh).",
        )
        ckpt.done("gtdbtk")
    else:
        comm("skipping gtdbtk classify_wf (already done; --resume)")

    summaries = sorted(
        glob.glob(os.path.join(gtdb_out, "classify", "gtdbtk.*.summary.tsv"))
    ) or sorted(glob.glob(os.path.join(gtdb_out, "gtdbtk.*.summary.tsv")))
    if dry_run():
        comm(
            "(dry run) would merge GTDB-Tk's summaries into %s"
            % os.path.join(args.output, "bin_taxonomy" + TSV_SUFFIX)
        )
        return
    if not summaries:
        error("GTDB-Tk produced no summary file under %s. Exiting." % gtdb_out)

    bin_taxonomy = os.path.join(args.output, "bin_taxonomy" + TSV_SUFFIX)
    n = 0
    with open(bin_taxonomy, "w") as out:
        for summary in summaries:
            with open(summary) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                try:
                    i_bin = header.index("user_genome")
                    i_tax = header.index("classification")
                except ValueError:
                    error("Unexpected GTDB-Tk summary columns in %s: %s" % (summary, header))
                for line in fh:
                    if not line.strip():
                        continue
                    cut = line.rstrip("\n").split("\t")
                    line_out = "%s\t%s" % (cut[i_bin], cut[i_tax])
                    out.write(line_out + "\n")
                    print(line_out)
                    n += 1
    comm("classified %d bins; consensus taxonomy of each bin is in %s" % (n, bin_taxonomy))


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not os.path.isdir(args.bins):
        error("%s does not exist. Exiting..." % args.bins)

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    validate_inputs(
        "classify_bins",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                args.bins,
                "bin folder (-b)",
            )
        ],
    )
    env = env_for("classify_bins", settings)
    rec = start_run("classify_bins", args, env, settings, inputs=[args.bins])

    ckpt = make_checkpoint(args.output)
    try:
        if args.gtdbtk:
            _classify_with_gtdbtk(args, settings, ckpt)
            announcement("BIN CLASSIFICATION PIPELINE FINISHED SUCCESSFULLY!!!")
            return 0

        blastdb = settings.db("BLASTDB")
        taxdump = settings.db("TAXDUMP")

        ok, problem = _blastdb.find(settings)
        if not ok:
            error(problem or "the configured BLAST database cannot be used")
        blastdb = _blastdb.db_path(settings)
        if not os.path.isfile(os.path.join(taxdump, "names.dmp")):
            error(
                "%s/names.dmp does not exist. Set TAXDUMP in metawrap2.toml to your NCBI_tax "
                "database, or download it (pub/taxonomy/taxdump.tar.gz)." % taxdump
            )

        announcement("ALIGN CONTIGS TO DATABASE WITH MEGABLAST")
        ensure_dir(args.output)
        comm("setting up output folder %s and merging contigs from all bins..." % args.output)
        all_contigs = os.path.join(args.output, "all_contigs.fa")
        if not dry_run():
            _concat_all_contigs(args.bins, all_contigs)
            if not (os.path.isfile(all_contigs) and os.path.getsize(all_contigs)):
                error(
                    "something went wrong with joining files in %s into %s"
                    % (args.bins, all_contigs)
                )

        raw = os.path.join(args.output, "megablast_out.raw" + TSV_SUFFIX)
        tab = os.path.join(args.output, "megablast_out" + TSV_SUFFIX)
        mapping = os.path.join(args.output, "mapping.taxids" + TSV_SUFFIX)
        if ckpt.todo("megablast"):
            if not dry_run() and os.path.isfile(raw) and os.path.getsize(raw):
                comm("megablast alignment already done. Skipping...")
            else:
                comm(
                    "aligning %s to %s database with MEGABLAST. This is the longest step - please "
                    "be patient. You may look at the classification progress in %s"
                    % (all_contigs, blastdb, raw)
                )
                run(
                    BLASTN.format(threads=args.threads, blastdb=blastdb, query=all_contigs),
                    env=env,
                    tool="blastn",
                    log_path=raw,
                    hint="Failed to run megablast.",
                )

            comm("removing unnecessary lines that lead to bad tax IDs (without a proper rank)")
            pruned = os.path.join(args.output, "megablast_out.pruned" + TSV_SUFFIX)
            if not dry_run():
                with open(pruned, "w") as out:
                    out.writelines(
                        line + "\n"
                        for line in prune_blast_hits.prune(os.path.join(taxdump, "nodes.dmp"), raw)
                    )
                _write_pruned_columns(pruned, tab)

                comm("making mapping file")
                _write_mapping(pruned, mapping)
            ckpt.done("megablast")
        else:
            comm("skipping MEGABLAST alignment (already done; --resume)")

        contig_taxonomy = os.path.join(args.output, "contig_taxonomy" + TSV_SUFFIX)
        if ckpt.todo("taxonomy"):
            announcement("GET TAXONOMY FROM MEGABLAST OUTPUT WITH TAXATOR-TK")
            comm("pulling out classifications with taxator")
            predictions = os.path.join(args.output, "predictions.gff3")
            run(
                [
                    "bash",
                    "-c",
                    TAXATOR.format(
                        taxdump=taxdump,
                        taxator_t=TAXATOR_LCA_THRESHOLD,
                        taxator_e=TAXATOR_MAX_EVALUE,
                        mapping=mapping,
                        input=tab,
                        output=predictions,
                    ),
                ],
                env=env,
                tool="taxator",
                hint="Failed to run taxator.",
            )

            comm("binning and consolidating classifications for each contig")
            binned = os.path.join(args.output, "binned_predictions.txt")
            run(
                [
                    "bash",
                    "-c",
                    BINNER.format(
                        taxdump=taxdump,
                        predictions=predictions,
                        genus_cutoff=BINNER_GENUS_CUTOFF,
                        output=binned,
                    ),
                ],
                env=env,
                tool="binner",
                cwd=args.output,
                hint="Failed to run binner.",
            )

            comm("pulling out full taxonomy path with taxknife")
            run(
                [
                    "bash",
                    "-c",
                    TAXKNIFE.format(taxdump=taxdump, binned=binned, output=contig_taxonomy),
                ],
                env=env,
                tool="taxknife",
                hint="Failed to extract full taxonomy path with taxknife.",
            )
            ckpt.done("taxonomy")
        else:
            comm("skipping taxator-tk classification (already done; --resume)")

        if ckpt.todo("consensus") and not dry_run():
            comm("finding consensus taxonomy for each bin")
            bin_taxonomy = os.path.join(args.output, "bin_taxonomy" + TSV_SUFFIX)
            with open(bin_taxonomy, "w") as out:
                for filename, consensus in _classify_helper.consensus_taxonomy(
                    contig_taxonomy, args.bins
                ):
                    line = filename + "\t" + consensus
                    out.write(line + "\n")
                    print(line)
            comm("you will find the consensus taxonomy of each bin in %s" % bin_taxonomy)
            ckpt.done("consensus")
        elif dry_run():
            comm(
                "(dry run) would write the per-bin consensus taxonomy to %s"
                % os.path.join(args.output, "bin_taxonomy" + TSV_SUFFIX)
            )
        else:
            comm("skipping per-bin consensus taxonomy (already done; --resume)")

        announcement("BIN CLASSIFICATION PIPELINE FINISHED SUCCESSFULLY!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
