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
import os
from typing import List

from ..config import load_settings
from ..io.seqio import iter_fasta
from ..scripts import classify_bins as _classify_helper
from ..scripts import prune_blast_hits
from ._common import (
    announcement, comm, ensure_dir, env_for, error, finish_run, make_checkpoint, run, start_run,
)

CONDA_ENV = "metawrap2-classify_bins"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
BLASTN = ("blastn -task megablast -num_threads {threads} -db {blastdb}/nt "
          "-outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length' "
          "-query {query}")
# taxator/binner/taxknife read stdin and are chained through pipes, so they run via bash -c.
TAXATOR  = "TAXATORTK_TAXONOMY_NCBI={taxdump} taxator -a megan-lca -t {taxator_t} -e {taxator_e} -g {mapping} < {input} > {output}"
BINNER   = "sort -k1,1 {predictions} | binner -n classification -i genus:{genus_cutoff} > {output}"
TAXKNIFE = "taxknife -f 2 --mode annotate -s path < {binned} | grep -v 'Could not' | cut -f1,2 > {output}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
TAXATOR_LCA_THRESHOLD = 0.3    # taxator -t (megan-lca support threshold)
TAXATOR_MAX_EVALUE    = 0.01   # taxator -e
BINNER_GENUS_CUTOFF   = 0.6    # binner -i genus:<cutoff>
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 classify_bins",
        usage="metawrap2 classify_bins [options] -b bin_folder -o output_dir",
        description="Assign a consensus taxonomy to each bin with MEGABLAST + taxator-tk.",
    )
    p.add_argument("-b", "--bins", required=True, help="folder with the bins to be classified (fasta)")
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument("-t", "--threads", type=int, default=1, help="number of threads (default 1)")
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _concat_all_contigs(bin_folder: str, all_contigs: str) -> None:
    """Merge every bin's contigs into one fasta (supports .gz bins)."""
    if os.path.isfile(all_contigs):
        os.remove(all_contigs)
    with open(all_contigs, "w") as out:
        for f in sorted(os.listdir(bin_folder)):
            for header, seq in iter_fasta(os.path.join(bin_folder, f)):
                out.write(">%s\n%s\n" % (header, seq))


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


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    if not os.path.isdir(args.bins):
        error("%s does not exist. Exiting..." % args.bins)

    settings = load_settings(args.config)
    env = env_for("classify_bins", settings)
    rec = start_run("classify_bins", args, env, settings, inputs=[args.bins])

    ckpt = make_checkpoint(args.output)
    try:
        blastdb = settings.db("BLASTDB")
        taxdump = settings.db("TAXDUMP")

        if not os.path.isfile(os.path.join(blastdb, "nt.00.nhd")):
            error("%s/nt.00.nhd does not exist. Set BLASTDB in metawrap2.toml to your NCBI_nt "
                  "database, or download it (blast/db/nt.*.tar.gz)." % blastdb)
        if not os.path.isfile(os.path.join(taxdump, "names.dmp")):
            error("%s/names.dmp does not exist. Set TAXDUMP in metawrap2.toml to your NCBI_tax "
                  "database, or download it (pub/taxonomy/taxdump.tar.gz)." % taxdump)

        announcement("ALIGN CONTIGS TO DATABASE WITH MEGABLAST")
        ensure_dir(args.output)
        comm("setting up output folder %s and merging contigs from all bins..." % args.output)
        all_contigs = os.path.join(args.output, "all_contigs.fa")
        _concat_all_contigs(args.bins, all_contigs)
        if not (os.path.isfile(all_contigs) and os.path.getsize(all_contigs)):
            error("something went wrong with joining files in %s into %s" % (args.bins, all_contigs))

        raw = os.path.join(args.output, "megablast_out.raw.tab")
        tab = os.path.join(args.output, "megablast_out.tab")
        mapping = os.path.join(args.output, "mapping.tax")
        if ckpt.todo("megablast"):
            if os.path.isfile(raw) and os.path.getsize(raw):
                comm("megablast alignment already done. Skipping...")
            else:
                comm("aligning %s to %s database with MEGABLAST. This is the longest step - please "
                     "be patient. You may look at the classification progress in %s"
                     % (all_contigs, blastdb, raw))
                run(BLASTN.format(threads=args.threads, blastdb=blastdb, query=all_contigs),
                    env=env, tool="blastn", log_path=raw, hint="Failed to run megablast.")

            comm("removing unnecessary lines that lead to bad tax IDs (without a proper rank)")
            pruned = os.path.join(args.output, "megablast_out.pruned.tab")
            with open(pruned, "w") as out:
                for line in prune_blast_hits.prune(os.path.join(taxdump, "nodes.dmp"), raw):
                    out.write(line + "\n")
            _write_pruned_columns(pruned, tab)

            comm("making mapping file")
            _write_mapping(pruned, mapping)
            ckpt.done("megablast")
        else:
            comm("skipping MEGABLAST alignment (already done; --resume)")

        contig_taxonomy = os.path.join(args.output, "contig_taxonomy.tab")
        if ckpt.todo("taxonomy"):
            announcement("GET TAXONOMY FROM MEGABLAST OUTPUT WITH TAXATOR-TK")
            comm("pulling out classifications with taxator")
            predictions = os.path.join(args.output, "predictions.gff3")
            run(["bash", "-c", TAXATOR.format(taxdump=taxdump, taxator_t=TAXATOR_LCA_THRESHOLD,
                                              taxator_e=TAXATOR_MAX_EVALUE, mapping=mapping,
                                              input=tab, output=predictions)],
                env=env, tool="taxator", hint="Failed to run taxator.")

            comm("binning and consolidating classifications for each contig")
            binned = os.path.join(args.output, "binned_predictions.txt")
            run(["bash", "-c", BINNER.format(predictions=predictions, genus_cutoff=BINNER_GENUS_CUTOFF,
                                             output=binned)],
                env=env, tool="binner", cwd=args.output, hint="Failed to run binner.")

            comm("pulling out full taxonomy path with taxknife")
            run(["bash", "-c", TAXKNIFE.format(binned=binned, output=contig_taxonomy)],
                env=env, tool="taxknife", hint="Failed to extract full taxonomy path with taxknife.")
            ckpt.done("taxonomy")
        else:
            comm("skipping taxator-tk classification (already done; --resume)")

        if ckpt.todo("consensus"):
            comm("finding consensus taxonomy for each bin")
            bin_taxonomy = os.path.join(args.output, "bin_taxonomy.tab")
            with open(bin_taxonomy, "w") as out:
                for filename, consensus in _classify_helper.consensus_taxonomy(contig_taxonomy, args.bins):
                    line = filename + "\t" + consensus
                    out.write(line + "\n")
                    print(line)
            comm("you will find the consensus taxonomy of each bin in %s" % bin_taxonomy)
            ckpt.done("consensus")
        else:
            comm("skipping per-bin consensus taxonomy (already done; --resume)")

        announcement("BIN CLASSIFICATION PIPELINE FINISHED SUCCESSFULLY!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
