#!/usr/bin/env python
# Select FASTA sequences from a file. The blobology module uses this only to pull
# a random subset of contigs ("-s r -n N"), which is what this implements.
#
# Replaces the original fastaqual_select.pl (Sujai Kumar, blaxterlab/blobology) for the
# options the blobology module actually uses.
import argparse
import random
import sys

from metawrap2.io.seqio import iter_fasta


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("-f", "--fastafile", required=True)
    ap.add_argument(
        "-s", "--sort", default="S", help="R = random order; anything else keeps input order"
    )
    ap.add_argument(
        "-n", "--numfasta", type=int, default=0, help="number of sequences (default all)"
    )
    ap.add_argument("-l", "--length", type=int, default=0, help="minimum sequence length")
    args = ap.parse_args()

    records = [(h, s) for h, s in iter_fasta(args.fastafile) if len(s) >= args.length]

    if args.sort.upper() == "R":
        random.shuffle(records)

    if args.numfasta:
        records = records[: args.numfasta]

    for header, seq in records:
        sys.stdout.write(">%s\n%s\n" % (header, seq))


if __name__ == "__main__":
    main()
