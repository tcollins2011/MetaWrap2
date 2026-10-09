"""The short-read aligners MetaWrap2 offers, and why minimap2 is the default.

Two modules map reads back to an assembly to measure per-contig coverage - `binning`, to build the
depth tables the binners read, and `blobology`, to put the abundance on a blobplot's y-axis. Neither
wants an alignment for its own sake. Both want coverage, which makes the aligner an implementation
detail of a coverage estimator, and means the same evidence settles the choice in both places. The
registry lives here so there is one answer rather than two that drift.

minimap2 is the default on measurement rather than novelty. Mapping one CAMI2 gut sample back to a
210k-contig assembly and binning the result with metaBAT2, minimap2 and bwa give bins of
indistinguishable quality - AMBER F1 0.4519 against 0.4493, recall 0.3118 against 0.3103, precision
0.8204 against 0.8133, all four fractionally in minimap2's favour and all four inside noise - for
1.37 CPU-hours against bwa's 7.77. A 5.7x saving with nothing traded away, and it keeps the BAM, so
CONCOCT is served too.

Note that minimap2 is *not* the arm that agrees most closely with bwa on the depth numbers
themselves: it ranks contigs almost identically (Spearman 0.984) but differs on magnitude (Pearson
0.946). That is the right thing to tolerate. bwa is the incumbent, not the truth, so agreement with
it is evidence about a replacement rather than the property being optimised - and the bins, which
are what the modules exist to produce, came out the same or better.

bwa-mem2 is offered but not chosen: it is bit-identical to bwa on the depth vector (Pearson 1.0000,
L1 0.000, the same 22,929 non-zero contigs) at 1.5x the speed - not the advertised 2-3x - but wants
23.6 GB of RAM against bwa's 9.8 and a 5.3 GB index against 1.7 GB. A reasonable trade to opt into
and a poor one to impose.

bowtie2 is not offered. Measured alongside these it cost what bwa cost (7.44 CPU-hours, of which
2.58 is index building alone) and agreed with bwa less closely than minimap2 did (Pearson 0.825), so
no input makes it the right answer. It was blobology's only aligner and a declared dependency of
binning that binning never actually called.

Caveat before treating the ranking as universal: this was one sample against CAMI2's pooled
gold-standard assembly, which is far less fragmented than a real metaSPAdes assembly, and
fragmentation is exactly what separates aligners.
"""

from __future__ import annotations

from typing import Dict

MINIMAP2_INDEX = "minimap2 -x sr -d {index} {assembly}"
MINIMAP2_PAIRED = "minimap2 -ax sr -t {threads} {index} {r1} {r2}"
MINIMAP2_SINGLE = "minimap2 -ax sr -t {threads} {index} {reads}"
BWA_INDEX = "bwa index {assembly}"
BWA_MEM_PAIRED = "bwa mem -v 1 -t {threads} {index} {r1} {r2}"
BWA_MEM_SINGLE = "bwa mem -t {threads} {index} {reads}"
BWA_MEM_INTERLEAVED = "bwa mem -v 1 -p -t {threads} {index} {reads}"
BWA_MEM2_INDEX = "bwa-mem2 index {assembly}"
BWA_MEM2_PAIRED = "bwa-mem2 mem -v 1 -t {threads} {index} {r1} {r2}"
BWA_MEM2_SINGLE = "bwa-mem2 mem -t {threads} {index} {reads}"
BWA_MEM2_INTERLEAVED = "bwa-mem2 mem -v 1 -p -t {threads} {index} {reads}"

#: ``--aligner`` choices. ``index_suffix`` is appended to the assembly copy to name the index the
#: map commands read; bwa and bwa-mem2 name their index files after the assembly itself, so for them
#: the prefix *is* the assembly and the suffix is empty.
ALIGNERS: Dict[str, Dict[str, str]] = {
    "minimap2": {
        "index_cmd": MINIMAP2_INDEX,
        "index_suffix": ".mmi",
        "paired": MINIMAP2_PAIRED,
        "single": MINIMAP2_SINGLE,
        # No -p equivalent to pass: the sr preset pairs consecutive records by itself, and does so
        # whether or not the read names carry /1 and /2 suffixes (checked both ways - the mates
        # come back flagged 0x1 with 0x40/0x80 set, and the suffixes are stripped from QNAME).
        "interleaved": MINIMAP2_SINGLE,
    },
    "bwa": {
        "index_cmd": BWA_INDEX,
        "index_suffix": "",
        "paired": BWA_MEM_PAIRED,
        "single": BWA_MEM_SINGLE,
        "interleaved": BWA_MEM_INTERLEAVED,
    },
    "bwa-mem2": {
        "index_cmd": BWA_MEM2_INDEX,
        "index_suffix": "",
        "paired": BWA_MEM2_PAIRED,
        "single": BWA_MEM2_SINGLE,
        "interleaved": BWA_MEM2_INTERLEAVED,
    },
}

DEFAULT_ALIGNER = "minimap2"

#: Every index file the three aligners write next to the assembly. Collected here because a module
#: that deletes spent alignments has to know what to sweep up, and the list is easy to get wrong:
#: minimap2's single ``.mmi`` was 479 MB for a 154 MB assembly, larger than the BAM beside it.
INDEX_SUFFIXES = (
    ".mmi",  # minimap2
    ".amb",  # bwa and bwa-mem2 from here down
    ".ann",
    ".bwt",
    ".pac",
    ".sa",
    ".0123",  # bwa-mem2 only from here down
    ".bwt.2bit.64",
)


def index_path(assembly: str, aligner: str) -> str:
    """Where *aligner* puts its index of *assembly* - a real file for minimap2, a prefix for bwa."""
    return assembly + ALIGNERS[aligner]["index_suffix"]
