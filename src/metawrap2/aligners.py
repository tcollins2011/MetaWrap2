"""Short-read aligners offered by MetaWrap2, and the default (minimap2).

`binning` and `blobology` both map reads back to an assembly to get per-contig coverage: binning
for the binners' depth tables, blobology for the blobplot y-axis. The aligner is only used to
estimate coverage, so one registry serves both modules.

minimap2 is the default based on a benchmark: one CAMI2 gut sample mapped to a 210k-contig assembly
and binned with metaBAT2. minimap2 and bwa gave bins of the same quality within noise (AMBER F1
0.4519 vs 0.4493, recall 0.3118 vs 0.3103, precision 0.8204 vs 0.8133), with minimap2 using 1.37
CPU-hours against bwa's 7.77 (5.7x less). It writes a BAM, so CONCOCT still works.

minimap2's depths rank contigs almost the same as bwa's (Spearman 0.984) but differ in magnitude
(Pearson 0.946). bwa is a reference point, not ground truth, and the resulting bins were the same
or better.

bwa-mem2 is available but not the default. Its depth vector is identical to bwa's (Pearson 1.0000,
L1 0.000, the same 22,929 non-zero contigs) at 1.5x the speed (the advertised figure is 2-3x), but
it needs 23.6 GB of RAM against bwa's 9.8 and a 5.3 GB index against 1.7 GB.

bowtie2 is not offered. In the same benchmark it cost about what bwa did (7.44 CPU-hours, 2.58 of
them building the index) and agreed with bwa less closely than minimap2 (Pearson 0.825). It was
previously blobology's only aligner and a declared, unused dependency of binning.

Limitation: the benchmark used one sample against CAMI2's pooled gold-standard assembly, which is
much less fragmented than a typical metaSPAdes assembly, and fragmentation is where aligners differ
most.
"""

from __future__ import annotations

from typing import Dict, Tuple

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
        # No -p needed: the sr preset pairs consecutive records, with or without /1 and /2 name
        # suffixes (mates come back flagged 0x1 with 0x40/0x80 set, suffixes stripped from QNAME).
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

#: Every index file the three aligners write next to the assembly, for modules that clean up after
#: alignment. minimap2's ``.mmi`` was 479 MB for a 154 MB assembly, larger than the BAM.
INDEX_SUFFIXES = (
    ".mmi",  # minimap2
    ".amb",  # bwa and bwa-mem2 from here down
    ".ann",
    ".bwt",
    ".pac",
    ".sa",
    ".0123",  # bwa-mem2 only from here down
    ".bwt.2bit.64",
    ".mw2index",  # MetaWrap2's record of the command that built the index (see INDEX_MARKER)
    ".mw2index.lock",
)

#: Written next to the reference after an index build succeeds, holding the exact index command.
#: An index is reused only when this matches the command about to run, so a half-written index or
#: one built with another minimap2 preset is rebuilt.
INDEX_MARKER = ".mw2index"


#: The files that make up each aligner's index. INDEX_SUFFIXES above is the
#: union, which is what a sweep needs; deciding whether an index is complete needs the per-aligner
#: set, because bwa writes five files under one prefix and minimap2 writes one.
INDEX_FILES: Dict[str, Tuple[str, ...]] = {
    "minimap2": (".mmi",),
    "bwa": (".amb", ".ann", ".bwt", ".pac", ".sa"),
    "bwa-mem2": (".amb", ".ann", ".pac", ".0123", ".bwt.2bit.64"),
}


def index_path(assembly: str, aligner: str) -> str:
    """Path of *aligner*'s index of *assembly*: a file for minimap2, a prefix for bwa."""
    return assembly + ALIGNERS[aligner]["index_suffix"]
