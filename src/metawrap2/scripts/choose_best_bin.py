"""Pick the best version of each bin from a reassembly ``.stats`` file.

``reassemble_bins`` produces up to three versions of every bin - the original (``.orig``) and
the two reassemblies (``.strict`` and ``.permissive``) - scores all of them with CheckM in one
go, and then has to choose one per bin. This module makes that choice.

The score deliberately weights contamination more heavily than completeness, because a
slightly less complete genome is usually more useful than a slightly more contaminated one::

    score = completeness + 5 * (100 - contamination)

Ties are broken by N50, preferring the less fragmented assembly. Versions failing the
requested completeness/contamination thresholds are not considered at all, so a bin whose
every version is poor simply drops out.

The ranking is identical to the original ``choose_best_bin.py``; this is a rewrite of
that script's structure (module-level ``sys.argv`` parsing, tab indentation, mutable
three-element lists) into something testable and readable.
"""

from __future__ import annotations

import sys
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

#: Contamination is penalised 5x relative to completeness (matches metawrap2.utils).
CONTAMINATION_WEIGHT = 5.0

#: Column positions in a MetaWrap2 ``.stats`` file (see metawrap2.checkm.STATS_HEADER).
_BIN, _COMPLETENESS, _CONTAMINATION, _N50 = 0, 1, 2, 5


class BinVersion(NamedTuple):
    """One scored candidate version of a bin, e.g. ``bin.3`` in its ``.strict`` form."""

    bin_name: str  # "bin.3"
    style: str  # "orig" | "strict" | "permissive"
    completeness: float
    contamination: float
    n50: int

    @property
    def full_name(self) -> str:
        """The name as it appears in the stats file and on disk, e.g. ``bin.3.strict``."""
        return "%s.%s" % (self.bin_name, self.style)

    @property
    def score(self) -> float:
        return self.completeness + CONTAMINATION_WEIGHT * (100.0 - self.contamination)

    def beats(self, other: BinVersion) -> bool:
        """True if this version should be preferred over *other* (N50 breaks a tie)."""
        if self.score != other.score:
            return self.score > other.score
        return self.n50 > other.n50


def parse_stats(lines: Iterable[str]) -> List[BinVersion]:
    """Parse ``.stats`` lines into :class:`BinVersion` records, skipping the header."""
    versions: List[BinVersion] = []
    for line in lines:
        if "contamination" in line:  # header row
            continue
        cut = line.rstrip("\n").split("\t")
        if len(cut) <= _N50:
            continue
        name = cut[_BIN]
        # "bin.3.strict" -> ("bin.3", "strict"); a name with no suffix keeps the whole name.
        bin_name, _, style = name.rpartition(".")
        if not bin_name:
            bin_name, style = name, ""
        try:
            versions.append(
                BinVersion(
                    bin_name=bin_name,
                    style=style,
                    completeness=float(cut[_COMPLETENESS]),
                    contamination=float(cut[_CONTAMINATION]),
                    n50=int(cut[_N50]),
                )
            )
        except ValueError:  # a row we cannot read is not a candidate
            continue
    return versions


def choose_best(
    versions: Iterable[BinVersion], min_completeness: float, max_contamination: float
) -> Dict[str, BinVersion]:
    """Return the winning version per bin, considering only versions passing the thresholds."""
    best: Dict[str, BinVersion] = {}
    for version in versions:
        if version.completeness < min_completeness or version.contamination > max_contamination:
            continue
        current = best.get(version.bin_name)
        if current is None or version.beats(current):
            best[version.bin_name] = version
    return best


def best_bin_names(stats_path: str, min_completeness: float, max_contamination: float) -> List[str]:
    """The chosen ``<bin>.<style>`` name for every bin that qualifies, sorted."""
    with open(stats_path) as fh:
        versions = parse_stats(fh)
    chosen = choose_best(versions, min_completeness, max_contamination)
    return sorted(v.full_name for v in chosen.values())


def summarize(
    stats_path: str, min_completeness: float, max_contamination: float
) -> Tuple[int, int, int]:
    """Count how many bins were won by (orig, strict, permissive), for the run summary."""
    with open(stats_path) as fh:
        versions = parse_stats(fh)
    chosen = choose_best(versions, min_completeness, max_contamination).values()
    styles = [v.style for v in chosen]
    return (styles.count("orig"), styles.count("strict"), styles.count("permissive"))


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        sys.stderr.write(
            "usage: python -m metawrap2.scripts.choose_best_bin "
            "<reassembled_bins.stats> <min completeness> <max contamination>\n"
        )
        return 2
    stats_path, min_completeness, max_contamination = argv[0], float(argv[1]), float(argv[2])
    for name in best_bin_names(stats_path, min_completeness, max_contamination):
        print(name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
