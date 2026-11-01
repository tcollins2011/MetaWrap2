"""`metawrap2 completion [bash|zsh]` - print a shell completion script.

Install with, for example, ``metawrap2 completion bash >> ~/.bashrc`` (or write it into a
completions directory). Completes module and utility names and the global flags.
"""

from __future__ import annotations

import argparse
from typing import List

from ..cliformat import make_parser

_BASH = """\
_metawrap2() {
    local cur prev words
    cur="${COMP_WORDS[COMP_CWORD]}"
    local subs="%(subs)s"
    local globals="--dry-run --force --resume --help -h --version -v"
    if [ "$COMP_CWORD" -le 1 ]; then
        COMPREPLY=( $(compgen -W "$subs $globals" -- "$cur") )
    else
        COMPREPLY=( $(compgen -o default -W "$globals" -- "$cur") )
    fi
}
complete -F _metawrap2 metawrap2
"""

_ZSH = """\
#compdef metawrap2
_metawrap2() {
    local -a subs
    subs=(%(subs)s)
    if (( CURRENT <= 2 )); then
        compadd -- $subs --dry-run --force --resume --help --version
    else
        _files
    fi
}
compdef _metawrap2 metawrap2
"""


#: The shells there is a script for.
SHELLS = ("bash", "zsh")


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 completion",
        usage="metawrap2 completion [bash|zsh]",
        description="Print a shell completion script for module and utility names and the "
        "global flags. Install with, for example, `metawrap2 completion bash >> ~/.bashrc`.",
    )
    p.add_argument(
        "shell",
        nargs="?",
        default="bash",
        choices=SHELLS,
        help="the shell to print a script for (default: %(default)s)",
    )
    return p.parse_args(argv)


def main(argv: List[str]) -> int:
    from ..cli import MODULES, UTILITIES

    shell = _parse_args(argv).shell
    subs = " ".join(sorted(list(MODULES) + list(UTILITIES)))
    if shell == "zsh":
        print(_ZSH % {"subs": subs})
    else:
        print(_BASH % {"subs": subs})
    return 0
