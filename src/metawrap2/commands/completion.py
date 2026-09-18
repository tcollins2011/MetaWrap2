"""`metawrap2 completion [bash|zsh]` - print a shell completion script.

Install with e.g.:  metawrap2 completion bash >> ~/.bashrc   (or into a completions dir)
Completes the module/utility names and the global flags.
"""

from __future__ import annotations

from typing import List

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


def main(argv: List[str]) -> int:
    from ..cli import MODULES, UTILITIES

    shell = argv[0] if argv else "bash"
    subs = " ".join(sorted(list(MODULES) + list(UTILITIES)))
    if shell == "zsh":
        print(_ZSH % {"subs": subs})
    else:
        print(_BASH % {"subs": subs})
    return 0
