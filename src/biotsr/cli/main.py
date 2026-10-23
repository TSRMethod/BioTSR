"""``main(argv)`` -- parse, dispatch, translate exceptions into exit codes.

This is the only place in the CLI that catches a broad exception, and it does so
for one reason: a scientist who mistyped a chain should see one line naming the
chain, not a forty-line traceback through the selection layer.  The distinction
it draws is the one the exit codes draw:

* a :class:`~biotsr.cli.errors.CLIError` and the package's own ``ValueError`` /
  ``KeyError`` / ``FileNotFoundError`` are **the request being wrong** -- printed
  as one message, exit 2;
* anything else is **the program being wrong** -- printed with its type, exit 1,
  and with a traceback under ``--debug``.

The package's error messages are printed verbatim.  They were written for a
human, they name the method and what it does support, and re-wording them here
would create a second set that drifts.
"""

from __future__ import annotations

import sys
from typing import Optional, Sequence

from .. import TSRError
from . import commands
from .console import NORMAL, QUIET, VERBOSE, Console
from .errors import EXIT_INTERNAL, EXIT_INTERRUPTED, EXIT_USAGE, CLIError
from .parser import build_parser

__all__ = ["main"]


def _subparser(parser, name):
    """The sub-parser argparse built for ``name``, so its help can be printed."""
    for action in parser._actions:
        choices = getattr(action, "choices", None)
        if isinstance(choices, dict) and name in choices:
            return choices[name]
    return parser

#: Exceptions that mean *the request was wrong*, not *biotsr is broken*.
#:
#: ``TSRError`` is the package's own hierarchy (``biotsr/domain/errors.py``) and is
#: the answer to the question this tuple used to answer by hand: every refusal
#: the package raises on purpose is under it, every one carries a message written
#: for the person who typed the command, and adding a new one does not require
#: editing the CLI.  The taxonomy used to live here, in the wrong layer,
#: as a list of five exception types the package never knew about.
#:
#: ``ValueError`` and ``KeyError`` remain, and not as legacy residue: the package
#: still raises bare ones at ~70 sites that were deliberately not swept
#: (a rewrite no test could verify).  They stay until those sites are converted,
#: and ``TSRError`` subclasses ``ValueError`` anyway, so the two overlap rather
#: than compete.
#:
#: The OS errors are here because a missing file or an unreadable directory is
#: the user's path being wrong, which is exit 2, not a defect in biotsr.
_USER_ERRORS = (
    TSRError,
    ValueError,
    KeyError,
    FileNotFoundError,
    IsADirectoryError,
    NotADirectoryError,
    PermissionError,
)

_DISPATCH = {
    "methods": commands.cmd_methods,
    "run": commands.cmd_run,
    "batch": commands.cmd_batch,
    "fetch": commands.cmd_fetch,
    "describe": commands.cmd_describe,
}


def _option_dests(parser, table=None):
    """Every option string of ``parser`` and its sub-parsers -> its ``dest``."""
    table = {} if table is None else table
    for action in parser._actions:
        for option in action.option_strings:
            table.setdefault(option, action.dest)
        choices = getattr(action, "choices", None)
        if isinstance(choices, dict):
            for sub in choices.values():
                if hasattr(sub, "_actions"):
                    _option_dests(sub, table)
    return table


def _explicit_dests(parser, argv) -> set:
    """The ``dest`` of every option the user actually typed.

    argparse fills in defaults, so a namespace alone cannot say whether
    ``--timeout 30`` was typed or is the default -- and a command that refuses
    an option it would ignore must know.  Long options may be abbreviated, as
    argparse allows; an abbreviation is resolved the way argparse resolves it.
    """
    table = _option_dests(parser)
    longs = [o for o in table if o.startswith("--")]
    found = set()
    for token in argv:
        if token == "--":
            break
        if token.startswith("--"):
            name = token.split("=", 1)[0]
            if name in table:
                found.add(table[name])
            else:
                matches = {table[o] for o in longs if o.startswith(name)}
                if len(matches) == 1:
                    found.update(matches)
        elif token.startswith("-") and len(token) >= 2 and not token[1].isdigit():
            if token[:2] in table:
                found.add(table[token[:2]])
    return found


class _SigtermAsInterrupt(object):
    """SIGTERM -- what a scheduler or ``kill`` sends -- is handled like Ctrl-C.

    Installed for the whole command, so cleanup that runs on Ctrl-C (removing
    a staging directory, writing a batch summary or a fetch report) also runs
    when a job is terminated.  ``SIGKILL`` and power loss cannot be caught;
    what they leave behind is recognised and cleaned up by a later run.
    """

    def __enter__(self):
        import signal

        self.previous = None
        try:
            self.previous = signal.signal(signal.SIGTERM, self._handler)
        except (ValueError, OSError, AttributeError):  # not the main thread
            self.previous = None
        return self

    @staticmethod
    def _handler(signum, frame):
        raise KeyboardInterrupt()

    def __exit__(self, *exc):
        if self.previous is not None:
            import signal

            signal.signal(signal.SIGTERM, self.previous)
        return False


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    args._argv = argv
    args._explicit = _explicit_dests(parser, argv)

    if not getattr(args, "command", None):
        parser.print_help()
        return EXIT_USAGE

    level = QUIET if getattr(args, "quiet", False) else (
        VERBOSE if getattr(args, "verbose", False) else NORMAL
    )
    console = Console(level)

    if args.command == "analyze":
        if not getattr(args, "analysis", None):
            # `print_help()`, not `parse_args([..., "--help"])`:
            # argparse's own --help exits 0 through SystemExit, and `biotsr analyze`
            # with no operation is a malformed command, not a request for help.
            # It has to leave through the same door as every other usage error --
            # exit 2, and a return rather than an exception, so an in-process
            # caller sees a code like everyone else.
            _subparser(parser, "analyze").print_help()
            console.error("`analyze` needs an operation")
            console.hint("operations: common-keys")
            return EXIT_USAGE
        handler = commands.cmd_analyze
    else:
        handler = _DISPATCH[args.command]

    try:
        with _SigtermAsInterrupt():
            return handler(args, console)
    except CLIError as exc:
        console.error(exc.message)
        console.hint(exc.hint)
        return exc.exit_code
    except _USER_ERRORS as exc:
        message = str(exc)
        if isinstance(exc, KeyError) and exc.args and isinstance(exc.args[0], str):
            # KeyError's str() adds quoting of its own; the package's KeyErrors
            # carry a written sentence, and printing it inside quotes is noise.
            message = exc.args[0]
        console.error(message)
        if getattr(args, "debug", False):
            import traceback

            traceback.print_exc()
        return EXIT_USAGE
    except KeyboardInterrupt:
        console.error("interrupted")
        return EXIT_INTERRUPTED
    except Exception as exc:  # noqa: BLE001 -- the internal-error boundary
        console.error("%s: %s" % (type(exc).__name__, exc))
        console.hint(
            "this is a defect in biotsr, not in the command. Re-run with --debug "
            "for a traceback."
        )
        if getattr(args, "debug", False):
            import traceback

            traceback.print_exc()
        return EXIT_INTERNAL
