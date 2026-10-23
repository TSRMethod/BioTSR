"""Exit codes and the one error type the CLI raises on purpose.

The taxonomy answers the only question a shell script actually asks: *did it
work, was it my fault, or did the program break?* -- plus the one extra
distinction a batch run needs, because "24 of 25 samples succeeded" is neither
success nor failure and a job scheduler must be able to tell it from both, and
the conventional code for a run stopped by Ctrl-C or SIGTERM.

``2`` is user error rather than ``1`` because ``argparse`` already exits ``2`` on a
malformed command line, and having a bad flag and a nonexistent chain report
different codes would be a distinction without a difference.
"""

from __future__ import annotations

__all__ = [
    "EXIT_OK",
    "EXIT_INTERNAL",
    "EXIT_USAGE",
    "EXIT_PARTIAL",
    "EXIT_INCOMPLETE",
    "EXIT_INTERRUPTED",
    "EXIT_CODES",
    "CLIError",
]

#: Everything asked for was produced.
EXIT_OK = 0
#: An unexpected failure inside the package: a bug, not a bad request.  The
#: message names the exception type; ``--debug`` prints the traceback.
EXIT_INTERNAL = 1
#: The request could not be honoured: bad path, unknown method, unknown chain,
#: illegal modifier combination, malformed manifest, output collision.  Also the
#: code ``argparse`` uses for a malformed command line.
EXIT_USAGE = 2
#: A batch run finished, some samples failed, and every failure is recorded in
#: the failure manifest.  Never returned by a single-structure run.
EXIT_PARTIAL = 3
#: Everything finished and nothing failed, but at least one result is
#: ``partial``: its selection did not get everything it asked for (a requested
#: residue, atom or molecule absent, a modified residue the method cannot
#: read, one the completeness rule rejected), or a resumed run kept such a
#: result.  The results are written and usable; the code exists so a script
#: can tell them from complete success (ADR-106).  A batch with a failure
#: returns :data:`EXIT_PARTIAL` instead.
EXIT_INCOMPLETE = 4
#: Stopped by Ctrl-C or SIGTERM (128 + SIGINT, the shell convention).  A batch
#: writes its summary first, marking each row interrupted or not started.
EXIT_INTERRUPTED = 130

EXIT_CODES = {
    EXIT_OK: "success",
    EXIT_INTERNAL: "unexpected internal error",
    EXIT_USAGE: "user, input or configuration error",
    EXIT_PARTIAL: "batch completed with failed samples",
    EXIT_INCOMPLETE: "completed, but some results are partial",
    EXIT_INTERRUPTED: "interrupted",
}


class CLIError(Exception):
    """A problem with the request, stated for the person who typed it.

    Raising this prints ``biotsr: error: <message>`` and exits :data:`EXIT_USAGE`.
    No traceback, because there is nothing in a traceback that helps somebody who
    misspelled a chain.  ``hint`` is printed on its own line when present and is
    for the next thing to try, not for restating the error.
    """

    def __init__(self, message: str, hint: str = "", exit_code: int = EXIT_USAGE):
        super().__init__(message)
        self.message = message
        self.hint = hint
        self.exit_code = exit_code
