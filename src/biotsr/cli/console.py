"""Terminal output: three levels, one stream discipline.

Results go to stdout, diagnostics go to stderr, and nothing interleaves them --
so ``biotsr methods --json | jq`` works and ``biotsr batch ... 2> errors.log`` keeps the
completion summary on the terminal.

Three levels, because two is not enough and four is a configuration file:

* ``--quiet`` -- nothing but errors.  For a job array where 5,000 identical
  success lines are noise.
* default -- one line per command and a completion summary.  In batch, one line
  per *failed* sample and none per successful one.
* ``--verbose`` -- one line per sample, and the package's own warnings
  (insertion-code skips, missing atoms, degenerate triangles).  Never per
  triangle: a 1,000-residue chain has 166 million candidates and no terminal is
  a reasonable place to put them.
"""

from __future__ import annotations

import sys

__all__ = ["Console", "QUIET", "NORMAL", "VERBOSE"]

QUIET = -1
NORMAL = 0
VERBOSE = 1


class Console(object):
    def __init__(self, level: int = NORMAL, stdout=None, stderr=None):
        self.level = level
        self._out = stdout if stdout is not None else sys.stdout
        self._err = stderr if stderr is not None else sys.stderr

    def out(self, text: str = "") -> None:
        """Result output.  Always written -- it is what the command produced."""
        self._out.write(text + "\n")

    def say(self, text: str = "") -> None:
        """Progress and summary.  Suppressed by ``--quiet``."""
        if self.level >= NORMAL:
            self._err.write(text + "\n")

    def detail(self, text: str = "") -> None:
        """Per-sample detail and package warnings.  ``--verbose`` only."""
        if self.level >= VERBOSE:
            self._err.write(text + "\n")

    def warn(self, text: str) -> None:
        """A recoverable problem.  Suppressed only by ``--quiet``."""
        if self.level >= NORMAL:
            self._err.write("biotsr: warning: " + text + "\n")

    def error(self, text: str) -> None:
        """Always written, at every level."""
        self._err.write("biotsr: error: " + text + "\n")

    def hint(self, text: str) -> None:
        if text:
            self._err.write("  " + text + "\n")
