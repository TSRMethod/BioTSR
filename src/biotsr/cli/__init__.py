"""The ``biotsr`` command line: a thin orchestration layer over the Python API.

What this package is allowed to do: read a request, call ``biotsr.api``, call
``biotsr.io``, report what happened.  What it must never do: geometry, ordering,
binning, key formulas, contact selection, SSE categories, or any default for a
scientific parameter.  Every number that reaches a key comes from a ``Profile``
or a ``Method``; the CLI's job is to name one, not to choose one.

A test asserts the first half of that (`tests/unit/test_cli_architecture.py`), and
the equivalence gate asserts the second: for every representative workflow, the
files the CLI writes are produced from a result that is identical, field for
field, to the one a direct ``analyze()`` call returns.

Entry points::

    biotsr ...            # the console script
    python -m biotsr ...  # the same thing, for a virtualenv without its bin on PATH
"""

from __future__ import annotations

from .errors import (
    EXIT_INTERNAL,
    EXIT_INCOMPLETE,
    EXIT_INTERRUPTED,
    EXIT_OK,
    EXIT_PARTIAL,
    EXIT_USAGE,
    CLIError,
)
from .main import main

__all__ = [
    "main",
    "CLIError",
    "EXIT_OK",
    "EXIT_INTERNAL",
    "EXIT_USAGE",
    "EXIT_PARTIAL",
    "EXIT_INCOMPLETE",
    "EXIT_INTERRUPTED",
]
