"""The package's exception hierarchy.

Five things a caller can distinguish, and no more (ADR-076).  The
taxonomy answers the only questions a program actually asks of a failure:

``InvalidRequest``
    You asked for something the package cannot do.  A bad chain, a cutoff of
    zero, a modifier the method does not support, a profile that does not belong
    to the method.  Fix the call.
``IncompatibleRepresentation``
    Both things are valid; they cannot be combined.  Two results whose integer
    keys mean different things, two frequency vectors drawn differently, two SSE
    category breakdowns named under different conventions.  Fix the *comparison*,
    not either input.
``MalformedResult``
    A serialized result on disk is corrupt, truncated, internally inconsistent,
    or assembled from files that do not belong together.  Fix the files.
``UnsupportedFormat``
    A serialized result is well-formed and this version of the package does not
    know how to read it -- a schema version from the future.  Upgrade.
``TSRError``
    The common base.  ``except TSRError`` means "the request was wrong"; anything
    outside this tree is a defect in the package.

**Every one of them also inherits ``ValueError``**, and that is deliberate rather
than untidy.  The package's refusal types were plain ``ValueError`` subclasses
before this hierarchy existed, the CLI's error boundary catches ``ValueError``,
and any user script that says ``except ValueError`` is correct code.
Re-parenting those types under a clean ``TSRError``-only base would silently stop
catching them.  The extra base costs one entry in an MRO and breaks nobody.

**What is deliberately NOT in the tree**: ``AttributeError``, ``TypeError`` from a
wrong-typed argument, ``AssertionError`` from an internal invariant, and every
other sign that the package itself is wrong.  That exclusion is what gives
``except TSRError`` its meaning, and it is what lets the CLI map "user error" to
exit 2 and "our bug" to exit 1 without maintaining a list of exception types of
its own.

This module lives in ``domain`` because ``domain`` is the package's dependency
sink: it imports nothing from ``biotsr``, so everything else may import it.  A
root-level ``biotsr/errors.py`` could not be imported by ``domain`` (which would
break the sink) nor by ``biotsr.analysis`` (whose architecture test allows exactly
one package-level dependency, and that dependency is ``domain``).
"""

from __future__ import annotations

__all__ = [
    "TSRError",
    "InvalidRequest",
    "IncompatibleRepresentation",
    "MalformedResult",
    "UnsupportedFormat",
]


class TSRError(ValueError):
    """Base for every refusal the package raises on purpose.

    Inherits ``ValueError`` so that code catching ``ValueError`` keeps working;
    see the module docstring.
    """


class InvalidRequest(TSRError):
    """The caller asked for something that cannot be done.

    Bad input, bad configuration, an argument the method does not take, a value
    outside its admissible range.  The message names what *is* accepted, because
    a refusal that does not is a worse error message than no refusal at all.
    """


class IncompatibleRepresentation(TSRError):
    """Two scientific representations cannot be combined.

    The parent of the four ``Incompatible*`` types.  Each input is individually
    valid; it is putting them together that is not defined.
    """


class MalformedResult(TSRError):
    """A serialized result is corrupt, incomplete, or internally inconsistent.

    Includes the case the format exists to make impossible: files from two
    different results placed side by side.  The message names the quantity that
    disagreed and both of its values.
    """


class UnsupportedFormat(TSRError):
    """A serialized result is well formed and this version cannot read it.

    Distinct from :class:`MalformedResult` because the remedy is different: a
    malformed file is broken, an unsupported one is merely newer.
    """
