"""The size-filter modifier: a predicate on ``max_distance``, hooked at
``CANDIDATE_FILTER``.

Boundary semantics were **characterized by execution**, not read off a comment
(``validation/scripts/legacy_size_filter_characterization.py``, evidence in
``validation/evidence/size_filter/``).  With the cutoff placed exactly on a
``max_distance`` value that really occurs in the 1UBQ golden, the legacy
implementation keeps 14,174 triangles and excludes every one of the 26 whose
``max_distance`` equals the cutoff; nudging the cutoff up by one ulp admits
exactly those 26.  The criterion is therefore **strict**::

    keep the candidate  <=>  max_distance < cutoff

The filter reads nothing but ``max_distance``, which is what makes it reusable
unchanged by the amino-acid atomic, drug, carbohydrate and nucleotide families --
there is no protein-specific logic here and none should be added.  All four now
use it, and none of them has a *historical* size filter of its own: the standard
package offers the shared capability with this boundary rule stated explicitly
rather than inherited by assumption (ADR-028, ADR-035).

**Optional and OFF BY DEFAULT, for every method** (ADR-035).  No profile carries
a filter unless one was asked for, and no family switches one on implicitly.  It
does not enter ``KeySpace`` -- it changes which triangles were drawn, not what a
key means -- and does enter ``SamplingProfile``, so two different cutoffs are
correctly refused as non-comparable frequency vectors while their keys stay
comparable.

**A rejected candidate is not a degenerate one.**  Legacy conflated them by
counting neither, dropping 56,126 of 70,300 triangles with no trace in any
output.  Here the two are separate counters with separate reasons, because they
mean different things: a degenerate triple has no definable angle, while a
size-filtered triple is perfectly well-defined and merely outside the range the
caller asked about.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .base import ModifierHook, ModifierName

__all__ = ["HOOK", "NAME", "Comparison", "SizeFilter", "LEGACY_DEFAULT_CUTOFF"]

NAME = ModifierName.SIZE_FILTER
HOOK = ModifierHook.CANDIDATE_FILTER

#: The legacy parameter's default. Not a real cutoff -- it is a sentinel large
#: enough that nothing is ever filtered, which is why legacy runs "without a size
#: filter" still pass a number.
LEGACY_DEFAULT_CUTOFF = 10000


class Comparison:
    """How the cutoff is compared against ``max_distance``."""

    #: The characterized legacy semantics: ``max_distance < cutoff``.
    STRICT_LESS_THAN = "strict_less_than"

    #: Not used by any legacy implementation; named so that a future method that
    #: genuinely needs it can say so explicitly rather than silently differing.
    LESS_THAN_OR_EQUAL = "less_than_or_equal"


@dataclass(frozen=True)
class SizeFilter:
    """Keep only triangles whose largest edge is below ``cutoff`` angstroms."""

    cutoff: float
    comparison: str = Comparison.STRICT_LESS_THAN

    def __post_init__(self) -> None:
        if not isinstance(self.cutoff, (int, float)) or isinstance(self.cutoff, bool):
            raise ValueError("size_filter cutoff must be a number, got %r" % (self.cutoff,))
        if self.cutoff <= 0:
            raise ValueError(
                "size_filter cutoff must be positive; %r would reject every "
                "triangle" % (self.cutoff,)
            )
        if self.comparison not in (
            Comparison.STRICT_LESS_THAN,
            Comparison.LESS_THAN_OR_EQUAL,
        ):
            raise ValueError("unknown size-filter comparison %r" % (self.comparison,))

    def accepts(self, max_distance: float) -> bool:
        if self.comparison == Comparison.STRICT_LESS_THAN:
            return max_distance < self.cutoff
        return max_distance <= self.cutoff

    def describe(self) -> str:
        symbol = "<" if self.comparison == Comparison.STRICT_LESS_THAN else "<="
        return "max_distance %s %r" % (symbol, self.cutoff)

    @classmethod
    def from_value(cls, value: Optional[float]) -> Optional["SizeFilter"]:
        """Build one from a user-supplied value, treating ``None`` as no filter."""
        if value is None:
            return None
        if isinstance(value, SizeFilter):
            return value
        return cls(cutoff=float(value))
