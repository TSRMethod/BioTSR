"""The candidate triangle value object.

A candidate is three *indices into a PointSet*, never coordinates: keeping this
layer index-only is what lets a future vectorized generator replace
enumeration without touching geometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

__all__ = ["TriangleCandidate"]


@dataclass(frozen=True, order=True)
class TriangleCandidate:
    """Three distinct point indices, always stored ascending (``i < j < k``).

    ``presentation`` is the order the three points are handed to canonical
    ordering, when that is not simply ``(i, j, k)``.  It exists because one
    ordering convention is position-sensitive: ``TSR-Package``'s two cross
    programs spell their all-distinct test as the chained ``l1 != l2 != l3``,
    which also fires when the **first and third** labels are equal, and then
    resolves by first position.  Those programs present a cross triple as
    ``(pair member, pair member, lone point)``, so the same three atoms reach the
    ordering rule in a different arrangement than ascending index order would
    give, and reach a different branch.

    It is ``None`` for every other convention and for every non-cross method,
    where ordering is presentation-invariant for all but exact-tie inputs.
    """

    i: int
    j: int
    k: int
    presentation: Optional[Tuple[int, int, int]] = None

    def __post_init__(self) -> None:
        if not (isinstance(self.i, int) and isinstance(self.j, int) and isinstance(self.k, int)):
            raise ValueError("candidate indices must be ints")
        if not (self.i < self.j < self.k):
            raise ValueError(
                "candidate indices must be strictly ascending: got (%r, %r, %r)"
                % (self.i, self.j, self.k)
            )
        if self.presentation is not None:
            if sorted(self.presentation) != [self.i, self.j, self.k]:
                raise ValueError(
                    "presentation %r must be a permutation of (%r, %r, %r)"
                    % (self.presentation, self.i, self.j, self.k)
                )

    def as_tuple(self) -> Tuple[int, int, int]:
        return (self.i, self.j, self.k)

    def vertices(self) -> Tuple[int, int, int]:
        """The three indices in the order the ordering rule should see them."""
        return self.presentation or (self.i, self.j, self.k)
