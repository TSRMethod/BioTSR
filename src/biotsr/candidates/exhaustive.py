"""Exhaustive candidate generation: all ``C(n, 3)`` index triples.

Direct, deterministic nested enumeration -- ``for i: for j>i: for k>j`` -- which is
both the legacy behaviour and the required one.  No set is built: order is fixed
by construction, so candidate order, triplet order and serialization order can
never depend on hash traversal (the direct fix for the confirmed nondeterminism
in the legacy frequency-table builder).

Yielding index triples rather than coordinates keeps this layer independent of
geometry, which is what makes a future vectorized generator a
drop-in replacement.
"""

from __future__ import annotations

from typing import Iterator

from ..domain.candidate import TriangleCandidate

__all__ = ["STRATEGY_NAME", "exhaustive", "count"]

STRATEGY_NAME = "exhaustive_c_n_3"


def exhaustive(n_points: int) -> Iterator[TriangleCandidate]:
    """Yield every ``i < j < k`` over ``range(n_points)``, in lexicographic order."""
    if n_points < 0:
        raise ValueError("n_points must be non-negative")
    for i in range(0, n_points - 2):
        for j in range(i + 1, n_points - 1):
            for k in range(j + 1, n_points):
                yield TriangleCandidate(i, j, k)


def count(n_points: int) -> int:
    """``C(n, 3)`` without enumerating."""
    if n_points < 3:
        return 0
    return n_points * (n_points - 1) * (n_points - 2) // 6
