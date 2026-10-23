"""Mixed-side candidate generation -- the one cross-specific strategy.

**Mixed-side cross candidate** (the canonical term, used everywhere in code and
documentation): a triple of points that takes at least one point from *each* of
the two sides.  With pools ``A`` and ``B`` the only allowed compositions are::

    A + A + B
    A + B + B

and the two monochromatic compositions ``A + A + A`` and ``B + B + B`` are not
cross triangles at all.  A cross key describes an *interaction*, so a triangle
entirely inside one side describes nothing about the other.

Counting.  With ``a = |A|`` and ``b = |B|``::

    a * C(b, 2)  +  b * C(a, 2)     ==     C(a+b, 3) - C(a, 3) - C(b, 3)

Both spellings are exact and a test asserts the identity over a grid of pool
sizes.  Every historical cross implementation emitted exactly this many
triangles -- ``Cross_TSR.py`` on 3EJ1 A/B wrote ``175 * C(169,2) + 169 * C(175,2)
= 5057325`` lines.

**This is a generator, not a filter.**  The historical ``cross_key.py`` walks
all ``C(a+b, 3)`` triples and ``continue``s on the monochromatic ones; that is
correct but wasteful, and it makes "did a monochromatic triple leak through?" a
property of a skipped branch rather than of the enumeration.  Here the
monochromatic triples are never proposed: for a same-side pair ``(i, j)`` the
third index is drawn only from the opposite side's index list, by binary search.
So no candidate is generated and discarded, and ``count()`` is exact without
enumerating.

The order is lexicographic ``i < j < k`` over the pooled point set
(``mixed_side``).  Order fixes the triplet file's layout, never which triangles
exist.  (Two historical programs enumerated the same set in two blocks; that
order is not reproduced.)
"""

from __future__ import annotations

from bisect import bisect_right
from typing import Iterator, Sequence, Tuple

from ..domain.candidate import TriangleCandidate

__all__ = [
    "STRATEGY_MIXED_SIDE",
    "COMPOSITION_RULE",
    "mixed_side",
    "count",
    "count_from_sides",
    "side_index_lists",
]

STRATEGY_MIXED_SIDE = "mixed_side_lexicographic"

#: The one canonical name for the rule, carried in ``SamplingProfile``.
COMPOSITION_RULE = "mixed_side_AAB_ABB"


def _c2(n: int) -> int:
    return n * (n - 1) // 2 if n >= 2 else 0


def _c3(n: int) -> int:
    return n * (n - 1) * (n - 2) // 6 if n >= 3 else 0


def count(a: int, b: int) -> int:
    """The exact number of mixed-side candidates for pools of size ``a`` and ``b``."""
    if a < 0 or b < 0:
        raise ValueError("pool sizes must be non-negative")
    return a * _c2(b) + b * _c2(a)


def count_via_complement(a: int, b: int) -> int:
    """The same number by the other spelling: ``C(a+b,3) - C(a,3) - C(b,3)``."""
    if a < 0 or b < 0:
        raise ValueError("pool sizes must be non-negative")
    return _c3(a + b) - _c3(a) - _c3(b)


def side_index_lists(sides: Sequence[str]) -> Tuple[str, Tuple[int, ...], str, Tuple[int, ...]]:
    """Split positions by side tag.  Exactly two distinct tags are required.

    Returns ``(tag_a, indices_a, tag_b, indices_b)`` with ``tag_a`` the tag of
    position 0, so the answer does not depend on tag spelling or on which side
    happens to sort first.
    """
    if not sides:
        raise ValueError("a cross candidate set needs a non-empty side tagging")
    distinct = []
    for tag in sides:
        if tag not in distinct:
            distinct.append(tag)
    if len(distinct) != 2:
        raise ValueError(
            "mixed-side candidates need exactly two sides; got %d (%s)"
            % (len(distinct), ", ".join(repr(t) for t in distinct) or "none")
        )
    tag_a, tag_b = distinct
    idx_a = tuple(i for i, t in enumerate(sides) if t == tag_a)
    idx_b = tuple(i for i, t in enumerate(sides) if t == tag_b)
    return tag_a, idx_a, tag_b, idx_b


def count_from_sides(sides: Sequence[str]) -> int:
    """``count`` for an already-tagged point set; ``0`` when a side is empty."""
    distinct = set(sides)
    if len(distinct) < 2:
        return 0
    _, idx_a, _, idx_b = side_index_lists(sides)
    return count(len(idx_a), len(idx_b))


def mixed_side(sides: Sequence[str]) -> Iterator[TriangleCandidate]:
    """Yield every mixed-side triple in lexicographic ``i < j < k`` order.

    Monochromatic triples are never proposed: when ``i`` and ``j`` share a side,
    ``k`` is drawn only from the opposite side's index list.
    """
    if len(set(sides)) < 2:
        return
    tag_a, idx_a, tag_b, idx_b = side_index_lists(sides)
    opposite = {tag_a: idx_b, tag_b: idx_a}

    n = len(sides)
    for i in range(0, n - 2):
        for j in range(i + 1, n - 1):
            if sides[i] == sides[j]:
                others = opposite[sides[i]]
                start = bisect_right(others, j)
                for k in others[start:]:
                    yield TriangleCandidate(i, j, k)
            else:
                # `i` and `j` already straddle the two sides, so every later
                # index completes a mixed-side triple whatever its tag.
                for k in range(j + 1, n):
                    yield TriangleCandidate(i, j, k)
