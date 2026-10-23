"""Canonical ordering -- the single shared component.

This logic was compared across ``TSR-Package/TSR.py``,
``TSR_AMINOACID_PACKAGE/analyzer.py``, ``CarbohydratesTSR/.../dual_mode.py`` and
``WebsiteTSR/key_triplets_protein.py`` decision by decision, and **every one
agrees**.  So this is one pure, stateless function
with no method flag (ADR-006): it receives three labels and three points and
returns a permutation.  It performs no selection, no candidate generation, no
geometry beyond the pairwise distances the tie-breaks need, and no encoding.

Two invariants a reimplementation must not "clean up":

1. The label->position map is keyed by **label**, so duplicate labels overwrite.
   That is safe only because every branch reads a label guaranteed unique within
   its triple; "fixing" the map into a multimap reads wrong positions.
2. The sort is **descending**.  An accidental ascending sort silently reverses
   every all-distinct triple.

The tie-breaks use ``>=``, not ``>``.  In cases 2 and 3 that makes the **lower
original position** win an exact distance tie; a naive ``>`` inverts those two
branches silently, which the legacy audit rated the highest-severity technical
risk of the whole rewrite.

**Two ordering conventions, still one function.**
-------------------------------------------------

The "every one agrees" survey covered the nine *non-cross* implementations.  The
two cross programs in ``TSR-Package`` -- ``Drug_TSR.py`` and ``Cross_TSR.py`` --
do **not** agree with the rest of the family.
Three differences, all read from source and all confirmed against frozen output:

* **String collation.**  ``load_atom_seq`` stores the lexicon's sequence numbers
  as *strings* (``atomSeq[row[2]] = row[1]``), and the ordering branch then calls
  ``max``/``min``/``sort`` on them.  So label 93 (CL) sorts *above* label 108 (N),
  because ``"93" > "108"``.  This is not a hypothesis: all 405 triangles
  containing 6G9M's ``ESW/CL1`` in the frozen ``Drug_TSR.py`` reference carry
  ``Label1 = 93`` with larger numeric labels below it.
* **Chained inequality.**  The all-distinct test is spelled ``l1 != l2 != l3``,
  which Python reads as ``(l1 != l2) and (l2 != l3)`` -- true whenever the
  *outer two* labels are equal and differ from the middle one.  That case then
  resolves by **first position**, not by distance.  Every other implementation in
  the family spells the test as three explicit pairwise comparisons.
* **Strict tie comparison.**  Cases 2 and 3 compare with ``>`` rather than ``>=``,
  inverting the branch on an exact distance tie.

Those are properties of that lineage, not of a family, so they are a *convention
value object* -- exactly like ``AngleConvention`` and ``SqrtConvention`` -- and
the function stays single.  ``OrderingConvention.version`` enters ``KeySpace``,
so keys produced under the two conventions can never be silently compared.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

from ..domain.point import Point
from .conventions import SQRT_POW_HALF, SqrtConvention

__all__ = [
    "OrderingCase",
    "CanonicalOrder",
    "OrderingConvention",
    "STANDARD_ORDERING",
    "canonical_order",
    "CANONICAL_ORDERING_VERSION",
]

#: Bumped only when the ordering *semantics* change.  Enters ``KeySpace``.
CANONICAL_ORDERING_VERSION = "legacy_four_case_v1"



class OrderingCase:
    """The four reachable branches, named."""

    ALL_LABELS_EQUAL = 1
    TOP_TWO_EQUAL = 2
    BOTTOM_TWO_EQUAL = 3
    ALL_DISTINCT = 4


@dataclass(frozen=True)
class OrderingConvention:
    """The four-case canonical ordering, named and versioned.

    There is one convention, :data:`STANDARD_ORDERING`; its ``version`` enters
    ``KeySpace`` so that a change to the ordering semantics would change every
    fingerprint rather than silently change keys.  (The historical TSR-Package
    cross programs spelled the ordering differently -- string collation of
    labels and a chained distinctness test; BioTSR does not reproduce that.)
    """

    name: str
    version: str

    def content_fingerprint(self) -> str:
        return self.version


#: The ordering every method uses.
STANDARD_ORDERING = OrderingConvention(
    name="legacy_four_case",
    version=CANONICAL_ORDERING_VERSION,
)


@dataclass(frozen=True)
class CanonicalOrder:
    """The result: a permutation of the input positions plus which case decided it.

    ``order[0]`` is the position of ``l1``, ``order[1]`` of ``l2``, ``order[2]``
    of ``l3``.  ``(l1, l2)`` is the *base* edge; ``l3`` is the apex.
    """

    order: Tuple[int, int, int]
    case: int


def _other_two(chosen: int) -> Tuple[int, int]:
    """The legacy ``indexFind`` helper: the two positions that are not ``chosen``,
    in ascending original order."""
    return tuple(p for p in (0, 1, 2) if p != chosen)  # type: ignore[return-value]


def canonical_order(
    labels: Sequence[int],
    points: Sequence[Point],
    sqrt: SqrtConvention = SQRT_POW_HALF,
    convention: OrderingConvention = STANDARD_ORDERING,
) -> CanonicalOrder:
    """Order one triple of labelled points into ``(l1, l2, l3)`` positions.

    ``labels`` and ``points`` are three-element sequences in the candidate's own
    order (the legacy ``i, j, k``).  The returned positions index back into them.

    ``sqrt`` must be the same square-root convention the geometry stage will use.
    Cases 1-3 break ties by comparing distances, and the two legacy spellings
    disagree by one ULP on some inputs, so a mismatched convention could in
    principle resolve a near-tie the other way.  Each legacy lineage uses one
    spelling throughout, and so does each profile here.

    ``convention`` must be :data:`STANDARD_ORDERING`.
    """
    if len(labels) != 3 or len(points) != 3:
        raise ValueError("canonical_order expects exactly three labels and three points")

    if convention is not STANDARD_ORDERING and convention.version != CANONICAL_ORDERING_VERSION:
        raise ValueError(
            "unsupported ordering convention %r; BioTSR implements %r"
            % (convention.version, CANONICAL_ORDERING_VERSION)
        )
    la, lb, lc = labels[0], labels[1], labels[2]

    # Legacy `keepLabelIndex`: keyed by label, so a duplicate label overwrites the
    # earlier position.  Preserved verbatim -- see invariant 1 in the module docstring.
    keep_label_index = {}
    keep_label_index[la] = 0
    keep_label_index[lb] = 1
    keep_label_index[lc] = 2

    # Invariant 2: descending, comparing the integer labels numerically.
    sorted_labels = sorted((la, lb, lc), reverse=True)

    def dist(p: int, q: int) -> float:
        return sqrt(points[p].squared_distance_to(points[q]))

    def wins(a: float, b: float) -> bool:
        """The cases 2/3 tie-break: ``>=`` gives an exact tie to the lower
        original position."""
        return a >= b

    all_distinct = la != lb and la != lc and lb != lc

    if sorted_labels[0] == sorted_labels[1] and sorted_labels[1] == sorted_labels[2]:
        # Case 1 -- all three labels equal.  The longest edge becomes the base;
        # the opposite vertex is the apex.  Ties resolve d(0,1) > d(0,2) > d(1,2).
        d01 = dist(0, 1)
        d02 = dist(0, 2)
        d12 = dist(1, 2)
        longest = max(d01, d02, d12)
        if d01 >= longest:
            order = (0, 1, 2)
        elif d02 >= longest:
            order = (0, 2, 1)
        else:
            order = (1, 2, 0)
        case = OrderingCase.ALL_LABELS_EQUAL

    elif all_distinct:
        # Case 4 -- all distinct.  Purely label-driven, geometry-independent.
        order = (
            keep_label_index[sorted_labels[0]],
            keep_label_index[sorted_labels[1]],
            keep_label_index[sorted_labels[2]],
        )
        case = OrderingCase.ALL_DISTINCT

    elif sorted_labels[0] == sorted_labels[1] and sorted_labels[1] != sorted_labels[2]:
        # Case 2 -- top two equal.  The unique smallest label is the apex (l3);
        # of the tied pair, the one FARTHER from the apex becomes l1.
        apex = keep_label_index[sorted_labels[2]]
        first, second = _other_two(apex)
        if wins(dist(first, apex), dist(second, apex)):
            order = (first, second, apex)
        else:
            order = (second, first, apex)
        case = OrderingCase.TOP_TWO_EQUAL

    else:
        # Case 3 -- bottom two equal.  The unique largest label is l1; of the tied
        # pair, the one FARTHER from l1 becomes l2.
        base = keep_label_index[sorted_labels[0]]
        first, second = _other_two(base)
        if wins(dist(base, first), dist(base, second)):
            order = (base, first, second)
        else:
            order = (base, second, first)
        case = OrderingCase.BOTTOM_TWO_EQUAL

    return CanonicalOrder(order=order, case=case)
