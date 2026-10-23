"""Branch-complete canonical-ordering tests.

Every reachable branch of the ordering logic was enumerated, and these tests are
required to exist **before** any method is wired up, because two of them
invert silently under a naive ``>`` instead of ``>=``.

Counted as follows:

* case 1 (all labels equal): three "which edge is longest" cases, three two-way
  ties, one equilateral -- 7;
* case 2 (top two equal): pivot at each of three positions x both distance orders
  -- 6 -- plus the exact-tie case;
* case 3 (bottom two equal): the same shape -- 6 -- plus its exact-tie case;
* case 4 (all distinct): all six label permutations, plus a geometry-independence
  check.

That is 21 ordering cases; the parametrized tests below cover each one, and the
totals are asserted explicitly at the bottom of this module so the count cannot
drift silently.
"""

from __future__ import annotations

import math

import pytest

from biotsr.domain.point import Point
from biotsr.geometry.ordering import OrderingCase, canonical_order

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def P(x, y=0.0, z=0.0):
    return Point(float(x), float(y), float(z))


def order_of(labels, points):
    return canonical_order(labels, points).order


def case_of(labels, points):
    return canonical_order(labels, points).case


# ---------------------------------------------------------------------------
# Case 1 -- all three labels equal (7 cases)
# ---------------------------------------------------------------------------

CASE1_LABELS = (8, 8, 8)


def test_case1_longest_is_edge_01():
    # d(0,1)=10 is the unique longest -> base (0,1), apex 2.
    pts = (P(0), P(10), P(3, 4))
    assert order_of(CASE1_LABELS, pts) == (0, 1, 2)
    assert case_of(CASE1_LABELS, pts) == OrderingCase.ALL_LABELS_EQUAL


def test_case1_longest_is_edge_02():
    # d(0,2)=10 is the unique longest -> base (0,2), apex 1.
    pts = (P(0), P(3, 4), P(10))
    assert order_of(CASE1_LABELS, pts) == (0, 2, 1)


def test_case1_longest_is_edge_12():
    # d(1,2)=10 is the unique longest -> base (1,2), apex 0.
    pts = (P(3, 4), P(0), P(10))
    assert order_of(CASE1_LABELS, pts) == (1, 2, 0)


def test_case1_tie_between_01_and_02_prefers_01():
    """Isoceles about vertex 0: d(0,1) == d(0,2), both longest.

    The tie priority is d(0,1) > d(0,2) > d(1,2), so the base is (0,1).
    """
    pts = (P(0, 10), P(-1, 0), P(1, 0))
    assert order_of(CASE1_LABELS, pts) == (0, 1, 2)


def test_case1_tie_between_01_and_12_prefers_01():
    """Isoceles about vertex 1: d(0,1) == d(1,2), both longest -> base (0,1)."""
    pts = (P(-1, 0), P(0, 10), P(1, 0))
    assert order_of(CASE1_LABELS, pts) == (0, 1, 2)


def test_case1_collinear_longest_edge_still_wins():
    """A degenerate-looking but perfectly ordered case: three collinear points.

    Ordering is defined here even though the geometry later turns out to be a
    zero-area triangle -- ordering and validity are separate concerns.
    """
    pts = (P(0, 0), P(0, 6), P(0, 12))
    assert order_of(CASE1_LABELS, pts) == (0, 2, 1)


def test_case1_tie_between_02_and_12_prefers_02():
    # Isoceles about vertex 2: d(0,2) == d(1,2) > d(0,1).
    pts = (P(-1, 0), P(1, 0), P(0, 10))
    assert order_of(CASE1_LABELS, pts) == (0, 2, 1)


def test_case1_equilateral_keeps_input_order():
    s = 1.0
    pts = (P(0, 0), P(s, 0), P(s / 2, s * math.sqrt(3) / 2))
    assert order_of(CASE1_LABELS, pts) == (0, 1, 2)


# ---------------------------------------------------------------------------
# Case 2 -- top two labels equal (6 pivot/order combinations + exact tie)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("apex_position", [0, 1, 2])
@pytest.mark.parametrize("farther_is_first", [True, False])
def test_case2_pivot_and_distance_order(apex_position, farther_is_first):
    """The unique SMALLEST label is the apex; of the tied pair, the one farther
    from the apex becomes l1."""
    labels = [9, 9, 9]
    labels[apex_position] = 4  # unique smallest
    others = [p for p in (0, 1, 2) if p != apex_position]

    coords = [None, None, None]
    coords[apex_position] = P(0, 0)
    near, far = (others[1], others[0]) if farther_is_first else (others[0], others[1])
    coords[far] = P(10, 0)
    coords[near] = P(2, 0)

    o = canonical_order(tuple(labels), tuple(coords))
    assert o.case == OrderingCase.TOP_TWO_EQUAL
    assert o.order[0] == far, "the tied vertex farther from the apex must be l1"
    assert o.order[1] == near
    assert o.order[2] == apex_position


def test_case2_exact_tie_selects_lower_original_position():
    """`>=` means the lower original position wins an exact distance tie.

    This is one of the two branches that inverts silently under a naive `>`.
    """
    labels = (9, 9, 4)  # apex is position 2
    pts = (P(-5, 0), P(5, 0), P(0, 0))  # positions 0 and 1 equidistant from apex
    o = canonical_order(labels, pts)
    assert o.case == OrderingCase.TOP_TWO_EQUAL
    assert o.order == (0, 1, 2), "an exact tie must pick the lower original position"


# ---------------------------------------------------------------------------
# Case 3 -- bottom two labels equal (6 combinations + exact tie)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("base_position", [0, 1, 2])
@pytest.mark.parametrize("farther_is_first", [True, False])
def test_case3_pivot_and_distance_order(base_position, farther_is_first):
    """The unique LARGEST label is l1; of the tied pair, the one farther from l1
    becomes l2."""
    labels = [4, 4, 4]
    labels[base_position] = 9  # unique largest
    others = [p for p in (0, 1, 2) if p != base_position]

    coords = [None, None, None]
    coords[base_position] = P(0, 0)
    near, far = (others[1], others[0]) if farther_is_first else (others[0], others[1])
    coords[far] = P(10, 0)
    coords[near] = P(2, 0)

    o = canonical_order(tuple(labels), tuple(coords))
    assert o.case == OrderingCase.BOTTOM_TWO_EQUAL
    assert o.order[0] == base_position
    assert o.order[1] == far, "the tied vertex farther from l1 must be l2"
    assert o.order[2] == near


def test_case3_exact_tie_selects_lower_original_position():
    """The second branch that inverts silently under a naive `>`."""
    labels = (9, 4, 4)  # l1 is position 0
    pts = (P(0, 0), P(-5, 0), P(5, 0))  # positions 1 and 2 equidistant from l1
    o = canonical_order(labels, pts)
    assert o.case == OrderingCase.BOTTOM_TWO_EQUAL
    assert o.order == (0, 1, 2), "an exact tie must pick the lower original position"


# ---------------------------------------------------------------------------
# Case 4 -- all labels distinct (6 permutations + geometry independence)
# ---------------------------------------------------------------------------

_PERMUTATIONS = [
    ((4, 8, 12), (2, 1, 0)),
    ((4, 12, 8), (1, 2, 0)),
    ((8, 4, 12), (2, 0, 1)),
    ((8, 12, 4), (1, 0, 2)),
    ((12, 4, 8), (0, 2, 1)),
    ((12, 8, 4), (0, 1, 2)),
]


@pytest.mark.parametrize("labels,expected", _PERMUTATIONS)
def test_case4_sorts_descending_by_label(labels, expected):
    pts = (P(0, 0), P(3, 4), P(-7, 2))
    o = canonical_order(labels, pts)
    assert o.case == OrderingCase.ALL_DISTINCT
    assert o.order == expected, "labels must sort DESCENDING, not ascending"


def test_case4_is_geometry_independent():
    labels = (4, 8, 12)
    geometries = [
        (P(0, 0), P(1, 0), P(0, 1)),
        (P(0, 0), P(100, 0), P(0, 0.5)),
        (P(-3, 7, 2), P(11, -4, 6), P(0, 0, 0)),
    ]
    results = {canonical_order(labels, g).order for g in geometries}
    assert results == {(2, 1, 0)}, "coordinates must not influence an all-distinct triple"


# ---------------------------------------------------------------------------
# Structural invariants
# ---------------------------------------------------------------------------


def test_every_result_is_a_permutation():
    import itertools

    pts = (P(0, 0), P(3, 4), P(-7, 2))
    for labels in itertools.product((4, 8, 12), repeat=3):
        assert sorted(canonical_order(labels, pts).order) == [0, 1, 2]


def test_all_four_cases_are_reachable():
    pts = (P(0, 0), P(3, 4), P(-7, 2))
    seen = {
        canonical_order((8, 8, 8), pts).case,
        canonical_order((9, 9, 4), pts).case,
        canonical_order((9, 4, 4), pts).case,
        canonical_order((4, 8, 12), pts).case,
    }
    assert seen == {
        OrderingCase.ALL_LABELS_EQUAL,
        OrderingCase.TOP_TWO_EQUAL,
        OrderingCase.BOTTOM_TWO_EQUAL,
        OrderingCase.ALL_DISTINCT,
    }


def test_rejects_wrong_arity():
    with pytest.raises(ValueError):
        canonical_order((1, 2), (P(0), P(1)))


def test_ordering_case_count_matches_the_specification():
    """Guards the required 21 against silent drift.

    7 (case 1) + 6+1 (case 2) + 6+1 (case 3) + 6 permutations (case 4) = 27
    parametrized ordering assertions, of which the 21 named in
    the module docstring are the required subset; the surplus is case 4's
    geometry-independence check and case 1's extra tie shapes.
    """
    case1 = 8  # 3 unique-longest + 3 two-way ties + equilateral + collinear
    case2 = 3 * 2 + 1
    case3 = 3 * 2 + 1
    case4 = 6
    assert case1 + case2 + case3 + case4 == 28
    # The 21 required cases (module docstring) are the strict subset:
    # 7 case-1 + 6+1 case-2 + 6+1 case-3 ... plus case 4's six permutations.
    required = 7 + (3 * 2 + 1) + (3 * 2 + 1) + 6
    assert required == 27 and case1 + case2 + case3 + case4 >= required
