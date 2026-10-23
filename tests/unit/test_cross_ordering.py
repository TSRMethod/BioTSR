"""The standard canonical ordering, on the triples where historical cross code differed.

A characterization of ``TSR-Package``'s two CROSS programs found three genuine
disagreements with the standard four-case ordering, all read from source:
string collation of labels, a chained ``l1 != l2 != l3`` test, and a strict
``>`` tie comparison.  BioTSR implements only the standard ordering.  The tests
below pin the STANDARD answer on exactly the triples that separated the two, and
record the historical answer as a documented constant (never computed), so a
future refactor cannot quietly drift towards the historical behaviour.
"""

from __future__ import annotations

import pytest

from biotsr.domain.point import Point
from biotsr.geometry.conventions import (
    MEDIAN_FROM_MIDPOINT,
    PI_CONVENTION,
    SQRT_MATH_SQRT,
)
from biotsr.geometry.ordering import (
    STANDARD_ORDERING,
    OrderingCase,
    canonical_order,
)
from biotsr.geometry.triangle import triangle_geometry


def P(x, y, z=0.0):
    return Point(x, y, z)


# ===========================================================================
# 1 -- labels are compared as integers
# ===========================================================================


def test_the_standard_convention_orders_labels_numerically():
    """93 < 106, so label 108 is l1 and 93 is last.

    Historical difference, recorded as a constant: `TSR-Package`'s cross
    programs compared labels as strings (``"93" > "108" > "106"``) and so gave
    ``HISTORICAL_STRING_COLLATED_LABELS`` = [93, 108, 106]; in frozen 6G9M
    output all 405 triangles containing ESW/CL1 (label 93) carried Label1 = 93.
    """
    HISTORICAL_STRING_COLLATED_LABELS = [93, 108, 106]  # documentation only
    order = canonical_order((106, 93, 108), (P(0, 0), P(3, 0), P(0, 4)),
                            SQRT_MATH_SQRT, STANDARD_ORDERING)
    labels = [(106, 93, 108)[i] for i in order.order]
    assert labels == [108, 106, 93]
    assert labels != HISTORICAL_STRING_COLLATED_LABELS


def test_the_ordering_case_follows_label_equality():
    for labels in ((93, 93, 106), (106, 93, 106), (93, 93, 93), (1, 20, 100)):
        result = canonical_order(labels, (P(0, 0), P(3, 0), P(0, 4)),
                                 SQRT_MATH_SQRT, STANDARD_ORDERING)
        distinct = len(set(labels))
        if distinct == 1:
            assert result.case == OrderingCase.ALL_LABELS_EQUAL
        elif distinct == 3:
            assert result.case == OrderingCase.ALL_DISTINCT
        else:
            assert result.case in (OrderingCase.TOP_TWO_EQUAL,
                                   OrderingCase.BOTTOM_TWO_EQUAL)


# ===========================================================================
# 2 -- first-equals-third triples are not "all distinct"
# ===========================================================================


def test_a_first_equals_third_triple_is_routed_to_the_distance_resolved_case():
    """The standard tests all three label pairs, so (5, 9, 5) is a
    bottom-two-equal triple resolved by distance.

    Historical difference, recorded as a constant: the chained test
    ``l1 != l2 != l3`` is ``(l1 != l2) and (l2 != l3)``, true for (5, 9, 5), so
    the historical cross programs treated it as all-distinct and resolved it by
    first position, giving ``HISTORICAL_CHAINED_ORDER`` = (1, 2, 0).
    """
    HISTORICAL_CHAINED_ORDER = (1, 2, 0)  # documentation only
    labels = (5, 9, 5)
    # Position 0 is FARTHER from the unique max than position 2 is, so the
    # distance-resolved branch puts position 0 second.
    points = (P(0, 0), P(10, 0), P(9, 1))
    standard = canonical_order(labels, points, SQRT_MATH_SQRT, STANDARD_ORDERING)
    assert standard.case == OrderingCase.BOTTOM_TWO_EQUAL
    assert standard.order == (1, 0, 2)
    assert standard.order != HISTORICAL_CHAINED_ORDER


def test_a_first_equals_third_triple_with_the_tied_labels_largest():
    """(9, 5, 9): the two 9s are the top pair and 5 (position 1) is the apex.
    The historical chained test called this all-distinct."""
    order = canonical_order((9, 5, 9), (P(0, 0), P(10, 0), P(0, 1)),
                            SQRT_MATH_SQRT, STANDARD_ORDERING)
    assert order.case == OrderingCase.TOP_TWO_EQUAL
    assert order.order[2] == 1
    assert sorted(order.order[:2]) == [0, 2]


def test_a_genuinely_all_distinct_triple():
    labels = (300, 100, 200)
    points = (P(0, 0), P(3, 0), P(0, 4))
    a = canonical_order(labels, points, SQRT_MATH_SQRT, STANDARD_ORDERING)
    assert a.case == OrderingCase.ALL_DISTINCT
    assert a.order == (0, 2, 1)


# ===========================================================================
# 3 -- an exact distance tie goes to the lower original position
# ===========================================================================


def test_an_exact_distance_tie_resolves_to_the_lower_position():
    """``>=`` gives the tie to the lower original position.  Cases 2 and 3 only.

    Historical difference, recorded as a constant: the historical cross programs
    compared with ``>`` and gave the tie to the higher position,
    ``HISTORICAL_STRICT_TIE_ORDER`` = (1, 0, 2).
    """
    HISTORICAL_STRICT_TIE_ORDER = (1, 0, 2)  # documentation only
    labels = (7, 7, 3)              # case 2: top two equal, apex is position 2
    points = (P(-1, 0), P(1, 0), P(0, 5))   # both tied points equidistant
    standard = canonical_order(labels, points, SQRT_MATH_SQRT, STANDARD_ORDERING)
    assert standard.case == OrderingCase.TOP_TWO_EQUAL
    assert standard.order == (0, 1, 2)
    assert standard.order != HISTORICAL_STRICT_TIE_ORDER


# ===========================================================================
# 4 -- swapping the case-1 base endpoints is key-neutral
# ===========================================================================


@pytest.mark.parametrize(
    "points",
    [
        (P(0, 0), P(1, 0), P(0.3, 9.0)),
        (P(0, 0), P(2, 1), P(-4, 7)),
        (P(1, 2, 3), P(4, 1, 0), P(-2, -3, 8)),
    ],
)
def test_swapping_the_case_one_base_endpoints_leaves_the_folded_theta_alone(points):
    """With all three labels equal, which base endpoint is ``l1`` swaps
    ``theta_raw`` for its supplement; the FOLDED theta -- and therefore the key --
    is identical under the standard geometry."""
    labels = (108, 108, 108)
    a = canonical_order(labels, points, SQRT_MATH_SQRT, STANDARD_ORDERING)
    assert a.case == OrderingCase.ALL_LABELS_EQUAL
    swapped = (a.order[1], a.order[0], a.order[2])

    def geometry(order):
        p = [points[i] for i in order]
        return triangle_geometry(p[0], p[1], p[2], PI_CONVENTION,
                                 sqrt=SQRT_MATH_SQRT,
                                 median=MEDIAN_FROM_MIDPOINT)

    ga, gb = geometry(a.order), geometry(swapped)
    assert ga.theta == pytest.approx(gb.theta, abs=1e-12)
    assert ga.max_distance == pytest.approx(gb.max_distance, abs=1e-12)
    assert ga.theta_raw == pytest.approx(180 - gb.theta_raw, abs=1e-9)


# ===========================================================================
# The median convention and the ordering version are key-space facts
# ===========================================================================


def test_the_median_convention_is_part_of_the_key_space():
    """Historical note: `TSR-Package`'s cross programs computed the median from
    the three side lengths (Apollonius); the standard uses the midpoint."""
    from biotsr.profiles.cross import STANDARD_ATOMIC_CROSS_V1

    assert STANDARD_ATOMIC_CROSS_V1.key_space().median_convention == "from_midpoint"


def test_the_standard_version_string_is_unchanged_by_the_cross_methods():
    """Standard keys from before the cross methods must stay comparable with today's."""
    from biotsr.geometry.ordering import CANONICAL_ORDERING_VERSION

    assert CANONICAL_ORDERING_VERSION == "legacy_four_case_v1"
    assert STANDARD_ORDERING.version == CANONICAL_ORDERING_VERSION


def test_the_standard_cross_profile_reads_double_precision_coordinates():
    """Historical note: Bio.PDB stores coordinates as float32, so the two
    `TSR-Package` cross programs worked from single-precision input.  The
    standard reads float64 and records it in the key space (ADR-053)."""
    from biotsr.profiles.cross import STANDARD_ATOMIC_CROSS_V1
    from biotsr.structure.backends import pdb_legacy

    assert (
        STANDARD_ATOMIC_CROSS_V1.coordinate_precision
        == pdb_legacy.COORD_PRECISION_DOUBLE
    )
    assert STANDARD_ATOMIC_CROSS_V1.key_space().numerical_precision == "float64"
