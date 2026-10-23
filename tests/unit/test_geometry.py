"""Geometry: theta, the fold, and the three degeneracy guards.

The degenerate cases below are enumerated from an audit of every degenerate
triangle shape (ADR-019).  The package implements the *representation* and
these tests; the downstream scientific policy
for a method that meets an invalid result is deliberately not decided here,
because the encoder only ever receives ``valid=True`` either way.
"""

from __future__ import annotations

import math

import pytest

from biotsr.domain.point import Point
from biotsr.geometry.conventions import DegeneracyReason, PI_CONVENTION
from biotsr.geometry.triangle import triangle_geometry


def P(x, y=0.0, z=0.0):
    return Point(float(x), float(y), float(z))


# ---------------------------------------------------------------------------
# Valid geometry
# ---------------------------------------------------------------------------


def test_known_geometry_right_isoceles():
    """A hand-checkable case.

    Base (0,0)-(2,0); midpoint (1,0); apex (1,1).
    s1 = d(l1, l3) = sqrt(2); s2 = base/2 = 1; s3 = median = 1.
    ratio = (2 - 1 - 1) / (2*1*1) = 0  ->  theta = 90 degrees exactly.
    """
    g = triangle_geometry(P(0, 0), P(2, 0), P(1, 1), PI_CONVENTION)
    assert g.valid
    assert g.theta_raw == pytest.approx(90.0, abs=1e-12)
    assert g.theta == pytest.approx(90.0, abs=1e-12)
    assert g.max_distance == pytest.approx(2.0)


def test_fold_boundary_at_exactly_90_is_not_folded():
    g = triangle_geometry(P(0, 0), P(2, 0), P(1, 1), PI_CONVENTION)
    assert g.theta_raw <= 90
    assert g.theta == g.theta_raw, "theta_raw == 90 must pass through unfolded"


def test_theta_raw_above_90_folds_to_the_supplement():
    # Apex pulled toward l1 so the median makes an obtuse angle at the midpoint.
    g = triangle_geometry(P(0, 0), P(10, 0), P(0.5, 1), PI_CONVENTION)
    assert g.valid
    assert g.theta_raw > 90
    assert g.theta == pytest.approx(abs(180 - g.theta_raw))
    assert 0 <= g.theta <= 90


def test_theta_raw_is_retained_alongside_the_folded_value():
    """No legacy triplet format stores the pre-fold angle; this one does."""
    g = triangle_geometry(P(0, 0), P(10, 0), P(0.5, 1), PI_CONVENTION)
    assert g.theta_raw is not None and g.theta is not None
    assert g.theta_raw != g.theta


def test_standard_convention_uses_math_pi():
    """The standard divisor is math.pi, so the pre-fold range is exactly [0, 180]."""
    assert PI_CONVENTION.radians_per_half_turn == math.pi
    assert PI_CONVENTION.degrees(math.pi) == pytest.approx(180.0, abs=1e-12)
    assert PI_CONVENTION.degrees(math.pi / 2) == pytest.approx(90.0, abs=1e-12)


def test_max_distance_is_not_necessarily_the_base_edge():
    """The trap worth restating: in cases 2-4 the base is chosen by LABEL order, so
    the longest edge may never enter the theta formula."""
    g = triangle_geometry(P(0, 0), P(1, 0), P(0, 20), PI_CONVENTION)
    base = P(0, 0).distance_to(P(1, 0))
    longest = P(1, 0).distance_to(P(0, 20))  # the l2-l3 edge, which theta never uses
    assert base == pytest.approx(1.0)
    assert g.max_distance == pytest.approx(longest)
    assert g.max_distance > base


# ---------------------------------------------------------------------------
# Near-collinear but VALID -- must not be rejected
# ---------------------------------------------------------------------------


def test_exactly_collinear_apex_outside_the_base_segment_is_valid():
    g = triangle_geometry(P(0, 0), P(2, 0), P(10, 0), PI_CONVENTION)
    assert g.valid, "a zero-area triangle is not automatically invalid"
    assert g.theta == pytest.approx(0.0, abs=1e-9)
    assert g.degeneracy_reason is None


def test_exactly_collinear_apex_between_the_base_endpoints_is_valid():
    """The apex lies strictly inside the base segment but not at its midpoint.

    The median points the other way, so theta_raw is 180 and folds to 0 -- valid,
    not an error.
    """
    g = triangle_geometry(P(0, 0), P(10, 0), P(3, 0), PI_CONVENTION)
    assert g.valid
    assert g.theta_raw == pytest.approx(180.0, abs=1e-9)
    assert g.theta == pytest.approx(0.0, abs=1e-9)


def test_near_collinear_at_realistic_magnitudes_is_valid_and_clamped_if_needed():
    """Realistic PDB coordinate magnitudes, an almost-collinear triple.

    This must produce a usable angle, not a rejection: 'nearly collinear' is a
    valid limiting geometry, distinct from 'mathematically degenerate'.
    """
    g = triangle_geometry(
        P(26.266, 25.413, 2.842),
        P(46.266, 25.413, 2.842),
        P(36.266, 25.413, 2.8420000000001),
        PI_CONVENTION,
    )
    assert g.valid
    assert g.degeneracy_reason is None
    assert 0 <= g.theta <= 90


# A real exactly-collinear triple at realistic PDB coordinate magnitudes whose
# law-of-cosines ratio lands at 1.0000000000072076 -- outside the acos domain by
# 7.2e-12, purely through floating-point rounding.  Found by an exhaustive search
# over collinear triples in the -60..60 A coordinate box; overshoots like this are
# common rather than exotic, because collinear triples are.
_OVERSHOOTING_COLLINEAR = (
    P(49.95199488290217, -2.1798641080232457, -53.58993872478085),
    P(52.80356517799265, 1.3133053241275765, -78.41028182741476),
    P(51.37790145986254, -0.4331306410922675, -66.00116720944837),
)
_KNOWN_EXCURSION = 7.207567875866516e-12


def test_acos_clamp_repairs_floating_point_excursion_and_is_counted():
    """A ratio outside [-1, 1] by rounding alone: clamped, counted, still valid."""
    g = triangle_geometry(*_OVERSHOOTING_COLLINEAR, convention=PI_CONVENTION)
    assert g.valid, "floating-point noise must not invalidate a real triangle"
    assert g.acos_clamped, "the clamp must fire, and must say so"
    assert g.degeneracy_reason is None
    assert g.theta == pytest.approx(0.0, abs=1e-9)


def test_clamp_is_refused_when_the_excursion_exceeds_the_tolerance():
    """Clamping must not disguise a genuinely invalid geometry.

    Same triple, but with a tolerance tighter than its known excursion: the result
    is reported as an acos domain error rather than silently pinned to the
    boundary.  This is the guard that keeps 'numerical safety' from quietly
    becoming 'accept anything'.
    """
    g = triangle_geometry(
        *_OVERSHOOTING_COLLINEAR,
        convention=PI_CONVENTION,
        acos_tolerance=_KNOWN_EXCURSION / 1000,
    )
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.ACOS_DOMAIN_ERROR
    assert g.theta is None


def test_clamp_disabled_reports_the_excursion_instead_of_repairing_it():
    g = triangle_geometry(*_OVERSHOOTING_COLLINEAR, convention=PI_CONVENTION, clamp_acos=False)
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.ACOS_DOMAIN_ERROR


def test_clamping_only_ever_affects_out_of_domain_input():
    """An ordinary triangle is bit-identical with and without the clamp.

    This is what makes clamping a numerical-safety measure rather than a
    scientific change.
    """
    ordinary = (P(0, 0), P(10, 0), P(3, 4))
    on = triangle_geometry(*ordinary, convention=PI_CONVENTION, clamp_acos=True)
    off = triangle_geometry(*ordinary, convention=PI_CONVENTION, clamp_acos=False)
    assert on.valid and off.valid
    assert on.theta == off.theta
    assert on.theta_raw == off.theta_raw
    assert not on.acos_clamped


# ---------------------------------------------------------------------------
# Mathematically degenerate -- must be reported, never forced to theta = 0
# ---------------------------------------------------------------------------


def test_zero_length_base_edge():
    g = triangle_geometry(P(1, 2, 3), P(1, 2, 3), P(4, 5, 6), PI_CONVENTION)
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.ZERO_BASE_EDGE
    assert g.theta is None and g.theta_raw is None
    assert g.max_distance > 0, "max_distance is still well defined"


def test_coincident_points_outside_the_base_pair():
    """The insidious case: without a pairwise check this completes with theta = 0
    and a nonzero max_distance, minting a plausible key for a zero-area triangle."""
    g = triangle_geometry(P(0, 0, 0), P(4, 5, 6), P(4, 5, 6), PI_CONVENTION)
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.COINCIDENT_POINTS
    assert g.theta is None


def test_coincident_l1_and_l3():
    g = triangle_geometry(P(4, 5, 6), P(0, 0, 0), P(4, 5, 6), PI_CONVENTION)
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.COINCIDENT_POINTS


def test_apex_exactly_at_the_base_midpoint():
    g = triangle_geometry(P(-2, 0, 0), P(2, 0, 0), P(0, 0, 0), PI_CONVENTION)
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.APEX_AT_BASE_MIDPOINT
    assert g.theta is None and g.theta_raw is None


def test_degenerate_results_never_carry_a_zero_theta():
    """Explicitly: an undefined angle is None, never 0.0.

    Forcing it to zero would mint a normal theta bin (class 1) and a normal key
    for a triangle that has no angle at all.
    """
    degenerate = [
        triangle_geometry(P(1, 1, 1), P(1, 1, 1), P(2, 2, 2), PI_CONVENTION),
        triangle_geometry(P(0, 0, 0), P(3, 3, 3), P(3, 3, 3), PI_CONVENTION),
        triangle_geometry(P(-1, 0, 0), P(1, 0, 0), P(0, 0, 0), PI_CONVENTION),
    ]
    for g in degenerate:
        assert not g.valid
        assert g.theta is None
        assert g.theta_raw is None
        assert g.degeneracy_reason is not None


# ===========================================================================
# Defined linear geometry vs undefined geometry
#
# The settled policy (ADR-033):
#
#   * a triangle whose denominator is valid and whose theta_raw is exactly 0 or
#     exactly 180 is VALID and is encoded;
#   * both fold to theta 0, so they share a theta class;
#   * mirror keys off theta_raw, so raw 0 and raw 180 have DIFFERENT orientation;
#   * theta that cannot be evaluated at all (d3 = 0, coincident points, any zero
#     denominator) is UNDEFINED: no key, triplet skipped, reason recorded;
#   * theta is never forced to zero to paper over an undefined case.
# ===========================================================================

from biotsr.binning.theta import THETA_SCHEME_29  # noqa: E402
from biotsr.geometry.conventions import SQRT_MATH_SQRT  # noqa: E402


def _collinear_raw_zero():
    """Apex on the far side of the base midpoint from p0, on the base line.

    theta is measured between the median (midpoint -> apex) and the half-base
    (midpoint -> p0), so putting the apex on the opposite ray gives an angle of
    180 there, and theta_raw = 180 - 180 = 0.
    """
    p0 = Point(0.0, 0.0, 0.0)
    p1 = Point(4.0, 0.0, 0.0)          # midpoint at x = 2
    p2 = Point(9.0, 0.0, 0.0)          # beyond p1, same line
    return p0, p1, p2


def _collinear_raw_180():
    """Apex on the same ray as p0 from the midpoint -- theta_raw = 180."""
    p0 = Point(0.0, 0.0, 0.0)
    p1 = Point(4.0, 0.0, 0.0)          # midpoint at x = 2
    p2 = Point(-5.0, 0.0, 0.0)         # beyond p0, same line
    return p0, p1, p2


def test_a_valid_collinear_triangle_with_theta_raw_zero_is_encoded():
    g = triangle_geometry(*_collinear_raw_zero(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    assert g.valid
    assert g.degeneracy_reason is None
    assert g.theta_raw == pytest.approx(0.0, abs=1e-9)
    assert g.theta == pytest.approx(0.0, abs=1e-9)


def test_a_valid_collinear_triangle_with_theta_raw_180_is_encoded():
    g = triangle_geometry(*_collinear_raw_180(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    assert g.valid
    assert g.degeneracy_reason is None
    assert g.theta_raw == pytest.approx(180.0, abs=1e-9)


def test_raw_zero_and_raw_180_share_the_folded_theta_and_the_theta_class():
    zero = triangle_geometry(*_collinear_raw_zero(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    one80 = triangle_geometry(*_collinear_raw_180(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    assert zero.theta == pytest.approx(one80.theta, abs=1e-9) == pytest.approx(0.0, abs=1e-9)
    assert THETA_SCHEME_29.classify(zero.theta) == THETA_SCHEME_29.classify(one80.theta) == 1


def test_raw_zero_and_raw_180_get_DIFFERENT_mirror_orientation():
    """Same folded theta class, opposite orientation -- the point of keeping raw."""
    from biotsr.modifiers.mirror import MIRROR_DISABLED, MIRROR_NEGATE_WHEN_PREFOLD_GT_90

    zero = triangle_geometry(*_collinear_raw_zero(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    one80 = triangle_geometry(*_collinear_raw_180(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)

    mirror = MIRROR_NEGATE_WHEN_PREFOLD_GT_90
    assert mirror.apply(500, zero.theta_raw) == 500      # raw <= 90 -> regular
    assert mirror.apply(500, one80.theta_raw) == -500    # raw > 90  -> mirrored
    assert mirror.category_for(zero.theta_raw) != mirror.category_for(one80.theta_raw)
    # ...and with mirror off, neither is ever negated.
    assert MIRROR_DISABLED.apply(500, zero.theta_raw) == 500
    assert MIRROR_DISABLED.apply(500, one80.theta_raw) == 500


def test_d3_zero_is_undefined_and_is_skipped_with_a_reason():
    """Apex exactly on the base midpoint: the median has zero length."""
    g = triangle_geometry(
        Point(0.0, 0.0, 0.0), Point(4.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
        PI_CONVENTION, sqrt=SQRT_MATH_SQRT,
    )
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.APEX_AT_BASE_MIDPOINT
    # Not forced to zero -- undefined means absent, not zero.
    assert g.theta_raw is None
    assert g.theta is None
    # max_distance is still well defined and still reported.
    assert g.max_distance == pytest.approx(4.0)


def test_coincident_base_points_are_undefined_and_skipped():
    g = triangle_geometry(
        Point(1.0, 2.0, 3.0), Point(1.0, 2.0, 3.0), Point(4.0, 0.0, 0.0),
        PI_CONVENTION, sqrt=SQRT_MATH_SQRT,
    )
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.ZERO_BASE_EDGE
    assert g.theta is None


def test_a_coincident_apex_is_undefined_and_skipped():
    """p2 == p0 would otherwise complete with a plausible-looking theta of 180."""
    g = triangle_geometry(
        Point(0.0, 0.0, 0.0), Point(4.0, 0.0, 0.0), Point(0.0, 0.0, 0.0),
        PI_CONVENTION, sqrt=SQRT_MATH_SQRT,
    )
    assert not g.valid
    assert g.degeneracy_reason == DegeneracyReason.COINCIDENT_POINTS
    assert g.theta is None


def test_near_collinear_but_valid_geometry_is_encoded_normally():
    """A sliver triangle is a triangle: tiny area, valid denominator, real theta."""
    g = triangle_geometry(
        Point(0.0, 0.0, 0.0), Point(4.0, 0.0, 0.0), Point(9.0, 1e-7, 0.0),
        PI_CONVENTION, sqrt=SQRT_MATH_SQRT,
    )
    assert g.valid
    assert g.degeneracy_reason is None
    assert 0.0 <= g.theta_raw <= 180.0
    assert 0.0 <= g.theta <= 90.0


def test_clamping_only_ever_repairs_a_tiny_excursion():
    """Clamping happens AFTER denominator validity, and only for rounding noise."""
    g = triangle_geometry(*_collinear_raw_180(), PI_CONVENTION, sqrt=SQRT_MATH_SQRT)
    assert g.valid
    # The exactly-collinear case is where an acos argument can round past -1.
    assert g.theta_raw == pytest.approx(180.0, abs=1e-9)
    # A genuinely out-of-domain ratio is refused rather than clamped: the guard
    # is the tolerance, and it is small.
    from biotsr.geometry.conventions import DEFAULT_ACOS_TOLERANCE

    assert DEFAULT_ACOS_TOLERANCE == 1e-9


def test_the_engine_encodes_valid_linear_triples_and_skips_undefined_ones():
    """End to end through the shared core, with the ledger balancing."""
    from biotsr.domain.point import LabeledPoint, PointSet
    from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
    from biotsr.methods.engine import run_core
    from biotsr.profiles.standard import STANDARD_PROTEIN_CALPHA_V1

    sid = StructureId("SYNTH")

    def point(x, name, seq, label):
        return LabeledPoint(
            identity=AtomIdentity(
                structure=sid, model=1, chain=ChainId("A"),
                residue_name=name, residue_seq=seq, atom_name="CA",
            ),
            coord=Point(x, 0.0, 0.0),
            label=label,
        )

    # Four collinear Ca points: every triple is zero-area but every denominator
    # is nonzero, so every triple must be encoded.
    points = PointSet([
        point(0.0, "ALA", 1, 21), point(4.0, "GLY", 2, 11),
        point(9.0, "LEU", 3, 8), point(15.0, "SER", 4, 16),
    ])
    out = run_core(points, STANDARD_PROTEIN_CALPHA_V1)
    assert out.counts.n_candidates == 4
    assert out.counts.n_degenerate == 0
    assert out.counts.n_triangles == 4
    assert all(t.theta == pytest.approx(0.0, abs=1e-9) for t in out.triplets)

    # Now put one point exactly on another's midpoint-with-a-third.
    degenerate = PointSet([
        point(0.0, "ALA", 1, 21), point(4.0, "GLY", 2, 11), point(2.0, "LEU", 3, 8),
    ])
    out = run_core(degenerate, STANDARD_PROTEIN_CALPHA_V1)
    assert out.counts.n_candidates == 1
    assert out.counts.n_triangles == 0
    assert out.counts.n_degenerate == 1
    assert out.skipped[0].reason == DegeneracyReason.APEX_AT_BASE_MIDPOINT
