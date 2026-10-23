"""Triangle geometry: distances, the median construction, theta, and validity.

The result type is ``GeometryResult``, which can represent geometric *failure*
explicitly rather than only ever an angle (ADR-019).  The encoder only ever
receives ``valid=True``.

A semantic trap worth restating: **the base edge is not generally the longest
edge.**  Only in canonical-ordering case 1 is it.  In cases 2-4 the base is
chosen by label ordering, so ``max_distance`` may be the one edge that never
enters the theta formula.  ``max_distance`` feeds the *distance* bin; the base
edge feeds *theta*.  The legacy variable names (``dist12Class``,
``dist12 = dist01``) actively mislead here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from ..domain.point import Point
from .conventions import (
    DEFAULT_ACOS_TOLERANCE,
    MEDIAN_FROM_MIDPOINT,
    SQRT_POW_HALF,
    AngleConvention,
    DegeneracyReason,
    MedianConvention,
    SqrtConvention,
)

__all__ = ["GeometryResult", "triangle_geometry"]


@dataclass(frozen=True)
class GeometryResult:
    """The outcome of computing one triangle's geometry.

    ``theta_raw`` is the pre-fold angle (may exceed 90 and, under a legacy
    divisor, may exceed 180).  ``theta`` is the folded angle in ``[0, 90]`` and
    is what the theta bin is computed from.  Both are retained because
    ``theta_raw`` is what the mirror convention keys its sign off, so the geometry
    never needs a second pass for it.

    ``max_distance`` is always populated -- it is well-defined even for a
    degenerate triple -- so a skipped candidate can still be reported usefully.

    ``median_length`` is ``s3``, the base-to-apex median the angle is computed
    from.  It is reported rather than discarded because the SSE triplet format is
    the one historical format that prints it (as ``height``), so reproducing that
    file must not mean recomputing a value the geometry already had.  Nothing in
    the key path reads it.
    """

    valid: bool
    theta_raw: Optional[float]
    theta: Optional[float]
    max_distance: float
    degeneracy_reason: Optional[str] = None
    acos_clamped: bool = False
    median_length: Optional[float] = None


def triangle_geometry(
    p0: Point,
    p1: Point,
    p2: Point,
    convention: AngleConvention,
    *,
    sqrt: SqrtConvention = SQRT_POW_HALF,
    clamp_acos: bool = True,
    acos_tolerance: float = DEFAULT_ACOS_TOLERANCE,
    median: MedianConvention = MEDIAN_FROM_MIDPOINT,
) -> GeometryResult:
    """Compute the geometry of one canonically ordered triangle.

    ``p0``/``p1`` are the base pair (``l1``, ``l2``); ``p2`` is the apex (``l3``).

    ``sqrt`` selects the family's square-root spelling (``** 0.5`` for Ca,
    ``math.sqrt`` for every atomic lineage).  It is not cosmetic: the two
    disagree by one ULP often enough to move a real triangle across a bin
    boundary.

    ``median`` names how ``s3`` is computed: from the base midpoint's
    coordinates, the one convention implemented.  ``s3`` is a denominator, so the
    convention is key-affecting and is recorded in ``KeySpace``.

    Three independent guards run before any ``acos`` call, because clamping alone
    catches none of them (ADR-019):

    1. a zero-length **base** edge (``s2 == 0`` -> zero denominator);
    2. any other coincident pair -- the insidious case, since it can otherwise
       complete with ``theta = 0`` and a nonzero ``max_distance``, producing a
       plausible-looking key for a zero-area triangle rather than raising;
    3. the apex sitting exactly on the base midpoint (``s3 == 0`` -> zero
       denominator).

    Only after those does the ``acos`` argument get clamped, and only for
    excursions small enough to be floating-point noise.  Clamping is numerical
    safety, not a scientific change: a mathematically valid triangle's ratio lies
    in ``[-1, 1]`` exactly, so a clamp can only ever repair rounding.  It is
    counted (``acos_clamped``) so the claim is auditable rather than assumed.
    """
    dist01 = sqrt(p0.squared_distance_to(p1))
    dist02 = sqrt(p0.squared_distance_to(p2))
    dist12 = sqrt(p1.squared_distance_to(p2))
    max_distance = max(dist01, dist02, dist12)

    # Guard 1 + 2: coincidence.  Checked pairwise and reported distinctly, since a
    # zero base edge and a coincident non-base pair fail for different reasons.
    if dist01 == 0.0:
        return GeometryResult(
            valid=False,
            theta_raw=None,
            theta=None,
            max_distance=max_distance,
            degeneracy_reason=DegeneracyReason.ZERO_BASE_EDGE,
        )
    if dist02 == 0.0 or dist12 == 0.0:
        return GeometryResult(
            valid=False,
            theta_raw=None,
            theta=None,
            max_distance=max_distance,
            degeneracy_reason=DegeneracyReason.COINCIDENT_POINTS,
        )

    s1 = dist02
    s2 = dist01 / 2
    if median is not MEDIAN_FROM_MIDPOINT and median.name != MEDIAN_FROM_MIDPOINT.name:
        raise ValueError(
            "unsupported median convention %r; BioTSR computes the median from "
            "the base midpoint (%r)" % (median.name, MEDIAN_FROM_MIDPOINT.name)
        )
    midpoint = p0.midpoint_to(p1)
    s3 = sqrt(
        (midpoint.x - p2.x) ** 2
        + (midpoint.y - p2.y) ** 2
        + (midpoint.z - p2.z) ** 2
    )

    # Guard 3: apex exactly at the base midpoint -- the median has zero length, so
    # the law-of-cosines denominator vanishes.
    if s3 == 0.0:
        return GeometryResult(
            valid=False,
            theta_raw=None,
            theta=None,
            max_distance=max_distance,
            degeneracy_reason=DegeneracyReason.APEX_AT_BASE_MIDPOINT,
        )

    ratio = (s1 ** 2 - s2 ** 2 - s3 ** 2) / (2 * s2 * s3)

    clamped = False
    if ratio < -1.0 or ratio > 1.0:
        excursion = abs(ratio) - 1.0
        if not clamp_acos or excursion > acos_tolerance:
            # Too far outside to be rounding: report it rather than disguise a
            # genuinely invalid geometry behind a clamp.
            return GeometryResult(
                valid=False,
                theta_raw=None,
                theta=None,
                max_distance=max_distance,
                degeneracy_reason=DegeneracyReason.ACOS_DOMAIN_ERROR,
            )
        ratio = 1.0 if ratio > 1.0 else -1.0
        clamped = True

    theta_raw = convention.degrees(math.acos(ratio))
    theta = theta_raw if theta_raw <= 90 else abs(180 - theta_raw)

    return GeometryResult(
        valid=True,
        theta_raw=theta_raw,
        theta=theta,
        max_distance=max_distance,
        degeneracy_reason=None,
        acos_clamped=clamped,
        median_length=s3,
    )
