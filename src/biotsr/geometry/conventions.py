"""Angle and degeneracy conventions.

The convention uses ``math.pi``.  The literal ``3.14`` some historical TSR
programs divided by appears nowhere in this package (ADR-002): there is no
flag, no environment variable and no profile that produces those angles.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

__all__ = [
    "AngleConvention",
    "SqrtConvention",
    "SQRT_POW_HALF",
    "SQRT_MATH_SQRT",
    "PI_CONVENTION",
    "MedianConvention",
    "MEDIAN_FROM_MIDPOINT",
    "DegeneracyReason",
    "DegeneratePolicy",
    "DEFAULT_ACOS_TOLERANCE",
]

#: How far outside ``[-1, 1]`` an ``acos`` argument may stray and still be treated
#: as floating-point noise rather than a structurally invalid geometry.  A
#: mathematically valid triangle always yields a ratio in ``[-1, 1]`` exactly, so
#: any excursion is numerical; a *large* excursion means something else is wrong
#: and must not be hidden by clamping.
DEFAULT_ACOS_TOLERANCE = 1e-9


@dataclass(frozen=True)
class SqrtConvention:
    """How a squared distance becomes a distance.

    The legacy family is not unanimous, and the difference is real rather than
    cosmetic:

    * ``TSR-Package`` (protein Ca, SSE) writes ``(...) ** 0.5``;
    * ``TSR_AMINOACID_PACKAGE`` and ``WebsiteTSR`` (every atomic lineage) write
      ``math.sqrt(...)``.

    ``** 0.5`` is evaluated by ``pow``, which IEEE-754 does not require to be
    correctly rounded; ``math.sqrt`` is the correctly-rounded hardware square
    root.  They disagree by one ULP on roughly 0.12% of inputs -- measured, not
    assumed -- and one ULP in a distance can move a triangle across a bin
    boundary and change its key.  So this is a key-affecting convention and it
    enters ``KeySpace`` alongside the angle convention.

    Discovered while reconciling the first byte comparison of
    amino-acid atomic output: 8 of 75 residues differed in the last digit of
    theta and nowhere else.
    """

    name: str
    function: Callable[[float], float]

    def __call__(self, value: float) -> float:
        return self.function(value)

    def content_fingerprint(self) -> str:
        return self.name


#: The Ca lineage's spelling.
SQRT_POW_HALF = SqrtConvention(name="pow_half", function=lambda v: v ** 0.5)

#: Every atomic lineage's spelling.
SQRT_MATH_SQRT = SqrtConvention(name="math_sqrt", function=math.sqrt)


@dataclass(frozen=True)
class AngleConvention:
    """How radians become degrees.

    ``radians_per_half_turn`` is the divisor in ``180 * acos(ratio) / D``.  The
    value BioTSR uses is ``math.pi``.
    With ``D = math.pi`` the pre-fold angle is exactly ``[0, 180]``; with the
    legacy divisor it overshoots to ``[0, 180.0913...]``, which is why the legacy
    fold needs ``abs(180 - theta)``.
    """

    name: str
    radians_per_half_turn: float

    def degrees(self, radians: float) -> float:
        """Convert exactly as the legacy expression does: multiply, then divide."""
        return 180 * radians / self.radians_per_half_turn

    def content_fingerprint(self) -> str:
        return "%s:%r" % (self.name, self.radians_per_half_turn)


#: The one standard convention.
PI_CONVENTION = AngleConvention(name="pi", radians_per_half_turn=math.pi)


@dataclass(frozen=True)
class MedianConvention:
    """How the apex-to-base-midpoint median ``s3`` is computed.

    Algebraically there is one answer; numerically there are two.  BioTSR uses
    the first:

    * ``from_midpoint`` -- build the base midpoint from coordinates and take the
      distance to the apex.  ``TSR-Package/TSR.py``, ``TSR_AMINOACID_PACKAGE``,
      ``CarbohydratesTSR``, ``WebsiteTSR`` and both ``cross_key.py`` scripts all
      spell it this way, so it is the standard.
    * ``from_side_lengths`` -- Apollonius' median formula over the three edge
      lengths, used only by two historical cross programs and not implemented
      here.

    The two disagree in the last bits, because one squares already-rounded edge
    lengths while the other never forms them.  ``s3`` is a denominator in the
    theta ratio, so the difference can move a triangle across a theta-bin
    boundary and change its key.  It is therefore key-affecting and enters
    ``KeySpace``, exactly like ``SqrtConvention``.
    """

    name: str

    def content_fingerprint(self) -> str:
        return self.name


#: The one median convention BioTSR implements.
MEDIAN_FROM_MIDPOINT = MedianConvention(name="from_midpoint")



class DegeneracyReason:
    """The named ways a triple can fail to define a usable triangle.

    Enumerated by the degenerate-triangle audit of the legacy family (ADR-019).

    **Zero area is not one of them (ADR-033).**  A triple can be exactly
    collinear, enclose no area at all, and still have a perfectly well-defined
    theta: ``theta_raw = 0`` when the apex lies beyond the far base point, and
    ``theta_raw = 180`` when it lies beyond the near one.  Both are VALID and are
    encoded; both fold to ``theta = 0``; and because mirror keys off the pre-fold
    angle they receive **opposite** orientations.  What makes a triple
    unusable is an undefined denominator, which is what every reason below is.
    """

    ZERO_BASE_EDGE = "zero_base_edge"
    APEX_AT_BASE_MIDPOINT = "apex_at_base_midpoint"
    COINCIDENT_POINTS = "coincident_points"
    ACOS_DOMAIN_ERROR = "acos_domain_error"


class DegeneratePolicy:
    """What a method does with an invalid ``GeometryResult``.

    Exactly one policy, and deliberately not a second: a mathematically
    **undefined** triple produces **no key** and is recorded with its reason
    (ADR-033).  Forcing an undefined theta to zero, or minting a normal theta bin
    for an undefined triangle, is explicitly not done -- ``theta`` is ``None``,
    not ``0.0``, because absent and zero are different facts.  A *defined* linear
    triple is not affected by this policy at all: it is valid and is encoded.  A
    dedicated scientific "degenerate key" category is deliberately not
    created.

    ``WebsiteTSR/key_triplets_protein.py`` is real, working legacy precedent for
    this behaviour (``try/except (ValueError, ZeroDivisionError): degenerate +=
    1; continue``, with the counter surfaced in its summary output); the other
    eight audited implementations have no protection at all and would propagate
    an uncaught exception.
    """

    SKIP_AND_RECORD = "skip_and_record"
