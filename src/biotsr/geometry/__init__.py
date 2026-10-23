"""Canonical ordering, triangle geometry, and the angle/degeneracy conventions."""

from .conventions import (
    AngleConvention,
    DegeneracyReason,
    DegeneratePolicy,
    PI_CONVENTION,
    SQRT_MATH_SQRT,
    SQRT_POW_HALF,
    SqrtConvention,
)
from .ordering import (
    CANONICAL_ORDERING_VERSION,
    CanonicalOrder,
    OrderingCase,
    canonical_order,
)
from .triangle import GeometryResult, triangle_geometry

__all__ = [
    "AngleConvention",
    "SQRT_MATH_SQRT",
    "SQRT_POW_HALF",
    "SqrtConvention",
    "CANONICAL_ORDERING_VERSION",
    "CanonicalOrder",
    "DegeneracyReason",
    "DegeneratePolicy",
    "GeometryResult",
    "OrderingCase",
    "PI_CONVENTION",
    "canonical_order",
    "triangle_geometry",
]
