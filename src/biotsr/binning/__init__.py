"""Bin schemes: continuous value -> class index."""

from .base import LadderBinScheme
from .distance import (
    ATOMIC_BOUNDARIES,
    CALPHA_BOUNDARIES,
    DISTANCE_SCHEME_ATOMIC_18,
    DISTANCE_SCHEME_CALPHA_35,
)
from .theta import THETA_BOUNDARIES_29, THETA_SCHEME_29

__all__ = [
    "ATOMIC_BOUNDARIES",
    "CALPHA_BOUNDARIES",
    "DISTANCE_SCHEME_ATOMIC_18",
    "DISTANCE_SCHEME_CALPHA_35",
    "LadderBinScheme",
    "THETA_BOUNDARIES_29",
    "THETA_SCHEME_29",
]
