"""Distance schemes: the 35-class C-alpha scale and the 18-class atomic scale.

The C-alpha scheme has a real 35th class for the tail:

    83.26 <= d < 132.45  ->  34
    d >= 132.45          ->  35

This is the standard for every Ca-based family (ADR-014), confirmed by direct
execution of ``WebsiteTSR/utils/distance_calpha_utils.py``.  (One historical
program merged the top two classes despite declaring 35; BioTSR does not
reproduce that defect.)

The **atomic** scheme is here too. Which scheme a method uses is
a property of its *family*, not a free choice (ADR-021):

    Protein Ca                      -> 29 theta, 35 distance
    every non-Ca atomic family      -> 29 theta, 18 distance

Substituting one for the other cannot silently produce comparable keys, because
the scheme's content fingerprint enters ``KeySpace`` -- the guard fires instead.
"""

from __future__ import annotations

from .base import LadderBinScheme

__all__ = [
    "CALPHA_BOUNDARIES",
    "DISTANCE_SCHEME_CALPHA_35",
    "ATOMIC_BOUNDARIES",
    "DISTANCE_SCHEME_ATOMIC_18",
]

#: The 34 cut points of the Ca scale.
CALPHA_BOUNDARIES = (
    3.83, 7.00, 9.00, 11.00, 14.00, 17.99, 21.25, 23.19, 24.8, 26.26,
    27.72, 28.9, 30.36, 31.62, 32.76, 33.84, 35.13, 36.26, 37.62, 38.73,
    40.12, 41.8, 43.41, 45.55, 47.46, 49.69, 52.65, 55.81, 60.2, 64.63,
    70.04, 76.15, 83.26, 132.45,
)

#: The standard: a real 35th class for the tail.
DISTANCE_SCHEME_CALPHA_35 = LadderBinScheme(
    name="distance_calpha_35",
    boundaries=CALPHA_BOUNDARIES,
    n_classes=35,
    overflow_class=35,
)

#: The atomic scale: seventeen 1 A cut points, everything from 17 A up in class 18.
#: Verbatim from ``TSR_AMINOACID_PACKAGE/.../utils/distance_utils.py``, confirmed
#: byte-identical (modulo line endings) to ``WebsiteTSR/utils/distance_utils.py``
#: -- one table shared by the whole atomic family.
ATOMIC_BOUNDARIES = tuple(float(n) for n in range(1, 18))

#: Unlike the Ca binner, this one is internally consistent: seventeen boundaries
#: give eighteen reachable classes, and it declares eighteen.
DISTANCE_SCHEME_ATOMIC_18 = LadderBinScheme(
    name="distance_atomic_18",
    boundaries=ATOMIC_BOUNDARIES,
    n_classes=18,
    overflow_class=18,
)


