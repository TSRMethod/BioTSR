"""The 29-class theta scheme.

One table, shared by every method in the legacy family and unchanged here.  The
scheme receives the **folded** theta (``[0, 90]``), never the raw pre-fold angle.

Out-of-range behaviour is closed deliberately: the legacy ladder raises
``UnboundLocalError`` at exactly ``theta == 90.0`` in its if/elif form (and
returns 29 in ``TSR.py``'s generator form).  This scheme is a total function
returning 29 there, which matches the generator form the historical golden outputs
were produced with.

Note a documented open domain question, carried forward unresolved and *not*
silently altered: the 72.69-79.2 degree bin is anomalously wide compared with its
neighbours, identically so in every legacy repo.
"""

from __future__ import annotations

from .base import LadderBinScheme

__all__ = ["THETA_BOUNDARIES_29", "THETA_SCHEME_29"]

THETA_BOUNDARIES_29 = (
    12.11, 17.32, 21.53, 25.21, 28.54, 31.64, 34.55, 37.34, 40.03, 42.64,
    45.17, 47.64, 50.05, 52.43, 54.77, 57.08, 59.38, 61.64, 63.87, 66.09,
    68.30, 70.5, 72.69, 79.2, 81.36, 83.51, 85.67, 87.80, 90.00,
)

THETA_SCHEME_29 = LadderBinScheme(
    name="theta29",
    boundaries=THETA_BOUNDARIES_29,
    n_classes=29,
    overflow_class=29,
)
