"""Cross TSR: two sides, one shared core.

This package contributes **only** what is genuinely cross-specific:

* side definitions and side-specific selection (``sides``);
* contact-pool construction under an explicit atom identity (``pools``);
* the mixed-side candidate strategy (``candidates``).

It contributes no geometry, no canonical ordering, no binning, no key formula and
no serialization -- those come from the same shared components every non-cross
method uses, and an architectural test asserts that no
``DrugProteinGeometry`` / ``CarbProteinGeometry`` / ``NucleotideProteinGeometry``
/ ``ProteinProteinGeometry`` class exists anywhere in the package (ADR-045).
"""

from .candidates import (
    COMPOSITION_RULE,
    STRATEGY_MIXED_SIDE,
    count,
    count_from_sides,
    count_via_complement,
    mixed_side,
)
from .pools import (
    SIDE_A,
    SIDE_B,
    ContactFilter,
    ContactPoolPolicy,
    CrossPools,
    PoolIdentity,
    PoolOrder,
    build_contact_pools,
)
from .runner import CrossRun, run_cross
from .sides import CROSS_SIDE_SELECTORS, CrossSideSpec, select_side

__all__ = [
    "COMPOSITION_RULE",
    "STRATEGY_MIXED_SIDE",
    "SIDE_A",
    "SIDE_B",
    "ContactFilter",
    "ContactPoolPolicy",
    "CrossPools",
    "CrossRun",
    "CrossSideSpec",
    "CROSS_SIDE_SELECTORS",
    "PoolIdentity",
    "PoolOrder",
    "build_contact_pools",
    "count",
    "count_from_sides",
    "count_via_complement",
    "mixed_side",
    "run_cross",
    "select_side",
]
