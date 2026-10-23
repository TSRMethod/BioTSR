"""Profiles for the cross (two-sided) families.

Every profile here is **configuration of the same components** the non-cross
families use: the same geometry function, the same canonical-ordering function,
the same 29-class theta scheme, the same key formula, the same encoder.  Nothing
in this module introduces an engine, and an architectural test asserts that no
``DrugProteinGeometry`` / ``CarbProteinGeometry`` / ``NucleotideProteinGeometry``
/ ``ProteinProteinGeometry`` definition exists anywhere in the package (ADR-045).

The historical cross programs used two other key spaces (a literal ``3.14``
divisor in one lineage; a 17-class distance scale, string collation in the
ordering and an edge-length median in the other).  None of their conventions
is adopted here: there is one standard cross key space, below.

The standard cross key space
----------------------------

One key space for every fully atomic cross method, whichever family:

* ``math.pi``;
* 29 theta classes;
* **18** atomic distance classes -- the same scheme every non-cross atomic family
  uses (ADR-021), so a standard cross key and a standard non-cross atomic key are
  built on the same scales;
* radix **112**, atoms labelled by element (ADR-113);
* the shared four-case ordering over integer labels;
* the midpoint median;
* mirror off, no size filter, undefined geometry skipped and recorded.

A cross method reaches it by naming ``STANDARD_ATOMIC_CROSS_V1``.  What differs
between the drug, carbohydrate, nucleotide and protein-protein cross methods is
*selection* -- which points each side offers and how the contact pool is built --
and selection lives in ``SamplingProfile``, never here.
"""

from __future__ import annotations

from dataclasses import replace

from ..binning.distance import DISTANCE_SCHEME_ATOMIC_18
from ..binning.theta import THETA_SCHEME_29
from ..geometry.conventions import (
    MEDIAN_FROM_MIDPOINT,
    PI_CONVENTION,
    SQRT_MATH_SQRT,
    DegeneratePolicy,
)
from ..geometry.ordering import STANDARD_ORDERING
from ..labels.elements import ELEMENT_LABELS_112
from ..modifiers.mirror import MIRROR_DISABLED
from ..structure.backends import pdb_legacy
from .standard import STANDARD_PROTEIN_CALPHA_V1, Profile

__all__ = [
    "STANDARD_ATOMIC_CROSS_V1",
    "STANDARD_CALPHA_CROSS_V1",
    "CROSS_PROFILES",
]


STANDARD_ATOMIC_CROSS_V1 = Profile(
    name="standard_atomic_cross_v1",
    angle_convention=PI_CONVENTION,
    theta_scheme=THETA_SCHEME_29,
    distance_scheme=DISTANCE_SCHEME_ATOMIC_18,
    label_scheme=ELEMENT_LABELS_112,
    parser_backend=pdb_legacy.BACKEND_NAME,
    sqrt_convention=SQRT_MATH_SQRT,
    median_convention=MEDIAN_FROM_MIDPOINT,
    ordering_convention=STANDARD_ORDERING,
    mirror_convention=MIRROR_DISABLED,
    size_filter=None,
    degenerate_policy=DegeneratePolicy.SKIP_AND_RECORD,
    clamp_acos=True,
    is_legacy_compatibility=False,
    description=(
        "Standard fully atomic cross TSR, shared by drug-protein, "
        "carbohydrate-protein, nucleotide/DNA/RNA-protein and atomic "
        "protein-protein cross: math.pi, 29 theta classes, the 18-class atomic "
        "distance scheme, the 112-radix element labels, the shared four-case "
        "ordering over integer labels, the midpoint median, no mirror signing, "
        "no size filter. Which points take part is a SamplingProfile concern: "
        "contact cutoff, side selectors and pool identity never enter this key "
        "space."
    ),
)


#: Ca protein-protein cross.  A Ca point **is a residue**, so it is labelled in
#: the residue label space at radix 30 on the 35-class Ca distance scale --
#: the same key space as ordinary protein Ca TSR and as ``aa_selection``, so a Ca
#: cross key and an ordinary Ca key are directly comparable.  It is emphatically
#: **not** the atomic 112 space: encoding a residue with atom labels would
#: make the two incomparable while looking identical (ADR-052).
#:
#: What makes it a *cross* is entirely in ``SamplingProfile``: which residues were
#: drawn (the interface, by heavy-atom contact) and how they were paired
#: (mixed-side AAB/ABB).  Nothing about the key formula changes.
STANDARD_CALPHA_CROSS_V1 = replace(
    STANDARD_PROTEIN_CALPHA_V1,
    name="standard_calpha_cross_v1",
    description=(
        "Standard Ca protein-protein cross: math.pi, 29 theta classes, the "
        "35-class Ca distance scheme, the residue labels at radix 30, the "
        "shared four-case ordering, the midpoint median, float64. Identical key "
        "space to ordinary protein Ca TSR -- only which residues take part "
        "differs, and that is a SamplingProfile concern: interface residues are "
        "identified by ANY non-hydrogen atom pair within the cutoff and are then "
        "represented by their Ca, which itself need not be within the cutoff."
    ),
)


#: Every cross profile, for the tests that assert family-wide invariants.
CROSS_PROFILES = (
    STANDARD_ATOMIC_CROSS_V1,
    STANDARD_CALPHA_CROSS_V1,
)
