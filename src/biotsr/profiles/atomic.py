"""Profiles for the atomic (non-Cα) families.

The two amino-acid atomic capabilities each select atoms the way their own
historical program did -- per-residue scope with that program's exact-atom-count
completeness rule and element token for ``amino_acid_atomic``, pooled
spec-correct selection for ``amino_acid_atom_selection`` -- because which atoms a
method selects is part of what the method *is* (ADR-004).  What the standard
changes is key semantics: ``math.pi`` rather than a literal divisor, and atom
labels taken from each atom's element field (ADR-113).

Every atomic profile shares **29 theta classes and 18 distance classes**, because
the distance-bin count is a property of the *family*, not of the profile
(ADR-021).  The Cα 35-class scheme is never reachable from any of them.

**AA grouping is not available here, by construction rather than by omission.**
Grouping requires the ``protein_amino_acid`` label domain; these profiles label
by element, whose domain is ``atom_element``, so
``Profile.with_modifiers(aa_grouping=True)`` raises ``IncompatibleModifier``
naming both domains (ADR-023). No per-family exclusion list is maintained.
"""

from __future__ import annotations

from ..binning.distance import DISTANCE_SCHEME_ATOMIC_18
from ..binning.theta import THETA_SCHEME_29
from ..geometry.conventions import SQRT_MATH_SQRT, DegeneratePolicy, PI_CONVENTION
from ..labels.elements import ELEMENT_LABELS_112
from ..modifiers.mirror import MIRROR_DISABLED
from ..selection.atomic import STANDARD_AMINOACID_RESIDUE_ADMISSION, WEBSITE_TSR_ADMISSION
from ..selection.completeness import EXACT_LEGACY_ATOM_COUNT, NO_COMPLETENESS_POLICY
from ..structure.backends import pdb_legacy
from .standard import Profile

__all__ = [
    "STANDARD_AMINOACID_ATOMIC_V1",
    "STANDARD_AMINOACID_ATOM_SELECTION_V1",
]


STANDARD_AMINOACID_ATOMIC_V1 = Profile(
    name="standard_aminoacid_atomic_v1",
    angle_convention=PI_CONVENTION,
    theta_scheme=THETA_SCHEME_29,
    distance_scheme=DISTANCE_SCHEME_ATOMIC_18,
    label_scheme=ELEMENT_LABELS_112,
    parser_backend=pdb_legacy.BACKEND_NAME,
    sqrt_convention=SQRT_MATH_SQRT,
    mirror_convention=MIRROR_DISABLED,
    size_filter=None,
    atom_admission=STANDARD_AMINOACID_RESIDUE_ADMISSION,
    completeness=EXACT_LEGACY_ATOM_COUNT,
    degenerate_policy=DegeneratePolicy.SKIP_AND_RECORD,
    clamp_acos=True,
    is_legacy_compatibility=False,
    description=(
        "Standard amino-acid atomic TSR: math.pi, the 18-class atomic "
        "distance scheme, 29-class theta, the 112-radix element labels, no mirror "
        "signing, no size filter. The completeness rule is this method's own "
        "historical exact-atom-count behaviour, which ADR-004 makes the standard for "
        "the methods that have one -- it is preserved, not relaxed. Atoms are "
        "admitted by the standard alternate-location rule shared by every "
        "method: blank or A kept, B and later excluded."
    ),
)


#: The pooled atom-selection standard. Identical to STANDARD_AMINOACID_ATOMIC_V1
#: except that it admits atoms by `WebsiteTSR`'s spec-correct rules and carries
#: no completeness policy, as that method always has.
STANDARD_AMINOACID_ATOM_SELECTION_V1 = Profile(
    name="standard_aminoacid_atom_selection_v1",
    angle_convention=PI_CONVENTION,
    theta_scheme=THETA_SCHEME_29,
    distance_scheme=DISTANCE_SCHEME_ATOMIC_18,
    label_scheme=ELEMENT_LABELS_112,
    parser_backend=pdb_legacy.BACKEND_NAME,
    sqrt_convention=SQRT_MATH_SQRT,
    mirror_convention=MIRROR_DISABLED,
    size_filter=None,
    atom_admission=WEBSITE_TSR_ADMISSION,
    completeness=NO_COMPLETENESS_POLICY,
    degenerate_policy=DegeneratePolicy.SKIP_AND_RECORD,
    clamp_acos=True,
    is_legacy_compatibility=False,
    description=(
        "Standard pooled amino-acid atom selection: math.pi, the 18-class "
        "atomic distance scheme, 29-class theta, the 112-radix element labels. Uses "
        "WebsiteTSR's spec-correct element and altLoc admission and no completeness "
        "policy, as this method's historical implementation did."
    ),
)


# ===========================================================================
# The remaining non-cross atomic families
#
# Every profile below is *configuration of the same components*: the same
# geometry, the same canonical ordering, the same 29-class theta scheme, the same
# 18-class atomic distance scheme, the same radix-112 element labels, the same key
# formula.  Nothing here introduces a family-specific engine, and a test asserts
# no `DrugGeometry` / `CarbohydrateGeometry` / `NucleotideGeometry` definition
# exists (ADR-027).
#
# What a family *does* get to choose is stated explicitly on every profile:
#   * which records are eligible          -> AtomAdmission
#   * (atom names are the full PDB field for every lineage, ADR-108)
#   * where its single pass stops         -> AtomAdmission.truncation
#   * whether entities must be complete   -> CompletenessPolicy
#   * which label table it uses           -> LabelScheme (one, shared)
# ===========================================================================

from ..labels.elements import ELEMENT_LABELS_112 as _STD_LABELS  # noqa: E402
from ..selection.atomic import STANDARD_HETERO_ADMISSION  # noqa: E402
from ..selection.completeness import (  # noqa: E402
    AVAILABLE_NUCLEOTIDE_ATOMS,
    EXACT_LEGACY_CARBOHYDRATE_ATOM_COUNT,
    EXACT_LEGACY_NUCLEOTIDE_ATOM_COUNT,
)

__all__ += [
    "STANDARD_AMINOACID_TYPE_V1",
    "STANDARD_DRUG_V1",
    "STANDARD_CARBOHYDRATE_V1",
    "STANDARD_CARBOHYDRATE_CHAIN_V1",
    "STANDARD_NUCLEOTIDE_V1",
    "STANDARD_NUCLEOTIDE_CHAIN_V1",
    "ATOMIC_PROFILES",
]


def _atomic_profile(
    name,
    *,
    angle_convention,
    label_scheme,
    atom_admission,
    completeness,
    is_legacy_compatibility,
    description,
):
    """Build one atomic profile. Only the five fields above ever differ."""
    return Profile(
        name=name,
        angle_convention=angle_convention,
        theta_scheme=THETA_SCHEME_29,
        distance_scheme=DISTANCE_SCHEME_ATOMIC_18,
        label_scheme=label_scheme,
        parser_backend=pdb_legacy.BACKEND_NAME,
        sqrt_convention=SQRT_MATH_SQRT,
        mirror_convention=MIRROR_DISABLED,
        size_filter=None,
        atom_admission=atom_admission,
        completeness=completeness,
        degenerate_policy=DegeneratePolicy.SKIP_AND_RECORD,
        clamp_acos=True,
        is_legacy_compatibility=is_legacy_compatibility,
        description=description,
    )


# --- amino-acid type / specific-residue selection (aa_ind, aa_ind_mix, ------
#     aa_ind_specific) ------------------------------------------------------
#
# Standard-only: no legacy implementation selects residues by amino-acid TYPE
# anywhere in the family.
# Because there is no lineage to pair with, this profile takes the **spec-correct**
# admission rather than the amino-acid analyzer's two off-spec column tokens: a
# new standard capability must not inherit a documented parsing defect. The two
# rules are verified to select identically on every pinned fixture, so the choice
# is measured rather than assumed.
STANDARD_AMINOACID_TYPE_V1 = _atomic_profile(
    "standard_aminoacid_type_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=WEBSITE_TSR_ADMISSION,
    completeness=EXACT_LEGACY_ATOM_COUNT,
    is_legacy_compatibility=False,
    description=(
        "Standard amino-acid selection by residue type or explicit residue "
        "identity (aa_ind / aa_ind_mix / aa_ind_specific): every admitted heavy "
        "atom of every matching residue, pooled into one calculation. 29 theta "
        "classes, 18 atomic distance classes, radix 112. Carries the family's "
        "exact-atom-count completeness rule (ADR-004), applied per residue "
        "instance. Standard-only: no historical implementation selects by type."
    ),
)

# --- drug -------------------------------------------------------------------

STANDARD_DRUG_V1 = _atomic_profile(
    "standard_drug_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=STANDARD_HETERO_ADMISSION,
    completeness=NO_COMPLETENESS_POLICY,
    is_legacy_compatibility=False,
    description=(
        "Standard drug/ligand TSR: the HETATM heavy atoms of one ligand copy, "
        "identified by (drug_chain, drug_name, drug_id), all C(n,3) triangles. "
        "29 theta classes, 18 atomic distance classes, radix 112. No completeness "
        "rule -- WebsiteTSR's ligand selector has none, and none is invented."
    ),
)

# --- carbohydrate -----------------------------------------------------------

STANDARD_CARBOHYDRATE_V1 = _atomic_profile(
    "standard_carbohydrate_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=STANDARD_HETERO_ADMISSION,
    completeness=EXACT_LEGACY_CARBOHYDRATE_ATOM_COUNT,
    is_legacy_compatibility=False,
    description=(
        "Standard carbohydrate monomer/residue TSR: the HETATM heavy atoms of one "
        "carbohydrate residue identified by (carb_chain, carb_name, carb_id). "
        "Carries the historical exact-atom-count rule from expected_atom_counts.json, "
        "which is this scope's own historical behaviour (ADR-004)."
    ),
)

STANDARD_CARBOHYDRATE_CHAIN_V1 = _atomic_profile(
    "standard_carbohydrate_chain_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=STANDARD_HETERO_ADMISSION,
    completeness=NO_COMPLETENESS_POLICY,
    is_legacy_compatibility=False,
    description=(
        "Standard carbohydrate polymer/chain TSR: every admitted carbohydrate "
        "residue of one chain, pooled into one inter-residue calculation. "
        "NO completeness rule, because the historical polymer path explicitly passes "
        "apply_ligand_atom_count_filter=False -- a characterized fact, not an "
        "omission."
    ),
)

# --- nucleotide -------------------------------------------------------------

STANDARD_NUCLEOTIDE_V1 = _atomic_profile(
    "standard_nucleotide_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=WEBSITE_TSR_ADMISSION,
    completeness=EXACT_LEGACY_NUCLEOTIDE_ATOM_COUNT,
    is_legacy_compatibility=False,
    description=(
        "Standard nucleotide monomer/residue TSR: the ATOM heavy atoms of one "
        "nucleotide identified by (nt_chain, nt_name, nt_id), with the historical "
        "exact-atom-count table (G 23, A 22, C 20, U 20, DG 22, DA 21, DC 19, "
        "DT 20) preserved as this method's standard (ADR-004). Uses the "
        "spec-correct element and altLoc rules; the historical analyzer has no altLoc "
        "test at all, which is characterized separately."
    ),
)

STANDARD_NUCLEOTIDE_CHAIN_V1 = _atomic_profile(
    "standard_nucleotide_chain_v1",
    angle_convention=PI_CONVENTION,
    label_scheme=_STD_LABELS,
    atom_admission=WEBSITE_TSR_ADMISSION,
    completeness=AVAILABLE_NUCLEOTIDE_ATOMS,
    is_legacy_compatibility=False,
    description=(
        "Standard DNA/RNA chain TSR: every admitted nucleotide of one strand, "
        "pooled into one inter-nucleotide calculation. STANDARD-ONLY -- no historical "
        "implementation pools nucleotide residues. Each of the eight standard "
        "nucleotides takes part with every heavy atom it has, at a chain end or "
        "inside the chain; a nucleotide whose atoms differ from the standard set "
        "is named in a warning (ADR-114). No atom count is required and none is "
        "reconstructed."
    ),
)

#: Every atomic profile, for the tests that assert family-wide invariants
#: (radix 112, 18 distance classes, no Ca scheme anywhere).
ATOMIC_PROFILES = (
    STANDARD_AMINOACID_ATOMIC_V1,
    STANDARD_AMINOACID_ATOM_SELECTION_V1,
    STANDARD_AMINOACID_TYPE_V1,
    STANDARD_DRUG_V1,
    STANDARD_CARBOHYDRATE_V1,
    STANDARD_CARBOHYDRATE_CHAIN_V1,
    STANDARD_NUCLEOTIDE_V1,
    STANDARD_NUCLEOTIDE_CHAIN_V1,
)
