"""The two amino-acid atomic method recipes.

Declarative configuration, no arithmetic and no engine — the same pattern the
``protein`` recipe follows, which is the point: adding an atomic family cost two
recipes and a selector, not a second TSR implementation.

``amino_acid_atomic``
    Per-residue scope.  Every admitted heavy atom of one residue, all ``C(n,3)``
    triangles.  Carries the exact-atom-count completeness rule, which ADR-004
    makes this method's standard because the behaviour is intentional.

``amino_acid_atom_selection``
    Pooled scope.  Exactly the atoms named in a ``RESIDUE_SEQNUM_ATOM`` list,
    across however many residues, triangulated together.  No completeness rule —
    asking for four atoms and then rejecting the residue for having four would be
    incoherent, and the legacy implementation has none either.

**Supported modifiers.**  Both declare ``mirror`` and ``size_filter``, per the
TSR Research Group's modifier applicability decisions.  Neither declares
``aa_grouping``: grouping is Cα-only, and it is refused twice over — the registry
rejects the undeclared modifier name, and the label domain would reject it anyway
(ADR-023), so removing the declaration alone could never let it slip through.

**A scope honesty note on ``size_filter``.**  No legacy amino-acid atomic
implementation has a size filter at all — not in ``TSR_AMINOACID_PACKAGE``, not in
``WebsiteTSR``.  So there is no historical boundary rule for this family to match
and no legacy reference for a size-filtered run can exist.  The TSR Research
Group decided that size filtering *applies* to this family, so it is offered on the
**standard profile only**, using the one boundary semantics that were actually
characterized anywhere in the family (strict ``<`` on ``max_distance``, ADR-024),
stated explicitly here rather than inherited by assumption.
"""

from __future__ import annotations

from ...candidates.exhaustive import STRATEGY_NAME
from ...profiles.atomic import (
    STANDARD_AMINOACID_ATOM_SELECTION_V1,
    STANDARD_AMINOACID_ATOMIC_V1,
)
from ...selection import atomic as atomic_selection
from ..registry import Method, register_method

__all__ = ["AMINO_ACID_ATOMIC_METHOD", "AMINO_ACID_ATOM_SELECTION_METHOD"]

_ATOMIC_MODIFIERS = frozenset({"mirror", "size_filter"})

AMINO_ACID_ATOMIC_METHOD = register_method(
    Method(
        name="amino_acid_atomic",
        family="atomic",
        selector=atomic_selection.SELECTOR_RESIDUE_ATOMS,
        scope=atomic_selection.SCOPE_PER_ENTITY,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_AMINOACID_ATOMIC_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_atom_count",
        composition_rule=None,
        description=(
            "Amino-acid atomic TSR: every admitted heavy atom of one residue, all "
            "C(n,3) triangles, the 112-radix element labels and the 18-class atomic "
            "distance scheme. One calculation per residue."
        ),
    )
)

AMINO_ACID_ATOM_SELECTION_METHOD = register_method(
    Method(
        name="amino_acid_atom_selection",
        family="atomic",
        selector=atomic_selection.SELECTOR_ATOM_SELECTION,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_AMINOACID_ATOM_SELECTION_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy=None,
        composition_rule=None,
        description=(
            "Amino-acid atom selection: exactly the atoms named as "
            "RESIDUE_SEQNUM_ATOM, pooled across residues into one calculation. "
            "Corresponds to the aa_atom_selection dataset column."
        ),
    )
)


# ---------------------------------------------------------------------------
# The dataset sheet's amino-acid selection modes (ADR-034)
#
# `aa_ind`, `aa_ind_mix` and `aa_ind_specific` are three *recipes* over one
# selector (`select_entity_atoms`), not three implementations. They differ only
# in how the target set is named:
#
#   aa_ind           one amino-acid TYPE               "TYR"
#   aa_ind_mix       two or more amino-acid TYPES      "TYR;PHE"
#   aa_ind_specific  explicit residue identities       "LYS_6;GLU18"
#
# All three use **every eligible non-H/D atom** of every admitted matching
# residue, pooled into one calculation (the dataset-sheet format's column
# definitions say every atom except hydrogen). Pooling is the reading `aa_ind_mix` forces: mixing two
# residue types is only meaningful if their atoms meet in the same triangle set.
# Per-residue scope remains available unchanged through `amino_acid_atomic` /
# `analyze_entities`.
#
# The `_2` columns (`aa_ind_2`, `aa_ind_mix_2`, `aa_ind_specific_2`) carry the
# same semantics against `aa_chain_2` and need no new method: they are these
# methods invoked with the other chain. No cross method is implemented here.
# ---------------------------------------------------------------------------

from ...profiles.atomic import STANDARD_AMINOACID_TYPE_V1  # noqa: E402

__all__ += [
    "AMINO_ACID_IND_METHOD",
    "AMINO_ACID_IND_MIX_METHOD",
    "AMINO_ACID_IND_SPECIFIC_METHOD",
]

AMINO_ACID_IND_METHOD = register_method(
    Method(
        name="aa_ind",
        family="atomic",
        selector=atomic_selection.SELECTOR_RESIDUE_TYPE_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_AMINOACID_TYPE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_atom_count",
        composition_rule=None,
        min_entity_types=1,
        max_entity_types=1,
        description=(
            "aa_ind: one specified amino-acid TYPE in aa_chain. Every eligible "
            "non-H/D atom of every admitted residue of that type, pooled into one "
            "calculation. Atomic label space, radix 112, 18 distance classes."
        ),
    )
)

AMINO_ACID_IND_MIX_METHOD = register_method(
    Method(
        name="aa_ind_mix",
        family="atomic",
        selector=atomic_selection.SELECTOR_RESIDUE_TYPE_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_AMINOACID_TYPE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_atom_count",
        composition_rule=None,
        min_entity_types=2,
        max_entity_types=None,
        description=(
            "aa_ind_mix: two or more specified amino-acid TYPES in aa_chain, "
            "pooled into one calculation. Identical to aa_ind in every respect "
            "except arity, which is validated so a one-type aa_ind_mix is a loud "
            "error rather than a silently different calculation."
        ),
    )
)

AMINO_ACID_IND_SPECIFIC_METHOD = register_method(
    Method(
        name="aa_ind_specific",
        family="atomic",
        selector=atomic_selection.SELECTOR_SPECIFIC_RESIDUE_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_AMINOACID_TYPE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_atom_count",
        composition_rule=None,
        description=(
            "aa_ind_specific: explicitly named amino-acid residues in aa_chain "
            "(RESIDUE_SEQNUM or RESIDUESEQNUM), every eligible non-H/D atom of "
            "each, pooled into one calculation. The residue NAME is matched as "
            "well as the number; a disagreement is reported, never silently "
            "resolved in favour of the structure."
        ),
    )
)
