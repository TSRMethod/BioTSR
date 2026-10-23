"""Carbohydrate TSR method recipes.

Declarative configuration only -- no ``CarbohydrateGeometry`` and no second TSR
engine (ADR-027).  Carbohydrate behaviour is read from
``CarbohydratesTSR/key_generator/key_generation_using_id_chain_res_seq/
key_triplets_generation_dual_mode.py`` and its ``expected_atom_counts.json``,
**never inferred from amino acids** -- and four of its rules really are different:

* it reads ``HETATM`` records, not ``ATOM``;
* it deliberately **does not stop at ``TER``** (its own comment: *"Keep scanning
  after TER because carbohydrate HETATM records can appear later in the same
  chain"*), so the peptide truncation rule would select nothing at all;
* its completeness table is a carbohydrate table (eight names, 10-14 atoms), and
  it is applied **per entity instance** so a chain-scoped run drops only the
  incomplete monomers;
* the **polymer scope has no completeness filter at all** -- the grouped caller
  passes ``apply_ligand_atom_count_filter=False``.

Three recipes:

``carbohydrate``
    Monomer / residue scope.  One carbohydrate residue, identified by
    ``(carb_chain, carb_name, carb_id)``.  Carries the exact-atom-count rule.

``carbohydrate_chain``
    Polymer / chain scope.  Every admitted carbohydrate residue of one chain,
    pooled into one inter-residue calculation.  **No completeness filter.**

``carbohydrate_atom_selection``
    The monomer narrowed to the atom names in ``carb_atom_selected``.

**Which residues count as carbohydrate.**  The legacy polymer path takes the
vocabulary from the *keys of the atom-count JSON* -- a scoping behaviour the legacy
audit flagged, because it silently diverges from the CSV's own ``carb_name`` column.
It is preserved, but it is now **declared on the method** (``default_entity_names``)
instead of being an implicit side effect of a configuration file, and it is
overridable per call.  Nothing about it is hidden.

**``carb_atom_selected`` has no carbohydrate-specific legacy implementation.**
``CarbohydratesTSR`` has no atom-narrowing path at all -- verified by a
repository-wide search, which finds no ``atom_selected`` anywhere in it.  The
capability is requested by the TSR Research Group and its *mechanism* is fully
characterized elsewhere in the family (``WebsiteTSR``'s ligand selector narrows a
HETATM entity by exactly this kind of name list), so it is offered on the
standard profile, and no historical reference for it exists or was fabricated.
"""

from __future__ import annotations

from ...candidates.exhaustive import STRATEGY_NAME
from ...profiles.atomic import STANDARD_CARBOHYDRATE_CHAIN_V1, STANDARD_CARBOHYDRATE_V1
from ...selection import atomic as atomic_selection
from ...selection.completeness import CARBOHYDRATE_RESIDUE_NAMES
from ..registry import Method, register_method

__all__ = [
    "CARBOHYDRATE_METHOD",
    "CARBOHYDRATE_CHAIN_METHOD",
    "CARBOHYDRATE_ATOM_SELECTION_METHOD",
]

_ATOMIC_MODIFIERS = frozenset({"mirror", "size_filter"})

CARBOHYDRATE_METHOD = register_method(
    Method(
        name="carbohydrate",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_CARBOHYDRATE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_carbohydrate_atom_count",
        composition_rule=None,
        description=(
            "Carbohydrate monomer/residue TSR: every admitted heavy HETATM of the "
            "carbohydrate residue identified by (carb_chain, carb_name, carb_id) "
            "-- or, with no id, of every residue of that name on carb_chain in "
            "one calculation (not the whole glycan: see carbohydrate_chain) -- "
            "all C(n,3) triangles, radix 112, 18 atomic distance classes. Carries "
            "the historical exact-atom-count rule."
        ),
    )
)

CARBOHYDRATE_CHAIN_METHOD = register_method(
    Method(
        name="carbohydrate_chain",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_CHAIN_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_CARBOHYDRATE_CHAIN_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy=None,
        composition_rule=None,
        default_entity_names=CARBOHYDRATE_RESIDUE_NAMES,
        description=(
            "Carbohydrate polymer/chain TSR: every carbohydrate residue of one "
            "chain pooled into a single inter-residue calculation. The residue "
            "vocabulary is declared here (BGC, BMA, FUC, FUL, GAL, MAN, NAG, NDG) "
            "and can be overridden per call. No completeness filter -- the historical "
            "grouped path has none."
        ),
    )
)

CARBOHYDRATE_ATOM_SELECTION_METHOD = register_method(
    Method(
        name="carbohydrate_atom_selection",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOM_SELECTION,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_CARBOHYDRATE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_carbohydrate_atom_count",
        composition_rule=None,
        description=(
            "Selected carbohydrate atoms: the monomer of (carb_chain, carb_name, "
            "carb_id) -- every copy when no id is given -- narrowed to the atom "
            "names in carb_atom_selected, applied to each copy."
        ),
    )
)
