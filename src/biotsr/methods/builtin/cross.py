"""Cross TSR method recipes -- declarative configuration only.

Four families, one framework.  Each recipe pins:

* the two side selectors and the admission rule each side reads records under;
* the contact-pool policy (cutoff, which sides it narrows, pool identity, order);
* the candidate strategy (always mixed-side -- the composition rule is the same
  for every cross family and is not a per-family choice);
* the default profile.

None of them pins geometry, ordering, binning or a key formula, because none of
them may: those are shared and identical across all four families.  There is no
``DrugProteinCross`` class, no ``CarbProteinGeometry``, and a test asserts it.

Scope and contact rule
----------------------

Decided by the TSR Research Group: **both participating chains are
user-specified**, for every cross method, and a chain the caller did not name is
never scanned.  Hydrogen is excluded everywhere -- contact detection, interface
identification, pools, vertices.  The contact cutoff narrows **both** sides
symmetrically and acts at the atom level: an atom enters iff it has at least one
non-hydrogen atom on the opposite specified side within the cutoff, and is never
admitted merely because another atom of its residue or of its ligand makes
contact (ADR-049, ADR-050, ADR-051).  The historical programs each did something
narrower or broader (a whole-structure protein scan, a one-sided drug filter);
none of those scopes is reproduced.

Modifiers
---------

``mirror`` and ``sse`` are **not** declared on any cross method.  Neither has been
defined or validated for cross calculations, so both are recorded as future
scope rather than offered untested (ADR-046).

``size_filter`` is declared on every cross method.  It is a different concept
from the contact cutoff and the two never interact: the cutoff decides pool
membership before any triangle exists, the size filter rejects a triangle after
its geometry is known.  It stays OFF by default.

``aa_grouping`` is not declared on any cross method.  A fully atomic cross labels
atoms by element, where the amino-acid grouping has no meaning at all; and
on ``protein_protein_calpha_cross``, where it *would* have a meaning, the
standard method is scoped to the ungrouped residue labels and grouping is not
part of this release, so it is not offered (ADR-052).
"""

from __future__ import annotations

from dataclasses import replace

from ...cross.candidates import COMPOSITION_RULE, STRATEGY_MIXED_SIDE
from ...cross import pools as cross_pools
from ...cross import sides as cross_sides
from ...profiles.cross import STANDARD_ATOMIC_CROSS_V1, STANDARD_CALPHA_CROSS_V1
from ...selection.atomic import STANDARD_HETERO_ADMISSION, WEBSITE_TSR_ADMISSION
from ...selection.base import CALPHA_MODIFIED_RESIDUES
from ..registry import Method, register_method

__all__ = [
    "CROSS_FAMILY",
    "SCOPE_CROSS",
    "CrossRecipe",
    "DRUG_PROTEIN_CROSS",
    "CARBOHYDRATE_PROTEIN_CROSS",
    "NUCLEOTIDE_PROTEIN_CROSS",
    "PROTEIN_PROTEIN_CROSS",
    "PROTEIN_PROTEIN_CALPHA_CROSS",
    "CROSS_METHODS",
    "cross_recipe",
]

#: The heavy atoms that decide Calpha interface residues: WebsiteTSR's residue
#: admission, reading ``HETATM`` records of the modified residues the Calpha
#: selector represents (ADR-107) -- so a supported modified residue decides
#: interface membership with all its heavy atoms, exactly like a regular one --
#: and of no other ``HETATM`` residue.
_CALPHA_INTERFACE_ADMISSION = replace(
    WEBSITE_TSR_ADMISSION,
    name="calpha_interface",
    hetero_residue_names=CALPHA_MODIFIED_RESIDUES,
)

CROSS_FAMILY = "cross"
SCOPE_CROSS = "cross_two_sided"

_CROSS_MODIFIERS = frozenset({"size_filter"})


class CrossRecipe:
    """The cross-specific half of a method: sides, admissions, pool policies.

    Kept beside the ``Method`` rather than inside it so the registry stays the
    one shared, family-agnostic type.  ``Method`` still carries everything
    provenance needs; this carries what the cross runner needs to build pools.

    ``entity_side`` says which side is the named entity (drug / carbohydrate /
    nucleotide) and which is the protein, so the public API can route the
    caller's arguments by **role** rather than by side letter.  ``None`` means
    both sides are protein.
    """

    __slots__ = (
        "name",
        "entity_side",
        "side_a_role",
        "side_b_role",
        # -- standard ------------------------------------------------------
        "side_a_selector",
        "side_b_selector",
        "side_a_admission",
        "side_b_admission",
        "side_a_entity_names",
        "side_b_entity_names",
        "default_cutoff",
        "standard_pool",
        #: Heavy-atom admission used to classify interface residues when the
        #: standard pool projects onto Ca. ``None`` for every atomic method,
        #: where the contacting atom IS the vertex and no second selection is
        #: needed.
        "interface_admission",
    )

    def __init__(self, **kw):
        for slot in self.__slots__:
            setattr(self, slot, kw.pop(slot, None))
        if kw:
            raise TypeError("unexpected CrossRecipe fields: %s" % ", ".join(sorted(kw)))


_RECIPES = {}


def cross_recipe(name) -> CrossRecipe:
    try:
        return _RECIPES[name]
    except KeyError:
        raise KeyError(
            "%r is not a cross method; cross methods are: %s"
            % (name, ", ".join(sorted(_RECIPES)))
        )


def _register(method: Method, recipe: CrossRecipe) -> Method:
    register_method(method)
    _RECIPES[method.name] = recipe
    return method


# ---------------------------------------------------------------------------
# drug - protein
# ---------------------------------------------------------------------------

DRUG_PROTEIN_CROSS = _register(
    Method(
        name="drug_protein_cross",
        family=CROSS_FAMILY,
        selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        scope=SCOPE_CROSS,
        candidate_strategy=STRATEGY_MIXED_SIDE,
        default_profile=STANDARD_ATOMIC_CROSS_V1,
        supported_modifiers=_CROSS_MODIFIERS,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        description=(
            "Drug-protein cross TSR: mixed-side triangles over a contact pool "
            "built from one ligand copy and the protein atoms near it. 29 theta "
            "classes, 18 atomic distance classes, radix 112."
        ),
    ),
    CrossRecipe(
        name="drug_protein_cross",
        entity_side="b",
        side_a_role="protein",
        side_b_role="drug",
        side_a_selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        side_b_selector=cross_sides.SIDE_SELECTOR_ENTITY,
        # `WEBSITE_TSR_ADMISSION` reads ATOM records and the recognised modified
        # amino acids (TPO, SEP, PTR; ADR-101) with the spec-correct element
        # field and excludes H and D, which is the standard heavy-atom rule.
        side_a_admission=WEBSITE_TSR_ADMISSION,
        side_b_admission=STANDARD_HETERO_ADMISSION,
        default_cutoff=3.5,
        # The standard does not keep `Drug_TSR.py`'s asymmetry (whole ligand,
        # contact-narrowed protein): contact filtering applies to BOTH sides, so
        # a ligand atom is retained only when it has a qualifying protein atom
        # within the cutoff. An explicit `drug_atom_selected` is applied FIRST
        # and the selected atoms are then contact-filtered (ADR-050, ADR-051).
        standard_pool=cross_pools.STANDARD_SYMMETRIC_ATOM_POOLS,
    ),
)


# ---------------------------------------------------------------------------
# carbohydrate - protein
# ---------------------------------------------------------------------------

CARBOHYDRATE_PROTEIN_CROSS = _register(
    Method(
        name="carbohydrate_protein_cross",
        family=CROSS_FAMILY,
        selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        scope=SCOPE_CROSS,
        candidate_strategy=STRATEGY_MIXED_SIDE,
        default_profile=STANDARD_ATOMIC_CROSS_V1,
        supported_modifiers=_CROSS_MODIFIERS,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        description=(
            "Carbohydrate-protein cross TSR: mixed-side triangles over a contact "
            "pool built from one carbohydrate residue (or a listed vocabulary of "
            "them) and the protein atoms near it."
        ),
    ),
    CrossRecipe(
        name="carbohydrate_protein_cross",
        entity_side="a",
        side_a_role="carbohydrate",
        side_b_role="protein",
        side_a_selector=cross_sides.SIDE_SELECTOR_ENTITY,
        side_b_selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        side_a_admission=cross_sides.STANDARD_CROSS_ENTITY_ADMISSION,
        side_b_admission=cross_sides.STANDARD_CROSS_PROTEIN_ADMISSION,
        side_b_entity_names=cross_sides.PROTEIN_RESIDUE_NAMES,
        default_cutoff=3.5,
        standard_pool=cross_pools.STANDARD_SYMMETRIC_ATOM_POOLS,
    ),
)


# ---------------------------------------------------------------------------
# nucleotide / DNA / RNA - protein
# ---------------------------------------------------------------------------

NUCLEOTIDE_PROTEIN_CROSS = _register(
    Method(
        name="nucleotide_protein_cross",
        family=CROSS_FAMILY,
        selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        scope=SCOPE_CROSS,
        candidate_strategy=STRATEGY_MIXED_SIDE,
        default_profile=STANDARD_ATOMIC_CROSS_V1,
        supported_modifiers=_CROSS_MODIFIERS,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        default_entity_names=cross_sides.NUCLEOTIDE_RESIDUE_NAMES,
        description=(
            "Nucleotide/DNA/RNA-protein cross TSR: mixed-side triangles over a "
            "contact pool built from one nucleotide residue (or a whole strand) "
            "and the protein atoms near it."
        ),
    ),
    CrossRecipe(
        name="nucleotide_protein_cross",
        entity_side="a",
        side_a_role="nucleotide",
        side_b_role="protein",
        side_a_selector=cross_sides.SIDE_SELECTOR_ENTITY,
        side_b_selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        side_a_admission=cross_sides.STANDARD_CROSS_ENTITY_ADMISSION,
        side_b_admission=cross_sides.STANDARD_CROSS_PROTEIN_ADMISSION,
        side_a_entity_names=cross_sides.NUCLEOTIDE_RESIDUE_NAMES,
        side_b_entity_names=cross_sides.PROTEIN_RESIDUE_NAMES,
        # The standard reference value for nucleotide/DNA/RNA-protein cross
        # is 3.5 A (the historical `drug_protein_3A.py` hard-coded 3.0). It is a
        # default, not a constant -- `contact_cutoff=` overrides it.
        default_cutoff=3.5,
        standard_pool=cross_pools.STANDARD_SYMMETRIC_ATOM_POOLS,
    ),
)


# ---------------------------------------------------------------------------
# protein - protein
# ---------------------------------------------------------------------------

PROTEIN_PROTEIN_CROSS = _register(
    Method(
        name="protein_protein_cross",
        family=CROSS_FAMILY,
        selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        scope=SCOPE_CROSS,
        candidate_strategy=STRATEGY_MIXED_SIDE,
        default_profile=STANDARD_ATOMIC_CROSS_V1,
        supported_modifiers=_CROSS_MODIFIERS,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        description=(
            "Protein-protein cross TSR: mixed-side triangles over a contact pool "
            "built from two protein chains of one structure -- aa_chain and "
            "aa_chain_2. Fully ATOMIC: both sides are generic atoms labelled by "
            "element (112-radix element labels), as Cross_TSR.py's atoms are and as "
            "aa_ind* / aa_atom_selection mean. Ca residue selections are a "
            "different key space (residue labels, radix 30) and are NOT mixed "
            "into this all-atom "
            "method; protein_protein_calpha_cross is the Ca method (ADR-045, "
            "ADR-052)."
        ),
    ),
    CrossRecipe(
        name="protein_protein_cross",
        entity_side=None,
        side_a_role="protein",
        side_b_role="protein",
        side_a_selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        side_b_selector=cross_sides.SIDE_SELECTOR_PROTEIN_CHAIN,
        side_a_admission=WEBSITE_TSR_ADMISSION,
        side_b_admission=WEBSITE_TSR_ADMISSION,
        default_cutoff=5.0,
        standard_pool=cross_pools.STANDARD_SYMMETRIC_ATOM_POOLS,
    ),
)


# ---------------------------------------------------------------------------
# protein - protein, Ca
# ---------------------------------------------------------------------------

PROTEIN_PROTEIN_CALPHA_CROSS = _register(
    Method(
        name="protein_protein_calpha_cross",
        family=CROSS_FAMILY,
        selector=cross_sides.SIDE_SELECTOR_CALPHA_CHAIN,
        scope=SCOPE_CROSS,
        candidate_strategy=STRATEGY_MIXED_SIDE,
        default_profile=STANDARD_CALPHA_CROSS_V1,
        supported_modifiers=_CROSS_MODIFIERS,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        description=(
            "Ca protein-protein cross TSR between two user-specified chains. Two "
            "steps, and they are not the same rule: a residue is an INTERFACE "
            "residue when at least one pair of non-hydrogen atoms -- side chains "
            "included -- from the two chains lies within the cutoff; each "
            "interface residue is then represented by its Ca, which itself need "
            "not be within the cutoff. Residue points, so the residue labels "
            "(20 regular amino acids and the supported modified residues) at "
            "radix 30 and the 35-class Ca distance scheme, NOT the atomic 112 "
            "space."
        ),
    ),
    CrossRecipe(
        name="protein_protein_calpha_cross",
        entity_side=None,
        side_a_role="protein",
        side_b_role="protein",
        side_a_selector=cross_sides.SIDE_SELECTOR_CALPHA_CHAIN,
        side_b_selector=cross_sides.SIDE_SELECTOR_CALPHA_CHAIN,
        # The Ca selector supplies the POINTS; this admission supplies the heavy
        # atoms that decide which residues are interface residues. Two selections
        # on each side, because the rule genuinely has two steps.
        side_a_admission=WEBSITE_TSR_ADMISSION,
        side_b_admission=WEBSITE_TSR_ADMISSION,
        # ATOM records and the supported modified residues: a residue that has a
        # Ca point decides interface membership with all its heavy atoms; a
        # modified residue without a residue label has neither and is reported
        # as not represented (ADR-101, ADR-107), never half-included.
        interface_admission=_CALPHA_INTERFACE_ADMISSION,
        default_cutoff=5.0,
        standard_pool=cross_pools.STANDARD_INTERFACE_CALPHA_POOLS,
    ),
)


CROSS_METHODS = (
    DRUG_PROTEIN_CROSS,
    CARBOHYDRATE_PROTEIN_CROSS,
    NUCLEOTIDE_PROTEIN_CROSS,
    PROTEIN_PROTEIN_CROSS,
    PROTEIN_PROTEIN_CALPHA_CROSS,
)
