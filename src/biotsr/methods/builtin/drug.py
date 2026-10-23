"""Drug / ligand TSR method recipes.

Declarative configuration only.  There is **no** ``DrugGeometry``,
``DrugOrdering`` or ``DrugEncoder``: drug TSR is the shared atomic core with a
different selector configuration, and a test asserts no such class exists
(ADR-027).

Two recipes:

``drug``
    Regular drug TSR.  Every admitted heavy atom of the ligand copy identified
    by ``(drug_chain, drug_name, drug_id)`` -- or of every copy of that name on
    the chain when no id is given (ADR-096) -- all ``C(n,3)`` triangles.

``drug_atom_selection``
    The same ligand narrowed to the atom names in ``drug_atom_selected``.

They are **the same method with and without a narrowing list** in
``WebsiteTSR``'s own implementation -- its ligand selector falls back to every
heavy atom when ``atom_selected`` is empty -- and they stay one selector here.
Two names exist because the historical dataset sheet drives them from two
different columns and a method name is what appears in provenance.

**Scope honesty.**

* ``TSR-Package/tsr_package/Drug_TSR.py`` is **not** regular drug TSR: it is
  drug-**protein cross** TSR (its own first line says so), enumerating mixed
  ``1 protein + 2 drug`` and ``1 drug + 2 protein`` triangles at ``dLen = 17``.
  It belongs to the cross methods and nothing here is derived from it.  (An
  earlier capability survey cited it for regular ligand TSR; the source supports
  only the cross reading.)
* No legacy drug implementation has a **size filter** or a **mirror** flag.
  Mirror applies to this family (a TSR Research Group decision), and the modifier is
  family-agnostic by construction (it reads only ``theta_raw``), so it is offered
  and its behaviour is identical to every other family's.  Size filtering is
  offered on the **standard profile only**, with the boundary rule stated
  explicitly (strict ``<`` on ``max_distance``, ADR-024) rather than inherited by
  assumption, and **no legacy reference for a size-filtered drug run exists or
  was fabricated** (ADR-028).
* ``aa_grouping`` is not declared, and would be refused by the label domain even
  if it were (ADR-023).
"""

from __future__ import annotations

from ...candidates.exhaustive import STRATEGY_NAME
from ...profiles.atomic import STANDARD_DRUG_V1
from ...selection import atomic as atomic_selection
from ..registry import Method, register_method

__all__ = ["DRUG_METHOD", "DRUG_ATOM_SELECTION_METHOD"]

_ATOMIC_MODIFIERS = frozenset({"mirror", "size_filter"})

DRUG_METHOD = register_method(
    Method(
        name="drug",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_DRUG_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy=None,  # WebsiteTSR's ligand selector has none.
        composition_rule=None,
        description=(
            "Regular drug/ligand TSR: every admitted heavy HETATM of the ligand "
            "identified by (drug_chain, drug_name, drug_id) -- or, with no id, of "
            "every copy of drug_name on drug_chain in one calculation -- all "
            "C(n,3) triangles, 29 theta classes, 18 atomic distance classes, "
            "radix 112."
        ),
    )
)

DRUG_ATOM_SELECTION_METHOD = register_method(
    Method(
        name="drug_atom_selection",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOM_SELECTION,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_DRUG_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy=None,
        composition_rule=None,
        description=(
            "Selected drug atoms: the ligand of (drug_chain, drug_name, drug_id) "
            "(every copy when no id is given) narrowed to the atom names in "
            "drug_atom_selected, applied to each copy, pooled into one "
            "calculation. A requested atom that is not present is reported in "
            "selection.atoms_missing, never silently dropped."
        ),
    )
)
