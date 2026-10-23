"""Nucleotide / DNA / RNA TSR method recipes.

Declarative configuration only -- no ``NucleotideGeometry`` (ADR-027).  Read from
``TSR_NUCLEOTIDE_PACKAGE/nucleotide_tsr_package/analyzer.py``, and **not** assumed
to match its amino-acid sibling, which mattered: two of its admission rules
genuinely differ.

===============================  ==================================  ============
Rule                             Nucleotide analyzer                 AA analyzer
===============================  ==================================  ============
record type                      ``ATOM`` only                       same
residue-**name** predicate       **absent** -- chain + seq only       present
altLoc predicate                 **absent**                           ``[16:20]``
                                                                      startswith B
H/D test                         ``line[77:80]``                     same
atom name                        ``line[13:16].rstrip()``            same
truncation                       ``ENDMDL`` / ``TER`` (chain) /       same
                                 ``MODEL > 1``
completeness                     exact counts, 8 residue names       23 names
                                 (``nucleotide_chain``: none)
===============================  ==================================  ============

The missing residue-name predicate is not cosmetic: two residues sharing a
sequence label in one chain would be merged by the legacy nucleotide reader and
kept apart by the amino-acid one.  It is **not** carried into the standard,
which matches name and sequence like every other entity family.

Three recipes:

``nucleotide``
    Monomer / residue scope, the only scope the legacy source implements.

``nucleotide_chain``
    DNA/RNA strand scope: every admitted nucleotide of one chain pooled into one
    inter-nucleotide calculation.  **Standard-only.**  No legacy
    program implements polymer nucleotide TSR; per-nucleotide atomic TSR does
    not mean polymer TSR.  So there is no historical behaviour to compare with,
    and no fabricated reference; the scope is composed from this family's own
    admission rule and the shared pooled candidate strategy that the
    carbohydrate polymer scope already uses, and is labelled standard-only
    everywhere it appears.  Like ``carbohydrate_chain`` it requires no atom
    count: each of the eight standard nucleotides takes part with every heavy
    atom it has, at a chain end or inside the chain, and each one whose atoms
    differ from the standard set is named in a warning (ADR-114).  The
    residue-specific methods below keep the exact count.

``nucleotide_atom_selection``
    One nucleotide narrowed to the atom names in ``nt_atom_selected``.  Like the
    carbohydrate case this has no nucleotide-specific legacy implementation; the
    mechanism is ``WebsiteTSR``'s, the capability is a TSR Research Group
    request, and it is standard-only.

DNA and RNA are distinguished by **residue name**, which is the only distinction
the legacy source makes: the completeness table lists ``DA/DC/DG/DT`` (DNA) and
``A/C/G/U`` (RNA) with different atom counts, and nothing else in the code
separates them.  ``DNA_RESIDUE_NAMES`` / ``RNA_RESIDUE_NAMES`` make that explicit
so a caller can scope a chain run to one or the other without inventing a
"strand type" concept the science does not have.
"""

from __future__ import annotations

from ...candidates.exhaustive import STRATEGY_NAME
from ...profiles.atomic import STANDARD_NUCLEOTIDE_CHAIN_V1, STANDARD_NUCLEOTIDE_V1
from ...selection import atomic as atomic_selection
from ...selection.completeness import NUCLEOTIDE_RESIDUE_NAMES
from ..registry import Method, register_method

__all__ = [
    "NUCLEOTIDE_METHOD",
    "NUCLEOTIDE_CHAIN_METHOD",
    "NUCLEOTIDE_ATOM_SELECTION_METHOD",
]

_ATOMIC_MODIFIERS = frozenset({"mirror", "size_filter"})

NUCLEOTIDE_METHOD = register_method(
    Method(
        name="nucleotide",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_NUCLEOTIDE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_nucleotide_atom_count",
        composition_rule=None,
        description=(
            "Nucleotide monomer/residue TSR: every admitted heavy ATOM of the "
            "nucleotide identified by (nt_chain, nt_name, nt_id) -- or, with no "
            "id, of every nucleotide of that name on nt_chain in one "
            "calculation -- all C(n,3) "
            "triangles, radix 112, 18 atomic distance classes, with the historical "
            "exact-atom-count table (G 23, A 22, C 20, U 20, DG 22, DA 21, "
            "DC 19, DT 20)."
        ),
    )
)

NUCLEOTIDE_CHAIN_METHOD = register_method(
    Method(
        name="nucleotide_chain",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_CHAIN_ATOMS,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_NUCLEOTIDE_CHAIN_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="available_nucleotide_atoms",
        composition_rule=None,
        default_entity_names=NUCLEOTIDE_RESIDUE_NAMES,
        description=(
            "DNA/RNA chain TSR: every admitted nucleotide of one strand pooled "
            "into a single inter-nucleotide calculation, each with every heavy "
            "atom it has; nucleotides whose atoms differ from the standard set "
            "are named. Pass "
            "entity_names=DNA_RESIDUE_NAMES or "
            "RNA_RESIDUE_NAMES to scope a mixed chain."
        ),
    )
)

NUCLEOTIDE_ATOM_SELECTION_METHOD = register_method(
    Method(
        name="nucleotide_atom_selection",
        family="atomic",
        selector=atomic_selection.SELECTOR_ENTITY_ATOM_SELECTION,
        scope=atomic_selection.SCOPE_POOLED,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_NUCLEOTIDE_V1,
        supported_modifiers=_ATOMIC_MODIFIERS,
        completeness_policy="exact_legacy_nucleotide_atom_count",
        composition_rule=None,
        description=(
            "Selected nucleotide atoms: the nucleotide of (nt_chain, nt_name, "
            "nt_id) -- every copy when no id is given -- narrowed to the atom "
            "names in nt_atom_selected, applied to each copy."
        ),
    )
)
