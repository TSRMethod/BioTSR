"""Which target arguments each method reads -- one table for the API and the CLI.

A method that is handed an argument it never reads has a choice between two
wrong answers: ignore it, and let a typo or a misunderstanding produce a result
that looks like what was asked for, or guess what was meant.  The package does
neither.  An argument the selected method does not read is refused, by name, in
the Python API and on the command line alike, and the refusal says what the
method does take.

Both layers read the tables below, so they cannot disagree.  The CLI adds
nothing scientific to them; it only spells the names with dashes.
"""

from __future__ import annotations

from typing import Any, Mapping, Tuple

from ..domain.errors import InvalidRequest

__all__ = [
    "SELECTION_ARGUMENTS",
    "CROSS_ARGUMENTS",
    "CROSS_ARGUMENTS_BY_METHOD",
    "TARGET_ARGUMENTS",
    "accepted_arguments",
    "required_arguments",
    "refuse_unused_arguments",
    "flag",
]

#: Which selection arguments each one-sided selector reads, and which of them it
#: requires: ``selector: (accepted, required)``.  A selector absent here takes no
#: selection argument at all.
SELECTION_ARGUMENTS: Mapping[str, Tuple[Tuple[str, ...], Tuple[str, ...]]] = {
    "protein_calpha": ((), ()),
    "protein_calpha_selection": (("residues",), ("residues",)),
    "amino_acid_residue_atoms": (
        ("residue", "residue_name"), ("residue", "residue_name"),
    ),
    "amino_acid_atom_selection": (("atoms",), ("atoms",)),
    "residue_type_atoms": (("residue_types",), ("residue_types",)),
    "specific_residue_atoms": (("residues",), ("residues",)),
    "entity_atoms": (("entity_name", "entity_id"), ("entity_name",)),
    "entity_atom_selection": (
        ("entity_name", "entity_id", "atoms"), ("entity_name", "atoms"),
    ),
    "entity_chain_atoms": (("entity_names",), ()),
}

#: Every cross-only argument, in the order ``analyze_cross`` documents them.
CROSS_ARGUMENTS: Tuple[str, ...] = (
    "chain_2",
    "contact_cutoff",
    "contact_filter",
    "no_contact_filter",
    "entity_name",
    "entity_id",
    "entity_names",
    "atoms",
    "residue_types",
    "residues",
    "atom_selection",
    "calpha_residues",
    "residue_types_2",
    "residues_2",
    "atom_selection_2",
    "calpha_residues_2",
)

_CROSS_CONTACT = ("chain_2", "contact_cutoff", "contact_filter", "no_contact_filter")
_CROSS_ENTITY = ("entity_name", "entity_id", "entity_names", "atoms")
_CROSS_ATOMIC_PROTEIN = ("residue_types", "residues", "atom_selection")

#: Which cross arguments each cross method reads.  Arguments routed by role: an
#: entity family (drug, carbohydrate, nucleotide) has one entity side and one
#: protein side; ``protein_protein_cross`` has two protein sides, the second
#: narrowed by the ``_2`` arguments; ``protein_protein_calpha_cross`` is a Cα
#: method, so its sides are narrowed only by residue-number selections, its pool
#: rule is symmetric by construction (no ``contact_filter``) and it is defined by
#: its contact cutoff (no ``no_contact_filter``).
CROSS_ARGUMENTS_BY_METHOD: Mapping[str, Tuple[str, ...]] = {
    "drug_protein_cross": _CROSS_CONTACT + _CROSS_ENTITY + _CROSS_ATOMIC_PROTEIN,
    "carbohydrate_protein_cross": (
        _CROSS_CONTACT + _CROSS_ENTITY + _CROSS_ATOMIC_PROTEIN
    ),
    "nucleotide_protein_cross": (
        _CROSS_CONTACT + _CROSS_ENTITY + _CROSS_ATOMIC_PROTEIN
    ),
    "protein_protein_cross": (
        _CROSS_CONTACT + _CROSS_ATOMIC_PROTEIN
        + ("residue_types_2", "residues_2", "atom_selection_2")
    ),
    "protein_protein_calpha_cross": (
        "chain_2", "contact_cutoff", "calpha_residues", "calpha_residues_2",
    ),
}

#: Every target argument any method reads.
TARGET_ARGUMENTS: Tuple[str, ...] = tuple(
    sorted(
        {name for accepted, _ in SELECTION_ARGUMENTS.values() for name in accepted}
        | set(CROSS_ARGUMENTS)
    )
)


def accepted_arguments(recipe) -> Tuple[str, ...]:
    """The target arguments ``recipe`` reads, in a stable order."""
    if recipe.family == "cross":
        return CROSS_ARGUMENTS_BY_METHOD[recipe.name]
    return SELECTION_ARGUMENTS.get(recipe.selector, ((), ()))[0]


def required_arguments(recipe) -> Tuple[str, ...]:
    if recipe.family == "cross":
        return ("chain_2",)
    return SELECTION_ARGUMENTS.get(recipe.selector, ((), ()))[1]


def flag(name: str) -> str:
    """The command-line spelling of an argument name."""
    return "--" + name.replace("_", "-")


def _is_set(value: Any) -> bool:
    if value is None or value is False:
        return False
    if isinstance(value, str) and not value.strip():
        return False
    if isinstance(value, (list, tuple)) and not value:
        return False
    return True


def refuse_unused_arguments(recipe, given: Mapping[str, Any], *, cli: bool = False,
                            prefix: str = "") -> None:
    """Raise :class:`InvalidRequest` naming every set argument ``recipe`` ignores.

    ``given`` maps argument names to values; unset values (``None``, ``False``,
    blank) are not arguments.  ``cli`` spells names as flags in the message.
    """
    accepted = set(accepted_arguments(recipe))
    unused = sorted(
        name for name, value in given.items()
        if name in TARGET_ARGUMENTS and name not in accepted and _is_set(value)
    )
    if not unused:
        return
    spell = flag if cli else (lambda n: "%s=" % n)
    takes = sorted(accepted)
    raise InvalidRequest(
        "%smethod %r does not use %s, so %s would be silently ignored; it is "
        "refused instead. %s"
        % (
            prefix,
            recipe.name,
            ", ".join(spell(n).rstrip("=") for n in unused),
            "it" if len(unused) == 1 else "they",
            ("It takes: %s." % ", ".join(spell(n).rstrip("=") for n in takes))
            if takes else
            "It takes no target argument beyond the chain.",
        )
    )
