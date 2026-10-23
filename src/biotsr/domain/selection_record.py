"""What a result's selection was asked for and what it actually selected.

Every calculation records, per side, a structured description of its selection
in the result's provenance (``provenance.selection_request``):

* ``selector`` and ``role`` -- how the side selects, and what it is (protein,
  drug, carbohydrate, nucleotide);
* ``chain`` -- the chain the caller named; a chain is never inferred;
* ``scope`` -- which of the documented selection scopes applied (see
  :data:`SCOPES`), so "this one NAG" and "every NAG of chain D" are never
  confused, even when the chain holds one copy;
* ``arguments`` -- the selection arguments as given, normalised (residue
  types and entity names upper-case, residue identifiers as PDB numbers);
* ``instances`` -- for the scopes that name molecules or residues, every
  instance that contributed atoms, as ``NAME_ID`` tokens in structure order,
  and ``excluded`` -- the requested instances that contributed nothing;
* ``not_found``, ``unsupported`` and ``declined`` -- present only when
  non-empty: what the request named that the structure does not contain,
  what it contains but the method cannot represent (a modified residue
  stored as ``HETATM`` other than TPO, SEP and PTR, or any of them in a
  Calpha method; ADR-101), and what is there but the method's completeness
  rule rejected (a residue or molecule missing atoms);
* ``atom_subset`` -- when named atoms narrow a molecule: the atoms, that the
  list applies to **each** selected instance, and which requested atoms each
  instance lacks.

The record is provenance.  It enters neither fingerprint and it never decides
whether two results may be compared: different structures, chains, residue
types, ligands, carbohydrates, residue numbers and atom subsets are what a
dataset is made of (ADR-054, ADR-094).  It is recorded so that every report can
say exactly what each sample used.

**A partial result** is one whose selection did not get everything it asked
for: :func:`selection_omissions` reads that from the record (and, for a result
written before these fields existed, from its warnings).  A partial result is
kept and is usable, but it is never presented as the complete selection: the
batch marks it ``partial`` and every analysis asks for it to be acknowledged
(ADR-106).
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Optional, Sequence

__all__ = [
    "SCOPES",
    "SCOPE_WHOLE_CHAIN",
    "SCOPE_CHAIN_VOCABULARY",
    "SCOPE_RESIDUE_LIST",
    "SCOPE_ONE_RESIDUE",
    "SCOPE_EACH_RESIDUE",
    "SCOPE_RESIDUE_TYPES",
    "SCOPE_NAMED_RESIDUES",
    "SCOPE_NAMED_ATOMS",
    "SCOPE_ONE_INSTANCE",
    "SCOPE_ALL_COPIES",
    "ARGUMENTS",
    "INCOMPLETE_WARNINGS",
    "STRUCTURE_NOTES",
    "OMISSION_KINDS",
    "side_record",
    "selection_record",
    "selection_sides",
    "selection_omissions",
]

#: Every atom of the chain the method admits (``protein``, an un-narrowed
#: cross protein side).
SCOPE_WHOLE_CHAIN = "whole_chain"
#: Every residue of the chain whose name is in a vocabulary
#: (``carbohydrate_chain``, ``nucleotide_chain``, ``entity_names=``): an entire
#: glycan or strand, different residue types included.
SCOPE_CHAIN_VOCABULARY = "chain_vocabulary"
#: Calpha atoms of listed residue numbers and ranges (``aa_selection``).
SCOPE_RESIDUE_LIST = "residue_list"
#: One amino-acid residue, named by its identifier (``amino_acid_atomic``).
SCOPE_ONE_RESIDUE = "one_residue"
#: One calculation per residue of the chain, optionally only residues of one
#: name (``amino_acid_atomic`` without ``residue``).
SCOPE_EACH_RESIDUE = "each_residue"
#: Every residue of the listed amino-acid types, pooled (``aa_ind``).
SCOPE_RESIDUE_TYPES = "residue_types"
#: The listed residue occurrences, each named by name and identifier.
SCOPE_NAMED_RESIDUES = "named_residues"
#: The listed atoms, each named by residue, identifier and atom.
SCOPE_NAMED_ATOMS = "named_atoms"
#: One molecule, named by its name and residue identifier.
SCOPE_ONE_INSTANCE = "one_instance"
#: Every molecule of one name on the chain, in one combined calculation
#: (no identifier given).  Triangles may span the copies, by design.
SCOPE_ALL_COPIES = "all_copies"

SCOPES = (
    SCOPE_WHOLE_CHAIN,
    SCOPE_CHAIN_VOCABULARY,
    SCOPE_RESIDUE_LIST,
    SCOPE_ONE_RESIDUE,
    SCOPE_EACH_RESIDUE,
    SCOPE_RESIDUE_TYPES,
    SCOPE_NAMED_RESIDUES,
    SCOPE_NAMED_ATOMS,
    SCOPE_ONE_INSTANCE,
    SCOPE_ALL_COPIES,
)

#: The selection arguments a record may carry.
ARGUMENTS = (
    "residue",
    "residue_name",
    "residues",
    "residue_types",
    "atoms",
    "atom_selection",
    "entity_name",
    "entity_id",
    "entity_names",
)

#: What a side records it did not use, each only when non-empty.
OMISSION_KINDS = ("not_found", "unsupported", "declined")

#: The warnings a result carries when its selection did not get everything it
#: asked for (``biotsr.api`` writes them).  Read only for results whose record
#: predates :data:`OMISSION_KINDS`; the structured record decides otherwise.
INCOMPLETE_WARNINGS = ("requested but not found:", "requested but excluded",
                       "not supported:", "not represented:")

#: The start of warnings that describe the structure file rather than the
#: selection: shown and kept with the result like the ones above, but never a
#: reason to call it partial -- every atom present was used (ADR-114).
STRUCTURE_NOTES = ("atom sets that differ from the standard nucleotide",)

_LISTS = ("residues", "residue_types", "atoms", "atom_selection", "entity_names")
_UPPER = ("residue_types", "entity_names", "entity_name", "residue_name")


def _names(value: Any) -> Sequence[str]:
    if value is None:
        return ()
    if isinstance(value, str):
        value = value.split(";")
    return tuple(str(v).strip() for v in value if str(v).strip())


def side_record(
    *,
    side: str,
    selector: str,
    chain: str,
    scope: str,
    role: str = "",
    instances: Iterable[str] = (),
    excluded: Iterable[str] = (),
    atom_subset_missing: Optional[Mapping[str, Sequence[str]]] = None,
    not_found: Iterable[str] = (),
    unsupported: Iterable[str] = (),
    declined: Iterable[str] = (),
    **arguments: Any,
) -> Mapping[str, Any]:
    """One side's selection record.

    ``arguments`` are the values actually used (parsed, resolved); unknown
    names and unknown scopes are refused, so a new selection argument cannot
    go unrecorded by accident.
    """
    if scope not in SCOPES:
        raise ValueError("unknown selection scope %r" % scope)
    unknown = sorted(set(arguments) - set(ARGUMENTS))
    if unknown:
        raise ValueError("selection argument(s) %s are not recordable"
                         % ", ".join(unknown))
    recorded = {}
    for name in ARGUMENTS:
        value = arguments.get(name)
        if name in _LISTS:
            names = _names(value)
            if names:
                if name in _UPPER:
                    names = tuple(n.upper() for n in names)
                recorded[name] = list(names)
        elif value not in (None, ""):
            recorded[name] = str(value).upper() if name in _UPPER else str(value)
    record = {
        "side": side,
        "role": role,
        "selector": selector,
        "chain": chain,
        "scope": scope,
        "arguments": recorded,
        "instances": list(instances),
        "excluded": list(excluded),
    }
    for name, values in (("not_found", not_found), ("unsupported", unsupported),
                         ("declined", declined)):
        values = list(values)
        if values:
            record[name] = values
    atoms = _names(arguments.get("atoms"))
    if atoms and scope in (SCOPE_ONE_INSTANCE, SCOPE_ALL_COPIES):
        record["atom_subset"] = {
            "atoms": list(atoms),
            "applies_to": "each_instance",
            "missing": {k: list(v) for k, v in
                        sorted((atom_subset_missing or {}).items())},
        }
    return record


def selection_record(sides: Iterable[Mapping[str, Any]]) -> Mapping[str, Any]:
    """The provenance record of a whole calculation: one entry per side."""
    return {"sides": [dict(s) for s in sides]}


def selection_sides(
        record: Optional[Mapping[str, Any]]) -> Optional[Sequence[Mapping[str, Any]]]:
    """The per-side entries of ``record``, or ``None`` when there is none."""
    if not record:
        return None
    return [dict(s) for s in record.get("sides", ())]


def selection_omissions(sides: Optional[Iterable[Mapping[str, Any]]],
                        warnings: Iterable[str] = ()) -> Dict[str, Any]:
    """What a calculation asked for and did not use; ``{}`` when it used it all.

    Per side (keyed by the side's name): ``not_found``, ``unsupported`` and
    ``declined`` from the record, and ``atoms_missing`` -- the requested atoms
    each selected instance lacks, when named atoms narrow a molecule.
    ``messages`` lists the result's warnings that say the same in words.  A
    result written before the structured fields existed has only those
    warnings, and they alone then make it partial: an older record is never
    read as complete because it cannot say otherwise.
    """
    out: Dict[str, Any] = {}
    for side in sides or ():
        entry: Dict[str, Any] = {}
        for kind in OMISSION_KINDS:
            if side.get(kind):
                entry[kind] = list(side[kind])
        missing = ((side.get("atom_subset") or {}).get("missing") or {})
        missing = {k: list(v) for k, v in sorted(missing.items()) if v}
        if missing:
            entry["atoms_missing"] = missing
        if entry:
            out[str(side.get("side") or "primary")] = entry
    messages = [w for w in warnings if str(w).startswith(INCOMPLETE_WARNINGS)]
    if messages:
        out["messages"] = messages
    return out
