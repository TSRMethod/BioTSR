"""Input checks a selection applies to the records it would use (ADR-108, ADR-109).

Each check **refuses** the calculation, naming the records, rather than
repairing, merging or dropping anything: a result computed from input BioTSR
cannot read reliably would look exactly like a correct one.  Each acts only on
the records a selection would actually use -- never on unrelated records
elsewhere in the file, so one damaged chain does not stop a calculation on
another.

* :func:`refuse_blank_elements` -- an atom with a blank element field (PDB
  columns 77-78) cannot be told apart from a hydrogen, and its label would be
  guessed from its name.  Not applied to Calpha points, whose ``HETATM`` CA may
  carry a blank element field (ADR-107).
* :func:`refuse_duplicate_atoms` -- two selected records for one atom (repeated
  or conflicting), or two distinct selected atoms at exactly the same
  coordinates.  No distance tolerance.
* :func:`refuse_unreadable_numbers` -- a record whose residue number cannot be
  read (hybrid-36) that the selection would include.
* :func:`refuse_unk` -- a residue named ``UNK`` inside the requested selection.
"""

from __future__ import annotations

from typing import Callable, Iterable, List, Optional, Sequence

from ..domain.errors import InvalidRequest

__all__ = [
    "UNK",
    "refuse_blank_elements",
    "refuse_duplicate_atoms",
    "refuse_coincident_points",
    "refuse_unreadable_numbers",
    "refuse_unk",
    "unk_residues",
]

#: The residue name the PDB uses for an unknown amino acid.
UNK = "UNK"

_SHOWN = 6


def _line(atom) -> str:
    return "line %d" % (atom.file_order + 1) if atom.file_order >= 0 else "a record"


def _more(items: Sequence[str]) -> str:
    shown = "; ".join(items[:_SHOWN])
    if len(items) > _SHOWN:
        shown += "; and %d more" % (len(items) - _SHOWN)
    return shown


def refuse_blank_elements(atoms: Iterable, context: str) -> None:
    """Refuse when any of ``atoms`` has a blank element field."""
    blank = [a for a in atoms if not (a.element or "").strip()]
    if not blank:
        return
    raise InvalidRequest(
        "%s: %d selected atom record(s) have a blank element field (PDB columns "
        "77-78): %s. Without it BioTSR cannot tell a hydrogen from a heavy atom, "
        "nor label the atom reliably, so the calculation is not computed. Write "
        "the element symbol of each atom into columns 77-78 (most structure "
        "tools can add it), or narrow the selection to leave these atoms out."
        % (context, len(blank),
           _more(["%s %s" % (_line(a), a.identity.describe()) for a in blank]))
    )


def _atom_key(atom):
    i = atom.identity
    # Where an atom is: model, chain, residue number and insertion code, atom
    # name. Neither the alternate-location indicator (after the blank-or-A rule,
    # a blank and an "A" record of one atom are two records of one atom) nor the
    # residue name (one residue number names one residue; two names for it are
    # a conflict, reported as one) is part of it.
    return (i.model, i.chain.value, i.residue_seq, i.insertion_code, i.atom_name)


_FIELDS = (("record type", 0, 6), ("altLoc", 16, 17), ("residue name", 17, 20),
           ("coordinates", 30, 54),
           ("occupancy", 54, 60), ("temperature factor", 60, 66),
           ("element", 76, 78), ("charge", 78, 80))


def _differences(a, b) -> List[str]:
    out = []
    for name, start, end in _FIELDS:
        if a.record_text[start:end].strip() != b.record_text[start:end].strip():
            out.append(name)
    return out


def refuse_duplicate_atoms(atoms: Sequence, context: str) -> None:
    """Refuse repeated, conflicting or coincident atom records.

    ``atoms`` are the parsed records a selection would use, after the model and
    alternate-location rules.  Three cases, all refused, none repaired:

    * the same atom written twice with every field equal (a repeated record);
    * one atom identity written twice with any field different (conflicting
      records -- BioTSR will not choose one);
    * two different atoms at exactly the same coordinates.
    """
    by_key = {}
    for atom in atoms:
        by_key.setdefault(_atom_key(atom), []).append(atom)
    repeated, conflicting = [], []
    for group in by_key.values():
        if len(group) < 2:
            continue
        first = group[0]
        for other in group[1:]:
            diff = _differences(first, other)
            text = "%s and %s (%s)" % (_line(first), _line(other),
                                       first.identity.describe())
            if diff:
                conflicting.append("%s differ in %s" % (text, ", ".join(diff)))
            else:
                repeated.append(text)
    if repeated or conflicting:
        parts = []
        if repeated:
            parts.append("repeated records of one atom: %s" % _more(repeated))
        if conflicting:
            parts.append("conflicting records of one atom: %s" % _more(conflicting))
        raise InvalidRequest(
            "%s: the selection contains %s. BioTSR does not merge, choose between "
            "or repair such records; correct the structure file so each atom is "
            "written once." % (context, "; ".join(parts))
        )
    by_coord = {}
    for atom in atoms:
        by_coord.setdefault(atom.coord.as_tuple(), []).append(atom)
    coincident = ["%s %s and %s %s" % (_line(g[0]), g[0].identity.describe(),
                                       _line(g[1]), g[1].identity.describe())
                  for g in by_coord.values() if len(g) > 1]
    if coincident:
        raise InvalidRequest(
            "%s: distinct selected atoms have exactly the same coordinates: %s. "
            "Two atoms cannot occupy one position, and BioTSR does not choose "
            "between them; correct the structure file." % (context, _more(coincident))
        )


def refuse_coincident_points(points: Sequence, context: str) -> None:
    """:func:`refuse_duplicate_atoms`' coordinate test, for labelled points
    gathered from two selections (the two sides of a cross)."""
    by_coord = {}
    for point in points:
        by_coord.setdefault(point.coord.as_tuple(), []).append(point)
    coincident = ["%s and %s" % (g[0].identity.describe(), g[1].identity.describe())
                  for g in by_coord.values() if len(g) > 1]
    if coincident:
        raise InvalidRequest(
            "%s: distinct selected atoms have exactly the same coordinates: %s. "
            "Two atoms cannot occupy one position, and BioTSR does not choose "
            "between them; correct the structure file." % (context, _more(coincident))
        )


def refuse_unreadable_numbers(structure, chain_id, cutoff: Optional[int],
                              accept: Callable, context: str) -> None:
    """Refuse when a record with an unreadable residue number would be selected.

    ``accept(record)`` says whether the selection's own record-type and name
    rules would take the record if its number could be read; ``cutoff`` is the
    selection's truncation line (records at or after it are not read).
    """
    found = [r for r in getattr(structure, "unreadable_numbers", ())
             if r.chain == chain_id
             and (cutoff is None or r.file_order < cutoff)
             and accept(r)]
    if not found:
        return
    raise InvalidRequest(
        "%s: %d record(s) the selection would include have a residue number "
        "BioTSR cannot read: %s. Hybrid-36 residue numbers (used past 9999 "
        "residues) are not decoded in BioTSR 1.0, and computing without those "
        "residues would report an incomplete selection as complete. Renumber "
        "the residues below 10000, or narrow the selection to leave them out."
        % (context, len(found), _more([r.describe() for r in found]))
    )


def unk_residues(records: Iterable) -> List[str]:
    """``NAME_ID`` tokens of the ``UNK`` residues among ``records``, in order."""
    out = []
    for atom in records:
        if atom.identity.residue_name != UNK:
            continue
        token = "UNK_%s" % atom.identity.residue.legacy_seq_label
        if token not in out:
            out.append(token)
    return out


def refuse_unk(tokens: Sequence[str], context: str) -> None:
    """Refuse a selection that includes ``UNK`` residues, naming them."""
    if not tokens:
        return
    raise InvalidRequest(
        "%s: the requested selection includes residue(s) named UNK (%s). UNK is "
        "an unknown residue: it has no residue label and no defined atom set, so "
        "no calculation is done for it, and BioTSR does not compute the rest of "
        "the selection without it. Leave it out by narrowing the selection "
        "(named residues, a residue range or residue types), or give the residue "
        "its real name in the structure file." % (context, _more(list(tokens)))
    )
