"""Protein C-alpha selection -- the first selector BioTSR implemented.

    Structure -> chain -> protein residues -> Ca points

Deliberately **not** implemented here, and not stubbed either: amino-acid-type
selection, grouping, mirror, SSE annotation, size filtering, atomic selection,
ligand/carbohydrate/nucleotide selection, and cross-entity pools.  Each lives in
its own module or modifier.  Candidate-strategy extensions (for example CA-kNN)
are not part of this release.

Admission rules, taken from ``TSR-Package/tsr_package/TSR.py`` with one
extension (ADR-107):

* ``ATOM`` records; and ``HETATM`` records of the supported modified residues
  (``base.CALPHA_MODIFIED_RESIDUES``: TPO, SEP, PTR, MSE, NIY, M3L, TRF, CSO,
  LLP, CME, CSD) whose ``CA`` atom's element field is carbon or blank -- never
  any other ``HETATM`` record, so a ligand or an ion with an atom named ``CA``
  (calcium) is not a point, and a supported residue's ``CA`` written with
  another element is skipped and recorded (ADR-107).  No other atom of the
  residue is required: the calculation uses the Ca alone;
* atom name ``CA``, read from columns 14-16 (``calpha_name``) as every
  historical Calpha program read it, so calcium (``CA`` from column 13) never
  qualifies and neither does ``CA1`` (ADR-108);
* altLoc blank or ``A``;
* chain equal to the requested chain;
* everything before the legacy truncation point (``selection/truncation.py``).

A residue named ``UNK`` inside the selection refuses the calculation: UNK has
no residue label, and computing the rest without it would report an incomplete
selection as complete (ADR-109).  So do two selected ``CA`` records of one
residue, two Calpha atoms at identical coordinates, and a record whose residue
number cannot be read (ADR-108).

Protein Ca TSR has **no completeness policy** -- ``TSR.py`` has none -- which is
part of why it is the right family to implement first: it is the one path that
never touches that decision.
"""

from __future__ import annotations

import re
from typing import List, Tuple

from ..domain.errors import InvalidRequest
from ..domain.identity import ChainId
from ..domain.point import LabeledPoint, PointSet
from ..domain.result import SkipRecord
from ..labels.schemes import LabelScheme, UnlabelledPointPolicy
from ..structure.model import Structure
from . import integrity
from .base import (
    CALPHA_MODIFIED_RESIDUES,
    AltLocPolicy,
    SelectionResult,
    SeqLabelPolicy,
)
from .truncation import LEGACY_SINGLE_PASS, legacy_truncation_line

__all__ = [
    "SELECTOR_NAME",
    "SELECTOR_CALPHA_SELECTION",
    "LABEL_BY_RESIDUE_IDENTITY",
    "SCOPE",
    "calpha_exclusions",
    "calpha_element_admissible",
    "select_calpha",
    "select_calpha_residues",
    "parse_residue_selection",
    "parse_residue_ranges",
    "check_range_insertion_codes",
]

SELECTOR_NAME = "protein_calpha"
#: ``aa_selection`` -- the Ca atoms of an explicitly listed set of residues.
SELECTOR_CALPHA_SELECTION = "protein_calpha_selection"
SCOPE = "pooled"

#: The standard reading of a protein Ca point: the residue's own identity, in the
#: residue label space (ADR-032, ADR-107).
LABEL_BY_RESIDUE_IDENTITY = "residue_identity"

#: Residue names no Calpha calculation includes: a selection containing one
#: is refused (ADR-109).
EXCLUDED_RESIDUE_NAMES = (integrity.UNK,)

_ACCEPTED_ALTLOC = ("", "A")


#: The element fields a Ca atom may carry (ADR-107): carbon, or none at all --
#: files written before the PDB element column existed leave it blank. An
#: explicit other element (``CA``, calcium) is never a Ca point.
_CALPHA_ELEMENTS = ("C", "")

#: Why a residue present in the chain has no Ca point (``calpha_exclusions``).
REASON_NOT_LABELLED = (
    "a modified residue stored as HETATM; the residue labels of the Calpha "
    "methods cover the 20 regular amino acids and only the modified residues "
    "%s" % ", ".join(CALPHA_MODIFIED_RESIDUES))
REASON_NOT_A_RESIDUE = "a HETATM record that is not an amino-acid residue the Calpha methods read"
REASON_UNK = "residue name UNK: no calculation is done for UNK"
REASON_NO_CA = "no CA atom"
REASON_CA_ALTLOC = ("its CA atom has only alternate-location indicator(s) %s; "
                    "blank or A is read")
REASON_CA_ELEMENT = "its CA atom's element field is %s, not carbon"


def calpha_element_admissible(atom) -> bool:
    """Whether a ``HETATM`` Ca atom's element field allows it as a point."""
    return atom.element.strip().upper() in _CALPHA_ELEMENTS


def _hetero_rejected(atom, chain_id) -> bool:
    """Whether a ``HETATM`` record is outside the Ca selection outright.

    Only records of a supported modified residue (``CALPHA_MODIFIED_RESIDUES``)
    of the requested chain go on to the Ca, alternate-location and element
    checks; every other ``HETATM`` record -- a ligand, an ion, a modified
    residue the labels do not cover -- is passed over, as before (ADR-107).
    """
    return atom.chain != chain_id or atom.residue_name not in CALPHA_MODIFIED_RESIDUES


def _element_skip(atom) -> SkipRecord:
    return SkipRecord(
        stage="selection",
        reason="calpha_element_not_carbon",
        detail="%s (element %r)" % (atom.identity.describe(), atom.element.strip()),
    )


def _calpha_record(record, chain_id) -> bool:
    """Whether a Calpha selection reads records like ``record`` at all."""
    return record.chain == chain_id and (
        not record.hetero or record.residue_name in CALPHA_MODIFIED_RESIDUES)


def _check_calpha_input(structure: Structure, chain: str, cutoff, chosen,
                        unk_tokens, refuse_unk: bool) -> None:
    """The checks a Calpha selection applies to what it would use (ADR-108,
    ADR-109): UNK, unreadable residue numbers, duplicate Calpha records."""
    context = "chain %s of %s" % (chain, structure.structure_id.value)
    if refuse_unk:
        integrity.refuse_unk(unk_tokens, context)
    chain_id = ChainId(chain)
    integrity.refuse_unreadable_numbers(
        structure, chain_id, cutoff, lambda r: _calpha_record(r, chain_id), context)
    integrity.refuse_duplicate_atoms(chosen, context)


def calpha_exclusions(structure: Structure, chain: str, admitted, requested=None):
    """``[(seq_label, residue_name, reason)]``: residues present in the chain's
    legacy read that have no Ca point, and why.

    With ``requested`` (residue identifiers, as ``aa_selection`` names them),
    every requested residue that is present but was not admitted; one that is
    absent from the chain is not listed -- the caller reports it as not found.
    Without it (a whole chain), the supported modified residues that were not
    admitted, and the ``HETATM`` amino acids the residue labels do not cover.
    ``admitted`` is the selection's ``residues_used``. File order.
    """
    from .atomic import hetero_amino_acid_residues

    chain_id = ChainId(chain)
    cutoff = legacy_truncation_line(structure, chain_id)
    present = {}
    order = []
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.chain != chain_id:
            continue
        seq = atom.identity.residue.legacy_seq_label
        entry = present.get(seq)
        if entry is None:
            entry = present[seq] = {"name": atom.residue_name, "hetero": atom.hetero,
                                    "ca": []}
            order.append(seq)
        if atom.calpha_name == "CA" and atom.residue_name == entry["name"]:
            entry["ca"].append(atom)
    amino_acids = {seq for seq, _ in
                   hetero_amino_acid_residues(structure, chain, LEGACY_SINGLE_PASS)}
    admitted = set(admitted)
    wanted = order if requested is None else [str(s) for s in requested]
    out = []
    for seq in wanted:
        if seq in admitted or seq not in present:
            continue
        entry = present[seq]
        name = entry["name"]
        if entry["hetero"] and name not in CALPHA_MODIFIED_RESIDUES:
            if requested is None and seq not in amino_acids:
                continue  # a ligand or an ion: not part of the chain's residues
            reason = REASON_NOT_LABELLED if seq in amino_acids else REASON_NOT_A_RESIDUE
        elif requested is None and name not in CALPHA_MODIFIED_RESIDUES:
            continue  # a regular residue without a Ca: as before, not named
        elif name in EXCLUDED_RESIDUE_NAMES:
            reason = REASON_UNK
        elif not entry["ca"]:
            reason = REASON_NO_CA
        else:
            readable = [a for a in entry["ca"] if a.identity.alt_loc in _ACCEPTED_ALTLOC]
            if not readable:
                reason = REASON_CA_ALTLOC % ", ".join(
                    sorted({a.identity.alt_loc for a in entry["ca"]}))
            elif entry["hetero"] and not any(calpha_element_admissible(a) for a in readable):
                reason = REASON_CA_ELEMENT % ", ".join(
                    sorted({repr(a.element.strip()) for a in readable}))
            else:
                continue
        out.append((seq, name, reason))
    return out


def select_calpha(
    structure: Structure,
    chain: str,
    label_scheme: LabelScheme,
    *,
    altloc_policy: str = AltLocPolicy.BLANK_OR_A,
    truncation: str = LEGACY_SINGLE_PASS,
    unlabelled_policy: str = UnlabelledPointPolicy.ERROR,
    seq_label_policy: str = SeqLabelPolicy.ACCEPT_ALL,
    refuse_unk: bool = True,
) -> SelectionResult:
    """Select the labelled C-alpha point cloud for one chain.

    Points come out in **file order**, which is what fixes deterministic
    candidate enumeration and therefore deterministic serialization.

    Every residue is admitted whatever its sequence label looks like
    (``SeqLabelPolicy.ACCEPT_ALL``), insertion-code residues included.  A
    ``UNK`` residue anywhere in the chain's read range refuses the selection
    (``refuse_unk=False`` only for reporting what a chain holds).
    """
    if altloc_policy != AltLocPolicy.BLANK_OR_A:
        raise ValueError(
            "unsupported altloc policy %r; BioTSR implements only %r"
            % (altloc_policy, AltLocPolicy.BLANK_OR_A)
        )
    if truncation != LEGACY_SINGLE_PASS:
        raise ValueError(
            "unsupported truncation rule %r; BioTSR implements only %r"
            % (truncation, LEGACY_SINGLE_PASS)
        )
    if seq_label_policy != SeqLabelPolicy.ACCEPT_ALL:
        raise ValueError(
            "unknown sequence-label policy %r; BioTSR implements only %r"
            % (seq_label_policy, SeqLabelPolicy.ACCEPT_ALL)
        )

    chain_id = ChainId(chain)
    cutoff = legacy_truncation_line(structure, chain_id)
    points: List[LabeledPoint] = []
    skipped: List[SkipRecord] = []
    residues_used: List[str] = []
    chains_seen: List[str] = []
    chosen = []

    unk = integrity.unk_residues(
        a for a in structure.atoms
        if (cutoff is None or a.file_order < cutoff) and _calpha_record(a, chain_id))
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.hetero and _hetero_rejected(atom, chain_id):
            continue
        if atom.calpha_name != "CA":
            continue
        if atom.identity.alt_loc not in _ACCEPTED_ALTLOC:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="altloc_rejected",
                    detail="%s (altLoc=%r, policy=%s)"
                    % (atom.identity.describe(), atom.identity.alt_loc, altloc_policy),
                )
            )
            continue
        if atom.chain != chain_id:
            continue
        if atom.hetero and not calpha_element_admissible(atom):
            skipped.append(_element_skip(atom))
            continue
        if atom.residue_name in EXCLUDED_RESIDUE_NAMES:
            # Only reachable with refuse_unk=False: reported, never a point.
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="excluded_residue_name",
                    detail=atom.identity.describe(),
                )
            )
            continue
        chosen.append(atom)

    _check_calpha_input(structure, chain, cutoff, chosen, unk, refuse_unk)
    for atom in chosen:
        label = label_scheme.label(atom.residue_name)
        if label is None:
            if unlabelled_policy == UnlabelledPointPolicy.ERROR:
                raise InvalidRequest(
                    "residue %s %s of chain %s cannot be labelled: the residue "
                    "labels %s have no code for %r (%s). A Calpha calculation is "
                    "not computed without a residue it selected; correct the "
                    "residue name in the structure file, or list the other "
                    "residues with aa_selection."
                    % (atom.residue_name, atom.identity.residue.legacy_seq_label,
                       chain, label_scheme.name, atom.residue_name,
                       atom.identity.describe())
                )
            skipped.append(
                SkipRecord(
                    stage="labelling",
                    reason="unlabelled_point",
                    detail="%s has no code in %s"
                    % (atom.identity.describe(), label_scheme.name),
                )
            )
            continue

        points.append(LabeledPoint(identity=atom.identity, coord=atom.coord, label=label))
        residues_used.append(atom.identity.residue.legacy_seq_label)
        if atom.chain.value not in chains_seen:
            chains_seen.append(atom.chain.value)

    return SelectionResult(
        points=PointSet(points),
        selector=SELECTOR_NAME,
        chains_requested=(chain,),
        chains_used=tuple(chains_seen),
        residues_used=tuple(residues_used),
        skipped=tuple(skipped),
        truncated_at_line=-1 if cutoff is None else cutoff,
    )


# ---------------------------------------------------------------------------
# aa_selection -- Ca atoms of explicitly listed residues
# ---------------------------------------------------------------------------


#: ``10-12`` or ``-5--3``: two plain residue numbers, optionally negative.
_RANGE = re.compile(r"^\s*(-?\d+)\s*-\s*(-?\d+)\s*$")


def _tokens(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        raw = [str(v) for v in value]
    else:
        raw = str(value).replace(",", ";").replace("\n", ";").split(";")
    return [t.strip() for t in raw if t.strip()]


def parse_residue_ranges(value) -> Tuple[Tuple[int, int], ...]:
    """The numeric ranges of an ``aa_selection`` cell, as ``(first, last)`` pairs."""
    out = []
    for token in _tokens(value):
        match = _RANGE.match(token)
        if match:
            out.append((int(match.group(1)), int(match.group(2))))
    return tuple(out)


def parse_residue_selection(value):
    """Expand an ``aa_selection`` cell into residue identifiers.

    The grammar, items separated by ``;`` (or ``,``)::

        48          one residue, by its PDB number
        52A         one residue with an insertion code
        -5          a negative residue number
        10-12       a range: plain numbers, ascending, both ends included,
                    expanding to 10, 11, 12 -- never to 11A (see below)
        -5--3       a range of negative numbers

    ``"10-12;48;52A"`` -> ``("10", "11", "12", "48", "52A")``.  Identifiers are
    normalised as :func:`biotsr.selection.atomic.parse_residue_id` does
    (``007`` is ``7``; an insertion code is kept as written and matched
    exactly).  Duplicates are dropped, preserving the first occurrence; the
    order of the list never matters -- points follow the structure file.

    Refused with a message naming the token, never kept verbatim or skipped: a
    range with an insertion code at either end (``27A-29``), a descending range
    (``12-10``), and anything else that is neither a residue identifier nor a
    numeric range.  A range spanning a residue that has an insertion code is
    checked against the structure by :func:`check_range_insertion_codes`.
    ``warnings`` is always empty.
    """
    from .atomic import parse_residue_id

    labels = []
    for token in _tokens(value):
        match = _RANGE.match(token)
        if match:
            first, last = int(match.group(1)), int(match.group(2))
            if first > last:
                raise InvalidRequest(
                    "residue range %r is descending: write it first-last, e.g. "
                    "%d-%d" % (token, last, first))
            expanded = [str(n) for n in range(first, last + 1)]
        elif "-" in token.lstrip("-"):
            raise InvalidRequest(
                "residue range %r is not supported: a range is two plain residue "
                "numbers, ascending (10-12, or -5--3 for negative numbers); a "
                "residue with an insertion code is listed on its own (52A), "
                "e.g. 10-12;12A" % token)
        else:
            expanded = [parse_residue_id(
                token, "residue %r in the residue selection" % token)]
        for label in expanded:
            if label not in labels:
                labels.append(label)
    return tuple(labels), ()


def check_range_insertion_codes(structure: Structure, chain: str, value,
                                labels) -> None:
    """Refuse a range that spans a residue with an insertion code.

    A numeric range names plain residue numbers, so ``50-53`` never includes
    ``52A``.  Rather than leave such a residue out silently, the selection is
    refused unless the residue is listed explicitly (``50-53;52A`` includes
    it); to leave it out, list the residues one by one (``50;51;52;53``).
    """
    ranges = parse_residue_ranges(value)
    if not ranges:
        return
    chain_id = ChainId(chain)
    listed = set(labels)
    spanned = []
    for atom in structure.atoms:
        residue = atom.identity.residue
        if atom.chain != chain_id or not residue.insertion_code:
            continue
        label = residue.legacy_seq_label
        if label in listed or label in spanned:
            continue
        for first, last in ranges:
            if first <= residue.seq <= last:
                spanned.append(label)
                break
    if spanned:
        raise InvalidRequest(
            "the residue range(s) %s of chain %s span residue(s) %s, which have "
            "an insertion code: a numeric range names plain residue numbers "
            "only. List them to include them (e.g. %s;%s), or list the residues "
            "one by one to leave them out"
            % (", ".join("%d-%d" % r for r in ranges), chain, ", ".join(spanned),
               "%d-%d" % ranges[0], ";".join(spanned)))


def select_calpha_residues(
    structure: Structure,
    chain: str,
    residue_seqs,
    label_scheme: LabelScheme,
    *,
    altloc_policy: str = AltLocPolicy.BLANK_OR_A,
    truncation: str = LEGACY_SINGLE_PASS,
    unlabelled_policy: str = UnlabelledPointPolicy.ERROR,
    label_key: str = LABEL_BY_RESIDUE_IDENTITY,
    refuse_unk: bool = True,
) -> SelectionResult:
    """Select the Ca atoms of an explicitly listed set of residues.

    This is the ``aa_selection`` column's selector, and it belongs to the **Ca
    family**: one point per residue, labelled by *residue identity* in the
    residue label space (ADR-032).  It admits records exactly as
    :func:`select_calpha` does, so a listed supported modified residue written
    as ``HETATM`` is selected too (ADR-107).  It is emphatically **not** an atomic method that
    happens to pick atoms named ``CA``.

    That distinction is load-bearing: the historical ``WebsiteTSR`` selector
    looked ``"CA"`` up in the 112-code atom lexicon, so every point carried the
    identical label 106 and the label axis carried no information.  BioTSR
    labels a Ca point by its residue's identity.

    Residue **names are not matched** -- the ``aa_selection`` column carries only
    numbers, so chain plus sequence label is the identity, which is what
    ``WebsiteTSR`` does too.  A requested residue with no Ca is reported.
    """
    if altloc_policy != AltLocPolicy.BLANK_OR_A:
        raise ValueError(
            "unsupported altloc policy %r; this package implements only %r"
            % (altloc_policy, AltLocPolicy.BLANK_OR_A)
        )
    if truncation != LEGACY_SINGLE_PASS:
        raise ValueError(
            "unsupported truncation rule %r; Ca selection implements only %r"
            % (truncation, LEGACY_SINGLE_PASS)
        )

    chain_id = ChainId(chain)
    cutoff = legacy_truncation_line(structure, chain_id)
    wanted = [str(s).strip() for s in residue_seqs]
    wanted_set = set(wanted)
    if not wanted_set:
        raise ValueError("aa_selection produced no residue sequence labels")

    points: List[LabeledPoint] = []
    skipped: List[SkipRecord] = []
    residues_used: List[str] = []
    seen = set()
    chosen = []

    unk = integrity.unk_residues(
        a for a in structure.atoms
        if (cutoff is None or a.file_order < cutoff) and _calpha_record(a, chain_id)
        and a.identity.residue.legacy_seq_label in wanted_set)
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.hetero and _hetero_rejected(atom, chain_id):
            continue
        if atom.calpha_name != "CA":
            continue
        if atom.chain != chain_id:
            continue
        # The requested-residue filter comes first, so an altLoc rejection is only
        # ever recorded for a residue the caller actually asked for. Recording
        # every altLoc B in the chain would bury the two that matter.
        seq_label = atom.identity.residue.legacy_seq_label
        if seq_label not in wanted_set:
            continue
        if atom.identity.alt_loc not in _ACCEPTED_ALTLOC:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="altloc_rejected",
                    detail="%s (altLoc=%r, policy=%s)"
                    % (atom.identity.describe(), atom.identity.alt_loc, altloc_policy),
                )
            )
            continue
        if atom.hetero and not calpha_element_admissible(atom):
            skipped.append(_element_skip(atom))
            continue
        if atom.residue_name in EXCLUDED_RESIDUE_NAMES:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="excluded_residue_name",
                    detail=atom.identity.describe(),
                )
            )
            continue
        chosen.append(atom)

    _check_calpha_input(structure, chain, cutoff, chosen, unk, refuse_unk)
    for atom in chosen:
        seq_label = atom.identity.residue.legacy_seq_label
        if label_key != LABEL_BY_RESIDUE_IDENTITY:
            raise ValueError(
                "unsupported label key %r; a Calpha point is labelled by its "
                "residue identity (%r)" % (label_key, LABEL_BY_RESIDUE_IDENTITY)
            )
        key_text = atom.residue_name
        label = label_scheme.label(key_text)
        if label is None:
            if unlabelled_policy == UnlabelledPointPolicy.ERROR:
                raise InvalidRequest(
                    "residue %s %s of chain %s cannot be labelled: the residue "
                    "labels %s have no code for %r (%s). A Calpha calculation is "
                    "not computed without a residue it selected; correct the "
                    "residue name in the structure file, or leave the residue "
                    "out of the selection."
                    % (key_text, seq_label, chain, label_scheme.name, key_text,
                       atom.identity.describe())
                )
            skipped.append(
                SkipRecord(
                    stage="labelling",
                    reason="unlabelled_point",
                    detail="%s has no code in %s"
                    % (atom.identity.describe(), label_scheme.name),
                )
            )
            continue

        seen.add(seq_label)
        points.append(LabeledPoint(identity=atom.identity, coord=atom.coord, label=label))
        residues_used.append(seq_label)

    missing = tuple(s for s in wanted if s not in seen)
    for label in missing:
        skipped.append(
            SkipRecord(
                stage="selection",
                reason="requested_residue_not_found",
                detail="residue %s was requested but has no admitted Ca in chain %s"
                % (label, chain),
            )
        )

    return SelectionResult(
        points=PointSet(points),
        selector=SELECTOR_CALPHA_SELECTION,
        chains_requested=(chain,),
        chains_used=(chain,) if points else (),
        residues_used=tuple(residues_used),
        atoms_missing=missing,
        skipped=tuple(skipped),
        truncated_at_line=-1 if cutoff is None else cutoff,
    )
