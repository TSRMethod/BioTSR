"""Turning ``HELIX``/``SHEET`` records into a residue-level annotation.

This is the only place the two PDB records are *interpreted*, and it is where
every convention the interpretation rests on is named rather than assumed:

``STRAND_IDENTITY_SHEET_ID_AND_SERIAL``
    What identifies the **individual β-strand** a sheet residue belongs to.  The
    scientific object the categories distinguish is the individual strand, not
    the whole β-sheet: a β-sheet contains one or more strands, and categories
    7-12, 15 and 16 are statements about strands (ADR-063, closed).

    A ``SHEET`` record does not name a strand with one field.  ``serNum`` is the
    strand's index *within its sheet*, so it restarts at 1 in every sheet and a
    structure with two sheets has two strands numbered 1; ``sheetID`` names the
    sheet, so on its own it merges every strand of a sheet into one group.  The
    standard therefore takes **both**: ``(sheetID, serNum)``, carried
    componentwise as :class:`~biotsr.sse.states.SSEElementId` so no concatenation
    can manufacture a collision.

    The historical ``SSE_TSR.py`` used ``serNum`` alone, which is exact inside
    one sheet and conflates strands across sheets; BioTSR does not reproduce it.

``HELIX_IDENTITY_BY_SERIAL``
    A ``HELIX`` record's ``serNum`` is unique across the file, not within some
    container, so it already identifies an individual helix and nothing needs
    disambiguating.  The alternative field (``helixID``) was measured against it
    over the whole frozen corpus and changes no category; helix handling is
    therefore unchanged, by evidence rather than by assumption.

``OVERLAP_FIRST_RECORD_WINS``
    What happens when two records of the same kind cover one residue.  Legacy
    guards each insertion with ``if pos not in dict``, so the **first record in
    file order** wins.

``HELIX_BEFORE_SHEET``
    What happens when a residue is in both a ``HELIX`` and a ``SHEET`` range.
    Legacy tests the helix dict first (``if ... elif ...``), so helix wins.

``RANGE_BY_ENDPOINTS_IN_CHAIN_ORDER``
    Which residues a record covers (ADR-110): its start and end residues are
    found in the chain by chain identifier, residue number **and insertion
    code**, and the record covers every residue of the chain from the start
    residue to the end residue in the chain's own order -- the order its
    records are written in the file, up to the chain's ``TER``.  So a helix
    from ``14C`` to ``14I`` covers ``14C``..``14I`` and not residue 14, and a
    record from ``183`` to ``185`` covers ``184A`` when the file places it
    there.  When a relevant endpoint cannot be resolved -- not found in the
    chain, found in two separate places, written with a residue number that
    cannot be read, on another chain, or ending before it starts -- the
    chain's SSE annotation is **refused**, naming the record: no numeric
    approximation, and no residue silently becomes coil.

    (The historical ``SSE_TSR.py`` read a record as the numeric range
    ``range(start, stop + 1)``, ignoring insertion codes, so a residue carrying
    one was never covered and a record could cover a residue outside it; until
    1.0.0rc5 BioTSR did the same.)

Nothing in this module knows about triangles, keys, geometry or ordering.  It
maps residues to states; :mod:`biotsr.sse.categories` maps three states to a
category; the engine does the rest with code it already had.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from ..domain.errors import InvalidRequest
from ..domain.identity import AtomIdentity, ChainId
from ..structure.model import SecondaryStructureRecord, Structure
from .categories import classify
from .states import HELIX, NONE, SHEET, SSEAssignment, SSEElementId, UNASSIGNED

__all__ = [
    "STRAND_IDENTITY_SHEET_ID_AND_SERIAL",
    "STRAND_IDENTITIES",
    "HELIX_IDENTITY_BY_SERIAL",
    "SSE_ASSIGNMENT_VERSION",
    "ResidueKey",
    "SSEAssignmentMap",
    "SSEAnnotator",
    "build_assignment_map",
]

#: **The standard.** An individual β-strand is ``(sheetID, serNum)``: the sheet
#: it belongs to and its index within that sheet.  Neither component identifies a
#: strand alone -- ``serNum`` restarts in every sheet, and ``sheetID`` names a
#: sheet, which contains one or more strands.
STRAND_IDENTITY_SHEET_ID_AND_SERIAL = "sheet_id_and_serial"

#: The defined convention.  ``sheetID`` alone is deliberately absent: it would
#: make every strand of one sheet the same group, which is the opposite of what
#: the categories distinguish (ADR-063).
STRAND_IDENTITIES = (STRAND_IDENTITY_SHEET_ID_AND_SERIAL,)

#: A helix is identified by its ``HELIX`` record's ``serNum``, under **both**
#: conventions.  That serial is unique across the file rather than within some
#: container, so unlike a strand's it needs no disambiguation -- measured, not
#: assumed: ``helixID`` changes no category anywhere in the frozen corpus.
HELIX_IDENTITY_BY_SERIAL = "helix_serial"

#: Bumped when the assignment rules change, so provenance records which set
#: produced an annotation.  ``v2`` is the corrected β-strand identity; ``v3``
#: matches range endpoints with their insertion codes, in chain order (ADR-110).
SSE_ASSIGNMENT_VERSION = "sse_assignment_v3"

#: How a residue is addressed: chain, sequence number, insertion code.  The
#: insertion code is a separate component, never concatenated.
ResidueKey = Tuple[str, int, str]


def _residue_key(identity: AtomIdentity) -> ResidueKey:
    return (identity.chain.value, identity.residue_seq, identity.insertion_code)


@dataclass(frozen=True)
class SSEAssignmentMap:
    """Every residue of one chain that a ``HELIX``/``SHEET`` record covers.

    Residues absent from ``by_residue`` are coil.  The map is deliberately
    *sparse*: it holds what the file asserts, and coil is the absence of an
    assertion, so nothing here silently manufactures an annotation.
    """

    chain: str
    by_residue: Mapping[ResidueKey, SSEAssignment]
    #: Which convention named the individual β-strands in this map.  It is
    #: recorded on the map, travels into the result's SSE **annotation
    #: provenance**, and enters neither ``KeySpace`` nor ``SamplingProfile``:
    #: it changes no key and draws no triangle (ADR-063).
    strand_identity: str = STRAND_IDENTITY_SHEET_ID_AND_SERIAL
    helix_identity: str = HELIX_IDENTITY_BY_SERIAL
    version: str = SSE_ASSIGNMENT_VERSION
    n_helix_records: int = 0
    n_sheet_records: int = 0
    warnings: Tuple[str, ...] = ()

    def __len__(self) -> int:
        return len(self.by_residue)

    @property
    def has_records(self) -> bool:
        """Did the structure assert *any* secondary structure for this chain?

        The distinction this exposes is the one the legacy three-state collapse
        destroys: with no records at all, every residue is coil and every
        triangle is ``17_3c`` -- an answer that looks like a finding and is
        actually a missing input.
        """
        return bool(self.n_helix_records or self.n_sheet_records)

    def assignment_for(self, identity: AtomIdentity) -> SSEAssignment:
        """This residue's assignment; coil when no record covers it."""
        return self.by_residue.get(_residue_key(identity), UNASSIGNED)

    def coverage(self, identities: Sequence[AtomIdentity]) -> Tuple[int, int, int]:
        """``(n_helix, n_sheet, n_coil)`` over the residues actually in play."""
        n_helix = n_sheet = n_coil = 0
        for identity in identities:
            state = self.assignment_for(identity).state
            if state == HELIX:
                n_helix += 1
            elif state == SHEET:
                n_sheet += 1
            else:
                n_coil += 1
        return (n_helix, n_sheet, n_coil)

    def to_dict(self) -> Mapping[str, object]:
        return {
            "chain": self.chain,
            "strand_identity": self.strand_identity,
            "helix_identity": self.helix_identity,
            "version": self.version,
            "n_helix_records": self.n_helix_records,
            "n_sheet_records": self.n_sheet_records,
            "n_residues_annotated": len(self.by_residue),
            "has_records": self.has_records,
        }


def build_assignment_map(
    structure: Structure,
    chain: str,
    *,
    strand_identity: str = STRAND_IDENTITY_SHEET_ID_AND_SERIAL,
) -> SSEAssignmentMap:
    """Build one chain's residue -> assignment map from the structure's records.

    ``strand_identity`` names which fields identify an individual β-strand;
    :data:`STRAND_IDENTITY_SHEET_ID_AND_SERIAL` is the one convention.
    """
    if strand_identity not in STRAND_IDENTITIES:
        raise ValueError(
            "unknown strand identity convention %r; the defined convention is %r "
            "-- (sheetID, serNum), which names an individual beta-strand. "
            "serNum alone conflates strands that share a serial across sheets, "
            "and sheetID alone merges every strand of a sheet (ADR-063)."
            % (strand_identity, STRAND_IDENTITY_SHEET_ID_AND_SERIAL)
        )

    chain_id = ChainId(chain)
    helix: Dict[ResidueKey, SSEAssignment] = {}
    sheet: Dict[ResidueKey, SSEAssignment] = {}
    warnings: List[str] = []
    n_helix_records = n_sheet_records = 0

    records = [r for r in structure.secondary_structure if r.chain == chain_id]
    order, position, split = _chain_residue_order(structure, chain_id)
    for record in records:
        if record.kind == HELIX:
            n_helix_records += 1
            # HELIX_IDENTITY_BY_SERIAL: a helix serial is file-unique already.
            element = SSEElementId(serial=record.serial)
        else:
            n_sheet_records += 1
            # An individual beta-strand: the sheet it belongs to is part of its
            # identity, because serNum restarts per sheet.
            element = SSEElementId(serial=record.serial, sheet_id=record.group_id)
        first, last = _resolve_endpoints(structure, chain_id, record, position, split)
        target = helix if record.kind == HELIX else sheet
        assignment = SSEAssignment(state=record.kind, group=element)
        for key in order[first:last + 1]:
            # OVERLAP_FIRST_RECORD_WINS.
            if key not in target:
                target[key] = assignment

    # HELIX_BEFORE_SHEET: helix assignments take precedence over sheet ones.
    by_residue: Dict[ResidueKey, SSEAssignment] = dict(sheet)
    overlaps = 0
    for key, assignment in helix.items():
        if key in by_residue:
            overlaps += 1
        by_residue[key] = assignment
    if overlaps:
        warnings.append(
            "%d residue(s) of chain %r fall in both a HELIX and a SHEET record; "
            "HELIX wins, which is the historical precedence (the reader tests the "
            "helix table first)" % (overlaps, chain_id.value)
        )

    if not (n_helix_records or n_sheet_records):
        warnings.append(
            "no HELIX or SHEET record covers chain %r, so every residue is coil "
            "and every triangle will be category 17_3c. That is the historical "
            "behaviour and it is deterministic, but it is an absence of input "
            "rather than a finding -- check that the structure carries secondary-"
            "structure annotation before reading anything into the result."
            % chain_id.value
        )

    return SSEAssignmentMap(
        chain=chain_id.value,
        by_residue=by_residue,
        strand_identity=strand_identity,
        n_helix_records=n_helix_records,
        n_sheet_records=n_sheet_records,
        warnings=tuple(warnings),
    )


def _chain_residue_order(structure: Structure, chain_id: ChainId):
    """The chain's residues in file order, up to where its Calpha read stops.

    Returns ``(order, position, split)``: residue keys in order of first
    appearance, each key's index, and the keys whose records appear in two
    separate places (whose position is therefore not one place).  The read
    stops where the Calpha selection's does -- the chain's ``TER``, or the end
    of model 1 -- so a water or ligand written after the chain is not part of
    its order.
    """
    cutoff = None
    for event in structure.truncation_events:
        stops = (event.kind == "ENDMDL"
                 or (event.kind == "MODEL" and event.model_number is not None
                     and event.model_number > 1)
                 or (event.kind == "TER" and event.chain == chain_id))
        if stops and (cutoff is None or event.file_order < cutoff):
            cutoff = event.file_order
    order: List[ResidueKey] = []
    position: Dict[ResidueKey, int] = {}
    split = set()
    previous = None
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.chain != chain_id:
            continue
        key = _residue_key(atom.identity)
        if key == previous:
            continue
        if key in position:
            split.add(key)
        else:
            position[key] = len(order)
            order.append(key)
        previous = key
    return order, position, split


def _resolve_endpoints(structure, chain_id, record, position, split):
    """``(first, last)`` indices of a record's endpoints in the chain's order,
    or a refusal naming the record and why it cannot be read (ADR-110)."""
    what = "%s record %r (line %d: %s to %s, chain %s)" % (
        record.kind, record.serial, record.file_order + 1,
        record.endpoint_text("start"), record.endpoint_text("end"),
        chain_id.value or "(blank)")

    def refuse(reason):
        raise InvalidRequest(
            "SSE annotation of chain %s of %s is refused: %s %s. BioTSR matches "
            "HELIX/SHEET endpoints by chain, residue number and insertion code "
            "and covers the residues between them in the chain's order; it "
            "does not fall back to a numeric range or call the residues coil. "
            "Correct the record, or run without SSE."
            % (chain_id.value or "(blank)", structure.structure_id.value, what,
               reason))

    if not record.readable:
        refuse("has a residue number BioTSR cannot read (hybrid-36 numbers are "
               "not decoded)")
    if not record.spans_one_chain:
        refuse("ends on another chain (%s)" % (record.end_chain.value or "(blank)"))
    out = []
    for which, seq, icode in (("start", record.start_seq, record.start_icode),
                              ("end", record.end_seq, record.end_icode)):
        key = (chain_id.value, seq, icode)
        if key not in position:
            refuse("names %s residue %s%s, which is not in the chain"
                   % (which, seq, icode))
        if key in split:
            refuse("names %s residue %s%s, whose records appear in two separate "
                   "places in the chain, so its position is not defined"
                   % (which, seq, icode))
        out.append(position[key])
    if out[1] < out[0]:
        refuse("ends before it starts in the chain's order")
    return out[0], out[1]


@dataclass(frozen=True)
class SSEAnnotator:
    """Maps a triangle's three vertices to their assignments and a category.

    This is the whole of what the engine calls.  It reads ``AtomIdentity`` and
    nothing else -- no coordinates, no labels, no key -- which is precisely why
    annotating a run cannot move a key, and why the annotator is indifferent to
    grouping (which changes labels), to mirror (which changes a sign off
    ``theta_raw``) and to size filtering (which removes candidates before this is
    ever reached).
    """

    assignments: SSEAssignmentMap

    def annotate(
        self, identities: Sequence[AtomIdentity]
    ) -> Tuple[Tuple[SSEAssignment, SSEAssignment, SSEAssignment], str]:
        """``((a0, a1, a2), category)`` for one canonically ordered triangle."""
        triple = (
            self.assignments.assignment_for(identities[0]),
            self.assignments.assignment_for(identities[1]),
            self.assignments.assignment_for(identities[2]),
        )
        return triple, classify(triple)
