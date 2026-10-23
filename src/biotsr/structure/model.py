"""The format-neutral structure model.

Parser internals never travel beyond this boundary -- a direct response to the
legacy pattern where fixed-width column slices leak into the geometry loop.
Everything downstream sees ``AtomRecord``s, not lines.

Two things beyond the atoms are modelled because a faithful legacy reader needs
them and they are genuinely structural facts, not parser trivia:

* ``TruncationEvent`` -- ``MODEL``/``ENDMDL``/``TER`` records, with their file
  order.  The legacy protein Ca reader stops at the first ``ENDMDL``, the first
  ``TER`` *for the requested chain*, or the first ``MODEL`` numbered above 1.
  That rule depends on which chain was asked for, so it is a **selection**
  decision, not a parsing one; the parser records the facts and
  ``selection/truncation.py`` applies the rule.
* ``file_order`` on every record -- the line index it was read from.  It fixes a
  deterministic, source-faithful ordering for points and hence for candidate
  enumeration and serialization.
* ``SecondaryStructureRecord`` -- the ``HELIX`` and ``SHEET`` records.  They are
  structural facts about the deposited model, exactly like ``TER``/``MODEL``, and
  like those they are parsed here and *interpreted* elsewhere: ``biotsr.sse`` turns
  them into a residue-level annotation, and which fields identify an individual
  β-strand is a named convention there, not a parser choice (ADR-063).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from ..domain.identity import AtomIdentity, ChainId, StructureId
from ..domain.point import Point

__all__ = [
    "AtomRecord",
    "TruncationEvent",
    "SecondaryStructureRecord",
    "Structure",
    "UnreadableNumberRecord",
]


@dataclass(frozen=True)
class AtomRecord:
    """One parsed atom."""

    identity: AtomIdentity
    coord: Point
    hetero: bool
    element: str = ""
    occupancy: Optional[float] = None
    temperature_factor: Optional[float] = None
    file_order: int = -1

    #: ``line[13:16].strip()`` -- the atom-name columns 14-16 only. Read by the
    #: Calpha test alone (``calpha_name == "CA"``), which therefore keeps the
    #: reading every historical Calpha program used: a Calpha is written
    #: `` CA `` (columns 14-15), while calcium writes ``CA`` from column 13 and
    #: so reads ``A`` here and never qualifies. ``identity.atom_name`` is the
    #: full PDB atom-name field, columns 13-16, which every atomic method reads
    #: (ADR-108).
    calpha_name: str = ""
    #: The record line as read, with the atom serial number (columns 7-11)
    #: blanked: what tells a repeated record from a conflicting one (ADR-108).
    record_text: str = ""

    @property
    def chain(self) -> ChainId:
        return self.identity.chain

    @property
    def atom_name(self) -> str:
        return self.identity.atom_name

    @property
    def residue_name(self) -> str:
        return self.identity.residue_name


@dataclass(frozen=True)
class TruncationEvent:
    """A ``MODEL`` / ``ENDMDL`` / ``TER`` record, with enough context to apply a
    truncation rule without re-reading the file."""

    kind: str  # "MODEL" | "ENDMDL" | "TER"
    file_order: int
    chain: Optional[ChainId] = None
    model_number: Optional[int] = None


@dataclass(frozen=True)
class SecondaryStructureRecord:
    """One ``HELIX`` or ``SHEET`` record, parsed but not interpreted.

    Both PDB records declare a *residue range on one chain*.  Everything the
    family's only SSE implementation reads is here, and so is everything it
    ignores, because "ignored" is a characterized behaviour and not an absence:

    * ``serial`` -- the ``serNum`` token, columns 8-10.  On a ``HELIX`` record
      that is the helix's serial number, unique across the file; on a ``SHEET``
      record it is the **strand's index within its sheet**, which restarts at 1
      in every sheet and so does not identify a strand on its own.
      ``TSR-Package/tsr_package/SSE_TSR.py`` uses this field alone as the group
      identity for both records.
    * ``group_id`` -- the spec's own ``helixID`` / ``sheetID`` (columns 12-14).
      For a ``SHEET`` record this is half of an individual strand's standard
      identity: ``biotsr.sse`` pairs it with ``serial`` so strands numbered alike in
      different sheets stay apart (ADR-063).  For a ``HELIX`` record nothing
      reads it -- the serial already suffices, which was measured.
    * ``end_chain`` -- ``endChainID``.  The legacy reader never looks at it and
      assumes a range never spans chains; carrying it lets that assumption be
      checked rather than inherited.
    * ``start_icode`` / ``end_icode`` -- the range endpoints' insertion codes,
      kept separate from the sequence numbers for the same reason
      ``ResidueIdentity`` does.  The legacy reader ignored them; ``biotsr.sse``
      matches each endpoint by chain, number and insertion code (ADR-110).
    * ``start_text`` / ``end_text`` -- an endpoint number as written, kept only
      when it cannot be read as an integer (a hybrid-36 number, for example).
      Such a record is kept rather than dropped, with ``start_seq`` /
      ``end_seq`` ``None``, so the chain it names refuses SSE annotation
      instead of losing the record without a trace (ADR-110).
    """

    kind: str  # "HELIX" | "SHEET"
    serial: str
    group_id: str
    chain: ChainId
    end_chain: ChainId
    start_seq: Optional[int]
    end_seq: Optional[int]
    start_icode: str = ""
    end_icode: str = ""
    file_order: int = -1
    start_text: str = ""
    end_text: str = ""

    @property
    def spans_one_chain(self) -> bool:
        return self.chain == self.end_chain

    @property
    def readable(self) -> bool:
        """Whether both endpoint numbers could be read."""
        return self.start_seq is not None and self.end_seq is not None

    def endpoint_text(self, which: str) -> str:
        """``"start"`` / ``"end"`` endpoint as the file writes it, e.g. ``"14C"``."""
        seq, icode, text = ((self.start_seq, self.start_icode, self.start_text)
                            if which == "start" else
                            (self.end_seq, self.end_icode, self.end_text))
        return ("%d%s" % (seq, icode)) if seq is not None else "%s%s" % (text, icode)


@dataclass(frozen=True)
class UnreadableNumberRecord:
    """An ``ATOM``/``HETATM`` record whose residue number cannot be read.

    Columns 23-26 must hold an integer; a hybrid-36 number (``A000`` for
    10000, used by some tools past 9999 residues) does not, and BioTSR 1.0
    does not decode it.  The record is kept here, out of ``Structure.atoms``,
    so a selection that would include it is **refused** rather than computed
    without it and reported complete (ADR-108).  Nothing reads its
    coordinates.
    """

    file_order: int
    hetero: bool
    chain: ChainId
    residue_name: str
    atom_name: str
    alt_loc: str
    seq_text: str
    element: str = ""

    def describe(self) -> str:
        return "line %d: %s %s %r chain %s atom %s" % (
            self.file_order + 1, "HETATM" if self.hetero else "ATOM",
            self.residue_name, self.seq_text, self.chain.value or "(blank)",
            self.atom_name)


@dataclass(frozen=True)
class Structure:
    """A parsed structure: atoms in file order, plus structural markers."""

    structure_id: StructureId
    atoms: Tuple[AtomRecord, ...]
    truncation_events: Tuple[TruncationEvent, ...] = ()
    #: ``HELIX`` / ``SHEET`` records in file order. Empty is a legitimate state
    #: and means "this file carries no secondary-structure annotation" -- which
    #: is NOT the same fact as "every residue is coil", and ``biotsr.sse`` warns
    #: rather than letting the two be confused.
    secondary_structure: Tuple["SecondaryStructureRecord", ...] = ()
    parser_backend: str = ""
    #: Records kept out of ``atoms`` because their residue number cannot be
    #: read (:class:`UnreadableNumberRecord`); a selection that would include
    #: one is refused.
    unreadable_numbers: Tuple["UnreadableNumberRecord", ...] = ()
    source_path: Optional[str] = None
    source_sha256: Optional[str] = None
    warnings: Tuple[str, ...] = ()

    def __len__(self) -> int:
        return len(self.atoms)

    def chains(self) -> Tuple[ChainId, ...]:
        seen = []
        for atom in self.atoms:
            if atom.chain not in seen:
                seen.append(atom.chain)
        return tuple(seen)
