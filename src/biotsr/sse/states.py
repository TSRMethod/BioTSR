"""The three secondary-structure states, and one residue's assignment.

The family's only SSE implementation --
``TSR-Package/tsr_package/SSE_TSR.py`` -- works with exactly three states and
names them ``HELIX``, ``SHEET`` and ``NONE``.  Those names are kept, because they
are what the 18 category names and the historical file header are written in;
``COIL`` is the same state under the name the literature uses, and is provided as
an alias so the vocabulary is readable either way.

**There is no fourth "unknown" state, and this package does not invent one.**  The
legacy lookup is a two-dict membership test with an unconditional ``else``::

    if s in helix_dict:   ...HELIX
    elif s in sheet_dict: ...SHEET
    else:                 'NONE_0'

so "this residue is in no ``HELIX``/``SHEET`` record" and "this file has no
annotation at all" collapse into the same state.  That collapse is **evidenced,
total and deterministic**, so it is reproduced rather than replaced -- but it is
also lossy, and the standard path closes the loss without changing a single
category: ``biotsr.sse.assignment`` reports coverage and warns when a structure
carries no records at all, so "everything is coil" can never be mistaken for
"nothing was annotated".  See ADR-065.

What a state alone does not say is **which** helix or **which** strand, and that
is the distinction categories ``1_3a1`` and ``3_3a3`` are made of.  That identity
is :class:`SSEElementId` -- a structured pair, never a concatenated string, so
two elements are the same element only when every component agrees (ADR-063).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

__all__ = [
    "HELIX",
    "SHEET",
    "COIL",
    "NONE",
    "SSE_STATES",
    "SSEElementId",
    "UNASSIGNED_GROUP",
    "SSEAssignment",
    "UNASSIGNED",
]

#: A residue inside a ``HELIX`` record's range.
HELIX = "HELIX"

#: A residue inside a ``SHEET`` record's range.
SHEET = "SHEET"

#: A residue inside neither.  ``NONE`` is the legacy spelling and the one the
#: category names are written in; ``COIL`` is the same state.
NONE = "NONE"
COIL = NONE

#: The three states, in the order the category table is organised by.
SSE_STATES = (HELIX, SHEET, NONE)


@dataclass(frozen=True, order=True)
class SSEElementId:
    """Which *individual* secondary-structure element a residue belongs to.

    An individual β-strand is what the categories distinguish (ADR-063), and a
    ``SHEET`` record does not identify one with a single field: ``serNum`` is the
    strand's index **within its sheet** and restarts at 1 in every sheet, so
    ``sheetID`` and ``serNum`` together name a strand and neither does alone.
    Both are carried here, as separate components:

    * ``sheet_id``  -- the ``SHEET`` record's ``sheetID``.  Empty for a helix
      (``HELIX`` serials are unique across the whole file, so there is nothing to
      disambiguate) and empty for a strand under the historical
      ``legacy_serial_only`` convention, which is exactly what makes that
      convention collapse strands across sheets.
    * ``serial``    -- the record's ``serNum``.

    Equality and hashing are **componentwise**, which is the whole point: no
    concatenation happens anywhere on the identity path, so a sheet id that ends
    in a digit can never merge with the serial beside it and manufacture a
    collision.  :attr:`text` exists only to put the identity in a file or a
    message, and nothing compares two of those.
    """

    serial: str
    sheet_id: str = ""

    @cached_property
    def text(self) -> str:
        """A display spelling: ``"1"`` for a helix, ``"BET/1"`` for a strand.

        Presentation only.  Two elements are the same element because their
        components agree, never because this string does.

        Cached, because an element is created once per ``HELIX``/``SHEET`` record
        and then shared by every residue it covers and every triangle those
        residues appear in -- a structure's triangles number in the millions and
        its elements in the dozens.  Caching makes them all reference the same
        string instead of formatting a fresh one each time, which is where the
        legacy plain-string group's memory behaviour came from.  The cache cannot
        go stale: both components are frozen.
        """
        return "%s/%s" % (self.sheet_id, self.serial) if self.sheet_id else self.serial

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.text


#: The element identity a coil residue carries.  Legacy writes the literal token
#: ``NONE_0``; the ``0`` is that constant.
UNASSIGNED_GROUP = SSEElementId(serial="0")


@dataclass(frozen=True)
class SSEAssignment:
    """One residue's secondary-structure state and the element it belongs to.

    ``group`` is what distinguishes "two vertices in the *same* helix" from "two
    vertices in *different* helices" -- the distinction that separates category
    ``1_3a1`` from ``2_3a2``, and, for β-strands, ``7_3b1`` from ``9_3b3``.  It is
    an :class:`SSEElementId` compared only for equality; which PDB fields supply
    its components is a named convention
    (``assignment.STRAND_IDENTITY_SHEET_ID_AND_SERIAL``),
    never a hard-coded choice.
    """

    state: str
    group: SSEElementId = UNASSIGNED_GROUP

    def __post_init__(self) -> None:
        if self.state not in SSE_STATES:
            raise ValueError(
                "unknown secondary-structure state %r; it is one of %s"
                % (self.state, ", ".join(SSE_STATES))
            )
        if not isinstance(self.group, SSEElementId):
            raise TypeError(
                "an SSE group is an SSEElementId, not %r. The identity of an "
                "individual helix or beta-strand is a structured pair; passing a "
                "bare string is how sheetID gets concatenated onto serNum and two "
                "different strands start comparing equal."
                % (type(self.group).__name__,)
            )
        if self.state is NONE and self.group != UNASSIGNED_GROUP:
            raise ValueError(
                "a coil residue has no group; got %r. Coil is the absence of a "
                "HELIX/SHEET record, so there is nothing for two coil residues to "
                "be 'the same' or 'different' instances of -- which is why there "
                "is one 17_3c category and not three." % (self.group,)
            )

    @property
    def token(self) -> str:
        """The ``STATE_GROUP`` spelling, e.g. ``HELIX_1`` / ``NONE_0``.

        Under the historical convention this is the legacy token byte for byte
        (``SHEET_2``); under the standard one a strand spells its sheet out
        (``SHEET_BET/2``), because that is what it is.
        """
        return "%s_%s" % (self.state, self.group.text)

    @property
    def is_annotated(self) -> bool:
        return self.state is not NONE


#: The assignment every residue starts from and coil residues keep.
UNASSIGNED = SSEAssignment(state=NONE, group=UNASSIGNED_GROUP)
