"""The PDB reader.

The one structure reader in this release.  It keeps the name ``pdb_legacy``
because that name is recorded in the provenance of every result; since
1.0.0rc6 it reads the atom-name and element fields at the PDB specification's
columns (ADR-108).

Field offsets used here (0-based Python slices):

===============  ===============  ==========================================
field            slice            note
===============  ===============  ==========================================
record name      ``[0:6]``
atom name        ``[12:16]``      PDB columns 13-16, the full field
                                  (``identity.atom_name``)
Calpha token     ``[13:16]``      columns 14-16, read by the Calpha test only
                                  (``calpha_name``); see below
altLoc           ``[16]``         case kept: ``a`` is not ``A``
residue name     ``[17:20]``
chain id         ``[21]``
residue seq      ``[22:26]``      an unreadable number (hybrid-36) keeps the
                                  record aside, see ``unreadable_numbers``
insertion code   ``[26]``         kept separate, never concatenated
x, y, z          ``[30:38]`` ...
element          ``[76:78]``      PDB columns 77-78 only, never the charge
                                  columns 79-80
===============  ===============  ==========================================

``HELIX`` / ``SHEET`` records are parsed too.  Their slices are
**identical in the legacy source and in the PDB v3.3 spec** -- unlike the atom
record, ``SSE_TSR.py`` reads every one of them at its correct offset:

===============  ==================  ==============================
field            HELIX               SHEET
===============  ==================  ==============================
serNum           ``[7:10]``          ``[7:10]`` (strand no. in sheet)
helixID/sheetID  ``[11:14]``         ``[11:14]``
initChainID      ``[19]``            ``[21]``
initSeqNum       ``[21:25]``         ``[22:26]``
initICode        ``[25]``            ``[26]``
endChainID       ``[31]``            ``[32]``
endSeqNum        ``[33:37]``         ``[33:37]``
endICode         ``[37]``            ``[37]``
===============  ==================  ==============================

The parser records them; it does not decide what they mean.  Which field carries
a strand's group identity, and what an unannotated residue is, are scientific
conventions and live in ``biotsr.sse``.

**On the atom-name field.**  The historical programs read the atom name from
columns 14-16 (``[13:16]``), which turns a name starting in column 13 -- a
two-letter element such as ``ZN`` or ``CL1`` -- into a different name (``N``,
``L1``).  Every atomic method now reads the full field, columns 13-16
(ADR-108).  The Calpha test does **not**: it keeps ``[13:16] == "CA"``, so a
calcium written ``CA`` from column 13 is never a Calpha, exactly as before.
``TSR-Package/tsr_package/TSR.py`` tests the two-character ``line[13:15]``;
for every atom name whose 16th column is blank the two agree, and this backend
records a warning whenever they disagree on an ``ATOM`` line, so a divergence
would be loud rather than silent.

**On residue sequence numbers.**  Legacy writes ``line[22:27].strip()`` -- the
number and the insertion code jammed together -- into its triplet output.  Here
they are separate fields, and ``ResidueIdentity.legacy_seq_label`` reconstructs
the legacy spelling for the legacy serializer only.  If a file ever contained a
sequence number whose stripped text does not round-trip through ``int`` (leading
zeros, for instance), that reconstruction would differ, so the parser checks and
warns rather than quietly producing different bytes.
"""

from __future__ import annotations

import hashlib
import os
from typing import List, Optional

from ...domain.identity import AtomIdentity, ChainId, StructureId
from ...domain.point import Point
from ..model import (
    AtomRecord,
    SecondaryStructureRecord,
    Structure,
    TruncationEvent,
    UnreadableNumberRecord,
)

__all__ = [
    "BACKEND_NAME",
    "COORD_PRECISION_DOUBLE",
    "backend_name_for",
    "parse_pdb_legacy",
    "read_pdb_legacy",
]

BACKEND_NAME = "pdb_legacy"

#: The coordinates as written in the file, read into a Python float. The standard.
COORD_PRECISION_DOUBLE = "float64"

_COORD_PRECISIONS = (COORD_PRECISION_DOUBLE,)


def backend_name_for(coordinate_precision: str) -> str:
    """The provenance name of this backend at a given coordinate precision."""
    if coordinate_precision == COORD_PRECISION_DOUBLE:
        return BACKEND_NAME
    return "%s+%s_coords" % (BACKEND_NAME, coordinate_precision)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _field(line: str, start: int, end: int) -> str:
    """A slice that tolerates short lines, as ``str`` slicing already does."""
    return line[start:end]


def _char(line: str, index: int) -> str:
    return line[index] if len(line) > index else " "


#: Column layout of the two secondary-structure records, as both the PDB v3.3
#: spec and ``SSE_TSR.py`` read them. ``(serial, group_id, init_chain,
#: init_seq, init_icode, end_chain, end_seq, end_icode)``, all 0-based.
_SSE_LAYOUT = {
    "HELIX": ((7, 10), (11, 14), 19, (21, 25), 25, 31, (33, 37), 37),
    "SHEET": ((7, 10), (11, 14), 21, (22, 26), 26, 32, (33, 37), 37),
}


def _int_or_none(text):
    try:
        return int(text)
    except ValueError:
        return None


def _parse_sse_record(line, order, warnings):
    """Parse one ``HELIX``/``SHEET`` line.

    A record whose endpoint number cannot be read (a hybrid-36 number, for
    example) is kept with that endpoint ``None`` and its text, and reported:
    the chain it names then refuses SSE annotation (ADR-110), so the record is
    never silently dropped and its residues never silently become coil.
    """
    kind = line[0:6].strip()
    layout = _SSE_LAYOUT[kind]
    (ser, gid, init_chain_col, init_seq_cols, init_icode_col,
     end_chain_col, end_seq_cols, end_icode_col) = layout

    serial = _field(line, *ser).strip()
    group_id = _field(line, *gid).strip()
    start_text = _field(line, *init_seq_cols).strip()
    end_text = _field(line, *end_seq_cols).strip()
    start_seq = _int_or_none(start_text)
    end_seq = _int_or_none(end_text)
    chain = ChainId(_char(line, init_chain_col).strip())
    if start_seq is None or end_seq is None:
        warnings.append(
            "line %d: %s record %r has a residue number BioTSR cannot read "
            "(%r..%r; hybrid-36 numbers are not decoded); SSE annotation of "
            "chain %s is refused" % (order + 1, kind, serial, start_text, end_text,
                                     chain.value or "(blank)")
        )
    return SecondaryStructureRecord(
        kind=kind,
        serial=serial,
        group_id=group_id,
        chain=chain,
        end_chain=ChainId(_char(line, end_chain_col).strip()),
        start_seq=start_seq,
        end_seq=end_seq,
        start_icode=_char(line, init_icode_col).strip(),
        end_icode=_char(line, end_icode_col).strip(),
        file_order=order,
        start_text=start_text if start_seq is None else "",
        end_text=end_text if end_seq is None else "",
    )


def parse_pdb_legacy(
    lines,
    structure_id: StructureId,
    *,
    source_path: Optional[str] = None,
    source_sha256: Optional[str] = None,
    coordinate_precision: str = COORD_PRECISION_DOUBLE,
) -> Structure:
    """Parse PDB text into a :class:`Structure`.

    Parsing is total and chain-agnostic: every ``ATOM``/``HETATM`` record in the
    file becomes an ``AtomRecord`` (or, when its residue number cannot be read,
    an ``UnreadableNumberRecord``), every ``MODEL``/``ENDMDL``/``TER`` record
    becomes a ``TruncationEvent``, and every ``HELIX``/``SHEET`` record becomes a
    ``SecondaryStructureRecord``.  No filtering, no truncation and no altLoc
    policy is applied here -- those are selection decisions.
    """
    if coordinate_precision not in _COORD_PRECISIONS:
        raise ValueError(
            "unknown coordinate_precision %r; BioTSR reads coordinates at %s"
            % (coordinate_precision, ", ".join(_COORD_PRECISIONS))
        )

    atoms: List[AtomRecord] = []
    events: List[TruncationEvent] = []
    sse_records: List[SecondaryStructureRecord] = []
    unreadable: List[UnreadableNumberRecord] = []
    warnings: List[str] = []
    current_model = 1

    for order, raw in enumerate(lines):
        line = raw.rstrip("\n").rstrip("\r")
        record = line[0:6]

        if line.startswith("MODEL"):
            token = _field(line, 10, 14).strip()
            try:
                current_model = int(token)
            except ValueError:
                warnings.append(
                    "line %d: MODEL record with unparseable serial %r; model number "
                    "left at %d" % (order + 1, token, current_model)
                )
                events.append(
                    TruncationEvent(kind="MODEL", file_order=order, model_number=None)
                )
                continue
            events.append(
                TruncationEvent(
                    kind="MODEL", file_order=order, model_number=current_model
                )
            )
            continue

        if line.startswith("ENDMDL"):
            events.append(TruncationEvent(kind="ENDMDL", file_order=order))
            continue

        if line.startswith("TER"):
            if len(line) > 21:
                chain_text = _char(line, 21).strip()
            else:
                chain_text = ""
                warnings.append(
                    "line %d: TER record shorter than the chain column; treated as "
                    "a blank chain id" % (order + 1)
                )
            events.append(
                TruncationEvent(
                    kind="TER", file_order=order, chain=ChainId(chain_text)
                )
            )
            continue

        if line.startswith("HELIX") or line.startswith("SHEET"):
            sse_records.append(_parse_sse_record(line, order, warnings))
            continue

        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        atom_name = _field(line, 12, 16).strip()
        calpha_name = _field(line, 13, 16).strip()
        if not atom_name:
            warnings.append(
                "line %d: %s record with a blank atom name; skipped"
                % (order + 1, record.strip())
            )
            continue

        # Prove, rather than assume, that the legacy two-character CA predicate and
        # the three-character atom-name field agree on this file.
        if line.startswith("ATOM"):
            if (_field(line, 13, 15).strip() == "CA") != (calpha_name == "CA"):
                warnings.append(
                    "line %d: legacy atom-name slices disagree -- [13:15]=%r vs "
                    "[13:16]=%r; Ca selection may diverge from legacy on this file"
                    % (order + 1, _field(line, 13, 15), _field(line, 13, 16))
                )

        seq_text = _field(line, 22, 26).strip()
        try:
            residue_seq = int(seq_text)
        except ValueError:
            # Kept aside, never silently dropped: a calculation whose selection
            # would include this record is refused (ADR-108).
            warnings.append(
                "line %d: residue number %r cannot be read (hybrid-36 numbers are "
                "not decoded); a calculation that would include this record is "
                "refused" % (order + 1, seq_text)
            )
            unreadable.append(UnreadableNumberRecord(
                file_order=order,
                hetero=line.startswith("HETATM"),
                chain=ChainId(_char(line, 21).strip()),
                residue_name=_field(line, 17, 20).strip(),
                atom_name=atom_name,
                alt_loc=_char(line, 16).strip(),
                seq_text=seq_text,
                element=_field(line, 76, 78).strip(),
            ))
            continue
        if str(residue_seq) != seq_text:
            warnings.append(
                "line %d: residue sequence text %r does not round-trip through int "
                "(%d); legacy output reconstructed from the integer would differ"
                % (order + 1, seq_text, residue_seq)
            )

        insertion_code = _char(line, 26).strip()
        alt_loc = _char(line, 16).strip()
        residue_name = _field(line, 17, 20).strip()
        chain_text = _char(line, 21).strip()

        try:
            x = float(_field(line, 30, 38))
            y = float(_field(line, 38, 46))
            z = float(_field(line, 46, 54))
        except ValueError:
            warnings.append(
                "line %d: unparseable coordinates; record skipped" % (order + 1)
            )
            continue

        identity = AtomIdentity(
            structure=structure_id,
            model=current_model,
            chain=ChainId(chain_text),
            residue_name=residue_name,
            residue_seq=residue_seq,
            insertion_code=insertion_code,
            atom_name=atom_name,
            alt_loc=alt_loc,
        )
        atoms.append(
            AtomRecord(
                identity=identity,
                coord=Point(x, y, z),
                hetero=line.startswith("HETATM"),
                element=_field(line, 76, 78).strip(),
                file_order=order,
                calpha_name=calpha_name,
                record_text=line[:6] + "     " + line[11:],
            )
        )

    return Structure(
        structure_id=structure_id,
        atoms=tuple(atoms),
        truncation_events=tuple(events),
        secondary_structure=tuple(sse_records),
        parser_backend=backend_name_for(coordinate_precision),
        unreadable_numbers=tuple(unreadable),
        source_path=source_path,
        source_sha256=source_sha256,
        warnings=tuple(warnings),
    )


def read_pdb_legacy(
    path: str,
    structure_id: Optional[StructureId] = None,
    *,
    coordinate_precision: str = COORD_PRECISION_DOUBLE,
) -> Structure:
    """Read a PDB file from disk.

    When ``structure_id`` is not supplied it is taken from the filename stem.
    **No four-character assumption is made**: the stem is used whole, whatever its
    length, so a ten-or-more-character identifier works unchanged.
    """
    if structure_id is None:
        stem = os.path.basename(path)
        if stem.lower().endswith(".pdb"):
            stem = stem[: -len(".pdb")]
        structure_id = StructureId.parse(stem)
    with open(path, "r") as fh:
        lines = fh.readlines()
    return parse_pdb_legacy(
        lines,
        structure_id,
        source_path=path,
        source_sha256=_sha256(path),
        coordinate_precision=coordinate_precision,
    )
