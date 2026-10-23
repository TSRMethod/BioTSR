"""Wide dataset sheets: one row per structure, methods chosen on the command line.

A *wide* sheet describes one structure per row and carries the columns of many
calculations at once -- several participant chains (``aa_chain``,
``drug_chain``, ``carb_chain``, ``nt_chain``, ``aa_chain_2``) and several
alternative selections (``aa_ind``, ``aa_ind_mix``, ``aa_ind_specific``,
``drug_name`` ...).  Which calculations to run is **not** in the sheet: it is
chosen when BioTSR is invoked (``--method``, repeatable), and a populated
column never means "run another method" (ADR-097).

**Method-aware resolution.**  A *layout* says, per method and participant,
which column supplies each canonical argument::

    carbohydrate, participant 1:  chain <- carb_chain, entity_name <- carb_name,
                                  entity_id <- carb_id
    aa_ind, participant 1:        chain <- aa_chain,   residue_types <- aa_ind
    aa_ind, participant 2:        chain <- aa_chain_2, residue_types <- aa_ind_2
    carbohydrate_protein_cross:   chain <- carb_chain, chain_2 <- aa_chain,
                                  entity_name <- carb_name, entity_id <- carb_id

so each method reads its own participant columns and nothing else, and one
sheet serves several methods whose chain columns differ.  The built-in ``tsr``
layout is the TSR group's dataset-sheet format (``biotsr.io.input_columns``);
``--map METHOD[.PARTICIPANT].FIELD=HEADER`` adds or overrides an entry, and the
``custom`` layout starts empty, so arbitrary headers work.

**Cross methods.**  A cross method's participant 1 covers both of its sides.
The ``tsr`` layout reads their chains and molecule names and identifiers; a
side can also be narrowed by one of several columns (the protein side of
``drug_protein_cross`` by ``aa_ind``, ``aa_ind_mix``, ``aa_ind_specific`` or
``aa_atom_selection``; the drug side by ``drug_atom_selected``), and the sheet
does not say which one a cross calculation means.  So BioTSR neither picks one
nor silently ignores them: when a row populates any of a side's candidate
columns, the side must be chosen on the command line --
``--cross-side drug_protein_cross.protein=aa_ind`` narrows by that column,
``--cross-side drug_protein_cross.protein=whole`` uses the whole side -- and
the run is refused until it is.  The columns a choice leaves unused are named
in the preview and the summary; a row whose chosen column is blank is
inapplicable (it does not describe that narrowed calculation), never computed
on the whole side.  :data:`CROSS_SIDE_COLUMNS` lists every side's candidates.

Participant 2 -- the ``_2`` columns describing a second protein chain -- runs
only when asked for (``--participant 2`` or ``all``).

A column no planned calculation reads is reported as unused; a ``comment*``
column is ignored.  Blank means unset, never a default.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from ..methods.arguments import accepted_arguments
from ..methods.registry import get_method
from .errors import CLIError
from .manifest import MANIFEST_ENCODING, _rows

__all__ = [
    "LAYOUT_TSR",
    "LAYOUT_CUSTOM",
    "LAYOUTS",
    "PARTICIPANTS",
    "PARTICIPANT_CHOICES",
    "SHEET_HEADER_ALIASES",
    "SheetMap",
    "SheetRequest",
    "CROSS_SIDE_COLUMNS",
    "WHOLE_SIDE",
    "builtin_layout",
    "parse_cross_side_options",
    "parse_map_options",
    "read_map_file",
    "read_sheet",
]

LAYOUT_TSR = "tsr"
LAYOUT_CUSTOM = "custom"
LAYOUTS = (LAYOUT_TSR, LAYOUT_CUSTOM)

PARTICIPANTS = ("1", "2")
PARTICIPANT_ALL = "all"
PARTICIPANT_CHOICES = PARTICIPANTS + (PARTICIPANT_ALL,)

#: Header spellings that are the same column.  ``car_atom_selected`` is the
#: header the TSR group's sheets carry; their column dictionary names it
#: ``carb_atom_selected``.  Both are read as the same column; a sheet carrying
#: both is refused.
SHEET_HEADER_ALIASES = {"car_atom_selected": "carb_atom_selected"}

#: Fields read from the sheet for every request, whatever the method.
GLOBAL_FIELDS = ("pdb_id", "pdb_path", "structure_id", "sample_id", "label")

#: ``--cross-side METHOD.SIDE=whole``: the side is used whole, unnarrowed.
WHOLE_SIDE = "whole"

_PROTEIN_SIDE = {"aa_ind": "residue_types", "aa_ind_mix": "residue_types",
                 "aa_ind_specific": "residues", "aa_atom_selection": "atom_selection"}
_PROTEIN_SIDE_2 = {"aa_ind_2": "residue_types_2", "aa_ind_mix_2": "residue_types_2",
                   "aa_ind_specific_2": "residues_2",
                   "aa_atom_selection_2": "atom_selection_2"}

#: The ``tsr`` layout's candidate narrowing columns of each cross side, as
#: ``{method: {side: {column: field}}}``, from the column dictionary
#: (``biotsr.io.input_columns``).  ``aa_selection`` narrows only the Calpha
#: cross: the atomic cross methods label atoms, not residues (ADR-045).
CROSS_SIDE_COLUMNS: Mapping[str, Mapping[str, Mapping[str, str]]] = {
    "drug_protein_cross": {"drug": {"drug_atom_selected": "atoms"},
                           "protein": _PROTEIN_SIDE},
    "carbohydrate_protein_cross": {"carbohydrate": {"carb_atom_selected": "atoms"},
                                   "protein": _PROTEIN_SIDE},
    "nucleotide_protein_cross": {"nucleotide": {"nt_atom_selected": "atoms"},
                                 "protein": _PROTEIN_SIDE},
    "protein_protein_cross": {"1": _PROTEIN_SIDE, "2": _PROTEIN_SIDE_2},
    "protein_protein_calpha_cross": {"1": {"aa_selection": "calpha_residues"},
                                     "2": {"aa_selection_2": "calpha_residues_2"}},
}


def builtin_layout(name: str) -> Tuple[Dict[str, str], Dict[str, Dict[str, Dict[str, str]]]]:
    """``(global fields, {method: {participant: {field: column}}})`` of a layout."""
    if name == LAYOUT_CUSTOM:
        return {}, {}
    if name != LAYOUT_TSR:
        raise CLIError("unknown sheet layout %r" % name,
                       "layouts: %s" % ", ".join(LAYOUTS))
    aa = ({"chain": "aa_chain"}, {"chain": "aa_chain_2"})

    def two(field_name=None, first=None, second=None):
        one = dict(aa[0])
        other = dict(aa[1])
        if field_name:
            one[field_name] = first
            other[field_name] = second
        return {"1": one, "2": other}

    def entity(chain, name, ident, atoms=None):
        out = {"chain": chain, "entity_name": name, "entity_id": ident}
        if atoms:
            out["atoms"] = atoms
        return {"1": out}

    def cross(chain, name=None, ident=None):
        out = {"chain": chain, "chain_2": "aa_chain"}
        if name:
            out.update(entity_name=name, entity_id=ident)
        return {"1": out}

    methods = {
        "protein": two(),
        "aa_selection": two("residues", "aa_selection", "aa_selection_2"),
        "amino_acid_atomic": two(),
        "amino_acid_atom_selection": two("atoms", "aa_atom_selection",
                                         "aa_atom_selection_2"),
        "aa_ind": two("residue_types", "aa_ind", "aa_ind_2"),
        "aa_ind_mix": two("residue_types", "aa_ind_mix", "aa_ind_mix_2"),
        "aa_ind_specific": two("residues", "aa_ind_specific", "aa_ind_specific_2"),
        "drug": entity("drug_chain", "drug_name", "drug_id"),
        "drug_atom_selection": entity("drug_chain", "drug_name", "drug_id",
                                      "drug_atom_selected"),
        "carbohydrate": entity("carb_chain", "carb_name", "carb_id"),
        "carbohydrate_atom_selection": entity("carb_chain", "carb_name", "carb_id",
                                              "carb_atom_selected"),
        "carbohydrate_chain": {"1": {"chain": "carb_chain"}},
        "nucleotide": entity("nt_chain", "nt_name", "nt_id"),
        "nucleotide_atom_selection": entity("nt_chain", "nt_name", "nt_id",
                                            "nt_atom_selected"),
        "nucleotide_chain": {"1": {"chain": "nt_chain"}},
        "drug_protein_cross": cross("drug_chain", "drug_name", "drug_id"),
        "carbohydrate_protein_cross": cross("carb_chain", "carb_name", "carb_id"),
        "nucleotide_protein_cross": cross("nt_chain", "nt_name", "nt_id"),
        "protein_protein_cross": {"1": {"chain": "aa_chain", "chain_2": "aa_chain_2"}},
        "protein_protein_calpha_cross": {"1": {"chain": "aa_chain",
                                               "chain_2": "aa_chain_2"}},
    }
    return {"pdb_id": "protein", "label": "aa_group"}, methods


def _fields_of(method: str) -> Tuple[str, ...]:
    recipe = get_method(method)
    accepted = tuple(accepted_arguments(recipe))
    if recipe.family == "cross":
        return ("chain",) + accepted
    return ("chain",) + accepted


def parse_map_options(values: Sequence[str]) -> List[Tuple[str, str, str, str, str]]:
    """``["carbohydrate.chain=glycan_chain", "aa_ind.2.chain=chain_b"]`` -> entries.

    Each entry is ``(method, participant, field, column, spelling)``.
    """
    out = []
    for item in values or ():
        text = str(item)
        target, sep, column = text.partition("=")
        parts = [p.strip() for p in target.split(".")]
        column = column.strip()
        if not sep or not column or len(parts) not in (2, 3) or not all(parts):
            raise CLIError(
                "--map %r is not METHOD[.PARTICIPANT].FIELD=HEADER" % text,
                "e.g. --map carbohydrate.chain=glycan_chain, or --map "
                "aa_ind.2.residue_types=second_types for participant 2",
            )
        if len(parts) == 2:
            method, field_name = parts
            participant = "1"
        else:
            method, participant, field_name = parts
        out.append((method, participant, field_name, column, "--map %s" % text))
    return out


def parse_cross_side_options(values: Sequence[str]) -> Dict[Tuple[str, str], str]:
    """``["drug_protein_cross.protein=aa_ind"]`` -> ``{(method, side): column}``.

    The column may be ``whole``.  Methods, sides and columns are checked
    against :data:`CROSS_SIDE_COLUMNS`; each side is chosen once.
    """
    out: Dict[Tuple[str, str], str] = {}
    for item in values or ():
        text = str(item)
        target, sep, column = text.partition("=")
        method, dot, side = target.strip().partition(".")
        column = SHEET_HEADER_ALIASES.get(column.strip(), column.strip())
        if not sep or not dot or not method or not side or not column:
            raise CLIError(
                "--cross-side %r is not METHOD.SIDE=COLUMN" % text,
                "e.g. --cross-side drug_protein_cross.protein=aa_ind, or "
                "--cross-side drug_protein_cross.protein=whole",
            )
        sides = CROSS_SIDE_COLUMNS.get(method)
        if sides is None:
            raise CLIError("--cross-side %s: %r is not a cross method"
                           % (text, method),
                           "cross methods: %s" % ", ".join(sorted(CROSS_SIDE_COLUMNS)))
        if side not in sides:
            raise CLIError("--cross-side %s: %s has no side %r"
                           % (text, method, side),
                           "its sides: %s" % ", ".join(sides))
        if column != WHOLE_SIDE and column not in sides[side]:
            raise CLIError(
                "--cross-side %s: column %r does not narrow the %s side of %s"
                % (text, column, side, method),
                "choose one of: %s, or %s" % (", ".join(sides[side]), WHOLE_SIDE))
        if (method, side) in out:
            raise CLIError("--cross-side chooses the %s side of %s twice"
                           % (side, method), "choose each side once")
        out[(method, side)] = column
    return out


def read_map_file(path: str) -> List[Tuple[str, str, str, str, str]]:
    """A mapping file: CSV with columns ``method,participant,field,column``.

    ``participant`` may be blank (participant 1).  Lines are reported by number.
    """
    if not os.path.isfile(path):
        raise CLIError("--map-file not found: %s" % path)
    out = []
    with open(path, "r", newline="", encoding=MANIFEST_ENCODING) as fh:
        reader = csv.reader(fh)
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration:
            raise CLIError("--map-file %s is empty" % path)
        wanted = ["method", "participant", "field", "column"]
        if header != wanted:
            raise CLIError(
                "--map-file %s must have the header %s; it has %s"
                % (path, ",".join(wanted), ",".join(header)),
            )
        for line, row in _rows(reader, path):
            if not any(c.strip() for c in row):
                continue
            if len(row) != 4:
                raise CLIError("--map-file %s line %d has %d cells, not 4"
                               % (path, line, len(row)))
            method, participant, field_name, column = (c.strip() for c in row)
            out.append((method, participant or "1", field_name, column,
                        "%s line %d" % (os.path.basename(path), line)))
    return out


@dataclass
class SheetMap:
    """The resolved mapping of one sheet: global fields and per-method participants."""

    global_fields: Dict[str, str]
    methods: Dict[str, Dict[str, Dict[str, str]]]
    #: Where each user mapping came from, for messages and the preview.
    user_entries: List[Tuple[str, str, str, str, str]] = field(default_factory=list)


@dataclass
class SheetRequest:
    """One planned request from one sheet row: a method, a participant, its values."""

    line: int
    where: str
    method: str
    participant: str
    values: Dict[str, Any]
    #: ``{field: column}`` actually read for this request.
    columns: Dict[str, str]
    metadata: Dict[str, str]
    #: Set when the row does not describe this calculation at all.
    inapplicable: str = ""
    #: Set when the row describes it but cannot run it as given.
    problem: str = ""
    #: The chain fields (``chain``, ``chain_2``) the method needs and the row
    #: leaves blank, with the columns they are read from; the batch applies an
    #: explicit ``--chain`` default before deciding the chain is not specified.
    chain_missing: Dict[str, str] = field(default_factory=dict)
    #: Populated columns that could narrow a cross side and were not applied,
    #: each with why (``[("aa_ind_mix", "the protein side uses aa_ind")]``).
    not_applied: List[Tuple[str, str]] = field(default_factory=list)


def _resolve_headers(header: Sequence[str], path: str) -> Dict[str, int]:
    """``{column: index}`` with aliases folded, refusing duplicates."""
    index: Dict[str, int] = {}
    spelled: Dict[str, str] = {}
    for i, raw in enumerate(header):
        name = str(raw or "").strip()
        if not name:
            continue
        canonical = SHEET_HEADER_ALIASES.get(name, name)
        if canonical in index:
            raise CLIError(
                "sheet %s has two columns for %r: %r and %r"
                % (path, canonical, spelled[canonical], name),
                "they are spellings of one column (%s); keep one"
                % ", ".join(sorted([canonical] + [a for a, c in
                                                  SHEET_HEADER_ALIASES.items()
                                                  if c == canonical])),
            )
        index[canonical] = i
        spelled[canonical] = name
    return index


def build_map(layout: str, header_index: Mapping[str, int], *, columns: Mapping[str, str],
              user_entries: Sequence[Tuple[str, str, str, str, str]],
              methods: Sequence[str]) -> SheetMap:
    """Combine a layout, ``--column`` globals and ``--map`` entries; refuse anything unclear."""
    global_fields, layout_methods = builtin_layout(layout)
    present = set(header_index)
    listing = ", ".join(sorted(present)) or "(none)"

    def column_of(name: str) -> str:
        return SHEET_HEADER_ALIASES.get(name, name)

    for field_name, source in sorted(columns.items()):
        if field_name not in GLOBAL_FIELDS:
            raise CLIError(
                "--column %s=%s: in a sheet, --column maps the fields every request "
                "shares (%s); a method's own fields are mapped with --map "
                "METHOD[.PARTICIPANT].FIELD=HEADER" % (field_name, source,
                                                        ", ".join(GLOBAL_FIELDS)),
            )
        if column_of(source) not in present:
            raise CLIError("--column %s=%s: the sheet has no column %r"
                           % (field_name, source, source), "its columns are: %s" % listing)
        global_fields[field_name] = column_of(source)
    if layout == LAYOUT_CUSTOM:
        for field_name in ("pdb_id", "pdb_path"):
            if field_name not in global_fields and field_name in present:
                global_fields[field_name] = field_name
    for field_name, source in list(global_fields.items()):
        if source not in present:
            del global_fields[field_name]

    merged = {m: {p: dict(f) for p, f in parts.items()}
              for m, parts in layout_methods.items()}
    seen = {}
    from ..api import list_methods

    known = set(list_methods())
    for method, participant, field_name, column, spelling in user_entries:
        if method not in known:
            raise CLIError("%s: unknown method %r" % (spelling, method),
                           "methods: %s" % ", ".join(sorted(known)))
        if participant not in PARTICIPANTS:
            raise CLIError("%s: participant %r is not 1 or 2" % (spelling, participant))
        recipe = get_method(method)
        if recipe.family == "cross" and participant != "1":
            raise CLIError(
                "%s: a cross method has one participant, covering both sides; "
                "map its second side with the *_2 fields" % spelling)
        allowed = _fields_of(method)
        if field_name not in allowed:
            raise CLIError(
                "%s: method %s does not read %r" % (spelling, method, field_name),
                "its fields: %s" % ", ".join(allowed))
        if column_of(column) not in present:
            raise CLIError("%s: the sheet has no column %r" % (spelling, column),
                           "its columns are: %s" % listing)
        key = (method, participant, field_name)
        if key in seen:
            raise CLIError(
                "%s and %s both map %s participant %s %s"
                % (seen[key], spelling, method, participant, field_name),
                "map each field once")
        seen[key] = spelling
        merged.setdefault(method, {}).setdefault(participant, {})[field_name] = \
            column_of(column)

    for method in methods:
        if not merged.get(method):
            raise CLIError(
                "the %s layout has no columns for method %s" % (layout, method),
                "map them with --map %s.chain=HEADER and the method's other "
                "fields (%s)" % (method, ", ".join(_fields_of(method)[1:]) or "none"),
            )
    return SheetMap(global_fields=global_fields, methods=merged,
                    user_entries=list(user_entries))


def _defining(method: str, fields: Mapping[str, str]) -> Tuple[str, ...]:
    """The fields whose presence says a row describes this calculation at all.

    A chain column is shared by many calculations, so it does not say a row
    wants this one -- unless the method reads nothing else.
    """
    recipe = get_method(method)
    if recipe.family == "cross":
        if "entity_name" in fields:
            return ("entity_name",)
        return ("chain_2",)
    if recipe.selector == "entity_atom_selection" and "atoms" in fields:
        # The atom subset is what makes a row an atom-selection calculation;
        # a row naming the molecule without one describes the whole molecule.
        return ("atoms",)
    specific = tuple(f for f in fields if f not in ("chain",))
    return specific or ("chain",)


def read_sheet(path: str, *, layout: str, methods: Sequence[str], participant: str,
               columns: Mapping[str, str], maps: Sequence[Tuple[str, str, str, str, str]],
               keep_columns: Sequence[str] = (), ignore_columns: Sequence[str] = (),
               cross_sides: Optional[Mapping[Tuple[str, str], str]] = None,
               ) -> Tuple[List[SheetRequest], SheetMap, List[str]]:
    """Plan every (row, method, participant) of a wide sheet.

    Returns ``(requests, sheet_map, unused_columns)``.  Rows that do not
    describe a calculation are returned with ``inapplicable`` set, never
    dropped; rows that describe one incompletely carry ``problem``.
    """
    if not os.path.isfile(path):
        raise CLIError("manifest not found: %s" % path)
    if not methods:
        raise CLIError(
            "a sheet names no methods: choose them with --method (repeatable)",
            "e.g. --method carbohydrate --method aa_ind; `biotsr methods` lists them",
        )
    if participant not in PARTICIPANT_CHOICES:
        raise CLIError("--participant %r is not one of %s"
                       % (participant, ", ".join(PARTICIPANT_CHOICES)))
    with open(path, "r", newline="", encoding=MANIFEST_ENCODING) as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise CLIError("sheet %s is empty" % path)
        index = _resolve_headers(header, path)
        #: The header as the sheet spells it, for messages.
        spelled = {SHEET_HEADER_ALIASES.get(str(h).strip(), str(h).strip()): str(h).strip()
                   for h in header if str(h or "").strip()}
        if "method" in index and "method" not in ignore_columns:
            raise CLIError(
                "sheet %s has a 'method' column, and --sheet takes its methods from "
                "--method" % path,
                "pass --ignore-column method to run the --method list on every "
                "row, or run the file as a one-method-per-row manifest (without "
                "--sheet and without --method)",
            )
        for option, names in (("--keep-column", keep_columns),
                              ("--ignore-column", ignore_columns)):
            for name in names:
                if SHEET_HEADER_ALIASES.get(name, name) not in index:
                    raise CLIError("%s %s: sheet %s has no column %r"
                                   % (option, name, path, name))
        sheet_map = build_map(layout, index, columns=columns, user_entries=maps,
                              methods=methods)
        if "pdb_id" not in sheet_map.global_fields and \
                "pdb_path" not in sheet_map.global_fields:
            raise CLIError(
                "sheet %s has no structure column" % path,
                "name it with --column pdb_id=HEADER (or pdb_path=HEADER)",
            )
        wanted_participants = (PARTICIPANTS if participant == PARTICIPANT_ALL
                               else (participant,))
        for method in methods:
            available = [p for p in wanted_participants
                         if sheet_map.methods.get(method, {}).get(p)]
            if not available:
                raise CLIError(
                    "method %s has no participant %s in the %s layout"
                    % (method, participant, layout),
                    "it has: %s" % ", ".join(sorted(sheet_map.methods.get(method, {})))
                    + ("; a cross method's participant 1 covers both of its sides"
                       if get_method(method).family == "cross" else ""),
                )
        cross_sides = dict(cross_sides or {})
        for (method, side), column in sorted(cross_sides.items()):
            if method not in methods:
                raise CLIError(
                    "--cross-side %s.%s=%s: %s is not among the --method choices"
                    % (method, side, column, method))
            if layout != LAYOUT_TSR:
                raise CLIError(
                    "--cross-side chooses among the tsr layout's columns; with "
                    "--sheet %s map the side's field with --map" % layout)
            if column != WHOLE_SIDE and column not in index:
                raise CLIError("--cross-side %s.%s=%s: sheet %s has no column %r"
                               % (method, side, column, path, column))
            mapped = [f for (m, _p, f, _c, _s) in sheet_map.user_entries
                      if m == method and f in CROSS_SIDE_COLUMNS[method][side].values()]
            if mapped:
                raise CLIError(
                    "the %s side of %s is chosen twice: by --cross-side and by "
                    "--map (%s)" % (side, method, ", ".join(mapped)),
                    "choose it once")
        #: ``{(method, side): {column: [lines]}}`` of unchosen populated columns.
        undecided: Dict[Tuple[str, str], Dict[str, List[int]]] = {}
        requests: List[SheetRequest] = []
        read_columns = set(sheet_map.global_fields.values())
        for line, raw in _rows(reader, path):
            if not any(str(cell).strip() for cell in raw):
                continue
            where = "%s line %d" % (os.path.basename(path), line)
            if len(raw) > len(header) and any(str(c).strip() for c in raw[len(header):]):
                raise CLIError(
                    "%s: the row has %d cells and the header %d" % (where, len(raw),
                                                                     len(header)),
                    "a value without a header would be silently dropped; separate "
                    "list items with ';' or quote the cell",
                )

            def cell(column: str) -> str:
                i = index.get(column)
                return str(raw[i]).strip() if i is not None and i < len(raw) else ""

            shared = {}
            for field_name, column in sheet_map.global_fields.items():
                if cell(column):
                    shared[field_name] = cell(column)
            metadata = {name: cell(SHEET_HEADER_ALIASES.get(name, name))
                        for name in keep_columns}
            for method in methods:
                for part in wanted_participants:
                    fields = sheet_map.methods.get(method, {}).get(part)
                    if not fields:
                        continue
                    values = dict(shared)
                    used = {}
                    for field_name, column in fields.items():
                        read_columns.add(column)
                        if cell(column):
                            values[field_name] = cell(column)
                            used[field_name] = column
                    request = SheetRequest(
                        line=line, where=where, method=method, participant=part,
                        values=values, columns=used, metadata=metadata,
                    )
                    defining = _defining(method, fields)
                    describes = any(values.get(f) for f in defining)
                    if describes and layout == LAYOUT_TSR \
                            and method in CROSS_SIDE_COLUMNS:
                        _choose_cross_sides(request, sheet_map, cross_sides, cell,
                                            spelled, undecided, read_columns)
                    if not describes:
                        request.inapplicable = (
                            "the row has no %s for %s%s"
                            % (" / ".join(spelled.get(fields[f], fields[f])
                                          for f in defining), method,
                               " (participant 2)" if part == "2" else ""))
                    else:
                        request.chain_missing = {
                            f: spelled.get(fields[f], fields[f])
                            for f in ("chain", "chain_2")
                            if f in fields and not values.get(f)}
                    requests.append(request)
    if undecided:
        lines = []
        for (method, side), found in sorted(undecided.items()):
            lines.append("%s of %s: %s" % (_side_name(side), method, "; ".join(
                "%s on line %s" % (spelled.get(c, c), _numbers(n))
                for c, n in sorted(found.items()))))
        first_method, first_side = sorted(undecided)[0]
        choices = ", ".join(CROSS_SIDE_COLUMNS[first_method][first_side])
        raise CLIError(
            "the sheet populates columns that could narrow a cross calculation, "
            "and BioTSR does not choose among them: %s" % " | ".join(lines),
            "choose each side on the command line: --cross-side %s.%s=COLUMN "
            "narrows it by that column (%s), and --cross-side %s.%s=whole uses "
            "the whole side; `biotsr batch --help` lists every side"
            % (first_method, first_side, choices, first_method, first_side),
        )
    ignored = set(SHEET_HEADER_ALIASES.get(n, n) for n in ignore_columns)
    kept = set(SHEET_HEADER_ALIASES.get(n, n) for n in keep_columns)
    unused = sorted(spelled.get(c, c) for c in index
                    if c not in read_columns and c not in ignored
                    and c not in kept and not c.startswith("comment"))
    if not requests:
        raise CLIError("sheet %s has a header but no data rows" % path)
    return requests, sheet_map, unused


def _numbers(lines: Sequence[int], limit: int = 8) -> str:
    shown = [str(n) for n in lines[:limit]]
    more = len(lines) - len(shown)
    return ",".join(shown) + (" and %d more" % more if more > 0 else "")


def _choose_cross_sides(request: SheetRequest, sheet_map: SheetMap,
                        cross_sides: Mapping[Tuple[str, str], str], cell,
                        spelled: Mapping[str, str],
                        undecided: Dict[Tuple[str, str], Dict[str, List[int]]],
                        read_columns: set) -> None:
    """Apply the ``--cross-side`` choices to one cross request (in place)."""
    method = request.method
    mapped_fields = {f for (m, _p, f, _c, _s) in sheet_map.user_entries if m == method}
    for side, candidates in CROSS_SIDE_COLUMNS[method].items():
        named = _side_name(side)
        present = [c for c in candidates if c in spelled]
        read_columns.update(present)
        populated = [c for c in present if cell(c)]
        choice = cross_sides.get((method, side))
        side_mapped = [f for f in candidates.values() if f in mapped_fields]
        if choice is None and side_mapped:
            # Chosen with --map: the mapped field was read like any other;
            # the side's other populated columns are named, not applied.
            used = set(request.columns.values())
            for column in populated:
                if column not in used:
                    request.not_applied.append((spelled.get(column, column),
                                                "%s is mapped with --map" % named))
            continue
        if choice is None:
            for column in populated:
                undecided.setdefault((method, side), {}).setdefault(
                    column, []).append(request.line)
            continue
        if choice == WHOLE_SIDE:
            for column in populated:
                request.not_applied.append((spelled.get(column, column),
                                            "%s is used whole" % named))
            continue
        field_name = candidates[choice]
        if cell(choice):
            request.values[field_name] = cell(choice)
            request.columns[field_name] = choice
        elif not request.inapplicable:
            request.inapplicable = (
                "the row leaves %s blank, chosen for %s of %s "
                "(--cross-side)" % (spelled.get(choice, choice), named, method))
        for column in populated:
            if column != choice:
                request.not_applied.append((spelled.get(column, column),
                                            "%s uses %s" % (named, choice)))


def _side_name(side: str) -> str:
    return "side %s" % side if side.isdigit() else "the %s side" % side
