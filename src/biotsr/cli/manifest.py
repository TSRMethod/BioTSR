"""The canonical dataset manifest: one row, one run.

A *wide* dataset sheet -- the TSR group's historical format, where one row is
one PDB entry carrying up to a dozen calculations, each with its own chain column (``aa_chain``,
``drug_chain``, ``carb_chain``, ``nt_chain``) -- cannot be a manifest: a row does
not describe a run, and ``chain`` is ambiguous *within* a row.
Such a sheet is read method-aware by ``biotsr.cli.sheet`` (``biotsr batch
--sheet``), where each method reads its own columns.

So the canonical manifest is long/tidy: **one row is exactly one ``biotsr run``**,
with the same argument names the command line uses, minus the leading dashes.
Anybody who can read ``biotsr run --help`` can write a manifest, and anybody reading
a manifest knows what command it stands for.

Aliases accept the column names used in the historical dataset sheets, because
reshaping a wide sheet into a long one should not also mean renaming every
column.  Any other column name can be mapped explicitly (``--column
label=group``), carried as metadata (``--keep-column``) or dropped
(``--ignore-column``).  Nothing is guessed: an unknown header is an error naming
the remedies, and two different aliases of the *same* canonical column populated
in one row is an error naming both rather than a silent precedence rule.

Blank means unset.  It never means a default -- defaults belong to the package,
or to an explicit command-line row default such as ``--method`` or ``--chain``.
"""

from __future__ import annotations

import csv
import difflib
import os
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from .errors import CLIError
from .spec import CROSS_ARGUMENTS, MODIFIER_ARGUMENTS, SELECTION_ARGUMENTS

#: How manifests are opened: UTF-8, with a byte-order mark accepted and removed.
#: Spreadsheet programs (Excel's "CSV UTF-8") write one, and without this the
#: first header would carry an invisible character and be "unrecognised".
MANIFEST_ENCODING = "utf-8-sig"

#: Columns earlier BioTSR releases read and this one does not, with the reason.
#: An unknown header is an error in any case; these get a message that says why.
_REMOVED_COLUMNS = {
    "profile": (
        "BioTSR runs every method under its standard profile, so there is no "
        "profile to choose; drop the column (or pass --ignore-column profile)"
    ),
}

__all__ = [
    "REQUIRED_COLUMNS",
    "UNIVERSAL_COLUMNS",
    "CANONICAL_COLUMNS",
    "COLUMN_ALIASES",
    "BOOLEAN_COLUMNS",
    "FLOAT_COLUMNS",
    "CHAIN_FROM_DEFAULT",
    "CHAIN_FROM_ROW",
    "ColumnMapping",
    "parse_column_options",
    "read_id_column",
    "read_manifest",
    "resolve_structure_path",
]

#: Universally meaningful, whatever the method.
UNIVERSAL_COLUMNS: Tuple[str, ...] = (
    "sample_id",
    "pdb_path",
    "pdb_id",
    "structure_id",
    "method",
    "chain",
    "label",
)

#: At least one structure column and a chain are always needed; ``method`` may
#: instead come from ``--method`` on the command line, and is checked there.
REQUIRED_COLUMNS: Tuple[str, ...] = ("chain",)

#: Columns whose value is a truth value.
BOOLEAN_COLUMNS: Tuple[str, ...] = ("aa_grouping", "mirror", "sse", "no_contact_filter")
#: Columns whose value is a number of angstroms.
FLOAT_COLUMNS: Tuple[str, ...] = ("size_filter", "contact_cutoff")


def _all_selection_columns() -> Tuple[str, ...]:
    names = set()
    for accepted, _required in SELECTION_ARGUMENTS.values():
        names.update(accepted)
    return tuple(sorted(names))


#: Every header the manifest reader accepts, canonical spelling.
CANONICAL_COLUMNS: Tuple[str, ...] = tuple(
    sorted(
        set(UNIVERSAL_COLUMNS)
        | set(MODIFIER_ARGUMENTS)
        | set(_all_selection_columns())
        | set(CROSS_ARGUMENTS)
    )
)

#: Historical dataset-sheet spellings that resolve to a canonical column.
#:
#: Only names that are unambiguous *as a single column* are here.  The four chain
#: columns all mean "the chain this calculation runs on" in their own row of a
#: long manifest, so they alias to ``chain``; in a wide sheet two of them are
#: populated at once, which is exactly the collision this reader refuses.
COLUMN_ALIASES: Mapping[str, str] = {
    "protein": "pdb_id",
    "pdb": "pdb_path",
    "aa_group": "label",
    "aa_chain": "chain",
    "drug_chain": "chain",
    "carb_chain": "chain",
    "nt_chain": "chain",
    "aa_chain_2": "chain_2",
    "aa_selection": "residues",
    "aa_ind_specific": "residues",
    "aa_ind": "residue_types",
    "aa_ind_mix": "residue_types",
    "aa_atom_selection": "atoms",
    "drug_atom_selected": "atoms",
    "carb_atom_selected": "atoms",
    # The header historical dataset sheets actually carry; their column
    # definitions name the same field carb_atom_selected. Both resolve, and the
    # data is never edited -- the same decision biotsr.io.input_columns makes.
    "car_atom_selected": "atoms",
    "nt_atom_selected": "atoms",
    "drug_name": "entity_name",
    "carb_name": "entity_name",
    "nt_name": "entity_name",
    "drug_id": "entity_id",
    "carb_id": "entity_id",
    "nt_id": "entity_id",
    "aa_selection_2": "residues_2",
    "aa_atom_selection_2": "atom_selection_2",
}

_TRUE = {"1", "true", "t", "yes", "y", "on"}
_FALSE = {"0", "false", "f", "no", "n", "off"}


def _canonical(header: str) -> Optional[str]:
    name = str(header or "").strip()
    if not name:
        return None
    if name in CANONICAL_COLUMNS:
        return name
    return COLUMN_ALIASES.get(name)


def _boolean(column: str, value: str, where: str) -> bool:
    text = value.strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise CLIError(
        "%s: column %r has value %r, which is not a truth value" % (where, column, value),
        "write one of: %s (true) or %s (false); leave it blank for 'not set'"
        % (", ".join(sorted(_TRUE)), ", ".join(sorted(_FALSE))),
    )


def _number(column: str, value: str, where: str) -> float:
    try:
        return float(value.strip())
    except ValueError:
        raise CLIError(
            "%s: column %r has value %r, which is not a number"
            % (where, column, value),
            "it is a distance in angstroms, e.g. 25.0",
        )


def resolve_structure_path(
    row: Mapping[str, Any],
    pdb_dir: str,
    where: str,
) -> Tuple[str, Optional[str]]:
    """``(path, structure_id)`` for one row.

    ``pdb_path`` wins when both are given, and ``pdb_id`` then still supplies the
    structure id -- which is how a dataset whose files are named
    ``download_4f2q.pdb`` keeps reporting ``4F2Q``.
    """
    path = row.get("pdb_path")
    pdb_id = row.get("pdb_id")
    structure_id = row.get("structure_id") or pdb_id or None
    if path:
        return str(path), structure_id
    if not pdb_id:
        raise CLIError(
            "%s: no structure" % where,
            "give either a pdb_path column or a pdb_id column (with --pdb-dir)",
        )
    if not pdb_dir:
        raise CLIError(
            "%s: the manifest identifies structures by pdb_id but no directory "
            "was given" % where,
            "pass --pdb-dir DIR, or use a pdb_path column instead",
        )
    for stem in (str(pdb_id), str(pdb_id).upper(), str(pdb_id).lower()):
        for ext in (".pdb", ".ent"):
            candidate = os.path.join(pdb_dir, stem + ext)
            if os.path.isfile(candidate):
                return candidate, structure_id
    raise CLIError(
        "%s: no structure file for %r in %s" % (where, pdb_id, pdb_dir),
        "looked for %s.pdb and %s.ent (as written, upper case and lower case)"
        % (pdb_id, pdb_id),
    )


class ColumnMapping(object):
    """How each header of one manifest is read.

    ``columns`` maps a canonical field to the header it is read from
    (``{"label": "group"}``, the Python spelling of ``--column label=group``).
    ``keep_columns`` are carried through as metadata -- into ``run.json`` and the
    batch summary -- and are never scientific.  ``ignore_columns`` are dropped.

    Everything else resolves exactly as it always has: a canonical header, or
    one of :data:`COLUMN_ALIASES`, or a ``comment*`` column, which is ignored.
    An unknown header is still an error: a mistyped ``size_fliter`` column must
    not be silently dropped, so dropping or keeping a column has to be asked for.
    """

    def __init__(
        self,
        columns: Optional[Mapping[str, str]] = None,
        keep_columns: Sequence[str] = (),
        ignore_columns: Sequence[str] = (),
    ):
        self.columns = dict(columns or {})
        self.keep_columns = tuple(keep_columns or ())
        self.ignore_columns = tuple(ignore_columns or ())

    def resolve(self, header: Sequence[str], path: str):
        """``(field_by_index, kept_by_index)`` for one header row, or raise.

        Every ambiguity is refused here, before a single row is read, with the
        remedy in the message.
        """
        names = [str(raw or "").strip() for raw in header]
        present = set(n for n in names if n)
        listing = ", ".join(n for n in names if n) or "(none)"
        repeated = sorted({n for n in names if n and names.count(n) > 1})
        if repeated:
            raise CLIError(
                "manifest %s has more than one column named %s"
                % (path, ", ".join(repr(n) for n in repeated)),
                "every header must be unique; rename or remove the duplicate",
            )

        for field_name, source in sorted(self.columns.items()):
            if field_name not in CANONICAL_COLUMNS:
                raise CLIError(
                    "--column %s=%s: %r is not a manifest field"
                    % (field_name, source, field_name),
                    "fields: %s%s" % (", ".join(CANONICAL_COLUMNS),
                                      _suggest(field_name, CANONICAL_COLUMNS)),
                )
            if source not in present:
                raise CLIError(
                    "--column %s=%s: manifest %s has no column %r"
                    % (field_name, source, path, source),
                    "its columns are: %s" % listing,
                )
        by_source: Dict[str, List[str]] = {}
        for field_name, source in self.columns.items():
            by_source.setdefault(source, []).append(field_name)
        for source, fields in sorted(by_source.items()):
            if len(fields) > 1:
                raise CLIError(
                    "column %r is mapped to more than one field: %s"
                    % (source, ", ".join(sorted(fields))),
                    "a column can be read as one field only",
                )
        for option, headers in (("--keep-column", self.keep_columns),
                                ("--ignore-column", self.ignore_columns)):
            for name in headers:
                if name not in present:
                    raise CLIError(
                        "%s %s: manifest %s has no column %r"
                        % (option, name, path, name),
                        "its columns are: %s" % listing,
                    )
                if name in by_source:
                    raise CLIError(
                        "column %r is both mapped (--column %s=%s) and passed to %s"
                        % (name, by_source[name][0], name, option),
                        "give each column one role",
                    )
        both = sorted(set(self.keep_columns) & set(self.ignore_columns))
        if both:
            raise CLIError(
                "column(s) %s are both kept and ignored" % ", ".join(both),
                "give each column one role",
            )
        for name in self.keep_columns:
            if _canonical(name) is not None:
                raise CLIError(
                    "--keep-column %s: %r is already a manifest column (%s)"
                    % (name, name, _canonical(name)),
                    "--keep-column carries a column BioTSR does not otherwise "
                    "read; pass --ignore-column %s to drop it instead" % name,
                )

        ignored = set(self.ignore_columns)
        mapped_sources = set(by_source)
        field_by_index: Dict[int, str] = {}
        kept_by_index: Dict[int, str] = {}
        seen: Dict[str, str] = {}
        unknown: List[str] = []
        for index, name in enumerate(names):
            if not name or name in ignored:
                continue
            if name in mapped_sources:
                field_by_index[index] = by_source[name][0]
                continue
            if name in self.keep_columns:
                kept_by_index[index] = name
                continue
            if name.startswith("comment"):
                continue
            canonical = _canonical(name)
            if canonical is None:
                if name in _REMOVED_COLUMNS:
                    raise CLIError(
                        "manifest %s has a %r column, which this version of "
                        "BioTSR does not read" % (path, name),
                        _REMOVED_COLUMNS[name],
                    )
                unknown.append(name)
                continue
            if canonical in self.columns:
                raise CLIError(
                    "--column %s=%s, but manifest %s also has column %r, which "
                    "is read as %s"
                    % (canonical, self.columns[canonical], path, name, canonical),
                    "drop the mapping, or pass --ignore-column %s" % name,
                )
            field_by_index[index] = canonical
            seen.setdefault(canonical, name)
        if unknown:
            raise CLIError(
                "manifest %s has unrecognised column(s): %s"
                % (path, ", ".join(sorted(unknown))),
                "for each one, either read it as a manifest field with --column "
                "FIELD=HEADER, carry it as metadata with --keep-column HEADER, or "
                "drop it with --ignore-column HEADER.%s Accepted columns: %s "
                "(aliases: %s); columns named comment* are ignored"
                % (
                    "".join(_suggest(u, CANONICAL_COLUMNS, prefix=" %s:" % u)
                            for u in sorted(unknown)),
                    ", ".join(CANONICAL_COLUMNS),
                    ", ".join(sorted(COLUMN_ALIASES)),
                ),
            )
        if not field_by_index:
            raise CLIError("manifest %s has no recognised columns" % path)
        return field_by_index, kept_by_index


def _suggest(name: str, choices: Sequence[str], prefix: str = "") -> str:
    """`` did you mean 'x'?`` -- or nothing, when nothing is close."""
    # A high cutoff on purpose: "size_fliter" should suggest size_filter, but a
    # class-label column called "group" must never be steered to aa_grouping.
    close = difflib.get_close_matches(name, list(choices), n=1, cutoff=0.8)
    if not close:
        return ""
    return "%s did you mean %r?" % (prefix or ";", close[0])


def parse_column_options(values: Sequence[str]) -> Dict[str, str]:
    """``["label=group", ...]`` -> ``{"label": "group"}``, refusing a malformed item."""
    out: Dict[str, str] = {}
    for item in values or ():
        field_name, sep, source = str(item).partition("=")
        field_name, source = field_name.strip(), source.strip()
        if not sep or not field_name or not source:
            raise CLIError(
                "--column %r is not FIELD=HEADER" % item,
                "e.g. --column label=group reads the label from the column "
                "named 'group'",
            )
        if field_name in out and out[field_name] != source:
            raise CLIError(
                "field %r is mapped to two columns: %r and %r"
                % (field_name, out[field_name], source),
                "a field is read from one column only",
            )
        out[field_name] = source
    return out


#: Where a row's chain came from, recorded per row in the batch summary.
CHAIN_FROM_ROW = "row"
CHAIN_FROM_DEFAULT = "default"


def read_manifest(
    path: str,
    *,
    defaults: Optional[Mapping[str, Any]] = None,
    columns: Optional[Mapping[str, str]] = None,
    keep_columns: Sequence[str] = (),
    ignore_columns: Sequence[str] = (),
) -> List[Dict[str, Any]]:
    """Parse a manifest into one normalised value-mapping per row.

    The mappings have canonical keys and typed values, with anything the row left
    blank simply absent -- the same shape ``argparse`` produces, so
    ``spec.build_spec`` validates a manifest row and a command line through one
    code path.

    ``defaults`` fill only what a row leaves blank or has no column for: a row's
    own value always wins.  Each row also carries ``_row`` (its line number),
    ``_where`` (for messages), ``_metadata`` (the kept columns, verbatim) and
    ``_chain_source`` (:data:`CHAIN_FROM_ROW` or :data:`CHAIN_FROM_DEFAULT`).
    """
    if not os.path.isfile(path):
        raise CLIError("manifest not found: %s" % path)
    defaults = dict(defaults or {})
    plan = ColumnMapping(columns, keep_columns, ignore_columns)

    with _open_manifest(path) as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise CLIError("manifest %s is empty" % path)
        except (csv.Error, UnicodeDecodeError) as exc:
            raise CLIError("manifest %s cannot be read as CSV: %s" % (path, exc),
                           "save it as UTF-8 comma-separated text")

        mapping, kept = plan.resolve(header, path)
        named = {i for i, h in enumerate(header) if str(h or "").strip()}

        rows: List[Dict[str, Any]] = []
        for line_number, raw_row in _rows(reader, path):
            if not any(str(cell).strip() for cell in raw_row):
                continue
            where = "%s line %d" % (os.path.basename(path), line_number)
            stray = [i for i, cell in enumerate(raw_row)
                     if i not in named and str(cell).strip()]
            if stray:
                first = stray[0]
                raise CLIError(
                    "%s: value %r is in column %d, which has no header "
                    "(the header has %d column%s)%s"
                    % (where, raw_row[first], first + 1, len(header),
                       "" if len(header) == 1 else "s",
                       "; the row has %d cells" % len(raw_row)
                       if len(raw_row) > len(header) else ""),
                    "a value without a header would be silently dropped, so it "
                    "is refused. If the cell itself contains commas (a list "
                    "such as '1,2,3'), separate the items with ';' or quote the "
                    "cell; otherwise add a header for that column or remove "
                    "the extra cells",
                )
            values: Dict[str, Any] = dict(defaults)
            filled: Dict[str, str] = {}
            for index, column in sorted(mapping.items()):
                if index >= len(raw_row):
                    continue
                text = str(raw_row[index]).strip()
                if not text:
                    continue
                if column in filled:
                    raise CLIError(
                        "%s: columns %r and %r both give %r"
                        % (where, filled[column], header[index], column),
                        "a manifest row is one run; split the wide sheet into "
                        "one row per calculation",
                    )
                filled[column] = header[index]
                if column in BOOLEAN_COLUMNS:
                    values[column] = _boolean(header[index], text, where)
                elif column in FLOAT_COLUMNS:
                    values[column] = _number(header[index], text, where)
                else:
                    values[column] = text
            values["_row"] = line_number
            values["_where"] = where
            values["_metadata"] = {
                name: (raw_row[index] if index < len(raw_row) else "")
                for index, name in sorted(kept.items())
            }
            values["_chain_source"] = (
                CHAIN_FROM_ROW if "chain" in filled
                else CHAIN_FROM_DEFAULT if values.get("chain") else ""
            )
            rows.append(values)

    if not rows:
        raise CLIError("manifest %s has a header but no data rows" % path)
    return rows


def _open_manifest(path: str):
    return open(path, "r", newline="", encoding=MANIFEST_ENCODING)


def _rows(reader, path: str):
    """``(line_number, row)`` pairs, with a CSV syntax error reported by line."""
    line_number = 1
    while True:
        try:
            row = next(reader)
        except StopIteration:
            return
        except (csv.Error, UnicodeDecodeError) as exc:
            raise CLIError(
                "manifest %s cannot be read as CSV after line %d: %s"
                % (path, line_number, exc),
                "save it as UTF-8 comma-separated text",
            )
        line_number = reader.line_num
        yield line_number, row


def read_id_column(path: str, column: Optional[str] = None) -> List[Tuple[int, str]]:
    """``[(line, pdb_id), ...]`` from one column of any CSV, and nothing else.

    For retrieval, which needs only identifiers: no method, no chain and no
    validation of the other columns, so a dataset sheet with labels and notes in
    it can be fetched as it is.  The column is ``column`` when given, otherwise
    the one header that means ``pdb_id`` (``pdb_id`` itself or an alias of it).
    Blank cells are skipped.
    """
    if not os.path.isfile(path):
        raise CLIError("manifest not found: %s" % path)
    with _open_manifest(path) as fh:
        reader = csv.reader(fh)
        try:
            header = [str(h or "").strip() for h in next(reader)]
        except StopIteration:
            raise CLIError("manifest %s is empty" % path)
        if column:
            if column not in header:
                raise CLIError(
                    "manifest %s has no column %r" % (path, column),
                    "its columns are: %s" % ", ".join(h for h in header if h),
                )
            index = header.index(column)
        else:
            candidates = [i for i, h in enumerate(header) if _canonical(h) == "pdb_id"]
            if not candidates:
                raise CLIError(
                    "manifest %s has no PDB id column" % path,
                    "name it with --column pdb_id=HEADER; recognised names are "
                    "pdb_id and %s. Its columns are: %s"
                    % (", ".join(sorted(a for a, c in COLUMN_ALIASES.items()
                                        if c == "pdb_id")),
                       ", ".join(h for h in header if h)),
                )
            if len(candidates) > 1:
                raise CLIError(
                    "manifest %s has more than one PDB id column: %s"
                    % (path, ", ".join(header[i] for i in candidates)),
                    "choose one with --column pdb_id=HEADER",
                )
            index = candidates[0]
        out: List[Tuple[int, str]] = []
        for line_number, raw_row in _rows(reader, path):
            if index < len(raw_row) and str(raw_row[index]).strip():
                out.append((line_number, str(raw_row[index]).strip()))
    return out
