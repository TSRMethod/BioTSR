"""The atomic label table: chemical element -> integer label, loaded as a package resource.

Every atomic method (amino-acid atomic, residue-type, drug, carbohydrate,
nucleotide and the atomic cross methods) labels a heavy atom by **its element**,
read from the structure's own element field -- PDB columns 77-78, never the
charge columns 79-80 -- case-insensitively (ADR-113).  The atom name is not
consulted: it stays the atom's identity in selections, records and outputs, but
it never chooses a label.

Why the element field
---------------------

The codes are element classes.  They come from the TSR group's historical atom
lexicon (``CarbohydratesTSR``, ADR-100), where every one of the 109 codes
belongs to exactly one element -- carbon 106, nitrogen 108, oxygen 109,
phosphorus 110, sulfur 111, calcium 99, holmium 37 ... -- and every standard
amino-acid and nucleotide atom name was listed on its element's code.  The
lexicon reached that code through the atom name, which made the label depend on
naming conventions: an ATP ``PA`` was read as protactinium (55), ``PB`` as lead
(56), a ligand ``CL2`` as carbon, a chlorine called ``CL35`` could not be
labelled at all, and which of two rows for one name came later decided some
labels.  The element field states the chemistry directly, so the label is now
taken from it.  No code was renumbered: each element keeps the code the
lexicon gave it.

The table
---------

``element_labels_v1.csv`` has two columns, ``element,label``: 110 element
symbols with 109 distinct labels, 4..112.  Chlorine is 93 (the lexicon carried
it only on its ``CL`` rows).  Hydrogen and deuterium share 107, the one
documented exception to one-label-per-element (:data:`SHARED_LABELS`); both are
excluded by every atomic admission before any atom is labelled, so neither ever
labels a point.  The superheavy elements the lexicon never listed (Ds, Rg, Cn,
Nh, Fl, Mc, Lv, Ts, Og) have no label.

The radix stays **112** -- the historical key space and the encoding radix of
every atomic family (ADR-032) -- and is not reduced to the number of labels: the
labels fit one window of 112 consecutive integers, which is what keeps the
encoding collision-free (ADR-107).

An atom whose element field is blank, or names no element in the table (``X``,
an unknown atom), cannot be labelled, and a calculation that selected it is
refused, naming it (ADR-108).  Nothing is repaired or inferred: there is no
element guessed from a name, no bond-length check and no correction.

The file is the one editable copy.  :func:`parse_element_labels` refuses a
malformed one -- a wrong header, a row with the wrong number of fields, a
symbol not written in chemical case, a label that is not a positive integer, a
symbol listed twice (conflicting or not), or two elements sharing a label
outside :data:`SHARED_LABELS` -- rather than loading part of it.
"""

from __future__ import annotations

import csv
import hashlib
import io
import os
import re
from typing import Dict, Tuple

from .schemes import LabelDomain, LabelScheme, LabelSource

__all__ = [
    "ELEMENT_LABELS_RESOURCE",
    "ATOMIC_LABEL_RADIX",
    "ELEMENT_RESOLVER",
    "SHARED_LABELS",
    "ElementTableError",
    "parse_element_labels",
    "load_element_labels",
    "element_label_scheme",
    "ELEMENT_LABELS_112",
]

#: The packaged element -> label table.
ELEMENT_LABELS_RESOURCE = "element_labels_v1.csv"

#: ``numOfLabels`` in every historical atomic implementation: the encoder
#: radix, not the number of labels.
ATOMIC_LABEL_RADIX = 112

#: How an atom becomes a label, recorded in the scheme's content fingerprint
#: and so in every atomic key space: the element field (PDB columns 77-78),
#: case-insensitive, looked up in the table; no atom-name lookup.
ELEMENT_RESOLVER = "element_field_v1"

#: Elements that share one label on purpose.  Hydrogen and deuterium keep the
#: historical lexicon's common code; neither labels a point, because every
#: atomic admission excludes them first.
SHARED_LABELS: Tuple[Tuple[str, ...], ...] = (("H", "D"),)

_HEADER = ["element", "label"]
_SYMBOL = re.compile(r"^[A-Z][a-z]?$")
_LABEL = re.compile(r"^[1-9][0-9]*$")

_RESOURCE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")


class ElementTableError(ValueError):
    """A malformed element label table."""


def _read_resource(name: str) -> str:
    try:  # Python 3.9+: the packaged-resource path, which is the one that matters.
        from importlib.resources import files  # type: ignore[attr-defined]

        return (files(__package__) / "resources" / name).read_text(encoding="utf-8")
    except Exception:
        # A source checkout with no package metadata still has to work.
        with open(os.path.join(_RESOURCE_DIR, name), encoding="utf-8") as fh:
            return fh.read()


def parse_element_labels(text: str, source: str = "element table") -> Dict[str, int]:
    """``{element symbol: label}`` from the text of an ``element,label`` table.

    Refuses, naming the line, anything that is not exactly the documented
    format; nothing is skipped or corrected.
    """
    rows = list(csv.reader(io.StringIO(text)))
    if not rows or rows[0] != _HEADER:
        raise ElementTableError(
            "%s: the first line must be exactly %r, got %r"
            % (source, ",".join(_HEADER), ",".join(rows[0]) if rows else ""))
    table: Dict[str, int] = {}
    seen_upper: Dict[str, Tuple[str, int, int]] = {}
    for line, row in enumerate(rows[1:], start=2):
        if len(row) != 2:
            raise ElementTableError(
                "%s line %d: expected 2 fields (element,label), got %d: %r"
                % (source, line, len(row), ",".join(row)))
        symbol, label = row
        if not _SYMBOL.match(symbol):
            raise ElementTableError(
                "%s line %d: %r is not an element symbol in chemical case "
                "(one capital letter, optionally one small letter)"
                % (source, line, symbol))
        if not _LABEL.match(label):
            raise ElementTableError(
                "%s line %d: label %r of %s is not a positive integer"
                % (source, line, label, symbol))
        code = int(label)
        previous = seen_upper.get(symbol.upper())
        if previous is not None:
            raise ElementTableError(
                "%s line %d: element %s is listed again (first on line %d, "
                "label %d; here label %d)%s; each element has exactly one row"
                % (source, line, symbol, previous[1], previous[2], code,
                   " -- a conflicting duplicate" if previous[2] != code else ""))
        seen_upper[symbol.upper()] = (symbol, line, code)
        table[symbol] = code
    if not table:
        raise ElementTableError("%s lists no elements" % source)
    by_label: Dict[int, list] = {}
    for symbol, code in table.items():
        by_label.setdefault(code, []).append(symbol)
    allowed = [set(group) for group in SHARED_LABELS]
    for code, symbols in sorted(by_label.items()):
        if len(symbols) > 1 and not any(set(symbols) <= group for group in allowed):
            raise ElementTableError(
                "%s: elements %s share label %d; only %s may share a label"
                % (source, ", ".join(sorted(symbols)), code,
                   "; ".join("/".join(g) for g in SHARED_LABELS)))
    for group in SHARED_LABELS:
        codes = {table[s] for s in group if s in table}
        if len(codes) > 1:
            raise ElementTableError(
                "%s: %s must share one label, got %s"
                % (source, "/".join(group), sorted(codes)))
    return table


def load_element_labels(
    resource: str = ELEMENT_LABELS_RESOURCE,
) -> Tuple[Dict[str, int], str]:
    """Return ``(element -> label, resource sha256)`` for a packaged table."""
    text = _read_resource(resource)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return parse_element_labels(text, resource), digest


def element_label_scheme(
    name: str = "element_labels_112",
    resource: str = ELEMENT_LABELS_RESOURCE,
    radix: int = ATOMIC_LABEL_RADIX,
) -> LabelScheme:
    """Build the atomic label scheme from an element table."""
    table, digest = load_element_labels(resource)
    distinct = len(set(table.values()))
    return LabelScheme(
        name=name,
        table=table,
        n_labels=radix,
        domain=LabelDomain.ATOM_ELEMENT,
        label_source=LabelSource.ELEMENT,
        resolver=ELEMENT_RESOLVER,
        description=(
            "Atomic labels by element from %s (sha256 %s): %d elements, %d "
            "distinct labels (H and D share one), encoded at radix %d. An atom "
            "is labelled by its element field (PDB columns 77-78, "
            "case-insensitive); its atom name is not consulted."
            % (resource, digest[:12], len(table), distinct, radix)
        ),
    )


#: **The atomic label scheme**, shared unchanged by every atomic method:
#: amino-acid atomic, residue-type, drug, carbohydrate, nucleotide and the
#: atomic cross methods.
ELEMENT_LABELS_112 = element_label_scheme()
