"""Output prefixes: one readable, deterministic name per calculation.

Every calculation's files share one prefix and sit side by side in one shared
output folder (ADR-098)::

    7LF9_D_carbohydrate_NAG_1.keys
    7LF9_D_carbohydrate_NAG_1.meta.json
    7LF9_D_carbohydrate_NAG_1.triplets
    7LF9_D_carbohydrate_NAG_1.run.json

The prefix names, in order: the **structure**, the **chain** (both chains for a
cross method, in argument order), the **method**, the **target** -- the molecule
and its residue identifier, ``allcopies`` when no identifier was given, residue
types, named residues or atoms, each cross side's narrowing prefixed by its
chain -- and the **modifiers and settings** that change the numbers (``mirror``,
``aagroup``, ``size25``, ``sse``, ``cutoff4``, ``filter-side_a_only``,
``nocontactfilter``).  The bin counts are not in the name; they are in the
metadata, with the full bin definitions.

The prefix is a function of the request alone, so the same request always gets
the same name, on any machine and in any batch.  Where the readable form cannot
represent the request exactly -- a character outside ``A-Z a-z 0-9 . _ + -``,
a lower-case chain identifier (which a case-insensitive file system cannot tell
from upper case), or a name longer than :data:`MAX_READABLE` characters -- a
deterministic request fingerprint is appended (``_h<10 hex digits>``), so two
different requests never share a name.  A caller may always choose the prefix
instead (``sample_id``).
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "SAFE_CHARACTERS",
    "MAX_READABLE",
    "request_prefix",
    "request_fingerprint",
    "safe_token",
]

SAFE_CHARACTERS = frozenset(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._+-"
)

#: Longer readable prefixes are shortened to this many characters, and the
#: request fingerprint is appended.
MAX_READABLE = 96

#: Arguments that narrow the first protein side of a cross method and the second.
_SIDE_1 = ("residue_types", "residues", "atom_selection", "calpha_residues")
_SIDE_2 = ("residue_types_2", "residues_2", "atom_selection_2", "calpha_residues_2")


def _items(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = value.replace(",", ";").split(";")
    return [str(v).strip() for v in value if str(v).strip()]


def safe_token(text: str) -> Tuple[str, bool]:
    """``(token, exact)``: ``text`` with unsafe characters replaced by ``-``.

    ``exact`` is false when anything was replaced, so the caller knows the
    token alone no longer identifies the request.
    """
    out = []
    exact = True
    for ch in str(text):
        if ch in SAFE_CHARACTERS:
            out.append(ch)
        elif ch == "'":
            # A prime, as in nucleotide atom names (C2'): readable as "p", and
            # inexact, so the fingerprint keeps C2' and C2P apart.
            exact = False
            out.append("p")
        else:
            exact = False
            if not out or out[-1] != "-":
                out.append("-")
    token = "".join(out).strip("-")
    if not token:
        return "x", False
    return token, exact


def _joined(values: Sequence[str]) -> str:
    return "+".join(values)


def _compact(name: str, rest: str) -> Tuple[str, bool]:
    """``NAME`` + ``rest``; exact only when the name has no digit to blur the join."""
    return name + rest, name.isalpha()


def _residue_tokens(values: Sequence[str]) -> Tuple[str, bool]:
    """``LYS_6;GLU18`` -> ``LYS6+GLU18`` (one residue, one token)."""
    tokens, exact = [], True
    for value in values:
        name, sep, rest = value.partition("_")
        if not sep:
            tokens.append(value)
            continue
        token, ok = _compact(name, rest)
        tokens.append(token)
        exact = exact and ok
    return _joined(tokens), exact


def _atom_tokens(values: Sequence[str]) -> Tuple[str, bool]:
    """``LYS_48_NZ;HIS_68_NE2`` -> ``LYS48.NZ+HIS68.NE2``."""
    tokens, exact = [], True
    for value in values:
        parts = value.split("_", 2)
        if len(parts) != 3:
            tokens.append(value)
            exact = False
            continue
        token, ok = _compact(parts[0], parts[1])
        tokens.append("%s.%s" % (token, parts[2]))
        exact = exact and ok
    return _joined(tokens), exact


def _modifier_tokens(modifiers: Mapping[str, Any], cross: Mapping[str, Any]) -> List[str]:
    out = []
    if modifiers.get("aa_grouping"):
        out.append("aagroup")
    if modifiers.get("mirror"):
        out.append("mirror")
    if modifiers.get("size_filter") is not None:
        out.append("size%g" % float(modifiers["size_filter"]))
    if modifiers.get("sse"):
        out.append("sse")
    if cross.get("contact_cutoff") is not None:
        out.append("cutoff%g" % float(cross["contact_cutoff"]))
    if cross.get("contact_filter"):
        out.append("filter-%s" % cross["contact_filter"])
    if cross.get("no_contact_filter"):
        out.append("nocontactfilter")
    return out


def request_fingerprint(request: Mapping[str, Any]) -> str:
    """Ten hex digits of SHA-256 over the canonical request."""
    blob = json.dumps(_canonical(request), sort_keys=True, separators=(",", ":"),
                      default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:10]


def _canonical(request: Mapping[str, Any]) -> Mapping[str, Any]:
    def clean(mapping):
        return {k: v for k, v in sorted((mapping or {}).items())
                if v is not None and v is not False and v != ""}
    return {
        "structure_id": str(request.get("structure_id") or ""),
        "method": request.get("method"),
        "chain": request.get("chain"),
        "chain_2": request.get("chain_2"),
        "sweep": bool(request.get("sweep")),
        "selection": clean(request.get("selection")),
        "modifiers": clean(request.get("modifiers")),
        "cross": clean(request.get("cross")),
    }


def _target_tokens(request: Mapping[str, Any]) -> List[Tuple[str, bool]]:
    """``[(token, exact)]``: what the calculation selected, in prefix order."""
    selection = dict(request.get("selection") or {})
    cross = dict(request.get("cross") or {})
    args = dict(selection)
    args.update(cross)
    out: List[Tuple[str, bool]] = []

    if request.get("sweep"):
        name = args.get("residue_name")
        out.append(("perresidue-%s" % str(name).upper() if name else "perresidue", True))
    elif args.get("residue") is not None and args.get("residue_name"):
        out.append(("%s_%s" % (str(args["residue_name"]).upper(), args["residue"]), True))

    if args.get("entity_name"):
        name = str(args["entity_name"]).upper()
        entity_id = args.get("entity_id")
        out.append(("%s_%s" % (name, entity_id) if entity_id not in (None, "")
                    else "%s_allcopies" % name, True))
        if args.get("atoms"):
            out.append(("atoms-" + _joined(_items(args["atoms"])), True))
    elif args.get("atoms"):
        # amino_acid_atom_selection: RESIDUE_ID_ATOM tokens.
        token, exact = _atom_tokens(_items(args["atoms"]))
        out.append(("atoms-" + token, exact))
    if args.get("entity_names"):
        out.append(("names-" + _joined([n.upper() for n in _items(args["entity_names"])]),
                    True))

    is_cross = request.get("chain_2") is not None
    if not is_cross:
        if args.get("residue_types"):
            out.append((_joined([t.upper() for t in _items(args["residue_types"])]), True))
        if args.get("residues"):
            values = _items(args["residues"])
            if request.get("method") == "aa_selection":
                out.append(("res-" + _joined(values), True))
            else:
                out.append(_residue_tokens(values))
        return out

    chain, chain_2 = str(request.get("chain")), str(request.get("chain_2"))
    same = chain == chain_2
    for names, marker in ((_SIDE_1, "s1" if same else chain),
                          (_SIDE_2, "s2" if same else chain_2)):
        for name in names:
            values = _items(args.get(name))
            if not values:
                continue
            exact = True
            if name.startswith("residue_types"):
                token = _joined([v.upper() for v in values])
            elif name.startswith("atom_selection"):
                token, exact = _atom_tokens(values)
                token = "atoms-" + token
            elif name.startswith("calpha_residues"):
                token = "res-" + _joined(values)
            else:
                token, exact = _residue_tokens(values)
            out.append(("%s-%s" % (marker, token), exact))
    return out


def request_prefix(request: Mapping[str, Any]) -> str:
    """The output prefix of one calculation request.

    ``request`` holds ``structure_id``, ``method``, ``chain``, ``chain_2``
    (``None`` for a one-sided method), ``selection``, ``cross`` and
    ``modifiers`` (the arguments as given) and ``sweep`` (true for the
    per-residue sweep of ``amino_acid_atomic``).
    """
    parts: List[Tuple[str, bool]] = []
    parts.append(safe_token(str(request.get("structure_id") or "structure")))
    for key in ("chain", "chain_2"):
        value = request.get(key)
        if value is None:
            continue
        token, exact = safe_token(str(value))
        # A case-insensitive file system cannot tell chain "a" from chain "A".
        parts.append((token, exact and token == token.upper()))
    parts.append(safe_token(str(request.get("method"))))
    for token, exact in _target_tokens(request):
        safe, ok = safe_token(token)
        parts.append((safe, ok and exact))
    for token in _modifier_tokens(dict(request.get("modifiers") or {}),
                                  dict(request.get("cross") or {})):
        parts.append(safe_token(token))
    readable = "_".join(token for token, _ in parts)
    exact = all(ok for _, ok in parts)
    if exact and len(readable) <= MAX_READABLE:
        return readable
    short = readable[:MAX_READABLE].rstrip("_+-.")
    return "%s_h%s" % (short, request_fingerprint(request))
