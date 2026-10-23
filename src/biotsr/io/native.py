"""The native serializers -- the standard output format.

One calculation is written as files sharing one **prefix**, side by side::

    7LF9_D_carbohydrate_NAG_1.keys        key frequencies
    7LF9_D_carbohydrate_NAG_1.meta.json   the metadata sidecar
    7LF9_D_carbohydrate_NAG_1.triplets    one line per triangle (optional)

The prefix is chosen by the caller (``save(..., prefix=)``); the command line
builds it from the request (``biotsr.io.naming``) and :func:`default_prefix`
builds one from a result.  The file names carry no bin counts: the metadata
sidecar records the full key space, including both bin schemes and their
content fingerprints, so nothing has to be recovered by parsing a filename.
(Release candidates up to 1.0.0rc3 wrote ``<id>_<chain>.keys_29_35`` inside one
directory per sample; that layout is not read by this version.)

**No four-character identifier assumption exists here.**  The structure id is
inserted whole.

**Both data files are tab-separated text with one header line**, UTF-8, ``\n``
line endings.  The key file reports **category and frequency explicitly**, one
row per ``(key, category)``::

    key      category   frequency
    123456   regular    10
    -123456  mirror     6

For a run with no category concept every row reads ``regular``, so this costs a
column and no ambiguity.  For a mirror run the sign already carries the
orientation, but the category is stated anyway -- because sign alone cannot: the
legacy ``key == 0`` edge case negates to itself, so a mirrored zero-key and an
unmirrored one are indistinguishable by sign.  The structured result remains the
source of truth; this file is an export of it.

Under AA grouping the ``label_name`` and ``residue`` columns differ -- the label
name is the group's representative, the residue is what is physically there.  Both
are written, because the legacy format keeps only the first and thereby misreports
which residue a grouped triangle contains.

Determinism is a requirement, not a nicety:

* the key file is sorted **numerically by key**, then alphabetically by category
  (the legacy file is in first-seen order, which is hash-independent but
  arbitrary);
* the triplet file is in candidate-enumeration order, which is fixed by the
  point set's file order.

The ``.meta.json`` sidecar carries the full ``KeySpace`` and
``SamplingProfile``, both fingerprints, the counts, the selection record and
the provenance of the run.
"""

from __future__ import annotations

import json
import os

from ..domain.result import TSRResult, TripletRecord

__all__ = [
    "ENCODING",
    "KEYS_SUFFIX",
    "TRIPLETS_SUFFIX",
    "METADATA_SUFFIX",
    "default_prefix",
    "keys_filename",
    "triplets_filename",
    "metadata_filename",
    "prefix_of",
    "format_native_triplet",
    "NativeTripletStreamWriter",
    "triplet_header",
    "KEY_COLUMNS",
    "keys_header",
    "write_keys",
    "write_triplets",
    "write_metadata",
    "save",
    "TRIPLET_COLUMNS",
    "SSE_ABSENT",
]

#: Every native file is read and written as UTF-8, explicitly.
#:
#: ``open(path, "w")`` with no ``encoding=`` resolves through
#: ``locale.getpreferredencoding(False)``, so the same result written on a
#: machine with ``LANG=C`` and on one with ``LANG=en_US.UTF-8`` could disagree in
#: bytes -- and a sidecar carrying a non-ASCII character (a structure path, which
#: is user-controlled and reaches ``source.path``) could raise
#: ``UnicodeEncodeError`` on write or come back mojibake on read.  A serialization
#: format whose bytes depend on the machine's locale is not a format.
#:
#: Byte-identical for everything the package has ever written: UTF-8 agrees with
#: ASCII, and every token these writers emit is ASCII.
ENCODING = "utf-8"

#: The three file kinds of one calculation, appended to its prefix.  Exported so
#: that anything recognising a native file -- the reader, and the command line
#: when it collects a folder -- uses the one spelling.
KEYS_SUFFIX = ".keys"
TRIPLETS_SUFFIX = ".triplets"
METADATA_SUFFIX = ".meta.json"

KEY_COLUMNS = ("key", "category", "frequency")

#: Per-vertex identity is written in full: structure, side, residue name,
#: chain, sequence label (number + insertion code) and atom name.
#:
#: ``structure`` and ``side`` are unconditional
#: rather than cross-only: a format whose columns depend on which method produced
#: it cannot be read without first knowing the method, which is the property the
#: legacy formats have and this one exists to avoid. A one-sided result writes
#: ``primary`` in every side column, which is the truth about it.
TRIPLET_COLUMNS = (
    "key",
    "label1",
    "label_name1",
    "structure1",
    "side1",
    "residue1",
    "chain1",
    "seq1",
    "atom1",
    "label2",
    "label_name2",
    "structure2",
    "side2",
    "residue2",
    "chain2",
    "seq2",
    "atom2",
    "label3",
    "label_name3",
    "structure3",
    "side3",
    "residue3",
    "chain3",
    "seq3",
    "atom3",
    "theta_class",
    "theta",
    "theta_raw",
    "distance_class",
    "max_distance",
    "x1",
    "y1",
    "z1",
    "x2",
    "y2",
    "z2",
    "x3",
    "y3",
    "z3",
    "category",
    #: The secondary-structure annotation, per vertex and then classified.
    #:
    #: Unconditional, for the same reason ``side`` is: a reader must not have to
    #: know which method wrote the file in order to know how wide it is. A run
    #: without SSE writes ``-`` in all seven, which is NOT the same token as
    #: ``NONE`` -- "this residue is in no HELIX/SHEET record" and "nobody looked"
    #: are different facts, and the one thing the legacy three-state collapse
    #: cannot express is the difference between them.
    "sse_state1",
    "sse_group1",
    "sse_state2",
    "sse_group2",
    "sse_state3",
    "sse_group3",
    "sse_category",
)

#: Written in every SSE column of a run that was not annotated.
SSE_ABSENT = "-"


def default_prefix(result: TSRResult) -> str:
    """A descriptive prefix for ``result``, by the command line's grammar.

    ``<structure>_<chain>[_<chain_2>]_<method>[_<target>][_<modifiers>]``,
    built from the result's own provenance (``biotsr.io.naming``).  The command
    line names a calculation from its request before computing it; for the
    same request the two agree.
    """
    from .naming import request_prefix

    return request_prefix(_request_of(result))


def _request_of(result: TSRResult):
    """The naming request a result describes, from its recorded provenance."""
    sides = (result.annotations.get("selection_request") or {}).get("sides") or []
    modifiers = {}
    for name in result.annotations.get("modifiers") or ():
        key = str(name).split("(")[0].split("=")[0]
        if key in ("mirror", "aa_grouping", "sse"):
            modifiers[key] = True
    if result.sampling.size_filter is not None:
        modifiers["size_filter"] = result.sampling.size_filter
    if result.annotations.get("sse"):
        modifiers["sse"] = True
    request = {
        "structure_id": result.structure_id.value,
        "method": result.method,
        "chain": result.chain,
        "chain_2": None,
        "selection": {},
        "cross": {},
        "modifiers": modifiers,
    }
    if result.sampling.cross_method is None:
        if sides:
            args = dict(sides[0].get("arguments") or {})
            if "atom_selection" in args:
                args["atoms"] = args.pop("atom_selection")
            request["selection"] = args
        return request
    # ``chain`` is the entity side's chain for an entity family (the caller's
    # argument order), and side a's for protein-protein methods.
    entity = [x for x in sides if x.get("selector") in
              ("cross_entity_atoms", "cross_entity_chain_atoms")]
    first = entity[0] if entity else (sides[0] if sides else {})
    others = [x for x in sides if x is not first]
    second = others[0] if others else {}
    request["chain"] = first.get("chain", result.chain)
    request["chain_2"] = second.get("chain")
    cross = {}
    for side, suffix in ((first, ""), (second, "_2")):
        args = side.get("arguments") or {}
        for key in ("entity_name", "entity_id", "atoms", "entity_names"):
            if args.get(key) and not suffix:
                cross[key] = args[key]
        if args.get("residue_types"):
            cross["residue_types" + suffix] = args["residue_types"]
        if args.get("atom_selection"):
            cross["atom_selection" + suffix] = args["atom_selection"]
        if args.get("residues"):
            key = ("calpha_residues" if side.get("scope") == "residue_list"
                   else "residues")
            cross[key + suffix] = args["residues"]
    if result.sampling.contact_cutoff is not None:
        # Named only when it is not the method's own reference cutoff, as the
        # command line names it only when the caller set it.
        from ..methods.builtin.cross import cross_recipe

        if result.sampling.contact_cutoff != cross_recipe(result.method).default_cutoff:
            cross["contact_cutoff"] = result.sampling.contact_cutoff
    if result.sampling.contact_filter not in (None, "both_sides"):
        if result.sampling.contact_filter == "none":
            cross["no_contact_filter"] = True
        else:
            cross["contact_filter"] = result.sampling.contact_filter
    request["cross"] = cross
    return request


def keys_filename(prefix: str) -> str:
    return prefix + KEYS_SUFFIX


def triplets_filename(prefix: str) -> str:
    return prefix + TRIPLETS_SUFFIX


def metadata_filename(prefix: str) -> str:
    return prefix + METADATA_SUFFIX


def prefix_of(path: str) -> str:
    """The prefix of a native file name (``.keys``, ``.triplets``, ``.meta.json``)."""
    base = os.path.basename(path)
    for suffix in (METADATA_SUFFIX, KEYS_SUFFIX, TRIPLETS_SUFFIX):
        if base.endswith(suffix) and len(base) > len(suffix):
            return base[: -len(suffix)]
    raise ValueError(
        "%r is not a native result file: its name does not end in %s, %s or %s"
        % (path, KEYS_SUFFIX, METADATA_SUFFIX, TRIPLETS_SUFFIX)
    )


def keys_header() -> str:
    return "\t".join(KEY_COLUMNS) + "\n"


def write_keys(result: TSRResult, path: str) -> str:
    """Write ``key<TAB>category<TAB>frequency``, one row per (key, category).

    Ordered numerically by key, then alphabetically by category, so the file is
    byte-stable across runs regardless of dict traversal order.
    """
    with open(path, "w", encoding=ENCODING) as fh:
        fh.write(keys_header())
        for key in sorted(result.keys):
            by_category = result.keys[key].by_category
            for category in sorted(by_category):
                fh.write("%d\t%s\t%d\n" % (key, category, by_category[category]))
    return path


def triplet_header() -> str:
    return "\t".join(TRIPLET_COLUMNS) + "\n"


def format_native_triplet(t: TripletRecord) -> str:
    """One native triplet line.

    Floats go through ``repr`` so every value round-trips exactly; the pre-fold
    ``theta_raw`` sits alongside the folded ``theta``, which no legacy format
    stores at all.
    """
    i0, i1, i2 = t.identities
    fields = [str(t.key)]
    for label, label_name, ident, side in zip(
        t.labels, t.label_names, (i0, i1, i2), t.sides
    ):
        fields.extend(
            [
                str(label),
                label_name,
                ident.structure.value,
                side,
                ident.residue_name,
                ident.chain.value,
                ident.residue.legacy_seq_label,
                ident.atom_name,
            ]
        )
    fields.extend(
        [
            str(t.theta_class),
            repr(t.theta),
            repr(t.theta_raw),
            str(t.distance_class),
            repr(t.max_distance),
        ]
    )
    for coord in t.coordinates:
        fields.extend(repr(v) for v in coord)
    fields.append(t.category)
    if t.sse is None:
        fields.extend([SSE_ABSENT] * 7)
    else:
        for state, group in zip(t.sse.states, t.sse.groups):
            fields.extend([state, group])
        fields.append(t.sse.category)
    return "\t".join(fields) + "\n"


class NativeTripletStreamWriter:
    """A ``triplet_sink`` that writes native lines without retaining records."""

    def __init__(self, path: str):
        self.path = path
        self._fh = open(path, "w", encoding=ENCODING)
        self._fh.write(triplet_header())
        self.n_written = 0

    def write(self, record: TripletRecord) -> None:
        self._fh.write(format_native_triplet(record))
        self.n_written += 1

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> "NativeTripletStreamWriter":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def write_triplets(result: TSRResult, path: str) -> str:
    """Write one header line plus one line per triangle, in enumeration order."""
    if result.triplets is None:
        raise ValueError("this result carries no triplet records")
    with open(path, "w", encoding=ENCODING) as fh:
        fh.write(triplet_header())
        for t in result.triplets:
            fh.write(format_native_triplet(t))
    return path


def write_metadata(result: TSRResult, path: str) -> str:
    with open(path, "w", encoding=ENCODING) as fh:
        json.dump(result.metadata(), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path


def save(result: TSRResult, out_dir: str, *, prefix: str = None,
         triplets: bool = True):
    """Write ``<prefix>.keys``, ``<prefix>.meta.json`` and optionally ``<prefix>.triplets``.

    ``prefix`` defaults to :func:`default_prefix`.  Returns the paths written,
    keyed ``keys``, ``metadata`` and ``triplets``.  Existing files of the same
    names are replaced; the command line, which shares one folder between many
    calculations, stages and installs files itself instead (``biotsr.cli.output``).
    """
    os.makedirs(out_dir, exist_ok=True)
    prefix = prefix or default_prefix(result)
    written = {}
    written["keys"] = write_keys(result, os.path.join(out_dir, keys_filename(prefix)))
    written["metadata"] = write_metadata(
        result, os.path.join(out_dir, metadata_filename(prefix)))
    if triplets and result.triplets is not None:
        written["triplets"] = write_triplets(
            result, os.path.join(out_dir, triplets_filename(prefix))
        )
    return written
