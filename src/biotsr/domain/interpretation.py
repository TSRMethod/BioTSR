"""Structure-interpretation semantics: which *reading* of a structure file produced a result.

A parser-backend audit constructed a case that neither scientific identity
catches (ADR-079).  On ``1KX5`` chain ``I`` with the
manganese ion named explicitly, two results have:

* the **same** ``KeySpace`` fingerprint,
* the **same** ``SamplingProfile`` fingerprint,
* and **35 of 35 triangles carrying different keys**,

because the two readings disagree about what the atom is *called* -- one reads
the PDB atom-name columns ``[12:16]``, the other ``[13:16]``, so ``MN`` becomes
``N`` and a manganese is labelled nitrogen.  The analysis layer combined them
and reported a Jaccard of 0.0 without a word.

Why this is **not** a ``KeySpace`` field
----------------------------------------

The obvious repair -- put ``parser_backend`` into ``KeySpace`` -- is wrong twice
over, and deliberately not done (ADR-079):

* **A backend name is an implementation fact, not a scientific semantic.** Two
  independently written readers that agree, byte for byte, about what every atom
  is called produce *scientifically identical* results.  Keying identity on the
  module that did the reading would declare them incomparable, which is false.
  What matters is the **interpretation**, not the implementation.
* **``KeySpace`` is the wrong layer, and changing it is not free.**
  ``KeySpace.from_dict`` refuses unknown fields by design, so adding one would
  make every result saved before this change unreadable, and would change every
  ``KeySpace`` fingerprint ever computed -- including those in frozen references.
  The interpretation does not change *what a key means*; identical keys under
  two interpretations denote the same geometry of the same label triple.  It
  changes *which atoms were there to be labelled*.

So this module follows the shape ADR-063 established for SSE strand identity:
**the semantic is recorded, kept out of both identities, and guarded at the
boundary where combining it would actually be wrong.**

The semantic, named
-------------------

``structure_interpretation`` is the answer to "which atom-name semantics did the
reader apply?".  It is derived from -- and strictly coarser than -- the recorded
``SourceInfo.parser_backend``:

======================================  =========================
``parser_backend``                      ``structure_interpretation``
======================================  =========================
``pdb_legacy``                          ``legacy_columns_v1``
``pdb_legacy+float32_coords``           ``legacy_columns_v1``
======================================  =========================

Both rows map to one interpretation on purpose: coordinate precision is already
a ``KeySpace`` field (``numerical_precision``), it is not an atom-naming
semantic, and two runs differing only there are comparable *as interpretations*.
That coarsening is the property that makes this durable -- a future second
implementation of the same reading would map here too and stay comparable.

BioTSR v1.0 ships exactly one reader, so every standard result carries
``legacy_columns_v1`` and nothing is refused in practice.  The guard exists
**before** the second interpretation does, which is the only time it can be
added without invalidating anything.
"""

from __future__ import annotations

from typing import Any, Mapping, Optional, Sequence

from .errors import IncompatibleRepresentation

__all__ = [
    "STRUCTURE_INTERPRETATION",
    "LEGACY_COLUMNS_V1",
    "SPEC_COLUMNS_V1",
    "INTERPRETATION_BY_BACKEND",
    "IncompatibleStructureInterpretation",
    "interpretation_of_backend",
    "structure_interpretation_of",
    "require_common_structure_interpretation",
]

#: The annotation/provenance key, defined once so a writer and a reader cannot
#: disagree about its spelling.
STRUCTURE_INTERPRETATION = "structure_interpretation"

#: The atom-name semantics BioTSR v1.0 reads: the legacy column slices.
LEGACY_COLUMNS_V1 = "legacy_columns_v1"

#: The spec-offset semantics.  **Named, not shipped** -- no backend in BioTSR
#: produces it (ADR-079).  It is defined here so that
#: the day a second interpretation exists, results already in the wild are
#: distinguishable from it rather than silently merged with it.
SPEC_COLUMNS_V1 = "spec_columns_v1"

#: Backend base name -> interpretation.  The key is the backend name with any
#: ``+suffix`` qualifier removed, because a qualifier records how the
#: coordinates were held, not how the atom names were read.
INTERPRETATION_BY_BACKEND: Mapping[str, str] = {
    "pdb_legacy": LEGACY_COLUMNS_V1,
    "pdb_strict": SPEC_COLUMNS_V1,
}


class IncompatibleStructureInterpretation(IncompatibleRepresentation):
    """Raised when results read under different atom-name semantics would be combined.

    The sibling of ``IncompatibleKeySpace``, ``IncompatibleSampling`` and
    ``IncompatibleSSEStratification``, for the axis the first two deliberately do
    not cover: two results can share both identity fingerprints exactly and still
    disagree about every key, because they disagree about what the atoms are
    called (ADR-079).
    """


def interpretation_of_backend(parser_backend: str) -> str:
    """Map a recorded ``parser_backend`` string to its interpretation semantics.

    Strips any ``+qualifier`` (``pdb_legacy+float32_coords``), which records
    coordinate handling rather than atom naming.  An unrecognized backend is an
    error rather than a default: guessing that an unknown reader means the same
    thing as the shipped one is exactly the silent merge this module exists to
    prevent.
    """
    base = (parser_backend or "").split("+", 1)[0].strip()
    if not base:
        raise IncompatibleStructureInterpretation(
            "this result records no parser backend, so which atom-name semantics "
            "produced it cannot be established. It was written by a build that "
            "predates structure-interpretation provenance."
        )
    try:
        return INTERPRETATION_BY_BACKEND[base]
    except KeyError:
        raise IncompatibleStructureInterpretation(
            "parser backend %r has no declared structure interpretation, so it "
            "cannot be compared with anything. Declare it in "
            "INTERPRETATION_BY_BACKEND, mapping it to an existing interpretation "
            "if it reads atom names identically to one, or to a new one if it "
            "does not." % parser_backend
        ) from None


def structure_interpretation_of(result: Any) -> str:
    """Which atom-name semantics produced this result.

    Reads ``result.source.parser_backend`` -- already recorded on every result
    and already serialized, so nothing new has to be written for a result to
    carry this.
    """
    source = getattr(result, "source", None)
    if source is None:
        raise IncompatibleStructureInterpretation(
            "this result carries no source information, so the structure "
            "interpretation that produced it cannot be established."
        )
    return interpretation_of_backend(getattr(source, "parser_backend", ""))


def _name(item: Any, index: int) -> str:
    for attr in ("sample_id", "structure_id"):
        value = getattr(item, attr, None)
        if value is not None:
            return getattr(value, "value", value)
    return "item %d" % index


def require_common_structure_interpretation(
    items: Sequence[Any],
    *,
    context: Optional[str] = None,
) -> str:
    """Refuse to combine results read under different atom-name semantics.

    Accepts results or :class:`~biotsr.analysis.samples.AnalysisSample`\\ s --
    anything that can state its interpretation, either through ``source`` or
    through an ``annotations[STRUCTURE_INTERPRETATION]`` entry.

    Returns the one interpretation they share, so a caller can record it.
    """
    if not items:
        raise ValueError(
            "no results: a structure interpretation cannot be established"
        )
    found = []
    for index, item in enumerate(items):
        annotations = getattr(item, "annotations", None) or {}
        stated = annotations.get(STRUCTURE_INTERPRETATION)
        found.append(
            (item, index, stated if stated else structure_interpretation_of(item))
        )
    first_item, first_index, first = found[0]
    for item, index, other in found[1:]:
        if other != first:
            raise IncompatibleStructureInterpretation(
                "%sthese results were read under different structure "
                "interpretations and their keys are NOT comparable: %s was read "
                "under %r and %s under %r. Two interpretations disagree about "
                "what an atom is CALLED -- the audited case labels a manganese "
                "ion nitrogen -- so the same physical atom receives a different "
                "label and every triangle it touches receives a different key. "
                "Their KeySpace and SamplingProfile fingerprints MATCH, because "
                "the interpretation is deliberately in neither identity "
                "(ADR-079); it is the atom labelling underneath them that "
                "differs. Re-read one of them under the other's interpretation."
                % (
                    "%s: " % context if context else "",
                    _name(first_item, first_index),
                    first,
                    _name(item, index),
                    other,
                )
            )
    return first
