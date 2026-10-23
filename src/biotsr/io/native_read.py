"""Reading a native TSR result back.  The inverse of ``biotsr.io.native``.

Without this module the native format would be **write-only**: a result could be
saved and never loaded, and every analysis of a previously computed structure
would begin by recomputing it.  This module closes that, and the one thing it
must never do is guess.

Three principles, each of which shows up as a refusal rather than a comment:

1. **The schema version decides what is readable.**  A sidecar states
   ``result_schema_version``; a reader that does not know that version says so
   and stops.  It never infers a shape from the fields that happen to be present,
   because a field's *absence* is exactly what distinguishes two schemas and
   guessing would silently read one format as another.
2. **Nothing is defaulted into existence.**  Every field of ``KeySpace`` and
   ``SamplingProfile`` is identity-bearing, and a fabricated value would produce
   a fingerprint -- a plausible, comparable, wrong fingerprint.  Missing means
   refused.
3. **The files must belong together.**  A saved result is two or three files
   related by a filename convention, and a filename convention is not a check.
   Five quantities the sidecar already carries -- both fingerprints,
   ``n_unique_keys``, ``n_triangles`` and ``category_totals`` -- are recomputed
   from the keys file and compared.  Pairing ``sample_A``'s keys with
   ``sample_B``'s sidecar fails on the first of them that differs, naming both
   values.

**The fidelity boundary is explicit** (ADR-075).  A loaded result carries every
key, every per-category frequency, both scientific identities with their exact
fingerprints, the profile provenance, the SSE annotation provenance, the
selection summary, the source provenance, the candidate ledger, the skip records
and the warnings.  It does **not** carry:

``triplets``
    The triplet file writes six of ``AtomIdentity``'s eight fields.  ``model`` is
    recoverable from the sidecar, but ``alt_loc`` is genuinely per-atom -- under
    ``blank_or_a`` one structure contributes both ``""`` and ``"A"`` atoms -- and
    is in no column.  Reconstructing it would be a guess, so triplets are not
    read and ``loaded.triplets is None``.
``key_order``
    First-seen key order, which no BioTSR output uses.
    Serializing it costs 408 KB against a 3.8 KB schema-1 sidecar on 1UBQ chain
    A -- 106x -- to enable re-writing historical files from a saved result, which
    is an operation nobody performs.  ``loaded.key_order == ()``, and the legacy
    writers refuse a loaded result **by name** rather than writing an empty file.

Both are capability distinctions, not silent losses: a loaded result says so in
``annotations["loaded_from"]``, and the operations that need the missing state
name what is missing and why.
"""

from __future__ import annotations

import json
import os
import re
from collections import Counter
from typing import Any, Dict, Mapping, Optional, Tuple

from ..domain.errors import MalformedResult, UnsupportedFormat
from ..domain.identity import StructureId
from ..domain.keyspace import KeySpace, SamplingProfile
from ..domain.result import (
    SUPPORTED_RESULT_SCHEMA_VERSIONS,
    CategoryCounts,
    ResultCounts,
    SelectionSummary,
    SkipRecord,
    SourceInfo,
    TSRResult,
)
from .native import ENCODING, KEY_COLUMNS, KEYS_SUFFIX, METADATA_SUFFIX

__all__ = [
    "METADATA_SUFFIX",
    "CANDIDATE_LAYOUT_MESSAGE",
    "LOADED_FROM",
    "read_metadata",
    "read_keys",
    "load_result",
    "sidecar_path_for",
]

#: Names written by the release candidates up to 1.0.0rc3: the bin counts in the
#: extension (``1UBQ_A.keys_29_35``) and the sidecar after the keys name.
_CANDIDATE_LAYOUT = re.compile(r"\.(keys|triplets)_\d+_\d+(\.meta\.json)?$")

#: Why such a file is not read, said once.
CANDIDATE_LAYOUT_MESSAGE = (
    "%s was written by a BioTSR 1.0 release candidate (1.0.0rc3 or earlier), "
    "in the per-sample directory layout with the bin counts in the file name. "
    "This version reads results written in its own layout only (PREFIX.keys, "
    "PREFIX.meta.json, result schema 2): the selection record, the "
    "alternate-location rule of amino_acid_atomic and the meaning of an omitted "
    "molecule identifier all changed, so an older result cannot be read as if "
    "it were current. Recompute it with this version; the same request "
    "reproduces every key that did not depend on those changes."
)

#: The annotation key under which a loaded result records that it was loaded,
#: where from, and what it therefore cannot do.  Present on loaded results and
#: absent on computed ones, so ``"loaded_from" in result.annotations`` is the
#: reliable test -- not a guess from a missing ``triplets``, which a
#: ``keep_triplets=False`` run also has.
LOADED_FROM = "loaded_from"


def sidecar_path_for(keys_path: str) -> str:
    """The sidecar that belongs to a keys file: ``PREFIX.keys`` -> ``PREFIX.meta.json``."""
    if not keys_path.endswith(KEYS_SUFFIX):
        raise ValueError("%s is not a native keys file (PREFIX%s)"
                         % (keys_path, KEYS_SUFFIX))
    return keys_path[: -len(KEYS_SUFFIX)] + METADATA_SUFFIX


class _Reconstructing(object):
    """Translate anything a malformed payload throws into ``MalformedResult``.

    The sub-record ``from_dict`` methods are written for clarity -- they index
    the fields they need -- so a truncated ``counts`` block raises ``KeyError``
    and a ``structure_id`` that is not a legal identifier raises ``ValueError``
    from ``StructureId``. Both are *correct diagnoses of a corrupt file* arriving
    under the wrong name: a caller who writes ``except MalformedResult`` would
    miss them, and the CLI would report a `KeyError` as though the request had
    been malformed rather than the file.

    Putting the translation at the boundary rather than inside each ``from_dict``
    keeps those methods usable on their own -- tests and callers construct them
    from hand-built dicts -- and keeps the contract stated in one place: **what
    this module raises for a bad file is** ``MalformedResult``.

    ``MalformedResult`` and ``UnsupportedFormat`` pass through unchanged; they
    already say what went wrong, in better words than a wrapper could.
    """

    def __init__(self, what: str, path: str):
        self.what = what
        self.path = path

    def __enter__(self) -> "_Reconstructing":
        return self

    def __exit__(self, exc_type, exc, _tb) -> bool:
        if exc is None or isinstance(exc, (MalformedResult, UnsupportedFormat)):
            return False
        if isinstance(exc, (KeyError, TypeError, ValueError, IndexError)):
            detail = exc.args[0] if isinstance(exc, KeyError) and exc.args else exc
            raise MalformedResult(
                "%s: cannot reconstruct the %s -- %s: %s"
                % (self.path, self.what, type(exc).__name__, detail)
            ) from None
        return False


# ---------------------------------------------------------------------------
# the sidecar
# ---------------------------------------------------------------------------


def _check_schema_version(payload: Mapping[str, Any], where: str) -> int:
    """Decide whether this build can read this sidecar at all.

    Four outcomes, and none of them is "read it anyway".
    """
    if "result_schema_version" not in payload:
        raise UnsupportedFormat(
            "%s declares no 'result_schema_version'. Sidecars written by "
            "pre-release development builds carry no version field, and this "
            "reader will not guess "
            "which shape an unversioned file has -- guessing is how one format "
            "gets silently read as another. Re-run the calculation to write a "
            "versioned result; nothing scientific is lost, because the keys are "
            "reproduced exactly by the same request." % where
        )
    version = payload["result_schema_version"]
    if isinstance(version, bool) or not isinstance(version, int):
        raise MalformedResult(
            "%s has a 'result_schema_version' of %r (%s); it must be an integer"
            % (where, version, type(version).__name__)
        )
    if version in SUPPORTED_RESULT_SCHEMA_VERSIONS:
        return version
    highest = max(SUPPORTED_RESULT_SCHEMA_VERSIONS)
    if version > highest:
        raise UnsupportedFormat(
            "%s was written by a newer biotsr: it declares result schema version %d "
            "and this build reads up to %d. Upgrade biotsr to read it. It is NOT "
            "read partially -- a schema version exists precisely so that a "
            "reader which does not understand a format refuses it instead of "
            "misinterpreting it." % (where, version, highest)
        )
    if version == 1:
        raise UnsupportedFormat(CANDIDATE_LAYOUT_MESSAGE % where)
    raise UnsupportedFormat(
        "%s declares result schema version %d, which this build does not "
        "support (supported: %s)"
        % (where, version, ", ".join(str(v) for v in SUPPORTED_RESULT_SCHEMA_VERSIONS))
    )


def read_metadata(path: str) -> Mapping[str, Any]:
    """Read and version-check one ``.meta.json`` sidecar.

    Returns the raw payload.  Reconstruction happens in :func:`load_result`,
    which also needs the keys file; this function exists on its own so that a
    caller can ask "what is this file?" without paying for the keys.
    """
    try:
        with open(path, "r", encoding=ENCODING) as fh:
            payload = json.load(fh)
    except ValueError as exc:
        raise MalformedResult("%s is not readable as JSON: %s" % (path, exc))
    if not isinstance(payload, dict):
        raise MalformedResult(
            "%s must contain a JSON object, got %s" % (path, type(payload).__name__)
        )
    _check_schema_version(payload, path)
    return payload


# ---------------------------------------------------------------------------
# the keys file
# ---------------------------------------------------------------------------


def read_keys(path: str) -> Tuple[Mapping[int, CategoryCounts], Mapping[str, int]]:
    """Read a native keys file into ``{key: CategoryCounts}`` and category totals.

    The file's own ordering guarantee is **verified, not assumed**.  The writer
    promises numerically ascending keys, then alphabetical categories; a file in
    any other order was not produced by this writer, and reading it as if it had
    been would mean accepting hand-edited or concatenated input as authentic.
    """
    keys: Dict[int, Dict[str, int]] = {}
    totals: "Counter[str]" = Counter()
    previous: Optional[Tuple[int, str]] = None

    with open(path, "r", encoding=ENCODING) as fh:
        header = fh.readline()
        if not header:
            raise MalformedResult("%s is empty; a keys file has a header" % path)
        columns = tuple(header.rstrip("\n").split("\t"))
        if columns != KEY_COLUMNS:
            raise MalformedResult(
                "%s has header %r; a native keys file's header is %r"
                % (path, "\t".join(columns), "\t".join(KEY_COLUMNS))
            )
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(KEY_COLUMNS):
                raise MalformedResult(
                    "%s line %d has %d field(s), expected %d: %r"
                    % (path, lineno, len(fields), len(KEY_COLUMNS), line.rstrip("\n"))
                )
            raw_key, category, raw_frequency = fields
            try:
                key = int(raw_key)
            except ValueError:
                raise MalformedResult(
                    "%s line %d: %r is not an integer key" % (path, lineno, raw_key)
                )
            if not category:
                raise MalformedResult(
                    "%s line %d: the category is empty. A native keys file always "
                    "states the category explicitly -- %r for a run with no "
                    "category concept -- because the key's sign alone cannot "
                    "carry it (key 0 negates to itself)."
                    % (path, lineno, "regular")
                )
            try:
                frequency = int(raw_frequency)
            except ValueError:
                raise MalformedResult(
                    "%s line %d: %r is not an integer frequency"
                    % (path, lineno, raw_frequency)
                )
            if frequency < 0:
                raise MalformedResult(
                    "%s line %d: frequency %d is negative" % (path, lineno, frequency)
                )
            current = (key, category)
            if previous is not None and current <= previous:
                raise MalformedResult(
                    "%s line %d breaks the format's ordering guarantee: "
                    "(%d, %r) does not follow (%d, %r). A native keys file is "
                    "sorted numerically by key, then alphabetically by category, "
                    "and a file in another order was not written by this package."
                    % (path, lineno, key, category, previous[0], previous[1])
                )
            previous = current
            keys.setdefault(key, {})[category] = frequency
            totals[category] += frequency

    return (
        {key: CategoryCounts(by_category=dict(by)) for key, by in keys.items()},
        dict(totals),
    )


# ---------------------------------------------------------------------------
# consistency between the files
# ---------------------------------------------------------------------------


def _agree(what: str, from_keys, from_sidecar, keys_path: str, meta_path: str) -> None:
    if from_keys != from_sidecar:
        raise MalformedResult(
            "these files do not belong to the same result: %s is %r in\n"
            "  %s\nbut %r in\n  %s\n"
            "A saved native result is two or three files related only by their "
            "names, so the reader checks the quantities they share rather than "
            "trusting the naming. Combining one sample's keys with another's "
            "sidecar is exactly what this refusal exists to catch."
            % (what, from_keys, keys_path, from_sidecar, meta_path)
        )


def _check_consistency(
    payload: Mapping[str, Any],
    keys: Mapping[int, CategoryCounts],
    totals: Mapping[str, int],
    key_space: KeySpace,
    sampling: SamplingProfile,
    counts: ResultCounts,
    keys_path: str,
    meta_path: str,
) -> None:
    """Six checks, in the order that makes a failure most informative.

    One more runs before this -- the schema version (and a release-candidate
    file name is refused before that) -- and three inside :func:`read_keys`,
    which validates the header, each row's shape and the file's sorted order.
    The file name is not checked against the result: the caller chooses the
    prefix, and the fingerprints and counts below are what pair a keys file
    with its sidecar.
    """
    # 1-2. The sidecar must agree with itself: a stated fingerprint that does not
    # recompute from the block beside it means the block was edited.
    for name, recomputed in (
        ("key_space_fingerprint", key_space.fingerprint),
        ("sampling_fingerprint", sampling.fingerprint),
    ):
        stated = payload.get(name)
        if stated is not None and stated != recomputed:
            raise MalformedResult(
                "%s: %s is %r but recomputing it from the %s block in the same "
                "file gives %r. The sidecar contradicts itself -- one of the two "
                "was edited, and the identity a result claims must be derivable "
                "from the semantics it records."
                % (meta_path, name, stated, name.rsplit("_", 1)[0], recomputed)
            )

    # 3-5. The sidecar must agree with the keys file.
    _agree("n_unique_keys", len(keys), counts.n_unique_keys, keys_path, meta_path)
    _agree("n_triangles", sum(totals.values()), counts.n_triangles, keys_path, meta_path)
    stated_totals = payload.get("category_totals")
    if stated_totals is not None:
        # Typed before it is compared: a hand-corrupted sidecar carrying, say, a
        # list here must leave through MalformedResult like everything else. An
        # AttributeError would leave through the CLI's INTERNAL boundary instead,
        # telling the user biotsr is broken when it is their file that is.
        if not isinstance(stated_totals, Mapping):
            raise MalformedResult(
                "%s: 'category_totals' must be a JSON object, got %s"
                % (meta_path, type(stated_totals).__name__)
            )
        _agree(
            "category_totals",
            dict(totals),
            dict(stated_totals),
            keys_path,
            meta_path,
        )

    # 6. The ledger the sidecar records must still balance.
    if not counts.ledger_balances:
        raise MalformedResult(
            "%s: the candidate ledger does not balance -- %d candidates against "
            "%d degenerate + %d size-filtered + %d triangles. Every candidate "
            "lands in exactly one of those three, so a result in which they do "
            "not sum is corrupt."
            % (
                meta_path,
                counts.n_candidates,
                counts.n_degenerate,
                counts.n_filtered_by_size,
                counts.n_triangles,
            )
        )


# ---------------------------------------------------------------------------
# the loader
# ---------------------------------------------------------------------------


def load_result(path: str, *, metadata_path: Optional[str] = None) -> TSRResult:
    """Load a previously saved native result.

    ``path`` is the keys file (``PREFIX.keys``) or the sidecar beside it
    (``PREFIX.meta.json``); either names the pair.  ``metadata_path``
    overrides where the sidecar is looked for, which is what makes the
    mismatched-files refusal testable without renaming anything on disk.

    The returned object is a real :class:`~biotsr.domain.result.TSRResult` and every
    analysis in :mod:`biotsr.analysis` consumes it unchanged, because every one of
    them reads ``frequencies``, ``key_space`` and ``sampling`` and nothing else.
    What it cannot do -- triplets and first-seen key order -- is stated in
    ``annotations[%r]`` and enforced by the operations that would need them.

    A file in the release-candidate layout (``1UBQ_A.keys_29_35``) is refused
    with :class:`~biotsr.domain.errors.UnsupportedFormat`, naming why.
    """
    if _CANDIDATE_LAYOUT.search(os.path.basename(path)):
        raise UnsupportedFormat(CANDIDATE_LAYOUT_MESSAGE % path)
    if path.endswith(METADATA_SUFFIX):
        keys_path = path[: -len(METADATA_SUFFIX)] + KEYS_SUFFIX
        default_meta = path
    elif path.endswith(KEYS_SUFFIX):
        keys_path = path
        default_meta = sidecar_path_for(path)
    else:
        raise ValueError(
            "%s is neither a native keys file (PREFIX%s) nor its sidecar "
            "(PREFIX%s)" % (path, KEYS_SUFFIX, METADATA_SUFFIX)
        )
    meta_path = metadata_path if metadata_path is not None else default_meta

    if not os.path.isfile(keys_path):
        raise FileNotFoundError(
            "no native keys file at %s. A saved native result is a keys file "
            "plus the '%s' sidecar beside it; point load_result at either."
            % (keys_path, METADATA_SUFFIX)
        )
    if not os.path.isfile(meta_path):
        raise FileNotFoundError(
            "no metadata sidecar at %s. The keys file alone cannot be read back: "
            "it carries integers and frequencies, and nothing about what those "
            "integers MEAN. The key space and the sampling profile live in the "
            "sidecar, and without them the numbers are not comparable with "
            "anything." % meta_path
        )

    payload = read_metadata(meta_path)
    key_space = KeySpace.from_dict(payload.get("key_space"))
    sampling = SamplingProfile.from_dict(payload.get("sampling"))

    for required in ("method", "structure_id", "chain", "counts", "source", "selection"):
        if required not in payload:
            raise MalformedResult(
                "%s is missing the required field %r" % (meta_path, required)
            )
    with _Reconstructing("candidate ledger", meta_path):
        counts = ResultCounts.from_dict(payload["counts"])

    keys, totals = read_keys(keys_path)
    _check_consistency(
        payload, keys, totals, key_space, sampling, counts, keys_path, meta_path
    )

    provenance = dict(payload.get("provenance") or {})
    annotations: Dict[str, Any] = {
        "profile": provenance.get("profile"),
        "base_profile": provenance.get("base_profile"),
        "legacy_compatibility_profile": provenance.get(
            "legacy_compatibility_profile", False
        ),
        "modifiers": tuple(provenance.get("modifiers") or ()),
        "target": provenance.get("target"),
        "truncated_at_line": provenance.get("truncated_at_line"),
        "selection_request": provenance.get("selection_request"),
        "sse": payload.get("sse"),
        "sse_coverage": payload.get("sse_coverage"),
        LOADED_FROM: {
            "keys": os.path.abspath(keys_path),
            "metadata": os.path.abspath(meta_path),
            "result_schema_version": payload["result_schema_version"],
            # Stated rather than implied. A caller asking "why is triplets None?"
            # gets the answer from the result itself.
            "triplets_loaded": False,
            "key_order_loaded": False,
        },
    }
    cross = provenance.get("cross")
    if cross:
        annotations.update(cross)

    with _Reconstructing("result", meta_path):
        return _build(payload, keys, key_space, sampling, counts, annotations)


def _build(payload, keys, key_space, sampling, counts, annotations) -> TSRResult:
    return TSRResult(
        method=payload["method"],
        structure_id=StructureId(payload["structure_id"]),
        chain=payload["chain"],
        keys=keys,
        # Empty, never fabricated: see the module docstring. `sorted(keys)` would
        # be a plausible-looking claim about an order nobody recorded.
        key_order=(),
        triplets=None,
        key_space=key_space,
        sampling=sampling,
        selection=SelectionSummary.from_dict(payload["selection"]),
        source=SourceInfo.from_dict(payload["source"]),
        counts=counts,
        skipped=tuple(SkipRecord.from_dict(s) for s in payload.get("skipped") or ()),
        warnings=tuple(payload.get("warnings") or ()),
        annotations=annotations,
    )


load_result.__doc__ = load_result.__doc__ % LOADED_FROM
