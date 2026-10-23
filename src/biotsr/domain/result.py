"""The structured result model.

Text files are an *export format*, never the source of truth.
Everything a consumer needs -- method, key space, sampling profile, provenance,
skip reasons, per-category frequencies, raw and folded angles -- rides on the
object, so nothing is ever recovered by parsing a filename.

Two structural requirements are built in rather than bolted on:

* ``CategoryCounts`` -- a key's frequency is always broken out by category.
  Methods with no category concept (ordinary protein Ca TSR, for example)
  use one implicit ``"regular"`` bucket, so this costs nothing today while
  leaving mirror ("regular"/"mirror") and SSE (up to 18 type buckets) expressible
  without a schema change.  Each key owns its own counter -- there is no shared
  mutable state between keys, which is the direct structural fix for the legacy
  shared-list aliasing that corrupts 40% of the historical SSE golden output.
* ``theta_raw`` -- the pre-fold angle is a first-class field on every triplet,
  closing a confirmed gap: no legacy triplet format stores it, forcing every
  consumer to reconstruct it from coordinates.

SSE adds ``TripletSSE`` and nothing else to this model.  SSE is an
annotation: it rides on the triplet and it names the bucket in the counter that
was already there.  It adds no field to ``KeySpace``, no field to
``SamplingProfile``, and no second result type -- which is what makes "the key is
unchanged by annotation" a structural property rather than a test result.
"""

from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence, Tuple

from .identity import AtomIdentity, StructureId
from .keyspace import KeySpace, SamplingProfile

__all__ = [
    "CategoryCounts",
    "TripletSSE",
    "DEFAULT_SIDE",
    "TripletRecord",
    "SkipRecord",
    "SourceInfo",
    "SelectionSummary",
    "ResultCounts",
    "TSRResult",
    "DEFAULT_CATEGORY",
    "RESULT_SCHEMA_VERSION",
    "SUPPORTED_RESULT_SCHEMA_VERSIONS",
    "PROVENANCE_KEYS",
    "LOADED_FROM",
]

#: The annotation key a loaded result carries, naming where it came from and
#: which parts of the original object were deliberately not reconstructed.
#: Defined here, beside the fields it talks about, and re-exported by
#: ``biotsr.io.native_read`` so that a reader and the domain cannot disagree about
#: the spelling.
LOADED_FROM = "loaded_from"

#: The version of the **native serialized result format**, written into every
#: metadata sidecar as ``result_schema_version``.
#:
#: Deliberately not ``biotsr.__version__``.  The package version moves for reasons
#: that have nothing to do with the file format -- a new method, a bug fix, a
#: documentation release -- and tying the two would either force a format
#: migration on every release or make the field useless as a compatibility
#: signal.  This integer changes if, and only if, what a reader must do to a
#: saved result changes.
#:
#: An integer rather than a string, because a reader must be able to say "that is
#: newer than anything I know" and act on it (ADR-074).
#:
#: 2 (1.0.0rc4): the selection record states each side's scope, instances and
#: per-instance atom subsets; every method records the one alternate-location
#: rule; an omitted molecule identifier means every copy on the chain; files
#: are named by a prefix (``PREFIX.keys``, ``PREFIX.meta.json``).  Version 1
#: was written by the 1.0 release candidates up to 1.0.0rc3 and is not read.
RESULT_SCHEMA_VERSION = 2

#: Every schema version this build can read.  A sidecar declaring a version
#: outside this set is refused -- above it as unsupported (upgrade biotsr), below it
#: as unknown.  A reader never guesses.
SUPPORTED_RESULT_SCHEMA_VERSIONS = (2,)

#: The annotation keys that constitute a result's *provenance block* in the
#: sidecar: how the run was asked for, as distinct from what it produced.
#:
#: Without them a saved result could not say which profile made it.  They are
#: serialized as one nested object rather than as loose top-level fields so that
#: "what this result IS" and "how it was REQUESTED" stay separable by a reader.
PROVENANCE_KEYS = (
    "profile",
    "base_profile",
    # Always false in BioTSR 1.0, which implements standard profiles only; kept
    # so the result schema (and every saved result) is unchanged.
    "legacy_compatibility_profile",
    "modifiers",
    "target",
    "truncated_at_line",
    # The structured selection record: per side, the chain, the scope, the
    # arguments and the instances selected (domain/selection_record.py).
    # Provenance only; it enters neither fingerprint and never decides
    # whether two results may be compared.
    "selection_request",
)

#: The extra annotation keys a cross result carries, gathered under
#: ``provenance.cross``.
CROSS_PROVENANCE_KEYS = ("cross_pools", "side_a", "side_b")

#: The single implicit bucket used by methods that have no category concept.
DEFAULT_CATEGORY = "regular"

#: The single implicit pool a one-sided method draws from. Matches
#: ``LabeledPoint.side``'s default, so a non-cross triplet reports it rather than
#: reporting nothing.
DEFAULT_SIDE = "primary"


@dataclass(frozen=True)
class CategoryCounts:
    """Per-key frequency broken out by category."""

    by_category: Mapping[str, int]

    @property
    def total(self) -> int:
        return sum(self.by_category.values())

    def to_dict(self) -> Mapping[str, Any]:
        return {"by_category": dict(self.by_category), "total": self.total}


@dataclass(frozen=True)
class TripletSSE:
    """One triangle's secondary-structure annotation.

    ``states`` and ``groups`` are in the same canonical vertex order as
    ``TripletRecord.labels``, so a vertex's residue, its label and its
    secondary-structure state are read off the same index.  ``category`` is the
    one of the 18 they classify to, and is exactly the string used as this
    triangle's key-frequency bucket.

    Both the per-vertex states and the category are kept.  The legacy triplet
    file stores only the category, which makes every historical row impossible to
    audit: there is no way to tell a misassigned residue from a misclassified
    triangle.  Storing both costs three strings and makes the classification
    checkable against its own inputs.
    """

    states: Tuple[str, str, str]
    groups: Tuple[str, str, str]
    category: str

    @property
    def tokens(self) -> Tuple[str, str, str]:
        """The legacy ``STATE_GROUP`` spellings, e.g. ``("SHEET_2", ...)``."""
        return tuple(
            "%s_%s" % (state, group) for state, group in zip(self.states, self.groups)
        )

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "states": list(self.states),
            "groups": list(self.groups),
            "category": self.category,
        }


@dataclass(frozen=True)
class TripletRecord:
    """One triangle, fully described.

    ``identities`` carries the *physical* identity of each vertex; ``labels``
    carries the integer code the key was built from; ``label_names`` carries the
    label scheme's own name for that code.

    Those last two differ under AA grouping, and the distinction matters: an ILE
    residue grouped with LEU has ``residue_name == "ILE"`` but label name
    ``"LEU"``, because LEU is the representative the grouped table lists first for
    code 8.  The legacy triplets format prints only the label name, which silently
    misreports which residue is present; the native format prints both, so the
    grouped output stays interpretable without losing the physical truth.
    """

    key: int
    labels: Tuple[int, int, int]
    identities: Tuple[AtomIdentity, AtomIdentity, AtomIdentity]
    label_names: Tuple[str, str, str]
    theta_raw: float
    theta: float
    theta_class: int
    distance_class: int
    max_distance: float
    coordinates: Tuple[
        Tuple[float, float, float],
        Tuple[float, float, float],
        Tuple[float, float, float],
    ]
    category: str = DEFAULT_CATEGORY
    #: Which pool each vertex came from, in the same canonical order as
    #: ``labels`` and ``identities``.
    #:
    #: For a one-sided method every vertex is ``"primary"``, the single default
    #: pool. For a cross method the three entries are the side tags, and a cross
    #: triangle always shows both of them -- which is what makes side membership
    #: readable off a triplet record without re-deriving it from the chain, and
    #: what a mixed-side composition assertion is checked against.
    sides: Tuple[str, str, str] = (DEFAULT_SIDE, DEFAULT_SIDE, DEFAULT_SIDE)
    #: The secondary-structure annotation, or ``None`` when SSE was not requested.
    #:
    #: When it is present, ``sse.category == category``: the SSE category *is*
    #: the frequency bucket, exactly as mirror's "regular"/"mirror" is when SSE
    #: is off. The mirror state is not lost when both are active -- it is carried
    #: by the key's sign, and ``theta_raw`` on this record is the value that
    #: decided it.
    sse: Optional["TripletSSE"] = None
    #: ``s3``, the base-to-apex median the angle was computed from. Carried
    #: because the SSE triplet format is the one historical format that prints it
    #: (its ``height`` column); no key depends on it.
    median_length: Optional[float] = None


@dataclass(frozen=True)
class SkipRecord:
    """A candidate, point or entity that did not contribute, and why.

    Skips are data.  Nothing is dropped silently -- the single most-criticised
    property of the legacy filters.
    """

    stage: str
    reason: str
    detail: str = ""
    indices: Optional[Tuple[int, ...]] = None

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "stage": self.stage,
            "reason": self.reason,
            "detail": self.detail,
            "indices": list(self.indices) if self.indices is not None else None,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "SkipRecord":
        indices = payload.get("indices")
        return cls(
            stage=payload["stage"],
            reason=payload["reason"],
            detail=payload.get("detail", ""),
            indices=tuple(indices) if indices is not None else None,
        )


@dataclass(frozen=True)
class SourceInfo:
    """Where the structure came from and how it was read."""

    structure_id: StructureId
    path: Optional[str]
    sha256: Optional[str]
    parser_backend: str
    model: int
    altloc_policy: str
    #: Which modified amino acids written as ``HETATM`` a protein selection
    #: read (``selection.base.ModifiedResiduePolicy``).  ``None`` only for a
    #: result written before the rule was recorded (ADR-105): *not recorded*,
    #: which is never read as the current rule.
    modified_residue_policy: Optional[str] = None
    #: How the structure's records were read and checked
    #: (``selection.base.InputRules``).  ``None`` for a result written by
    #: 1.0.0rc5 or earlier: *not recorded* (ADR-108).
    input_rules: Optional[str] = None

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "structure_id": self.structure_id.value,
            "path": self.path,
            "sha256": self.sha256,
            "parser_backend": self.parser_backend,
            "model": self.model,
            "altloc_policy": self.altloc_policy,
            "modified_residue_policy": self.modified_residue_policy,
            "input_rules": self.input_rules,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "SourceInfo":
        return cls(
            structure_id=StructureId(payload["structure_id"]),
            path=payload.get("path"),
            sha256=payload.get("sha256"),
            parser_backend=payload["parser_backend"],
            model=payload["model"],
            altloc_policy=payload["altloc_policy"],
            modified_residue_policy=payload.get("modified_residue_policy"),
            input_rules=payload.get("input_rules"),
        )


@dataclass(frozen=True)
class SelectionSummary:
    """What selection actually used -- computed, never echoed from the request."""

    selector: str
    chains_requested: Tuple[str, ...]
    chains_used: Tuple[str, ...]
    n_points: int
    residues_used: Tuple[str, ...] = ()
    atoms_missing: Tuple[str, ...] = ()

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "selector": self.selector,
            "chains_requested": list(self.chains_requested),
            "chains_used": list(self.chains_used),
            "n_points": self.n_points,
            "n_residues_used": len(self.residues_used),
            # The tuple itself, alongside -- not instead of -- its length.
            # ``n_residues_used`` predates this and every existing consumer reads
            # it; ``residues_used`` is what makes the summary reconstructible.
            "residues_used": list(self.residues_used),
            "atoms_missing": list(self.atoms_missing),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "SelectionSummary":
        """Rebuild from :meth:`to_dict`.

        ``n_residues_used`` is derived and is deliberately *not* read back: it is
        ``len(residues_used)`` by construction, and reading both would let a
        hand-edited sidecar assert a count its own list contradicts.
        """
        return cls(
            selector=payload["selector"],
            chains_requested=tuple(payload.get("chains_requested") or ()),
            chains_used=tuple(payload.get("chains_used") or ()),
            n_points=payload["n_points"],
            residues_used=tuple(payload.get("residues_used") or ()),
            atoms_missing=tuple(payload.get("atoms_missing") or ()),
        )


@dataclass(frozen=True)
class ResultCounts:
    """The candidate ledger.

    The four ways a candidate can end up, kept as separate counters because they
    mean different things and legacy conflated them by counting none of them:

    * ``n_degenerate``     -- no definable angle (geometry invalid);
    * ``n_filtered_by_size`` -- perfectly valid, but outside the caller's cutoff;
    * ``n_triangles``      -- encoded into a key;
    * (the rest is nothing: every candidate lands in exactly one of the above.)

    The ledger balances exactly, and a test asserts it:
    ``n_candidates == n_degenerate + n_filtered_by_size + n_triangles``.
    """

    n_points: int
    n_candidates: int
    n_triangles: int
    n_unique_keys: int
    n_valid_geometry: int = 0
    n_degenerate: int = 0
    n_acos_clamped: int = 0
    n_filtered_by_size: int = 0

    @property
    def ledger_balances(self) -> bool:
        return (
            self.n_candidates
            == self.n_degenerate + self.n_filtered_by_size + self.n_triangles
        )

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "n_points": self.n_points,
            "n_candidates": self.n_candidates,
            "n_valid_geometry": self.n_valid_geometry,
            "n_triangles": self.n_triangles,
            "n_unique_keys": self.n_unique_keys,
            "n_degenerate": self.n_degenerate,
            "n_acos_clamped": self.n_acos_clamped,
            "n_filtered_by_size": self.n_filtered_by_size,
            "ledger_balances": self.ledger_balances,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ResultCounts":
        """Rebuild from :meth:`to_dict`.

        ``ledger_balances`` is derived and is not read back, for the same reason
        ``n_residues_used`` is not: it is a property of the other six numbers, and
        accepting it as input would let a file assert that it balances when it
        does not.  The reader recomputes it and the loader checks it.
        """
        return cls(
            n_points=payload["n_points"],
            n_candidates=payload["n_candidates"],
            n_triangles=payload["n_triangles"],
            n_unique_keys=payload["n_unique_keys"],
            n_valid_geometry=payload.get("n_valid_geometry", 0),
            n_degenerate=payload.get("n_degenerate", 0),
            n_acos_clamped=payload.get("n_acos_clamped", 0),
            n_filtered_by_size=payload.get("n_filtered_by_size", 0),
        )


@dataclass(frozen=True)
class TSRResult:
    """The structured result of one TSR calculation."""

    method: str
    structure_id: StructureId
    chain: str
    keys: Mapping[int, CategoryCounts]
    key_order: Tuple[int, ...]
    triplets: Optional[Sequence[TripletRecord]]
    key_space: KeySpace
    sampling: SamplingProfile
    selection: SelectionSummary
    source: SourceInfo
    counts: ResultCounts
    skipped: Sequence[SkipRecord] = ()
    warnings: Sequence[str] = ()
    annotations: Mapping[str, Any] = field(default_factory=dict)

    @property
    def frequencies(self) -> Mapping[int, int]:
        """Flat ``key -> total frequency`` view, for the common single-category case."""
        return {k: v.total for k, v in self.keys.items()}

    @property
    def n_triangles(self) -> int:
        return self.counts.n_triangles

    @property
    def n_unique_keys(self) -> int:
        return self.counts.n_unique_keys

    def category_totals(self) -> Mapping[str, int]:
        totals: "Counter[str]" = Counter()
        for counts in self.keys.values():
            for category, n in counts.by_category.items():
                totals[category] += n
        return dict(totals)

    def provenance(self) -> Mapping[str, Any]:
        """How this run was *asked for*, as distinct from what it produced.

        The profile that ran, the profile it was derived from, whether that
        profile reproduces a historical program, which modifiers were active, the
        entity a per-entity run targeted, and -- for a cross result -- the two
        sides and the contact pool.

        Without this block a saved result could not say which convention
        produced it, and that must never be inferred from a filename.
        """
        block = {}
        for key in PROVENANCE_KEYS:
            value = self.annotations.get(key)
            block[key] = list(value) if isinstance(value, tuple) else value
        cross = {
            key: copy.deepcopy(self.annotations[key])
            for key in CROSS_PROVENANCE_KEYS
            if key in self.annotations
        }
        block["cross"] = cross or None
        return block

    def metadata(self) -> Mapping[str, Any]:
        """The self-describing sidecar payload.

        Every nested block is a **copy**.  ``metadata()`` used to embed
        ``annotations["sse"]`` by reference, so mutating the returned payload
        mutated the result's own SSE provenance -- harmless while nothing read a
        sidecar back, and a live corruption hazard now that something does.
        """
        return {
            "result_schema_version": RESULT_SCHEMA_VERSION,
            "method": self.method,
            "structure_id": self.structure_id.value,
            "chain": self.chain,
            "key_space": dict(self.key_space.to_dict()),
            "key_space_fingerprint": self.key_space.fingerprint,
            "sampling": dict(self.sampling.to_dict()),
            "sampling_fingerprint": self.sampling.fingerprint,
            "selection": dict(self.selection.to_dict()),
            "source": dict(self.source.to_dict()),
            "counts": dict(self.counts.to_dict()),
            "category_totals": dict(self.category_totals()),
            "provenance": dict(self.provenance()),
            "sse": copy.deepcopy(self.annotations.get("sse")),
            "sse_coverage": copy.deepcopy(self.annotations.get("sse_coverage")),
            "skipped": [dict(s.to_dict()) for s in self.skipped],
            "warnings": list(self.warnings),
        }
