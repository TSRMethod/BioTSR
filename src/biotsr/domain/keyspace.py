"""The two scientific identities: ``KeySpace`` and ``SamplingProfile``.

``KeySpace`` answers "can these integers be compared at all?"; ``SamplingProfile``
answers "were these counts produced the same way?"  (ADR-003).  Keeping them
separate is the design's central safety property, so they are separate objects
with separate fingerprints and are never merged into one version string.

Each fingerprint hashes every scheme's name, content fingerprint and radix and
every convention's name, so two schemes that share a name but differ in content
get different fingerprints, and renaming a scheme also changes the fingerprint.
The content fingerprint answers a hazard that is real: the legacy
``TSR-Package`` Ca distance binner and ``WebsiteTSR``'s share every cut point,
both describe themselves as 35-class, and still produce different keys for any
triangle with ``max_distance >= 132.45`` -- because only their overflow behaviour
differs (ADR-014).  A human-assigned version string would have been identical
for both.  ``SchemeRef`` therefore carries the scheme's *content* fingerprint,
which each scheme computes over its boundaries, class count, interval semantics
and overflow behaviour together.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, fields
from typing import Any, Mapping, Optional

from .errors import IncompatibleRepresentation, MalformedResult

__all__ = ["SchemeRef", "KeySpace", "SamplingProfile", "IncompatibleKeySpace", "IncompatibleSampling"]


class IncompatibleKeySpace(IncompatibleRepresentation):
    """Raised when results whose keys mean different things would be combined."""


class IncompatibleSampling(IncompatibleRepresentation):
    """Raised when results whose counts are not comparable would be combined."""


def _require_mapping(payload: Any, what: str) -> Mapping[str, Any]:
    if not isinstance(payload, Mapping):
        raise MalformedResult(
            "%s must be a JSON object, got %s" % (what, type(payload).__name__)
        )
    return payload


def _require_str(payload: Mapping[str, Any], name: str, what: str) -> str:
    """Pull a required identity-bearing string, with no default of any kind.

    Every field this is used for changes what a key MEANS.  Defaulting a missing
    one would manufacture a key space that was never measured, and its
    fingerprint would be a plausible-looking hash of a fiction -- which is worse
    than an error, because it compares equal to something.
    """
    if name not in payload:
        raise MalformedResult(
            "%s is missing the required field %r; every field of a key space is "
            "identity-bearing and none of them has a default" % (what, name)
        )
    value = payload[name]
    if not isinstance(value, str):
        raise MalformedResult(
            "%s field %r must be a string, got %r (%s)"
            % (what, name, value, type(value).__name__)
        )
    return value


def _digest(payload: Mapping[str, Any]) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SchemeRef:
    """A reference to a bin or label scheme: readable name + content fingerprint.

    ``radix`` is the value the key encoder actually multiplies by.  It is stated
    explicitly and is *not* inferred from a table's size -- the grouped residue
    scheme has 23 distinct codes but radix 30, because its codes span 4..33
    (ADR-107), and the legacy encoder used radix 20 for 13 grouped codes
    (ADR-015); this field is what makes such facts visible.
    """

    name: str
    content_fingerprint: str
    radix: int

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "name": self.name,
            "content_fingerprint": self.content_fingerprint,
            "radix": self.radix,
        }

    @classmethod
    def from_dict(cls, payload: Any, what: str = "scheme") -> "SchemeRef":
        """Rebuild a reference from :meth:`to_dict`.

        Note what is *not* rebuilt: the scheme itself.  A ``SchemeRef`` is a name,
        a content fingerprint and a radix -- three scalars -- and that is the
        whole of what enters ``KeySpace``.  There is no ``LadderBinScheme`` or
        ``LabelScheme`` to reconstruct here, and that is precisely why a key
        space round-trips exactly: the serialized form is not a lossy projection
        of a live object, it IS the object's contribution to the identity.
        """
        payload = _require_mapping(payload, what)
        radix = payload.get("radix")
        if not isinstance(radix, int) or isinstance(radix, bool):
            raise MalformedResult(
                "%s field 'radix' must be an integer, got %r (%s)"
                % (what, radix, type(radix).__name__)
            )
        if radix < 1:
            raise MalformedResult("%s radix must be positive, got %d" % (what, radix))
        return cls(
            name=_require_str(payload, "name", what),
            content_fingerprint=_require_str(payload, "content_fingerprint", what),
            radix=radix,
        )


@dataclass(frozen=True)
class KeySpace:
    """Everything that changes what an integer key *means*.

    Every field here is key-defining.  Nothing that merely changes which
    triangles were drawn belongs in this object -- that is ``SamplingProfile``.
    """

    angle_convention: str
    #: How a squared distance becomes a distance. Key-affecting, because the two
    #: legacy spellings (``** 0.5`` and ``math.sqrt``) disagree by one ULP often
    #: enough to move a triangle across a bin boundary -- see
    #: ``geometry/conventions.py::SqrtConvention``.
    sqrt_convention: str
    #: How the apex-to-base-midpoint median is computed. Key-affecting for the
    #: same reason: ``TSR-Package``'s two cross programs derive it from
    #: Apollonius' formula over already-rounded edge lengths rather than from the
    #: midpoint coordinate, and it is a denominator in the theta ratio
    #: (``geometry/conventions.py::MedianConvention``, ADR-039).
    median_convention: str
    #: The floating-point precision the coordinates were held at when the key was
    #: computed: ``"float64"`` (the standard) or ``"float32"``. This is
    #: key-affecting and is **not** mere provenance. Measured, not assumed: on
    #: 1BRS A/D, 2 of 1,230,714 triangles land on opposite sides of a theta-bin
    #: boundary under the two precisions (37.339998 vs the 37.34 cut, 52.430001
    #: vs 52.43), because ``Bio.PDB`` stores ``Atom.coord`` as float32 and the
    #: two ``TSR-Package`` cross programs read through it. Two runs that differ
    #: only here can disagree on a key, so they do not share a key space
    #: (ADR-053).
    numerical_precision: str
    theta_scheme: SchemeRef
    distance_scheme: SchemeRef
    label_scheme: SchemeRef
    canonical_ordering_version: str
    key_formula_version: str
    mirror_convention: str

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "angle_convention": self.angle_convention,
            "sqrt_convention": self.sqrt_convention,
            "median_convention": self.median_convention,
            "numerical_precision": self.numerical_precision,
            "theta_scheme": dict(self.theta_scheme.to_dict()),
            "distance_scheme": dict(self.distance_scheme.to_dict()),
            "label_scheme": dict(self.label_scheme.to_dict()),
            "canonical_ordering_version": self.canonical_ordering_version,
            "key_formula_version": self.key_formula_version,
            "mirror_convention": self.mirror_convention,
        }

    #: The scalar (non-``SchemeRef``) fields, in declaration order.  Derived from
    #: the dataclass rather than restated, so a field added to ``KeySpace``
    #: without a corresponding reader change fails loudly at load time instead of
    #: being silently dropped.
    _SCHEME_FIELDS = ("theta_scheme", "distance_scheme", "label_scheme")

    @classmethod
    def from_dict(cls, payload: Any) -> "KeySpace":
        """Rebuild a key space from :meth:`to_dict`.  The exact inverse.

        Three properties this guarantees, each gated:

        * **every identity-bearing field is required.**  There is no default for
          any of them, because a defaulted field would produce a key space that
          describes a calculation nobody ran, carrying a fingerprint that
          compares equal to other fictions;
        * **unknown fields are refused.**  A payload carrying a field this
          version does not know is a payload from a different schema, and
          silently ignoring it would mean reading a key space while discarding
          part of what defines it;
        * ``fingerprint(from_dict(ks.to_dict())) == fingerprint(ks)``, bit for
          bit, because ``to_dict`` is what the fingerprint is computed over.
        """
        payload = _require_mapping(payload, "key_space")
        names = tuple(f.name for f in fields(cls))
        unknown = sorted(set(payload) - set(names))
        if unknown:
            raise MalformedResult(
                "key_space carries %d field(s) this version of biotsr does not "
                "know: %s. A key space is only meaningful if every field that "
                "defines it is understood, so this is refused rather than read "
                "partially." % (len(unknown), ", ".join(repr(u) for u in unknown))
            )
        values = {}
        for name in names:
            if name in cls._SCHEME_FIELDS:
                if name not in payload:
                    raise MalformedResult(
                        "key_space is missing the required field %r" % name
                    )
                values[name] = SchemeRef.from_dict(
                    payload[name], "key_space.%s" % name
                )
            else:
                values[name] = _require_str(payload, name, "key_space")
        return cls(**values)

    @property
    def fingerprint(self) -> str:
        """Stable hash over every key-defining semantic, canonically serialized."""
        return _digest(dict(self.to_dict()))

    def require_compatible(self, other: "KeySpace") -> None:
        if self.fingerprint != other.fingerprint:
            diffs = []
            a, b = dict(self.to_dict()), dict(other.to_dict())
            for field_name in sorted(a):
                if a[field_name] != b[field_name]:
                    diffs.append("%s (%r vs %r)" % (field_name, a[field_name], b[field_name]))
            raise IncompatibleKeySpace(
                "results do not share a key space; their integer keys are not "
                "comparable. Differences: " + "; ".join(diffs)
            )


@dataclass(frozen=True)
class SamplingProfile:
    """Everything that changes *which* triangles were drawn.

    None of these change what a key means, so none of them belong in
    ``KeySpace``.  Equally, nothing here is a runtime setting: worker count,
    output directory, output format, cache location and log level cannot change
    a key or a count and appear in neither identity (ADR-003).
    """

    selector: str
    scope: str
    candidate_strategy: str
    completeness_policy: Optional[str] = None
    composition_rule: Optional[str] = None
    size_filter: Optional[float] = None
    degenerate_policy: str = "skip_and_record"
    #: Which residues a Ca selector admitted by the shape of their sequence
    #: label; ``"accept_all"`` in every BioTSR calculation. A rule here would
    #: change WHICH triangles exist, never what a key means, so it belongs here
    #: and not in ``KeySpace``.
    seq_label_policy: str = "accept_all"
    # -- cross-only sampling facts -------------------------------------------
    # All of these change WHICH points and triangles were drawn and none of them
    # changes what an integer key means, so this is where they belong. In
    # particular, changing `contact_cutoff` must leave `KeySpace` untouched --
    # a test asserts it -- because a 3.0 A pool and a 3.5 A pool over the same
    # structure produce keys that are directly comparable, just sampled
    # differently (ADR-037).
    cross_method: Optional[str] = None
    #: The chain the caller named for each side. Both are user-specified for
    #: every cross method (ADR-049).
    side_a_chain: Optional[str] = None
    side_b_chain: Optional[str] = None
    side_a_selector: Optional[str] = None
    side_b_selector: Optional[str] = None
    contact_cutoff: Optional[float] = None
    contact_filter: Optional[str] = None
    #: The contact rule spelled out scientifically rather than as a mechanism
    #: flag: ``"symmetric_atom_contact_filter"`` (the standard -- an atom is kept
    #: iff it has at least one opposite-side atom within the cutoff, applied to
    #: both sides) versus e.g. ``"protein_side_only_contact_filter"`` (historical
    #: ``Drug_TSR.py``) or ``"no_contact_filter"``.
    contact_filter_rule: Optional[str] = None
    #: How a contact becomes a TSR vertex: ``"atom_level_contact_pool"`` -- the
    #: contacting heavy atom itself is the vertex -- or
    #: ``"heavy_atom_interface_to_ca_projection"``, where a residue is classified
    #: as an interface residue by ANY heavy-atom pair within the cutoff and is
    #: then represented by its Ca, which itself need not be within the cutoff.
    contact_projection: Optional[str] = None
    pool_identity: Optional[str] = None
    pool_order: Optional[str] = None

    def to_dict(self) -> Mapping[str, Any]:
        return dict(asdict(self))

    #: The two fields that carry a physical quantity rather than a name.  They
    #: are normalised to ``float`` on the way in, for the reason stated at
    #: length in :meth:`from_dict`.
    _NUMERIC_FIELDS = ("size_filter", "contact_cutoff")

    @classmethod
    def from_dict(cls, payload: Any) -> "SamplingProfile":
        """Rebuild a sampling profile from :meth:`to_dict`.  The exact inverse.

        **The numeric normalisation is the load-bearing part.**  ``fingerprint``
        hashes through ``json.dumps``, which spells the integer ``3`` as ``3`` and
        the float ``3.0`` as ``3.0``.  Two runs with a 3 A contact cutoff would
        therefore carry different sampling fingerprints, and be refused as
        incomparable, purely because one caller wrote ``3`` and the other
        ``3.0`` -- a difference with no science in it whatsoever.

        The package already forecloses this on the write side: ``SizeFilter``
        calls ``float()`` on its cutoff and ``ContactPoolPolicy`` calls
        ``float()`` on its own, the latter with a comment naming this exact
        hazard.  A reader that did not do the same would re-open it for anyone
        who hand-writes or hand-edits a sidecar, so it is done here too, on the
        same two fields and no others.

        ``None`` stays ``None``: "no filter" and "a filter of 0.0" are different
        facts, and the second is rejected by ``SizeFilter`` anyway.
        """
        payload = _require_mapping(payload, "sampling")
        names = tuple(f.name for f in fields(cls))
        unknown = sorted(set(payload) - set(names))
        if unknown:
            raise MalformedResult(
                "sampling carries %d field(s) this version of biotsr does not know: "
                "%s. A sampling profile is only meaningful if every field that "
                "defines it is understood, so this is refused rather than read "
                "partially." % (len(unknown), ", ".join(repr(u) for u in unknown))
            )
        values = {}
        for f in fields(cls):
            name = f.name
            if name not in payload:
                raise MalformedResult(
                    "sampling is missing the required field %r; a sampling "
                    "profile is serialized in full, including its nulls, so a "
                    "missing field means a truncated or foreign payload" % name
                )
            value = payload[name]
            if value is None:
                values[name] = None
                continue
            if name in cls._NUMERIC_FIELDS:
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    raise MalformedResult(
                        "sampling field %r must be a number or null, got %r (%s)"
                        % (name, value, type(value).__name__)
                    )
                if not value > 0:
                    raise MalformedResult(
                        "sampling field %r must be positive, got %r" % (name, value)
                    )
                values[name] = float(value)
                continue
            if not isinstance(value, str):
                raise MalformedResult(
                    "sampling field %r must be a string or null, got %r (%s)"
                    % (name, value, type(value).__name__)
                )
            values[name] = value
        return cls(**values)

    @property
    def fingerprint(self) -> str:
        return _digest(dict(self.to_dict()))

    def require_compatible(self, other: "SamplingProfile") -> None:
        if self.fingerprint != other.fingerprint:
            a, b = dict(self.to_dict()), dict(other.to_dict())
            diffs = [
                "%s (%r vs %r)" % (f, a[f], b[f]) for f in sorted(a) if a[f] != b[f]
            ]
            raise IncompatibleSampling(
                "results share a key space but were not sampled the same way; "
                "their frequency vectors are not comparable. Differences: "
                + "; ".join(diffs)
            )
