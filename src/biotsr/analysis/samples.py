"""The unit of analysis, and the comparison rule that guards every comparison.

An ``AnalysisSample`` is one key-frequency vector plus the scientific identities
it was produced under.  Nothing above this module knows which engine made the
keys: a Ca ``TSRResult``, an atomic one, a mirror-signed one and a cross one all
reduce to the same things -- an id, ``key -> count``, the ``KeySpace`` and the
``SamplingProfile``, and the recorded selection.  That is what makes the
analysis layer method-agnostic.

**One comparison rule, used by every analysis** (ADR-094).  Feature matrices,
common keys, occurrence tables, similarity and clustering all apply
:func:`check_comparable`:

1. **Key definitions -- always.**  The samples must share one ``KeySpace``:
   the same number must mean the same triangle -- the same point and label
   conventions, the same theta and distance bins, the same encoding and
   mirror convention.  This is decided from the scientific key-space
   definition, never from method names: different methods that share a key
   space (``aa_ind`` and ``drug`` both label atoms by element at radix 112
   in 29 x 18 bins) are comparable; the same method with a different modifier
   (``mirror``, ``aa_grouping``) is not.  No option relaxes this.
2. **Processing -- checked separately, and acknowledged when it differs.**
   The settings that change how triangles are drawn from whatever was
   selected (:data:`PROCESSING_FIELDS`: size filter, contact cutoff and
   contact rule, completeness rule, candidate strategy, degenerate-triangle
   policy, pool rules, and the alternate-location rule) must agree, unless
   the caller passes ``allow_mixed_processing=True``
   (``--allow-mixed-processing``).  Then the comparison proceeds, the report
   records every sample's settings, and a warning names each differing
   setting.  The opt-in never relaxes step 1.

3. **Partial samples -- acknowledged, separately.**  A sample whose
   selection did not get everything it asked for (a requested residue, atom
   or molecule absent from the structure, a modified residue the method
   cannot read, one the completeness rule rejected) is not the selection
   that was requested.  It is never dropped and never passed off as
   complete: the comparison is refused, naming each partial sample and what
   it lacks, unless the caller passes ``allow_partial=True``
   (``--allow-partial``).  Then it proceeds, the report records each
   sample's omissions, and a warning names them (ADR-106).  This
   acknowledgement is independent of ``allow_mixed_processing``: neither
   implies the other, and neither relaxes step 1.

**What never blocks a comparison** is what was selected: the structure, the
chain(s), amino-acid types (LYS-only versus TYR-only), ligands,
carbohydrates, nucleotides, residue numbers, atom subsets, a single copy
versus all copies, and the method name itself.  Those are recorded per
sample (``selection``, ``method``, :data:`TARGET_FIELDS`) and reported.  Two
samples need not share any key: an empty intersection is a valid result.
A *complete* LYS-only sample and a *complete* TYR-only sample compare
without any acknowledgement; only a selection that fell short of its own
request needs one.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

from ..domain.errors import TSRError
from ..domain.keyspace import (
    IncompatibleKeySpace,
    IncompatibleSampling,
    KeySpace,
    SamplingProfile,
)
from ..domain.result import TSRResult
from ..domain.selection_record import selection_omissions, selection_sides
from ..domain.interpretation import (
    STRUCTURE_INTERPRETATION,
    IncompatibleStructureInterpretation,
    require_common_structure_interpretation,
    structure_interpretation_of,
)

__all__ = [
    "require_common_structure_interpretation",
    "AnalysisSample",
    "IncompatibleProcessing",
    "PartialSamples",
    "PROCESSING_FIELDS",
    "TARGET_FIELDS",
    "label_from_name",
    "require_common_key_space",
    "require_common_sampling",
    "require_common_processing",
    "check_comparable",
    "describe_omissions",
    "processing_policy",
    "processing_fingerprint",
    "sample_provenance",
    "sampling_variation",
]


class IncompatibleProcessing(IncompatibleSampling):
    """Raised when samples with one key space were processed differently.

    A subclass of ``IncompatibleSampling``, so code that catches that keeps
    working.  ``allow_mixed_processing=True`` turns the refusal into a recorded
    warning; nothing turns a key-space refusal into one.
    """


class PartialSamples(TSRError):
    """Raised when an analysis would use samples whose selection fell short.

    Not an incompatibility: a partial sample has the same key definitions and
    may be compared deliberately -- ``allow_partial=True`` (``--allow-partial``)
    proceeds and records what each one lacks.  It is a separate
    acknowledgement from ``allow_mixed_processing``, and neither relaxes the
    key-space rule.  ``partial`` maps each partial sample's id to its
    omissions (:func:`biotsr.domain.selection_record.selection_omissions`).
    """

    def __init__(self, message: str, partial: Mapping[str, Any]):
        super().__init__(message)
        self.partial = dict(partial)


#: ``SamplingProfile`` fields that are **processing policy**: how triangles are
#: drawn from whatever was selected.  A difference here must be acknowledged.
PROCESSING_FIELDS = (
    "candidate_strategy",
    "completeness_policy",
    "composition_rule",
    "size_filter",
    "degenerate_policy",
    "seq_label_policy",
    "contact_cutoff",
    "contact_filter",
    "contact_filter_rule",
    "contact_projection",
    "pool_identity",
    "pool_order",
)

#: ``SamplingProfile`` fields that identify **what was selected**: the
#: selector, its scope, the cross method and each side's chain and free-text
#: description (which names the ligand or residues).  Recorded, never a reason
#: to refuse.  A test keeps the two tuples a partition of the profile's fields.
TARGET_FIELDS = (
    "selector",
    "scope",
    "cross_method",
    "side_a_chain",
    "side_b_chain",
    "side_a_selector",
    "side_b_selector",
)


@dataclass(frozen=True)
class AnalysisSample:
    """One key-frequency vector, carrying both identities.

    ``counts`` is sparse by construction: a key absent from the mapping has
    frequency zero.  That is the same fact the legacy dense CSV writes as a
    literal ``0`` column, and no analysis in this layer distinguishes the two
    (ADR-058).

    Because they are the same fact, a zero passed in explicitly is **dropped** at
    construction.  Keeping it would break the equivalence in a way that shows: the
    presence-based measures define presence as ``value > 0`` -- legacy spells it
    ``a_jac[a_jac > 0] = 1`` -- while ``key in counts`` would call a stored zero
    present, so ``binary_jaccard`` and ``sarika`` would disagree with
    ``generalized_jaccard`` about the same pair.  Canonicalising once here is why
    every measure can read the mapping directly.
    """

    sample_id: str
    counts: Mapping[int, int]
    key_space: KeySpace
    sampling: SamplingProfile
    label: Optional[str] = None
    annotations: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for key, value in self.counts.items():
            if not isinstance(key, int) or isinstance(key, bool):
                raise TypeError(
                    "sample %r: key %r is %s, not int -- keys are integer "
                    "coordinates and are never compared as strings"
                    % (self.sample_id, key, type(key).__name__)
                )
            if not isinstance(value, int) or isinstance(value, bool):
                raise TypeError(
                    "sample %r: frequency for key %d is %r, not int"
                    % (self.sample_id, key, value)
                )
            if value < 0:
                raise ValueError(
                    "sample %r: frequency for key %d is negative (%d)"
                    % (self.sample_id, key, value)
                )
        if not all(self.counts.values()):
            object.__setattr__(
                self, "counts", {k: v for k, v in self.counts.items() if v}
            )

    @classmethod
    def from_result(
        cls,
        result: TSRResult,
        sample_id: Optional[str] = None,
        label: Optional[str] = None,
    ) -> "AnalysisSample":
        """Reduce any ``TSRResult`` to a sample, whatever method produced it.

        The per-category breakdown is summed: every evidenced analysis operates
        on the total frequency of a key, and none of them splits mirror from
        regular or one SSE bucket from another.  The breakdown is not discarded --
        it stays on the originating result -- it is simply not what these
        measures are defined over.
        """
        name = sample_id
        if name is None:
            name = "%s_%s" % (result.structure_id.value, result.chain) if result.chain else result.structure_id.value
        annotations = {
            "method": result.method,
            "structure_id": result.structure_id.value,
            "chain": result.chain,
            # The alternate-location rule that admitted the atoms: part of the
            # processing record, since it decides which atoms exist.
            "altloc_policy": result.source.altloc_policy,
            # Which modified amino acids were read as protein residues: part
            # of the processing record too (ADR-105).  ``None`` for a result
            # written before it was recorded -- reported as not recorded,
            # never taken to be the current rule.
            "modified_residue_policy": result.source.modified_residue_policy,
            # How the records were read and checked (ADR-108); ``None`` for a
            # result written by 1.0.0rc5 or earlier.
            "input_rules": getattr(result.source, "input_rules", None),
        }
        sides = selection_sides(result.annotations.get("selection_request"))
        if sides is not None:
            annotations["selection"] = sides
        if result.annotations.get("target") is not None:
            annotations["target"] = result.annotations.get("target")
        if result.warnings:
            annotations["result_warnings"] = list(result.warnings)
        # Which reading of the structure file produced these keys.  Carried
        # here, not in either identity, because two results can share both
        # fingerprints exactly and still disagree about every key when the
        # atom-name semantics differ (ADR-079).  Derived, never stored, so a
        # result written before this existed still answers.
        try:
            annotations[STRUCTURE_INTERPRETATION] = structure_interpretation_of(result)
        except IncompatibleStructureInterpretation:
            # A loaded result whose provenance cannot state an interpretation
            # is not silently labelled with the shipped one.
            pass
        return cls(
            sample_id=name,
            counts=dict(result.frequencies),
            key_space=result.key_space,
            sampling=result.sampling,
            label=label,
            annotations=annotations,
        )

    @property
    def omissions(self) -> Mapping[str, Any]:
        """What this sample's selection asked for and did not get; ``{}`` if nothing.

        Read from the structured selection record, and from the result's
        warnings for a result written before that record listed omissions
        (:func:`biotsr.domain.selection_record.selection_omissions`).  A
        sample built by hand from counts records no selection and has none.
        """
        return selection_omissions(self.annotations.get("selection"),
                                   self.annotations.get("result_warnings") or ())

    @property
    def is_partial(self) -> bool:
        return bool(self.omissions)

    @property
    def n_keys(self) -> int:
        """Number of distinct keys present, i.e. non-zero columns."""
        return len(self.counts)

    @property
    def total_frequency(self) -> int:
        return sum(self.counts.values())

    def provenance(self) -> Mapping[str, Any]:
        return {
            "sample_id": self.sample_id,
            "label": self.label,
            "key_space_fingerprint": self.key_space.fingerprint,
            "sampling_fingerprint": self.sampling.fingerprint,
            "n_keys": self.n_keys,
            "total_frequency": self.total_frequency,
            "annotations": dict(self.annotations),
        }


def label_from_name(name: str, token: int, separator: str = "_") -> str:
    """The evidenced class-label rule: the *n*-th separator token of the id.

    Every legacy consumer that needs a class recovers it this way and they
    disagree only about which token -- ``clustermap.py`` takes token 3 of
    ``PDB_chain_seq_RES_SITE``, ``clustermap_chain_hexamer.py`` takes token 2 of
    ``PDB_chain_RES_seq_x``, and the ``ai_ml`` trainers take the last one.  The
    index is therefore a parameter of the dataset's naming convention and never
    a default.  Nothing in BioTSR calls this for you: a class label is whatever
    the caller supplies, and a table matched by sample id is usually the safer
    source (``docs/tutorials/08_similarity_clustering.md``).
    """
    parts = name.split(separator)
    if token < 0:
        token += len(parts)
    if not 0 <= token < len(parts):
        raise ValueError(
            "sample id %r has %d %r-separated tokens; token %d does not exist"
            % (name, len(parts), separator, token)
        )
    return parts[token]


def require_common_key_space(samples: Sequence[AnalysisSample]) -> KeySpace:
    """Refuse to compare key coordinates that do not mean the same thing.

    Raises ``IncompatibleKeySpace`` naming every differing field, so the caller
    learns *which* semantic diverged rather than only that something did.
    """
    if not samples:
        raise ValueError("no samples: a key space cannot be established")
    first = samples[0]
    for other in samples[1:]:
        try:
            first.key_space.require_compatible(other.key_space)
        except IncompatibleKeySpace as exc:
            raise IncompatibleKeySpace(
                "samples %r and %r cannot be compared. %s"
                % (first.sample_id, other.sample_id, exc)
            ) from None
    return first.key_space


def require_common_sampling(samples: Sequence[AnalysisSample]) -> SamplingProfile:
    """Opt-in: demand identical sampling as well as an identical key space.

    Not the default.  Only an analysis whose own specification says the counts
    must have been drawn the same way should call this (ADR-055).
    """
    if not samples:
        raise ValueError("no samples: a sampling profile cannot be established")
    first = samples[0]
    for other in samples[1:]:
        try:
            first.sampling.require_compatible(other.sampling)
        except IncompatibleSampling as exc:
            raise IncompatibleSampling(
                "samples %r and %r were not sampled the same way. %s"
                % (first.sample_id, other.sample_id, exc)
            ) from None
    return first.sampling


def processing_policy(sample: AnalysisSample) -> Mapping[str, Any]:
    """The processing record two samples must share to be compared without opt-in.

    Every :data:`PROCESSING_FIELDS` value of the sample's ``SamplingProfile``,
    and the alternate-location and modified-residue rules its result recorded
    (``None`` for a sample built by hand from counts, which records neither,
    and for a result written before the modified-residue rule was recorded).
    """
    sampling = sample.sampling.to_dict()
    return {
        "sampling": {name: sampling.get(name) for name in PROCESSING_FIELDS},
        "altloc_policy": sample.annotations.get("altloc_policy"),
        "modified_residue_policy": sample.annotations.get("modified_residue_policy"),
        "input_rules": sample.annotations.get("input_rules"),
    }


def processing_fingerprint(sample: AnalysisSample) -> str:
    """SHA-256 over :func:`processing_policy`, canonically serialized."""
    blob = json.dumps(processing_policy(sample), sort_keys=True,
                      separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _flat_processing(sample: AnalysisSample) -> Mapping[str, Any]:
    policy = processing_policy(sample)
    out = dict(policy["sampling"])
    out["altloc_policy"] = policy["altloc_policy"]
    out["modified_residue_policy"] = policy["modified_residue_policy"]
    out["input_rules"] = policy["input_rules"]
    return out


def _shown(value: Any) -> str:
    return "not recorded" if value is None else repr(value)


def _differences(a: AnalysisSample, b: AnalysisSample) -> Tuple[str, ...]:
    fa, fb = _flat_processing(a), _flat_processing(b)
    return tuple(
        "%s (%s vs %s)" % (name, _shown(fa[name]), _shown(fb[name]))
        for name in sorted(fa) if fa[name] != fb[name]
    )


def _varying(samples: Sequence[AnalysisSample]) -> Tuple[str, ...]:
    seen: dict = {}
    for sample in samples:
        for name, value in _flat_processing(sample).items():
            seen.setdefault(name, set()).add(json.dumps(value, default=str))
    return tuple(sorted(name for name, values in seen.items() if len(values) > 1))


def require_common_processing(samples: Sequence[AnalysisSample]) -> None:
    """Demand that every sample was processed by the same settings.

    Samples that differ only in what was selected -- the structure, chains,
    residue types, molecules, residue numbers, atoms, or the method that
    selected them -- pass.  A different size filter, contact cutoff, contact
    rule, completeness rule, candidate strategy, pool rule or
    alternate-location rule raises :class:`IncompatibleProcessing` naming the
    settings.
    """
    if not samples:
        raise ValueError("no samples: a processing policy cannot be established")
    first = samples[0]
    reference = processing_fingerprint(first)
    for other in samples[1:]:
        if processing_fingerprint(other) == reference:
            continue
        raise IncompatibleProcessing(
            "samples %r and %r share key definitions but were processed with "
            "different settings: %s. What was selected may differ freely "
            "(structures, chains, residue types, molecules, residue numbers, "
            "atoms); a processing difference changes how the frequencies were "
            "drawn, so it must be acknowledged: allow_mixed_processing=True "
            "(--allow-mixed-processing), which records the difference in the "
            "report." % (first.sample_id, other.sample_id,
                         "; ".join(_differences(first, other)))
        )


def describe_omissions(omissions: Mapping[str, Any]) -> str:
    """One line saying what a partial sample lacks, for messages and warnings."""
    parts = []
    for side, entry in omissions.items():
        if side == "messages":
            continue
        words = []
        for kind, label in (("not_found", "not found"),
                            ("unsupported", "not supported or not represented"),
                            ("declined", "rejected by the completeness rule")):
            if entry.get(kind):
                words.append("%s: %s" % (label, ", ".join(entry[kind])))
        if entry.get("atoms_missing"):
            words.append("atoms missing: %s" % "; ".join(
                "%s lacks %s" % (k, ", ".join(v))
                for k, v in entry["atoms_missing"].items()))
        if words:
            parts.append("side %s: %s" % (side, "; ".join(words))
                         if side != "primary" else "; ".join(words))
    if not parts:
        parts = list(omissions.get("messages") or ())
    return "; ".join(parts)


def _partial_samples(samples: Sequence[AnalysisSample]) -> Mapping[str, Any]:
    return {s.sample_id: dict(s.omissions) for s in samples if s.is_partial}


def check_comparable(
    samples: Sequence[AnalysisSample], *, allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> Tuple[Mapping[str, Any], Tuple[str, ...]]:
    """The comparison rule every analysis applies; returns ``(summary, warnings)``.

    1. A common key space is required, always (``IncompatibleKeySpace``).
    2. A common processing policy is required (:class:`IncompatibleProcessing`)
       unless ``allow_mixed_processing=True``; then the summary says the set is
       mixed and a warning names every setting that differs.
    3. Partial samples are refused (:class:`PartialSamples`) unless
       ``allow_partial=True``; then the summary lists each one's omissions and
       a warning names them.  Independent of step 2.

    Nothing about what was selected is checked.  ``summary`` records the key
    space, the processing fingerprints and which samples share each, the
    differing settings, and -- for the reader -- the methods and selection
    scopes in the set and how many distinct full sampling fingerprints it holds
    (those include the selector and scope, and for cross methods the chains
    and side descriptions, so they can differ while the processing agrees).
    """
    samples = list(samples)
    key_space = require_common_key_space(samples)
    groups: dict = {}
    for sample in samples:
        groups.setdefault(processing_fingerprint(sample), []).append(sample.sample_id)
    varying = _varying(samples)
    warnings = []
    if len(groups) > 1:
        if not allow_mixed_processing:
            require_common_processing(samples)
        warnings.append(
            "the samples were not all processed with the same settings (%d "
            "processing policies; differing in: %s). Mixed processing was "
            "explicitly allowed, so the comparison is reported, but frequencies "
            "drawn under different settings answer different questions."
            % (len(groups), ", ".join(varying))
        )
    partial = _partial_samples(samples)
    if partial:
        listed = "; ".join("%s (%s)" % (sid, describe_omissions(om))
                           for sid, om in partial.items())
        if not allow_partial:
            raise PartialSamples(
                "%d of the %d samples are partial -- their selections did not get "
                "everything they asked for: %s. A partial sample is not the "
                "selection that was requested, so using it must be acknowledged: "
                "allow_partial=True (--allow-partial), which records what each "
                "one lacks in the report. That is separate from mixed processing "
                "and never relaxes the key definitions."
                % (len(partial), len(samples), listed), partial)
        warnings.append(
            "%d of the %d samples are partial and were explicitly allowed: %s. "
            "They are compared as computed, not as the complete selections that "
            "were requested." % (len(partial), len(samples), listed))
    scopes = sorted({side.get("scope") for s in samples
                     for side in (s.annotations.get("selection") or ())
                     if side.get("scope")})
    summary = {
        "key_space_fingerprint": key_space.fingerprint,
        "processing_uniform": len(groups) <= 1,
        "mixed_processing_allowed": bool(allow_mixed_processing),
        "n_processing_policies": len(groups),
        "samples_by_processing_fingerprint": {k: list(v) for k, v in sorted(groups.items())},
        "processing_differences": list(varying),
        "methods": sorted({str(s.annotations.get("method"))
                           for s in samples if s.annotations.get("method")}),
        "selection_scopes": scopes,
        "n_sampling_fingerprints": len({s.sampling.fingerprint for s in samples}),
        "partial_allowed": bool(allow_partial),
        "partial_samples": {k: dict(v) for k, v in partial.items()},
    }
    return summary, tuple(warnings)


def sample_provenance(sample: AnalysisSample) -> Mapping[str, Any]:
    """What an analysis report records about one sample.

    Its identities (the key-space, full sampling and processing fingerprints),
    its processing settings, the sampling settings that are set, the selection
    it actually used (per side: chain, scope, arguments, instances), and the
    warnings its own result carried.
    """
    settings = {k: v for k, v in sorted(sample.sampling.to_dict().items())
                if v is not None}
    return {
        "sample_id": sample.sample_id,
        "label": sample.label,
        "method": sample.annotations.get("method"),
        "structure_id": sample.annotations.get("structure_id"),
        "key_space_fingerprint": sample.key_space.fingerprint,
        "sampling_fingerprint": sample.sampling.fingerprint,
        "processing_fingerprint": processing_fingerprint(sample),
        "processing": processing_policy(sample),
        "sampling": settings,
        "selection": sample.annotations.get("selection"),
        "target": sample.annotations.get("target"),
        "chains": {k: getattr(sample.sampling, k)
                   for k in ("side_a_chain", "side_b_chain")
                   if getattr(sample.sampling, k) is not None}
        or ({"chain": sample.annotations["chain"]}
            if sample.annotations.get("chain") else None),
        "n_keys": sample.n_keys,
        "total_frequency": sample.total_frequency,
        "result_warnings": list(sample.annotations.get("result_warnings") or ()),
        # What the selection asked for and did not get; ``{}`` when complete.
        "partial": sample.is_partial,
        "omissions": dict(sample.omissions),
    }


def sampling_variation(samples: Iterable[AnalysisSample]) -> Mapping[str, Any]:
    """Describe how sampling differs across a set, without judging it.

    This is what an analysis result exposes so a reader can decide whether a
    comparison across, say, two contact cutoffs answers their question.
    """
    fingerprints: dict = {}
    fields: dict = {}
    for sample in samples:
        payload = dict(sample.sampling.to_dict())
        fingerprints.setdefault(sample.sampling.fingerprint, []).append(sample.sample_id)
        for name, value in payload.items():
            fields.setdefault(name, set()).add(repr(value))
    varying = sorted(name for name, seen in fields.items() if len(seen) > 1)
    return {
        "n_distinct_profiles": len(fingerprints),
        "uniform": len(fingerprints) <= 1,
        "varying_fields": tuple(varying),
        "samples_by_fingerprint": {k: tuple(v) for k, v in fingerprints.items()},
    }
