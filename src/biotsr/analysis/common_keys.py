"""Common keys -- four distinct operations that legacy code gives one name.

What "common keys" means is taken from source, not guessed.  Source says it means
four different things, implemented by four programs that live in four sibling
directories and share a word:

1. **Global intersection** (``common_keys.py``) -- the keys present in *every*
   sample, plus, per sample, its total key count and the sum of its frequencies
   over that intersection.  Presence only; frequency never decides membership.
2. **Group-wise intersection** (``common_keys_chain.py``, and
   ``common_keys_residue.py``, which is the same code with a different grouping
   token) -- definition 1 applied within each group separately.
3. **Occurrence percentage** (``common_keys_percent.py``,
   ``common_keys_percent_res.py``) -- for every key, in how many samples of a
   group it occurs at all, as a percentage.  This is a *per-key* table, not a
   subset: it answers "how common is this key", where definition 1 answers
   "which keys are universal".
4. **One-vs-all occurrence percentage** (``common_keys_percent_one_vs_all.py``)
   -- definition 3 for every group at once, joined into one wide table whose
   missing cells are ``0``.

They are kept under four names because collapsing them behind one
``common_keys(mode=...)`` would let a caller who wanted a subset silently get a
percentage table.  Every one of them is presence-based; none of them weights by
frequency (ADR-062).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional, Sequence, Tuple

from .samples import AnalysisSample, check_comparable, sample_provenance

__all__ = [
    "SampleCommonKeyRow",
    "CommonKeyIntersection",
    "KeyOccurrenceRow",
    "KeyOccurrenceTable",
    "OneVsAllOccurrence",
    "common_key_intersection",
    "common_key_intersection_by_group",
    "key_occurrence_percent",
    "key_occurrence_percent_one_vs_all",
]


@dataclass(frozen=True)
class SampleCommonKeyRow:
    """One sample's line of the intersection report."""

    sample_id: str
    total_keys: int
    n_common_keys: int
    sum_common_key_frequency: int


@dataclass(frozen=True)
class CommonKeyIntersection:
    """Definition 1/2: the keys every sample in the set has.

    ``keys`` is numerically ascending.  Legacy returned ``list(set(...) & set(...))``,
    whose order is hash-dependent; the *membership* is what it defines and that is
    reproduced exactly.
    """

    keys: Tuple[int, ...]
    rows: Tuple[SampleCommonKeyRow, ...]
    n_samples: int
    group: Optional[str] = None
    #: Per-sample identities, in sample order (:func:`sample_provenance`).
    samples: Tuple[Mapping[str, Any], ...] = ()
    #: How the set compares -- key space, processing, methods, selection
    #: scopes (:func:`biotsr.analysis.samples.check_comparable`).
    comparison: Mapping[str, Any] = None  # type: ignore[assignment]
    warnings: Tuple[str, ...] = ()

    @property
    def n_common_keys(self) -> int:
        return len(self.keys)

    @property
    def is_empty(self) -> bool:
        return not self.keys

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "group": self.group,
            "n_samples": self.n_samples,
            "n_common_keys": self.n_common_keys,
            "key_space_fingerprint": (
                self.samples[0]["key_space_fingerprint"] if self.samples else None
            ),
            "comparison": dict(self.comparison or {}),
            "samples": [dict(p) for p in self.samples],
            "warnings": list(self.warnings),
            "keys": list(self.keys),
            "rows": [
                {
                    "fileName": r.sample_id,
                    "total_keys": r.total_keys,
                    "common_keys": r.n_common_keys,
                    "sum_common_keys_freq": r.sum_common_key_frequency,
                }
                for r in self.rows
            ],
        }


def _guard(samples: Sequence[AnalysisSample], allow_mixed_processing: bool,
           allow_partial: bool = False):
    """The comparison rule every analysis applies (:func:`check_comparable`).

    Key space always; processing settings unless
    ``allow_mixed_processing=True``.  What was selected -- structures, chains,
    residue types, molecules, residue numbers, atoms -- is recorded and never
    checked.  Returns ``(summary, warnings)`` for the report.
    """
    return check_comparable(samples, allow_mixed_processing=allow_mixed_processing,
                            allow_partial=allow_partial)


def common_key_intersection(
    samples: Sequence[AnalysisSample],
    group: Optional[str] = None,
    *,
    allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> CommonKeyIntersection:
    """Keys present in every sample, with the per-sample report columns.

    The samples must share a key space, always, and processing settings,
    unless ``allow_mixed_processing=True`` -- in which case the report says so
    and names what differs.  Each sample's fingerprints, settings and
    selection are recorded in the result either way.

    An empty intersection is a result, not an error: ``1UBQ``'s 75 residue-level
    amino-acid samples share no key at all, and the report still counts each
    sample's totals.
    """
    samples = list(samples)
    if not samples:
        raise ValueError("common keys over an empty sample set is undefined")
    summary, warnings = _guard(samples, allow_mixed_processing, allow_partial)
    common: set = set(samples[0].counts)
    for sample in samples[1:]:
        common &= set(sample.counts)
        if not common:
            break
    keys = tuple(sorted(common))
    rows = tuple(
        SampleCommonKeyRow(
            sample_id=s.sample_id,
            total_keys=len(s.counts),
            n_common_keys=len(keys),
            sum_common_key_frequency=sum(s.counts[k] for k in keys),
        )
        for s in samples
    )
    return CommonKeyIntersection(
        keys=keys, rows=rows, n_samples=len(samples), group=group,
        samples=tuple(sample_provenance(s) for s in samples),
        comparison=summary, warnings=warnings,
    )


def common_key_intersection_by_group(
    samples: Sequence[AnalysisSample],
    group_of: Callable[[AnalysisSample], str],
    *,
    allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> Mapping[str, CommonKeyIntersection]:
    """Definition 2: definition 1 inside each group, groups in sorted order.

    ``group_of`` is supplied by the caller because the legacy grouping key is a
    filename convention, not a property of the science -- ``common_keys_chain.py``
    takes the trailing underscore token and calls it a chain, while
    ``common_keys_residue.py`` takes the same token and calls it a residue.
    """
    samples = list(samples)
    if not samples:
        raise ValueError("common keys over an empty sample set is undefined")
    _guard(samples, allow_mixed_processing, allow_partial)
    buckets: dict = {}
    for sample in samples:
        buckets.setdefault(group_of(sample), []).append(sample)
    return {
        name: common_key_intersection(
            members, group=name, allow_mixed_processing=allow_mixed_processing,
            allow_partial=allow_partial)
        for name, members in sorted(buckets.items())
    }


@dataclass(frozen=True)
class KeyOccurrenceRow:
    key: int
    occurrence: int
    total_samples: int

    @property
    def percent(self) -> float:
        """``100 * occurrence / total_samples``, unrounded.

        The three legacy programs round this to 0, 4 and 3 decimals
        respectively; rounding is an export concern and lives in the exporters,
        not in the number.
        """
        return (100.0 * self.occurrence) / self.total_samples


@dataclass(frozen=True)
class KeyOccurrenceTable:
    """Definition 3: how many samples each key occurs in, as a percentage."""

    rows: Tuple[KeyOccurrenceRow, ...]
    total_samples: int
    group: Optional[str] = None
    #: Per-sample identities, the comparison summary and any warning, exactly
    #: as :class:`CommonKeyIntersection` records them.
    samples: Tuple[Mapping[str, Any], ...] = ()
    comparison: Mapping[str, Any] = None  # type: ignore[assignment]
    warnings: Tuple[str, ...] = ()

    def percents(self) -> Mapping[int, float]:
        return {r.key: r.percent for r in self.rows}

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "group": self.group,
            "total_samples": self.total_samples,
            "comparison": dict(self.comparison or {}),
            "samples": [dict(p) for p in self.samples],
            "warnings": list(self.warnings),
            "rows": [
                {
                    "key": r.key,
                    "key_occurrence": r.occurrence,
                    "total_files": r.total_samples,
                    "key_percent": r.percent,
                }
                for r in self.rows
            ],
        }


def key_occurrence_percent(
    samples: Sequence[AnalysisSample],
    group: Optional[str] = None,
    *,
    allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> KeyOccurrenceTable:
    """Definition 3.  Presence only -- a key counts once per sample, never by frequency.

    Same comparison rule as :func:`common_key_intersection`.
    """
    samples = list(samples)
    if not samples:
        raise ValueError("key occurrence over an empty sample set is undefined")
    summary, warnings = _guard(samples, allow_mixed_processing, allow_partial)
    occurrence: dict = {}
    for sample in samples:
        for key in sample.counts:
            occurrence[key] = occurrence.get(key, 0) + 1
    total = len(samples)
    rows = tuple(
        KeyOccurrenceRow(key=k, occurrence=occurrence[k], total_samples=total)
        for k in sorted(occurrence)
    )
    return KeyOccurrenceTable(
        rows=rows, total_samples=total, group=group,
        samples=tuple(sample_provenance(s) for s in samples),
        comparison=summary, warnings=warnings,
    )


@dataclass(frozen=True)
class OneVsAllOccurrence:
    """Definition 4: one row per key, one percentage column per group.

    A key absent from a group is ``0.0``, which is what the legacy
    ``fillna(0)`` produces and is the correct reading: it occurred in none of
    that group's samples.
    """

    groups: Tuple[str, ...]
    keys: Tuple[int, ...]
    percents: Mapping[int, Tuple[float, ...]]
    totals: Mapping[str, int]
    #: Per-sample identities, the comparison summary and any warning, exactly
    #: as :class:`CommonKeyIntersection` records them.
    samples: Tuple[Mapping[str, Any], ...] = ()
    comparison: Mapping[str, Any] = None  # type: ignore[assignment]
    warnings: Tuple[str, ...] = ()

    def to_dict(self) -> Mapping[str, Any]:
        return {
            "groups": list(self.groups),
            "totals": dict(self.totals),
            "comparison": dict(self.comparison or {}),
            "samples": [dict(p) for p in self.samples],
            "warnings": list(self.warnings),
            "rows": [
                dict(
                    [("key", k)]
                    + [
                        ("%s_key%%" % g, self.percents[k][i])
                        for i, g in enumerate(self.groups)
                    ]
                )
                for k in self.keys
            ],
        }


def key_occurrence_percent_one_vs_all(
    samples: Sequence[AnalysisSample],
    group_of: Callable[[AnalysisSample], str],
    *,
    allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> OneVsAllOccurrence:
    """Definition 4, groups in sorted order and keys numerically ascending.

    Same comparison rule as :func:`common_key_intersection`.
    """
    samples = list(samples)
    if not samples:
        raise ValueError("key occurrence over an empty sample set is undefined")
    summary, warnings = _guard(samples, allow_mixed_processing, allow_partial)
    buckets: dict = {}
    for sample in samples:
        buckets.setdefault(group_of(sample), []).append(sample)
    groups = tuple(sorted(buckets))
    tables = {g: key_occurrence_percent(buckets[g], group=g,
                                        allow_mixed_processing=allow_mixed_processing,
                                        allow_partial=allow_partial)
              for g in groups}
    keys: set = set()
    for table in tables.values():
        keys.update(r.key for r in table.rows)
    ordered = tuple(sorted(keys))
    per_group = {g: tables[g].percents() for g in groups}
    percents = {
        k: tuple(per_group[g].get(k, 0.0) for g in groups) for k in ordered
    }
    return OneVsAllOccurrence(
        groups=groups,
        keys=ordered,
        percents=percents,
        totals={g: tables[g].total_samples for g in groups},
        samples=tuple(sample_provenance(s) for s in samples),
        comparison=summary,
        warnings=warnings,
    )
