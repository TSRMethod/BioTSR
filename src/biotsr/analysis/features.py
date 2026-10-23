"""The shared feature matrix: samples x key coordinates, deterministically ordered.

Two column orderings are offered, and they are *the same mathematical object
under a permutation*, because every measure built on top is a sum over all
columns (ADR-056):

* ``numeric_ascending`` -- the standard.  Integer keys sorted as integers, so
  ``-713374251 < 9 < 10``.  No legacy program does this; it is adopted because
  the alternatives are respectively non-reproducible and mis-sorted, and because
  it changes no number that any evidenced analysis produces (ADR-056).
* ``first_encounter`` -- ``hsp70_actin``'s insertion order, driven by the order
  of the feature-selection file.

Historical programs also sorted keys as decimal strings, and some used the
``list(set(...))``
hash order of ``Kinases-and-Phosphatases-Clustering`` and ``Nucleotide_Analysis``.
the latter is not reproducible (it changes with ``PYTHONHASHSEED``), and neither
is offered.

Storage is sparse and results are dense-identical.  A row is a mapping; a missing
key *is* a zero, and ``row()`` materialises the dense vector on demand (ADR-058).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

from ..domain.keyspace import KeySpace
from .samples import AnalysisSample, check_comparable, sample_provenance

__all__ = [
    "FEATURE_ORDERINGS",
    "NUMERIC_ASCENDING",
    "FIRST_ENCOUNTER",
    "FeatureMatrix",
    "feature_matrix",
]

#: Standard: integer keys ordered as integers.
NUMERIC_ASCENDING = "numeric_ascending"
#: Order of first appearance across samples, in sample order.
FIRST_ENCOUNTER = "first_encounter"

FEATURE_ORDERINGS = (NUMERIC_ASCENDING, FIRST_ENCOUNTER)


def _ordered_columns(samples: Sequence[AnalysisSample], ordering: str) -> Tuple[int, ...]:
    if ordering == NUMERIC_ASCENDING:
        union: set = set()
        for sample in samples:
            union.update(sample.counts)
        return tuple(sorted(union))
    if ordering == FIRST_ENCOUNTER:
        seen: dict = {}
        for sample in samples:
            for key in sample.counts:
                if key not in seen:
                    seen[key] = None
        return tuple(seen)
    raise ValueError(
        "unknown feature ordering %r; choose one of %r. A hash order "
        "(list(set(...))) is not offered: it is not reproducible."
        % (ordering, FEATURE_ORDERINGS)
    )


@dataclass(frozen=True)
class FeatureMatrix:
    """Samples x keys, with both identities preserved per sample.

    ``columns`` is the union of keys over all samples, in ``ordering``.  Rows are
    stored sparsely; ``row(i)`` and ``dense()`` expand them with explicit zeros.
    """

    sample_ids: Tuple[str, ...]
    columns: Tuple[int, ...]
    rows: Tuple[Mapping[int, int], ...]
    key_space: KeySpace
    sampling_fingerprints: Tuple[str, ...]
    ordering: str
    labels: Tuple[Optional[str], ...]
    #: How the rows compare: key space, processing, methods, selection scopes
    #: (``check_comparable``'s summary).
    comparison: Mapping[str, Any]
    #: Per row, in row order: identities, settings, selection and warnings
    #: (:func:`biotsr.analysis.samples.sample_provenance`).
    samples: Tuple[Mapping[str, Any], ...] = ()
    #: Set when mixed processing was explicitly allowed.
    warnings: Tuple[str, ...] = ()

    # -- shape ---------------------------------------------------------------
    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def n_columns(self) -> int:
        return len(self.columns)

    @property
    def shape(self) -> Tuple[int, int]:
        return (self.n_samples, self.n_columns)

    @property
    def nnz(self) -> int:
        return sum(len(r) for r in self.rows)

    @property
    def density(self) -> float:
        cells = self.n_samples * self.n_columns
        return (self.nnz / cells) if cells else 0.0

    # -- access --------------------------------------------------------------
    def index_of(self, sample_id: str) -> int:
        try:
            return self.sample_ids.index(sample_id)
        except ValueError:
            raise KeyError("no sample %r in this matrix" % sample_id) from None

    def value(self, i: int, key: int) -> int:
        """The frequency of ``key`` in row ``i``; a key not present is zero."""
        return self.rows[i].get(key, 0)

    def row(self, i: int) -> Tuple[int, ...]:
        """Row ``i`` as a dense integer vector in column order."""
        row = self.rows[i]
        return tuple(row.get(k, 0) for k in self.columns)

    def dense(self) -> Tuple[Tuple[int, ...], ...]:
        return tuple(self.row(i) for i in range(self.n_samples))

    def row_sum(self, i: int) -> int:
        return sum(self.rows[i].values())

    def column_occurrence(self) -> Mapping[int, int]:
        """How many samples each key appears in (non-zero), in column order."""
        counts = {k: 0 for k in self.columns}
        for row in self.rows:
            for key in row:
                counts[key] += 1
        return counts

    # -- export --------------------------------------------------------------
    def to_header_csv(self, name_column: str = "sample_id") -> str:
        """``sample_id,<key>,<key>,...`` then one row of frequencies per sample.

        Columns are the integer keys in :attr:`ordering`; a cell is a key's
        frequency in that sample, ``0`` when the key is absent.
        """
        lines = ["%s,%s\r\n" % (name_column, ",".join(str(k) for k in self.columns))]
        for i, name in enumerate(self.sample_ids):
            lines.append("%s,%s\r\n" % (name, ",".join(str(v) for v in self.row(i))))
        return "".join(lines)

    def provenance(self) -> Mapping[str, Any]:
        return {
            "n_samples": self.n_samples,
            "n_columns": self.n_columns,
            "ordering": self.ordering,
            "nnz": self.nnz,
            "density": self.density,
            "key_space_fingerprint": self.key_space.fingerprint,
            "comparison": dict(self.comparison),
            "sample_ids": list(self.sample_ids),
            "labels": list(self.labels),
            "samples": [dict(p) for p in self.samples],
            "warnings": list(self.warnings),
        }


def feature_matrix(
    samples: Sequence[AnalysisSample],
    ordering: str = NUMERIC_ASCENDING,
    *,
    allow_mixed_processing: bool = False,
    allow_partial: bool = False,
) -> FeatureMatrix:
    """Build the shared feature matrix over ``samples``.

    The same comparison rule as every other analysis
    (:func:`biotsr.analysis.samples.check_comparable`): a common ``KeySpace`` is
    required, always -- the columns are key coordinates (ADR-054) -- and common
    processing settings are required unless ``allow_mixed_processing=True``, in
    which case the matrix records the difference as a warning.  A partial
    sample -- one whose selection did not get everything it asked for -- is
    refused unless ``allow_partial=True``, which records its omissions; that
    is a separate acknowledgement (ADR-106).  Samples from
    different methods sharing the key space, structures, chains, residue types,
    molecules, copies or residues compare freely, and need not share any key.
    """
    samples = list(samples)
    if not samples:
        raise ValueError("a feature matrix needs at least one sample")
    seen: dict = {}
    for sample in samples:
        if sample.sample_id in seen:
            raise ValueError("duplicate sample id %r" % sample.sample_id)
        seen[sample.sample_id] = None
    summary, warnings = check_comparable(
        samples, allow_mixed_processing=allow_mixed_processing,
        allow_partial=allow_partial)
    columns = _ordered_columns(samples, ordering)
    return FeatureMatrix(
        sample_ids=tuple(s.sample_id for s in samples),
        columns=columns,
        rows=tuple(dict(s.counts) for s in samples),
        key_space=samples[0].key_space,
        sampling_fingerprints=tuple(s.sampling.fingerprint for s in samples),
        ordering=ordering,
        labels=tuple(s.label for s in samples),
        comparison=summary,
        samples=tuple(sample_provenance(s) for s in samples),
        warnings=warnings,
    )
