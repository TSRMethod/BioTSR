"""Pairwise similarity and distance over key-frequency vectors.

Five measures are evidenced, by two independent implementations that agree
line for line: the Python ``hsp70_actin/lib/jacOnlyFromVect*.py`` family and the
C kernel ``hsp70_actin/src/c/calc_coeffs.h``.  The more recent
``jaccard_parallel_code`` keeps only the generalized one, and it is the only
measure for which an executable byte-exact reference output exists.  Formulas
below are quoted from source, not from a paper.

    generalized_jaccard   J = sum_k min(a_k,b_k) / sum_k max(a_k,b_k)
    binary_jaccard        J = |A and B| / |A or B|         (presence only)
    wu                    J = sum_k min(a_k,b_k) / min( sum_k max(a_k,b_k),
                                                        max(sum a, sum b) )
    sarika                J = sum_{k: a_k>0 and b_k>0} (a_k + b_k) / (sum a + sum b)
    cosine                J = a.b / (|a| |b|)

Every one of them is defined here as a **similarity**; the distance is ``1 - J``
and the conversion is never implicit.  Legacy stored the distance and called the
file ``jaccard_similarity``, and source settles which it is:
``jaccard_mpi.c`` writes ``1 - J`` unless ``-s`` is given, and the shipped
``generalised.csv`` has a ``0.000`` diagonal, so the historical
"jaccard_similarity" outputs are distances (ADR-057).

**Input.**  A vector is a mapping ``key -> count``.  Counts must be finite and
non-negative (``int`` or ``float``; ``bool`` is refused); anything else raises
``ValueError`` naming the key.  An explicit zero is the same fact as an absent
key and is dropped before any measure sees it, so presence (``binary_jaccard``)
always means ``count > 0``.

**Empty input.**  One rule for every measure: two empty vectors (no key with a
positive count) have ``J = 1.0`` -- they are identical -- and so ``D = 0.0``;
an empty vector against a non-empty one has ``J = 0.0`` (``D = 1.0``).  This is
``jaccard_mpi.c``'s convention (it returns 1 when the denominator is zero).  A
matrix over samples that include an empty one carries a warning saying so,
because "identical because both are empty" is rarely what an analysis wants.

**Scale.**  Every value is a fraction in ``[0, 1]`` unless the matrix was
explicitly :meth:`SimilarityMatrix.scaled`; a scaled matrix records its
``scale`` (100 for percent), its diagonal and conversions follow it
(``D = scale - J``), and every export names it.  Nothing rescales implicitly.
"""

from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from typing import Any, Callable, Mapping, Optional, Sequence, Tuple

from .features import FeatureMatrix

__all__ = [
    "SIMILARITY",
    "DISTANCE",
    "MEASURES",
    "MeasureSpec",
    "SimilarityMatrix",
    "similarity_matrix",
    "distance_matrix",
    "generalized_jaccard",
    "binary_jaccard",
    "wu",
    "sarika",
    "cosine",
]

SIMILARITY = "similarity"
DISTANCE = "distance"


@dataclass(frozen=True)
class MeasureSpec:
    """What a measure is, written down next to the code that computes it."""

    name: str
    formula: str
    input_domain: str
    empty_input: str
    source: str
    has_frozen_reference: bool


def _counts(vector: Mapping[int, Any], which: str) -> Mapping[int, Any]:
    """``vector`` validated, with explicit zeros dropped.

    Negative, non-finite and non-numeric counts are refused rather than
    computed with: ``min``/``max`` over a negative count produce a "similarity"
    outside ``[0, 1]`` that nothing downstream would notice.
    """
    if not isinstance(vector, Mapping):
        raise TypeError("%s must be a mapping key -> count, not %s"
                        % (which, type(vector).__name__))
    clean = {}
    for key, value in vector.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("%s: count for key %r is %r, not a number"
                             % (which, key, value))
        if value != value or value in (math.inf, -math.inf):
            raise ValueError("%s: count for key %r is not finite (%r)"
                             % (which, key, value))
        if value < 0:
            raise ValueError("%s: count for key %r is negative (%r)"
                             % (which, key, value))
        if value:
            clean[key] = value
    return clean


def _sum_min(a: Mapping[int, int], b: Mapping[int, int]) -> int:
    if len(b) < len(a):
        a, b = b, a
    total = 0
    for key, va in a.items():
        vb = b.get(key)
        if vb is not None:
            total += va if va < vb else vb
    return total


def _shared(a: Mapping[int, int], b: Mapping[int, int]):
    if len(b) < len(a):
        a, b = b, a
    for key, va in a.items():
        vb = b.get(key)
        if vb is not None:
            yield va, vb


def _generalized_jaccard(a: Mapping[int, int], b: Mapping[int, int]) -> float:
    """``sum min / sum max`` over the union of keys, weighted by frequency."""
    smin = _sum_min(a, b)
    smax = sum(a.values()) + sum(b.values()) - smin
    # Integer counts give a correctly rounded ratio <= 1; float counts can
    # round a hair above it, and a similarity is never allowed above 1.
    return min(1.0, smin / smax) if smax > 0 else 1.0


def _binary_jaccard(a: Mapping[int, int], b: Mapping[int, int]) -> float:
    """``|A and B| / |A or B|`` on the presence sets; frequencies ignored."""
    inter = 0
    small, large = (a, b) if len(a) <= len(b) else (b, a)
    for key in small:
        if key in large:
            inter += 1
    union = len(a) + len(b) - inter
    return (inter / union) if union > 0 else 1.0


def _wu(a: Mapping[int, int], b: Mapping[int, int]) -> float:
    """``sum min / min(sum max, max(sum a, sum b))``.

    Written exactly as source writes it.  For non-negative counts
    ``sum max >= max(sum a, sum b)`` always, so the inner ``min`` never binds and
    the measure is ``sum min / max(sum a, sum b)``; a unit test pins that
    identity rather than leaving it as a remark.
    """
    sa = sum(a.values())
    sb = sum(b.values())
    smin = _sum_min(a, b)
    smax = sa + sb - smin
    den = smax if smax < (sa if sa > sb else sb) else (sa if sa > sb else sb)
    return min(1.0, smin / den) if den > 0 else 1.0


def _sarika(a: Mapping[int, int], b: Mapping[int, int]) -> float:
    """``sum over co-present keys of (a_k + b_k) / (sum a + sum b)``."""
    num = 0
    for va, vb in _shared(a, b):
        num += va + vb
    den = sum(a.values()) + sum(b.values())
    return min(1.0, num / den) if den > 0 else 1.0


def _cosine(a: Mapping[int, int], b: Mapping[int, int]) -> float:
    """``a.b / (|a| |b|)``.

    Evidenced only by the Python family; the C kernel's cosine slot is a
    disabled stub (``result = 1.0`` with the norm computation commented out), so
    it always wrote ``100`` and carries no information.  That stub is documented,
    not reproduced.
    """
    dot = 0
    for va, vb in _shared(a, b):
        dot += va * vb
    na = math.sqrt(sum(v * v for v in a.values()))
    nb = math.sqrt(sum(v * v for v in b.values()))
    if na == 0.0 or nb == 0.0:
        return 1.0 if (na == 0.0 and nb == 0.0) else 0.0
    # Two identical vectors can round to 1 + 2**-52; a similarity is never
    # allowed above 1, or its distance would be negative.
    return min(1.0, dot / (na * nb))


MEASURES: Mapping[str, MeasureSpec] = {
    "generalized_jaccard": MeasureSpec(
        name="generalized_jaccard",
        formula="J(a,b) = sum_k min(a_k,b_k) / sum_k max(a_k,b_k)",
        input_domain="integer frequencies",
        empty_input="J = 1.0 when sum max == 0",
        source="legacy/jaccard_parallel_code/src/jaccard_mpi.c; "
        "hsp70_actin/src/c/calc_coeffs.h; hsp70_actin/lib/jacOnlyFromVect.py",
        has_frozen_reference=True,
    ),
    "binary_jaccard": MeasureSpec(
        name="binary_jaccard",
        formula="J(a,b) = |A and B| / |A or B| on presence sets",
        input_domain="presence/absence",
        empty_input="J = 1.0 when the union is empty",
        source="hsp70_actin/src/c/calc_coeffs.h ('normal'); "
        "hsp70_actin/lib/jacOnlyFromVect.py",
        has_frozen_reference=False,
    ),
    "wu": MeasureSpec(
        name="wu",
        formula="J(a,b) = sum_k min(a_k,b_k) / min(sum_k max(a_k,b_k), max(sum a, sum b))",
        input_domain="integer frequencies",
        empty_input="J = 1.0 when the denominator is zero",
        source="hsp70_actin/src/c/calc_coeffs.h; hsp70_actin/lib/jacOnlyFromVect.py",
        has_frozen_reference=False,
    ),
    "sarika": MeasureSpec(
        name="sarika",
        formula="J(a,b) = sum_{k: a_k>0 and b_k>0} (a_k+b_k) / (sum a + sum b)",
        input_domain="integer frequencies",
        empty_input="J = 1.0 when both vectors are empty",
        source="hsp70_actin/src/c/calc_coeffs.h; hsp70_actin/lib/jacOnlyFromVect.py",
        has_frozen_reference=False,
    ),
    "cosine": MeasureSpec(
        name="cosine",
        formula="J(a,b) = a.b / (|a| |b|)",
        input_domain="integer frequencies",
        empty_input="J = 1.0 when both norms are zero, else 0.0",
        source="hsp70_actin/lib/jacOnlyFromVect.py "
        "(the C kernel's cosine is a disabled stub that always writes 100)",
        has_frozen_reference=False,
    ),
}

#: The unchecked kernels.  ``similarity_matrix`` uses these directly: its rows
#: come from ``AnalysisSample``, which already holds validated positive integers.
_FUNCTIONS: Mapping[str, Callable[[Mapping[int, int], Mapping[int, int]], float]] = {
    "generalized_jaccard": _generalized_jaccard,
    "binary_jaccard": _binary_jaccard,
    "wu": _wu,
    "sarika": _sarika,
    "cosine": _cosine,
}


def _public(name: str):
    kernel = _FUNCTIONS[name]

    def measure(a: Mapping[int, Any], b: Mapping[int, Any]) -> float:
        return kernel(_counts(a, "first vector"), _counts(b, "second vector"))

    measure.__name__ = name
    measure.__qualname__ = name
    measure.__doc__ = (
        "%s\n\n``MEASURES[%r].formula``: %s.  A similarity in ``[0, 1]``; "
        "counts are validated (finite, non-negative) and explicit zeros dropped; "
        "two empty vectors give 1.0." % (kernel.__doc__.strip().splitlines()[0],
                                         name, MEASURES[name].formula)
    )
    return measure


generalized_jaccard = _public("generalized_jaccard")
binary_jaccard = _public("binary_jaccard")
wu = _public("wu")
sarika = _public("sarika")
cosine = _public("cosine")


def _condensed_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return n * i - (i * (i + 1)) // 2 + (j - i - 1)


@dataclass(frozen=True)
class SimilarityMatrix:
    """A symmetric all-against-all matrix that knows whether it is a distance.

    Stored condensed (upper triangle, row-major) -- the same layout
    ``jaccard_mpi.c`` writes -- because the diagonal is a definition, not a
    computed value: ``scale`` for a similarity and ``0`` for a distance.

    ``scale`` is ``1.0`` (values are fractions in ``[0, 1]``) unless the matrix
    was made with :meth:`scaled`; ``scaled(100)`` gives percent, and the
    conversions then follow it: ``D = scale - J``.
    """

    sample_ids: Tuple[str, ...]
    measure: str
    kind: str
    values: Tuple[float, ...]
    key_space_fingerprint: str
    comparison: Mapping[str, Any]
    feature_ordering: str
    scale: float = 1.0
    #: Carried from the feature matrix, plus any warning the measure adds.
    warnings: Tuple[str, ...] = ()
    #: Per sample, in row order (``sample_provenance``).
    samples: Tuple[Mapping[str, Any], ...] = ()
    labels: Tuple[Optional[str], ...] = ()

    def __post_init__(self) -> None:
        n = len(self.sample_ids)
        expected = n * (n - 1) // 2
        if len(self.values) != expected:
            raise ValueError(
                "condensed length %d does not match %d samples (expected %d)"
                % (len(self.values), n, expected)
            )
        if self.kind not in (SIMILARITY, DISTANCE):
            raise ValueError("kind must be %r or %r" % (SIMILARITY, DISTANCE))
        if len(set(self.sample_ids)) != n:
            raise ValueError("sample ids must be unique")
        if (isinstance(self.scale, bool) or not isinstance(self.scale, (int, float))
                or not 0 < self.scale < math.inf):
            raise ValueError("scale must be a positive finite number, got %r"
                             % (self.scale,))
        for v in self.values:
            if not (0.0 <= v <= self.scale):
                raise ValueError(
                    "%s %s value %r is outside [0, %g]: not a %s on this scale"
                    % (self.measure, self.kind, v, self.scale, self.kind)
                )
        if self.labels and len(self.labels) != n:
            raise ValueError("%d labels for %d samples" % (len(self.labels), n))

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def diagonal(self) -> float:
        return float(self.scale) if self.kind == SIMILARITY else 0.0

    @property
    def units(self) -> str:
        """``"fraction"`` (0 to 1), ``"percent"`` (scale 100), or ``"x<scale>"``."""
        if self.scale == 1.0:
            return "fraction"
        if self.scale == 100.0:
            return "percent"
        return "x%g" % self.scale

    @property
    def value_label(self) -> str:
        """What one cell is, in words: ``"generalized Jaccard distance (0-1)"``."""
        return "%s %s (%s)" % (
            self.measure.replace("_", " ").replace("jaccard", "Jaccard"),
            self.kind,
            "0-1" if self.scale == 1.0
            else "percent, 0-100" if self.scale == 100.0
            else "0-%g" % self.scale,
        )

    def value(self, i: int, j: int) -> float:
        if i == j:
            return self.diagonal
        return self.values[_condensed_index(self.n_samples, i, j)]

    def square(self) -> Tuple[Tuple[float, ...], ...]:
        n = self.n_samples
        return tuple(tuple(self.value(i, j) for j in range(n)) for i in range(n))

    def condensed(self) -> Tuple[float, ...]:
        return self.values

    def _with(self, kind: str, values, scale: float) -> "SimilarityMatrix":
        return SimilarityMatrix(
            sample_ids=self.sample_ids,
            measure=self.measure,
            kind=kind,
            values=tuple(values),
            key_space_fingerprint=self.key_space_fingerprint,
            comparison=self.comparison,
            feature_ordering=self.feature_ordering,
            scale=scale,
            warnings=self.warnings,
            samples=self.samples,
            labels=self.labels,
        )

    def to_distance(self) -> "SimilarityMatrix":
        """Explicit ``D = scale - J`` (``1 - J`` unscaled).  Never implicit."""
        if self.kind == DISTANCE:
            return self
        return self._with(DISTANCE, (self.scale - v for v in self.values), self.scale)

    def to_similarity(self) -> "SimilarityMatrix":
        """Explicit ``J = scale - D``."""
        if self.kind == SIMILARITY:
            return self
        return self._with(SIMILARITY, (self.scale - v for v in self.values), self.scale)

    def scaled(self, factor: float) -> "SimilarityMatrix":
        """Every value, the diagonal and the scale multiplied by ``factor``.

        ``scaled(100)`` is the percent form the historical clustermap plotted.
        The result records ``scale``, so a percent matrix is never mistaken for
        a fraction: its diagonal, its conversions and its exports all say so.
        """
        if (isinstance(factor, bool) or not isinstance(factor, (int, float))
                or not 0 < factor < math.inf):
            raise ValueError("factor must be a positive finite number, got %r"
                             % (factor,))
        return self._with(self.kind, (factor * v for v in self.values),
                          self.scale * factor)

    def unscaled(self) -> "SimilarityMatrix":
        """Back to fractions in ``[0, 1]``."""
        if self.scale == 1.0:
            return self
        return self._with(self.kind, (v / self.scale for v in self.values), 1.0)

    @property
    def corner_label(self) -> str:
        """The first header cell of :meth:`to_csv`: ``generalized_jaccard_distance``."""
        return "%s_%s%s" % (self.measure, self.kind,
                            "" if self.scale == 1.0 else "_x%g" % self.scale)

    def to_csv(self) -> str:
        """The square matrix, labelled on both axes.

        The first header cell names what the cells are
        (``generalized_jaccard_distance``, ``..._x100`` when scaled), then one
        column per sample id; each row starts with its sample id.  Values are
        written with ``repr``, so they read back exactly.
        """
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow([self.corner_label] + list(self.sample_ids))
        for i, name in enumerate(self.sample_ids):
            writer.writerow([name] + [repr(float(v)) for v in
                                      (self.value(i, j) for j in range(self.n_samples))])
        return out.getvalue()

    def to_long_csv(self) -> str:
        """One row per unordered pair: ``sample_id_1,sample_id_2,measure,kind,scale,value``."""
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["sample_id_1", "sample_id_2", "measure", "kind", "scale", "value"])
        n = self.n_samples
        for i in range(n):
            for j in range(i + 1, n):
                writer.writerow([self.sample_ids[i], self.sample_ids[j], self.measure,
                                 self.kind, repr(float(self.scale)),
                                 repr(float(self.value(i, j)))])
        return out.getvalue()

    def provenance(self) -> Mapping[str, Any]:
        spec = MEASURES.get(self.measure)
        return {
            "measure": self.measure,
            "kind": self.kind,
            "formula": spec.formula if spec else None,
            "empty_input": spec.empty_input if spec else None,
            "conversion": "distance = %g - similarity" % self.scale,
            "scale": self.scale,
            "units": self.units,
            "value_label": self.value_label,
            "n_samples": self.n_samples,
            "sample_ids": list(self.sample_ids),
            "labels": list(self.labels) if self.labels else None,
            "diagonal": self.diagonal,
            "key_space_fingerprint": self.key_space_fingerprint,
            "feature_ordering": self.feature_ordering,
            "comparison": dict(self.comparison),
            "samples": [dict(p) for p in self.samples],
            "warnings": list(self.warnings),
        }


def similarity_matrix(matrix: FeatureMatrix, measure: str = "generalized_jaccard") -> SimilarityMatrix:
    """All-against-all similarity over the rows of ``matrix``, as fractions in ``[0, 1]``.

    Time and memory are quadratic in the number of samples: ``n (n - 1) / 2``
    pair evaluations, each linear in the two samples' distinct keys, and a
    condensed result of that many floats.  It runs in one process.
    """
    try:
        fn = _FUNCTIONS[measure]
    except KeyError:
        raise ValueError(
            "unknown measure %r; evidenced measures are %r"
            % (measure, tuple(sorted(_FUNCTIONS)))
        ) from None
    rows = matrix.rows
    n = len(rows)
    values = []
    for i in range(n):
        a = rows[i]
        for j in range(i + 1, n):
            values.append(fn(a, rows[j]))
    warnings = tuple(matrix.warnings)
    empty = [matrix.sample_ids[i] for i in range(n) if not rows[i]]
    if empty:
        warnings += (
            "%d sample(s) have no keys (%s). By definition two empty samples have "
            "similarity 1 (distance 0) and an empty sample against any other has "
            "similarity 0 (distance 1), so empty samples cluster together."
            % (len(empty), ", ".join(empty[:10]) + (" ..." if len(empty) > 10 else "")),
        )
    return SimilarityMatrix(
        sample_ids=matrix.sample_ids,
        measure=measure,
        kind=SIMILARITY,
        values=tuple(values),
        key_space_fingerprint=matrix.key_space.fingerprint,
        comparison=matrix.comparison,
        feature_ordering=matrix.ordering,
        warnings=warnings,
        samples=tuple(getattr(matrix, "samples", ()) or ()),
        labels=tuple(matrix.labels) if any(l is not None for l in matrix.labels) else (),
    )


def distance_matrix(matrix: FeatureMatrix, measure: str = "generalized_jaccard") -> SimilarityMatrix:
    """All-against-all ``1 - J``; the form every evidenced clustering path consumes."""
    return similarity_matrix(matrix, measure).to_distance()
