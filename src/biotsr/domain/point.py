"""Geometric point value objects and the labelled point set the engine consumes.

``PointSet`` is the boundary between selection and the mathematics: everything
downstream of it sees coordinates, integer labels and identities -- never a PDB
line.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence, Tuple

from .identity import AtomIdentity

__all__ = ["Point", "LabeledPoint", "PointSet"]


@dataclass(frozen=True)
class Point:
    """An immutable 3-D coordinate in angstroms."""

    x: float
    y: float
    z: float

    def squared_distance_to(self, other: "Point") -> float:
        """The sum of squared differences -- the square root is the caller's choice.

        Separated because the legacy family is not unanimous about how it takes
        that root, and the difference is observable: ``TSR-Package`` (the Ca
        lineage) writes ``... ** 0.5`` while both atomic lineages write
        ``math.sqrt(...)``.  Those are not the same function -- ``** 0.5`` goes
        through ``pow``, which is not required to be correctly rounded, while
        ``math.sqrt`` maps to the correctly-rounded hardware instruction -- and
        they disagree by one ULP on roughly 0.12% of values.  One ULP in a
        distance is enough to move a triangle across a bin boundary, so the
        choice is a key-affecting convention, not an implementation detail.

        The expression order matches the legacy ``calcDist``
        (``(x1-x2)**2 + (y2-y1)**2 + (z2-z1)**2``).  Squaring makes the sign
        irrelevant and IEEE-754 negation is exact, so the mixed orderings across
        the legacy files are bit-identical; the legacy spelling is kept so the
        equivalence is visible rather than assumed.
        """
        return (
            (self.x - other.x) ** 2
            + (other.y - self.y) ** 2
            + (other.z - self.z) ** 2
        )

    def distance_to(self, other: "Point") -> float:
        """Euclidean distance using the ``** 0.5`` spelling (the Ca convention).

        Geometry and ordering take an explicit square-root function instead of
        calling this, so a family can pin its own spelling; this stays as the
        readable default for callers that do not care.
        """
        return self.squared_distance_to(other) ** 0.5

    def midpoint_to(self, other: "Point") -> "Point":
        return Point(
            (self.x + other.x) / 2,
            (self.y + other.y) / 2,
            (self.z + other.z) / 2,
        )

    def as_tuple(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass(frozen=True)
class LabeledPoint:
    """One participating point: full identity, coordinate, and integer label.

    ``side`` tags which pool a point came from.  A one-sided method uses the single
    default pool; a cross method builds mixed pools through the same field,
    without changing this type or anything downstream of it.
    """

    identity: AtomIdentity
    coord: Point
    label: int
    side: str = "primary"

    def __post_init__(self) -> None:
        if not isinstance(self.label, int):
            raise ValueError("label must be an int")


class PointSet(Sequence[LabeledPoint]):
    """An ordered, immutable sequence of labelled points.

    Order is significant and deterministic: it is the enumeration order the
    candidate strategy sees, and therefore the order that fixes triplet
    serialization order.
    """

    __slots__ = ("_points",)

    def __init__(self, points: Sequence[LabeledPoint]):
        self._points: Tuple[LabeledPoint, ...] = tuple(points)

    def __len__(self) -> int:
        return len(self._points)

    def __getitem__(self, index):  # type: ignore[override]
        if isinstance(index, slice):
            return PointSet(self._points[index])
        return self._points[index]

    def __iter__(self) -> Iterator[LabeledPoint]:
        return iter(self._points)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PointSet):
            return NotImplemented
        return self._points == other._points

    def __hash__(self) -> int:
        return hash(self._points)

    def __repr__(self) -> str:  # pragma: no cover - diagnostic
        return "PointSet(n=%d)" % len(self._points)

    @property
    def coords(self) -> Tuple[Point, ...]:
        return tuple(p.coord for p in self._points)

    @property
    def labels(self) -> Tuple[int, ...]:
        return tuple(p.label for p in self._points)

    @property
    def identities(self) -> Tuple[AtomIdentity, ...]:
        return tuple(p.identity for p in self._points)
