"""The shared TSR core: point set + profile -> structured result.

This is the one place the stages are wired together, and it contains **no
science of its own** -- no ordering rules, no geometry, no bin tables, no key
arithmetic, no modifier logic.  Every scientific decision lives in the component
it names.  That is what lets a future method reuse this loop unchanged by
supplying a different selector, label scheme, candidate strategy or bin scheme.

**Modifiers do not branch here.**  The loop never asks "is grouping on?" or "is
mirror on?": the profile has already resolved the label scheme, the size filter
and the mirror convention (``profiles/standard.py::Profile.with_modifiers``), and
this loop simply uses them.  That is what makes all eight combinations one code
path, and why no ``GroupedMirrorTSR``-style class exists or is needed.

**SSE does not branch here either, and does not fork the loop.**  It adds
one optional collaborator, ``sse_annotator``, consulted at the ANNOTATION hook
*after* the key is complete.  Everything above that point -- selection, candidate
enumeration, ordering, geometry, the size filter, binning, encoding and the
mirror sign -- is the same code executing on the same values whether SSE is on or
off, which is why the ordinary key is invariant to annotation by construction
rather than by agreement.  There is no ``SSETSR`` class, no second engine, and no
duplicated arithmetic anywhere in ``biotsr.sse``.
"""

from __future__ import annotations

from collections import Counter
from typing import Callable, Dict, Iterable, List, Optional, Tuple

from ..candidates.exhaustive import exhaustive
from ..domain.candidate import TriangleCandidate
from ..domain.point import PointSet
from ..domain.result import (
    CategoryCounts,
    ResultCounts,
    SkipRecord,
    TripletRecord,
    TripletSSE,
)
from ..encoding.key import encode
from ..geometry.ordering import canonical_order
from ..geometry.triangle import triangle_geometry
from ..profiles.standard import Profile

__all__ = ["CoreOutput", "run_core"]


class CoreOutput:
    """Everything the core produced, before provenance is attached."""

    __slots__ = (
        "key_counts",
        "key_order",
        "triplets",
        "skipped",
        "counts",
    )

    def __init__(
        self,
        key_counts: Dict[int, "Counter[str]"],
        key_order: Tuple[int, ...],
        triplets: Optional[List[TripletRecord]],
        skipped: List[SkipRecord],
        counts: ResultCounts,
    ):
        self.key_counts = key_counts
        self.key_order = key_order
        self.triplets = triplets
        self.skipped = skipped
        self.counts = counts

    def as_category_counts(self) -> Dict[int, CategoryCounts]:
        return {
            key: CategoryCounts(by_category=dict(counter))
            for key, counter in self.key_counts.items()
        }


def run_core(
    points: PointSet,
    profile: Profile,
    *,
    candidates: Optional[Iterable[TriangleCandidate]] = None,
    candidate_count: Optional[int] = None,
    keep_triplets: bool = True,
    triplet_sink: Optional[Callable[[TripletRecord], None]] = None,
    max_degenerate_reported: int = 1000,
    sse_annotator=None,
) -> CoreOutput:
    """Run canonical ordering -> geometry -> binning -> encoding over candidates.

    ``key_order`` records the order in which distinct keys were first seen.  It
    is kept on the in-memory result for callers who want it; the native
    serializer ignores it and sorts numerically instead.

    ``sse_annotator`` is the ANNOTATION hook.  When it is ``None`` the category is
    the mirror category, exactly as before.  When it is supplied, the category is
    the triangle's SSE category instead -- the key, including its mirror sign, is
    already fully determined by the time it is consulted, and the annotator is
    handed identities only, so it *cannot* influence one.  Under mirror the two do
    not compete for the bucket name: the mirror state is carried by the key's
    sign, and the bucket names the secondary structure.

    ``triplet_sink`` lets a serializer consume records as they are produced.  It
    is a consumption mode, not a second result model: the record handed to the
    sink is the same ``TripletRecord`` that would otherwise be retained, and the
    counts are identical either way.  With ``keep_triplets=False`` and a sink,
    memory stays flat over runs producing millions of triangles -- which the
    larger acceptance fixtures do.
    """
    coords = points.coords
    labels = points.labels
    identities = points.identities
    point_sides = tuple(p.side for p in points)

    formula = profile.key_formula()
    theta_scheme = profile.theta_scheme
    distance_scheme = profile.distance_scheme
    convention = profile.angle_convention
    sqrt = profile.sqrt_convention
    median = profile.median_convention
    ordering_convention = profile.ordering_convention

    if candidates is None:
        candidates = exhaustive(len(points))
    if candidate_count is None:
        n = len(points)
        candidate_count = 0 if n < 3 else n * (n - 1) * (n - 2) // 6

    key_counts: Dict[int, "Counter[str]"] = {}
    key_order: List[int] = []
    triplets: Optional[List[TripletRecord]] = [] if keep_triplets else None
    build_record = keep_triplets or triplet_sink is not None
    skipped: List[SkipRecord] = []

    size_filter = profile.size_filter
    mirror = profile.mirror_convention
    # The scheme's own name for each code. Precomputed once: under AA grouping the
    # reverse lookup is many-to-one, and doing it per triangle would be needless
    # work on a multi-million-candidate run.
    code_names = {}
    # The first spelling of each code wins: the residue name for the Ca
    # schemes, the element symbol for the atomic scheme (hydrogen before
    # deuterium, which share 107 and never label a point).
    for entity_name, code in profile.label_scheme.table.items():
        code_names.setdefault(code, entity_name)

    n_triangles = 0
    n_valid_geometry = 0
    n_degenerate = 0
    n_clamped = 0
    n_filtered_by_size = 0

    for candidate in candidates:
        idx = candidate.vertices()
        triple_labels = (labels[idx[0]], labels[idx[1]], labels[idx[2]])
        triple_points = (coords[idx[0]], coords[idx[1]], coords[idx[2]])

        ordering = canonical_order(
            triple_labels, triple_points, sqrt, ordering_convention
        )
        o0, o1, o2 = ordering.order
        p0, p1, p2 = triple_points[o0], triple_points[o1], triple_points[o2]

        geometry = triangle_geometry(
            p0,
            p1,
            p2,
            convention,
            sqrt=sqrt,
            clamp_acos=profile.clamp_acos,
            acos_tolerance=profile.acos_tolerance,
            median=median,
        )

        if not geometry.valid:
            n_degenerate += 1
            if len(skipped) < max_degenerate_reported:
                skipped.append(
                    SkipRecord(
                        stage="geometry",
                        reason=geometry.degeneracy_reason or "unknown",
                        detail="; ".join(
                            identities[idx[p]].describe() for p in (o0, o1, o2)
                        ),
                        indices=(idx[o0], idx[o1], idx[o2]),
                    )
                )
            continue

        n_valid_geometry += 1
        if geometry.acos_clamped:
            n_clamped += 1

        # CANDIDATE_FILTER hook. Valid geometry that is simply outside the
        # caller's cutoff -- a different outcome from degeneracy, counted
        # separately, and rejected here so nothing downstream does needless work.
        if size_filter is not None and not size_filter.accepts(geometry.max_distance):
            n_filtered_by_size += 1
            continue

        theta_class = theta_scheme.classify(geometry.theta)
        distance_class = distance_scheme.classify(geometry.max_distance)
        ordered_labels = (
            triple_labels[o0],
            triple_labels[o1],
            triple_labels[o2],
        )
        key = encode(ordered_labels, theta_class, distance_class, formula)

        # KEY_SIGN hook. The magnitude is already fully determined; mirror only
        # decides the sign, off the PRE-FOLD angle, and names the category.
        category = mirror.category_for(geometry.theta_raw)
        key = mirror.apply(key, geometry.theta_raw)

        # ANNOTATION hook. The key is final -- sign and all -- before this runs,
        # and the annotator sees identities only. It names the frequency bucket;
        # it can move nothing.
        sse_record = None
        if sse_annotator is not None:
            ordered_identities = (
                identities[idx[o0]],
                identities[idx[o1]],
                identities[idx[o2]],
            )
            triple, category = sse_annotator.annotate(ordered_identities)
            sse_record = TripletSSE(
                states=(triple[0].state, triple[1].state, triple[2].state),
                # `.text` is the display spelling of a structured element
                # identity; the identity itself never leaves `biotsr.sse`, where
                # it is compared componentwise (ADR-063).
                groups=(
                    triple[0].group.text,
                    triple[1].group.text,
                    triple[2].group.text,
                ),
                category=category,
            )

        counter = key_counts.get(key)
        if counter is None:
            counter = Counter()
            key_counts[key] = counter
            key_order.append(key)
        counter[category] += 1
        n_triangles += 1

        if build_record:
            record = TripletRecord(
                key=key,
                labels=ordered_labels,
                identities=(
                    identities[idx[o0]],
                    identities[idx[o1]],
                    identities[idx[o2]],
                ),
                label_names=(
                    code_names.get(ordered_labels[0], ""),
                    code_names.get(ordered_labels[1], ""),
                    code_names.get(ordered_labels[2], ""),
                ),
                theta_raw=geometry.theta_raw,
                theta=geometry.theta,
                theta_class=theta_class,
                distance_class=distance_class,
                max_distance=geometry.max_distance,
                coordinates=(
                    p0.as_tuple(),
                    p1.as_tuple(),
                    p2.as_tuple(),
                ),
                category=category,
                sides=(
                    point_sides[idx[o0]],
                    point_sides[idx[o1]],
                    point_sides[idx[o2]],
                ),
                sse=sse_record,
                median_length=geometry.median_length,
            )
            if triplets is not None:
                triplets.append(record)
            if triplet_sink is not None:
                triplet_sink(record)

    if n_degenerate > max_degenerate_reported:
        skipped.append(
            SkipRecord(
                stage="geometry",
                reason="degenerate_records_truncated",
                detail="%d degenerate candidates total; only the first %d are listed"
                % (n_degenerate, max_degenerate_reported),
            )
        )

    counts = ResultCounts(
        n_points=len(points),
        n_candidates=candidate_count,
        n_triangles=n_triangles,
        n_unique_keys=len(key_counts),
        n_valid_geometry=n_valid_geometry,
        n_degenerate=n_degenerate,
        n_acos_clamped=n_clamped,
        n_filtered_by_size=n_filtered_by_size,
    )
    return CoreOutput(key_counts, tuple(key_order), triplets, skipped, counts)
