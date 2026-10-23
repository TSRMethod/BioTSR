"""The Ca distance top class end to end, and the standard theta bounds.

Two gaps these tests close:

1. The Ca distance scheme was unit-tested, but **no real acceptance chain
   reaches 132.45 A**, so the standard 35th distance class was never exercised
   end to end.  A synthetic triangle closes that -- as a regression test, never
   a production special case.
2. Standard-profile theta bounds were implied but not asserted as invariants.
"""

from __future__ import annotations

import math

import pytest

import biotsr
from biotsr.candidates.exhaustive import exhaustive
from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
from biotsr.domain.point import LabeledPoint, Point, PointSet
from biotsr.encoding.key import encode
from biotsr.geometry.ordering import canonical_order
from biotsr.geometry.triangle import triangle_geometry
from biotsr.methods.engine import run_core

# ---------------------------------------------------------------------------
# 1. End-to-end: the 35th distance class above 132.45 A
# ---------------------------------------------------------------------------


def _point(seq, residue, coord):
    return LabeledPoint(
        identity=AtomIdentity(
            structure=StructureId("SYNTH"),
            model=1,
            chain=ChainId("A"),
            residue_name=residue,
            residue_seq=seq,
            atom_name="CA",
        ),
        coord=coord,
        label=biotsr.STANDARD_PROTEIN_CALPHA_V1.label_scheme.table[residue],
    )


@pytest.fixture(scope="module")
def long_edge_points():
    """Three Ca points whose longest edge is ~133 A -- past the 132.45 A cut.

    Distinct labels put this in canonical-ordering case 4, so the ordering is
    label-driven and the base edge is NOT the longest one. The 133 A edge is
    therefore the distance-bin input while playing no part in theta -- exactly the
    configuration that makes this a real end-to-end exercise rather than a
    restatement of the unit test.
    """
    return PointSet(
        [
            _point(1, "MET", Point(0.0, 0.0, 0.0)),
            _point(2, "ILE", Point(133.0, 0.0, 0.0)),
            _point(3, "GLN", Point(60.0, 40.0, 0.0)),
        ]
    )


def test_the_synthetic_triangle_really_exceeds_the_top_cut(long_edge_points):
    coords = long_edge_points.coords
    longest = max(
        coords[0].distance_to(coords[1]),
        coords[0].distance_to(coords[2]),
        coords[1].distance_to(coords[2]),
    )
    assert longest == pytest.approx(133.0)
    assert longest > 132.45


def _single_triangle(points, profile):
    candidate = next(exhaustive(len(points)))
    idx = candidate.as_tuple()
    labels = tuple(points.labels[i] for i in idx)
    pts = tuple(points.coords[i] for i in idx)
    order = canonical_order(labels, pts)
    o0, o1, o2 = order.order
    geometry = triangle_geometry(
        pts[o0], pts[o1], pts[o2], profile.angle_convention
    )
    theta_class = profile.theta_scheme.classify(geometry.theta)
    distance_class = profile.distance_scheme.classify(geometry.max_distance)
    ordered_labels = (labels[o0], labels[o1], labels[o2])
    key = encode(ordered_labels, theta_class, distance_class, profile.key_formula())
    return geometry, theta_class, distance_class, key, ordered_labels


def test_end_to_end_the_long_triangle_reaches_the_35th_class(long_edge_points):
    """candidate -> ordering -> geometry -> binning -> encoding, standard profile."""
    s_geom, s_tclass, s_dclass, s_key, s_labels = _single_triangle(
        long_edge_points, biotsr.STANDARD_PROTEIN_CALPHA_V1
    )
    assert s_geom.valid
    assert s_dclass == 35, "the standard scheme has a real 35th class"
    assert s_geom.max_distance == pytest.approx(133.0)


def test_the_top_class_lands_exactly_in_the_distance_component(long_edge_points):
    """The distance digit's positional weight is the number of theta classes.

    Re-encoding the same triangle with the distance class one lower (34, the
    class below the top cut) must change the key by exactly ``n_theta`` (29),
    so the 35th class is a genuine, separately encoded key component.
    """
    s_geom, s_tclass, s_dclass, s_key, s_labels = _single_triangle(
        long_edge_points, biotsr.STANDARD_PROTEIN_CALPHA_V1
    )
    spec = biotsr.STANDARD_PROTEIN_CALPHA_V1.key_formula()
    assert encode(s_labels, s_tclass, s_dclass, spec) == s_key
    one_lower = encode(s_labels, s_tclass, s_dclass - 1, spec)
    n_theta = biotsr.STANDARD_PROTEIN_CALPHA_V1.theta_scheme.n_classes
    assert s_key - one_lower == n_theta * 1 == 29


def test_the_top_class_survives_the_whole_engine(long_edge_points):
    standard = run_core(long_edge_points, biotsr.STANDARD_PROTEIN_CALPHA_V1)
    _, _, _, s_key, _ = _single_triangle(
        long_edge_points, biotsr.STANDARD_PROTEIN_CALPHA_V1
    )
    assert standard.counts.n_triangles == 1
    assert standard.triplets[0].distance_class == 35
    assert standard.triplets[0].key == s_key


def test_a_size_filter_can_exclude_the_long_triangle(long_edge_points):
    """The two mechanisms are independent: the 35th class exists, and a caller may
    still choose not to look at triangles that reach it."""
    profile = biotsr.STANDARD_PROTEIN_CALPHA_V1.with_modifiers(size_filter=100.0)
    out = run_core(long_edge_points, profile)
    assert out.counts.n_candidates == 1
    assert out.counts.n_valid_geometry == 1
    assert out.counts.n_filtered_by_size == 1
    assert out.counts.n_degenerate == 0, "filtered is not degenerate"
    assert out.counts.n_triangles == 0
    assert out.counts.ledger_balances


# ---------------------------------------------------------------------------
# 2. Standard theta invariants
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ubq_standard(ubq_pdb):
    return biotsr.analyze(ubq_pdb, method="protein", chain="A")


def test_standard_theta_raw_stays_within_zero_to_180(ubq_standard):
    for t in ubq_standard.triplets:
        assert 0.0 <= t.theta_raw <= 180.0


def test_standard_folded_theta_stays_within_zero_to_90(ubq_standard):
    for t in ubq_standard.triplets:
        assert 0.0 <= t.theta <= 90.0


def test_the_standard_bound_is_exact_not_approximate():
    """math.pi makes the pre-fold range exactly [0, 180]; the divisor is the reason."""
    convention = biotsr.STANDARD_PROTEIN_CALPHA_V1.angle_convention
    assert convention.degrees(0.0) == 0.0
    assert convention.degrees(math.pi) == pytest.approx(180.0, abs=1e-12)
    assert convention.degrees(math.pi) <= 180.0 + 1e-12


def test_the_fold_rule_holds_under_the_standard_profile(ubq_standard):
    for t in ubq_standard.triplets[:2000]:
        expected = t.theta_raw if t.theta_raw <= 90 else abs(180 - t.theta_raw)
        assert t.theta == pytest.approx(expected, abs=1e-12)
