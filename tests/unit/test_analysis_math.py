"""Jaccard, distances, scale and clustering, checked against hand calculations.

Every expected number here is worked out by hand in the comment beside it, not
taken from the code under test.
"""

from __future__ import annotations

import csv
import io
import math
from fractions import Fraction

import pytest

from biotsr.analysis import (
    binary_jaccard,
    cluster,
    clustermap_spec,
    distance_matrix,
    feature_matrix,
    generalized_jaccard,
    similarity_matrix,
)
from biotsr.analysis.similarity import DISTANCE, SIMILARITY, SimilarityMatrix

from .test_analysis_features import sample


# ---------------------------------------------------------------------------
# The two Jaccard measures on hand-calculable vectors
# ---------------------------------------------------------------------------

A = {1: 2, 2: 1}
B = {1: 1, 3: 3}


def test_generalized_jaccard_by_hand():
    # keys 1,2,3:  min = (1, 0, 0) -> 1;  max = (2, 1, 3) -> 6;  J = 1/6
    assert generalized_jaccard(A, B) == pytest.approx(1 / 6)


def test_binary_jaccard_by_hand():
    # presence {1,2} and {1,3}: intersection {1}, union {1,2,3}: J = 1/3
    assert binary_jaccard(A, B) == pytest.approx(1 / 3)


def test_frequencies_matter_to_generalized_and_not_to_binary():
    # {1:1} vs {1:5}: binary sees the same set (J = 1); generalized 1/5
    assert binary_jaccard({1: 1}, {1: 5}) == 1.0
    assert generalized_jaccard({1: 1}, {1: 5}) == pytest.approx(0.2)


@pytest.mark.parametrize("measure", [generalized_jaccard, binary_jaccard])
def test_identical_vectors_have_similarity_one(measure):
    assert measure(A, dict(A)) == 1.0


@pytest.mark.parametrize("measure", [generalized_jaccard, binary_jaccard])
def test_disjoint_vectors_have_similarity_zero(measure):
    assert measure({1: 4}, {2: 4}) == 0.0


@pytest.mark.parametrize("measure", [generalized_jaccard, binary_jaccard])
def test_two_empty_vectors_are_identical_and_one_empty_vector_is_disjoint(measure):
    assert measure({}, {}) == 1.0
    assert measure({}, A) == 0.0
    assert measure(A, {}) == 0.0


@pytest.mark.parametrize("measure", [generalized_jaccard, binary_jaccard])
def test_an_explicit_zero_is_an_absent_key(measure):
    # {1:2, 2:0} is {1:2}; {3:0} is empty
    assert measure({1: 2, 2: 0}, {1: 2}) == 1.0
    assert measure({3: 0}, {}) == 1.0
    assert measure({3: 0}, {3: 0, 4: 0}) == 1.0
    assert measure({1: 1, 2: 0}, {2: 5}) == 0.0


@pytest.mark.parametrize("measure", [generalized_jaccard, binary_jaccard])
@pytest.mark.parametrize("bad", [-1, -0.5, math.nan, math.inf, -math.inf, True, "3", None])
def test_invalid_counts_are_refused(measure, bad):
    with pytest.raises(ValueError, match="key 7"):
        measure({7: bad}, {1: 1})
    with pytest.raises(ValueError, match="key 7"):
        measure({1: 1}, {7: bad})


@pytest.mark.parametrize("a, b", [
    ({0: 0.22, 1: 0.42, 2: 0.03}, {0: 0.22, 1: 0.42, 2: 0.03}),
    ({0: 0.8, 1: 0.4, 2: 0.2}, {0: 0.5, 1: 0.7, 2: 0.7}),
])
def test_float_counts_never_give_a_similarity_above_one(a, b):
    """Rounding of float sums can land a hair above 1; it is capped there."""
    from biotsr.analysis import sarika, wu

    for measure in (generalized_jaccard, sarika, wu):
        assert 0.0 <= measure(a, b) <= 1.0


def test_a_sample_refuses_negative_and_non_integer_counts():
    with pytest.raises(ValueError, match="negative"):
        sample("a", {1: -1})
    with pytest.raises(TypeError):
        sample("a", {1: 1.5})
    with pytest.raises(TypeError):
        sample("a", {1: True})


@pytest.mark.parametrize("measure", ["generalized_jaccard", "binary_jaccard"])
def test_the_matrix_agrees_with_the_pairwise_function(measure):
    fm = feature_matrix([sample("a", A), sample("b", B), sample("c", {9: 1})])
    sm = similarity_matrix(fm, measure)
    fn = {"generalized_jaccard": generalized_jaccard, "binary_jaccard": binary_jaccard}[measure]
    rows = [A, B, {9: 1}]
    for i in range(3):
        for j in range(3):
            expected = 1.0 if i == j else fn(rows[i], rows[j])
            assert sm.value(i, j) == expected


# ---------------------------------------------------------------------------
# Matrix properties: symmetry, diagonal, range, D = 1 - J
# ---------------------------------------------------------------------------


def _five():
    return feature_matrix([
        sample("s1", {1: 3, 2: 1}), sample("s2", {1: 1, 2: 1, 3: 2}),
        sample("s3", {4: 2}), sample("s4", {1: 3, 2: 1}), sample("s5", {2: 7, 4: 1}),
    ])


@pytest.mark.parametrize("measure", ["generalized_jaccard", "binary_jaccard"])
def test_symmetry_diagonal_range_and_complement(measure):
    sm = similarity_matrix(_five(), measure)
    dm = sm.to_distance()
    sq, dq = sm.square(), dm.square()
    for i in range(5):
        assert sq[i][i] == 1.0 and dq[i][i] == 0.0
        for j in range(5):
            assert sq[i][j] == sq[j][i] and dq[i][j] == dq[j][i]
            assert 0.0 <= sq[i][j] <= 1.0 and 0.0 <= dq[i][j] <= 1.0
            assert dq[i][j] == 1.0 - sq[i][j]
    assert dm.kind == DISTANCE and sm.kind == SIMILARITY


def test_generalized_jaccard_matches_exact_rational_arithmetic():
    """Against fractions, not floats: s1 = {1:3, 2:1}, s2 = {1:1, 2:1, 3:2}:
    min = 1 + 1 + 0 = 2; max = 3 + 1 + 2 = 6; J = 1/3."""
    sm = similarity_matrix(_five())
    assert Fraction(sm.value(0, 1)).limit_denominator(100) == Fraction(1, 3)
    assert sm.value(0, 3) == 1.0  # s1 and s4 are identical
    assert sm.value(0, 2) == 0.0  # s1 and s3 are disjoint


# ---------------------------------------------------------------------------
# Scale: a percent matrix can never pass for a fraction
# ---------------------------------------------------------------------------


def test_percent_is_explicit_and_follows_through_every_conversion():
    dm = distance_matrix(_five())
    pct = dm.scaled(100)
    assert pct.units == "percent" and pct.scale == 100.0
    assert pct.value(0, 1) == pytest.approx(100 * dm.value(0, 1))
    # D% = 100 - J%, never 1 - J%
    assert pct.to_similarity().value(0, 1) == pytest.approx(100 - pct.value(0, 1))
    assert pct.to_similarity().diagonal == 100.0
    assert pct.unscaled().value(0, 1) == pytest.approx(dm.value(0, 1))
    assert pct.corner_label == "generalized_jaccard_distance_x100"
    assert pct.provenance()["conversion"] == "distance = 100 - similarity"


def test_a_value_outside_the_scale_is_refused():
    """A percent value labelled as a fraction is caught at construction."""
    with pytest.raises(ValueError, match="outside"):
        SimilarityMatrix(sample_ids=("a", "b"), measure="generalized_jaccard",
                         kind=DISTANCE, values=(83.4,), key_space_fingerprint="fp",
                         comparison={}, feature_ordering="numeric_ascending")
    with pytest.raises(ValueError, match="outside"):
        SimilarityMatrix(sample_ids=("a", "b"), measure="generalized_jaccard",
                         kind=DISTANCE, values=(-1e-9,), key_space_fingerprint="fp",
                         comparison={}, feature_ordering="numeric_ascending")


def test_an_empty_sample_is_computed_by_the_rule_and_warned_about():
    fm = feature_matrix([sample("a", {}), sample("b", {}), sample("c", {1: 1})])
    sm = similarity_matrix(fm)
    assert sm.value(0, 1) == 1.0 and sm.value(0, 2) == 0.0
    assert any("no keys" in w for w in sm.warnings)
    assert any("no keys" in w for w in sm.to_distance().warnings)


# ---------------------------------------------------------------------------
# Export: sample ids on both axes, kind and scale named
# ---------------------------------------------------------------------------


def test_the_square_csv_reads_back_exactly_with_its_labels():
    dm = distance_matrix(_five())
    rows = list(csv.reader(io.StringIO(dm.to_csv())))
    assert rows[0] == ["generalized_jaccard_distance", "s1", "s2", "s3", "s4", "s5"]
    assert [r[0] for r in rows[1:]] == ["s1", "s2", "s3", "s4", "s5"]
    assert [[float(v) for v in r[1:]] for r in rows[1:]] == [list(r) for r in dm.square()]


def test_the_long_csv_names_measure_kind_and_scale_on_every_row():
    sm = similarity_matrix(_five(), "binary_jaccard").scaled(100)
    rows = list(csv.DictReader(io.StringIO(sm.to_long_csv())))
    assert len(rows) == 10
    assert {(r["measure"], r["kind"], r["scale"]) for r in rows} == {
        ("binary_jaccard", "similarity", "100.0")}


# ---------------------------------------------------------------------------
# Clustering: alignment, ties, determinism
# ---------------------------------------------------------------------------


def _blocks():
    """Two tight groups {p1,p2,p3} and {q1,q2}, built from counts."""
    return [
        sample("q1", {10: 5, 11: 5}, label="Q"), sample("p1", {1: 5, 2: 5, 3: 1}, label="P"),
        sample("p2", {1: 5, 2: 4, 3: 1}, label="P"), sample("q2", {10: 5, 11: 4}, label="Q"),
        sample("p3", {1: 4, 2: 5, 3: 1}, label="P"),
    ]


def test_ids_rows_leaves_clusters_and_labels_stay_aligned():
    samples = _blocks()
    dm = distance_matrix(feature_matrix(samples))
    spec = clustermap_spec(dm, [s.label for s in samples], n_clusters=2)
    a = spec.assignment
    by_id = dict(zip(a.sample_ids, a.labels))
    assert by_id["p1"] == by_id["p2"] == by_id["p3"] != by_id["q1"] == by_id["q2"]
    assert spec.adjusted_rand == 1.0
    # the assignment CSV is in input order and its leaf positions are the plot order
    rows = list(csv.DictReader(io.StringIO(spec.to_assignments_csv())))
    assert [r["sample_id"] for r in rows] == list(dm.sample_ids)
    assert [r["class_label"] for r in rows] == ["Q", "P", "P", "Q", "P"]
    order = sorted(rows, key=lambda r: int(r["leaf_position"]))
    assert tuple(r["sample_id"] for r in order) == spec.ordered_sample_ids()
    assert spec.to_row_names_csv().splitlines()[1:] == [r["sample_id"] for r in order]
    # palette slots follow the rows they describe
    assert spec.row_class_slots() == (1, 0, 0, 1, 0)  # sorted classes: P=0, Q=1
    cluster_slot = dict(zip(a.sample_ids, spec.row_cluster_slots()))
    assert cluster_slot["p1"] == cluster_slot["p3"] != cluster_slot["q1"]


def test_reordering_the_input_does_not_change_the_partition():
    samples = _blocks()
    first = cluster(distance_matrix(feature_matrix(samples)), 2)
    second = cluster(distance_matrix(feature_matrix(list(reversed(samples)))), 2)

    def parts(assignment):
        return {frozenset(v) for v in assignment.clusters().values()}

    assert parts(first) == parts(second)


def test_all_tied_distances_are_deterministic_and_say_how_many_clusters_formed():
    # four samples, every pairwise distance 1.0 (pairwise disjoint)
    samples = [sample("s%d" % i, {i: 1}) for i in range(4)]
    dm = distance_matrix(feature_matrix(samples))
    assert set(dm.values) == {1.0}
    runs = [cluster(dm, 2) for _ in range(3)]
    assert len({(r.labels, r.leaves, r.linkage.rows) for r in runs}) == 1
    formed = runs[0].n_clusters
    assert formed in (1, 2, 3, 4)
    if formed != 2:
        assert any("requested and %d formed" % formed in w for w in runs[0].warnings)


def test_class_labels_are_supplied_never_inferred_from_names():
    samples = [sample("1ABC_A_LYS", {1: 1}), sample("2DEF_A_ARG", {2: 1})]
    dm = distance_matrix(feature_matrix(samples))
    spec = clustermap_spec(dm, n_clusters=2)
    assert spec.has_classes is False and spec.unique_classes == ()
    with pytest.raises(TypeError, match="sequence in the matrix's sample order"):
        clustermap_spec(dm, {"1ABC_A_LYS": "LYS"}, n_clusters=2)


def test_one_to_n_clusters_when_merge_heights_are_distinct():
    """Heights 1/11 (p1-p2; p1-p3 also at 1/11 but p3 then joins the pair at
    their average, 1.5/11), 0.1 (q1-q2), then the two blocks: every t from 1
    to 5 is a cut of its own, so exactly t clusters form, with no warning."""
    samples = _blocks()
    dm = distance_matrix(feature_matrix(samples))
    assert dm.value(1, 2) == dm.value(1, 4) == pytest.approx(1 / 11)
    for t in range(1, 6):
        assignment = cluster(dm, t)
        assert assignment.n_clusters == t
        assert not any("requested" in w for w in assignment.warnings)
