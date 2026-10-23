"""The five evidenced measures: formulas, empty input, symmetry, and the sim/dist split."""

from __future__ import annotations

import math

import pytest

from biotsr.analysis import (
    DISTANCE,
    MEASURES,
    SIMILARITY,
    binary_jaccard,
    cosine,
    distance_matrix,
    feature_matrix,
    generalized_jaccard,
    sarika,
    similarity_matrix,
    wu,
)

from .test_analysis_features import key_space, sample, sampling


# --- formulas, computed by hand ---------------------------------------------


def test_generalized_jaccard_is_sum_min_over_sum_max():
    a, b = {1: 3, 2: 1}, {1: 1, 3: 4}
    # min: key1 -> 1; max: key1 -> 3, key2 -> 1, key3 -> 4  => 1/8
    assert generalized_jaccard(a, b) == pytest.approx(1 / 8)


def test_binary_jaccard_is_intersection_over_union_of_presence():
    a, b = {1: 3, 2: 1}, {1: 1, 3: 4}
    assert binary_jaccard(a, b) == pytest.approx(1 / 3)


def test_binary_jaccard_ignores_how_often_a_key_occurs():
    assert binary_jaccard({1: 1, 2: 1}, {1: 999}) == binary_jaccard({1: 5, 2: 5}, {1: 1})


def test_wu_divides_by_min_of_sum_max_and_the_larger_row_sum():
    a, b = {1: 3, 2: 1}, {1: 1, 3: 4}
    # sum min = 1; sum max = 8; max(sum a, sum b) = max(4, 5) = 5 -> 1/5
    assert wu(a, b) == pytest.approx(1 / 5)


def test_wu_s_inner_min_never_binds_for_non_negative_counts():
    """``sum max >= max(sum a, sum b)`` always, so wu is ``sum min / max(sum a, sum b)``."""
    cases = [({1: 3, 2: 1}, {1: 1, 3: 4}), ({5: 2}, {5: 2}), ({1: 1}, {2: 1}), ({}, {1: 1})]
    for a, b in cases:
        sa, sb = sum(a.values()), sum(b.values())
        denominator = max(sa, sb)
        smin = sum(min(a.get(k, 0), b.get(k, 0)) for k in set(a) | set(b))
        expected = (smin / denominator) if denominator else 1.0
        assert wu(a, b) == pytest.approx(expected)


def test_sarika_sums_both_frequencies_over_co_present_keys():
    a, b = {1: 3, 2: 1}, {1: 1, 3: 4}
    # co-present: key1 only -> (3 + 1) / (4 + 5) = 4/9
    assert sarika(a, b) == pytest.approx(4 / 9)


def test_cosine_is_the_normalised_dot_product():
    a, b = {1: 3, 2: 4}, {1: 3, 2: 4}
    assert cosine(a, b) == pytest.approx(1.0)
    assert cosine({1: 1}, {2: 1}) == pytest.approx(0.0)
    assert cosine({1: 1, 2: 1}, {1: 1}) == pytest.approx(1 / math.sqrt(2))


# --- identity, symmetry, range ----------------------------------------------


@pytest.mark.parametrize("measure", sorted(MEASURES))
def test_every_measure_is_symmetric(measure):
    from biotsr.analysis.similarity import _FUNCTIONS

    fn = _FUNCTIONS[measure]
    a, b = {1: 3, 2: 1, 5: 7}, {1: 1, 3: 4, 5: 2}
    assert fn(a, b) == pytest.approx(fn(b, a))


@pytest.mark.parametrize("measure", sorted(MEASURES))
def test_every_measure_calls_a_vector_identical_to_itself(measure):
    from biotsr.analysis.similarity import _FUNCTIONS

    a = {1: 3, 2: 1, 5: 7}
    assert _FUNCTIONS[measure](a, dict(a)) == pytest.approx(1.0)


@pytest.mark.parametrize("measure", sorted(MEASURES))
def test_every_measure_stays_inside_the_unit_interval(measure):
    from biotsr.analysis.similarity import _FUNCTIONS

    fn = _FUNCTIONS[measure]
    for a, b in [({1: 1}, {2: 1}), ({1: 5, 2: 1}, {1: 1}), ({}, {}), ({1: 9}, {})]:
        assert 0.0 <= fn(a, b) <= 1.0


# --- empty input ------------------------------------------------------------


@pytest.mark.parametrize("measure", sorted(MEASURES))
def test_two_empty_vectors_are_called_identical(measure):
    """``jaccard_mpi.c``'s rule: ``J = 1.0`` when the denominator is zero."""
    from biotsr.analysis.similarity import _FUNCTIONS

    assert _FUNCTIONS[measure]({}, {}) == 1.0


def test_one_empty_vector_is_maximally_dissimilar():
    assert generalized_jaccard({}, {1: 5}) == 0.0
    assert binary_jaccard({}, {1: 5}) == 0.0
    assert sarika({}, {1: 5}) == 0.0
    assert cosine({}, {1: 5}) == 0.0


# --- matrices ---------------------------------------------------------------


def _matrix():
    return feature_matrix(
        [sample("a", {1: 3, 2: 1}), sample("b", {1: 1, 3: 4}), sample("c", {1: 3, 2: 1})]
    )


def test_a_similarity_matrix_knows_it_is_a_similarity():
    sm = similarity_matrix(_matrix())
    assert sm.kind == SIMILARITY
    assert sm.diagonal == 1.0
    assert sm.value(0, 2) == pytest.approx(1.0)


def test_a_distance_matrix_knows_it_is_a_distance():
    dm = distance_matrix(_matrix())
    assert dm.kind == DISTANCE
    assert dm.diagonal == 0.0
    assert dm.value(0, 2) == pytest.approx(0.0)


def test_the_conversion_between_them_is_explicit_and_round_trips():
    sm = similarity_matrix(_matrix())
    assert sm.to_distance().to_similarity().values == pytest.approx(sm.values)
    assert sm.to_distance().value(0, 1) == pytest.approx(1 - sm.value(0, 1))
    assert sm.to_similarity() is sm


def test_a_matrix_is_symmetric_and_sample_ids_are_preserved():
    dm = distance_matrix(_matrix())
    square = dm.square()
    assert dm.sample_ids == ("a", "b", "c")
    for i in range(3):
        for j in range(3):
            assert square[i][j] == pytest.approx(square[j][i])


def test_an_unknown_measure_names_the_evidenced_ones():
    with pytest.raises(ValueError, match="evidenced measures"):
        similarity_matrix(_matrix(), "euclidean")


def test_every_measure_records_its_formula_and_whether_a_reference_exists():
    assert MEASURES["generalized_jaccard"].has_frozen_reference is True
    for name in ("binary_jaccard", "wu", "sarika", "cosine"):
        assert MEASURES[name].has_frozen_reference is False
        assert MEASURES[name].formula
        assert MEASURES[name].source


def test_provenance_carries_the_key_space_and_the_sampling_variation():
    dm = distance_matrix(_matrix())
    p = dm.provenance()
    assert p["measure"] == "generalized_jaccard"
    assert p["kind"] == DISTANCE
    assert p["key_space_fingerprint"]
    assert p["comparison"]["processing_uniform"] is True
    assert p["scale"] == 1.0 and p["units"] == "fraction"
    assert p["conversion"] == "distance = 1 - similarity"


# --- the removed legacy exports --------------------------------------------


def test_the_legacy_matrix_csv_and_quantization_are_gone():
    import biotsr.analysis as analysis

    assert not hasattr(analysis, "to_legacy_matrix_csv")
    assert not hasattr(distance_matrix(_matrix()), "quantized")


def test_scaling_multiplies_every_value_and_records_the_scale():
    dm = distance_matrix(_matrix())
    scaled = dm.scaled(100.0)
    assert scaled.values == pytest.approx(tuple(100.0 * v for v in dm.values))
    assert scaled.kind == DISTANCE
    assert (scaled.scale, scaled.units) == (100.0, "percent")
    assert scaled.diagonal == 0.0
    assert scaled.to_similarity().diagonal == 100.0
    assert "percent" in scaled.value_label


def test_a_condensed_matrix_of_the_wrong_length_is_refused():
    from biotsr.analysis.similarity import SimilarityMatrix

    with pytest.raises(ValueError, match="does not match"):
        SimilarityMatrix(
            sample_ids=("a", "b", "c"),
            measure="generalized_jaccard",
            kind=DISTANCE,
            values=(0.1, 0.2),
            key_space_fingerprint="x",
            comparison={},
            feature_ordering="numeric_ascending",
        )
