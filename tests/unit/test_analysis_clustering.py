"""Clustering: the computation, its provenance, and its independence from any figure."""

from __future__ import annotations

import pytest

from biotsr.analysis import (
    AVERAGE,
    adjusted_rand_index,
    average_linkage,
    cluster,
    clustermap_spec,
    distance_matrix,
    feature_matrix,
    flat_clusters,
    leaf_order,
    similarity_matrix,
)
from biotsr.analysis.similarity import DISTANCE, SimilarityMatrix

from .test_analysis_features import sample


def _dm(values, n):
    return SimilarityMatrix(
        sample_ids=tuple("s%d" % i for i in range(n)),
        measure="generalized_jaccard",
        kind=DISTANCE,
        values=tuple(values),
        key_space_fingerprint="fp",
        comparison={"uniform": True},
        feature_ordering="numeric_ascending",
    )


def test_a_similarity_matrix_is_refused_rather_than_silently_inverted():
    fm = feature_matrix([sample("a", {1: 1}), sample("b", {2: 1})])
    with pytest.raises(ValueError, match="needs a distance matrix"):
        average_linkage(similarity_matrix(fm))


def test_fewer_than_two_samples_cannot_be_clustered():
    with pytest.raises(ValueError, match="at least 2 samples"):
        average_linkage(_dm([], 1))


def test_average_linkage_on_two_obvious_pairs():
    #   s0,s1 close (0.1); s2,s3 close (0.1); the two pairs far apart (0.9)
    values = [0.1, 0.9, 0.9, 0.9, 0.9, 0.1]
    link = average_linkage(_dm(values, 4))
    assert link.method == AVERAGE
    assert len(link.rows) == 3
    assert link.heights[0] == pytest.approx(0.1)
    assert link.heights[1] == pytest.approx(0.1)
    assert link.heights[2] == pytest.approx(0.9)
    assert link.is_monotonic


def test_the_linkage_is_computed_on_the_matrix_as_given_and_records_its_scale():
    """No hidden x100: a percent matrix gives percent heights, and says so."""
    values = [0.1, 0.9, 0.9, 0.9, 0.9, 0.1]
    plain = average_linkage(_dm(values, 4))
    percent = average_linkage(_dm(values, 4).scaled(100.0))
    assert plain.scale == 1.0 and percent.scale == 100.0
    assert percent.heights == pytest.approx(tuple(100.0 * h for h in plain.heights))
    assert max(plain.heights) <= 1.0
    assert flat_clusters(plain, 2) == flat_clusters(percent, 2)
    assert leaf_order(plain) == leaf_order(percent)


def test_flat_clusters_recovers_the_two_pairs():
    link = average_linkage(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4))
    labels = flat_clusters(link, 2)
    assert labels[0] == labels[1]
    assert labels[2] == labels[3]
    assert labels[0] != labels[2]
    assert sorted(set(labels)) == [1, 2]


def test_flat_clusters_labels_are_one_based():
    link = average_linkage(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4))
    assert min(flat_clusters(link, 4)) == 1


def test_tied_heights_yield_fewer_clusters_than_requested():
    """Two merges tie at 0.1, so a threshold either keeps both or neither.

    ``maxclust`` cuts at a *height*: the lowest threshold leaving at most t
    clusters.  For t = 2 and t = 3 that is 0.1 (two clusters; no threshold
    gives three); for t >= 4 it is below every merge, so every sample stands
    alone.
    """
    link = average_linkage(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4))
    assert flat_clusters(link, 1) == (1, 1, 1, 1)
    for t in (2, 3):
        assert flat_clusters(link, t) == (1, 1, 2, 2)
    for t in (4, 99):
        assert sorted(flat_clusters(link, t)) == [1, 2, 3, 4]


def test_maxclust_finds_the_lowest_cut_with_distinct_heights():
    """Merge heights 0.1, 0.2, 0.75: t clusters for every t, as defined.

    scipy 1.11.1's binary search never evaluated the lowest candidate cut and
    returned two clusters for t = 3 and t = 4 here; BioTSR follows the
    definition (as later scipy releases do), so asking for three gives three
    (ADR-060, amended in 1.0.0rc3).
    """
    link = average_linkage(_dm([0.1, 0.8, 0.9, 0.7, 0.6, 0.2], 4))
    assert link.heights == pytest.approx((0.1, 0.2, 0.75))
    assert flat_clusters(link, 1) == (1, 1, 1, 1)
    assert flat_clusters(link, 2) == (1, 1, 2, 2)
    for t in (3, 4):
        assert len(set(flat_clusters(link, t))) == t


def test_a_requested_count_is_reached_when_the_search_can_find_it():
    """The 1BNA reference asks for 4 over 11 samples and gets exactly 4."""
    link = average_linkage(
        _dm([0.05, 0.40, 0.42, 0.44, 0.41, 0.43, 0.45, 0.06, 0.47, 0.48, 0.07,
             0.50, 0.52, 0.54, 0.51], 6)
    )
    assert len(set(flat_clusters(link, 3))) == 3


def test_asking_for_fewer_than_one_cluster_is_refused():
    link = average_linkage(_dm([0.5], 2))
    with pytest.raises(ValueError, match="at least 1"):
        flat_clusters(link, 0)


def test_two_samples_asked_for_one_cluster_get_one_cluster():
    """scipy 1.11.1 answered ``[1, 2]`` here by reading one element past its
    array (ADR-060); scipy 1.18.1 and BioTSR answer ``[1, 1]``."""
    link = average_linkage(_dm([0.5], 2))
    assert flat_clusters(link, 1) == (1, 1)
    assert flat_clusters(link, 2) == (1, 2)


def test_leaf_order_is_a_permutation_of_every_sample():
    link = average_linkage(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4))
    assert sorted(leaf_order(link)) == [0, 1, 2, 3]


def test_leaf_order_keeps_a_tight_pair_adjacent():
    link = average_linkage(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4))
    order = list(leaf_order(link))
    assert abs(order.index(0) - order.index(1)) == 1
    assert abs(order.index(2) - order.index(3)) == 1


def test_a_cluster_assignment_exists_with_no_figure_anywhere_near_it():
    assignment = cluster(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4), 2)
    assert assignment.n_clusters == 2
    assert assignment.ordered_sample_ids()
    assert set(assignment.clusters()) == {1, 2}
    assert len(assignment.clusters()[1]) == 2


def test_provenance_records_everything_a_rerun_needs():
    assignment = cluster(_dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4).scaled(100.0), 2)
    p = assignment.provenance()
    assert p["linkage_method"] == "average"
    assert p["criterion"] == "maxclust"
    assert p["metric"] == "precomputed:generalized_jaccard:distance"
    assert p["scale"] == 100.0
    assert p["n_clusters_requested"] == 2
    assert p["n_clusters_formed"] == 2
    assert p["deterministic"] is True
    assert p["random_seed"] is None
    assert p["key_space_fingerprint"] == "fp"
    assert p["sample_order"] == ["s0", "s1", "s2", "s3"]
    assert p["leaf_order"] == [assignment.sample_ids[i] for i in assignment.leaves]
    assert p["cluster_labels"] == list(assignment.labels)
    assert "never takes it from known class labels" in p["cut_rule_choice"]


def test_clustering_is_deterministic_across_repeated_calls():
    values = [0.1, 0.9, 0.9, 0.9, 0.9, 0.1]
    first = cluster(_dm(values, 4), 2)
    second = cluster(_dm(values, 4), 2)
    assert first.labels == second.labels
    assert first.leaves == second.leaves
    assert first.linkage.rows == second.linkage.rows


# --- adjusted Rand index ----------------------------------------------------


def test_ari_is_one_for_a_perfect_relabelling():
    assert adjusted_rand_index(["a", "a", "b", "b"], [2, 2, 1, 1]) == pytest.approx(1.0)


def test_ari_is_one_when_both_labellings_are_all_singletons():
    assert adjusted_rand_index([0, 1, 2], [5, 6, 7]) == pytest.approx(1.0)


def test_ari_is_one_when_both_labellings_are_a_single_cluster():
    assert adjusted_rand_index([0, 0, 0], [1, 1, 1]) == pytest.approx(1.0)


def test_ari_is_zero_when_one_side_is_uninformative():
    assert adjusted_rand_index(["a", "a", "b", "b"], [1, 1, 1, 1]) == pytest.approx(0.0)


def test_ari_can_be_negative():
    assert adjusted_rand_index(["a", "a", "b", "b"], [1, 2, 1, 2]) < 0.0


def test_ari_refuses_mismatched_lengths():
    with pytest.raises(ValueError, match="different lengths"):
        adjusted_rand_index([1, 2], [1])


# --- the clustermap spec ----------------------------------------------------


def test_the_clustermap_spec_is_data_and_computes_no_number_by_drawing():
    dm = _dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4)
    spec = clustermap_spec(dm, ["x", "x", "y", "y"], n_clusters=2)
    assert spec.assignment.n_clusters == 2
    assert spec.adjusted_rand == pytest.approx(1.0)
    assert spec.unique_classes == ("x", "y")
    assert spec.row_class_slots() == (0, 0, 1, 1)
    assert set(spec.row_cluster_slots()) == {0, 1}
    assert spec.figsize == (10, 8)
    assert spec.provenance()["visualization_affects_results"] is False
    assert spec.provenance()["plotted_matrix"] == "generalized Jaccard distance (0-1)"


def test_the_number_of_clusters_is_never_inferred():
    """Not even from the class count: that would make agreement a construction.

    Without a cutting rule the spec is the tree and the ordered heat map, with
    no partition and no ARI; with one, the rule given is the rule used."""
    dm = _dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4)
    uncut = clustermap_spec(dm, ["x", "x", "y", "y"])
    assert uncut.assignment is None
    assert uncut.adjusted_rand is None
    assert uncut.provenance()["flat_partition"] is None
    assert clustermap_spec(dm, ["x", "x", "y", "y"],
                           n_clusters=3).assignment.n_clusters_requested == 3


def test_an_unlabelled_clustermap_has_no_class_stripe_and_no_ari():
    dm = _dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4)
    spec = clustermap_spec(dm, n_clusters=2)
    assert spec.has_classes is False
    assert spec.adjusted_rand is None
    assert spec.to_assignments_csv().splitlines()[0] == "sample_id,cluster,leaf_position"


def test_the_clustermap_refuses_a_similarity_matrix():
    fm = feature_matrix([sample("a", {1: 1}), sample("b", {2: 1})])
    with pytest.raises(ValueError, match="plots a distance matrix"):
        clustermap_spec(similarity_matrix(fm), ["x", "y"], n_clusters=2)


def test_the_clustermap_refuses_a_label_count_that_does_not_match():
    dm = _dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4)
    with pytest.raises(ValueError, match="3 class labels for 4 samples"):
        clustermap_spec(dm, ["x", "y", "z"], n_clusters=2)


def test_the_order_and_grouped_reports_are_both_reachable():
    dm = _dm([0.1, 0.9, 0.9, 0.9, 0.9, 0.1], 4)
    spec = clustermap_spec(dm, ["x", "x", "y", "y"], n_clusters=2)
    flat = spec.to_row_names_csv()
    assert flat.splitlines()[0] == "sample_id"
    assert flat.splitlines()[1:] == list(spec.ordered_sample_ids())
    grouped = spec.to_grouped_csv()
    assert grouped.startswith("Cluster 1:\n")
    assert grouped.count("Cluster ") == 2
