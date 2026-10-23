"""Feature matrices: ordering, zero filling, negative keys, and the two identities."""

from __future__ import annotations

import pytest

from biotsr.analysis import (
    FIRST_ENCOUNTER,
    NUMERIC_ASCENDING,
    AnalysisSample,
    feature_matrix,
    label_from_name,
    require_common_key_space,
    require_common_sampling,
    sampling_variation,
)
from biotsr.domain.keyspace import (
    IncompatibleKeySpace,
    IncompatibleSampling,
    KeySpace,
    SamplingProfile,
    SchemeRef,
)


def _scheme(name, radix=20, fingerprint=None):
    return SchemeRef(name=name, content_fingerprint=fingerprint or ("fp-" + name), radix=radix)


def key_space(label="labels", theta="theta29"):
    return KeySpace(
        angle_convention="standard",
        sqrt_convention="math_sqrt",
        median_convention="midpoint",
        numerical_precision="float64",
        theta_scheme=_scheme(theta, 29),
        distance_scheme=_scheme("dist18", 18),
        label_scheme=_scheme(label),
        canonical_ordering_version="v1",
        key_formula_version="v1",
        mirror_convention="disabled",
    )


def sampling(cutoff=None, selector="calpha"):
    return SamplingProfile(
        selector=selector,
        scope="chain",
        candidate_strategy="exhaustive",
        contact_cutoff=cutoff,
    )


def sample(name, counts, ks=None, sp=None, label=None):
    return AnalysisSample(
        sample_id=name,
        counts=counts,
        key_space=ks or key_space(),
        sampling=sp or sampling(),
        label=label,
    )


# --- ordering ---------------------------------------------------------------


def test_numeric_ascending_sorts_keys_as_integers_not_strings():
    a = sample("a", {10: 1, 9: 2, -5: 3})
    b = sample("b", {100: 4, 2: 5})
    fm = feature_matrix([a, b], ordering=NUMERIC_ASCENDING)
    assert fm.columns == (-5, 2, 9, 10, 100)


def test_the_string_sorted_ordering_is_no_longer_offered():
    """Only the numeric and first-encounter orderings exist."""
    import biotsr.analysis as analysis

    assert analysis.FEATURE_ORDERINGS == (NUMERIC_ASCENDING, FIRST_ENCOUNTER)
    assert not hasattr(analysis, "LEGACY_STRING_SORTED")
    with pytest.raises(ValueError, match="unknown feature ordering"):
        feature_matrix([sample("a", {1: 1})], ordering="legacy_string_sorted")


def test_first_encounter_follows_the_order_keys_are_met():
    a = sample("a", dict([(9, 1), (10, 1)]))
    b = sample("b", dict([(3, 1), (9, 2)]))
    fm = feature_matrix([a, b], ordering=FIRST_ENCOUNTER)
    assert fm.columns == (9, 10, 3)


def test_the_hash_order_variant_is_not_offered():
    with pytest.raises(ValueError, match="not reproducible"):
        feature_matrix([sample("a", {1: 1})], ordering="legacy_hash_order")


def test_ordering_is_a_permutation_and_nothing_else():
    """Both orderings hold the same multiset of values per row."""
    samples = [sample("a", {-3: 2, 7: 1, 40: 5}), sample("b", {7: 9, 12: 1})]
    rows = {}
    for ordering in (NUMERIC_ASCENDING, FIRST_ENCOUNTER):
        fm = feature_matrix(samples, ordering=ordering)
        rows[ordering] = [sorted(fm.row(i)) for i in range(fm.n_samples)]
    assert rows[NUMERIC_ASCENDING] == rows[FIRST_ENCOUNTER]


# --- zeros, negatives, shape ------------------------------------------------


def test_a_missing_key_is_an_explicit_zero():
    fm = feature_matrix([sample("a", {1: 4}), sample("b", {2: 7})])
    assert fm.columns == (1, 2)
    assert fm.row(0) == (4, 0)
    assert fm.row(1) == (0, 7)
    assert fm.value(0, 2) == 0


def test_negative_mirror_keys_are_ordinary_columns():
    fm = feature_matrix([sample("a", {-713374251: 3, 713491179: 1})])
    assert fm.columns == (-713374251, 713491179)
    assert fm.row(0) == (3, 1)


def test_sparse_storage_and_dense_view_agree():
    fm = feature_matrix([sample("a", {1: 1, 5: 2}), sample("b", {5: 3})])
    assert fm.nnz == 3
    assert fm.shape == (2, 2)
    assert fm.dense() == ((1, 2), (0, 3))
    assert fm.density == pytest.approx(3 / 4)


def test_column_occurrence_counts_samples_not_frequencies():
    fm = feature_matrix([sample("a", {1: 50}), sample("b", {1: 1, 2: 1})])
    assert dict(fm.column_occurrence()) == {1: 2, 2: 1}


def test_row_order_is_the_caller_s_order():
    fm = feature_matrix([sample("z", {1: 1}), sample("a", {1: 1})])
    assert fm.sample_ids == ("z", "a")


def test_duplicate_sample_ids_are_refused():
    with pytest.raises(ValueError, match="duplicate sample id"):
        feature_matrix([sample("a", {1: 1}), sample("a", {2: 1})])


def test_an_empty_sample_set_is_refused():
    with pytest.raises(ValueError):
        feature_matrix([])


# --- typing guards ----------------------------------------------------------


def test_string_keys_are_refused_at_construction():
    """Legacy read keys as text and intersected strings; this layer never does."""
    with pytest.raises(TypeError, match="never compared as strings"):
        AnalysisSample("a", {"713491179": 1}, key_space(), sampling())


def test_negative_frequencies_are_refused():
    with pytest.raises(ValueError, match="negative"):
        AnalysisSample("a", {1: -1}, key_space(), sampling())


# --- the two identities -----------------------------------------------------


def test_a_differing_key_space_is_refused_and_the_difference_is_named():
    a = sample("a", {1: 1})
    b = sample("b", {1: 1}, ks=key_space(label="grouped"))
    with pytest.raises(IncompatibleKeySpace) as exc:
        feature_matrix([a, b])
    message = str(exc.value)
    assert "'a'" in message and "'b'" in message
    assert "label_scheme" in message


def test_a_processing_difference_is_refused_unless_explicitly_allowed():
    """A contact cutoff change keeps every key's meaning but draws the
    frequencies differently, so a matrix over both is refused by default --
    the same rule common keys apply -- and recorded when allowed."""
    a = sample("a", {1: 1}, sp=sampling(cutoff=3.0))
    b = sample("b", {1: 2}, sp=sampling(cutoff=3.5))
    with pytest.raises(IncompatibleSampling, match="contact_cutoff") as exc:
        feature_matrix([a, b])
    assert "allow_mixed_processing=True" in str(exc.value)
    fm = feature_matrix([a, b], allow_mixed_processing=True)
    assert fm.n_samples == 2
    assert fm.comparison["processing_uniform"] is False
    assert fm.comparison["mixed_processing_allowed"] is True
    assert fm.comparison["n_processing_policies"] == 2
    assert fm.warnings and "contact_cutoff" in fm.warnings[0]
    assert len(set(fm.sampling_fingerprints)) == 2
    assert [p["sampling"]["contact_cutoff"] for p in fm.samples] == [3.0, 3.5]
    assert fm.provenance()["warnings"] == list(fm.warnings)


def test_the_opt_in_never_relaxes_the_key_space():
    a = sample("a", {1: 1}, sp=sampling(cutoff=3.0))
    b = sample("b", {1: 2}, ks=key_space(label="grouped"), sp=sampling(cutoff=3.5))
    with pytest.raises(IncompatibleKeySpace):
        feature_matrix([a, b], allow_mixed_processing=True)


def test_require_common_key_space_returns_the_shared_space():
    a, b = sample("a", {1: 1}), sample("b", {2: 1})
    assert require_common_key_space([a, b]).fingerprint == a.key_space.fingerprint


def test_require_common_sampling_returns_the_shared_profile():
    a, b = sample("a", {1: 1}), sample("b", {2: 1})
    assert require_common_sampling([a, b]).fingerprint == a.sampling.fingerprint


def test_sampling_variation_describes_without_judging():
    a = sample("a", {1: 1}, sp=sampling(cutoff=3.0, selector="calpha"))
    b = sample("b", {1: 1}, sp=sampling(cutoff=3.0, selector="atoms"))
    summary = sampling_variation([a, b])
    assert summary["varying_fields"] == ("selector",)
    assert summary["uniform"] is False


# --- labels -----------------------------------------------------------------


def test_label_from_name_takes_the_token_the_dataset_convention_says():
    assert label_from_name("1BNA_A_10_DG", 3) == "DG"
    assert label_from_name("1BNA_A_10_DG", 2) == "10"
    assert label_from_name("6ORP_B_614_NAG_NCI", 3) == "NAG"
    assert label_from_name("1UBQ_A_57_SER", -1) == "SER"


def test_label_from_name_refuses_a_token_that_does_not_exist():
    with pytest.raises(ValueError, match="token 7 does not exist"):
        label_from_name("1BNA_A_10_DG", 7)


# --- export -----------------------------------------------------------------


def test_the_headerless_legacy_csv_export_is_gone():
    fm = feature_matrix([sample("x", {7: 1})])
    assert not hasattr(fm, "to_legacy_csv")


def test_header_csv_carries_the_key_names_and_uses_crlf():
    fm = feature_matrix([sample("x", {7: 1})])
    assert fm.to_header_csv() == "sample_id,7\r\nx,1\r\n"


def test_an_explicitly_stored_zero_is_dropped_at_construction():
    """A zero count and an absent key are the same fact, so only one is stored."""
    s = AnalysisSample("a", {1: 3, 2: 0, 3: 1}, key_space(), sampling())
    assert dict(s.counts) == {1: 3, 3: 1}
    assert s.n_keys == 2
    assert s.total_frequency == 4


def test_a_column_that_is_zero_everywhere_is_not_a_column():
    fm = feature_matrix([sample("a", {1: 1, 9: 0}), sample("b", {2: 1, 9: 0})])
    assert fm.columns == (1, 2)


def test_presence_based_measures_agree_with_the_frequency_ones_about_presence():
    """The reason zeros are canonicalised: ``key in counts`` must mean ``value > 0``.

    Legacy spells presence ``a_jac[a_jac > 0] = 1``. If a stored zero counted as
    present, ``binary_jaccard`` and ``sarika`` would see a key that
    ``generalized_jaccard`` treats as absent.
    """
    from biotsr.analysis import binary_jaccard, generalized_jaccard, sarika

    with_zero = AnalysisSample("a", {1: 3, 2: 0}, key_space(), sampling()).counts
    without = AnalysisSample("a", {1: 3}, key_space(), sampling()).counts
    other = AnalysisSample("b", {2: 5}, key_space(), sampling()).counts
    for fn in (binary_jaccard, generalized_jaccard, sarika):
        assert fn(with_zero, other) == fn(without, other)
