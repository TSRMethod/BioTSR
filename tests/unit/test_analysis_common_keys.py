"""The four common-key definitions, kept apart on purpose."""

from __future__ import annotations

import pytest

from biotsr.analysis import (
    AnalysisSample,
    common_key_intersection,
    common_key_intersection_by_group,
    key_occurrence_percent,
    key_occurrence_percent_one_vs_all,
)
from biotsr.analysis.samples import IncompatibleProcessing, processing_fingerprint
from biotsr.domain.keyspace import (
    IncompatibleKeySpace,
    IncompatibleSampling,
    SamplingProfile,
)
from biotsr.domain.selection_record import (
    SCOPE_ALL_COPIES,
    SCOPE_ONE_INSTANCE,
    SCOPE_RESIDUE_TYPES,
    SCOPE_WHOLE_CHAIN,
    selection_record,
    selection_sides,
    side_record,
)

from .test_analysis_features import key_space, sample, sampling


def test_global_intersection_is_membership_in_every_sample():
    s = [
        sample("a", {1: 5, 2: 1, 3: 1}),
        sample("b", {1: 2, 2: 9}),
        sample("c", {1: 1, 2: 4, 9: 1}),
    ]
    report = common_key_intersection(s)
    assert report.keys == (1, 2)
    assert report.n_common_keys == 2


def test_the_intersection_ignores_frequency_entirely():
    """A key present once counts exactly as much as a key present a thousand times."""
    rare = common_key_intersection([sample("a", {1: 1}), sample("b", {1: 1})]).keys
    common = common_key_intersection([sample("a", {1: 1000}), sample("b", {1: 999})]).keys
    assert rare == common == (1,)


def test_the_per_sample_report_carries_the_historical_columns():
    s = [sample("a", {1: 5, 2: 1, 3: 1}), sample("b", {1: 2, 2: 9})]
    rows = common_key_intersection(s).rows
    assert [(r.sample_id, r.total_keys, r.n_common_keys, r.sum_common_key_frequency)
            for r in rows] == [("a", 3, 2, 6), ("b", 2, 2, 11)]


def test_an_empty_intersection_is_a_result_not_an_error():
    """1UBQ's 75 residue samples share no key at all; the report still balances."""
    report = common_key_intersection([sample("a", {1: 1}), sample("b", {2: 1})])
    assert report.keys == ()
    assert report.is_empty
    assert [r.sum_common_key_frequency for r in report.rows] == [0, 0]
    assert [r.total_keys for r in report.rows] == [1, 1]


def test_a_single_sample_is_its_own_intersection():
    report = common_key_intersection([sample("a", {4: 2, 1: 1})])
    assert report.keys == (1, 4)


def test_intersection_keys_come_back_numerically_ordered():
    s = [sample("a", {10: 1, 9: 1, -3: 1}), sample("b", {10: 1, 9: 1, -3: 1})]
    assert common_key_intersection(s).keys == (-3, 9, 10)


def test_group_wise_intersection_is_computed_inside_each_group():
    s = [
        sample("x_DG", {1: 1, 2: 1}),
        sample("y_DG", {1: 1, 3: 1}),
        sample("z_DA", {7: 1}),
    ]
    by_group = common_key_intersection_by_group(s, lambda t: t.sample_id.split("_")[-1])
    assert sorted(by_group) == ["DA", "DG"]
    assert by_group["DG"].keys == (1,)
    assert by_group["DA"].keys == (7,)
    assert by_group["DG"].group == "DG"


def test_occurrence_percent_counts_samples_not_frequencies():
    s = [sample("a", {1: 100, 2: 1}), sample("b", {1: 1}), sample("c", {3: 1})]
    table = key_occurrence_percent(s)
    percents = table.percents()
    assert percents[1] == pytest.approx(200 / 3)
    assert percents[2] == pytest.approx(100 / 3)
    assert percents[3] == pytest.approx(100 / 3)
    assert table.total_samples == 3


def test_occurrence_percent_is_unrounded_and_rounding_belongs_to_export():
    """Three legacy programs round this to 0, 4 and 3 decimals; the number does not."""
    s = [sample("a", {1: 1}), sample("b", {}), sample("c", {})]
    value = key_occurrence_percent(s).percents()[1]
    assert value == pytest.approx(100 / 3)
    assert round(value, 4) == 33.3333
    assert round(value, 3) == 33.333


def test_one_vs_all_fills_absent_groups_with_zero():
    s = [
        sample("a_DG", {1: 1}),
        sample("b_DG", {1: 1, 2: 1}),
        sample("c_DA", {2: 1}),
    ]
    table = key_occurrence_percent_one_vs_all(s, lambda t: t.sample_id.split("_")[-1])
    assert table.groups == ("DA", "DG")
    assert table.keys == (1, 2)
    assert table.percents[1] == (0.0, 100.0)
    assert table.percents[2] == (100.0, 50.0)
    assert dict(table.totals) == {"DA": 1, "DG": 2}


def test_every_common_key_operation_refuses_a_mixed_key_space():
    a = sample("a", {1: 1})
    b = sample("b", {1: 1}, ks=key_space(label="grouped"))
    for call in (
        lambda: common_key_intersection([a, b]),
        lambda: common_key_intersection_by_group([a, b], lambda s: "g"),
        lambda: key_occurrence_percent([a, b]),
        lambda: key_occurrence_percent_one_vs_all([a, b], lambda s: "g"),
    ):
        with pytest.raises(IncompatibleKeySpace):
            call()


def test_every_common_key_operation_refuses_an_empty_set():
    for call in (
        lambda: common_key_intersection([]),
        lambda: common_key_intersection_by_group([], lambda s: "g"),
        lambda: key_occurrence_percent([]),
        lambda: key_occurrence_percent_one_vs_all([], lambda s: "g"),
    ):
        with pytest.raises(ValueError):
            call()


def test_the_four_definitions_answer_different_questions():
    """A key in 2 of 3 samples is 66.7% common and is *not* in the intersection."""
    s = [sample("a", {1: 1, 2: 1}), sample("b", {1: 1, 2: 1}), sample("c", {1: 1})]
    assert common_key_intersection(s).keys == (1,)
    assert key_occurrence_percent(s).percents()[2] == pytest.approx(200 / 3)


# ---------------------------------------------------------------------------
# The comparison rule: key definitions always; processing on explicit opt-in;
# selections never (ADR-094)
# ---------------------------------------------------------------------------

#: The four operations, each called with or without the opt-in.
OPERATIONS = {
    "intersection": lambda s, **kw: common_key_intersection(s, **kw),
    "intersection_by_group": lambda s, **kw: common_key_intersection_by_group(
        s, lambda t: "g", **kw),
    "occurrence_percent": lambda s, **kw: key_occurrence_percent(s, **kw),
    "occurrence_one_vs_all": lambda s, **kw: key_occurrence_percent_one_vs_all(
        s, lambda t: "g", **kw),
}


def cross_sampling(chain_a, chain_b, entity="ESW_404", cutoff=4.0):
    """A drug-protein cross sampling profile, as a cross run records it."""
    return SamplingProfile(
        selector="cross",
        scope="pooled",
        candidate_strategy="exhaustive_c_n_3",
        cross_method="drug_protein_cross",
        side_a_chain=chain_a,
        side_b_chain=chain_b,
        side_a_selector="protein cross_protein_atoms chain=%s" % chain_a,
        side_b_selector="drug cross_entity_atoms chain=%s entity=%s" % (chain_b, entity),
        contact_cutoff=cutoff,
        contact_filter_rule="symmetric_atom_contact_filter",
    )


def cross_sample(name, counts, chain_a="A", chain_b="A", entity="ESW", entity_id="404",
                 cutoff=4.0, atoms=(), residue_types=(), entity_names=(),
                 protein_selector="cross_protein_chain_atoms", altloc="blank_or_a"):
    """A drug-protein cross sample carrying the selection record a real run writes."""
    record = selection_record([
        side_record(side="a", role="protein", selector=protein_selector,
                    chain=chain_a,
                    scope=SCOPE_RESIDUE_TYPES if residue_types else SCOPE_WHOLE_CHAIN,
                    residue_types=residue_types),
        side_record(side="b", role="drug", selector="cross_entity_atoms",
                    chain=chain_b,
                    scope=SCOPE_ONE_INSTANCE if entity_id else SCOPE_ALL_COPIES,
                    entity_name=entity or None,
                    entity_id=entity_id if entity else None, atoms=atoms,
                    entity_names=entity_names),
    ])
    return AnalysisSample(
        sample_id=name, counts=counts, key_space=key_space(),
        sampling=cross_sampling(chain_a, chain_b, "%s_%s" % (entity, entity_id), cutoff),
        annotations={"method": "drug_protein_cross", "altloc_policy": altloc,
                     "selection": selection_sides(record)},
    )


def _mixed_pair():
    """Same key space, different processing (a size filter on one side)."""
    import dataclasses

    a = sample("a", {1: 3, 2: 1})
    b = sample("b", {1: 1, 3: 2}, sp=dataclasses.replace(sampling(), size_filter=25.0))
    return a, b


@pytest.mark.parametrize("name", sorted(OPERATIONS))
def test_mixed_processing_is_refused_by_default(name):
    a, b = _mixed_pair()
    with pytest.raises(IncompatibleProcessing) as exc:
        OPERATIONS[name]([a, b])
    message = str(exc.value)
    assert "size_filter" in message
    assert "allow_mixed_processing=True" in message
    # Callers catching the domain exception keep working.
    assert isinstance(exc.value, IncompatibleSampling)


@pytest.mark.parametrize("name", sorted(OPERATIONS))
def test_mixed_processing_is_allowed_with_the_explicit_opt_in(name):
    a, b = _mixed_pair()
    OPERATIONS[name]([a, b], allow_mixed_processing=True)  # does not raise


def test_an_allowed_mixed_intersection_warns_and_records_each_sample():
    a, b = _mixed_pair()
    report = common_key_intersection([a, b], allow_mixed_processing=True)
    assert report.keys == (1,)
    assert len(report.warnings) == 1
    assert "size_filter" in report.warnings[0]
    assert "explicitly allowed" in report.warnings[0]
    assert report.comparison["processing_uniform"] is False
    assert report.comparison["mixed_processing_allowed"] is True
    assert report.comparison["n_processing_policies"] == 2
    assert report.comparison["processing_differences"] == ["size_filter"]
    assert [p["sample_id"] for p in report.samples] == ["a", "b"]
    for provenance, s in zip(report.samples, (a, b)):
        assert provenance["key_space_fingerprint"] == s.key_space.fingerprint
        assert provenance["sampling_fingerprint"] == s.sampling.fingerprint
        assert provenance["processing_fingerprint"] == processing_fingerprint(s)
        assert provenance["processing"]["sampling"]["size_filter"] == s.sampling.size_filter
    assert report.samples[0]["processing_fingerprint"] != (
        report.samples[1]["processing_fingerprint"])
    payload = report.to_dict()
    assert payload["warnings"] == list(report.warnings)
    assert payload["comparison"]["processing_uniform"] is False
    assert len(payload["samples"]) == 2


def test_an_allowed_mixed_group_report_warns_inside_each_group():
    a, b = _mixed_pair()
    by_group = common_key_intersection_by_group(
        [a, b], lambda t: "g", allow_mixed_processing=True)
    assert by_group["g"].warnings and "size_filter" in by_group["g"].warnings[0]
    assert len(by_group["g"].samples) == 2


def test_uniform_processing_carries_no_warning_and_still_records_provenance():
    a, b = sample("a", {1: 1}), sample("b", {1: 2})
    report = common_key_intersection([a, b])
    assert report.warnings == ()
    assert report.comparison["processing_uniform"] is True
    assert report.comparison["mixed_processing_allowed"] is False
    assert report.comparison["n_processing_policies"] == 1
    assert [p["sample_id"] for p in report.samples] == ["a", "b"]


@pytest.mark.parametrize("name", sorted(OPERATIONS))
@pytest.mark.parametrize("change", [
    {"chain_a": "B", "chain_b": "C"},
    {"entity": "ATP", "entity_id": "501"},
    {"entity_id": None},
    {"atoms": ("C1", "C2", "N1")},
    {"residue_types": ("TYR",)},
    {"protein_selector": "cross_residue_type_atoms", "residue_types": ("LYS",)},
    {"entity": "", "entity_names": ("NAG", "MAN")},
])
def test_a_selection_difference_never_blocks_a_comparison(name, change):
    """Other chains, another ligand, all copies instead of one, an atom subset,
    a residue-type narrowing, another selector, a vocabulary: each is what was
    selected, not how triangles were drawn from it.  Compared without the
    opt-in and without a warning; recorded per sample."""
    a = cross_sample("a", {1: 1, 2: 1})
    b = cross_sample("b", {1: 2}, **change)
    assert processing_fingerprint(a) == processing_fingerprint(b)
    result = OPERATIONS[name]([a, b])
    assert getattr(result, "warnings", ()) == () or all(
        not r.warnings for r in result.values())


def test_a_selection_difference_is_recorded_side_by_side():
    a = cross_sample("a", {1: 1, 2: 1}, "A", "A")
    b = cross_sample("b", {1: 2}, "B", "C", entity="BEN", entity_id="1")
    report = common_key_intersection([a, b])
    assert report.warnings == ()
    assert report.comparison["processing_uniform"] is True
    assert report.comparison["n_sampling_fingerprints"] == 2
    assert report.samples[0]["chains"] == {"side_a_chain": "A", "side_b_chain": "A"}
    assert report.samples[1]["chains"] == {"side_a_chain": "B", "side_b_chain": "C"}
    entity = report.samples[1]["selection"][1]
    assert (entity["chain"], entity["scope"], entity["arguments"]) == (
        "C", SCOPE_ONE_INSTANCE, {"entity_name": "BEN", "entity_id": "1"})


@pytest.mark.parametrize("name", sorted(OPERATIONS))
@pytest.mark.parametrize("change, field", [
    ({"cutoff": 5.0}, "contact_cutoff"),
    ({"altloc": None}, "altloc_policy"),
])
def test_a_processing_difference_is_refused_and_named(name, change, field):
    """A different contact cutoff or a different (or unrecorded)
    alternate-location rule draws the frequencies differently: refused
    without the opt-in, named, and allowed with it."""
    a = cross_sample("a", {1: 1})
    b = cross_sample("b", {1: 1}, "B", "B", **change)
    with pytest.raises(IncompatibleProcessing, match=field):
        OPERATIONS[name]([a, b])
    OPERATIONS[name]([a, b], allow_mixed_processing=True)


@pytest.mark.parametrize("name", sorted(OPERATIONS))
def test_a_key_space_mismatch_is_refused_even_with_the_opt_in(name):
    a = sample("a", {1: 1})
    b = sample("b", {1: 1}, ks=key_space(label="grouped"))
    with pytest.raises(IncompatibleKeySpace):
        OPERATIONS[name]([a, b], allow_mixed_processing=True)


def test_the_key_space_is_checked_before_processing():
    """Both differ: the key-space refusal wins, whatever the opt-in says."""
    import dataclasses

    a = sample("a", {1: 1})
    b = sample("b", {1: 1}, ks=key_space(label="grouped"),
               sp=dataclasses.replace(sampling(), size_filter=25.0))
    for allow in (False, True):
        with pytest.raises(IncompatibleKeySpace):
            common_key_intersection([a, b], allow_mixed_processing=allow)


@pytest.mark.parametrize("name", ["occurrence_percent", "occurrence_one_vs_all"])
def test_an_allowed_mixed_occurrence_table_records_warnings_and_provenance(name):
    """Every common-key operation's result records per-sample identities
    (``.samples``), a ``.comparison`` summary and ``.warnings``."""
    a, b = _mixed_pair()
    table = OPERATIONS[name]([a, b], allow_mixed_processing=True)
    assert [p["sample_id"] for p in table.samples] == ["a", "b"]
    assert table.comparison["processing_uniform"] is False
    assert table.warnings and "size_filter" in table.warnings[0]
