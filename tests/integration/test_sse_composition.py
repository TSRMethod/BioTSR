"""SSE against every supported Ca configuration.

Two claims are under test, and they are the whole of the SSE layer's science:

1. **Key invariance.** Annotating a run changes no key, no sign, no count and
   neither identity fingerprint. Asserted for regular Ca, selected Ca, grouping,
   mirror, size filtering and every valid combination of them -- not once on a
   representative case.
2. **Frequency conservation.** The SSE categories partition the retained
   triangles exactly: the category totals sum to ``n_triangles``, and each key's
   category vector sums to that key's own frequency. There is no missing-SSE
   escape hatch that could make this inexact, because coil is a state rather
   than an absence (ADR-065).

The β-strand identity convention (ADR-063) is a fact about the annotation: it is
recorded in the annotation's provenance and in neither identity fingerprint.
"""

from __future__ import annotations

import itertools
import os

import pytest

import biotsr
from biotsr.sse import SSE_CATEGORIES

CHAIN = "A"

#: The three modifiers, as the eight on/off combinations the package promises.
MODIFIER_COMBINATIONS = [
    {"aa_grouping": g, "mirror": m, **({"size_filter": 25.0} if s else {})}
    for g, m, s in itertools.product((False, True), repeat=3)
]


def _label(combo):
    on = [k for k, v in combo.items() if v]
    return "+".join(on) if on else "plain"


@pytest.fixture(scope="module")
def ubq(ubq_pdb):
    return ubq_pdb


# ---------------------------------------------------------------------------
# 1. Key invariance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "combo", MODIFIER_COMBINATIONS, ids=[_label(c) for c in MODIFIER_COMBINATIONS]
)
def test_sse_annotation_changes_no_key_for_any_modifier_combination(ubq, combo):
    plain = biotsr.analyze(ubq, method="protein", chain=CHAIN, **combo)
    annotated = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True, **combo)

    assert annotated.frequencies == plain.frequencies
    assert annotated.key_order == plain.key_order
    assert annotated.counts.to_dict() == plain.counts.to_dict()
    assert annotated.key_space.fingerprint == plain.key_space.fingerprint
    assert annotated.sampling.fingerprint == plain.sampling.fingerprint


@pytest.mark.parametrize(
    "combo", MODIFIER_COMBINATIONS, ids=[_label(c) for c in MODIFIER_COMBINATIONS]
)
def test_every_triplet_is_geometrically_identical_with_and_without_sse(ubq, combo):
    """Not just the aggregate: the per-triangle record matches field for field."""
    plain = biotsr.analyze(ubq, method="protein", chain=CHAIN, **combo)
    annotated = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True, **combo)
    assert len(plain.triplets) == len(annotated.triplets)
    for a, b in zip(plain.triplets, annotated.triplets):
        assert a.key == b.key
        assert a.labels == b.labels
        assert a.identities == b.identities
        assert a.theta_raw == b.theta_raw
        assert a.theta == b.theta
        assert a.theta_class == b.theta_class
        assert a.distance_class == b.distance_class
        assert a.max_distance == b.max_distance
        assert a.coordinates == b.coordinates
        assert a.sse is None and b.sse is not None


def test_key_invariance_holds_for_selected_calpha(ubq):
    """``aa_selection`` stays a Ca method: radix 30, never atomic radix 112."""
    selection = "1-7;23-34;40-45;56-59"
    plain = biotsr.analyze(
        ubq, method="aa_selection", chain=CHAIN, residues=selection
    )
    annotated = biotsr.analyze(
        ubq, method="aa_selection", chain=CHAIN, residues=selection, sse=True
    )
    assert annotated.frequencies == plain.frequencies
    assert annotated.key_space.fingerprint == plain.key_space.fingerprint
    assert annotated.key_space.label_scheme.radix == 30
    assert annotated.key_space.distance_scheme.radix == 35


@pytest.mark.parametrize(
    "combo", MODIFIER_COMBINATIONS, ids=[_label(c) for c in MODIFIER_COMBINATIONS]
)
def test_key_invariance_holds_for_selected_calpha_under_every_modifier(ubq, combo):
    selection = "1-7;23-34;56-59;64-72"
    kwargs = dict(method="aa_selection", chain=CHAIN, residues=selection, **combo)
    plain = biotsr.analyze(ubq, **kwargs)
    annotated = biotsr.analyze(ubq, sse=True, **kwargs)
    assert annotated.frequencies == plain.frequencies
    assert annotated.sampling.fingerprint == plain.sampling.fingerprint


# ---------------------------------------------------------------------------
# 2. Frequency conservation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "combo", MODIFIER_COMBINATIONS, ids=[_label(c) for c in MODIFIER_COMBINATIONS]
)
def test_sse_categories_partition_the_retained_triangles(ubq, combo):
    result = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True, **combo)
    totals = result.category_totals()
    assert set(totals) <= set(SSE_CATEGORIES)
    assert sum(totals.values()) == result.n_triangles
    for key, counts in result.keys.items():
        assert counts.total == sum(counts.by_category.values())
    assert result.counts.ledger_balances


def test_each_triplets_category_is_its_keys_bucket(ubq):
    import collections

    result = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True)
    rebuilt = collections.defaultdict(collections.Counter)
    for t in result.triplets:
        assert t.category == t.sse.category
        rebuilt[t.key][t.sse.category] += 1
    assert {k: dict(v) for k, v in rebuilt.items()} == {
        k: dict(v.by_category) for k, v in result.keys.items()
    }


# ---------------------------------------------------------------------------
# Mirror: the sign is theta_raw's, never the category's
# ---------------------------------------------------------------------------


def test_the_mirror_sign_is_decided_by_theta_raw_and_not_by_sse(ubq):
    result = biotsr.analyze(ubq, method="protein", chain=CHAIN, mirror=True, sse=True)
    assert any(k < 0 for k in result.keys), "1UBQ must produce mirrored keys"
    for t in result.triplets:
        assert (t.key < 0) == (t.theta_raw > 90.0)


def test_negative_and_positive_keys_are_both_sse_stratified(ubq):
    result = biotsr.analyze(ubq, method="protein", chain=CHAIN, mirror=True, sse=True)
    negative = {k: v for k, v in result.keys.items() if k < 0}
    positive = {k: v for k, v in result.keys.items() if k > 0}
    assert negative and positive
    for bucket in (negative, positive):
        categories = set()
        for counts in bucket.values():
            categories |= set(counts.by_category)
        assert categories <= set(SSE_CATEGORIES)
        assert len(categories) > 1


def test_mirror_and_plain_split_the_same_triangles_the_same_way(ubq):
    """Turning SSE on must not move a triangle between the +k and -k buckets."""
    plain = biotsr.analyze(ubq, method="protein", chain=CHAIN, mirror=True)
    annotated = biotsr.analyze(ubq, method="protein", chain=CHAIN, mirror=True, sse=True)
    assert plain.frequencies == annotated.frequencies


# ---------------------------------------------------------------------------
# Grouping: label space collapses, SSE identity does not
# ---------------------------------------------------------------------------


def test_grouping_changes_the_label_space_but_not_the_sse_annotation(ubq):
    ungrouped = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True)
    grouped = biotsr.analyze(
        ubq, method="protein", chain=CHAIN, aa_grouping=True, sse=True
    )
    assert ungrouped.key_space.label_scheme.name == "amino_acid_30"
    assert grouped.key_space.label_scheme.name == "amino_acid_group_23"
    assert grouped.key_space.fingerprint != ungrouped.key_space.fingerprint
    # Same triangles, same count, same SSE distribution: grouping relabels
    # vertices, it does not re-annotate them.
    assert grouped.n_triangles == ungrouped.n_triangles
    assert grouped.category_totals() == ungrouped.category_totals()


def test_grouping_keeps_the_original_residues_secondary_structure(ubq):
    """The grouped label may print as another residue; the SSE state is the real one."""
    grouped = biotsr.analyze(
        ubq, method="protein", chain=CHAIN, aa_grouping=True, sse=True
    )
    st = biotsr.load_structure(ubq)
    from biotsr.sse import build_assignment_map

    truth = build_assignment_map(st, CHAIN)
    renamed = 0
    for t in grouped.triplets:
        for identity, state, group in zip(t.identities, t.sse.states, t.sse.groups):
            expected = truth.assignment_for(identity)
            # `groups` carries the element identity's display spelling; the
            # structured identity itself stays inside `biotsr.sse` (ADR-063).
            assert (state, group) == (expected.state, expected.group.text)
        renamed += sum(
            1
            for identity, name in zip(t.identities, t.label_names)
            if identity.residue_name != name
        )
    assert renamed > 0, "grouping must actually rename some vertices on this fixture"


# ---------------------------------------------------------------------------
# Size filter: it decides membership before annotation ever runs
# ---------------------------------------------------------------------------


def test_the_size_filter_decides_membership_and_sse_does_not_alter_it(ubq):
    unfiltered = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True)
    filtered = biotsr.analyze(
        ubq, method="protein", chain=CHAIN, size_filter=20.0, sse=True
    )
    assert filtered.counts.n_filtered_by_size > 0
    assert filtered.n_triangles < unfiltered.n_triangles
    plain = biotsr.analyze(ubq, method="protein", chain=CHAIN, size_filter=20.0)
    assert plain.counts.to_dict() == filtered.counts.to_dict()
    assert sum(filtered.category_totals().values()) == filtered.n_triangles


def test_the_size_filter_is_off_by_default(ubq):
    result = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True)
    assert result.sampling.size_filter is None
    assert result.counts.n_filtered_by_size == 0


def test_every_retained_triangle_is_strictly_below_the_cutoff(ubq):
    """The historical Ca rule is ``max_distance < cutoff``, not ``<=``."""
    result = biotsr.analyze(
        ubq, method="protein", chain=CHAIN, size_filter=20.0, sse=True
    )
    assert all(t.max_distance < 20.0 for t in result.triplets)


# ---------------------------------------------------------------------------
# Structures chosen for the secondary structure they contain
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fixture,chain",
    [("barnase_path", "A"), ("p53_path", "A"), ("nucleosome_path", "A")],
)
def test_sse_over_helix_rich_and_sheet_rich_structures(request, fixture, chain):
    path = request.getfixturevalue(fixture)
    plain = biotsr.analyze(path, method="protein", chain=chain)
    annotated = biotsr.analyze(path, method="protein", chain=chain, sse=True)
    assert annotated.frequencies == plain.frequencies
    totals = annotated.category_totals()
    assert sum(totals.values()) == annotated.n_triangles
    assert len(totals) >= 10, "a real structure should reach most of the 18"


def test_a_structure_without_secondary_structure_records_is_all_coil_and_says_so(
    bdna_path,
):
    """1BNA is B-DNA. Every triangle is ``17_3c`` -- and the result warns why.

    This is the one case where the three-state collapse is lossy, and it is
    closed by reporting rather than by inventing a fourth state (ADR-065).
    """
    result = biotsr.analyze(bdna_path, method="protein", chain="A", sse=True)
    if result.n_triangles:
        assert set(result.category_totals()) == {"17_3c"}
    assert result.annotations["sse"]["has_records"] is False
    assert any("no HELIX or SHEET record" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------


def test_sse_is_refused_for_atomic_and_entity_families(ubq, trypsin_path):
    with pytest.raises(ValueError, match="does not support"):
        biotsr.analyze(
            ubq,
            method="amino_acid_atomic",
            chain="A",
            residue="15",
            residue_name="LEU",
            sse=True,
        )
    with pytest.raises(ValueError, match="does not support"):
        biotsr.analyze(
            trypsin_path,
            method="drug",
            chain="A",
            entity_name="BEN",
            entity_id="1",
            sse=True,
        )


def test_sse_is_refused_through_the_per_entity_entry_point_too(ubq):
    """``analyze_entities`` routes through the same validation, not around it."""
    with pytest.raises(ValueError, match="does not support 'sse'"):
        biotsr.analyze_entities(
            ubq, method="amino_acid_atomic", chain="A", sse=True
        )


# ---------------------------------------------------------------------------
# 3. The strand-identity convention is provenance, not identity
# ---------------------------------------------------------------------------


def test_the_convention_is_recorded_on_the_result_and_in_neither_identity(ubq):
    """Provenance, not identity -- the placement ADR-063 argues for.

    A `SamplingProfile` field would break key invariance; a `KeySpace`
    field would make two comparable results refuse to be compared. The
    annotation's own provenance is where a fact about the annotation belongs.
    """
    import dataclasses

    from biotsr.domain.keyspace import KeySpace, SamplingProfile
    from biotsr.sse import STRAND_IDENTITY_SHEET_ID_AND_SERIAL

    annotated = biotsr.analyze(ubq, method="protein", chain=CHAIN, sse=True)
    assert annotated.annotations["sse"]["strand_identity"] == "sheet_id_and_serial"
    assert (
        annotated.annotations["sse"]["strand_identity"]
        == STRAND_IDENTITY_SHEET_ID_AND_SERIAL
    )
    assert annotated.annotations["sse"]["version"] == "sse_assignment_v3"

    fields = {f.name for f in dataclasses.fields(KeySpace)} | {
        f.name for f in dataclasses.fields(SamplingProfile)
    }
    assert not any("strand" in f or "sse" in f for f in fields)


def test_an_unannotated_run_records_no_convention(ubq):
    """No annotation, nothing to have a convention about."""
    plain = biotsr.analyze(ubq, method="protein", chain=CHAIN)
    assert plain.annotations["sse"] is None
