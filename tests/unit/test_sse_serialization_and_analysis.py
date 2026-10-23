"""SSE through the serializers, and SSE meeting the analysis layer."""

from __future__ import annotations

import json
import os

import pytest

import biotsr
from biotsr.analysis.features import feature_matrix
from biotsr.analysis.samples import AnalysisSample
from biotsr.analysis.similarity import generalized_jaccard
from biotsr.io import native
from biotsr.sse import (
    SSE_CATEGORIES,
    IncompatibleSSEStratification,
    category_totals_in_order,
    check_conservation,
    require_common_strand_identity,
    sse_annotation_summary,
    strand_identity_of,
    stratified_key_counts,
)


@pytest.fixture(scope="module")
def annotated(ubq_pdb):
    return biotsr.analyze(ubq_pdb, method="protein", chain="A", sse=True)


@pytest.fixture(scope="module")
def plain(ubq_pdb):
    return biotsr.analyze(ubq_pdb, method="protein", chain="A")


# ---------------------------------------------------------------------------
# Native serialization
# ---------------------------------------------------------------------------


def test_the_native_key_file_is_already_sse_stratified(annotated, tmp_path):
    """No new format: the key file has been ``key/category/frequency`` since 2B."""
    path = native.write_keys(annotated, str(tmp_path / "keys.tsv"))
    rows = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    assert {r[1] for r in rows} <= set(SSE_CATEGORIES)
    assert sum(int(r[2]) for r in rows) == annotated.n_triangles


def test_the_native_triplet_width_does_not_depend_on_whether_sse_was_requested(
    annotated, plain
):
    width = len(native.TRIPLET_COLUMNS)
    for result in (annotated, plain):
        line = native.format_native_triplet(result.triplets[0])
        assert len(line.rstrip("\n").split("\t")) == width
    assert len(native.triplet_header().rstrip("\n").split("\t")) == width


def test_an_unannotated_run_writes_the_absent_marker_not_coil(plain):
    """``-`` and ``NONE`` are different facts and must not share a spelling."""
    fields = native.format_native_triplet(plain.triplets[0]).rstrip("\n").split("\t")
    assert fields[-7:] == [native.SSE_ABSENT] * 7
    assert native.SSE_ABSENT not in ("NONE", "COIL")


def test_the_annotated_triplet_carries_per_vertex_states_and_the_category(annotated):
    for t in annotated.triplets[:200]:
        fields = native.format_native_triplet(t).rstrip("\n").split("\t")
        states = fields[-7], fields[-5], fields[-3]
        groups = fields[-6], fields[-4], fields[-2]
        assert states == t.sse.states
        assert groups == t.sse.groups
        assert fields[-1] == t.sse.category == t.category


def test_the_metadata_sidecar_records_the_assignment_and_the_coverage(
    annotated, tmp_path
):
    path = native.write_metadata(annotated, str(tmp_path / "m.json"))
    meta = json.load(open(path))
    assert meta["sse"]["has_records"] is True
    # The SSE strand-identity convention travels in ANNOTATION provenance --
    # not in KeySpace, not in SamplingProfile (ADR-063).
    assert meta["sse"]["strand_identity"] == "sheet_id_and_serial"
    assert meta["sse"]["helix_identity"] == "helix_serial"
    assert sum(meta["sse_coverage"].values()) == annotated.selection.n_points
    assert sum(meta["category_totals"].values()) == annotated.n_triangles


def test_an_unannotated_results_metadata_says_so(plain, tmp_path):
    meta = json.load(open(native.write_metadata(plain, str(tmp_path / "m.json"))))
    assert meta["sse"] is None
    assert meta["sse_coverage"] is None


def test_save_writes_a_readable_triple(annotated, tmp_path):
    written = native.save(annotated, str(tmp_path))
    assert set(written) == {"keys", "metadata", "triplets"}
    for path in written.values():
        assert os.path.getsize(path) > 0


# ---------------------------------------------------------------------------
# Stratification: reshaping, not analysis
# ---------------------------------------------------------------------------


def test_stratified_counts_reconcile_with_the_flat_frequencies(annotated):
    strata = stratified_key_counts(annotated)
    rebuilt = {}
    for per_key in strata.values():
        for key, n in per_key.items():
            rebuilt[key] = rebuilt.get(key, 0) + n
    assert rebuilt == dict(annotated.frequencies)


def test_stratified_counts_hold_no_zeros(annotated):
    for per_key in stratified_key_counts(annotated).values():
        assert all(n > 0 for n in per_key.values())


def test_the_category_vector_is_dense_and_in_canonical_order(annotated):
    ordered = category_totals_in_order(annotated)
    assert [c for c, _ in ordered] == list(SSE_CATEGORIES)
    assert sum(n for _, n in ordered) == annotated.n_triangles


def test_conservation_is_checkable_and_holds(annotated):
    check_conservation(annotated)


def test_the_summary_distinguishes_annotated_from_unannotated(annotated, plain):
    assert sse_annotation_summary(annotated)["annotated"] is True
    assert sse_annotation_summary(plain)["annotated"] is False


# ---------------------------------------------------------------------------
# Interaction with the analysis layer
# ---------------------------------------------------------------------------


def test_an_sse_result_reduces_to_the_same_analysis_sample_as_a_plain_one(
    annotated, plain
):
    """Annotation must be invisible to every evidenced analysis operation."""
    a = AnalysisSample.from_result(annotated, sample_id="ubq")
    b = AnalysisSample.from_result(plain, sample_id="ubq")
    assert a.counts == b.counts
    assert a.key_space.fingerprint == b.key_space.fingerprint
    assert a.sampling.fingerprint == b.sampling.fingerprint
    assert a.total_frequency == b.total_frequency == annotated.n_triangles


def test_sse_and_plain_results_share_a_feature_matrix(annotated, plain):
    matrix = feature_matrix(
        [
            AnalysisSample.from_result(annotated, sample_id="annotated"),
            AnalysisSample.from_result(plain, sample_id="plain"),
        ]
    )
    assert matrix.rows[0] == matrix.rows[1]


def test_an_sse_result_is_identical_to_its_plain_twin_under_every_measure(
    annotated, plain
):
    a = AnalysisSample.from_result(annotated, sample_id="a")
    b = AnalysisSample.from_result(plain, sample_id="b")
    assert generalized_jaccard(a.counts, b.counts) == pytest.approx(1.0)


def test_the_analysis_layer_is_not_given_an_sse_feature_representation():
    """None is defined, deliberately -- no source establishes one (ADR-067).

    If a later release adds one it must arrive with its own evidence, and this test
    is what will fail to remind it.
    """
    import biotsr.analysis as analysis

    surface = " ".join(dir(analysis))
    assert "sse" not in surface.lower()


def test_the_stratification_module_offers_no_measures():
    from biotsr.sse import stratification

    for banned in ("jaccard", "distance", "similarity", "cluster", "linkage"):
        assert banned not in " ".join(dir(stratification)).lower()


# ---------------------------------------------------------------------------
# Streaming: the memory-flat path must carry the annotation too
# ---------------------------------------------------------------------------


def test_the_triplet_sink_receives_annotated_records(ubq_pdb):
    """``keep_triplets=False`` plus a sink is how a 1.2M-triangle run stays flat.

    If SSE only rode on retained records, the exact runs that need streaming
    would be the ones that silently lost the annotation.
    """
    import collections

    received = []
    result = biotsr.analyze(
        ubq_pdb,
        method="protein",
        chain="A",
        sse=True,
        keep_triplets=False,
        triplet_sink=received.append,
    )
    assert result.triplets is None
    assert len(received) == result.n_triangles
    assert all(r.sse is not None for r in received)
    assert (
        dict(collections.Counter(r.sse.category for r in received))
        == result.category_totals()
    )


def test_a_point_set_too_small_for_a_triangle_annotates_cleanly(ubq_pdb):
    """No triangles, no categories, and conservation still holds vacuously."""
    result = biotsr.analyze(
        ubq_pdb, method="aa_selection", chain="A", residues="1-2", sse=True
    )
    assert result.n_triangles == 0
    assert result.category_totals() == {}
    check_conservation(result)


# ---------------------------------------------------------------------------
# The strand-identity convention: preserved in provenance, and un-mixable
# ---------------------------------------------------------------------------
#
# The convention is in NEITHER identity (ADR-063), which is what makes ordinary
# key analysis across conventions correct -- and what makes `KeySpace` and
# `SamplingProfile` structurally unable to catch a bad *category* comparison.
# The package computes only the standard ``sheet_id_and_serial`` convention; a
# result carrying another one (e.g. read from an older output) is simulated here
# by relabelling a copy of a standard result's annotation, which is exactly the
# situation the guard exists for.
#
# 1TUP chain A has two β-sheets, `A` (4 strands) and `B` (7), numbering their
# strands from 1, so serials 1-4 collide across them. Five residues are enough --
# ARG 110 and LEU 111 in sheet `A` strand 1, CYS 124 in sheet `B` strand 1,
# LYS 132 in sheet `B` strand 2, SER 94 coil -- and it runs in one triangle
# enumeration over 5 points.

COLLIDING_RESIDUES = "110-111;124;132;94"
OTHER_CONVENTION = "legacy_serial_only"


@pytest.fixture(scope="module")
def tup_standard(p53_path):
    return biotsr.analyze(
        p53_path,
        method="aa_selection",
        chain="A",
        residues=COLLIDING_RESIDUES,
        sse=True,
    )


@pytest.fixture(scope="module")
def tup_other(tup_standard):
    """A copy of the standard result whose annotation names another convention."""
    import dataclasses

    sse = dict(tup_standard.annotations["sse"])
    sse["strand_identity"] = OTHER_CONVENTION
    annotations = dict(tup_standard.annotations)
    annotations["sse"] = sse
    return dataclasses.replace(tup_standard, annotations=annotations)


def test_the_standard_separates_same_serial_strands_of_different_sheets(tup_standard):
    """Sheet A strand 1 and sheet B strand 1 are two strands under the standard."""
    by_residues = {
        frozenset(i.residue_seq for i in t.identities): t.category
        for t in tup_standard.triplets
    }
    # two residues in A/1, one in B/1: two distinct strands
    assert by_residues[frozenset({110, 111, 124})] == "8_3b2"
    # A/1, B/1, B/2: three distinct strands
    assert by_residues[frozenset({110, 124, 132})] == "9_3b3"


def test_the_convention_survives_into_the_metadata_sidecar_on_disk(
    tup_standard, tmp_path
):
    """Serialized provenance: the categories and the convention that named them
    travel in the same file, so a sidecar can never be read without it."""
    written = native.save(tup_standard, str(tmp_path))
    meta = json.load(open(written["metadata"]))
    assert meta["sse"]["strand_identity"] == "sheet_id_and_serial"
    assert meta["sse"]["helix_identity"] == "helix_serial"
    assert meta["sse"]["version"] == "sse_assignment_v3"
    # ...in the same payload as the thing it qualifies.
    assert sum(meta["category_totals"].values()) == tup_standard.n_triangles


def test_the_convention_is_recorded_and_is_recoverable_in_one_call(
    tup_standard, tup_other
):
    assert strand_identity_of(tup_standard) == "sheet_id_and_serial"
    assert sse_annotation_summary(tup_standard)["assignment"]["strand_identity"] == (
        "sheet_id_and_serial"
    )
    assert strand_identity_of(tup_other) == OTHER_CONVENTION


def test_asking_an_unannotated_result_for_its_convention_refuses(plain):
    with pytest.raises(IncompatibleSSEStratification, match="no SSE annotation"):
        strand_identity_of(plain)


def test_neither_identity_can_catch_a_mixed_convention(tup_standard):
    """The reason a dedicated guard is needed at all.

    This is not a gap being papered over -- it is the deliberate consequence of
    keeping the convention out of both identities, and it is what makes the
    ordinary-key comparison below legitimate.
    """
    assert "strand_identity" not in json.dumps(tup_standard.key_space.to_dict())
    assert "strand_identity" not in json.dumps(tup_standard.sampling.to_dict())


def test_stratified_analysis_across_conventions_is_refused_by_name(
    tup_other, tup_standard
):
    with pytest.raises(IncompatibleSSEStratification) as excinfo:
        require_common_strand_identity([tup_standard, tup_other])
    message = str(excinfo.value)
    assert "legacy_serial_only" in message and "sheet_id_and_serial" in message
    assert "1TUP" in message
    # The refusal must say what IS still comparable, or it teaches the wrong
    # lesson -- that these two results cannot be compared at all.
    assert "result.frequencies" in message


def test_the_guard_passes_and_returns_the_shared_convention(tup_standard, annotated):
    """Same convention, different structures: comparable, and it says which."""
    assert require_common_strand_identity([tup_standard]) == "sheet_id_and_serial"
    assert require_common_strand_identity([tup_standard, annotated]) == (
        "sheet_id_and_serial"
    )


def test_the_guard_refuses_an_unannotated_result_in_the_set(tup_standard, plain):
    with pytest.raises(IncompatibleSSEStratification, match="no SSE annotation"):
        require_common_strand_identity([tup_standard, plain])


def test_the_guard_needs_at_least_one_result(annotated):
    with pytest.raises(ValueError, match="no results"):
        require_common_strand_identity([])


def test_a_context_label_is_carried_into_the_refusal(tup_other, tup_standard):
    with pytest.raises(IncompatibleSSEStratification, match="category feature matrix"):
        require_common_strand_identity(
            [tup_standard, tup_other], context="category feature matrix"
        )


def test_ordinary_key_analysis_across_conventions_stays_allowed(
    tup_other, tup_standard
):
    """The other half, and it must not regress into over-strictness.

    The convention changes no key, so a key-frequency comparison across the two
    is correct science. `AnalysisSample.from_result` makes it unavoidable rather
    than merely permitted: it sums the categories away, so nothing downstream
    can even see which convention produced the counts.
    """
    a = AnalysisSample.from_result(tup_other, sample_id="other")
    b = AnalysisSample.from_result(tup_standard, sample_id="standard")

    assert a.counts == b.counts
    assert generalized_jaccard(a.counts, b.counts) == pytest.approx(1.0)
    assert "strand_identity" not in json.dumps(a.provenance())
    assert not any("sse" in str(v).lower() for v in a.annotations.values())

    from biotsr.analysis.samples import require_common_key_space, require_common_sampling

    require_common_key_space([a, b])
    require_common_sampling([a, b])
    matrix = feature_matrix([a, b])
    assert matrix.rows[0] == matrix.rows[1]
