"""One analysis layer, every method family -- and a hard stop where the keys differ.

The point of the analysis layer is that ``biotsr.analysis`` does not know, and
cannot ask, which engine produced a result.  These tests feed it real results
from unrelated families -- ordinary protein Ca, Ca with AA grouping,
mirror-signed, size-filtered, amino-acid atomic, drug, carbohydrate, nucleotide
and protein-protein cross -- and assert two things: that compatible ones flow
through the *same* code with no branch on ``result.method``, and that
incompatible key spaces are refused rather than merged.
"""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.analysis import (
    AnalysisSample,
    average_linkage,
    cluster,
    common_key_intersection,
    distance_matrix,
    feature_matrix,
    key_occurrence_percent,
    similarity_matrix,
)
from biotsr.domain.keyspace import IncompatibleKeySpace, IncompatibleSampling

from ..conftest import structure as tracked_structure

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = HERE  # the package project is the repository root


def _fixture(pdb_id):
    return tracked_structure(pdb_id)


@pytest.fixture(scope="module")
def brs_path():
    """1BRS -- barnase/barstar, the protein-protein cross fixture."""
    return _fixture("1BRS")


@pytest.fixture(scope="module")
def ac0_cross_path():
    """2AC0 -- the nucleotide-protein cross fixture."""
    return _fixture("2AC0")


def _sample(result, name, label=None):
    return AnalysisSample.from_result(result, sample_id=name, label=label)


# ---------------------------------------------------------------------------
# Same key space, different sampling rule: refused unless explicitly allowed
# ---------------------------------------------------------------------------


def test_the_same_chain_at_two_contact_cutoffs_needs_the_opt_in(brs_path):
    """A cutoff change keeps every key's meaning but draws the frequencies
    differently: refused by default, compared and recorded when allowed."""
    samples = []
    for cutoff in (3.5, 5.0):
        result = biotsr.analyze_cross(
            brs_path, method="protein_protein_cross", chain="A", chain_2="D",
            contact_cutoff=cutoff, keep_triplets=False,
        )
        samples.append(_sample(result, "1BRS_A_D_cut%.1f" % cutoff))
    with pytest.raises(IncompatibleSampling, match="contact_cutoff"):
        feature_matrix(samples)
    fm = feature_matrix(samples, allow_mixed_processing=True)
    assert fm.n_samples == 2
    assert fm.comparison["processing_uniform"] is False
    assert any("contact_cutoff" in w for w in fm.warnings)
    dm = distance_matrix(fm)
    assert 0.0 <= dm.value(0, 1) <= 1.0
    assert dm.warnings == fm.warnings
    assert cluster(dm, 1).warnings[: len(fm.warnings)] == fm.warnings


def test_size_filtered_and_unfiltered_calpha_need_the_opt_in(ubq_pdb):
    """Size filtering changes which triangles were drawn, never what a key means."""
    plain = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    filtered = biotsr.analyze(
        ubq_pdb, method="protein", chain="A", size_filter=25.0, keep_triplets=False
    )
    assert plain.key_space.fingerprint == filtered.key_space.fingerprint
    pair = [_sample(plain, "1UBQ_plain"), _sample(filtered, "1UBQ_sf25")]
    with pytest.raises(IncompatibleSampling, match="size_filter"):
        feature_matrix(pair)
    fm = feature_matrix(pair, allow_mixed_processing=True)
    assert fm.comparison["processing_uniform"] is False
    assert similarity_matrix(fm).value(0, 1) > 0.0


# ---------------------------------------------------------------------------
# Different key space: the case that must be REFUSED
# ---------------------------------------------------------------------------


def test_grouped_and_ungrouped_calpha_cannot_be_compared(ubq_pdb):
    """AA grouping changes the label scheme and its radix, so key 5 means two things."""
    plain = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    grouped = biotsr.analyze(
        ubq_pdb, method="protein", chain="A", aa_grouping=True, keep_triplets=False
    )
    assert plain.key_space.fingerprint != grouped.key_space.fingerprint
    with pytest.raises(IncompatibleKeySpace, match="label_scheme"):
        feature_matrix([_sample(plain, "plain"), _sample(grouped, "grouped")])


def test_mirror_signed_and_unsigned_calpha_cannot_be_compared(ubq_pdb):
    plain = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    mirrored = biotsr.analyze(
        ubq_pdb, method="protein", chain="A", mirror=True, keep_triplets=False
    )
    assert plain.key_space.fingerprint != mirrored.key_space.fingerprint
    with pytest.raises(IncompatibleKeySpace, match="mirror_convention"):
        feature_matrix([_sample(plain, "plain"), _sample(mirrored, "mirror")])


def test_calpha_and_atomic_cannot_be_compared(ubq_pdb):
    """Residue labels and atom labels are different alphabets in the same integers."""
    calpha = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    atomic = biotsr.analyze(
        ubq_pdb, method="amino_acid_atomic", chain="A",
        residue="57", residue_name="SER", keep_triplets=False,
    )
    with pytest.raises(IncompatibleKeySpace):
        feature_matrix([_sample(calpha, "ca"), _sample(atomic, "atomic")])


def test_every_analysis_entry_point_refuses_a_mixed_key_space(ubq_pdb):
    plain = _sample(
        biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False), "plain"
    )
    grouped = _sample(
        biotsr.analyze(
            ubq_pdb, method="protein", chain="A", aa_grouping=True, keep_triplets=False
        ),
        "grouped",
    )
    for call in (
        lambda: feature_matrix([plain, grouped]),
        lambda: common_key_intersection([plain, grouped]),
        lambda: key_occurrence_percent([plain, grouped]),
    ):
        with pytest.raises(IncompatibleKeySpace):
            call()


# ---------------------------------------------------------------------------
# One implementation, several families
# ---------------------------------------------------------------------------


def test_atomic_families_share_one_analysis_path(mapk_path, lf9_path, bdna_path):
    """Drug, carbohydrate, nucleotide and amino-acid atomic all reduce to samples.

    Each family keeps its own key space, so each gets its own matrix -- but the
    *code* that builds those matrices, the similarity and the clustering is one
    implementation with no ``if method ==`` anywhere in it.
    """
    families = {
        "drug": biotsr.analyze(
            mapk_path, method="drug", chain="A", entity_name="ESW", entity_id="404",
            keep_triplets=False,
        ),
        "carbohydrate": biotsr.analyze(
            lf9_path, method="carbohydrate", chain="D", entity_name="NAG", entity_id="1",
            keep_triplets=False,
        ),
        "nucleotide": biotsr.analyze(
            bdna_path, method="nucleotide", chain="A", entity_name="DG", entity_id="2",
            keep_triplets=False,
        ),
    }
    for name, result in families.items():
        sample = _sample(result, name)
        assert sample.n_keys > 0, name
        fm = feature_matrix([sample])
        assert fm.shape == (1, sample.n_keys)
        assert fm.provenance()["key_space_fingerprint"] == result.key_space.fingerprint
        assert fm.provenance()["comparison"]["processing_uniform"] is True
        # the layer records the method as provenance and never branches on it
        assert sample.annotations["method"] == result.method


def test_several_nucleotide_residues_cluster_through_the_generic_path(bdna_path):
    """A real multi-sample run: one structure, several residues, one analysis."""
    targets = [("2", "DG"), ("3", "DC"), ("4", "DG"), ("5", "DA"), ("7", "DT")]
    samples = []
    for seq, name in targets:
        result = biotsr.analyze(
            bdna_path, method="nucleotide", chain="A", entity_name=name,
            entity_id=seq, keep_triplets=False,
        )
        samples.append(_sample(result, "1BNA_A_%s_%s" % (seq, name), label=name))
    fm = feature_matrix(samples)
    assert fm.n_samples == 5
    dm = distance_matrix(fm)
    assignment = cluster(dm, 2)
    assert assignment.n_clusters == 2
    assert sorted(assignment.ordered_sample_ids()) == sorted(fm.sample_ids)
    assert assignment.provenance()["linkage_method"] == "average"


def test_a_cross_family_result_flows_through_unchanged(brs_path, ac0_cross_path):
    """Cross results are samples like any other; only their key space differs."""
    pp = biotsr.analyze_cross(
        brs_path, method="protein_protein_cross", chain="A", chain_2="D",
        keep_triplets=False,
    )
    np_ = biotsr.analyze_cross(
        ac0_cross_path, method="nucleotide_protein_cross", chain="E", chain_2="B",
        entity_name="DG", entity_id="8", keep_triplets=False,
    )
    for result, name in ((pp, "1BRS_A_D"), (np_, "1TUP_A_E")):
        fm = feature_matrix([_sample(result, name)])
        assert fm.n_samples == 1
        assert fm.nnz == len(result.frequencies)
    if pp.key_space.fingerprint != np_.key_space.fingerprint:
        with pytest.raises(IncompatibleKeySpace):
            feature_matrix([_sample(pp, "pp"), _sample(np_, "np")])


def test_analysis_never_needs_the_triplets_or_the_structure(ubq_pdb):
    """``keep_triplets=False`` loses nothing the analysis layer uses."""
    result = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    assert result.triplets is None
    sample = _sample(result, "1UBQ_A")
    fm = feature_matrix([sample, _sample(result, "1UBQ_A_copy")])
    assert similarity_matrix(fm).value(0, 1) == pytest.approx(1.0)
    assert average_linkage(distance_matrix(fm)).heights == pytest.approx((0.0,))


def test_a_sample_can_be_built_from_a_bare_mapping_with_no_engine_at_all(ubq_pdb):
    """Nothing in the layer requires a TSRResult; the identities are what it needs."""
    result = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    from_result = AnalysisSample.from_result(result, sample_id="x")
    hand_made = AnalysisSample(
        sample_id="x",
        counts=dict(result.frequencies),
        key_space=result.key_space,
        sampling=result.sampling,
    )
    assert from_result.counts == hand_made.counts
    assert feature_matrix([hand_made]).columns == feature_matrix([from_result]).columns
