"""The round trip, once per representative family.

A round trip that works on protein Cα and nowhere else is not a round trip.  This
module runs one representative calculation per family the package supports,
saves each, loads it back, and asserts the same eight properties every time.

The list is deliberately the *whole* list of supported families rather
than a sample: the families differ in exactly the ways a serializer can get
wrong.  Entity-scoped methods put a target in the filename; cross methods put a
pipe in the chain and fill all 21 sampling fields; SSE fills eighteen category
buckets where every other method fills one; grouping changes the label radix and
mirror signs keys.
"""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.io import native, native_read

from ..conftest import structure as tracked_structure

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = HERE  # the package project is the repository root


def fixture(pdb_id):
    return tracked_structure(pdb_id)


#: One case per family.  ``cross`` marks the ones that go through
#: ``analyze_cross`` rather than ``analyze``; everything else is identical.
#:
#: The cross cutoffs are narrowed exactly as ``test_cross_families`` narrows
#: them: 1BRS A/D at its historical 5 A emits 1.2 million triangles, and this
#: module asserts round-trip fidelity, not scale.
FAMILIES = {
    "protein_calpha": dict(
        source="1UBQ", method="protein", chain="A",
    ),
    "protein_calpha_grouped": dict(
        source="1UBQ", method="protein", chain="A", aa_grouping=True,
    ),
    "protein_calpha_mirror": dict(
        source="1UBQ", method="protein", chain="A", mirror=True,
    ),
    "protein_calpha_size_filtered": dict(
        source="1UBQ", method="protein", chain="A", size_filter=14.0,
    ),
    "protein_calpha_sse": dict(
        source="1TUP", method="protein", chain="A", sse=True,
    ),
    "amino_acid_atomic": dict(
        source="1UBQ", method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU",
    ),
    "aa_selection": dict(
        source="1UBQ", method="aa_selection", chain="A", residues="35-40;54;60-64",
    ),
    "drug": dict(
        source="3PTB", method="drug", chain="A", entity_name="BEN", entity_id="1",
    ),
    "carbohydrate": dict(
        source="5WT9", method="carbohydrate", chain="A", entity_name="NAG",
        entity_id="1",
    ),
    "nucleotide": dict(
        source="1BNA", method="nucleotide", chain="A", entity_name="DC",
        entity_id="3",
    ),
    "drug_protein_cross": dict(
        cross=True, source="3PTB", method="drug_protein_cross", chain="A",
        chain_2="A", entity_name="BEN", entity_id="1",
    ),
    "carbohydrate_protein_cross": dict(
        cross=True, source="7LF9", method="carbohydrate_protein_cross", chain="C",
        chain_2="A", entity_name="NAG", entity_id="1",
    ),
    "nucleotide_protein_cross": dict(
        cross=True, source="2AC0", method="nucleotide_protein_cross", chain="E",
        chain_2="B", entity_name="DG", entity_id="8",
    ),
    "protein_protein_cross_all_atom": dict(
        cross=True, source="1BRS", method="protein_protein_cross", chain="A",
        chain_2="D", contact_cutoff=3.2,
    ),
    "protein_protein_cross_calpha_interface": dict(
        cross=True, source="1BRS", method="protein_protein_calpha_cross", chain="A",
        chain_2="D", contact_cutoff=3.2,
    ),
}


def run(name):
    spec = dict(FAMILIES[name])
    is_cross = spec.pop("cross", False)
    source = fixture(spec.pop("source"))
    call = biotsr.analyze_cross if is_cross else biotsr.analyze
    return call(source, keep_triplets=False, **spec)


@pytest.fixture(scope="module")
def computed():
    """Every family computed once; the round-trip assertions are cheap, these are not."""
    return {name: run(name) for name in FAMILIES}


def _round_trip(result, tmp_path, name):
    written = native.save(result, str(tmp_path / name), triplets=False)
    return native_read.load_result(written["keys"]), written


@pytest.mark.parametrize("name", list(FAMILIES))
def test_every_key_and_category_survives(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert {k: dict(v.by_category) for k, v in loaded.keys.items()} == {
        k: dict(v.by_category) for k, v in original.keys.items()
    }


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_key_space_fingerprint_is_unchanged(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert loaded.key_space.fingerprint == original.key_space.fingerprint
    assert loaded.key_space == original.key_space


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_sampling_fingerprint_is_unchanged(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert loaded.sampling.fingerprint == original.sampling.fingerprint
    assert loaded.sampling == original.sampling


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_profile_provenance_is_unchanged(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert loaded.annotations["profile"] == original.annotations["profile"]
    assert (
        loaded.annotations["legacy_compatibility_profile"]
        == original.annotations["legacy_compatibility_profile"]
    )


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_ledger_and_the_selection_are_unchanged(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert loaded.counts == original.counts
    assert loaded.selection == original.selection
    assert loaded.source == original.source


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_category_totals_are_unchanged(name, computed, tmp_path):
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    assert loaded.category_totals() == original.category_totals()


@pytest.mark.parametrize("name", list(FAMILIES))
def test_re_serialization_is_byte_identical(name, computed, tmp_path):
    original = computed[name]
    once = native.save(original, str(tmp_path / ("%s_a" % name)), triplets=False)
    loaded = native_read.load_result(once["keys"])
    twice = native.save(loaded, str(tmp_path / ("%s_b" % name)), triplets=False)
    for role in ("keys", "metadata"):
        with open(once[role], "rb") as fh:
            first = fh.read()
        with open(twice[role], "rb") as fh:
            second = fh.read()
        assert first == second, "%s: %s differs after a round trip" % (name, role)


@pytest.mark.parametrize("name", list(FAMILIES))
def test_the_loaded_result_reduces_to_the_same_analysis_sample(name, computed, tmp_path):
    """The operational test: what the analysis layer sees is identical.

    Every evidenced analysis in the package is defined over an
    ``AnalysisSample``, which is a key-frequency vector plus the two identities.
    If those three agree, every measure, every feature matrix, every common-key
    operation and every clustering agrees -- without this module having to
    re-test each one per family.
    """
    from biotsr.analysis import AnalysisSample

    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, name)
    a = AnalysisSample.from_result(original)
    b = AnalysisSample.from_result(loaded)
    assert dict(a.counts) == dict(b.counts)
    assert a.key_space.fingerprint == b.key_space.fingerprint
    assert a.sampling.fingerprint == b.sampling.fingerprint
    assert a.sample_id == b.sample_id


# ---------------------------------------------------------------------------
# family-specific properties the shared assertions cannot reach
# ---------------------------------------------------------------------------


def test_an_entity_scoped_result_keeps_its_target_and_its_prefix(computed, tmp_path):
    """The residue is in the prefix, and the prefix must round-trip.

    ``1UBQ_A_amino_acid_atomic_LEU_15`` -- structure, chain, method, residue.  A
    loaded result that lost its selection record would re-serialize under a
    name without the residue and collide with every other residue of the chain.
    """
    original = computed["amino_acid_atomic"]
    loaded, written = _round_trip(original, tmp_path, "entity")
    assert os.path.basename(written["keys"]) == "1UBQ_A_amino_acid_atomic_LEU_15.keys"
    assert loaded.annotations["target"] == original.annotations["target"]
    assert native.default_prefix(loaded) == native.default_prefix(original)


def test_an_sse_result_keeps_its_strand_identity_convention(computed, tmp_path):
    """ADR-063: the convention is in NEITHER fingerprint, so nothing else carries it."""
    from biotsr.sse import require_common_strand_identity, strand_identity_of

    original = computed["protein_calpha_sse"]
    loaded, _ = _round_trip(original, tmp_path, "sse")
    assert strand_identity_of(loaded) == strand_identity_of(original)
    # The SSE strand-identity guard still works on a loaded result, which is the point:
    # the safeguard must not be weakened by the result having come off disk.
    require_common_strand_identity([original, loaded])


def test_an_sse_result_keeps_its_eighteen_category_buckets(computed, tmp_path):
    from biotsr.sse import category_totals_in_order, check_conservation

    original = computed["protein_calpha_sse"]
    loaded, _ = _round_trip(original, tmp_path, "sse_cats")
    assert category_totals_in_order(loaded) == category_totals_in_order(original)
    assert len(category_totals_in_order(loaded)) == 18
    check_conservation(loaded)


def test_an_sse_result_keeps_its_record_availability(computed, tmp_path):
    """``has_records`` distinguishes "all coil" from "nobody looked"."""
    original = computed["protein_calpha_sse"]
    loaded, _ = _round_trip(original, tmp_path, "sse_records")
    assert loaded.annotations["sse"] == original.annotations["sse"]
    assert loaded.annotations["sse"]["has_records"] is True
    assert loaded.annotations["sse_coverage"] == original.annotations["sse_coverage"]


def test_a_mirror_result_keeps_its_negative_keys_and_both_buckets(computed, tmp_path):
    original = computed["protein_calpha_mirror"]
    loaded, _ = _round_trip(original, tmp_path, "mirror")
    assert any(k < 0 for k in loaded.keys)
    assert set(loaded.category_totals()) == {"regular", "mirror"}
    assert loaded.category_totals() == original.category_totals()


def test_a_cross_result_keeps_its_two_sides(computed, tmp_path):
    original = computed["protein_protein_cross_all_atom"]
    loaded, _ = _round_trip(original, tmp_path, "cross")
    assert loaded.chain == original.chain == "A|D"
    assert loaded.sampling.side_a_chain == "A"
    assert loaded.sampling.side_b_chain == "D"
    assert loaded.sampling.contact_cutoff == 3.2
    assert loaded.annotations["cross_pools"] == original.annotations["cross_pools"]


def test_a_size_filtered_result_keeps_its_cutoff_as_a_float(computed, tmp_path):
    original = computed["protein_calpha_size_filtered"]
    loaded, _ = _round_trip(original, tmp_path, "sized")
    assert loaded.sampling.size_filter == 14.0
    assert isinstance(loaded.sampling.size_filter, float)
    assert loaded.counts.n_filtered_by_size == original.counts.n_filtered_by_size


@pytest.mark.parametrize("name", [n for n in FAMILIES if FAMILIES[n].get("cross")])
def test_a_cross_result_keeps_its_numerical_precision(name, computed, tmp_path):
    """``numerical_precision`` is key-affecting (ADR-053) and must survive."""
    original = computed[name]
    loaded, _ = _round_trip(original, tmp_path, "precision_" + name)
    assert original.key_space.numerical_precision == "float64"
    assert loaded.key_space.numerical_precision == (
        original.key_space.numerical_precision
    )


def test_a_grouped_result_keeps_its_grouped_scheme(computed, tmp_path):
    original = computed["protein_calpha_grouped"]
    loaded, _ = _round_trip(original, tmp_path, "grouped")
    assert loaded.key_space.label_scheme.name == "amino_acid_group_23"
    assert loaded.key_space.label_scheme.radix == 30
    assert loaded.key_space.fingerprint == original.key_space.fingerprint


def test_results_from_incompatible_key_spaces_are_still_refused(computed, tmp_path):
    """Loading must not launder an incompatibility."""
    from biotsr.analysis import AnalysisSample, require_common_key_space

    grouped, _ = _round_trip(computed["protein_calpha_grouped"], tmp_path, "g")
    plain, _ = _round_trip(computed["protein_calpha"], tmp_path, "p")
    with pytest.raises(biotsr.IncompatibleKeySpace):
        require_common_key_space(
            [
                AnalysisSample.from_result(plain, sample_id="plain"),
                AnalysisSample.from_result(grouped, sample_id="grouped"),
            ]
        )
