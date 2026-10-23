"""Amino-acid atomic TSR end to end: modifiers, standard gate, determinism.

The architectural claim under test is that the atomic families reuse the
protein C-alpha core (geometry, encoding and modifiers) rather than becoming a
second TSR implementation.  Concretely: the
same canonical ordering, the same ``GeometryResult``, the same encoder, the same
mirror and size-filter modifiers, the same result model and the same
serialization determinism — configured differently.
"""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.api import analyze_entities, list_residue_targets
from biotsr.domain.result import DEFAULT_CATEGORY
from biotsr.io import native
from biotsr.modifiers.mirror import MIRROR_CATEGORY
from biotsr.profiles.atomic import STANDARD_AMINOACID_ATOMIC_V1

from ..conftest import FIXTURE_PDB, _require

CUTOFF = 4.0
#: 1UBQ chain A: six side-chain atoms across six residues.
SELECTION = "LEU_8_CD1;LEU_43_CD2;ILE_44_CD1;LYS_48_NZ;HIS_68_NE2;VAL_70_CG1"


@pytest.fixture(scope="module")
def ubq_path():
    path = os.path.join(FIXTURE_PDB, "1UBQ.pdb")
    _require(path)
    return path


@pytest.fixture(scope="module")
def leu15(ubq_path):
    return biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU",
    )


# ---------------------------------------------------------------------------
# Regular amino-acid atomic TSR
# ---------------------------------------------------------------------------


def test_one_residue_end_to_end(leu15):
    assert leu15.counts.n_points == 8
    assert leu15.counts.n_candidates == 56  # C(8,3)
    assert leu15.counts.n_valid_geometry == 56
    assert leu15.counts.n_degenerate == 0
    assert leu15.counts.n_filtered_by_size == 0
    assert leu15.counts.n_triangles == 56
    assert leu15.counts.ledger_balances
    assert leu15.method == "amino_acid_atomic"


def test_the_standard_profile_uses_the_atomic_regime(leu15):
    assert leu15.key_space.distance_scheme.radix == 18
    assert leu15.key_space.theta_scheme.radix == 29
    assert leu15.key_space.label_scheme.radix == 112
    assert leu15.key_space.angle_convention.startswith("pi:")
    assert leu15.key_space.sqrt_convention == "math_sqrt"


def test_a_whole_chain_runs_per_residue(ubq_path):
    results = analyze_entities(ubq_path, method="amino_acid_atomic", chain="A")
    assert len(results) == 75, "76 residues offered, Gly76 declined by completeness"
    assert "15_LEU" in results
    assert "76_GLY" not in results
    assert all(r.counts.ledger_balances for r in results.values())


def test_rejected_residues_are_inspectable_on_request(ubq_path):
    results = analyze_entities(
        ubq_path, method="amino_acid_atomic", chain="A",
        keep_triplets=False, include_rejected=True,
    )
    assert len(results) == 76
    gly76 = results["76_GLY"]
    assert gly76.counts.n_triangles == 0
    assert [s.reason for s in gly76.skipped if s.stage == "completeness"] == [
        "atom_count_mismatch"
    ]


def test_targets_are_offered_before_they_are_admitted(ubq_path):
    targets = list_residue_targets(ubq_path, chain="A")
    assert len(targets) == 76
    assert ("76", "GLY") in targets, "offered; the policy declines it later"


def test_a_missing_target_is_an_error_not_a_silent_empty_result(ubq_path):
    with pytest.raises(ValueError, match="needs a target residue"):
        biotsr.analyze(ubq_path, method="amino_acid_atomic", chain="A")


def test_the_calpha_method_refuses_atomic_arguments(ubq_path):
    with pytest.raises(biotsr.InvalidRequest, match="does not use residue"):
        biotsr.analyze(ubq_path, method="protein", chain="A", residue="15")


# ---------------------------------------------------------------------------
# Atom selection
# ---------------------------------------------------------------------------


def test_atom_selection_end_to_end(ubq_path):
    result = biotsr.analyze(
        ubq_path, method="amino_acid_atom_selection", chain="A", atoms=SELECTION,
    )
    assert result.counts.n_points == 6
    assert result.counts.n_triangles == 20  # C(6,3)
    assert result.selection.atoms_missing == ()
    assert result.sampling.scope == "pooled"
    assert result.sampling.completeness_policy is None


def test_atom_selection_spans_several_residues_in_one_calculation(ubq_path):
    result = biotsr.analyze(
        ubq_path, method="amino_acid_atom_selection", chain="A", atoms=SELECTION,
    )
    residues = {t.identities[0].residue_seq for t in result.triplets}
    residues |= {t.identities[1].residue_seq for t in result.triplets}
    assert len(residues) > 1, "the whole point of the pooled scope"


def test_atom_selection_requires_a_selection(ubq_path):
    with pytest.raises(ValueError, match="needs an atom selection"):
        biotsr.analyze(ubq_path, method="amino_acid_atom_selection", chain="A")


def test_an_entirely_unusable_selection_is_refused_loudly(ubq_path):
    with pytest.raises(ValueError, match="'nonsense' is not RESIDUE_ID_ATOM"):
        biotsr.analyze(
            ubq_path, method="amino_acid_atom_selection", chain="A", atoms="nonsense",
        )


def test_a_partly_malformed_selection_is_refused_naming_the_token(ubq_path):
    """Since 1.0.0rc4 a malformed token is refused, not skipped with a warning:
    a selection that silently loses an atom is a different calculation."""
    with pytest.raises(ValueError, match="'junk' is not RESIDUE_ID_ATOM"):
        biotsr.analyze(
            ubq_path, method="amino_acid_atom_selection", chain="A",
            atoms="LEU_8_CD1;junk;LEU_43_CD2;ILE_44_CD1",
        )
    result = biotsr.analyze(
        ubq_path, method="amino_acid_atom_selection", chain="A",
        atoms="LEU_8_CD1;LEU_43_CD2;ILE_44_CD1",
    )
    assert result.counts.n_points == 3


# ---------------------------------------------------------------------------
# Mirror
# ---------------------------------------------------------------------------


def test_mirror_off_never_produces_a_negative_key(ubq_path):
    results = analyze_entities(ubq_path, method="amino_acid_atomic", chain="A")
    for label, r in results.items():
        assert all(k >= 0 for k in r.keys), label
        assert set(r.category_totals()) == {DEFAULT_CATEGORY}


def test_mirror_on_signs_by_the_pre_fold_angle(ubq_path):
    results = analyze_entities(
        ubq_path, method="amino_acid_atomic", chain="A", mirror=True
    )
    checked = 0
    for r in results.values():
        for t in r.triplets:
            if t.theta_raw > 90:
                assert t.key < 0 and t.category == MIRROR_CATEGORY
            else:
                assert t.key >= 0 and t.category == DEFAULT_CATEGORY
            checked += 1
    assert checked == 4860


def test_mirror_is_a_pure_sign_multiplier(ubq_path, leu15):
    mirrored = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", mirror=True,
    )
    assert [abs(t.key) for t in mirrored.triplets] == [t.key for t in leu15.triplets]


def test_mirror_changes_the_key_space_not_the_triangle_count(ubq_path, leu15):
    mirrored = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", mirror=True,
    )
    assert mirrored.counts.n_triangles == leu15.counts.n_triangles
    assert mirrored.key_space.fingerprint != leu15.key_space.fingerprint
    assert mirrored.sampling.fingerprint == leu15.sampling.fingerprint


def test_both_angles_are_retained(leu15):
    for t in leu15.triplets:
        assert 0 <= t.theta <= 90
        assert 0 <= t.theta_raw <= 180
        expected = t.theta_raw if t.theta_raw <= 90 else abs(180 - t.theta_raw)
        assert t.theta == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
# Size filtering
# ---------------------------------------------------------------------------


def test_the_standard_profile_carries_no_size_filter_of_its_own():
    """Recorded as a scope fact, not assumed away.

    Neither ``TSR_AMINOACID_PACKAGE`` nor ``WebsiteTSR`` has a size filter, so
    the modifier is offered on the standard profile only, and the profile
    itself applies none unless it is asked for.
    """
    assert STANDARD_AMINOACID_ATOMIC_V1.size_filter is None


def test_size_filter_below_at_and_above_the_boundary(ubq_path, leu15):
    """Strict `<`: a triangle exactly at the cutoff is excluded."""
    distances = sorted({t.max_distance for t in leu15.triplets})
    exact = distances[len(distances) // 2]

    at = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", size_filter=exact,
    )
    assert all(t.max_distance < exact for t in at.triplets)
    assert not any(t.max_distance == exact for t in at.triplets)

    above = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", size_filter=exact * (1 + 1e-12),
    )
    assert above.counts.n_triangles > at.counts.n_triangles

    below = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", size_filter=distances[0],
    )
    assert below.counts.n_triangles == 0


def test_size_filter_diagnostics_are_separate_from_degeneracy(ubq_path):
    filtered = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", size_filter=CUTOFF,
    )
    assert filtered.counts.n_candidates == 56
    assert filtered.counts.n_valid_geometry == 56
    assert filtered.counts.n_degenerate == 0
    assert filtered.counts.n_filtered_by_size > 0
    assert filtered.counts.ledger_balances
    assert "size_filtered" not in {s.reason for s in filtered.skipped}


def test_size_filter_composes_with_mirror(ubq_path):
    both = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", mirror=True, size_filter=CUTOFF,
    )
    only_filter = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", size_filter=CUTOFF,
    )
    assert both.counts.n_triangles == only_filter.counts.n_triangles
    assert any(k < 0 for k in both.keys)
    assert all(k >= 0 for k in only_filter.keys)
    assert list(both.annotations["modifiers"]) == ["mirror", "size_filter"]


def test_size_filtering_works_for_atom_selection_too(ubq_path):
    result = biotsr.analyze(
        ubq_path, method="amino_acid_atom_selection", chain="A",
        atoms=SELECTION, size_filter=10.0,
    )
    assert result.counts.n_filtered_by_size > 0
    assert all(t.max_distance < 10.0 for t in result.triplets)
    assert result.counts.ledger_balances


# ---------------------------------------------------------------------------
# Key-space separation
# ---------------------------------------------------------------------------


def test_an_atomic_result_never_shares_a_key_space_with_a_calpha_one(ubq_path, leu15):
    from biotsr.domain.keyspace import IncompatibleKeySpace

    calpha = biotsr.analyze(ubq_path, method="protein", chain="A")
    with pytest.raises(IncompatibleKeySpace):
        leu15.key_space.require_compatible(calpha.key_space)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def test_entity_scoped_native_prefix(leu15):
    assert native.default_prefix(leu15) == "1UBQ_A_amino_acid_atomic_LEU_15"


def test_the_metadata_distinguishes_atomic_from_calpha_output(ubq_path, leu15):
    """The file names no longer carry bin counts; the key space does."""
    calpha = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert native.default_prefix(calpha) == "1UBQ_A_protein"
    assert calpha.metadata()["key_space"]["distance_scheme"]["radix"] == 35
    assert leu15.metadata()["key_space"]["distance_scheme"]["radix"] == 18


def test_native_triplets_identify_the_three_physical_atoms(leu15, tmp_path):
    path = str(tmp_path / "t")
    native.write_triplets(leu15, path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    for n in (1, 2, 3):
        assert row[header.index("chain%d" % n)] == "A"
        assert row[header.index("residue%d" % n)] == "LEU"
        assert row[header.index("seq%d" % n)] == "15"
        assert row[header.index("atom%d" % n)] in {
            "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2",
        }
    for column in ("theta_raw", "theta", "theta_class", "distance_class",
                   "max_distance", "category", "key"):
        assert column in header


def test_mirror_output_carries_the_category_alongside_the_signed_key(ubq_path, tmp_path):
    mirrored = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU", mirror=True,
    )
    path = str(tmp_path / "k")
    native.write_keys(mirrored, path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]
    assert header == ["key", "category", "frequency"]
    assert {r[1] for r in rows} == {DEFAULT_CATEGORY, MIRROR_CATEGORY}
    for key, category, _freq in rows:
        assert (int(key) < 0) == (category == MIRROR_CATEGORY)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"mirror": True},
        {"size_filter": CUTOFF},
        {"mirror": True, "size_filter": CUTOFF},
    ],
)
def test_repeated_runs_are_byte_identical(ubq_path, kwargs):
    def once():
        r = biotsr.analyze(
            ubq_path, method="amino_acid_atomic", chain="A",
            residue="15", residue_name="LEU", **kwargs
        )
        return (
            tuple(r.key_order),
            "".join(native.format_native_triplet(t) for t in r.triplets),
            dict(r.frequencies),
        )

    assert once() == once()


def test_streaming_matches_the_retained_run(ubq_path):
    seen = []
    streamed = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU",
        keep_triplets=False, triplet_sink=seen.append,
    )
    retained = biotsr.analyze(
        ubq_path, method="amino_acid_atomic", chain="A",
        residue="15", residue_name="LEU",
    )
    assert streamed.triplets is None
    assert [t.key for t in seen] == [t.key for t in retained.triplets]
    assert streamed.frequencies == retained.frequencies
