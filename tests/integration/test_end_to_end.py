"""End-to-end behaviour of the public API and the shared core loop."""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.methods.engine import run_core
from biotsr.methods.registry import get_method, list_methods

from ..conftest import FIXTURE_PDB, _require


@pytest.fixture(scope="module")
def ubq_path():
    path = os.path.join(FIXTURE_PDB, "1UBQ.pdb")
    _require(path)
    return path


# ---------------------------------------------------------------------------
# The one-call path
# ---------------------------------------------------------------------------


def test_simple_things_are_one_call(ubq_path):
    result = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert result.n_triangles == 70300
    assert result.n_unique_keys == 45349  # standard profile
    assert result.method == "protein"
    assert result.chain == "A"


def test_the_default_profile_is_the_standard_one(ubq_path):
    result = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert result.annotations["profile"] == "standard_protein_calpha_v1"
    assert result.annotations["legacy_compatibility_profile"] is False
    assert result.key_space.angle_convention.startswith("pi")


def test_the_default_path_uses_the_standard_arithmetic(ubq_path):
    """Every default path uses math.pi and the standard 35-class Ca scheme."""
    import math

    result = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert result.key_space.angle_convention == "pi:%r" % math.pi
    assert result.key_space.sqrt_convention == "pow_half"
    assert result.key_space.distance_scheme.name == "distance_calpha_35"
    assert biotsr.STANDARD_PROTEIN_CALPHA_V1.angle_convention.radians_per_half_turn is math.pi


def test_result_is_self_describing(ubq_path):
    result = biotsr.analyze(ubq_path, method="protein", chain="A")
    meta = result.metadata()
    assert meta["structure_id"] == "1UBQ"
    assert meta["chain"] == "A"
    assert meta["source"]["parser_backend"] == "pdb_legacy"
    assert meta["source"]["sha256"]
    assert meta["selection"]["n_points"] == 76
    assert meta["counts"]["n_candidates"] == 70300
    assert meta["key_space_fingerprint"] != meta["sampling_fingerprint"]


def test_structure_id_can_be_overridden_with_a_long_identifier(ubq_path):
    result = biotsr.analyze(
        ubq_path, method="protein", chain="A", structure_id="PDB_0000001UBQ"
    )
    assert result.structure_id.value == "PDB_0000001UBQ"
    from biotsr.io import native

    assert native.default_prefix(result) == "PDB_0000001UBQ_A_protein"
    assert native.keys_filename(native.default_prefix(result)) == (
        "PDB_0000001UBQ_A_protein.keys")


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_repeated_runs_are_identical(ubq_path):
    a = biotsr.analyze(ubq_path, method="protein", chain="A")
    b = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert a.frequencies == b.frequencies
    assert a.key_order == b.key_order
    assert [t.key for t in a.triplets] == [t.key for t in b.triplets]


def test_streaming_and_retained_runs_agree(ubq_path):
    retained = biotsr.analyze(ubq_path, method="protein", chain="A")
    streamed_keys = []

    def sink(record):
        streamed_keys.append(record.key)

    streamed = biotsr.analyze(
        ubq_path, method="protein", chain="A", keep_triplets=False, triplet_sink=sink
    )
    assert streamed.triplets is None
    assert streamed.frequencies == retained.frequencies
    assert streamed.counts.to_dict() == retained.counts.to_dict()
    assert streamed_keys == [t.key for t in retained.triplets]


def test_triplet_order_follows_candidate_enumeration(ubq_path):
    result = biotsr.analyze(ubq_path, method="protein", chain="A")
    assert len(result.triplets) == result.counts.n_triangles


# ---------------------------------------------------------------------------
# The registry and its guard rails
# ---------------------------------------------------------------------------


def test_the_registered_methods_are_exactly_the_implemented_ones():
    """The one-sided methods, plus five cross methods.

    Four atomic ones migrated from legacy implementations, plus
    ``protein_protein_calpha_cross``, which has no legacy counterpart at all: it
    is the approved interface-residue-to-Ca method (ADR-052), standard only.

    **SSE** is an annotation, not a method (ADR-064), and is absent rather than
    stubbed -- a registry entry would advertise a method that does not exist.
    """
    assert list_methods() == (
        "aa_ind",
        "aa_ind_mix",
        "aa_ind_specific",
        "aa_selection",
        "amino_acid_atom_selection",
        "amino_acid_atomic",
        "carbohydrate",
        "carbohydrate_atom_selection",
        "carbohydrate_chain",
        "carbohydrate_protein_cross",
        "drug",
        "drug_atom_selection",
        "drug_protein_cross",
        "nucleotide",
        "nucleotide_atom_selection",
        "nucleotide_chain",
        "nucleotide_protein_cross",
        "protein",
        "protein_protein_calpha_cross",
        "protein_protein_cross",
    )


def test_only_cross_methods_declare_a_composition_rule():
    """A composition rule is what makes a method two-sided.

    Before any cross method existed, *no* method declared one. The other half of
    that guarantee still holds: a one-sided method
    still declares none, so the field cannot drift into meaning something else.
    """
    from biotsr.cross.candidates import COMPOSITION_RULE

    for name in list_methods():
        recipe = get_method(name)
        if name.endswith("_cross"):
            assert recipe.composition_rule == COMPOSITION_RULE
            assert recipe.family == "cross"
        else:
            assert recipe.composition_rule is None
            assert recipe.family != "cross"


def test_the_protein_recipe_declares_no_completeness_policy():
    """TSR.py has no completeness filter, so none is invented."""
    recipe = get_method("protein")
    assert recipe.completeness_policy is None
    assert recipe.composition_rule is None


def test_the_protein_recipe_declares_the_three_calpha_modifiers_and_sse():
    """Three modifiers, plus the SSE annotation.

    ``sse`` sits in the same declaration for one reason only: so that asking for
    it on a family where the 18 categories are undefined is refused by name. It
    is not a modifier in the sense the other three are -- it changes neither
    identity (``modifiers/base.py``).
    """
    recipe = get_method("protein")
    assert recipe.supported_modifiers == frozenset(
        {"aa_grouping", "mirror", "size_filter", "sse"}
    )
    assert get_method("aa_selection").supported_modifiers == frozenset(
        {"aa_grouping", "mirror", "size_filter", "sse"}
    )


def test_unsupported_modifiers_fail_loudly_rather_than_being_ignored(ubq_path):
    """Asking for something a method does not declare names what is supported."""
    for modifier in ("annotate", "labels", "scope", "completeness"):
        with pytest.raises(ValueError, match="does not support"):
            biotsr.analyze(ubq_path, method="protein", chain="A", **{modifier: True})


def test_sse_is_refused_by_every_family_it_is_undefined_for(ubq_path):
    """ADR-036's boundary, asserted rather than documented.

    An atom inside a residue has no secondary structure of its own, and a ligand,
    a glycan and a nucleotide have none at all. Only the two Ca methods take
    ``sse``; everything else refuses it and says what it does support.
    """
    for method in sorted(biotsr.list_methods()):
        supported = get_method(method).supported_modifiers
        if method in ("protein", "aa_selection"):
            assert "sse" in supported, method
        else:
            assert "sse" not in supported, method
            with pytest.raises(ValueError, match="does not support"):
                get_method(method).validate_modifiers({"sse": True})


def test_unknown_method_names_the_registered_ones(ubq_path):
    with pytest.raises(KeyError, match="registered methods: .*protein"):
        biotsr.analyze(ubq_path, method="not_a_method", chain="A")


def test_a_cross_method_is_refused_by_analyze_naming_analyze_cross(ubq_path):
    """Until 1.0.0rc5 a cross method passed to ``analyze`` fell through to the
    Cα selector with the atom lexicon and failed on the first residue name;
    the old test above matched that message by accident."""
    with pytest.raises(biotsr.InvalidRequest, match="analyze_cross"):
        biotsr.analyze(ubq_path, method="drug_protein_cross", chain="A")


def test_unknown_parser_backend_is_refused(ubq_path):
    with pytest.raises(ValueError, match="pdb_legacy"):
        biotsr.load_structure(ubq_path, backend="pdb_strict")


# ---------------------------------------------------------------------------
# The core loop in isolation
# ---------------------------------------------------------------------------


def test_core_handles_a_point_set_too_small_for_a_triangle():
    from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
    from biotsr.domain.point import LabeledPoint, Point, PointSet

    def point(i):
        return LabeledPoint(
            identity=AtomIdentity(
                structure=StructureId("X"),
                model=1,
                chain=ChainId("A"),
                residue_name="LEU",
                residue_seq=i,
                atom_name="CA",
            ),
            coord=Point(float(i), 0.0, 0.0),
            label=8,
        )

    for n in (0, 1, 2):
        out = run_core(PointSet([point(i) for i in range(n)]), biotsr.STANDARD_PROTEIN_CALPHA_V1)
        assert out.counts.n_triangles == 0
        assert out.counts.n_unique_keys == 0
        assert out.triplets == []


def test_degenerate_candidates_are_skipped_and_recorded_not_encoded():
    """Three collinear points with the apex exactly at the base midpoint.

    The triangle is mathematically degenerate, so it must produce no key at all
    and appear in ``skipped`` with its reason -- never a theta of 0 and a normal
    key.
    """
    from biotsr.domain.identity import AtomIdentity, ChainId, StructureId
    from biotsr.domain.point import LabeledPoint, Point, PointSet

    def point(i, x):
        return LabeledPoint(
            identity=AtomIdentity(
                structure=StructureId("X"),
                model=1,
                chain=ChainId("A"),
                residue_name="LEU",
                residue_seq=i,
                atom_name="CA",
            ),
            coord=Point(x, 0.0, 0.0),
            label=8,
        )

    # Labels all equal -> case 1 -> the longest edge (-2, 2) is the base, and the
    # third point sits exactly on its midpoint.
    points = PointSet([point(0, -2.0), point(1, 2.0), point(2, 0.0)])
    out = run_core(points, biotsr.STANDARD_PROTEIN_CALPHA_V1)

    assert out.counts.n_candidates == 1
    assert out.counts.n_triangles == 0
    assert out.counts.n_degenerate == 1
    assert out.counts.n_unique_keys == 0
    assert out.key_counts == {}
    assert len(out.skipped) == 1
    assert out.skipped[0].stage == "geometry"
    assert out.skipped[0].reason == "apex_at_base_midpoint"
    assert out.skipped[0].indices is not None
