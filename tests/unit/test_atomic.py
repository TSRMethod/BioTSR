"""The atomic foundation: element labels, 18-bin scheme, admission, completeness, selection."""

from __future__ import annotations

import os

import pytest

import biotsr
from biotsr.binning.distance import (
    ATOMIC_BOUNDARIES,
    DISTANCE_SCHEME_ATOMIC_18,
    DISTANCE_SCHEME_CALPHA_35,
)
from biotsr.domain.identity import ChainId
from biotsr.labels.elements import (
    ATOMIC_LABEL_RADIX,
    ELEMENT_LABELS_112,
    element_label_scheme,
    load_element_labels,
)
from biotsr.labels.schemes import LabelDomain
from biotsr.modifiers.base import IncompatibleModifier
from biotsr.profiles.atomic import (
    STANDARD_AMINOACID_ATOM_SELECTION_V1,
    STANDARD_AMINOACID_ATOMIC_V1,
)
from biotsr.selection.atomic import (
    STANDARD_AMINOACID_RESIDUE_ADMISSION,
    WEBSITE_TSR_ADMISSION,
    AtomSelectionRequest,
    parse_atom_selection,
    residue_targets,
    select_atom_selection,
    select_residue_atoms,
)
from biotsr.selection.completeness import (
    EXACT_LEGACY_ATOM_COUNT,
    LEGACY_AMINOACID_ATOM_COUNTS,
    NO_COMPLETENESS_POLICY,
    ExactAtomCountPolicy,
)
from biotsr.structure.backends.pdb_legacy import read_pdb_legacy

from ..conftest import FIXTURE_PDB, WORKSPACE, _require

EPS = 1e-9
ACCEPT_PDB = os.path.join(WORKSPACE, "tests", "data", "structures")


@pytest.fixture(scope="module")
def ubq():
    path = os.path.join(FIXTURE_PDB, "1UBQ.pdb")
    _require(path)
    return read_pdb_legacy(path)


@pytest.fixture(scope="module")
def mapk():
    path = os.path.join(ACCEPT_PDB, "6G9M.pdb")
    _require(path)
    return read_pdb_legacy(path)


# ===========================================================================
# The 18-class atomic distance scheme
# ===========================================================================


def test_the_atomic_scheme_has_seventeen_cut_points_and_eighteen_classes():
    assert len(ATOMIC_BOUNDARIES) == 17
    assert ATOMIC_BOUNDARIES == tuple(float(n) for n in range(1, 18))
    assert DISTANCE_SCHEME_ATOMIC_18.n_classes == 18
    assert DISTANCE_SCHEME_ATOMIC_18.max_emitted_class == 18


def test_the_atomic_scheme_has_no_unreachable_class():
    """Unlike the Ca binner, this one is internally consistent."""
    assert not DISTANCE_SCHEME_ATOMIC_18.allow_unreachable_classes
    assert DISTANCE_SCHEME_ATOMIC_18.max_emitted_class == DISTANCE_SCHEME_ATOMIC_18.n_classes


@pytest.mark.parametrize("index", range(len(ATOMIC_BOUNDARIES)))
def test_every_atomic_boundary_is_lower_inclusive_upper_exclusive(index):
    cut = ATOMIC_BOUNDARIES[index]
    assert DISTANCE_SCHEME_ATOMIC_18.classify(cut - EPS) == index + 1
    assert DISTANCE_SCHEME_ATOMIC_18.classify(cut) == index + 2


@pytest.mark.parametrize(
    "value,expected",
    [
        (0.0, 1), (0.5, 1), (0.999999, 1),
        (1.0, 2), (1.5, 2),
        (8.5, 9),
        (16.999999, 17),
        (17.0, 18), (17.5, 18), (100.0, 18), (1e9, 18),
    ],
)
def test_atomic_classification(value, expected):
    assert DISTANCE_SCHEME_ATOMIC_18.classify(value) == expected


def test_the_atomic_and_calpha_schemes_are_different_key_spaces():
    """ADR-021: 18 for every atomic family, 35 for Ca, never interchanged."""
    assert (
        DISTANCE_SCHEME_ATOMIC_18.content_fingerprint()
        != DISTANCE_SCHEME_CALPHA_35.content_fingerprint()
    )
    assert DISTANCE_SCHEME_ATOMIC_18.boundaries != DISTANCE_SCHEME_CALPHA_35.boundaries


# ===========================================================================
# The element labels
# ===========================================================================


def test_the_radix_is_not_the_label_count():
    """112 is the encoder radix; the alphabet is 109 labels. Measured, not assumed."""
    assert ELEMENT_LABELS_112.n_labels == ATOMIC_LABEL_RADIX == 112
    assert ELEMENT_LABELS_112.n_distinct_codes == 109
    assert ELEMENT_LABELS_112.code_range == (4, 112)


def test_the_radix_wider_than_the_alphabet_needs_no_exception():
    """The window rule (ADR-107): 109 labels over 4..112 fit one radix-112 window,
    so the scheme is valid as it stands -- the three unused radix slots are
    visible in its description and label count, not hidden behind an override."""
    assert "109 distinct labels" in ELEMENT_LABELS_112.description
    assert "radix 112" in ELEMENT_LABELS_112.description
    assert not hasattr(ELEMENT_LABELS_112, "allow_radix_mismatch")


def test_the_encoding_is_still_injective():
    """109 consecutive labels span less than one radix, so no digit can carry."""
    assert ELEMENT_LABELS_112.codes_are_consecutive
    assert ELEMENT_LABELS_112.encoding_is_injective


def test_the_labels_are_in_the_atom_domain_not_the_amino_acid_one():
    assert ELEMENT_LABELS_112.domain == LabelDomain.ATOM_ELEMENT


def test_the_table_loads_once_with_its_digest():
    table, digest = load_element_labels()
    assert len(digest) == 64
    assert table == dict(ELEMENT_LABELS_112.table)
    assert len(table) == 110 and len(set(table.values())) == 109


@pytest.mark.parametrize(
    "atom_name,element,code",
    [
        ("N", "N", 108), ("CA", "C", 106), ("C", "C", 106), ("O", "O", 109),
        ("CB", "C", 106), ("OD1", "O", 109), ("NZ", "N", 108), ("SD", "S", 111),
        ("P", "P", 110),
    ],
)
def test_amino_acid_heavy_atoms_resolve_to_element_level_codes(atom_name, element, code):
    assert ELEMENT_LABELS_112.resolve_atom(atom_name, element) == code


def test_cysteine_sg_is_labelled_sulfur_and_seaborgium_keeps_72():
    """The label follows the element field: cysteine's SG (element S) is sulfur,
    a record whose element is Sg is Seaborgium."""
    assert ELEMENT_LABELS_112.resolve_atom("SG", "S") == 111
    assert ELEMENT_LABELS_112.resolve_atom("SG", "SG") == 72
    for sibling in ("SG1", "SG2", "SG5", "SD"):
        assert ELEMENT_LABELS_112.resolve_atom(sibling, "S") == 111


def test_an_unknown_element_returns_none_rather_than_guessing():
    assert ELEMENT_LABELS_112.resolve_atom("C1", "X") is None
    assert ELEMENT_LABELS_112.resolve_atom("C1", "") is None
    assert ELEMENT_LABELS_112.label("NOT_AN_ATOM") is None


def test_a_second_scheme_instance_fingerprints_identically():
    assert (
        element_label_scheme().content_fingerprint()
        == ELEMENT_LABELS_112.content_fingerprint()
    )


# ===========================================================================
# Completeness
# ===========================================================================


def test_the_completeness_table_has_the_twenty_three_legacy_entries():
    assert len(LEGACY_AMINOACID_ATOM_COUNTS) == 23
    for phosphorylated in ("TPO", "SEP", "PTR"):
        assert phosphorylated in LEGACY_AMINOACID_ATOM_COUNTS


@pytest.mark.parametrize(
    "residue,count", sorted(LEGACY_AMINOACID_ATOM_COUNTS.items())
)
def test_a_residue_with_exactly_the_expected_count_is_admitted(residue, count):
    assert EXACT_LEGACY_ATOM_COUNT.admit(residue, count).admitted


def test_the_rule_is_exact_equality_not_a_minimum():
    """One atom too many is rejected exactly as firmly as one too few.

    This is what excludes 1UBQ's C-terminal Gly76, which has an OXT and so five
    atoms against the table's four.  Reading the rule as "at least N" would
    silently admit 76 residues where history admits 75.
    """
    assert EXACT_LEGACY_ATOM_COUNT.admit("GLY", 4).admitted
    assert not EXACT_LEGACY_ATOM_COUNT.admit("GLY", 5).admitted
    assert not EXACT_LEGACY_ATOM_COUNT.admit("GLY", 3).admitted


def test_rejection_carries_a_reason_and_a_detail():
    verdict = EXACT_LEGACY_ATOM_COUNT.admit("GLY", 5)
    assert verdict.reason == "atom_count_mismatch"
    assert "expected exactly 4" in verdict.detail and "got 5" in verdict.detail


def test_an_unlisted_residue_type_is_rejected_whatever_its_count():
    for count in (0, 4, 8, 20):
        verdict = EXACT_LEGACY_ATOM_COUNT.admit("MSE", count)
        assert not verdict.admitted
        assert verdict.reason == "unlisted_entity_type"


def test_the_no_completeness_policy_admits_everything():
    assert NO_COMPLETENESS_POLICY.admit("GLY", 0).admitted
    assert NO_COMPLETENESS_POLICY.admit("ANYTHING", 999).admitted


def test_a_policy_is_a_value_object_other_families_can_supply():
    custom = ExactAtomCountPolicy(expected={"XYZ": 3}, name="custom")
    assert custom.admit("XYZ", 3).admitted
    assert not custom.admit("GLY", 4).admitted


# ===========================================================================
# Atom admission
# ===========================================================================


def test_the_two_lineages_share_the_altloc_rule_and_the_element_column(ubq):
    """One alternate-location rule for every method (blank or A, ADR-095), and
    since 1.0.0rc6 one element column too, 77-78 (ADR-108)."""
    for admission in (STANDARD_AMINOACID_RESIDUE_ADMISSION, WEBSITE_TSR_ADMISSION):
        assert admission.accepted_alt_locs == ("", "A")
        assert not hasattr(admission, "use_legacy_tokens")
    assert not hasattr(STANDARD_AMINOACID_RESIDUE_ADMISSION, "reject_alt_loc_prefix")


def test_hetatm_is_never_admitted_by_either(ubq):
    hetero = [a for a in ubq.atoms if a.hetero]
    assert hetero, "1UBQ has HETATM water records"
    for policy in (STANDARD_AMINOACID_RESIDUE_ADMISSION, WEBSITE_TSR_ADMISSION):
        ok, reason = policy.accepts(hetero[0])
        assert not ok and reason == "hetatm_excluded"


def test_the_off_by_one_element_token_is_no_longer_captured(ubq):
    """Columns 78-80 read into the charge field; nothing reads them now
    (ADR-108), so the token is gone rather than kept unused."""
    ca = [a for a in ubq.atoms if a.atom_name == "CA" and not a.hetero][0]
    assert ca.element == "C"
    assert not hasattr(ca, "legacy_element_token")


def test_the_historical_altloc_token_is_no_longer_captured(ubq):
    """`line[16:20]` (altLoc plus residue name) fed the amino-acid analyzer's
    "not B" rule; the standard reads the parsed altLoc field instead, so the
    token is gone rather than kept unused."""
    ca = [a for a in ubq.atoms if a.atom_name == "CA" and not a.hetero][0]
    assert ca.identity.alt_loc == ""
    assert not hasattr(ca, "legacy_altloc_residue_token")


# ===========================================================================
# Per-residue selection
# ===========================================================================


def test_residue_targets_lists_every_residue_of_the_chain(ubq):
    targets = residue_targets(ubq, "A")
    assert len(targets) == 76
    assert targets[0] == ("1", "MET")
    assert targets[-1] == ("76", "GLY")


def test_a_complete_residue_is_selected_with_its_expected_atoms(ubq):
    result = select_residue_atoms(
        ubq, "A", "15", "LEU", ELEMENT_LABELS_112,
        completeness=EXACT_LEGACY_ATOM_COUNT,
    )
    assert len(result.points) == LEGACY_AMINOACID_ATOM_COUNTS["LEU"] == 8
    assert [p.identity.atom_name for p in result.points] == [
        "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2",
    ]


def test_every_selected_point_carries_full_physical_identity(ubq):
    result = select_residue_atoms(
        ubq, "A", "15", "LEU", ELEMENT_LABELS_112,
        completeness=EXACT_LEGACY_ATOM_COUNT,
    )
    for p in result.points:
        i = p.identity
        assert i.structure.value == "1UBQ"
        assert i.chain == ChainId("A")
        assert i.residue_name == "LEU"
        assert i.residue_seq == 15
        assert i.model == 1
        assert i.atom_name
    # Full identity, so two atoms are never conflated by name alone.
    assert len({p.identity.as_tuple() for p in result.points}) == len(result.points)


def test_the_terminal_residue_is_rejected_and_says_why(ubq):
    """1UBQ's Gly76 carries an OXT: five atoms against the table's four."""
    result = select_residue_atoms(
        ubq, "A", "76", "GLY", ELEMENT_LABELS_112,
        completeness=EXACT_LEGACY_ATOM_COUNT,
    )
    assert len(result.points) == 0
    reasons = [s for s in result.skipped if s.stage == "completeness"]
    assert len(reasons) == 1
    assert reasons[0].reason == "atom_count_mismatch"
    assert "expected exactly 4 atoms, got 5" in reasons[0].detail


def test_without_a_completeness_policy_the_same_residue_is_admitted(ubq):
    """Proving the exclusion is the policy's doing, not the selector's."""
    result = select_residue_atoms(
        ubq, "A", "76", "GLY", ELEMENT_LABELS_112,
        completeness=NO_COMPLETENESS_POLICY,
    )
    assert len(result.points) == 5
    assert "OXT" in [p.identity.atom_name for p in result.points]


def test_hydrogens_are_excluded_when_the_structure_has_them(mapk):
    """6G9M is an X-ray structure with no hydrogens; assert the rule holds anyway
    by checking no admitted atom is one."""
    result = select_residue_atoms(
        mapk, "A", "40", "CYS", ELEMENT_LABELS_112,
        completeness=EXACT_LEGACY_ATOM_COUNT,
    )
    assert len(result.points) == LEGACY_AMINOACID_ATOM_COUNTS["CYS"] == 6
    assert all(p.identity.atom_name[0] not in ("H", "D") for p in result.points)


def test_the_residue_name_must_match_not_just_the_number(ubq):
    """A sequence number alone does not identify a residue."""
    wrong = select_residue_atoms(
        ubq, "A", "15", "ALA", ELEMENT_LABELS_112,
        completeness=NO_COMPLETENESS_POLICY,
    )
    assert len(wrong.points) == 0


# ===========================================================================
# Atom selection
# ===========================================================================


def test_parse_atom_selection_reads_the_acceptance_sheet_format():
    requests, warnings = parse_atom_selection("LYS_48_NZ;HIS_68_NE2")
    assert warnings == ()
    assert requests == (
        AtomSelectionRequest("LYS", "48", "NZ"),
        AtomSelectionRequest("HIS", "68", "NE2"),
    )


def test_parse_atom_selection_splits_from_the_left_exactly_twice():
    """So an atom name containing an underscore survives."""
    requests, _ = parse_atom_selection("LYS_48_N_Z")
    assert requests == (AtomSelectionRequest("LYS", "48", "N_Z"),)


def test_parse_atom_selection_drops_duplicates_keeping_order():
    requests, _ = parse_atom_selection("A_1_X;B_2_Y;A_1_X")
    assert [r.token() for r in requests] == ["A_1_X", "B_2_Y"]


def test_parse_atom_selection_refuses_malformed_tokens_rather_than_dropping_them():
    """A token that is not RESIDUE_ID_ATOM is refused by name: dropping it would
    change the triangle set, and guessing could change it differently."""
    from biotsr.domain.errors import InvalidRequest

    for bad in ("nonsense", "LYS_X_NZ", "LYS_48", "LYS_4 8_NZ", "LYSINE_48_NZ"):
        with pytest.raises(InvalidRequest, match=bad.split("_")[0]):
            parse_atom_selection("LYS_48_NZ;%s;HIS_68_NE2" % bad)
    requests, warnings = parse_atom_selection("LYS_48_NZ;HIS_68_NE2;LYS_52A_NZ;GLY_-3_CA")
    assert [r.residue_seq for r in requests] == ["48", "68", "52A", "-3"]
    assert warnings == ()


def test_parse_atom_selection_accepts_a_list_and_empty_input():
    requests, _ = parse_atom_selection(["LYS_48_NZ", "HIS_68_NE2"])
    assert len(requests) == 2
    assert parse_atom_selection(None) == ((), ())
    assert parse_atom_selection("") == ((), ())


def test_atom_selection_pools_across_residues(mapk):
    requests, _ = parse_atom_selection("TYR_30_OH;TYR_30_CZ;LYS_207_NZ")
    result = select_atom_selection(mapk, "A", requests, ELEMENT_LABELS_112)
    assert len(result.points) == 3
    residues = {p.identity.residue_seq for p in result.points}
    assert residues == {30, 207}, "one pooled calculation spanning two residues"


def test_atom_selection_returns_structure_order_not_request_order(mapk):
    requests, _ = parse_atom_selection("LYS_207_NZ;TYR_30_OH")
    result = select_atom_selection(mapk, "A", requests, ELEMENT_LABELS_112)
    assert [p.identity.residue_seq for p in result.points] == [30, 207]


def test_a_requested_atom_that_is_absent_is_reported(mapk):
    requests, _ = parse_atom_selection("LYS_207_NZ;LYS_207_NOPE")
    result = select_atom_selection(mapk, "A", requests, ELEMENT_LABELS_112)
    assert len(result.points) == 1
    assert result.atoms_missing == ("LYS_207_NOPE",)
    assert any(s.reason == "requested_atom_not_found" for s in result.skipped)


def test_a_residue_name_mismatch_is_reported_not_silently_accepted(mapk):
    """The CSV and the structure disagreeing is worth saying out loud."""
    requests, _ = parse_atom_selection("ALA_207_NZ")
    result = select_atom_selection(mapk, "A", requests, ELEMENT_LABELS_112)
    assert len(result.points) == 0
    assert any(s.reason == "residue_name_mismatch" for s in result.skipped)


# ===========================================================================
# Profiles and identity
# ===========================================================================


def test_each_standard_atomic_profile_pins_its_selection_rules():
    """Which atoms a profile admits is part of the standard, not an accident.

    The residue-complete method uses the amino-acid admission rule with the exact
    per-residue atom-count table; the atom-selection method uses the website
    admission rule with no completeness filter.
    """
    assert STANDARD_AMINOACID_ATOMIC_V1.atom_admission.name == STANDARD_AMINOACID_RESIDUE_ADMISSION.name
    assert STANDARD_AMINOACID_ATOMIC_V1.completeness.name == EXACT_LEGACY_ATOM_COUNT.name
    assert STANDARD_AMINOACID_ATOM_SELECTION_V1.atom_admission.name == WEBSITE_TSR_ADMISSION.name
    assert STANDARD_AMINOACID_ATOM_SELECTION_V1.completeness.name == NO_COMPLETENESS_POLICY.name
    for profile in (STANDARD_AMINOACID_ATOMIC_V1, STANDARD_AMINOACID_ATOM_SELECTION_V1):
        assert profile.label_scheme is ELEMENT_LABELS_112
        assert profile.angle_convention.name == "pi"
        assert not profile.is_legacy_compatibility


def test_each_method_defaults_to_the_matching_standard_profile():
    assert (
        biotsr.methods.get_method("amino_acid_atomic").default_profile
        is STANDARD_AMINOACID_ATOMIC_V1
    )
    assert (
        biotsr.methods.get_method("amino_acid_atom_selection").default_profile
        is STANDARD_AMINOACID_ATOM_SELECTION_V1
    )


def test_the_atomic_profiles_use_the_18_class_scheme_never_the_calpha_one():
    for profile in (
        STANDARD_AMINOACID_ATOMIC_V1,
        STANDARD_AMINOACID_ATOM_SELECTION_V1,
    ):
        assert profile.distance_scheme is DISTANCE_SCHEME_ATOMIC_18
        assert profile.key_formula().n_distance == 18
        assert profile.key_formula().n_theta == 29
        assert profile.key_formula().n_labels == 112


def test_the_atomic_key_space_differs_from_protein_calpha():
    atomic = STANDARD_AMINOACID_ATOMIC_V1.key_space()
    calpha = biotsr.STANDARD_PROTEIN_CALPHA_V1.key_space()
    assert atomic.fingerprint != calpha.fingerprint
    assert atomic.distance_scheme.radix != calpha.distance_scheme.radix
    assert atomic.label_scheme.radix != calpha.label_scheme.radix


def test_changing_the_atomic_label_scheme_changes_the_fingerprint():
    from biotsr.labels.schemes import LabelScheme

    altered = LabelScheme(
        name="element_labels_112",
        table={"N": 1, "C": 2},
        n_labels=112,
        domain=LabelDomain.ATOM_ELEMENT,
    )
    assert altered.content_fingerprint() != ELEMENT_LABELS_112.content_fingerprint()


def test_changing_the_distance_semantics_changes_the_fingerprint():
    import dataclasses

    swapped = dataclasses.replace(
        STANDARD_AMINOACID_ATOMIC_V1, distance_scheme=DISTANCE_SCHEME_CALPHA_35
    )
    assert swapped.key_space().fingerprint != STANDARD_AMINOACID_ATOMIC_V1.key_space().fingerprint


def test_aa_grouping_is_rejected_for_the_atomic_family():
    """Ca-only, enforced by label domain rather than by a per-family list."""
    with pytest.raises(IncompatibleModifier, match="Protein Ca label domain"):
        STANDARD_AMINOACID_ATOMIC_V1.with_modifiers(aa_grouping=True)


def test_the_atomic_methods_do_not_even_declare_grouping():
    """Refused twice over: the registry rejects the name, and the domain would
    reject it anyway."""
    for name in ("amino_acid_atomic", "amino_acid_atom_selection"):
        recipe = biotsr.methods.get_method(name)
        assert "aa_grouping" not in recipe.supported_modifiers
        assert recipe.supported_modifiers == frozenset({"mirror", "size_filter"})


# ===========================================================================
# The sqrt convention -- discovered by characterizing the two lineages
# ===========================================================================


def test_the_two_lineages_take_the_square_root_differently():
    from biotsr.geometry.conventions import SQRT_MATH_SQRT, SQRT_POW_HALF

    assert biotsr.STANDARD_PROTEIN_CALPHA_V1.sqrt_convention is SQRT_POW_HALF
    assert STANDARD_AMINOACID_ATOMIC_V1.sqrt_convention is SQRT_MATH_SQRT


def test_the_two_spellings_really_do_disagree():
    """Not a theoretical distinction: they differ on real values.

    ``** 0.5`` goes through ``pow``, which IEEE-754 does not require to be
    correctly rounded; ``math.sqrt`` is the correctly-rounded hardware root.
    """
    import random

    from biotsr.geometry.conventions import SQRT_MATH_SQRT, SQRT_POW_HALF

    random.seed(1)
    disagreements = sum(
        1
        for _ in range(20000)
        if SQRT_POW_HALF(v := random.uniform(0.01, 2000.0)) != SQRT_MATH_SQRT(v)
    )
    assert disagreements > 0, "the two spellings must actually differ somewhere"


def test_the_sqrt_convention_enters_the_key_space():
    """One ULP in a distance can move a triangle across a bin boundary, so this
    is key-affecting and cannot be silently mixed."""
    import dataclasses

    from biotsr.geometry.conventions import SQRT_POW_HALF

    swapped = dataclasses.replace(
        STANDARD_AMINOACID_ATOMIC_V1, sqrt_convention=SQRT_POW_HALF
    )
    assert (
        swapped.key_space().fingerprint
        != STANDARD_AMINOACID_ATOMIC_V1.key_space().fingerprint
    )
    assert STANDARD_AMINOACID_ATOMIC_V1.key_space().sqrt_convention == "math_sqrt"
