"""The label-space rule, AA selection semantics, and size filtering.

**The settled label-space rule.**

======================================  ======
family                                  radix
======================================  ======
every C-alpha method                    **30**
C-alpha + AA grouping                   **30**
every non-C-alpha atomic method         **112**
======================================  ======

The rule is about the *family*, not the candidate or residue-selection strategy:
an exhaustive whole-chain run, a residue-selected run, a size-filtered run and a
mirrored run all stay in the C-alpha residue label space. The ungrouped and the
grouped residue labels share radix 30 and are still two key spaces (ADR-107).

**The five amino-acid selection modes of the acceptance sheet.**

**Size filtering is optional and off by default.**
"""

from __future__ import annotations

import pytest

import biotsr
from biotsr.labels.schemes import LabelDomain
from biotsr.methods.registry import get_method, list_methods
from biotsr.profiles.atomic import ATOMIC_PROFILES
from biotsr.selection.atomic import (
    ResidueRef,
    parse_atom_names,
    parse_residue_refs,
    parse_residue_types,
)
from biotsr.selection.protein import parse_residue_selection

# The methods whose points are protein residues, whatever strategy chose them.
CALPHA_METHODS = ("protein", "aa_selection")
# Ca protein-protein cross is a residue method too -- an
# interface residue is represented by its Ca -- so it labels at radix 30 on the
# 35-class Ca distance scale, exactly like the other two. It is listed separately
# only because `aa_grouping` is not offered on it: the method is scoped to the
# ungrouped residue labels and grouping is left to a future mode (ADR-052).
CALPHA_CROSS_METHODS = ("protein_protein_calpha_cross",)
ALL_CALPHA_METHODS = CALPHA_METHODS + CALPHA_CROSS_METHODS
# Everything else that produces atomic keys.
ATOMIC_METHODS = (
    "amino_acid_atomic",
    "amino_acid_atom_selection",
    "aa_ind",
    "aa_ind_mix",
    "aa_ind_specific",
    "drug",
    "drug_atom_selection",
    "carbohydrate",
    "carbohydrate_chain",
    "carbohydrate_atom_selection",
    "nucleotide",
    "nucleotide_chain",
    "nucleotide_atom_selection",
)
# Every cross method is FULLY ATOMIC too -- both sides are generic
# atoms by element at radix 112, and a Ca residue side is never mixed in
# (ADR-040). They are listed separately because they differ from the one-sided
# atomic methods in exactly one declared respect: they do NOT offer `mirror`.
CROSS_METHODS = (
    "drug_protein_cross",
    "carbohydrate_protein_cross",
    "nucleotide_protein_cross",
    "protein_protein_cross",
)
ALL_ATOMIC_METHODS = ATOMIC_METHODS + CROSS_METHODS


# ===========================================================================
# The label-space rule
# ===========================================================================


def test_every_registered_method_is_in_exactly_one_of_the_two_families():
    assert set(ALL_CALPHA_METHODS) | set(ALL_ATOMIC_METHODS) == set(list_methods())
    assert not set(ALL_CALPHA_METHODS) & set(ALL_ATOMIC_METHODS)


@pytest.mark.parametrize("name", ALL_CALPHA_METHODS)
def test_calpha_methods_use_radix_30(name):
    profile = get_method(name).default_profile
    assert profile.label_scheme.name == "amino_acid_30"
    assert profile.label_scheme.n_labels == 30
    assert profile.label_scheme.domain == LabelDomain.PROTEIN_AMINO_ACID
    assert profile.distance_scheme.n_classes == 35


@pytest.mark.parametrize("name", CALPHA_METHODS)
def test_grouped_calpha_methods_use_radix_30(name):
    profile = get_method(name).default_profile.with_modifiers(aa_grouping=True)
    assert profile.label_scheme.name == "amino_acid_group_23"
    assert profile.label_scheme.n_labels == 30
    assert profile.label_scheme.domain == LabelDomain.PROTEIN_AMINO_ACID


@pytest.mark.parametrize("name", ALL_ATOMIC_METHODS)
def test_every_non_calpha_atomic_method_uses_radix_112(name):
    profile = get_method(name).default_profile
    assert profile.label_scheme.n_labels == 112
    assert profile.label_scheme.domain == LabelDomain.ATOM_ELEMENT


@pytest.mark.parametrize("name", ALL_ATOMIC_METHODS)
def test_every_non_calpha_atomic_method_uses_18_distance_bins(name):
    assert get_method(name).default_profile.distance_scheme.n_classes == 18


@pytest.mark.parametrize("name", ALL_ATOMIC_METHODS)
def test_no_atomic_method_can_reach_the_calpha_label_schemes(name):
    profile = get_method(name).default_profile
    assert profile.label_scheme.name == "element_labels_112"
    assert "amino_acid" not in profile.label_scheme.name
    assert profile.label_scheme.domain != LabelDomain.PROTEIN_AMINO_ACID


@pytest.mark.parametrize("name", ALL_ATOMIC_METHODS)
def test_aa_grouping_is_refused_for_every_atomic_method(name):
    """Refused twice over: the registry rejects the name, and the domain rejects
    the scheme even if the declaration were removed (ADR-023)."""
    from biotsr.modifiers.aa_grouping import IncompatibleModifier

    recipe = get_method(name)
    assert "aa_grouping" not in recipe.supported_modifiers
    with pytest.raises(ValueError, match="aa_grouping"):
        recipe.validate_modifiers({"aa_grouping": True})
    with pytest.raises(IncompatibleModifier):
        recipe.default_profile.with_modifiers(aa_grouping=True, method=name)


@pytest.mark.parametrize("profile", ATOMIC_PROFILES, ids=lambda p: p.name)
def test_every_atomic_profile_is_112_over_18(profile):
    assert profile.label_scheme.n_labels == 112
    assert profile.distance_scheme.n_classes == 18
    assert profile.theta_scheme.n_classes == 29


def test_the_three_label_spaces_have_three_different_keyspace_fingerprints():
    calpha = get_method("protein").default_profile
    grouped = calpha.with_modifiers(aa_grouping=True)
    atomic = get_method("drug").default_profile
    prints = {
        calpha.key_space().fingerprint,
        grouped.key_space().fingerprint,
        atomic.key_space().fingerprint,
    }
    assert len(prints) == 3
    assert calpha.key_space().label_scheme.radix == 30
    assert grouped.key_space().label_scheme.radix == 30
    assert atomic.key_space().label_scheme.radix == 112


def test_every_atomic_family_shares_one_label_scheme_fingerprint():
    """Drug, carbohydrate, nucleotide and AA atomic keys are mutually comparable."""
    prints = {
        get_method(name).default_profile.key_space().label_scheme.content_fingerprint
        for name in ATOMIC_METHODS
    }
    assert len(prints) == 1


def test_calpha_label_space_survives_every_modifier_and_strategy():
    """Candidate/selection strategy never changes the label space (A4)."""
    base = get_method("aa_selection").default_profile
    for kwargs in (
        {},
        {"mirror": True},
        {"size_filter": 25.0},
        {"mirror": True, "size_filter": 25.0},
    ):
        assert base.with_modifiers(**kwargs).label_scheme.n_labels == 30


# ===========================================================================
# Amino-acid selection semantics
# ===========================================================================


def test_residue_type_parsing():
    assert parse_residue_types("TYR")[0] == ("TYR",)
    assert parse_residue_types("TYR;TPO")[0] == ("TYR", "TPO")
    assert parse_residue_types("tyr ; tpo ;")[0] == ("TYR", "TPO")
    assert parse_residue_types(["LEU", "ILE", "LEU"])[0] == ("LEU", "ILE")
    assert parse_residue_types("")[0] == ()


def test_residue_ref_parsing_accepts_both_spellings_dataset_sheets_use():
    refs, warnings = parse_residue_refs("LYS_6;GLU18")
    assert refs == (ResidueRef("LYS", "6"), ResidueRef("GLU", "18"))
    assert warnings == ()


def test_residue_ref_parsing_keeps_an_insertion_code_with_the_sequence():
    refs, _ = parse_residue_refs("SER52A;SER_52A")
    assert refs == (ResidueRef("SER", "52A"),)


def test_a_malformed_residue_ref_is_refused_not_silently_dropped():
    from biotsr.domain.errors import InvalidRequest

    with pytest.raises(InvalidRequest, match="NOTARESIDUE"):
        parse_residue_refs("LYS_6;NOTARESIDUE;TYR12")


def test_both_residue_spellings_of_the_dataset_sheets_are_read_explicitly():
    """``LYS_6`` and ``GLU18`` share one cell in the TSR sheets; both are
    supported, normalised to one identity, and one residue named twice is one."""
    refs, warnings = parse_residue_refs("LYS_6;LYS11;LYS_11;GLU18;LYS11A;5CM_3")
    assert refs == (ResidueRef("LYS", "6"), ResidueRef("LYS", "11"),
                    ResidueRef("GLU", "18"), ResidueRef("LYS", "11A"),
                    ResidueRef("5CM", "3"))
    assert warnings == ()


def test_a_compact_spelling_that_cannot_be_split_is_refused():
    """A name containing a digit gives no boundary in NAMEID; never guessed."""
    from biotsr.domain.errors import InvalidRequest

    for bad in ("5CM3", "LYS", "LYS_A", "LYS_1 4", "LYS_11AB"):
        with pytest.raises(InvalidRequest):
            parse_residue_refs(bad)


def test_residue_selection_expands_ranges_like_dataset_sheets_mean_them():
    labels, warnings = parse_residue_selection("10-12;48;70-71")
    assert labels == ("10", "11", "12", "48", "70", "71")
    assert warnings == ()


def test_an_insertion_code_range_is_refused_naming_it():
    """Until 1.0.0rc4 such a token was kept verbatim and then matched nothing;
    a range is defined only between plain numbers, so it is refused."""
    from biotsr.domain.errors import InvalidRequest

    with pytest.raises(InvalidRequest, match="'27A-27B' is not supported"):
        parse_residue_selection("27A-27B;48")


def test_atom_name_parsing_preserves_case_exactly():
    """Case folding here would undo the namespace separation (ADR-031)."""
    assert parse_atom_names("N5;C6;O7") == ("N5", "C6", "O7")
    assert parse_atom_names("C2';C4;C5;O2") == ("C2'", "C4", "C5", "O2")
    assert parse_atom_names("SG") == ("SG",)
    assert parse_atom_names(None) == ()


def test_aa_ind_takes_exactly_one_type(mapk_path):
    with pytest.raises(ValueError, match="at most 1"):
        biotsr.analyze(mapk_path, method="aa_ind", chain="A", residue_types="TYR;PHE")


def test_aa_ind_mix_takes_two_or_more_types(mapk_path):
    with pytest.raises(ValueError, match="at least 2"):
        biotsr.analyze(mapk_path, method="aa_ind_mix", chain="A", residue_types="TYR")


def test_aa_ind_uses_every_eligible_non_hydrogen_atom(mapk_path):
    """One TRP is 14 atoms; the selection must be a whole number of them.

    TRP is chosen because 6G9M chain A has only three of them -- the semantics are
    identical for a commoner residue, and C(42,3) keeps the suite fast.
    """
    result = biotsr.analyze(
        mapk_path, method="aa_ind", chain="A", residue_types="TRP",
        keep_triplets=True,
    )
    assert result.counts.n_points == 42
    assert result.counts.n_points == 14 * len(result.selection.residues_used)
    assert all(r.endswith("_TRP") for r in result.selection.residues_used)
    names = {i.atom_name for t in result.triplets for i in t.identities}
    assert {"NE1", "CZ2", "CZ3", "CH2"} <= names
    assert not any(n.startswith("H") for n in names)


def test_aa_ind_mix_pools_both_types_into_one_calculation(mapk_path):
    trp = biotsr.analyze(mapk_path, method="aa_ind", chain="A",
                      residue_types="TRP", keep_triplets=False)
    met = biotsr.analyze(mapk_path, method="aa_ind", chain="A",
                      residue_types="MET", keep_triplets=False)
    mix = biotsr.analyze(mapk_path, method="aa_ind_mix", chain="A",
                      residue_types="TRP;MET", keep_triplets=False)
    assert mix.counts.n_points == trp.counts.n_points + met.counts.n_points
    # Pooled, not concatenated: the mix has triangles neither single run can have.
    assert mix.counts.n_triangles > trp.counts.n_triangles + met.counts.n_triangles
    assert set(r.split("_")[1] for r in mix.selection.residues_used) == {"TRP", "MET"}


def test_aa_ind_specific_selects_exactly_the_named_residues(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="aa_ind_specific", chain="A", residues="LYS_207;SER266"
    )
    assert set(result.selection.residues_used) == {"207_LYS", "266_SER"}
    assert result.counts.n_points == 15          # LYS (9) and SER (6) heavy atoms
    assert result.selection.atoms_missing == ()


def test_aa_ind_specific_reports_a_residue_it_could_not_find(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="aa_ind_specific", chain="A", residues="LYS_207;TRP_9999"
    )
    assert result.selection.atoms_missing == ("TRP_9999",)
    assert any(s.reason == "requested_entity_not_found" for s in result.skipped)


def test_aa_selection_is_calpha_only(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="aa_selection", chain="A",
        residues="26-28;192;207",
    )
    assert result.counts.n_points == 5           # one Ca per residue
    names = {i.atom_name for t in result.triplets for i in t.identities}
    assert names == {"CA"}
    # ...and labelled by RESIDUE identity, not by the atom name.
    assert result.key_space.label_scheme.radix == 30
    labels = {i.residue_name: l for t in result.triplets
              for i, l in zip(t.identities, t.labels)}
    assert len(set(labels.values())) > 1        # not one repeated CA code


def test_aa_atom_selection_pools_named_atoms(mapk_path):
    result = biotsr.analyze(
        mapk_path, method="amino_acid_atom_selection", chain="A",
        atoms="TYR_30_OH;TYR_30_CZ;ALA_35_CB;LYS_207_NZ;TRP_192_NE1",
    )
    assert result.counts.n_points == 5
    assert result.counts.n_candidates == 10      # C(5, 3)


def test_the_2_columns_need_no_new_method(mapk_path):
    """`aa_ind_2` etc. are the same methods against `aa_chain_2` (A5)."""
    a = biotsr.analyze(mapk_path, method="aa_ind", chain="A",
                    residue_types="TRP", keep_triplets=False)
    # A second chain is simply a second call. No cross method is involved.
    assert a.sampling.composition_rule is None
    assert a.chain == "A"


# ===========================================================================
# Size filtering is optional and OFF BY DEFAULT
# ===========================================================================


@pytest.mark.parametrize("name", list(CALPHA_METHODS) + list(ALL_ATOMIC_METHODS))
def test_size_filtering_is_off_by_default_for_every_method(name):
    profile = get_method(name).default_profile
    assert profile.size_filter is None
    assert "size_filter" not in profile.active_modifiers()


@pytest.mark.parametrize("name", list(CALPHA_METHODS) + list(ALL_ATOMIC_METHODS))
def test_size_filtering_is_available_to_every_method(name):
    recipe = get_method(name)
    assert "size_filter" in recipe.supported_modifiers
    filtered = recipe.default_profile.with_modifiers(size_filter=8.0)
    assert filtered.size_filter is not None
    assert filtered.size_filter.cutoff == 8.0


@pytest.mark.parametrize("name", list(CALPHA_METHODS) + list(ALL_ATOMIC_METHODS))
def test_size_filtering_never_changes_the_key_space(name):
    """It changes which triangles were drawn, not what a key means."""
    base = get_method(name).default_profile
    assert (
        base.with_modifiers(size_filter=8.0).key_space().fingerprint
        == base.key_space().fingerprint
    )
    assert (
        base.with_modifiers(size_filter=8.0).name != base.name
    )


def test_the_size_filter_boundary_is_strictly_less_than():
    from biotsr.modifiers.size_filter import SizeFilter

    f = SizeFilter.from_value(10.0)
    assert f.accepts(9.999999)
    assert not f.accepts(10.0)
    assert not f.accepts(10.000001)


@pytest.mark.parametrize("name", ATOMIC_METHODS)
def test_mirror_is_available_to_every_one_sided_atomic_method(name):
    recipe = get_method(name)
    assert "mirror" in recipe.supported_modifiers
    assert recipe.default_profile.mirror_convention.enabled is False


@pytest.mark.parametrize("name", CROSS_METHODS)
def test_mirror_is_not_offered_to_any_cross_method(name):
    """Mirror is NOT extended to cross by analogy.

    No legacy cross implementation contains a mirror flag, and the modifier being
    family-agnostic is a reason it *could* be offered, not evidence that it
    *should* be. Offering it would ship a capability with no reference output and
    no scientific confirmation, so it is refused and recorded as unconfirmed
    future scope (ADR-047). The refusal names what the method does support.
    """
    recipe = get_method(name)
    assert "mirror" not in recipe.supported_modifiers
    assert recipe.default_profile.mirror_convention.enabled is False
    with pytest.raises(ValueError) as excinfo:
        recipe.validate_modifiers({"mirror": True})
    assert "mirror" in str(excinfo.value)
    assert "size_filter" in str(excinfo.value)
