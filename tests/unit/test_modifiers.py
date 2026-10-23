"""The three Cα modifiers, in isolation.

Composition is tested separately (``tests/integration/test_modifier_composition.py``);
this module tests each mechanism on its own terms.
"""

from __future__ import annotations

import pytest

import biotsr
from biotsr.domain.result import DEFAULT_CATEGORY
from biotsr.labels.schemes import (
    AMINO_ACID_30,
    AMINO_ACID_GROUP_23,
    LabelDomain,
    LabelScheme,
)
from biotsr.modifiers import aa_grouping as grouping_mod
from biotsr.modifiers.base import IncompatibleModifier, ModifierHook
from biotsr.modifiers.mirror import (
    MIRROR_CATEGORY,
    MIRROR_DISABLED,
    MIRROR_NEGATE_WHEN_PREFOLD_GT_90,
)
from biotsr.modifiers.size_filter import Comparison, SizeFilter

STANDARD = biotsr.STANDARD_PROTEIN_CALPHA_V1


# ===========================================================================
# AA grouping
# ===========================================================================


def test_default_mode_uses_the_residue_label_scheme_at_radix_30():
    assert STANDARD.label_scheme is AMINO_ACID_30
    assert STANDARD.label_scheme.n_labels == 30
    assert STANDARD.label_scheme.n_distinct_codes == 30
    assert len(STANDARD.label_scheme.table) == 31, "MSE shares MET's code"
    assert STANDARD.aa_grouping is False


def test_grouping_maps_31_residues_onto_23_labels():
    assert len(AMINO_ACID_GROUP_23.table) == 31, "every residue of the base scheme"
    assert AMINO_ACID_GROUP_23.n_distinct_codes == 23
    regular = [r for r, c in AMINO_ACID_30.table.items() if c <= 23 and r != "MSE"]
    assert len({AMINO_ACID_GROUP_23.label(r) for r in regular}) == 13


def test_the_standard_grouped_radix_is_30_the_span_of_its_codes():
    """ADR-107: the radix follows the codes' span (4..33), not their count (23)."""
    grouped = STANDARD.with_modifiers(aa_grouping=True)
    assert grouped.label_scheme is AMINO_ACID_GROUP_23
    assert grouped.label_scheme.n_labels == 30
    assert grouped.label_scheme.code_range == (4, 33)
    assert grouped.key_formula().n_labels == 30
    assert grouped.key_space().label_scheme.radix == 30


def test_the_defective_grouped_scheme_is_not_reachable_from_the_standard_profile():
    """The historical 13-group/radix-20 scheme, and the 13-group/radix-13 one it
    was replaced by, are gone; only the 23-group radix-30 scheme exists."""
    import biotsr.labels.schemes as schemes

    grouped = STANDARD.with_modifiers(aa_grouping=True)
    assert grouped.label_scheme is AMINO_ACID_GROUP_23
    assert not grouped.is_legacy_compatibility
    for gone in ("AMINO_ACID_GROUP_13_LEGACY_RADIX_20", "AMINO_ACID_GROUP_13",
                 "AMINO_ACID_20"):
        assert not hasattr(schemes, gone)


def test_a_radix_narrower_than_the_codes_is_refused():
    with pytest.raises(ValueError, match="radix"):
        LabelScheme(name="bad", table={"A": 1, "B": 22}, n_labels=20)


@pytest.mark.parametrize(
    "residues,shared",
    [
        (("ASP", "GLU"), True),      # acidic
        (("ASN", "GLN"), True),      # amide
        (("ALA", "VAL"), True),      # small hydrophobic
        (("LEU", "ILE"), True),      # branched hydrophobic
        (("PHE", "TRP"), True),      # aromatic
        (("SER", "THR"), True),      # hydroxyl
        (("LYS", "ARG"), True),      # basic
        (("MET", "MSE"), True),      # selenomethionine with methionine
        (("GLY", "PRO"), False),     # each its own group
        (("CYS", "MET"), False),     # sulphur, but not grouped together
        (("TYR", "PHE"), False),     # TYR is separated from the other aromatics
        (("HIS", "LYS"), False),     # HIS is its own group
        (("SEP", "SER"), False),     # a phosphorylated residue is its own group
        (("TPO", "THR"), False),
        (("PTR", "TYR"), False),
        (("SEP", "TPO"), False),
        (("CME", "CSO"), False),
        (("CSD", "CYS"), False),
        (("M3L", "LYS"), False),
        (("LLP", "LYS"), False),
        (("NIY", "TYR"), False),
        (("TRF", "TRP"), False),
    ],
)
def test_grouping_map_boundaries(residues, shared):
    a, b = (AMINO_ACID_GROUP_23.label(r) for r in residues)
    assert (a == b) is shared


def test_every_residue_the_base_scheme_labels_has_a_group():
    assert set(AMINO_ACID_GROUP_23.table) == set(AMINO_ACID_30.table)


def test_both_residue_schemes_fit_one_radix_30_window():
    """Collision-free by the window rule: every code within 4..33 at radix 30.
    The grouped codes are not consecutive (17..23 are unused), and that is
    harmless: what matters is the span."""
    assert AMINO_ACID_30.code_range == (4, 33)
    assert AMINO_ACID_30.codes_are_consecutive
    assert AMINO_ACID_GROUP_23.code_range == (4, 33)
    assert not AMINO_ACID_GROUP_23.codes_are_consecutive
    assert AMINO_ACID_GROUP_23.encoding_is_injective
    assert AMINO_ACID_30.encoding_is_injective


def test_grouping_changes_the_key_space_fingerprint_at_the_same_radix():
    """The two residue schemes share radix 30 and remain different key spaces:
    the scheme name and content fingerprint differ."""
    plain = STANDARD.key_space()
    grouped = STANDARD.with_modifiers(aa_grouping=True).key_space()
    assert plain.fingerprint != grouped.fingerprint
    assert plain.label_scheme.radix == grouped.label_scheme.radix == 30
    assert plain.label_scheme.name != grouped.label_scheme.name
    assert plain.label_scheme.content_fingerprint != grouped.label_scheme.content_fingerprint


def test_changing_the_grouping_map_changes_the_fingerprint():
    altered = dict(AMINO_ACID_GROUP_23.table)
    altered["CYS"] = altered["MET"]  # merge two groups -> 22 distinct codes
    other = LabelScheme(
        name="amino_acid_group_23",
        table=altered,
        n_labels=30,
        domain=LabelDomain.PROTEIN_AMINO_ACID,
    )
    assert other.content_fingerprint() != AMINO_ACID_GROUP_23.content_fingerprint()


def test_changing_only_the_radix_changes_the_fingerprint():
    radix_40 = LabelScheme(
        name="amino_acid_group_23",
        table=dict(AMINO_ACID_GROUP_23.table),
        n_labels=40,
        domain=LabelDomain.PROTEIN_AMINO_ACID,
    )
    assert AMINO_ACID_GROUP_23.table == radix_40.table
    assert AMINO_ACID_GROUP_23.content_fingerprint() != radix_40.content_fingerprint()


def test_grouping_is_rejected_outside_the_protein_amino_acid_domain():
    """AA grouping is Ca-only, checked against the label DOMAIN.

    Atomic/drug/carbohydrate/nucleotide points carry residue metadata but are
    labelled by element, so grouping is meaningless there and must
    fail clearly rather than be silently ignored.
    """
    atomic = LabelScheme(
        name="element_labels_112",
        table={"C": 1, "N": 2, "O": 3},
        n_labels=3,
        domain=LabelDomain.ATOM_ELEMENT,
    )
    with pytest.raises(IncompatibleModifier, match="Protein Ca label domain"):
        grouping_mod.check_applicable(atomic)


def test_the_rejection_names_the_required_domain_and_the_offending_one():
    atomic = LabelScheme(
        name="element_labels_112",
        table={"C": 1},
        n_labels=1,
        domain=LabelDomain.ATOM_ELEMENT,
    )
    with pytest.raises(IncompatibleModifier) as exc:
        grouping_mod.apply_grouping(atomic, AMINO_ACID_GROUP_23)
    message = str(exc.value)
    assert LabelDomain.PROTEIN_AMINO_ACID in message
    assert LabelDomain.ATOM_ELEMENT in message
    assert "element_labels_112" in message


def test_a_grouped_scheme_missing_a_residue_is_refused():
    partial = LabelScheme(
        name="partial_groups",
        table={"LEU": 4, "ILE": 4},
        n_labels=1,
        domain=LabelDomain.PROTEIN_AMINO_ACID,
    )
    with pytest.raises(IncompatibleModifier, match="no group for"):
        grouping_mod.apply_grouping(AMINO_ACID_30, partial)


def test_a_grouped_scheme_without_the_modified_residues_is_refused():
    """The historical 20-residue grouping cannot serve the residue labels: the
    base scheme labels eleven modified residues it has no group for."""
    regular_only = LabelScheme(
        name="groups_of_the_20_regular_residues",
        table={r: c for r, c in AMINO_ACID_GROUP_23.table.items() if c <= 16 and r != "MSE"},
        n_labels=13,
        domain=LabelDomain.PROTEIN_AMINO_ACID,
    )
    with pytest.raises(IncompatibleModifier, match="no group for 11 residue"):
        grouping_mod.apply_grouping(AMINO_ACID_30, regular_only)


def test_grouping_hook_is_the_label_stage():
    assert grouping_mod.HOOK == ModifierHook.LABELS


# ===========================================================================
# Mirror
# ===========================================================================


def test_mirror_is_off_by_default():
    assert STANDARD.mirror_convention is MIRROR_DISABLED
    assert STANDARD.mirror_convention.enabled is False


def test_mirror_disabled_never_negates_under_any_angle():
    for theta_raw in (0.0, 45.0, 89.9, 90.0, 90.1, 179.9, 180.0):
        assert MIRROR_DISABLED.apply(12345, theta_raw) == 12345
        assert MIRROR_DISABLED.category_for(theta_raw) == DEFAULT_CATEGORY


@pytest.mark.parametrize("theta_raw", [0.0, 1.0, 45.0, 89.999, 90.0])
def test_mirror_on_and_at_or_below_90_is_positive(theta_raw):
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(12345, theta_raw) == 12345
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.category_for(theta_raw) == DEFAULT_CATEGORY


@pytest.mark.parametrize("theta_raw", [90.0000001, 91.0, 137.49, 180.0])
def test_mirror_on_and_above_90_is_negative(theta_raw):
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(12345, theta_raw) == -12345
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.category_for(theta_raw) == MIRROR_CATEGORY


def test_exactly_90_is_the_boundary_and_is_not_mirrored():
    """Strictly greater than, never `>=`."""
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(7, 90.0) == 7
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(7, 90.0 + 1e-12) == -7


def test_the_zero_key_edge_case_is_preserved_and_the_category_closes_it():
    """key == 0 negates to itself, so sign alone cannot say it was mirrored.

    The category field carries the information the sign cannot -- which is why
    category is reported explicitly rather than inferred from the sign.
    """
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(0, 137.0) == 0
    assert MIRROR_NEGATE_WHEN_PREFOLD_GT_90.category_for(137.0) == MIRROR_CATEGORY


def test_mirror_is_a_pure_sign_multiplier():
    for key in (1, 12345, 987654321):
        for theta_raw in (10.0, 100.0):
            signed = MIRROR_NEGATE_WHEN_PREFOLD_GT_90.apply(key, theta_raw)
            assert abs(signed) == key


def test_mirror_changes_the_key_space_fingerprint():
    plain = STANDARD.key_space()
    mirrored = STANDARD.with_modifiers(mirror=True).key_space()
    assert plain.fingerprint != mirrored.fingerprint
    assert mirrored.mirror_convention == "negate_when_prefold_angle_gt_90"


def test_mirror_does_not_change_the_sampling_profile():
    """Mirror changes what a key means, not which triangles were drawn."""
    from biotsr.api import _sampling_profile
    from biotsr.methods.registry import get_method

    recipe = get_method("protein")
    plain = _sampling_profile(recipe, STANDARD)
    mirrored = _sampling_profile(recipe, STANDARD.with_modifiers(mirror=True))
    assert plain.fingerprint == mirrored.fingerprint


# ===========================================================================
# Size filter
# ===========================================================================


def test_no_filter_by_default():
    assert STANDARD.size_filter is None
    assert SizeFilter.from_value(None) is None


def test_the_characterized_criterion_is_strict_less_than():
    f = SizeFilter(cutoff=25.0)
    assert f.comparison == Comparison.STRICT_LESS_THAN
    assert f.describe() == "max_distance < 25.0"


@pytest.mark.parametrize(
    "max_distance,accepted",
    [
        (0.0, True),
        (24.999999, True),
        (25.0, False),      # exactly at the cutoff: EXCLUDED, per the characterization
        (25.000001, False),
        (1000.0, False),
    ],
)
def test_size_filter_boundary_matches_the_characterized_legacy_semantics(
    max_distance, accepted
):
    assert SizeFilter(cutoff=25.0).accepts(max_distance) is accepted


def test_the_inclusive_variant_exists_but_is_not_the_legacy_default():
    inclusive = SizeFilter(cutoff=25.0, comparison=Comparison.LESS_THAN_OR_EQUAL)
    assert inclusive.accepts(25.0) is True
    assert SizeFilter(cutoff=25.0).accepts(25.0) is False


def test_a_nonsensical_cutoff_is_refused():
    for bad in (0, -1.0, -25):
        with pytest.raises(ValueError, match="positive"):
            SizeFilter(cutoff=bad)
    for bad in ("25", None, True):
        with pytest.raises(ValueError):
            SizeFilter(cutoff=bad)


def test_unknown_comparison_is_refused():
    with pytest.raises(ValueError, match="comparison"):
        SizeFilter(cutoff=25.0, comparison="roughly_less_than")


def test_size_filter_changes_sampling_not_key_space():
    from biotsr.api import _sampling_profile
    from biotsr.methods.registry import get_method

    recipe = get_method("protein")
    filtered_profile = STANDARD.with_modifiers(size_filter=25.0)

    assert STANDARD.key_space().fingerprint == filtered_profile.key_space().fingerprint
    plain = _sampling_profile(recipe, STANDARD)
    filtered = _sampling_profile(recipe, filtered_profile)
    assert plain.fingerprint != filtered.fingerprint
    assert filtered.size_filter == 25.0


def test_size_filter_reads_only_max_distance():
    """What makes it reusable by every future family without change."""
    import inspect

    source = inspect.getsource(SizeFilter.accepts)
    for protein_specific in ("residue", "chain", "amino", "calpha", "label"):
        assert protein_specific not in source.lower()
