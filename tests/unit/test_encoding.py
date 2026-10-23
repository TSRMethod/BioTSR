"""Key encoding: hand-calculated keys, collision-free labels, and radix behaviour."""

from __future__ import annotations

import itertools

import pytest

from biotsr.encoding.key import KEY_FORMULA_VERSION, KeyFormulaSpec, encode
from biotsr.labels.schemes import (
    AMINO_ACID_30,
    AMINO_ACID_GROUP_23,
    LabelScheme,
    label_triple_collisions,
    ordered_label_triples,
)

CALPHA_SPEC = KeyFormulaSpec(n_labels=30, n_theta=29, n_distance=35)
#: The radix-20 formula, kept to state the arithmetic of the historical
#: 20-residue keys and the collision the old radix produces with 30 codes.
RADIX_20_SPEC = KeyFormulaSpec(n_labels=20, n_theta=29, n_distance=35)


def test_hand_calculated_key_at_radix_30():
    """MET(13) ILE(10) GLN(6), theta class 29, distance class 2:

        35*29*900*(13-1) = 10,962,000
      + 35*29*30*(10-1)  =    274,050
      + 35*29*(6-1)      =      5,075
      + 29*(2-1)         =         29
      + (29-1)           =         28
      -------------------------------
                           11,241,182
    """
    assert encode((13, 10, 6), 29, 2, CALPHA_SPEC) == 11241182


def test_second_hand_calculated_key_at_radix_30():
    """MET(13) GLN(6) PHE(4), theta class 5, distance class 4 -> 11,117,386."""
    assert encode((13, 6, 4), 5, 4, CALPHA_SPEC) == 11117386


def test_the_same_labels_at_radix_20_give_the_historical_keys():
    """The 1UBQ golden's first two rows, under the historical radix 20: the
    formula is unchanged, only ``L`` differs."""
    assert encode((13, 10, 6), 29, 2, RADIX_20_SPEC) == 5059832
    assert encode((13, 6, 4), 5, 4, RADIX_20_SPEC) == 4976636


def test_hand_calculated_modified_residue_keys_are_distinct():
    """The pair that collided at radix 20, in the synthetic geometry where both
    fall in theta class 26 and distance class 3:

      SEP-PHE-PHE (25,4,4): 1015*(900*24 + 30*3 + 3) + 29*2 + 25 = 22,018,478
      TPO-TPO-PHE (24,24,4): 1015*(900*23 + 30*23 + 3) + 29*2 + 25 = 21,713,978
    """
    assert encode((25, 4, 4), 26, 3, CALPHA_SPEC) == 22018478
    assert encode((24, 24, 4), 26, 3, CALPHA_SPEC) == 21713978


def test_formula_is_written_out_longhand():
    labels, tclass, dclass = (33, 17, 9), 12, 30
    d_len, d_theta, n = 35, 29, 30
    expected = (
        d_len * d_theta * (n ** 2) * (labels[0] - 1)
        + d_len * d_theta * n * (labels[1] - 1)
        + d_len * d_theta * (labels[2] - 1)
        + d_theta * (dclass - 1)
        + (tclass - 1)
    )
    assert encode(labels, tclass, dclass, CALPHA_SPEC) == expected


def test_theta_is_the_lowest_order_digit():
    a = encode((13, 10, 6), 1, 2, CALPHA_SPEC)
    b = encode((13, 10, 6), 2, 2, CALPHA_SPEC)
    assert b - a == 1


def test_distance_is_the_next_digit_up():
    a = encode((13, 10, 6), 1, 1, CALPHA_SPEC)
    b = encode((13, 10, 6), 1, 2, CALPHA_SPEC)
    assert b - a == 29


def test_label_order_matters():
    assert encode((13, 10, 6), 5, 5, CALPHA_SPEC) != encode(
        (6, 10, 13), 5, 5, CALPHA_SPEC
    )


@pytest.mark.parametrize("scheme, triples", [(AMINO_ACID_30, 4960),
                                             (AMINO_ACID_GROUP_23, 2300)])
def test_every_ordered_label_triple_has_its_own_label_component(scheme, triples):
    """Exhaustive: every triple the canonical ordering can emit (c1 >= c2 >= c3)
    -- 4,960 for the 30 residue codes, 2,300 for the 23 grouped codes -- has a
    label component no other triple shares, at radix 30."""
    assert len(ordered_label_triples(scheme)) == triples
    assert label_triple_collisions(scheme) == {}
    components = {encode(t, 1, 1, CALPHA_SPEC) for t in ordered_label_triples(scheme)}
    assert len(components) == triples


def test_the_bin_part_never_reaches_the_label_digits():
    """Every (theta, distance) pair gives a distinct bin part below 35*29 =
    1,015, the step of one label component: with distinct label components
    (above), full keys are distinct for every triple and every bin pair. The
    smallest and largest keys are the boundary classes of the extreme codes."""
    base = encode((4, 4, 4), 1, 1, CALPHA_SPEC)
    parts = {encode((4, 4, 4), t, d, CALPHA_SPEC) - base
             for t in range(1, 30) for d in range(1, 36)}
    assert parts == set(range(1015))
    assert base == 2834895
    assert encode((33, 33, 33), 29, 35, CALPHA_SPEC) == 30239894
    assert (encode((33, 33, 33), 29, 35, CALPHA_SPEC)
            - encode((33, 33, 33), 1, 1, CALPHA_SPEC)) == 1014


def test_the_radix_20_collision_is_gone_in_every_bin_pair():
    """At radix 20, SEP-PHE-PHE (25,4,4) and TPO-TPO-PHE (24,24,4) share the
    label component 9663 (400*24 + 20*3 + 3 = 400*23 + 20*23 + 3), so they
    collide for every bin pair. At radix 30 they never do."""
    for theta_class, distance_class in itertools.product(range(1, 30), range(1, 36)):
        assert (encode((25, 4, 4), theta_class, distance_class, RADIX_20_SPEC)
                == encode((24, 24, 4), theta_class, distance_class, RADIX_20_SPEC))
        assert (encode((25, 4, 4), theta_class, distance_class, CALPHA_SPEC)
                != encode((24, 24, 4), theta_class, distance_class, CALPHA_SPEC))


def test_intentional_grouping_still_shares_codes():
    """Grouping is deliberate: MSE shares MET's code in both tables, and the
    grouped table puts SER and THR (and ASP and GLU, ...) in one group, so their
    triangles share keys on purpose."""
    assert AMINO_ACID_30.label("MSE") == AMINO_ACID_30.label("MET") == 13
    assert AMINO_ACID_GROUP_23.label("MSE") == AMINO_ACID_GROUP_23.label("MET") == 14
    assert AMINO_ACID_GROUP_23.label("SER") == AMINO_ACID_GROUP_23.label("THR") == 12
    assert AMINO_ACID_GROUP_23.label("SEP") != AMINO_ACID_GROUP_23.label("SER")
    ser = AMINO_ACID_GROUP_23.label("SER")
    asp = AMINO_ACID_GROUP_23.label("ASP")
    assert encode((ser, ser, asp), 7, 7, CALPHA_SPEC) == encode(
        (AMINO_ACID_GROUP_23.label("THR"), ser, AMINO_ACID_GROUP_23.label("GLU")), 7, 7,
        CALPHA_SPEC)


def test_changing_a_radix_remaps_every_key():
    """There is no partial compatibility between radix triples."""
    labels, tclass, dclass = (13, 10, 6), 7, 7
    assert encode(labels, tclass, dclass, CALPHA_SPEC) != encode(
        labels, tclass, dclass, RADIX_20_SPEC
    )


def test_out_of_range_classes_are_rejected():
    with pytest.raises(ValueError, match="theta_class"):
        encode((13, 10, 6), 30, 2, CALPHA_SPEC)
    with pytest.raises(ValueError, match="theta_class"):
        encode((13, 10, 6), 0, 2, CALPHA_SPEC)
    with pytest.raises(ValueError, match="distance_class"):
        encode((13, 10, 6), 1, 36, CALPHA_SPEC)
    with pytest.raises(ValueError, match="distance_class"):
        encode((13, 10, 6), 1, 0, CALPHA_SPEC)


def test_wrong_label_arity_is_rejected():
    with pytest.raises(ValueError, match="three ordered labels"):
        encode((13, 10), 1, 1, CALPHA_SPEC)


def test_spec_rejects_non_positive_radices():
    for kwargs in (
        {"n_labels": 0, "n_theta": 29, "n_distance": 35},
        {"n_labels": 20, "n_theta": -1, "n_distance": 35},
        {"n_labels": 20, "n_theta": 29, "n_distance": 0},
    ):
        with pytest.raises(ValueError):
            KeyFormulaSpec(**kwargs)


def test_formula_version_is_stated():
    assert KEY_FORMULA_VERSION == "mixed_radix_v1"
    assert CALPHA_SPEC.version == KEY_FORMULA_VERSION


def test_the_encoder_refuses_a_label_outside_the_schemes_codes():
    """With the scheme's code range (every profile supplies it), a label the
    scheme cannot emit is refused rather than carried into the next digit."""
    spec = KeyFormulaSpec(n_labels=30, n_theta=29, n_distance=35, label_range=(4, 33))
    assert encode((33, 4, 4), 1, 1, spec) == encode((33, 4, 4), 1, 1, CALPHA_SPEC)
    for bad in ((34, 4, 4), (13, 3, 4), (13, 10, 0)):
        with pytest.raises(ValueError, match=r"outside the label scheme's codes 4\.\.33"):
            encode(bad, 1, 1, spec)


def test_a_label_range_wider_than_the_radix_is_refused():
    with pytest.raises(ValueError, match="does not fit in one radix-20 window"):
        KeyFormulaSpec(n_labels=20, n_theta=29, n_distance=35, label_range=(4, 33))
    with pytest.raises(ValueError):
        KeyFormulaSpec(n_labels=30, n_theta=29, n_distance=35, label_range=(0, 29))


def test_every_profile_states_its_label_range():
    import biotsr

    for name in biotsr.list_methods():
        profile = biotsr.get_method(name).default_profile
        assert profile.key_formula().label_range == profile.label_scheme.code_range
        grouped = profile.grouped_label_scheme
        assert grouped.code_range == (4, 33)


def test_label_scheme_refuses_codes_outside_one_radix_window():
    """The supplied tables at the old radices are refused before any key is
    computed: 4..33 does not fit in 20 slots, nor in 13."""
    with pytest.raises(ValueError,
                       match=r"codes 4\.\.33 \(a span of 30\) but declares radix 20"):
        LabelScheme(name="residues_at_20", table=dict(AMINO_ACID_30.table), n_labels=20)
    with pytest.raises(ValueError,
                       match=r"codes 4\.\.33 \(a span of 30\) but declares radix 13"):
        LabelScheme(name="groups_at_13", table=dict(AMINO_ACID_GROUP_23.table), n_labels=13)


def test_the_code_count_is_not_the_test():
    """The rule this release replaces compared the number of distinct codes with
    the radix. The grouped table has 23 codes, so radix 23 passed that check --
    yet its codes span 4..33 (17..23 are unused), so 98 label components
    collide at radix 23. The window rule refuses it; radix 30 is collision-free.
    """
    table = dict(AMINO_ACID_GROUP_23.table)
    assert len(set(table.values())) == 23
    with pytest.raises(ValueError, match="span of 30.*radix 23.*is not the test"):
        LabelScheme(name="groups_at_23", table=table, n_labels=23)
    # What radix 23 would have done, with the encoder itself:
    spec_23 = KeyFormulaSpec(n_labels=23, n_theta=1, n_distance=1)
    components = {}
    for triple in ordered_label_triples(AMINO_ACID_GROUP_23):
        components.setdefault(encode(triple, 1, 1, spec_23), []).append(triple)
    assert sum(1 for v in components.values() if len(v) > 1) == 98
    assert label_triple_collisions(AMINO_ACID_GROUP_23) == {}


def test_a_radix_wider_than_the_codes_is_accepted():
    """The element labels declare radix 112 for 109 codes over 4..112: they fit
    in one window, so the encoding is collision-free and the scheme is valid."""
    scheme = LabelScheme(name="wide", table={"A": 4, "B": 5}, n_labels=20)
    assert scheme.encoding_is_injective
    assert label_triple_collisions(scheme) == {}


def test_a_code_below_one_is_refused():
    with pytest.raises(ValueError, match="integer >= 1"):
        LabelScheme(name="zero", table={"A": 0, "B": 1}, n_labels=2)


def test_the_old_radix_20_collisions_are_what_the_window_rule_prevents():
    """Recorded so the reason for radix 30 stays visible: the supplied residue
    codes at radix 20 give 255 label components shared by two or more triples,
    every one involving an added code; the 1,540 triples of the 20 regular codes
    are unaffected."""
    spec_20 = KeyFormulaSpec(n_labels=20, n_theta=1, n_distance=1)
    triples = ordered_label_triples(AMINO_ACID_30)
    by_component = {}
    for triple in triples:
        by_component.setdefault(encode(triple, 1, 1, spec_20), []).append(triple)
    groups = [v for v in by_component.values() if len(v) > 1]
    assert len(groups) == 255
    assert all(any(c >= 24 for t in g for c in t) for g in groups)
    regular = [t for t in triples if max(t) <= 23]
    assert len(regular) == 1540
    assert all(len(by_component[encode(t, 1, 1, spec_20)]) == 1 for t in regular)
