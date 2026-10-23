"""Cross KeySpace vs SamplingProfile, label space, and the 18-class distance scale."""

from __future__ import annotations

from dataclasses import replace

import pytest

from biotsr.binning.distance import DISTANCE_SCHEME_ATOMIC_18
from biotsr.cross.candidates import COMPOSITION_RULE, STRATEGY_MIXED_SIDE
from biotsr.cross.pools import (
    ContactFilter,
    ContactPoolPolicy,
    PoolIdentity,
    PoolOrder,
    PoolProjection,
)
from biotsr.domain.keyspace import SamplingProfile
from biotsr.labels.schemes import LabelDomain
from biotsr.methods.registry import get_method
from biotsr.profiles.atomic import ATOMIC_PROFILES
from biotsr.profiles.cross import (
    CROSS_PROFILES,
    STANDARD_ATOMIC_CROSS_V1,
    STANDARD_CALPHA_CROSS_V1,
)

CROSS_METHODS = (
    "drug_protein_cross",
    "carbohydrate_protein_cross",
    "nucleotide_protein_cross",
    "protein_protein_cross",
)


# ===========================================================================
# Label space
# ===========================================================================


@pytest.mark.parametrize("name", CROSS_METHODS)
def test_every_cross_method_defaults_to_radix_112(name):
    """Fully atomic cross means the 112-radix element labels. Not 20, not 13,
    not 109 -- and 109 in particular is the number of *distinct labels* the
    table uses, never the radix the encoder multiplies by."""
    scheme = get_method(name).default_profile.label_scheme
    assert scheme.n_labels == 112
    assert scheme.domain == LabelDomain.ATOM_ELEMENT


#: Every cross profile except the Ca one. The Ca profile is a residue profile by
#: design: an interface residue is represented by its Ca, so it labels in the
#: 20-code amino-acid space -- the same key space as ordinary protein Ca TSR, and
#: deliberately NOT the atomic 112 space (ADR-052).
ATOMIC_CROSS_PROFILES = [
    p for p in CROSS_PROFILES if p.name != "standard_calpha_cross_v1"
]


@pytest.mark.parametrize("profile", ATOMIC_CROSS_PROFILES)
def test_no_atomic_cross_profile_can_reach_a_calpha_label_scheme(profile):
    assert profile.label_scheme.n_labels == 112
    assert "amino_acid" not in profile.label_scheme.name


def test_the_calpha_cross_profile_is_the_ordinary_calpha_key_space():
    """A Ca cross key and an ordinary protein Ca key must be directly comparable:
    same labels, same distance scale, same everything a key means.  Only which
    residues were drawn differs, and that is a SamplingProfile fact."""
    from biotsr.profiles.cross import STANDARD_CALPHA_CROSS_V1
    from biotsr.profiles.standard import STANDARD_PROTEIN_CALPHA_V1

    cross = STANDARD_CALPHA_CROSS_V1.key_space()
    plain = STANDARD_PROTEIN_CALPHA_V1.key_space()
    assert cross.label_scheme.radix == 30
    assert cross.distance_scheme.radix == 35
    cross.require_compatible(plain)


def test_the_standard_cross_key_space_matches_the_standard_atomic_one():
    """A standard cross key and a standard non-cross atomic key are built on the
    same scales: 29 theta classes, 18 distance classes, radix 112.  Only the
    *sampling* differs, which is the whole point of the KeySpace/SamplingProfile
    split."""
    from biotsr.profiles.atomic import STANDARD_DRUG_V1

    a = STANDARD_ATOMIC_CROSS_V1.key_space()
    b = STANDARD_DRUG_V1.key_space()
    assert a.theta_scheme == b.theta_scheme
    assert a.distance_scheme == b.distance_scheme
    assert a.label_scheme == b.label_scheme
    assert a.fingerprint == b.fingerprint


# ===========================================================================
# The standard atomic cross distance scale: 18 classes
# ===========================================================================


#: Historical constant, documentation only: `TSR-Package`'s `Drug_TSR.py` /
#: `Cross_TSR.py` binned cross distances on a 17-class ladder (``dLen = 17``,
#: sixteen 1 A cut points plus an overflow). The standard does not adopt it.
HISTORICAL_TSR_PACKAGE_CROSS_DISTANCE_CLASSES = 17


def test_the_standard_cross_profile_uses_the_eighteen_class_atomic_scale():
    """The standard uses the atomic family's 18-class scheme (ADR-021)."""
    scheme = STANDARD_ATOMIC_CROSS_V1.distance_scheme
    assert scheme is DISTANCE_SCHEME_ATOMIC_18
    assert scheme.n_classes == 18
    assert scheme.n_classes != HISTORICAL_TSR_PACKAGE_CROSS_DISTANCE_CLASSES
    # The two scales part company above 17 A: the standard keeps a class there.
    assert scheme.classify(16.5) == 17
    assert scheme.classify(17.5) == 18


def test_every_atomic_profile_uses_the_eighteen_class_scale():
    for profile in tuple(ATOMIC_PROFILES) + tuple(ATOMIC_CROSS_PROFILES):
        assert profile.distance_scheme is DISTANCE_SCHEME_ATOMIC_18, profile.name


# ===========================================================================
# KeySpace vs SamplingProfile
# ===========================================================================


def sampling(**kw):
    base = dict(
        selector="cross_entity_atoms|cross_protein_chain_atoms",
        scope="cross_two_sided",
        candidate_strategy=STRATEGY_MIXED_SIDE,
        composition_rule=COMPOSITION_RULE,
        cross_method="nucleotide_protein_cross",
        contact_cutoff=3.0,
        contact_filter=ContactFilter.BOTH_SIDES,
        pool_identity=PoolIdentity.FULL_ATOM_IDENTITY,
        pool_order=PoolOrder.BLOCKED_BY_SIDE,
    )
    base.update(kw)
    return SamplingProfile(**base)


@pytest.mark.parametrize(
    "field,value",
    [
        ("contact_cutoff", 4.5),
        ("contact_filter", ContactFilter.SIDE_A_ONLY),
        ("contact_projection", PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA),
        ("side_a_selector", "cross_entity_atoms chain=E entity=DG_8"),
        ("cross_method", "carbohydrate_protein_cross"),
        ("composition_rule", None),
    ],
)
def test_every_cross_sampling_fact_moves_the_sampling_fingerprint(field, value):
    base = sampling()
    assert replace(base, **{field: value}).fingerprint != base.fingerprint


def test_changing_the_contact_cutoff_does_not_redefine_the_key_space():
    """A 3.0 A pool and a 3.5 A pool over the same structure produce keys that
    are directly comparable -- they were just sampled differently.  The cutoff is
    therefore a SamplingProfile fact and must never enter KeySpace (ADR-043)."""
    fingerprint = STANDARD_ATOMIC_CROSS_V1.key_space().fingerprint
    for cutoff in (2.0, 3.0, 3.5, 5.0, 12.0):
        policy = ContactPoolPolicy(name="p", cutoff=cutoff)
        assert policy.cutoff == cutoff
        # building a policy cannot touch the profile at all
        assert STANDARD_ATOMIC_CROSS_V1.key_space().fingerprint == fingerprint
    a, b = sampling(contact_cutoff=3.0), sampling(contact_cutoff=3.5)
    assert a.fingerprint != b.fingerprint  # sampling differs ...
    # ... and the key space is untouched, since neither carries one.
    assert "contact_cutoff" not in STANDARD_ATOMIC_CROSS_V1.key_space().to_dict()


def test_the_pool_identity_is_sampling_not_key_semantics():
    """Which atoms were pooled is a sampling fact, recorded in SamplingProfile;
    it never enters what a key means."""
    assert "pool_identity" not in STANDARD_ATOMIC_CROSS_V1.key_space().to_dict()
    assert "pool_order" not in STANDARD_ATOMIC_CROSS_V1.key_space().to_dict()
    assert sampling().to_dict()["pool_identity"] == PoolIdentity.FULL_ATOM_IDENTITY


def test_the_key_space_carries_only_key_defining_semantics():
    keys = set(STANDARD_ATOMIC_CROSS_V1.key_space().to_dict())
    assert keys == {
        "angle_convention",
        "sqrt_convention",
        "median_convention",
        "numerical_precision",
        "theta_scheme",
        "distance_scheme",
        "label_scheme",
        "canonical_ordering_version",
        "key_formula_version",
        "mirror_convention",
    }


@pytest.mark.parametrize(
    "field",
    ["cross_method", "side_a_selector", "side_b_selector", "contact_cutoff",
     "contact_filter", "pool_identity", "pool_order", "composition_rule"],
)
def test_the_sampling_profile_carries_every_cross_specific_fact(field):
    assert field in sampling().to_dict()


# ===========================================================================
# The two cross key spaces are mutually incomparable
# ===========================================================================


def test_the_two_cross_key_spaces_are_different():
    prints = {p.name: p.key_space().fingerprint for p in CROSS_PROFILES}
    assert len(prints) == 2
    assert len(set(prints.values())) == len(prints)


def test_the_guard_names_the_differences_rather_than_just_refusing():
    from biotsr.domain.keyspace import IncompatibleKeySpace

    with pytest.raises(IncompatibleKeySpace) as excinfo:
        STANDARD_ATOMIC_CROSS_V1.key_space().require_compatible(
            STANDARD_CALPHA_CROSS_V1.key_space()
        )
    message = str(excinfo.value)
    assert "distance_scheme" in message
    assert "label_scheme" in message


def test_the_standard_cross_profiles_use_math_pi():
    """Historical note, documentation only: the `cross_key.py` lineage divided
    by 3.14, while `TSR-Package`'s cross programs used ``math.pi``.  The standard
    uses ``math.pi`` for every cross profile."""
    for profile in CROSS_PROFILES:
        assert profile.angle_convention.name == "pi", profile.name
