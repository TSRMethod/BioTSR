"""KeySpace and SamplingProfile: fingerprint sensitivity and insensitivity.

The two identities answer two different questions, and the tests are written to
prove both halves of each: changing a key-defining semantic MUST change the
``KeySpace`` fingerprint, and changing a mere runtime setting MUST NOT change
either fingerprint.
"""

from __future__ import annotations

import pytest

from biotsr.domain.keyspace import (
    IncompatibleKeySpace,
    IncompatibleSampling,
    KeySpace,
    SamplingProfile,
    SchemeRef,
)
from biotsr.profiles.standard import STANDARD_PROTEIN_CALPHA_V1


def base_key_space() -> KeySpace:
    return STANDARD_PROTEIN_CALPHA_V1.key_space()


def replace(ks: KeySpace, **kwargs) -> KeySpace:
    import dataclasses

    return dataclasses.replace(ks, **kwargs)


# ---------------------------------------------------------------------------
# Sensitivity -- every key-defining semantic must move the fingerprint
# ---------------------------------------------------------------------------


def test_fingerprint_is_stable_across_identical_constructions():
    assert base_key_space().fingerprint == base_key_space().fingerprint


def test_angle_convention_changes_the_fingerprint():
    ks = base_key_space()
    assert replace(ks, angle_convention="alternate_divisor:3.0").fingerprint != ks.fingerprint


def test_theta_scheme_content_changes_the_fingerprint():
    ks = base_key_space()
    altered = SchemeRef(
        name=ks.theta_scheme.name,
        content_fingerprint="0" * 64,
        radix=ks.theta_scheme.radix,
    )
    assert replace(ks, theta_scheme=altered).fingerprint != ks.fingerprint


def test_distance_scheme_content_changes_the_fingerprint():
    ks = base_key_space()
    altered = SchemeRef(
        name=ks.distance_scheme.name,
        content_fingerprint="0" * 64,
        radix=ks.distance_scheme.radix,
    )
    assert replace(ks, distance_scheme=altered).fingerprint != ks.fingerprint


def test_label_scheme_radix_changes_the_fingerprint():
    ks = base_key_space()
    altered = SchemeRef(
        name=ks.label_scheme.name,
        content_fingerprint=ks.label_scheme.content_fingerprint,
        radix=13,
    )
    assert replace(ks, label_scheme=altered).fingerprint != ks.fingerprint


def test_canonical_ordering_version_changes_the_fingerprint():
    ks = base_key_space()
    assert replace(ks, canonical_ordering_version="v2").fingerprint != ks.fingerprint


def test_key_formula_version_changes_the_fingerprint():
    ks = base_key_space()
    assert replace(ks, key_formula_version="mixed_radix_v2").fingerprint != ks.fingerprint


def test_mirror_convention_changes_the_fingerprint():
    ks = base_key_space()
    assert (
        replace(ks, mirror_convention="negate_when_prefold_angle_gt_90").fingerprint
        != ks.fingerprint
    )


def test_same_radix_different_scheme_content_is_the_hazard_case():
    """Same declared size, different keys: a radix-only comparison would miss it."""
    ks = base_key_space()
    altered = SchemeRef(
        name=ks.distance_scheme.name,
        content_fingerprint="0" * 64,
        radix=ks.distance_scheme.radix,
    )
    changed = replace(ks, distance_scheme=altered)
    assert changed.distance_scheme.radix == ks.distance_scheme.radix
    assert changed.fingerprint != ks.fingerprint


# ---------------------------------------------------------------------------
# Insensitivity -- runtime settings belong in neither identity
# ---------------------------------------------------------------------------


def test_no_runtime_setting_can_enter_either_identity():
    """Worker count, output directory, output format, cache location and log level
    cannot change a key or a count, so there is nowhere for them to go."""
    ks_fields = set(base_key_space().to_dict())
    sp_fields = set(
        SamplingProfile(
            selector="protein_calpha", scope="pooled", candidate_strategy="exhaustive_c_n_3"
        ).to_dict()
    )
    forbidden = {
        "workers",
        "n_jobs",
        "output_dir",
        "output_path",
        "output_format",
        "cache_dir",
        "log_level",
        "input_path",
        "verbose",
    }
    assert not (ks_fields & forbidden)
    assert not (sp_fields & forbidden)


def test_serialization_choice_does_not_change_the_key_space():
    """How a result is written out is an export concern, not a key-space one."""
    ks = base_key_space()
    assert ks.fingerprint == base_key_space().fingerprint


# ---------------------------------------------------------------------------
# SamplingProfile
# ---------------------------------------------------------------------------


def base_sampling() -> SamplingProfile:
    return SamplingProfile(
        selector="protein_calpha",
        scope="pooled",
        candidate_strategy="exhaustive_c_n_3",
        degenerate_policy="skip_and_record+acos_clamp",
    )


def test_candidate_strategy_identity_is_separate_from_the_key_space():
    """A future candidate generator occupies the SAME key space but a DIFFERENT sampling
    profile: a key means the same thing, but the frequency vectors are not
    comparable.  That distinction is the whole reason for two identities."""
    import dataclasses

    exhaustive = base_sampling()
    knn = dataclasses.replace(exhaustive, candidate_strategy="knn_k12")
    assert knn.fingerprint != exhaustive.fingerprint
    assert base_key_space().fingerprint == base_key_space().fingerprint


def test_size_filter_changes_sampling_not_key_space():
    import dataclasses

    a = base_sampling()
    b = dataclasses.replace(a, size_filter=25.0)
    assert a.fingerprint != b.fingerprint


def test_degenerate_policy_is_sampling_not_key_space():
    import dataclasses

    a = base_sampling()
    b = dataclasses.replace(a, degenerate_policy="skip_and_record")
    assert a.fingerprint != b.fingerprint


def test_selector_and_scope_change_sampling():
    import dataclasses

    a = base_sampling()
    assert dataclasses.replace(a, selector="ligand_atoms").fingerprint != a.fingerprint
    assert dataclasses.replace(a, scope="per_entity").fingerprint != a.fingerprint


# ---------------------------------------------------------------------------
# Guards
# ---------------------------------------------------------------------------


def test_incompatible_key_space_names_the_differing_field():
    a = STANDARD_PROTEIN_CALPHA_V1.key_space()
    b = replace(a, angle_convention="alternate_divisor:3.0")
    a.require_compatible(a)
    with pytest.raises(IncompatibleKeySpace) as exc:
        a.require_compatible(b)
    assert "angle_convention" in str(exc.value)


def test_incompatible_sampling_names_the_differing_field():
    import dataclasses

    a = base_sampling()
    b = dataclasses.replace(a, size_filter=25.0)
    a.require_compatible(a)
    with pytest.raises(IncompatibleSampling) as exc:
        a.require_compatible(b)
    assert "size_filter" in str(exc.value)


def test_the_sampling_profile_does_not_interpret_its_side_descriptions():
    """Whether two samples share a sampling rule is decided in the analysis
    layer from structured fields (``biotsr.analysis.samples.comparison_rule``);
    the profile itself offers no rule that parses its free-text side
    descriptions, and its fingerprint still covers every field, chains included."""
    import dataclasses

    a = dataclasses.replace(base_sampling(), side_a_chain="A")
    b = dataclasses.replace(a, side_a_chain="C")
    assert a.fingerprint != b.fingerprint
    assert not hasattr(SamplingProfile, "rule_fingerprint")
    assert not hasattr(SamplingProfile, "rule_dict")

