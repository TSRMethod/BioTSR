"""The native result round trip, and where its fidelity stops.

Three questions, asked separately because they have different answers:

1. do the two scientific identities survive a round trip **exactly**? (yes, and
   that is the whole reason a loaded result can be compared with anything);
2. does the whole result survive? (everything the format records -- and two
   things it deliberately does not);
3. does a reader refuse what it must refuse? (an unknown schema, a corrupt file,
   and files from two different results placed side by side).
"""

from __future__ import annotations

import copy
import json
import os
import shutil
from dataclasses import replace

import pytest

import biotsr
from biotsr.domain.errors import (
    IncompatibleRepresentation,
    MalformedResult,
    TSRError,
    UnsupportedFormat,
)
from biotsr.domain.keyspace import KeySpace, SamplingProfile, SchemeRef
from biotsr.domain.result import RESULT_SCHEMA_VERSION
from biotsr.io import native, native_read


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def save_and_load(result, tmp_path, name="out", triplets=False):
    out = str(tmp_path / name)
    written = native.save(result, out, triplets=triplets)
    return native_read.load_result(written["keys"]), written


@pytest.fixture(scope="module")
def ubq(ubq_pdb):
    """The session fixture, renamed locally so every test here reads the same."""
    return ubq_pdb


@pytest.fixture(scope="module")
def result(ubq):
    return biotsr.analyze(ubq, method="protein", chain="A", keep_triplets=False)


# ---------------------------------------------------------------------------
# KeySpace
# ---------------------------------------------------------------------------


def test_key_space_round_trips_to_an_equal_object(result):
    restored = KeySpace.from_dict(dict(result.key_space.to_dict()))
    assert restored == result.key_space


def test_key_space_round_trips_to_an_identical_fingerprint(result):
    restored = KeySpace.from_dict(dict(result.key_space.to_dict()))
    assert restored.fingerprint == result.key_space.fingerprint


def test_key_space_round_trip_survives_json(result):
    """Through the encoding the sidecar actually uses, not just through a dict."""
    payload = json.loads(json.dumps(dict(result.key_space.to_dict())))
    assert KeySpace.from_dict(payload).fingerprint == result.key_space.fingerprint


def test_a_round_tripped_key_space_is_compatible_with_its_original(result):
    KeySpace.from_dict(dict(result.key_space.to_dict())).require_compatible(
        result.key_space
    )


@pytest.mark.parametrize(
    "missing",
    [
        "angle_convention",
        "sqrt_convention",
        "median_convention",
        "numerical_precision",
        "canonical_ordering_version",
        "key_formula_version",
        "mirror_convention",
        "theta_scheme",
        "distance_scheme",
        "label_scheme",
    ],
)
def test_no_key_space_field_has_a_default(result, missing):
    """Every field is identity-bearing, so a missing one is refused, not guessed.

    A default here would build a key space describing a calculation nobody ran,
    and hand it a fingerprint that compares equal to other results built from the
    same default -- the one failure mode ``KeySpace`` exists to make impossible.
    """
    payload = dict(result.key_space.to_dict())
    del payload[missing]
    with pytest.raises(MalformedResult) as exc:
        KeySpace.from_dict(payload)
    assert missing in str(exc.value)


def test_an_unknown_key_space_field_is_refused(result):
    payload = dict(result.key_space.to_dict())
    payload["quantum_convention"] = "many_worlds"
    with pytest.raises(MalformedResult) as exc:
        KeySpace.from_dict(payload)
    assert "quantum_convention" in str(exc.value)


def test_a_scheme_ref_round_trips(result):
    ref = result.key_space.theta_scheme
    assert SchemeRef.from_dict(dict(ref.to_dict())) == ref


@pytest.mark.parametrize("bad", ["29", 29.0, None, True, 0, -1])
def test_a_scheme_radix_must_be_a_positive_integer(result, bad):
    payload = dict(result.key_space.theta_scheme.to_dict())
    payload["radix"] = bad
    with pytest.raises(MalformedResult):
        SchemeRef.from_dict(payload)


# ---------------------------------------------------------------------------
# SamplingProfile
# ---------------------------------------------------------------------------


def test_sampling_profile_round_trips_exactly(result):
    restored = SamplingProfile.from_dict(dict(result.sampling.to_dict()))
    assert restored == result.sampling
    assert restored.fingerprint == result.sampling.fingerprint


def test_sampling_profile_round_trips_with_every_optional_field_set():
    """The cross case: none of the 21 fields is None, so all 21 are exercised."""
    full = SamplingProfile(
        selector="a|b",
        scope="pooled",
        candidate_strategy="mixed_side",
        completeness_policy="none",
        composition_rule="mixed",
        size_filter=25.0,
        degenerate_policy="skip_and_record+acos_clamp",
        seq_label_policy="accept_all",
        cross_method="drug_protein_cross",
        side_a_chain="A",
        side_b_chain="B",
        side_a_selector="entity",
        side_b_selector="calpha",
        contact_cutoff=3.5,
        contact_filter="both_sides",
        contact_filter_rule="symmetric_atom_contact_filter",
        contact_projection="atom_level_contact_pool",
        pool_identity="physical",
        pool_order="first_contact",
    )
    restored = SamplingProfile.from_dict(json.loads(json.dumps(dict(full.to_dict()))))
    assert restored == full
    assert restored.fingerprint == full.fingerprint


@pytest.mark.parametrize("field_name", ["size_filter", "contact_cutoff"])
def test_an_integer_and_a_float_cutoff_give_the_same_fingerprint(result, field_name):
    """``3`` and ``3.0`` are the same science and must hash the same.

    ``fingerprint`` goes through ``json.dumps``, which spells them differently.
    Without normalisation two runs at a 3 A cutoff would be refused as
    incomparable because one caller wrote an int.
    """
    payload = dict(result.sampling.to_dict())
    payload[field_name] = 3
    as_int = SamplingProfile.from_dict(payload)
    payload[field_name] = 3.0
    as_float = SamplingProfile.from_dict(payload)
    assert as_int.fingerprint == as_float.fingerprint
    assert isinstance(getattr(as_int, field_name), float)


@pytest.mark.parametrize("bad", [0, -1.0, "3.5", True])
def test_a_non_positive_or_non_numeric_cutoff_is_refused(result, bad):
    payload = dict(result.sampling.to_dict())
    payload["contact_cutoff"] = bad
    with pytest.raises(MalformedResult):
        SamplingProfile.from_dict(payload)


def test_a_missing_sampling_field_is_refused(result):
    payload = dict(result.sampling.to_dict())
    del payload["seq_label_policy"]
    with pytest.raises(MalformedResult) as exc:
        SamplingProfile.from_dict(payload)
    assert "seq_label_policy" in str(exc.value)


def test_an_unknown_sampling_field_is_refused(result):
    payload = dict(result.sampling.to_dict())
    payload["worker_count"] = "8"
    with pytest.raises(MalformedResult):
        SamplingProfile.from_dict(payload)


# ---------------------------------------------------------------------------
# the whole result
# ---------------------------------------------------------------------------


def test_the_sidecar_declares_a_schema_version(result):
    assert result.metadata()["result_schema_version"] == RESULT_SCHEMA_VERSION


def test_every_key_and_frequency_survives(result, tmp_path):
    loaded, _ = save_and_load(result, tmp_path)
    assert dict(loaded.frequencies) == dict(result.frequencies)
    assert {k: dict(v.by_category) for k, v in loaded.keys.items()} == {
        k: dict(v.by_category) for k, v in result.keys.items()
    }


def test_both_identities_survive_with_their_fingerprints(result, tmp_path):
    loaded, _ = save_and_load(result, tmp_path)
    assert loaded.key_space == result.key_space
    assert loaded.sampling == result.sampling
    assert loaded.key_space.fingerprint == result.key_space.fingerprint
    assert loaded.sampling.fingerprint == result.sampling.fingerprint


def test_the_scalar_provenance_survives(result, tmp_path):
    loaded, _ = save_and_load(result, tmp_path)
    assert loaded.method == result.method
    assert loaded.structure_id == result.structure_id
    assert loaded.chain == result.chain
    assert loaded.counts == result.counts
    assert loaded.selection == result.selection
    assert loaded.source == result.source
    assert tuple(loaded.warnings) == tuple(result.warnings)
    assert tuple(loaded.skipped) == tuple(result.skipped)


def test_the_profile_survives(result, tmp_path):
    """A saved result must say which convention created it.

    The first sidecar format recorded no profile at all, so this was not answerable
    from a saved result by any means other than the filename.
    """
    loaded, _ = save_and_load(result, tmp_path)
    for key in ("profile", "base_profile", "legacy_compatibility_profile"):
        assert loaded.annotations[key] == result.annotations[key]
    assert tuple(loaded.annotations["modifiers"]) == tuple(
        result.annotations["modifiers"]
    )


def test_a_loaded_result_says_it_was_loaded(result, tmp_path):
    loaded, written = save_and_load(result, tmp_path)
    marker = loaded.annotations[native_read.LOADED_FROM]
    assert marker["keys"] == os.path.abspath(written["keys"])
    assert marker["result_schema_version"] == RESULT_SCHEMA_VERSION
    assert marker["triplets_loaded"] is False
    assert marker["key_order_loaded"] is False
    assert native_read.LOADED_FROM not in result.annotations


# -- the two documented fidelity boundaries ---------------------------------


def test_triplets_are_not_reconstructed(result, ubq, tmp_path):
    """ADR-075: ``alt_loc`` is per-atom and is not a column, so it would be a guess."""
    full = biotsr.analyze(ubq, method="protein", chain="A", keep_triplets=True)
    loaded, _ = save_and_load(full, tmp_path, name="with_triplets", triplets=True)
    assert loaded.triplets is None


def test_key_order_is_empty_rather_than_fabricated(result, tmp_path):
    loaded, _ = save_and_load(result, tmp_path)
    assert loaded.key_order == ()
    assert sorted(loaded.keys) != []  # the keys themselves are all there


# ---------------------------------------------------------------------------
# byte-stable re-serialization
# ---------------------------------------------------------------------------


def test_re_serializing_a_loaded_result_is_byte_identical(result, tmp_path):
    """``serialize(deserialize(serialize(r))) == serialize(r)``.

    Claimed for the keys file and the sidecar, both of which are canonical: the
    keys file is sorted and the JSON has sorted keys at every level.  It is not
    claimed for the triplet file, which is not read back at all.
    """
    once = native.save(result, str(tmp_path / "a"), triplets=False)
    loaded = native_read.load_result(once["keys"])
    twice = native.save(loaded, str(tmp_path / "b"), triplets=False)
    for role in ("keys", "metadata"):
        with open(once[role], "rb") as fh:
            first = fh.read()
        with open(twice[role], "rb") as fh:
            second = fh.read()
        assert first == second, "%s differs after a round trip" % role


def test_a_loaded_result_loads_again(result, tmp_path):
    once = native.save(result, str(tmp_path / "a"), triplets=False)
    twice = native.save(
        native_read.load_result(once["keys"]), str(tmp_path / "b"), triplets=False
    )
    again = native_read.load_result(twice["keys"])
    assert dict(again.frequencies) == dict(result.frequencies)
    assert again.key_space.fingerprint == result.key_space.fingerprint


# ---------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------


def _rewrite_sidecar(written, mutate):
    with open(written["metadata"]) as fh:
        payload = json.load(fh)
    mutate(payload)
    with open(written["metadata"], "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return written


def test_a_sidecar_with_no_schema_version_is_refused(result, tmp_path):
    """Every sidecar written before the schema was versioned is in this state,
    and is refused deliberately."""
    written = native.save(result, str(tmp_path / "old"), triplets=False)
    _rewrite_sidecar(written, lambda p: p.pop("result_schema_version"))
    with pytest.raises(UnsupportedFormat) as exc:
        native_read.load_result(written["keys"])
    assert "result_schema_version" in str(exc.value)


def test_a_future_schema_version_is_refused_and_says_to_upgrade(result, tmp_path):
    written = native.save(result, str(tmp_path / "future"), triplets=False)
    _rewrite_sidecar(
        written,
        lambda p: p.__setitem__("result_schema_version", RESULT_SCHEMA_VERSION + 7),
    )
    with pytest.raises(UnsupportedFormat) as exc:
        native_read.load_result(written["keys"])
    message = str(exc.value)
    assert "newer biotsr" in message and "Upgrade" in message


@pytest.mark.parametrize("bad", ["1", 1.0, None, [1]])
def test_a_non_integer_schema_version_is_malformed(result, tmp_path, bad):
    written = native.save(result, str(tmp_path / ("bad%s" % type(bad).__name__)),
                          triplets=False)
    _rewrite_sidecar(written, lambda p: p.__setitem__("result_schema_version", bad))
    with pytest.raises((MalformedResult, UnsupportedFormat)):
        native_read.load_result(written["keys"])


def test_a_keys_file_with_no_sidecar_is_refused(result, tmp_path):
    written = native.save(result, str(tmp_path / "lonely"), triplets=False)
    os.remove(written["metadata"])
    with pytest.raises(FileNotFoundError) as exc:
        native_read.load_result(written["keys"])
    assert "key space" in str(exc.value)


def test_a_truncated_keys_file_is_refused(result, tmp_path):
    """The sidecar's counts no longer match what the file contains."""
    written = native.save(result, str(tmp_path / "cut"), triplets=False)
    with open(written["keys"]) as fh:
        head = fh.readlines()[:50]
    with open(written["keys"], "w") as fh:
        fh.writelines(head)
    with pytest.raises(MalformedResult) as exc:
        native_read.load_result(written["keys"])
    assert "n_unique_keys" in str(exc.value)


def test_a_keys_file_with_a_wrong_header_is_refused(result, tmp_path):
    written = native.save(result, str(tmp_path / "hdr"), triplets=False)
    with open(written["keys"]) as fh:
        rest = fh.readlines()[1:]
    with open(written["keys"], "w") as fh:
        fh.write("key\tfrequency\n")
        fh.writelines(rest)
    with pytest.raises(MalformedResult) as exc:
        native_read.load_result(written["keys"])
    assert "header" in str(exc.value)


def test_a_keys_file_out_of_order_is_refused(result, tmp_path):
    """The format promises sorted order; a file in another order is not ours."""
    written = native.save(result, str(tmp_path / "shuffled"), triplets=False)
    with open(written["keys"]) as fh:
        lines = fh.readlines()
    lines[1], lines[2] = lines[2], lines[1]
    with open(written["keys"], "w") as fh:
        fh.writelines(lines)
    with pytest.raises(MalformedResult) as exc:
        native_read.load_result(written["keys"])
    assert "ordering guarantee" in str(exc.value)


def test_an_edited_key_space_is_caught_by_its_own_fingerprint(result, tmp_path):
    """The sidecar must not be able to claim an identity it does not derive."""
    written = native.save(result, str(tmp_path / "edited"), triplets=False)

    def mutate(payload):
        payload["key_space"]["mirror_convention"] = "negate_when_prefold_gt_90"

    _rewrite_sidecar(written, mutate)
    with pytest.raises(MalformedResult) as exc:
        native_read.load_result(written["keys"])
    assert "contradicts itself" in str(exc.value)


def test_two_different_results_cannot_be_combined(ubq, tmp_path):
    """``sample_A.keys`` + ``sample_B.meta.json`` must fail loudly."""
    a = biotsr.analyze(ubq, method="protein", chain="A", keep_triplets=False)
    b = biotsr.analyze(
        ubq, method="protein", chain="A", size_filter=12.0, keep_triplets=False
    )
    wa = native.save(a, str(tmp_path / "a"), triplets=False)
    wb = native.save(b, str(tmp_path / "b"), triplets=False)
    assert a.counts.n_triangles != b.counts.n_triangles
    with pytest.raises(MalformedResult) as exc:
        native_read.load_result(wa["keys"], metadata_path=wb["metadata"])
    assert "do not belong to the same result" in str(exc.value)


@pytest.mark.parametrize(
    "name,mutate",
    [
        ("counts_missing_a_field", lambda p: p["counts"].pop("n_points")),
        ("illegal_structure_id", lambda p: p.__setitem__("structure_id", "a b")),
        ("non_string_structure_id", lambda p: p.__setitem__("structure_id", 7)),
        ("source_missing_a_field", lambda p: p["source"].pop("parser_backend")),
        ("selection_missing_a_field", lambda p: p["selection"].pop("n_points")),
        ("malformed_skip_record", lambda p: p.__setitem__("skipped", [{"stage": "x"}])),
        ("counts_not_an_object", lambda p: p.__setitem__("counts", [1, 2])),
        (
            "category_totals_not_an_object",
            lambda p: p.__setitem__("category_totals", ["regular", 1]),
        ),
    ],
)
def test_a_corrupt_sidecar_always_raises_malformed_result(
    result, tmp_path, name, mutate
):
    """The contract is a type, not just a message.

    Without translation at the reader's boundary these arrive as ``KeyError``,
    ``TypeError`` or ``ValueError`` -- correct diagnoses of a corrupt file under
    the wrong name. A caller writing ``except MalformedResult`` would miss them.
    """
    written = native.save(result, str(tmp_path / name), triplets=False)
    _rewrite_sidecar(written, mutate)
    with pytest.raises(MalformedResult):
        native_read.load_result(written["keys"])


def test_a_file_of_the_candidate_layout_is_refused_by_name(result, tmp_path):
    """``1UBQ_A.keys_29_35`` (1.0.0rc3 and earlier) is never read as current."""
    from biotsr.domain.errors import UnsupportedFormat

    written = native.save(result, str(tmp_path / "old"), triplets=False)
    old = os.path.join(os.path.dirname(written["keys"]), "1UBQ_A.keys_29_35")
    os.rename(written["keys"], old)
    os.rename(written["metadata"], old + ".meta.json")
    with pytest.raises(UnsupportedFormat, match="release candidate"):
        native_read.load_result(old)


def test_a_schema_1_sidecar_is_refused_as_a_release_candidate_result(result, tmp_path):
    from biotsr.domain.errors import UnsupportedFormat

    written = native.save(result, str(tmp_path / "v1"), triplets=False)
    _rewrite_sidecar(written, lambda p: p.update(result_schema_version=1))
    with pytest.raises(UnsupportedFormat, match="1.0.0rc3 or earlier"):
        native_read.load_result(written["keys"])


# ---------------------------------------------------------------------------
# the error hierarchy
# ---------------------------------------------------------------------------


def test_every_package_refusal_is_a_tsr_error():
    for kind in (
        biotsr.IncompatibleKeySpace,
        biotsr.IncompatibleSampling,
        biotsr.IncompatibleModifier,
        biotsr.IncompatibleSSEStratification,
        biotsr.MalformedResult,
        biotsr.UnsupportedFormat,
        biotsr.InvalidRequest,
    ):
        assert issubclass(kind, TSRError)


def test_every_package_refusal_is_still_a_value_error():
    """Existing code that writes ``except ValueError`` must keep working (ADR-076)."""
    for kind in (
        biotsr.IncompatibleKeySpace,
        biotsr.IncompatibleSampling,
        biotsr.IncompatibleModifier,
        biotsr.IncompatibleSSEStratification,
        biotsr.MalformedResult,
        biotsr.UnsupportedFormat,
        biotsr.InvalidRequest,
        TSRError,
    ):
        assert issubclass(kind, ValueError)


def test_the_four_incompatibility_types_share_a_parent():
    for kind in (
        biotsr.IncompatibleKeySpace,
        biotsr.IncompatibleSampling,
        biotsr.IncompatibleModifier,
        biotsr.IncompatibleSSEStratification,
    ):
        assert issubclass(kind, IncompatibleRepresentation)


def test_a_programming_error_is_not_a_tsr_error():
    """The exclusion is what gives ``except TSRError`` its meaning."""
    for kind in (AttributeError, TypeError, AssertionError, KeyError):
        assert not issubclass(kind, TSRError)


# ---------------------------------------------------------------------------
# mutation safety
# ---------------------------------------------------------------------------


def test_metadata_does_not_share_its_sse_block_with_the_result(ubq):
    """Mutating the metadata payload must not reach the result (round-trip path)."""
    annotated = biotsr.analyze(
        ubq, method="protein", chain="A", sse=True, keep_triplets=False
    )
    before = copy.deepcopy(annotated.annotations["sse"])
    payload = annotated.metadata()
    payload["sse"]["strand_identity"] = "tampered"
    assert annotated.annotations["sse"] == before


def test_metadata_does_not_share_its_provenance_with_the_result(result):
    payload = result.metadata()
    payload["provenance"]["profile"] = "tampered"
    assert result.annotations["profile"] != "tampered"


def test_two_results_do_not_share_a_category_counter(ubq):
    a = biotsr.analyze(ubq, method="protein", chain="A", keep_triplets=False)
    b = biotsr.analyze(ubq, method="protein", chain="A", keep_triplets=False)
    key = sorted(a.keys)[0]
    assert a.keys[key] is not b.keys[key]
    assert a.keys[key].by_category is not b.keys[key].by_category


def test_a_loaded_result_shares_nothing_with_its_sidecar(result, tmp_path):
    loaded, written = save_and_load(result, tmp_path)
    with open(written["metadata"]) as fh:
        payload = json.load(fh)
    payload["sampling"]["size_filter"] = 1.0
    assert loaded.sampling.size_filter is None


# ---------------------------------------------------------------------------
# configuration validation reaches the Python API, not only the CLI
# ---------------------------------------------------------------------------


def test_a_profile_from_another_method_is_refused(ubq):
    """The CLI could not express this; the Python API could, and must refuse it too."""
    with pytest.raises(biotsr.InvalidRequest) as exc:
        biotsr.analyze(
            ubq, method="drug", chain="A", entity_name="BEN",
            profile=biotsr.STANDARD_PROTEIN_CALPHA_V1,
        )
    message = str(exc.value)
    assert "standard_protein_calpha_v1" in message
    assert "is not the standard profile of method 'drug'" in message
    assert "Omit profile=" in message


def test_a_derived_profile_is_refused_even_in_the_right_domain(ubq):
    """Every method runs under its own standard profile; a derived one is refused."""
    derived = replace(biotsr.STANDARD_PROTEIN_CALPHA_V1, name="derived_for_a_study")
    with pytest.raises(biotsr.InvalidRequest) as exc:
        biotsr.analyze(ubq, method="protein", chain="A", profile=derived,
                       keep_triplets=False)
    message = str(exc.value)
    assert "derived_for_a_study" in message
    assert "is not the standard profile of method 'protein'" in message


def test_a_chain_the_structure_does_not_contain_is_refused(ubq_pdb):
    """An unknown chain is refused. It used to return an empty result with a
    fingerprint, no warning and no skip record."""
    with pytest.raises(biotsr.InvalidRequest) as exc:
        biotsr.analyze(ubq_pdb, method="protein", chain="ZZ")
    message = str(exc.value)
    assert "'ZZ'" in message and "1UBQ" in message
    # The command-line contract's two phrases are preserved: moving
    # WHERE a request is refused must not change WHAT the user is told.
    assert "no atoms" in message and "chains in this structure" in message


def test_a_chain_that_exists_but_selects_nothing_is_not_an_error(bdna_path):
    """The distinction the check is built around.

    ``protein`` over 1BNA's DNA chain A selects no Ca, and that is a *fact about
    the structure* -- `biotsr/cli/spec.py::_require_points` says so, and says why a
    library must not call it an error. Only a chain that is not there is refused.
    """
    result = biotsr.analyze(bdna_path, method="protein", chain="A",
                         keep_triplets=False)
    assert result.counts.n_points == 0
    assert result.n_triangles == 0


@pytest.mark.parametrize(
    "call",
    [
        lambda p: biotsr.analyze(p, method="protein", chain="ZZ"),
        lambda p: biotsr.analyze_entities(p, method="amino_acid_atomic", chain="ZZ"),
        lambda p: biotsr.list_residue_targets(p, chain="ZZ"),
        lambda p: biotsr.list_entity_targets(p, chain="ZZ", method="drug"),
    ],
)
def test_every_entry_point_refuses_an_absent_chain(ubq_pdb, call):
    with pytest.raises(biotsr.InvalidRequest):
        call(ubq_pdb)


def test_a_cross_run_refuses_an_absent_chain_on_either_side(barnase_path):
    for kwargs in (
        dict(chain="ZZ", chain_2="D"),
        dict(chain="A", chain_2="ZZ"),
    ):
        with pytest.raises(biotsr.InvalidRequest):
            biotsr.analyze_cross(
                barnase_path, method="protein_protein_cross",
                contact_cutoff=3.2, keep_triplets=False, **kwargs
            )


def test_the_methods_own_default_is_always_accepted():
    from biotsr.api import _require_standard_profile
    from biotsr.methods.registry import get_method

    for name in biotsr.list_methods():
        recipe = get_method(name)
        _require_standard_profile(recipe, None)
        _require_standard_profile(recipe, recipe.default_profile)
        # An equal (unchanged) copy is the same profile.
        _require_standard_profile(recipe, replace(recipe.default_profile))


# ---------------------------------------------------------------------------
# the analysis row label is a label; the persisted identity is an identity
# ---------------------------------------------------------------------------
#
# `biotsr analyze --from-results` names each row after the sample DIRECTORY, so that
# two saved runs over the same chain -- which share a keys filename -- can be
# told apart in a report.  That is a presentation decision about one report, and
# it must not reach the result: a `TSRResult` says which calculation it is, and
# moving or renaming a directory on disk changes no calculation.
#
# Nothing enforces the separation structurally -- `AnalysisSample.from_result`
# takes a `sample_id` and could as easily have been written to stamp it back --
# so it is pinned here.


def _persisted_identity(result):
    """Everything a saved result asserts about *which calculation it is*.

    Deliberately not a fingerprint comparison alone: a fingerprint would still
    match if `structure_id`, `chain`, `method` or the profile were rewritten,
    because none of those enter either identity.  They are exactly the fields a
    path-derived label could plausibly have overwritten, so they are compared
    field by field.
    """
    return {
        "structure_id": result.structure_id,
        "chain": result.chain,
        "method": result.method,
        "key_space": result.key_space,
        "key_space_fingerprint": result.key_space.fingerprint,
        "sampling": result.sampling,
        "sampling_fingerprint": result.sampling.fingerprint,
        "counts": result.counts,
        "selection": result.selection,
        "source": result.source,
        "frequencies": dict(result.frequencies),
        "categories": {k: dict(v.by_category) for k, v in result.keys.items()},
        "profile": result.annotations["profile"],
        "base_profile": result.annotations["base_profile"],
        "legacy_compatibility": result.annotations["legacy_compatibility_profile"],
        "modifiers": tuple(result.annotations["modifiers"]),
        "target": result.annotations["target"],
        "sse": result.annotations["sse"],
    }


def _move_directory(keys_path, tmp_path, new_name):
    """Move the containing directory somewhere else, under a different name.

    Returns the keys file's new path.  The keys FILE keeps its name, because the
    format's filename is part of the contract and renaming it is refused
    (`test_a_renamed_keys_file_is_caught_before_the_counts`).  The directory is
    not part of the contract, and this is what proves it.
    """
    old_dir = os.path.dirname(keys_path)
    new_dir = str(tmp_path / "relocated" / new_name)
    os.makedirs(os.path.dirname(new_dir), exist_ok=True)
    shutil.move(old_dir, new_dir)
    return os.path.join(new_dir, os.path.basename(keys_path))


def test_moving_and_renaming_a_result_directory_changes_no_persisted_identity(
    result, tmp_path
):
    """The load-bearing assertion: a directory is not an identity."""
    written = native.save(result, str(tmp_path / "original_name"), triplets=False)
    before = _persisted_identity(native_read.load_result(written["keys"]))

    moved = _move_directory(written["keys"], tmp_path, "a totally different name")
    after = _persisted_identity(native_read.load_result(moved))

    assert after == before
    # Spelled out for the fields a path-derived label would have overwritten, so
    # a failure names which one moved rather than printing two large dicts.
    assert after["structure_id"].value == "1UBQ"
    assert after["chain"] == "A"
    assert after["method"] == "protein"
    assert after["profile"] == "standard_protein_calpha_v1"


def test_a_relocated_result_keeps_its_source_provenance(result, tmp_path):
    """Where the *structure* came from is unrelated to where the result now sits."""
    written = native.save(result, str(tmp_path / "before"), triplets=False)
    original = native_read.load_result(written["keys"])
    relocated = native_read.load_result(
        _move_directory(written["keys"], tmp_path, "after")
    )
    assert relocated.source == original.source
    assert relocated.source.sha256 == result.source.sha256
    assert relocated.source.path == result.source.path


def test_only_loaded_from_tracks_the_new_location(result, tmp_path):
    """The one field that *should* change, and the fact that it is the only one.

    `loaded_from` is provenance about the read, not about the science, which is
    why it is kept out of `_persisted_identity` and out of `metadata()`.
    """
    written = native.save(result, str(tmp_path / "here"), triplets=False)
    before = native_read.load_result(written["keys"])
    moved_path = _move_directory(written["keys"], tmp_path, "there")
    after = native_read.load_result(moved_path)

    assert before.annotations[native_read.LOADED_FROM]["keys"] != (
        after.annotations[native_read.LOADED_FROM]["keys"]
    )
    assert after.annotations[native_read.LOADED_FROM]["keys"] == (
        os.path.abspath(moved_path)
    )
    # And it does not leak into what gets re-serialized: a relocated result
    # writes the same sidecar bytes as the original.
    a = native.save(before, str(tmp_path / "reser_a"), triplets=False)
    b = native.save(after, str(tmp_path / "reser_b"), triplets=False)
    for role in ("keys", "metadata"):
        with open(a[role], "rb") as fh:
            first = fh.read()
        with open(b[role], "rb") as fh:
            second = fh.read()
        assert first == second, "%s differs after relocation" % role


def test_the_prefix_names_the_analysis_row_and_nothing_else(ubq_pdb, tmp_path):
    """The completion record's prefix labels the row; it is never written back."""
    from biotsr.analysis import AnalysisSample
    from biotsr.cli import EXIT_OK, main
    from biotsr.cli.commands import _saved_results
    from biotsr.cli.console import QUIET, Console

    root = str(tmp_path / "root")
    assert main(["run", "-m", "protein", "--pdb", ubq_pdb, "-c", "A", "-o", root,
                 "--sample-id", "run_one", "--no-triplets", "-q"]) == EXIT_OK
    [(label, loaded, _, keys)] = _saved_results(root, Console(QUIET))
    assert label == "run_one" and os.path.basename(keys) == "run_one.keys"

    sample = AnalysisSample.from_result(loaded, sample_id=label)
    assert sample.sample_id == "run_one"          # the row is labelled
    assert loaded.structure_id.value == "1UBQ"    # the result is not renamed
    assert loaded.chain == "A"
    # With no override the id comes from the RESULT's own identity, not the path.
    assert AnalysisSample.from_result(loaded).sample_id == "1UBQ_A"


def test_the_cli_analysis_is_unchanged_by_relocating_the_whole_folder(
    ubq_pdb, trypsin_path, tmp_path
):
    """End to end: the same science, the same numbers, the same row labels."""
    from biotsr.cli import EXIT_OK, main

    original = str(tmp_path / "runs")
    for path in (ubq_pdb, trypsin_path):
        assert main(["run", "-m", "protein", "--pdb", path, "-c", "A",
                     "-o", original, "--no-triplets", "-q"]) == EXIT_OK
    assert main(["analyze", "common-keys", "--from-results", original,
                 "-o", str(tmp_path / "an_before"), "-q"]) == EXIT_OK

    relocated = str(tmp_path / "moved elsewhere" / "runs")
    shutil.move(original, relocated)
    assert main(["analyze", "common-keys", "--from-results", relocated,
                 "-o", str(tmp_path / "an_after"), "-q"]) == EXIT_OK

    def report(where):
        path = tmp_path / where / "analysis" / "common_keys" / "all_samples.json"
        with open(str(path)) as fh:
            return json.load(fh)

    before, after = report("an_before"), report("an_after")
    assert after["keys"] == before["keys"]
    assert after["rows"] == before["rows"]
    assert [r["fileName"] for r in after["rows"]] == ["1UBQ_A_protein", "3PTB_A_protein"]

    # The results on disk still assert their own identity, unchanged.
    for prefix in ("1UBQ_A_protein", "3PTB_A_protein"):
        loaded = native_read.load_result(os.path.join(relocated, prefix + ".keys"))
        assert loaded.method == "protein"
        assert loaded.chain == "A"
        assert loaded.annotations["profile"] == "standard_protein_calpha_v1"


def test_a_relocated_result_still_refuses_a_foreign_sidecar(ubq_pdb, tmp_path):
    """Relocation must not weaken the multi-file consistency checks."""
    a = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    b = biotsr.analyze(ubq_pdb, method="protein", chain="A", size_filter=12.0,
                    keep_triplets=False)
    wa = native.save(a, str(tmp_path / "a"), triplets=False)
    wb = native.save(b, str(tmp_path / "b"), triplets=False)
    moved = _move_directory(wa["keys"], tmp_path, "relocated a")
    with pytest.raises(MalformedResult):
        native_read.load_result(moved, metadata_path=wb["metadata"])


def test_building_an_analysis_sample_never_mutates_the_result(result, tmp_path):
    """The structural reason a row label cannot reach an identity.

    `AnalysisSample.from_result` takes a `sample_id` override, and the CLI uses
    it for every `--from-results` row.  If that override were stamped back onto
    the result -- or if the sample aliased the result's mappings -- the label
    would have become part of the identity.  It does neither.
    """
    from biotsr.analysis import AnalysisSample

    written = native.save(result, str(tmp_path / "mutation"), triplets=False)
    loaded = native_read.load_result(written["keys"])
    before = _persisted_identity(loaded)

    sample = AnalysisSample.from_result(loaded, sample_id="a label from a path")
    assert sample.sample_id == "a label from a path"
    assert _persisted_identity(loaded) == before

    # Nor does mutating the sample's own counts reach back into the result.
    sample.counts[999999] = 1234
    assert 999999 not in loaded.frequencies
    assert _persisted_identity(loaded) == before


def test_the_cli_loader_yields_samples_matching_an_independent_reload(
    ubq_pdb, trypsin_path, tmp_path
):
    """`_samples_from_saved` relabels; it must not alter what it loaded.

    Each sample the CLI builds is compared against a *separate, direct* load of
    the same file, so anything the CLI path did to the science on the way through
    shows up as a difference against the file itself rather than against another
    copy of the CLI's own answer.
    """
    from biotsr.analysis import AnalysisSample
    from biotsr.cli import EXIT_OK, main
    from biotsr.cli.commands import _samples_from_saved
    from biotsr.cli.console import QUIET, Console

    root = str(tmp_path / "runs")
    for path in (ubq_pdb, trypsin_path):
        assert main(["run", "-m", "protein", "--pdb", path, "-c", "A",
                     "-o", root, "--no-triplets", "-q"]) == EXIT_OK

    class _Args(object):
        from_results = root

    produced = _samples_from_saved(_Args(), Console(QUIET))
    assert {s.sample_id for s in produced} == {"1UBQ_A_protein", "3PTB_A_protein"}

    for sample in produced:
        fresh = native_read.load_result(os.path.join(root, sample.sample_id + ".keys"))
        reference = AnalysisSample.from_result(fresh)

        # The science matches the file, independently read ...
        assert dict(sample.counts) == dict(reference.counts)
        assert sample.key_space == reference.key_space
        assert sample.sampling == reference.sampling
        assert sample.key_space.fingerprint == reference.key_space.fingerprint
        assert sample.sampling.fingerprint == reference.sampling.fingerprint
        assert sample.annotations["method"] == fresh.method == "protein"
        # ... while the id the file itself implies is the structural one, and the
        # directory label the CLI chose is only ever on the sample.
        assert reference.sample_id == "%s_%s" % (fresh.structure_id.value,
                                                 fresh.chain)
        assert sample.sample_id != reference.sample_id
        assert fresh.structure_id.value in ("1UBQ", "3PTB")
