#!/usr/bin/env python3
"""API-hardening gates: the hardening promises, each checked as a promise.

API hardening adds no science.  What it adds is a set of *guarantees* about the
package as an artefact, and a gate here exists for each one -- not one per test,
and not one per function.  The categories are the five things API hardening promised:

``round_trip``
    A saved result can be read back, and what comes back is the same science.
    Every representative family, both identities, byte-stable re-serialization.
``refusal``
    The reader refuses what it must refuse: an unversioned or future schema, a
    corrupt file, a file pair from two different results, an edited sidecar.
``analysis``
    Every established analysis operates on a loaded result and produces exactly
    what it produces on the computed one.  No new metric; the same ones.
``surface``
    The public API is what it says it is: documented imports resolve, the error
    hierarchy is coherent, validation is in the package rather than in the CLI,
    and scientific defaults have one source of truth.
``robustness``
    Paths with spaces and non-ASCII characters, read-only inputs, output outside
    the repository, no dependency on the checkout layout, and no dependency
    beyond the one the package declares.

Usage::

    python3 validation/scripts/gates_api_hardening.py [--json PATH]
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import copy
import io
import json
import os
import shutil
import stat
import sys
import tempfile
from dataclasses import fields, replace

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
FIXTURES = os.path.join(WORKSPACE, "tests", "data", "structures")
EVIDENCE = os.path.join(WORKSPACE, "validation", "evidence", "api_hardening")

sys.path.insert(0, SRC)

import biotsr  # noqa: E402
from biotsr.analysis import (  # noqa: E402
    AnalysisSample,
    common_key_intersection,
    distance_matrix,
    feature_matrix,
    generalized_jaccard,
    similarity_matrix,
)
from biotsr.cli import EXIT_OK, EXIT_USAGE, main as cli_main  # noqa: E402
from biotsr.domain.errors import (  # noqa: E402
    IncompatibleRepresentation,
    InvalidRequest,
    MalformedResult,
    TSRError,
    UnsupportedFormat,
)
from biotsr.domain.keyspace import KeySpace, SamplingProfile  # noqa: E402
from biotsr.domain.result import RESULT_SCHEMA_VERSION  # noqa: E402
from biotsr.io import native, native_read  # noqa: E402
from biotsr.methods.registry import get_method  # noqa: E402

RESULTS = []


def check(category, name, detail=""):
    """Decorator: run the body, record pass/fail, never let one gate stop the rest."""

    def decorate(fn):
        try:
            note = fn()
            RESULTS.append(
                {"category": category, "name": name, "passed": True,
                 "detail": note or detail}
            )
        except AssertionError as exc:
            RESULTS.append(
                {"category": category, "name": name, "passed": False,
                 "detail": str(exc)[:300]}
            )
        except Exception as exc:  # noqa: BLE001
            RESULTS.append(
                {"category": category, "name": name, "passed": False,
                 "detail": "%s: %s" % (type(exc).__name__, str(exc)[:260])}
            )
        return fn

    return decorate


def cli(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli_main(argv)
    return code, out.getvalue(), err.getvalue()


def fixture(name):
    path = os.path.join(FIXTURES, "%s.pdb" % name)
    if os.path.exists(path):
        return path
    raise AssertionError("pinned fixture %s.pdb is missing" % name)


TMP = tempfile.mkdtemp(prefix="tsr-api-hardening-gates-")


def workdir(name):
    path = os.path.join(TMP, name)
    os.makedirs(path, exist_ok=True)
    return path


def save_load(result, name, triplets=False):
    written = native.save(result, workdir(name), triplets=triplets)
    return native_read.load_result(written["keys"]), written


def read_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


# ===========================================================================
# the representative families, computed once
# ===========================================================================

#: The representative families.  Each one differs from the others in a way a
#: serializer can get wrong: an entity in the filename, a pipe in the chain, all
#: 21 sampling fields populated, eighteen category buckets instead of one, a
#: different label scheme.
FAMILY_CASES = [
    ("protein_calpha", False,
     dict(source="1UBQ", method="protein", chain="A")),
    ("protein_calpha_grouped", False,
     dict(source="1UBQ", method="protein", chain="A", aa_grouping=True)),
    ("protein_calpha_mirror", False,
     dict(source="1UBQ", method="protein", chain="A", mirror=True)),
    ("protein_calpha_size_filtered", False,
     dict(source="1UBQ", method="protein", chain="A", size_filter=14.0)),
    ("protein_calpha_sse", False,
     dict(source="1TUP", method="protein", chain="A", sse=True)),
    ("amino_acid_atomic", False,
     dict(source="1UBQ", method="amino_acid_atomic", chain="A",
          residue="15", residue_name="LEU")),
    ("drug", False,
     dict(source="3PTB", method="drug", chain="A", entity_name="BEN",
          entity_id="1")),
    ("carbohydrate", False,
     dict(source="5WT9", method="carbohydrate", chain="A", entity_name="NAG",
          entity_id="1")),
    ("nucleotide", False,
     dict(source="1BNA", method="nucleotide", chain="A", entity_name="DC",
          entity_id="3")),
    ("drug_protein_cross", True,
     dict(source="3PTB", method="drug_protein_cross", chain="A", chain_2="A",
          entity_name="BEN", entity_id="1")),
    ("protein_protein_cross_all_atom", True,
     dict(source="1BRS", method="protein_protein_cross", chain="A", chain_2="D",
          contact_cutoff=3.2)),
    ("protein_protein_cross_calpha_interface", True,
     dict(source="1BRS", method="protein_protein_calpha_cross", chain="A",
          chain_2="D", contact_cutoff=3.2)),
]


def _compute(spec, is_cross):
    spec = dict(spec)
    source = fixture(spec.pop("source"))
    call = biotsr.analyze_cross if is_cross else biotsr.analyze
    return call(source, keep_triplets=False, **spec)


COMPUTED = {}
for _name, _is_cross, _spec in FAMILY_CASES:
    COMPUTED[_name] = _compute(_spec, _is_cross)


# ===========================================================================
# round_trip -- a saved result is the same science when it comes back
# ===========================================================================

for _case_name in list(COMPUTED):

    def _make(case_name):
        @check("round_trip", "keys_and_categories_survive__%s" % case_name)
        def _():
            original = COMPUTED[case_name]
            loaded, _w = save_load(original, "rt_%s" % case_name)
            assert {k: dict(v.by_category) for k, v in loaded.keys.items()} == {
                k: dict(v.by_category) for k, v in original.keys.items()
            }, "the key/category table changed"
            return "%d keys" % len(loaded.keys)

        @check("round_trip", "both_fingerprints_survive__%s" % case_name)
        def _():
            original = COMPUTED[case_name]
            loaded, _w = save_load(original, "fp_%s" % case_name)
            assert loaded.key_space.fingerprint == original.key_space.fingerprint, (
                "key space fingerprint changed"
            )
            assert loaded.sampling.fingerprint == original.sampling.fingerprint, (
                "sampling fingerprint changed"
            )
            return original.key_space.fingerprint[:16]

        @check("round_trip", "re_serialization_is_byte_identical__%s" % case_name)
        def _():
            original = COMPUTED[case_name]
            once = native.save(original, workdir("b1_%s" % case_name), triplets=False)
            loaded = native_read.load_result(once["keys"])
            twice = native.save(loaded, workdir("b2_%s" % case_name), triplets=False)
            for role in ("keys", "metadata"):
                assert read_bytes(once[role]) == read_bytes(twice[role]), (
                    "%s differs after a round trip" % role
                )
            return "keys + sidecar identical"

    _make(_case_name)


@check("round_trip", "the_sidecar_declares_an_explicit_schema_version")
def _():
    payload = COMPUTED["protein_calpha"].metadata()
    assert payload["result_schema_version"] == RESULT_SCHEMA_VERSION
    assert isinstance(payload["result_schema_version"], int)
    return "result_schema_version = %d" % RESULT_SCHEMA_VERSION


@check("round_trip", "the_schema_version_is_not_the_package_version")
def _():
    """Tying them would force a format migration on every release (ADR-074)."""
    assert str(RESULT_SCHEMA_VERSION) != biotsr.__version__
    return "%d vs %s" % (RESULT_SCHEMA_VERSION, biotsr.__version__)


@check("round_trip", "key_space_round_trips_through_json_exactly")
def _():
    ks = COMPUTED["protein_calpha"].key_space
    restored = KeySpace.from_dict(json.loads(json.dumps(dict(ks.to_dict()))))
    assert restored == ks and restored.fingerprint == ks.fingerprint
    return "%d fields" % len(fields(KeySpace))


@check("round_trip", "sampling_profile_round_trips_through_json_exactly")
def _():
    sp = COMPUTED["drug_protein_cross"].sampling
    restored = SamplingProfile.from_dict(json.loads(json.dumps(dict(sp.to_dict()))))
    assert restored == sp and restored.fingerprint == sp.fingerprint
    populated = sum(1 for f in fields(SamplingProfile)
                    if getattr(sp, f.name) is not None)
    return "%d of %d fields populated" % (populated, len(fields(SamplingProfile)))


@check("round_trip", "an_integer_and_a_float_cutoff_fingerprint_identically")
def _():
    """``3`` and ``3.0`` are the same science; ``json.dumps`` spells them apart."""
    base = dict(COMPUTED["protein_calpha"].sampling.to_dict())
    base["contact_cutoff"] = 3
    as_int = SamplingProfile.from_dict(base)
    base["contact_cutoff"] = 3.0
    as_float = SamplingProfile.from_dict(base)
    assert as_int.fingerprint == as_float.fingerprint
    assert isinstance(as_int.contact_cutoff, float)
    return "both -> %r" % as_int.contact_cutoff


@check("round_trip", "the_profile_that_produced_a_result_survives")
def _():
    """Not inferred from a filename: recorded in the sidecar."""
    for name in ("protein_calpha", "protein_calpha_grouped", "drug_protein_cross"):
        original = COMPUTED[name]
        loaded, _w = save_load(original, "prov_%s" % name)
        assert loaded.annotations["profile"] == original.annotations["profile"], name
        assert (
            loaded.annotations["legacy_compatibility_profile"]
            == original.annotations["legacy_compatibility_profile"]
        ), name
    return ", ".join(sorted({COMPUTED[n].annotations["profile"] for n in (
        "protein_calpha", "protein_calpha_grouped", "drug_protein_cross")}))


@check("round_trip", "sse_annotation_provenance_survives_in_full")
def _():
    from biotsr.sse import category_totals_in_order, check_conservation, strand_identity_of

    original = COMPUTED["protein_calpha_sse"]
    loaded, _w = save_load(original, "sse_prov")
    assert loaded.annotations["sse"] == original.annotations["sse"], "assignment map"
    assert loaded.annotations["sse_coverage"] == original.annotations["sse_coverage"]
    assert strand_identity_of(loaded) == strand_identity_of(original)
    assert category_totals_in_order(loaded) == category_totals_in_order(original)
    assert loaded.annotations["sse"]["has_records"] is True
    check_conservation(loaded)
    return "%s, %d categories" % (
        strand_identity_of(loaded),
        sum(1 for _c, n in category_totals_in_order(loaded) if n),
    )


@check("round_trip", "an_entity_scoped_result_keeps_its_target_in_the_prefix")
def _():
    """The default prefix names the molecule or residue and its identifier
    (ADR-098), and is derived from the result's own provenance, so a loaded
    result names itself exactly as the computed one did."""
    expected = {
        "amino_acid_atomic": "1UBQ_A_amino_acid_atomic_LEU_15",
        "drug": "3PTB_A_drug_BEN_1",
        "carbohydrate": "5WT9_A_carbohydrate_NAG_1",
        "nucleotide": "1BNA_A_nucleotide_DC_3",
        "drug_protein_cross": "3PTB_A_A_drug_protein_cross_BEN_1",
    }
    for case, prefix in sorted(expected.items()):
        original = COMPUTED[case]
        loaded, written = save_load(original, "entity_%s" % case)
        assert os.path.basename(written["keys"]) == prefix + native.KEYS_SUFFIX, (
            case, os.path.basename(written["keys"]))
        assert os.path.basename(written["metadata"]) == prefix + native.METADATA_SUFFIX
        assert native.default_prefix(loaded) == native.default_prefix(original), case
    return ", ".join(expected[c] for c in sorted(expected))


@check("round_trip", "a_cross_result_keeps_both_sides_and_its_pool")
def _():
    original = COMPUTED["protein_protein_cross_all_atom"]
    loaded, _w = save_load(original, "cross_sides")
    assert loaded.chain == "A|D"
    assert loaded.sampling.side_a_chain == "A" and loaded.sampling.side_b_chain == "D"
    assert loaded.sampling.contact_cutoff == 3.2
    assert loaded.annotations["cross_pools"] == original.annotations["cross_pools"]
    return "A|D at 3.2 A"


@check("round_trip", "a_loaded_result_can_be_saved_and_loaded_again")
def _():
    original = COMPUTED["protein_calpha"]
    once = native.save(original, workdir("twice_a"), triplets=False)
    twice = native.save(
        native_read.load_result(once["keys"]), workdir("twice_b"), triplets=False
    )
    again = native_read.load_result(twice["keys"])
    assert dict(again.frequencies) == dict(original.frequencies)
    assert again.key_space.fingerprint == original.key_space.fingerprint
    return "stable under repetition"


@check("round_trip", "triplets_are_declared_absent_rather_than_guessed")
def _():
    """ADR-075: ``alt_loc`` is per-atom and is in no column, so it would be a guess."""
    full = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                       keep_triplets=True)
    loaded, _w = save_load(full, "with_triplets", triplets=True)
    assert loaded.triplets is None
    assert loaded.annotations[native_read.LOADED_FROM]["triplets_loaded"] is False
    return "stated on the result, not silently None"


@check("round_trip", "key_order_is_empty_rather_than_fabricated")
def _():
    loaded, _w = save_load(COMPUTED["protein_calpha"], "ko")
    assert loaded.key_order == ()
    assert loaded.annotations[native_read.LOADED_FROM]["key_order_loaded"] is False
    return "()"


# ===========================================================================
# refusal -- the reader never guesses
# ===========================================================================


def _sidecar_variant(name, mutate):
    """A saved result whose sidecar has been edited, in its own directory."""
    written = native.save(COMPUTED["protein_calpha"], workdir(name), triplets=False)
    with open(written["metadata"]) as fh:
        payload = json.load(fh)
    mutate(payload)
    with open(written["metadata"], "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return written


@check("refusal", "an_unversioned_sidecar_is_refused")
def _():
    """Every sidecar written before result-schema versioning is in exactly this state."""
    written = _sidecar_variant("no_version", lambda p: p.pop("result_schema_version"))
    try:
        native_read.load_result(written["keys"])
    except UnsupportedFormat as exc:
        assert "result_schema_version" in str(exc)
        return "refused, not guessed"
    raise AssertionError("an unversioned sidecar was read")


@check("refusal", "a_future_schema_version_is_refused_and_says_to_upgrade")
def _():
    written = _sidecar_variant(
        "future_version",
        lambda p: p.__setitem__("result_schema_version", RESULT_SCHEMA_VERSION + 7),
    )
    try:
        native_read.load_result(written["keys"])
    except UnsupportedFormat as exc:
        message = str(exc)
        assert "newer biotsr" in message and "Upgrade" in message, message
        return "v%d refused by a v%d reader" % (
            RESULT_SCHEMA_VERSION + 7, RESULT_SCHEMA_VERSION,
        )
    raise AssertionError("a future schema version was read")


@check("refusal", "a_non_integer_schema_version_is_refused")
def _():
    for i, bad in enumerate(("1", 1.5, None, [1])):
        written = _sidecar_variant(
            "badver%d" % i, lambda p, b=bad: p.__setitem__("result_schema_version", b)
        )
        try:
            native_read.load_result(written["keys"])
        except (MalformedResult, UnsupportedFormat):
            continue
        raise AssertionError("schema version %r was accepted" % (bad,))
    return "4 malformed spellings refused"


@check("refusal", "a_keys_file_with_no_sidecar_is_refused")
def _():
    written = native.save(COMPUTED["protein_calpha"], workdir("lonely"), triplets=False)
    os.remove(written["metadata"])
    try:
        native_read.load_result(written["keys"])
    except FileNotFoundError as exc:
        assert "key space" in str(exc)
        return "the keys alone do not say what they mean"
    raise AssertionError("a keys file was read without its sidecar")


@check("refusal", "a_truncated_keys_file_is_refused")
def _():
    written = native.save(COMPUTED["protein_calpha"], workdir("cut"), triplets=False)
    with open(written["keys"]) as fh:
        head = fh.readlines()[:50]
    with open(written["keys"], "w") as fh:
        fh.writelines(head)
    try:
        native_read.load_result(written["keys"])
    except MalformedResult as exc:
        assert "n_unique_keys" in str(exc)
        return "caught by the sidecar's own count"
    raise AssertionError("a truncated keys file was read")


@check("refusal", "a_keys_file_out_of_order_is_refused")
def _():
    written = native.save(COMPUTED["protein_calpha"], workdir("shuffled"),
                          triplets=False)
    with open(written["keys"]) as fh:
        lines = fh.readlines()
    lines[1], lines[2] = lines[2], lines[1]
    with open(written["keys"], "w") as fh:
        fh.writelines(lines)
    try:
        native_read.load_result(written["keys"])
    except MalformedResult as exc:
        assert "ordering guarantee" in str(exc)
        return "sorted order is a contract, not a coincidence"
    raise AssertionError("an out-of-order keys file was read")


@check("refusal", "a_keys_file_with_a_foreign_header_is_refused")
def _():
    written = native.save(COMPUTED["protein_calpha"], workdir("hdr"), triplets=False)
    with open(written["keys"]) as fh:
        rest = fh.readlines()[1:]
    with open(written["keys"], "w") as fh:
        fh.write("key\tfrequency\n")
        fh.writelines(rest)
    try:
        native_read.load_result(written["keys"])
    except MalformedResult as exc:
        assert "header" in str(exc)
        return "the legacy two-column shape is not silently accepted"
    raise AssertionError("a foreign header was read")


@check("refusal", "an_edited_key_space_is_caught_by_its_own_fingerprint")
def _():
    def mutate(payload):
        payload["key_space"]["mirror_convention"] = "negate_when_prefold_gt_90"

    written = _sidecar_variant("edited_ks", mutate)
    try:
        native_read.load_result(written["keys"])
    except MalformedResult as exc:
        assert "contradicts itself" in str(exc)
        return "a sidecar cannot claim an identity it does not derive"
    raise AssertionError("an edited key space was read")


@check("refusal", "an_edited_sampling_profile_is_caught_by_its_own_fingerprint")
def _():
    def mutate(payload):
        payload["sampling"]["size_filter"] = 12.0

    written = _sidecar_variant("edited_sp", mutate)
    try:
        native_read.load_result(written["keys"])
    except MalformedResult as exc:
        assert "contradicts itself" in str(exc)
        return "sampling too, not only the key space"
    raise AssertionError("an edited sampling profile was read")


@check("refusal", "two_different_results_cannot_be_combined")
def _():
    """A mismatched pair: ``sample_A.keys`` beside ``sample_B.meta.json``."""
    a = COMPUTED["protein_calpha"]
    b = COMPUTED["protein_calpha_size_filtered"]
    wa = native.save(a, workdir("pair_a"), triplets=False)
    wb = native.save(b, workdir("pair_b"), triplets=False)
    assert a.counts.n_triangles != b.counts.n_triangles
    try:
        native_read.load_result(wa["keys"], metadata_path=wb["metadata"])
    except MalformedResult as exc:
        assert "do not belong to the same result" in str(exc)
        return "caught on n_triangles"
    raise AssertionError("two results' files were combined into one")


@check("refusal", "a_release_candidate_file_is_refused_by_name")
def _():
    """``1UBQ_A.keys_29_35`` (1.0.0rc3 and earlier) is never read as current,
    whether the keys file or its sidecar is named, and even beside a sidecar
    that would otherwise be readable."""
    written = native.save(COMPUTED["protein_calpha"], workdir("old_layout"),
                          triplets=False)
    old = os.path.join(os.path.dirname(written["keys"]), "1UBQ_A.keys_29_35")
    os.rename(written["keys"], old)
    os.rename(written["metadata"], old + ".meta.json")
    for path, extra in ((old, {}), (old + ".meta.json", {}),
                        (old, {"metadata_path": old + ".meta.json"})):
        try:
            native_read.load_result(path, **extra)
        except UnsupportedFormat as exc:
            assert "release candidate" in str(exc), str(exc)[:200]
            continue
        raise AssertionError("%s was read" % os.path.basename(path))
    return "keys, sidecar and an explicit pairing all refused"


@check("refusal", "a_schema_1_sidecar_is_refused_as_a_release_candidate_result")
def _():
    written = _sidecar_variant("schema_1",
                               lambda p: p.__setitem__("result_schema_version", 1))
    try:
        native_read.load_result(written["keys"])
    except UnsupportedFormat as exc:
        assert "1.0.0rc3 or earlier" in str(exc), str(exc)[:200]
        return "schema 1 refused; this build reads schema %d" % RESULT_SCHEMA_VERSION
    raise AssertionError("a schema 1 sidecar was read")


@check("refusal", "a_renamed_keys_file_is_never_paired_by_guessing")
def _():
    """A prefix carries no facts the reader trusts: a keys file renamed away
    from its sidecar is refused (no sidecar of its prefix), and pairing it
    explicitly with another result's sidecar fails on the shared quantities."""
    written = native.save(COMPUTED["protein_calpha"], workdir("renamed"),
                          triplets=False)
    moved = os.path.join(os.path.dirname(written["keys"]), "9XYZ_Z_protein.keys")
    os.rename(written["keys"], moved)
    try:
        native_read.load_result(moved)
    except FileNotFoundError as exc:
        assert "9XYZ_Z_protein.meta.json" in str(exc), str(exc)[:200]
    else:
        raise AssertionError("a keys file was read without the sidecar of its prefix")
    other = native.save(COMPUTED["protein_calpha_mirror"], workdir("renamed_other"),
                        triplets=False)
    try:
        native_read.load_result(moved, metadata_path=other["metadata"])
    except MalformedResult as exc:
        assert "do not belong to the same result" in str(exc), str(exc)[:200]
        return "no sidecar of its prefix; another result's sidecar refused"
    raise AssertionError("a renamed keys file was paired with another result's sidecar")


@check("refusal", "every_corrupt_sidecar_raises_MalformedResult_and_not_something_else")
def _():
    """The contract is a type, not just a message.

    The sub-record ``from_dict`` methods index the fields they need, so a
    truncated ``counts`` block naturally raises ``KeyError`` and an illegal
    ``structure_id`` raises ``ValueError`` from ``StructureId``. Both are correct
    diagnoses arriving under the wrong name: a caller writing
    ``except MalformedResult`` would miss them, and the CLI would blame the
    request rather than the file. The reader translates at its boundary.
    """
    cases = {
        "counts missing a field": lambda p: p["counts"].pop("n_points"),
        "an illegal structure_id": lambda p: p.__setitem__(
            "structure_id", "with space"
        ),
        "a non-string structure_id": lambda p: p.__setitem__("structure_id", 7),
        "source missing a field": lambda p: p["source"].pop("parser_backend"),
        "selection missing a field": lambda p: p["selection"].pop("n_points"),
        "a malformed skip record": lambda p: p.__setitem__(
            "skipped", [{"stage": "x"}]
        ),
        "counts not an object": lambda p: p.__setitem__("counts", [1, 2]),
        "category_totals not an object": lambda p: p.__setitem__(
            "category_totals", ["regular", 1]
        ),
    }
    for i, (name, mutate) in enumerate(sorted(cases.items())):
        written = _sidecar_variant("corrupt%d" % i, mutate)
        try:
            native_read.load_result(written["keys"])
        except MalformedResult:
            continue
        except Exception as exc:  # noqa: BLE001
            raise AssertionError(
                "%s raised %s, not MalformedResult" % (name, type(exc).__name__)
            )
        raise AssertionError("%s was accepted" % name)
    return "%d corruptions, all MalformedResult" % len(cases)


@check("refusal", "a_key_space_field_is_never_defaulted")
def _():
    ks = dict(COMPUTED["protein_calpha"].key_space.to_dict())
    for name in list(ks):
        payload = dict(ks)
        del payload[name]
        try:
            KeySpace.from_dict(payload)
        except MalformedResult:
            continue
        raise AssertionError("key_space field %r was defaulted" % name)
    return "%d fields, none defaultable" % len(ks)


@check("refusal", "a_sampling_field_is_never_defaulted")
def _():
    sp = dict(COMPUTED["protein_calpha"].sampling.to_dict())
    for name in list(sp):
        payload = dict(sp)
        del payload[name]
        try:
            SamplingProfile.from_dict(payload)
        except MalformedResult:
            continue
        raise AssertionError("sampling field %r was defaulted" % name)
    return "%d fields, none defaultable" % len(sp)


@check("refusal", "an_unknown_identity_field_is_refused_not_ignored")
def _():
    """Silently dropping a field would mean reading a key space partially."""
    for cls, block in ((KeySpace, "key_space"), (SamplingProfile, "sampling")):
        payload = dict(getattr(COMPUTED["protein_calpha"], block).to_dict())
        payload["a_field_from_the_future"] = "x"
        try:
            cls.from_dict(payload)
        except MalformedResult:
            continue
        raise AssertionError("%s accepted an unknown field" % cls.__name__)
    return "both identities"


# ===========================================================================
# analysis -- established measures, now from saved results
# ===========================================================================

_ANALYSIS_PAIR = ("1UBQ", "3PTB")


def _analysis_samples():
    live, loaded = [], []
    for pdb in _ANALYSIS_PAIR:
        result = biotsr.analyze(fixture(pdb), method="protein", chain="A",
                             keep_triplets=False)
        _l, _w = save_load(result, "an_%s" % pdb)
        live.append(AnalysisSample.from_result(result))
        loaded.append(AnalysisSample.from_result(_l))
    return live, loaded


LIVE_SAMPLES, LOADED_SAMPLES = _analysis_samples()


@check("analysis", "common_keys_from_loaded_results_are_identical")
def _():
    a = common_key_intersection(LIVE_SAMPLES)
    b = common_key_intersection(LOADED_SAMPLES)
    assert a.keys == b.keys, "the intersection changed"
    assert [(r.sample_id, r.total_keys, r.n_common_keys, r.sum_common_key_frequency)
            for r in a.rows] == [
        (r.sample_id, r.total_keys, r.n_common_keys, r.sum_common_key_frequency)
        for r in b.rows
    ]
    return "%d common keys over %d samples" % (b.n_common_keys, b.n_samples)


@check("analysis", "feature_vectors_from_loaded_results_are_identical")
def _():
    a, b = feature_matrix(LIVE_SAMPLES), feature_matrix(LOADED_SAMPLES)
    assert a.columns == b.columns, "column order changed"
    assert a.rows == b.rows, "a row changed"
    assert a.key_space.fingerprint == b.key_space.fingerprint
    return "%d columns" % len(b.columns)


@check("analysis", "generalized_jaccard_from_loaded_results_is_identical")
def _():
    a = generalized_jaccard(LIVE_SAMPLES[0].counts, LIVE_SAMPLES[1].counts)
    b = generalized_jaccard(LOADED_SAMPLES[0].counts, LOADED_SAMPLES[1].counts)
    assert a == b, "%r vs %r" % (a, b)
    return "%.12f" % b


@check("analysis", "every_similarity_measure_from_loaded_results_is_identical")
def _():
    from biotsr.analysis import MEASURES

    a, b = feature_matrix(LIVE_SAMPLES), feature_matrix(LOADED_SAMPLES)
    for measure in sorted(MEASURES):
        assert (
            similarity_matrix(a, measure).values == similarity_matrix(b, measure).values
        ), measure
    return "%d measures" % len(MEASURES)


@check("analysis", "the_distance_matrix_from_loaded_results_is_identical")
def _():
    a, b = feature_matrix(LIVE_SAMPLES), feature_matrix(LOADED_SAMPLES)
    assert distance_matrix(a).values == distance_matrix(b).values
    return "identical"


@check("analysis", "clustering_from_loaded_results_is_identical")
def _():
    """Four samples, so average linkage and maxclust have something to do.

    The class labels are supplied rather than derived: `label_from_name` takes a
    token index because the legacy consumers disagree about which one it is, and
    a gate must not invent a convention either.
    """
    from biotsr.analysis import clustermap_spec

    cases = [("1UBQ", "A", "small"), ("3PTB", "A", "large"),
             ("1BRS", "A", "large"), ("1BRS", "B", "large")]
    live, saved = [], []
    for pdb, chain, label in cases:
        result = biotsr.analyze(fixture(pdb), method="protein", chain=chain,
                             keep_triplets=False)
        loaded_result, _w = save_load(result, "cl_%s_%s" % (pdb, chain))
        sid = "%s_%s" % (pdb, chain)
        live.append(AnalysisSample.from_result(result, sample_id=sid, label=label))
        loaded = AnalysisSample.from_result(loaded_result, sample_id=sid, label=label)
        saved.append(loaded)
    # Two clusters, chosen here: BioTSR never infers the count.
    a = clustermap_spec(distance_matrix(feature_matrix(live)),
                        [s.label for s in live], n_clusters=2)
    b = clustermap_spec(distance_matrix(feature_matrix(saved)),
                        [s.label for s in saved], n_clusters=2)
    assert a.assignment.leaves == b.assignment.leaves, "leaf order changed"
    assert a.assignment.labels == b.assignment.labels, "cluster assignment changed"
    assert a.adjusted_rand == b.adjusted_rand, "the ARI changed"
    return "%d leaves, ARI %.6f" % (len(b.assignment.leaves), b.adjusted_rand)


@check("analysis", "the_uncut_tree_from_loaded_results_is_identical")
def _():
    """The default hierarchical workflow needs no cluster count (ADR-099): the
    tree, its leaf order and the clustermap without a partition are the same
    on loaded results, and no partition is invented."""
    from biotsr.analysis import clustermap_spec, hierarchy

    cases = [("1UBQ", "A", "small"), ("3PTB", "A", "large"),
             ("1BRS", "A", "large"), ("1BRS", "B", "large")]
    live, saved = [], []
    for pdb, chain, label in cases:
        result = biotsr.analyze(fixture(pdb), method="protein", chain=chain,
                                keep_triplets=False)
        loaded_result, _w = save_load(result, "tree_%s_%s" % (pdb, chain))
        sid = "%s_%s" % (pdb, chain)
        live.append(AnalysisSample.from_result(result, sample_id=sid, label=label))
        saved.append(AnalysisSample.from_result(loaded_result, sample_id=sid,
                                                label=label))
    da = distance_matrix(feature_matrix(live))
    db = distance_matrix(feature_matrix(saved))
    ta, tb = hierarchy(da), hierarchy(db)
    assert ta.to_linkage_csv() == tb.to_linkage_csv(), "the merges changed"
    assert ta.ordered_sample_ids() == tb.ordered_sample_ids(), "leaf order changed"
    assert tb.provenance()["flat_partition"] is None
    spec_a = clustermap_spec(da, [s.label for s in live])
    spec_b = clustermap_spec(db, [s.label for s in saved])
    for spec in (spec_a, spec_b):
        assert spec.assignment is None and spec.adjusted_rand is None
        assert spec.cut_description() is None
        try:
            spec.to_assignments_csv()
        except ValueError:
            pass
        else:
            raise AssertionError("an uncut clustermap produced cluster assignments")
    assert spec_a.ordered_sample_ids() == spec_b.ordered_sample_ids()
    assert spec_a.row_class_slots() == spec_b.row_class_slots()
    assert spec_b.ordered_sample_ids() == tb.ordered_sample_ids()
    return "%d leaves, %d merges, no partition" % (
        tb.n_samples, len(tb.linkage.rows))


@check("analysis", "loading_does_not_launder_an_incompatible_key_space")
def _():
    from biotsr.analysis import require_common_key_space

    plain, _w = save_load(COMPUTED["protein_calpha"], "inc_plain")
    grouped, _w2 = save_load(COMPUTED["protein_calpha_grouped"], "inc_grouped")
    try:
        require_common_key_space(
            [
                AnalysisSample.from_result(plain, sample_id="plain"),
                AnalysisSample.from_result(grouped, sample_id="grouped"),
            ]
        )
    except biotsr.IncompatibleKeySpace:
        return "refused with the differing field named"
    raise AssertionError("two key spaces were combined after loading")


@check("analysis", "the_cli_analyses_saved_results_identically_to_recomputing_them")
def _():
    out = workdir("cli_saved")
    for pdb in _ANALYSIS_PAIR:
        code, _o, err = cli(
            ["run", "-m", "protein", "--pdb", fixture(pdb), "-c", "A",
             "-o", os.path.join(out, "runs"), "--no-triplets", "-q"]
        )
        assert code == EXIT_OK, err
    code, _o, err = cli(
        ["analyze", "common-keys", "--from-results", os.path.join(out, "runs"),
         "-o", os.path.join(out, "from_saved"), "-q"]
    )
    assert code == EXIT_OK, err

    # The manifest's sample ids are set to the output prefixes `biotsr run`
    # chose, because that is what `--from-results` reports: each calculation's
    # completion record carries its prefix, and the prefix is the sample id.
    # Aligning them is what makes a byte comparison meaningful rather than
    # accidental.
    manifest = os.path.join(out, "m.csv")
    with open(manifest, "w") as fh:
        fh.write("sample_id,pdb_path,chain,method\n")
        for pdb in _ANALYSIS_PAIR:
            fh.write("%s_A_protein,%s,A,protein\n" % (pdb, fixture(pdb)))
    code, _o, err = cli(
        ["analyze", "common-keys", "--manifest", manifest,
         "-o", os.path.join(out, "recomputed"), "-q"]
    )
    assert code == EXIT_OK, err

    for name in ("all_samples.json", "all_samples.csv"):
        a = read_bytes(os.path.join(out, "from_saved", "analysis", "common_keys", name))
        b = read_bytes(os.path.join(out, "recomputed", "analysis", "common_keys", name))
        assert a == b, "%s differs between the two paths" % name
    return "both reports byte-identical"


@check("analysis", "an_api_loaded_result_equals_a_cli_loaded_result")
def _():
    """API == CLI, the extended form: the same request, saved and loaded both ways.

    The command-line gates established API-compute == CLI-compute.  This adds a save/load
    step on each side, and this asserts the two loaded objects are the same
    result -- not that each is self-consistent, which would prove nothing about
    the other.

    Every scientifically meaningful field is compared.  ``source.path`` is not:
    the two sides are handed the same file by two different code paths and one
    of them may spell it relatively, which is provenance about the invocation
    rather than about the science.  The structure's SHA-256 IS compared, and that
    is the field that says they read the same bytes.
    """
    out = workdir("api_vs_cli_loaded")
    pdb = fixture("1UBQ")

    api_result = biotsr.analyze(pdb, method="protein", chain="A", mirror=True,
                             keep_triplets=False)
    api_written = native.save(api_result, os.path.join(out, "api"), triplets=False)
    api_loaded = native_read.load_result(api_written["keys"])

    code, stdout, err = cli(
        ["run", "-m", "protein", "--pdb", pdb, "-c", "A", "--mirror",
         "-o", os.path.join(out, "cli"), "--no-triplets", "-q"]
    )
    assert code == EXIT_OK, err
    with open(stdout.strip(), encoding="utf-8") as fh:
        record = json.load(fh)
    assert record["complete"] is True
    [entry] = record["results"]
    cli_keys = os.path.join(out, "cli", entry["files"]["keys"])
    # The CLI names the request exactly as the API names the result.
    assert os.path.basename(cli_keys) == os.path.basename(api_written["keys"]), (
        os.path.basename(cli_keys), os.path.basename(api_written["keys"]))
    cli_loaded = native_read.load_result(cli_keys)

    assert cli_loaded.method == api_loaded.method
    assert cli_loaded.structure_id == api_loaded.structure_id
    assert cli_loaded.chain == api_loaded.chain
    assert {k: dict(v.by_category) for k, v in cli_loaded.keys.items()} == {
        k: dict(v.by_category) for k, v in api_loaded.keys.items()
    }, "the key/category tables differ"
    assert cli_loaded.key_space == api_loaded.key_space
    assert cli_loaded.sampling == api_loaded.sampling
    assert cli_loaded.key_space.fingerprint == api_loaded.key_space.fingerprint
    assert cli_loaded.sampling.fingerprint == api_loaded.sampling.fingerprint
    assert cli_loaded.counts == api_loaded.counts
    assert cli_loaded.selection == api_loaded.selection
    assert cli_loaded.source.sha256 == api_loaded.source.sha256
    assert cli_loaded.source.parser_backend == api_loaded.source.parser_backend
    assert cli_loaded.category_totals() == api_loaded.category_totals()
    for key in ("profile", "base_profile", "legacy_compatibility_profile"):
        assert cli_loaded.annotations[key] == api_loaded.annotations[key], key
    assert tuple(cli_loaded.annotations["modifiers"]) == tuple(
        api_loaded.annotations["modifiers"]
    )

    # And the bytes: the two saved key files must be identical, not merely
    # equivalent once parsed.
    assert read_bytes(api_written["keys"]) == read_bytes(cli_keys), (
        "the two key files differ byte for byte")
    return "%d keys, keyspace %s" % (
        len(api_loaded.keys), api_loaded.key_space.fingerprint[:16],
    )


@check("analysis", "the_cli_refuses_a_dataset_named_two_ways_at_once")
def _():
    code, _o, err = cli(
        ["analyze", "common-keys", "--manifest", "m.csv",
         "--from-results", TMP, "-o", workdir("both")]
    )
    assert code == EXIT_USAGE, "exit %d" % code
    assert "both name a dataset" in err, err
    return "exit 2"


@check("analysis", "a_relocated_folder_keeps_every_identity_and_every_row_label")
def _():
    """The row label is the calculation's prefix, read from its completion
    record -- not the folder it sits in.

    Moving or renaming the output folder changes no calculation: the label and
    everything the serialized result asserts about itself survive intact.
    Renaming a data file away from its record is not a move but a damaged
    calculation, and is refused.
    """
    from biotsr.cli.commands import _saved_results
    from biotsr.cli.console import QUIET, Console
    from biotsr.cli.errors import CLIError

    root = os.path.join(TMP, "relocate_root")
    code, _o, err = cli(
        ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A", "-o", root,
         "--sample-id", "first_name", "--no-triplets", "-q"])
    assert code == EXIT_OK, err
    [(sample_id, before, _label, keys)] = _saved_results(root, Console(QUIET))
    assert sample_id == "first_name" and os.path.basename(keys) == "first_name.keys"

    moved_root = os.path.join(workdir("relocate_new_root"), "a totally different name")
    shutil.move(root, moved_root)
    [(sample_id, after, _label, moved_keys)] = _saved_results(moved_root,
                                                               Console(QUIET))
    # The row label stays the prefix, wherever the folder goes ...
    assert sample_id == "first_name", sample_id
    assert os.path.dirname(moved_keys) == moved_root
    # ... and nothing the result asserts about itself moves either.
    for name in ("structure_id", "chain", "method", "counts", "selection",
                 "source", "key_space", "sampling"):
        assert getattr(after, name) == getattr(before, name), name
    assert after.key_space.fingerprint == before.key_space.fingerprint
    assert after.sampling.fingerprint == before.sampling.fingerprint
    assert dict(after.frequencies) == dict(before.frequencies)
    for key in ("profile", "base_profile", "legacy_compatibility_profile",
                "target"):
        assert after.annotations[key] == before.annotations[key], key
    assert after.structure_id.value == "1UBQ" and after.chain == "A"

    # Re-serializing a relocated result reproduces the original bytes: the new
    # location reached `loaded_from` and nothing that gets written.
    a = native.save(before, workdir("relocate_ser_a"), triplets=False)
    b = native.save(after, workdir("relocate_ser_b"), triplets=False)
    for role in ("keys", "metadata"):
        assert read_bytes(a[role]) == read_bytes(b[role]), role

    # A data file renamed away from its completion record is refused.
    os.rename(moved_keys, os.path.join(moved_root, "renamed.keys"))
    try:
        _saved_results(moved_root, Console(QUIET))
    except CLIError as exc:
        assert "not complete" in exc.message, exc.message
        return "identity and label intact after the move; a renamed file refused"
    raise AssertionError("a calculation with a renamed keys file was read")


@check("analysis", "the_cli_refuses_two_saved_results_that_would_share_an_id")
def _():
    """Two output folders of one tree can hold calculations with one prefix.

    `--from-results` reads a tree recursively and reports each calculation
    under its prefix.  A user who gathers two `--output` folders under one root
    can therefore have two `1UBQ_A_protein` calculations -- here, by giving two
    different requests the same --sample-id -- and the report would have two
    rows under one name that cannot be matched back to what produced them.
    """
    merged = workdir("merged_roots")
    for i, extra in enumerate((["--sample-id", "1UBQ_A_protein"],
                               ["--size-filter", "14", "--sample-id",
                                "1UBQ_A_protein"])):
        code, _o, err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", os.path.join(merged, "run%d" % i), "--no-triplets", "-q"] + extra
        )
        assert code == EXIT_OK, err
    code, _o, err = cli(
        ["analyze", "common-keys", "--from-results", merged,
         "-o", workdir("merged_out")]
    )
    assert code == EXIT_USAGE, "exit %d" % code
    assert "would both be reported as" in err, err
    return "exit 2, naming both files"


@check("analysis", "the_cli_refuses_a_results_directory_with_nothing_in_it")
def _():
    empty = workdir("empty_results")
    code, _o, err = cli(
        ["analyze", "common-keys", "--from-results", empty, "-o", workdir("nothing")]
    )
    assert code == EXIT_USAGE, "exit %d" % code
    assert "no completed calculations" in err and ".run.json" in err, err
    return "exit 2, naming what a completed calculation is"


# ===========================================================================
# surface -- the public API is what it says it is
# ===========================================================================


@check("surface", "every_name_in_tsr___all___resolves")
def _():
    missing = [n for n in biotsr.__all__ if not hasattr(biotsr, n)]
    assert missing == [], missing
    return "%d names" % len(biotsr.__all__)


@check("surface", "the_documented_import_paths_work")
def _():
    namespace = {}
    exec(
        "from biotsr import analyze, analyze_cross, load_result, TSRResult, KeySpace\n"
        "from biotsr.analysis import AnalysisSample, feature_matrix, distance_matrix\n"
        "from biotsr.sse import require_common_strand_identity\n"
        "from biotsr.io import save, load_result as io_load_result\n"
        "from biotsr.profiles import STANDARD_PROTEIN_CALPHA_V1, Profile\n"
        "from biotsr.cli import main\n",
        namespace,
    )
    return "6 import paths"


@check("surface", "the_loader_is_reachable_from_the_top_level")
def _():
    assert biotsr.load_result is native_read.load_result
    assert "load_result" in biotsr.__all__
    return "biotsr.load_result"


@check("surface", "sse_is_exported_alongside_analysis")
def _():
    assert "sse" in biotsr.__all__ and "analysis" in biotsr.__all__
    return "both sub-packages named"


@check("surface", "every_serializer_module_is_declared_by_tsr_io")
def _():
    from biotsr import io as tsr_io

    for name in ("native", "native_read"):
        assert name in tsr_io.__all__, name
        assert hasattr(tsr_io, name), name
    leftovers = [n for n in tsr_io.__all__ if "legacy" in n]
    assert leftovers == [], "removed legacy writers still declared: %s" % leftovers
    return "%d names, native format only" % len(tsr_io.__all__)


@check("surface", "every_package_refusal_is_a_tsr_error")
def _():
    kinds = (
        biotsr.IncompatibleKeySpace, biotsr.IncompatibleSampling, biotsr.IncompatibleModifier,
        biotsr.IncompatibleSSEStratification, biotsr.MalformedResult, biotsr.UnsupportedFormat,
        biotsr.InvalidRequest,
    )
    for kind in kinds:
        assert issubclass(kind, TSRError), kind.__name__
    return "%d types under TSRError" % len(kinds)


@check("surface", "every_package_refusal_is_still_a_value_error")
def _():
    """Existing code that writes ``except ValueError`` must keep working."""
    for kind in (
        biotsr.IncompatibleKeySpace, biotsr.IncompatibleSampling, biotsr.IncompatibleModifier,
        biotsr.IncompatibleSSEStratification, biotsr.MalformedResult, biotsr.UnsupportedFormat,
        biotsr.InvalidRequest, TSRError,
    ):
        assert issubclass(kind, ValueError), kind.__name__
    return "backward compatible"


@check("surface", "the_four_incompatibility_types_share_a_parent")
def _():
    for kind in (biotsr.IncompatibleKeySpace, biotsr.IncompatibleSampling,
                 biotsr.IncompatibleModifier, biotsr.IncompatibleSSEStratification):
        assert issubclass(kind, IncompatibleRepresentation), kind.__name__
    return "IncompatibleRepresentation"


@check("surface", "the_comparison_rule_has_one_entry_point_and_one_opt_in")
def _():
    """ADR-094: one ``check_comparable``; a processing refusal is an
    ``IncompatibleSampling`` and never an ``IncompatibleKeySpace``; every
    analysis entry point takes ``allow_mixed_processing`` and none still takes
    the rc3 ``allow_mixed_sampling``."""
    import inspect

    from biotsr import analysis

    assert biotsr.IncompatibleProcessing is analysis.IncompatibleProcessing
    assert issubclass(analysis.IncompatibleProcessing, biotsr.IncompatibleSampling)
    assert not issubclass(analysis.IncompatibleProcessing, biotsr.IncompatibleKeySpace)
    entry_points = (analysis.check_comparable, analysis.feature_matrix,
                    analysis.common_key_intersection, analysis.key_occurrence_percent)
    for fn in entry_points:
        parameters = inspect.signature(fn).parameters
        assert "allow_mixed_processing" in parameters, fn.__name__
        assert "allow_mixed_sampling" not in parameters, fn.__name__
    return "%d entry points take allow_mixed_processing" % len(entry_points)


@check("surface", "a_programming_error_is_not_a_tsr_error")
def _():
    """The exclusion is what gives ``except TSRError`` its meaning."""
    for kind in (AttributeError, TypeError, AssertionError, KeyError, RuntimeError):
        assert not issubclass(kind, TSRError), kind.__name__
    return "5 stdlib types stay outside"


@check("surface", "the_cli_no_longer_keeps_its_own_error_taxonomy")
def _():
    """It caught five package exception types by name; now it catches the base."""
    # ``biotsr.cli.main`` names two things: the module, and the function that
    # ``biotsr/cli/__init__.py`` re-exports over it. ``import biotsr.cli.main as x``
    # binds the attribute, i.e. the function. Ask the module registry instead.
    import importlib

    cli_module = importlib.import_module("biotsr.cli.main")

    assert TSRError in cli_module._USER_ERRORS
    for name in ("IncompatibleKeySpace", "IncompatibleSampling",
                 "IncompatibleModifier"):
        assert name not in [t.__name__ for t in cli_module._USER_ERRORS], name
    return "%d entries" % len(cli_module._USER_ERRORS)


@check("surface", "the_python_api_validates_a_profile_against_its_method")
def _():
    """The CLI could not express this; Python could."""
    try:
        biotsr.analyze(fixture("1UBQ"), method="drug", chain="A", entity_name="BEN",
                    profile=biotsr.STANDARD_PROTEIN_CALPHA_V1)
    except InvalidRequest as exc:
        message = str(exc)
        assert "standard_protein_calpha_v1" in message, message
        assert "standard_drug_v1" in message, message
        return "refused, naming the given and the standard profile"
    raise AssertionError("a residue profile ran a method that produces atoms")


@check("surface", "a_derived_or_named_profile_is_refused")
def _():
    """Every method runs under its own standard profile and nothing else.

    A profile built by ``replace()`` -- even one in the right label domain -- or
    a profile given by name defines keys nobody has defined for the method, so
    both are refused before any structure is read.
    """
    derived = replace(biotsr.STANDARD_PROTEIN_CALPHA_V1, name="derived_for_a_study")
    refused = 0
    for profile in (derived, "standard_protein_calpha_v1"):
        try:
            biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                           profile=profile, keep_triplets=False)
        except InvalidRequest:
            refused += 1
            continue
        raise AssertionError("profile %r was accepted" % (profile,))
    out = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                         profile=biotsr.STANDARD_PROTEIN_CALPHA_V1,
                         keep_triplets=False)
    assert out.annotations["profile"] == "standard_protein_calpha_v1"
    return "%d refused; the standard profile itself is accepted" % refused


@check("surface", "every_method_accepts_its_own_default_profile")
def _():
    from biotsr.api import _require_standard_profile

    for name in biotsr.list_methods():
        recipe = get_method(name)
        _require_standard_profile(recipe, recipe.default_profile)
        _require_standard_profile(recipe, None)
    return "%d methods" % len(biotsr.list_methods())


@check("surface", "the_python_api_validates_every_scientific_argument_the_cli_does")
def _():
    """Direct Python users get the same refusals, from the package."""
    cases = []

    def refuses(what, fn):
        try:
            fn()
        except (TSRError, ValueError, KeyError):
            cases.append(what)
            return
        raise AssertionError("%s was accepted" % what)

    refuses("a non-positive size filter",
            lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                                size_filter=0))
    refuses("a non-positive contact cutoff",
            lambda: biotsr.analyze_cross(fixture("3PTB"), method="drug_protein_cross",
                                      chain="A", chain_2="A", entity_name="BEN",
                                      entity_id="1", contact_cutoff=0))
    refuses("an unknown method",
            lambda: biotsr.analyze(fixture("1UBQ"), method="no_such_method", chain="A"))
    refuses("an unknown chain",
            lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="ZZ"))
    refuses("a modifier the method does not support",
            lambda: biotsr.analyze(fixture("3PTB"), method="drug", chain="A",
                                entity_name="BEN", entity_id="1", sse=True))
    refuses("aa_grouping outside the amino-acid domain",
            lambda: biotsr.analyze(fixture("3PTB"), method="drug", chain="A",
                                entity_name="BEN", entity_id="1", aa_grouping=True))
    refuses("a profile given by name",
            lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                                profile="no_such_profile"))
    refuses("an argument the method does not read",
            lambda: biotsr.analyze(fixture("1UBQ"), method="protein", chain="A",
                                entity_name="BEN"))
    refuses("a profile from the wrong label domain",
            lambda: biotsr.analyze(fixture("1UBQ"), method="drug", chain="A",
                                entity_name="BEN",
                                profile=biotsr.STANDARD_PROTEIN_CALPHA_V1))
    refuses("an unknown parser backend",
            lambda: biotsr.load_structure(fixture("1UBQ"), backend="mmcif"))
    return "%d refusals, all from the package" % len(cases)


@check("surface", "defaults_have_one_source_of_truth")
def _():
    """No profile is built outside the package's profile definitions.

    Both are already asserted by the command-line tests; this gate states the
    property itself so that the API-hardening promise is checkable in one place.
    """
    cli_dir = os.path.join(SRC, "biotsr", "cli")
    offenders = []
    for name in sorted(os.listdir(cli_dir)):
        if not name.endswith(".py"):
            continue
        with open(os.path.join(cli_dir, name)) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == "Profile":
                    offenders.append("%s builds a Profile" % name)
    assert offenders == [], offenders
    # The profile each method runs under is the one its recipe declares; the
    # CLI's method listing reports that same object's name, not a restatement.
    code, text, _err = cli(["methods", "--json"])
    assert code == EXIT_OK
    listed = {r["name"]: r["profile"] for r in json.loads(text)}
    for method in biotsr.list_methods():
        assert listed[method] == get_method(method).default_profile.name, method
    return "%d methods resolve their profile through their own recipe" % len(
        biotsr.list_methods()
    )


@check("surface", "a_result_cannot_be_mutated_through_its_own_metadata")
def _():
    """A result's metadata is a copy -- checked directly on the round-trip path."""
    annotated = COMPUTED["protein_calpha_sse"]
    before = copy.deepcopy(annotated.annotations["sse"])
    payload = annotated.metadata()
    payload["sse"]["strand_identity"] = "tampered"
    payload["provenance"]["profile"] = "tampered"
    assert annotated.annotations["sse"] == before, "the SSE block was shared"
    assert annotated.annotations["profile"] != "tampered", "provenance was shared"
    return "metadata() copies what it embeds"


@check("surface", "two_results_never_share_a_category_counter")
def _():
    a = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A", keep_triplets=False)
    b = biotsr.analyze(fixture("1UBQ"), method="protein", chain="A", keep_triplets=False)
    key = sorted(a.keys)[0]
    assert a.keys[key] is not b.keys[key]
    assert a.keys[key].by_category is not b.keys[key].by_category
    return "the mirror-convention validation's aliasing fix still holds"


@check("surface", "a_loaded_result_shares_nothing_with_the_file_it_came_from")
def _():
    loaded, written = save_load(COMPUTED["protein_calpha"], "no_share")
    with open(written["metadata"]) as fh:
        payload = json.load(fh)
    payload["sampling"]["size_filter"] = 1.0
    payload["key_space"]["mirror_convention"] = "x"
    assert loaded.sampling.size_filter is None
    assert loaded.key_space.mirror_convention != "x"
    return "the loader copies, it does not alias"


# ===========================================================================
# robustness -- filesystem and packaging
# ===========================================================================


@check("robustness", "a_path_with_spaces_round_trips")
def _():
    out = os.path.join(TMP, "a directory with spaces", "and another")
    written = native.save(COMPUTED["protein_calpha"], out, triplets=False)
    loaded = native_read.load_result(written["keys"])
    assert loaded.key_space.fingerprint == (
        COMPUTED["protein_calpha"].key_space.fingerprint
    )
    return "saved and loaded"


@check("robustness", "a_path_with_non_ascii_characters_round_trips")
def _():
    out = os.path.join(TMP, "étude-使用者", "résultats")
    written = native.save(COMPUTED["protein_calpha"], out, triplets=False)
    loaded = native_read.load_result(written["keys"])
    assert loaded.counts == COMPUTED["protein_calpha"].counts
    return os.path.basename(os.path.dirname(written["keys"]))


@check("robustness", "a_non_ascii_source_path_survives_into_the_sidecar")
def _():
    """``source.path`` is user-controlled and reaches the JSON; UTF-8 is explicit."""
    staging = os.path.join(TMP, "étude-input")
    os.makedirs(staging, exist_ok=True)
    pdb = os.path.join(staging, "1UBQ.pdb")
    shutil.copy(fixture("1UBQ"), pdb)
    result = biotsr.analyze(pdb, method="protein", chain="A", keep_triplets=False)
    written = native.save(result, workdir("nonascii_src"), triplets=False)
    loaded = native_read.load_result(written["keys"])
    assert loaded.source.path == result.source.path
    assert "é" in loaded.source.path
    # ``json.dump`` escapes non-ASCII by default, so the sidecar's BYTES stay
    # ASCII (``\u00e9``) while its VALUE is the character. That is the stronger
    # position and is kept deliberately: an ASCII-only sidecar is byte-stable on
    # every machine whatever the locale. What must hold is that the value
    # survives -- asserted above -- and that the reader decodes it as UTF-8
    # rather than as whatever the locale happens to be.
    with open(written["metadata"], "rb") as fh:
        raw = fh.read()
    raw.decode("ascii")  # raises if the sidecar is not ASCII-safe
    assert b"\\u00e9" in raw, "the non-ASCII path was not escaped"
    return "value preserved, bytes ASCII-stable"


@check("robustness", "every_native_file_declares_its_encoding")
def _():
    """A format whose bytes depend on the machine's locale is not a format."""
    offenders = []
    for name in ("native.py", "native_read.py"):
        path = os.path.join(SRC, "biotsr", "io", name)
        with open(path) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id == "open":
                    if not any(kw.arg == "encoding" for kw in node.keywords):
                        offenders.append("%s line %d" % (name, node.lineno))
    assert offenders == [], offenders
    return "every open() in the native format"


@check("robustness", "a_read_only_input_is_readable")
def _():
    staging = workdir("readonly")
    pdb = os.path.join(staging, "ro.pdb")
    shutil.copy(fixture("1UBQ"), pdb)
    os.chmod(pdb, stat.S_IRUSR)
    try:
        result = biotsr.analyze(pdb, method="protein", chain="A", keep_triplets=False)
        written = native.save(result, workdir("readonly_out"), triplets=False)
        os.chmod(written["keys"], stat.S_IRUSR)
        os.chmod(written["metadata"], stat.S_IRUSR)
        loaded = native_read.load_result(written["keys"])
        assert loaded.counts.n_triangles == result.counts.n_triangles
    finally:
        os.chmod(pdb, stat.S_IRUSR | stat.S_IWUSR)
    return "input and saved result both read-only"


@check("robustness", "a_relative_path_resolves_against_the_working_directory")
def _():
    here = os.getcwd()
    work = workdir("relative")
    try:
        os.chdir(work)
        shutil.copy(fixture("1UBQ"), "in.pdb")
        result = biotsr.analyze("in.pdb", method="protein", chain="A",
                             keep_triplets=False)
        written = native.save(result, "out", triplets=False)
        loaded = native_read.load_result(written["keys"])
        assert loaded.counts == result.counts
    finally:
        os.chdir(here)
    return "relative in, relative out"


@check("robustness", "loading_works_from_an_unrelated_working_directory")
def _():
    written = native.save(COMPUTED["protein_calpha"], workdir("abs_load"),
                          triplets=False)
    here = os.getcwd()
    decoy = workdir("decoy")
    os.makedirs(os.path.join(decoy, "validation"), exist_ok=True)
    try:
        os.chdir(decoy)
        loaded = native_read.load_result(os.path.abspath(written["keys"]))
        assert loaded.counts == COMPUTED["protein_calpha"].counts
    finally:
        os.chdir(here)
    return "no dependency on the checkout layout"


@check("robustness", "an_output_root_outside_the_repository_works")
def _():
    outside = tempfile.mkdtemp(prefix="tsr-outside-")
    try:
        code, out, err = cli(
            ["run", "-m", "protein", "--pdb", fixture("1UBQ"), "-c", "A",
             "-o", os.path.join(outside, "results"), "--no-triplets", "-q"]
        )
        assert code == EXIT_OK, err
        record_path = out.strip()
        results = os.path.join(outside, "results")
        assert os.path.dirname(record_path) == results, record_path
        with open(record_path, encoding="utf-8") as fh:
            record = json.load(fh)
        assert record["complete"] is True
        loaded = native_read.load_result(
            os.path.join(results, record["results"][0]["files"]["keys"]))
        assert loaded.method == "protein"
    finally:
        shutil.rmtree(outside, ignore_errors=True)
    return "written and read outside the checkout"


@check("robustness", "loading_a_directory_or_a_missing_file_fails_cleanly")
def _():
    lonely = native.save(COMPUTED["protein_calpha"], workdir("only_keys"),
                         triplets=False)
    os.remove(lonely["metadata"])
    cases = (
        # A directory is not a result file: refused by name, before any I/O.
        (TMP, ValueError, "neither a native keys file"),
        (os.path.join(TMP, "nope.keys"), FileNotFoundError, "no native keys file"),
        (os.path.join(TMP, "nope.meta.json"), FileNotFoundError, "no native keys file"),
        (lonely["keys"], FileNotFoundError, "no metadata sidecar"),
    )
    for bad, kind, words in cases:
        try:
            native_read.load_result(bad)
        except kind as exc:
            assert words in str(exc), (bad, str(exc)[:160])
            continue
        raise AssertionError("%r was read" % bad)
    return "%d bad paths, each refused with a message" % len(cases)


@check("robustness", "the_package_declares_no_dependency_beyond_joblib")
def _():
    pyproject = os.path.join(WORKSPACE, "pyproject.toml")
    with open(pyproject) as fh:
        text = fh.read()
    assert '\ndependencies = ["joblib>=1.3"]\n' in text, (
        "the runtime dependencies changed (ADR-086 allows joblib alone)")
    return "dependencies = joblib>=1.3"


@check("robustness", "the_reader_imports_nothing_outside_the_standard_library")
def _():
    path = os.path.join(SRC, "biotsr", "io", "native_read.py")
    with open(path) as fh:
        tree = ast.parse(fh.read())
    allowed = {"__future__", "json", "os", "re", "collections", "typing"}
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] not in allowed:
                    offenders.append(alias.name)
        elif isinstance(node, ast.ImportFrom) and not node.level:
            if (node.module or "").split(".")[0] not in allowed:
                offenders.append(node.module)
    assert offenders == [], offenders
    return "stdlib + relative imports only"


@check("robustness", "the_domain_layer_is_still_a_dependency_sink")
def _():
    """``domain/errors.py`` is new; it must not have broken the layering."""
    domain = os.path.join(SRC, "biotsr", "domain")
    offenders = []
    for name in sorted(os.listdir(domain)):
        if not name.endswith(".py"):
            continue
        with open(os.path.join(domain, name)) as fh:
            tree = ast.parse(fh.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level and node.level > 1:
                offenders.append("%s: level %d" % (name, node.level))
    assert offenders == [], offenders
    return "%d modules, none reaches upward" % len(
        [n for n in os.listdir(domain) if n.endswith(".py")]
    )


#: The 1.0 methods (README, "Method families"); a deferred candidate strategy
#: that became registered would show up as an extra name.
RELEASED_METHODS = frozenset((
    "aa_ind", "aa_ind_mix", "aa_ind_specific", "aa_selection",
    "amino_acid_atom_selection", "amino_acid_atomic", "carbohydrate",
    "carbohydrate_atom_selection", "carbohydrate_chain", "carbohydrate_protein_cross",
    "drug", "drug_atom_selection", "drug_protein_cross", "nucleotide",
    "nucleotide_atom_selection", "nucleotide_chain", "nucleotide_protein_cross",
    "protein", "protein_protein_calpha_cross", "protein_protein_cross",
))


@check("robustness", "no_deferred_candidate_method_became_reachable")
def _():
    names = set(biotsr.list_methods())
    assert names == RELEASED_METHODS, sorted(names ^ RELEASED_METHODS)
    return "%d methods, none deferred" % len(names)


# ===========================================================================


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", dest="json_path", default=None)
    args = parser.parse_args()

    by_category = {}
    for record in RESULTS:
        bucket = by_category.setdefault(record["category"], {"total": 0, "failed": 0})
        bucket["total"] += 1
        bucket["failed"] += not record["passed"]

    current = None
    for record in RESULTS:
        if record["category"] != current:
            current = record["category"]
            print("\n== %s ==" % current)
        print(
            "  %s %-62s %s"
            % ("PASS" if record["passed"] else "FAIL", record["name"], record["detail"])
        )

    failed = sum(b["failed"] for b in by_category.values())
    print()
    print("%d gates, %d failed  %r" % (len(RESULTS), failed, by_category))

    payload = {
        "n_gates": len(RESULTS),
        "n_failed": failed,
        "result_schema_version": RESULT_SCHEMA_VERSION,
        "by_category": by_category,
        "gates": RESULTS,
    }
    path = args.json_path or os.path.join(EVIDENCE, "GATE_RESULTS.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("-> %s" % os.path.relpath(path, WORKSPACE))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
