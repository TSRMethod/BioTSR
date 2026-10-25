"""Scientific identities that must agree on every platform.

Keys are integers, and every KeySpace and SamplingProfile fingerprint is a hash
of configuration, so a standard calculation must give the same identities, the
same distinct-key count and the same key-frequency map on macOS arm64, Linux
x86_64 and Linux aarch64 alike.  (Native triplet files print raw floating-point
geometry, which a host's C library may round in the last digit; they are
compared byte for byte only on macOS arm64, in
``tests/integration/test_standard_references.py``.)
"""

from __future__ import annotations

import hashlib
import os

import pytest

import biotsr

from ..conftest import FIXTURE_PDB, WORKSPACE, _require

F = lambda name: os.path.join(FIXTURE_PDB, name)     # noqa: E731

#: KeySpace fingerprint, SamplingProfile fingerprint, unique keys.  Measured equal
#: on macOS arm64, Linux x86_64 and Linux aarch64 (docs/PLATFORM_VALIDATION.md)
#: up to 1.0.0rc6.  The atomic key-space fingerprints changed in 1.0.0rc7 with the
#: element labels (ADR-113; 1187735d.../dbcde351... before); the keys, the
#: sampling fingerprints and the frequency maps below did not.
IDENTITIES = {
    "calpha_standard": (
        lambda: biotsr.analyze(F("1UBQ.pdb"), method="protein", chain="A", keep_triplets=False),
        "da39d73510c54a4c07529e17f34449789584fad4ea79b020b1badd3b5550d30b",
        "866885436e9a02bbbd52eecfc8372be2191d80ea737a0da509a43d0f3d56511d", 45349),
    "sse_standard": (
        lambda: biotsr.analyze(F("1UBQ.pdb"), method="protein", chain="A", sse=True,
                               keep_triplets=False),
        "da39d73510c54a4c07529e17f34449789584fad4ea79b020b1badd3b5550d30b",
        "866885436e9a02bbbd52eecfc8372be2191d80ea737a0da509a43d0f3d56511d", 45349),
    "aa_atomic_standard_mirror": (
        lambda: biotsr.analyze(F("1UBQ.pdb"), method="amino_acid_atomic", chain="A",
                               residue="1", residue_name="MET", mirror=True,
                               keep_triplets=False),
        "a4ed7b5ac78e690f5f5773d3c5934064d58d14d5983e13b4718329e1658694c7",
        "9f3b52e7dca90beb1c59520707839ed8e1c2a5de4a55e2ba32eb98532c1be726", 54),
    "carbohydrate_standard": (
        lambda: biotsr.analyze(F("5WT9.pdb"), method="carbohydrate", chain="A",
                               entity_name="NAG", entity_id="1", keep_triplets=False),
        "1ecbef1b83344df6d739edf4dc570056aa7e8ea9ce5c43237a746ab16e611ea4",
        "0ede534e4b6324c67d332a2d14a3eee9106f34d7b6daa75d4afb2ffed09ce8da", 247),
    "drug_standard": (
        lambda: biotsr.analyze(F("6G9M.pdb"), method="drug", chain="A", entity_name="ESW",
                               entity_id="404", keep_triplets=False),
        "1ecbef1b83344df6d739edf4dc570056aa7e8ea9ce5c43237a746ab16e611ea4",
        "158a8f785b6e5da0174f96a7adb5acd7873ae66f20cc2d07eb853af097d1cdab", 1880),
    "nucleotide_standard": (
        lambda: biotsr.analyze(F("1BNA.pdb"), method="nucleotide", chain="A", entity_name="DG",
                               entity_id="10", keep_triplets=False),
        "1ecbef1b83344df6d739edf4dc570056aa7e8ea9ce5c43237a746ab16e611ea4",
        "4828568b93a83af28f1aeb60b7b443f29dc6ab94ecd813346dc1aa9f51744380", 977),
}

#: sha256 of ``key<TAB>category<TAB>frequency`` rows, sorted -- the standard
#: results' whole key-frequency map, which no historical file protects.
STANDARD_FREQUENCIES = {
    # The two C-alpha maps at radix 30 (ADR-107), derived from the radix-20 ones
    # (c60bf3352a87..., 6fe4ac13b6d1...) by re-encoding each key; the counts are
    # the same, row for row.
    "calpha_standard": "233009c7b151632dbb690f456fd87403d3bd1e49157cbfa1f4220db36236ce45",
    "sse_standard": "6f7cc6763e6390cb1a05cfc580123bc603b45bf0c486ff7e558e315d6549cd64",
    "aa_atomic_standard_mirror": "9126f57be54b1d666a169ac2cb7fc52bf4946414dc47d5a771c6c8f60a110263",
    "carbohydrate_standard": "ca50ab48a72d1fc2a2c2d81013e5a82c94c2090bbd526b3508478190cd4556c9",
    "drug_standard": "3b825983c1d9284af920ff51a4aa1caafb3f61c5ecdd21b438788578c0ec5f6e",
    "nucleotide_standard": "0c7c4d40e8bfbfad1a44e13b2b500e18684d7244d7ed6e000bfe993ef2e5d4c0",
}


@pytest.mark.parametrize("name", sorted(IDENTITIES))
def test_scientific_identities_are_the_same_on_every_platform(name):
    run, keyspace, sampling, n_unique = IDENTITIES[name]
    result = run()
    assert result.key_space.fingerprint == keyspace
    assert result.sampling.fingerprint == sampling
    assert result.n_unique_keys == n_unique
    if name in STANDARD_FREQUENCIES:
        rows = "".join("%d\t%s\t%d\n" % (k, c, n) for k in sorted(result.keys)
                       for c, n in sorted(result.keys[k].by_category.items()))
        assert hashlib.sha256(rows.encode()).hexdigest() == STANDARD_FREQUENCIES[name]




# ---------------------------------------------------------------------------
# CI wiring
# ---------------------------------------------------------------------------


def _read(*parts):
    with open(_require(os.path.join(WORKSPACE, *parts)), encoding="utf-8") as fh:
        return fh.read()


def test_ci_runs_the_suite_on_linux_and_on_macos_arm64():
    ci = _read(".github", "workflows", "ci.yml")
    assert "ubuntu-24.04" in ci and "macos-15" in ci
    for gone in ("frozen_text_bundle.py", "platform_envelope.py", "frozen-text",
                 "BIOTSR_FROZEN_TEXT_DIR", "fetch-depth: 0"):
        assert gone not in ci, gone


def test_no_temporary_diagnostic_workflow_remains():
    workflows = _require(os.path.join(WORKSPACE, ".github", "workflows"))
    assert sorted(os.listdir(workflows)) == [
        "ci.yml", "publish.yml", "release-validation.yml"]
