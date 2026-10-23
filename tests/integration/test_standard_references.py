"""Valid, unambiguous standard calculations keep the results they had.

``tests/data/reference/standard_baseline.json`` was computed by the last build
that still carried the historical compatibility code (1.0.0rc1), before any of
it was removed, for 42 standard calculations covering all 20 methods and the
supported modifiers.  Each case here recomputes one of
them with the current package and requires the same keys file (by SHA-256), the
same counts, both identity fingerprints, the same residues used and the same
category totals.  The reference is never regenerated to match changed code; a
difference is a finding.

**Protein C-alpha cases since 1.0.0rc5 (ADR-107).**  The residue labels moved to
radix 30, which re-maps every C-alpha key, and the C-alpha selectors now read the
supported modified residues.  ``standard_baseline.json`` is kept exactly as it was;
the C-alpha cases' current values are in ``calpha_base30.json``, which was derived
from the rc1 files by an explicit remap rather than regenerated (see its
description).  Each C-alpha case is checked three ways:

1. **The old selection, through the remap.**  The triangles all of whose
   vertices the previous selection read -- decided by **residue identity** (one
   of the 20 regular amino acids), never by the numeric label, since MSE shares
   MET's 13 -- are decoded and re-encoded at the old radix; they must
   reproduce the frozen rc1 keys file byte for byte.
2. **The added triangles, separately.**  The triangles with a newly read
   modified residue are counted: on 6G9M chain A, CME 161 adds exactly
   C(341, 2) = 57,970.  Their keys were established outside this package (an
   independent calculation of the same triangles) when ``calpha_base30.json``
   was derived; they are never mapped back to the old radix.
3. **The whole updated output, against the combined expectation.**  The keys
   (and, on macOS arm64, triplet) file must equal ``calpha_base30.json``'s
   hash: the remapped old triangles plus the added ones.

**Atomic cases since 1.0.0rc7 (ADR-113, ADR-114).**  Atomic labels come from
each atom's element field, so every atomic key space has a new fingerprint, and
``nucleotide_chain`` includes each nucleotide with the atoms it has.  The
current values are in ``atomic_element_labels.json`` (see its description), and
each atomic case is checked against rc1 as well:

1. every keys file, count, residue list and sampling fingerprint is rc1's own,
   except ``nucleotide_chain_1BNA_A`` (a test below reads the reference file
   to make sure of it);
2. a triplet file is rc1's own, or differs only in how a label name is spelled
   (``Cl`` for rc1's ``CL``): spelled back, it is rc1's file byte for byte;
3. ``nucleotide_chain_1BNA_A``: the triangles without a DC 1 atom reproduce the
   rc1 keys file, and DC 1's 16 atoms add exactly C(243, 3) - C(227, 3) =
   438,216 triangles.

Triplet files carry raw floating-point geometry and are compared by SHA-256 only
on the platform that produced them (macOS arm64); everywhere else the keys, which
are integers, carry the comparison.  Cases over a few hundred thousand triangles
run only with ``-m slow`` (release validation).
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import sys

import pytest

import biotsr
from biotsr.io import native
from biotsr.labels import STANDARD_AMINO_ACIDS

from ..conftest import REFERENCE_DATA, STRUCTURES

_BASELINE = os.path.join(REFERENCE_DATA, "standard_baseline.json")
with open(_BASELINE) as _fh:
    _CASES = json.load(_fh)["cases"]
#: The frozen rc1 file's own SHA-256: it is never edited (ADR-107 added a
#: second file instead).
_BASELINE_SHA256 = "30ef29ea0bf0bbe6db5ede7aeba861251ddd15161eb4519ac307890d0a4071bb"
_CALPHA_BASE30 = os.path.join(REFERENCE_DATA, "calpha_base30.json")
with open(_CALPHA_BASE30) as _fh:
    _CALPHA = {c["request"]["id"]: c for c in json.load(_fh)["cases"]}
_ATOMIC_ELEMENT = os.path.join(REFERENCE_DATA, "atomic_element_labels.json")
with open(_ATOMIC_ELEMENT) as _fh:
    _ATOMIC = {c["request"]["id"]: c for c in json.load(_fh)["cases"]}
#: The one atomic case whose selection changed (ADR-114), and what it added.
_NUCLEOTIDE_CHAIN = "nucleotide_chain_1BNA_A"

#: The positional constants of the key formula, for the remap to the rc1 radix.
_THETA, _DISTANCE, _LOWEST_CODE = 29, 35, 4
_BIN = _THETA * _DISTANCE

#: Above this many triangles a case is part of the slow set.
_SLOW_ABOVE = 250_000

_ON_REFERENCE_PLATFORM = sys.platform == "darwin" and platform.machine() == "arm64"


def _params():
    for case in _CASES:
        expected = case["expected"]
        size = expected.get("n_triangles", expected.get("n_triangles_total", 0))
        marks = [pytest.mark.slow] if size > _SLOW_ABOVE else []
        yield pytest.param(case, id=case["request"]["id"], marks=marks)


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _check(result, expected, tmp_path, want_triplets, rc1=None):
    keys = tmp_path / native.keys_filename("reference")
    native.write_keys(result, str(keys))
    assert _sha(keys) == expected["keys_sha256"]
    assert result.counts.n_points == expected["n_points"]
    assert result.counts.n_triangles == expected["n_triangles"]
    assert result.counts.n_unique_keys == expected["n_unique_keys"]
    assert result.key_space.fingerprint == expected["key_space_fingerprint"]
    assert result.sampling.fingerprint == expected["sampling_fingerprint"]
    assert list(result.selection.residues_used) == expected["residues_used"]
    assert dict(result.category_totals()) == expected["category_totals"]
    if want_triplets and "triplets_sha256_macos_arm64" in expected:
        triplets = tmp_path / native.triplets_filename("reference")
        native.write_triplets(result, str(triplets))
        assert _sha(triplets) == expected["triplets_sha256_macos_arm64"]
        if rc1 is not None:
            assert _respelled_sha(triplets, expected.get("rc1_label_names", {})) == \
                rc1["triplets_sha256_macos_arm64"]


def _respelled_sha(path, names):
    """SHA-256 of a triplet file with each label name spelled the rc1 way."""
    with open(str(path), encoding=native.ENCODING) as fh:
        header = fh.readline()
        columns = header.rstrip("\n").split("\t")
        at = [columns.index("label_name%d" % i) for i in (1, 2, 3)]
        lines = [header]
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            for i in at:
                fields[i] = names.get(fields[i], fields[i])
            lines.append("\t".join(fields) + "\n")
    return hashlib.sha256("".join(lines).encode(native.ENCODING)).hexdigest()


def _to_rc1_radix(key, radix):
    """Decode a radix-30 key and re-encode the same labels and bins at ``radix``."""
    sign = -1 if key < 0 else 1
    magnitude = abs(key)
    bins = magnitude % _BIN
    component = magnitude // _BIN - (_LOWEST_CODE - 1) * (30 * 30 + 30 + 1)
    e1, rest = divmod(component, 30 * 30)
    e2, e3 = divmod(rest, 30)
    labels = [e + _LOWEST_CODE for e in (e1, e2, e3)]
    assert all(_LOWEST_CODE <= c < _LOWEST_CODE + radix for c in labels), labels
    c1, c2, c3 = labels
    return sign * (_BIN * (radix * radix * (c1 - 1) + radix * (c2 - 1) + (c3 - 1)) + bins)


def in_previous_selection(triplet):
    """Whether every vertex of a C-alpha triangle is one the previous selection
    read: a regular amino acid, by residue identity.

    Never by label: MSE shares MET's code 13, which lies inside the old radix-20
    window, so a label test would map an MSE triangle onto an old MET key.
    """
    return all(i.residue_name in STANDARD_AMINO_ACIDS for i in triplet.identities)


def _keys_sha(rows):
    """SHA-256 of the keys file the ``(key, category) -> count`` rows give."""
    text = native.keys_header() + "".join(
        "%d\t%s\t%d\n" % (k, c, rows[(k, c)]) for k, c in sorted(rows))
    return hashlib.sha256(text.encode(native.ENCODING)).hexdigest()


def _rc1_keys_sha(rows, radix):
    """SHA-256 of the keys file the ``(key, category) -> count`` rows give at
    the rc1 radix."""
    remapped = {}
    for (key, category), count in rows.items():
        old = (_to_rc1_radix(key, radix), category)
        assert old not in remapped, "two radix-30 keys remapped onto one rc1 key"
        remapped[old] = count
    text = native.keys_header() + "".join(
        "%d\t%s\t%d\n" % (k, c, remapped[(k, c)]) for k, c in sorted(remapped))
    return hashlib.sha256(text.encode(native.ENCODING)).hexdigest()


@pytest.mark.parametrize("case", list(_params()))
def test_the_standard_result_is_unchanged(case, tmp_path):
    request = dict(case["request"])
    case_id = request.pop("id")
    api = request.pop("api")
    pdb = os.path.join(STRUCTURES, request.pop("pdb") + ".pdb")
    rc1 = case["expected"]
    calpha = _CALPHA.get(case_id)
    atomic = _ATOMIC.get(case_id)
    expected = calpha["expected"] if calpha else atomic["expected"]
    # Triplet bytes are compared on the reference platform for the smaller
    # cases; a keys-only run is the same calculation and is much lighter.
    want_triplets = _ON_REFERENCE_PLATFORM and expected.get("n_triangles", 0) <= 120_000
    # C-alpha: the triangles of the previous selection, by key and category,
    # for the remap to the rc1 file; the others are counted.
    standard_rows, added = {}, [0]

    def tally(t):
        if not in_previous_selection(t):
            added[0] += 1
            return
        row = (t.key, t.category)
        standard_rows[row] = standard_rows.get(row, 0) + 1

    def tally_nucleotides(t):
        if any(i.residue_name == "DC" and i.residue.legacy_seq_label == "1"
               for i in t.identities):
            added[0] += 1
            return
        row = (t.key, t.category)
        standard_rows[row] = standard_rows.get(row, 0) + 1

    sink = tally if calpha else tally_nucleotides if case_id == _NUCLEOTIDE_CHAIN else None
    if api in ("analyze", "analyze_cross"):
        run = biotsr.analyze if api == "analyze" else biotsr.analyze_cross
        result = run(pdb, keep_triplets=want_triplets, triplet_sink=sink, **request)
        _check(result, expected, tmp_path, want_triplets, rc1=None if calpha else rc1)
        assert result.key_space.fingerprint != rc1["key_space_fingerprint"] or not atomic
        if calpha:
            assert _rc1_keys_sha(standard_rows, calpha["rc1_radix"]) == rc1["keys_sha256"]
            assert added[0] == calpha.get("added_triangles", 0)
            assert result.sampling.fingerprint == rc1["sampling_fingerprint"]
            assert result.key_space.fingerprint != rc1["key_space_fingerprint"]
        if case_id == _NUCLEOTIDE_CHAIN:
            # Atomic labels are the same codes as rc1's for these atoms, so the
            # old selection's triangles give rc1's keys file unchanged.
            assert _keys_sha(standard_rows) == rc1["keys_sha256"]
            assert added[0] == expected["added_triangles"] == 438216
    else:
        results = biotsr.analyze_entities(pdb, keep_triplets=False, **request)
        assert sorted(results) == expected["entities"]
        for name, result in results.items():
            path = tmp_path / ("%s.keys" % name)
            native.write_keys(result, str(path))
            assert _sha(path) == expected["entity_keys_sha256"][name], name
            assert result.key_space.fingerprint == expected["key_space_fingerprint"]
            assert result.sampling.fingerprint == expected["sampling_fingerprint"]
        assert sum(r.counts.n_triangles for r in results.values()) == \
            expected["n_triangles_total"]


def test_the_rc1_reference_file_is_unchanged():
    with open(_BASELINE, "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == _BASELINE_SHA256


def test_every_calpha_case_has_a_base30_expectation_and_nothing_else_does():
    calpha_methods = {"protein", "aa_selection", "protein_protein_calpha_cross"}
    ids = {c["request"]["id"] for c in _CASES
           if c["request"]["method"] in calpha_methods}
    assert set(_CALPHA) == ids and len(ids) == 17
    for case_id, case in _CALPHA.items():
        assert case["request"] == next(c["request"] for c in _CASES
                                       if c["request"]["id"] == case_id)


def test_every_atomic_case_has_a_current_expectation_that_differs_from_rc1_only_as_approved():
    """The reference file itself, read against rc1: nothing in it moved that
    ADR-113 and ADR-114 do not account for."""
    calpha_methods = {"protein", "aa_selection", "protein_protein_calpha_cross"}
    rc1 = {c["request"]["id"]: c for c in _CASES}
    assert set(_ATOMIC) == {i for i, c in rc1.items()
                            if c["request"]["method"] not in calpha_methods}
    assert len(_ATOMIC) == 25
    for case_id, case in _ATOMIC.items():
        assert case["request"] == rc1[case_id]["request"]
        old, new = rc1[case_id]["expected"], case["expected"]
        changed = {k for k in set(old) | set(new) if old.get(k) != new.get(k)}
        changed -= {"rc1_label_names", "added_triangles", "rc1_selection_excludes"}
        allowed = {"key_space_fingerprint", "triplets_sha256_macos_arm64"}
        if case_id == _NUCLEOTIDE_CHAIN:
            allowed |= {"keys_sha256", "n_points", "n_triangles", "n_unique_keys",
                        "residues_used", "category_totals", "sampling_fingerprint"}
            assert new["n_points"] == old["n_points"] + 16 == 243
            assert new["residues_used"] == ["1_DC"] + old["residues_used"]
        assert "key_space_fingerprint" in changed, case_id
        assert changed <= allowed, (case_id, changed - allowed)
        if "triplets_sha256_macos_arm64" in changed and "triplets_sha256_macos_arm64" in new:
            assert new["rc1_label_names"], case_id
    # One key space for every atomic case without a modifier that changes it.
    plain = {c["expected"]["key_space_fingerprint"] for i, c in _ATOMIC.items()
             if not c["request"].get("mirror")}
    assert len(plain) == 1


def test_every_method_is_covered():
    methods = {c["request"]["method"] for c in _CASES}
    assert methods == set(biotsr.list_methods())


def test_every_supported_modifier_is_covered_somewhere():
    seen = set()
    for case in _CASES:
        request = case["request"]
        recipe = biotsr.get_method(request["method"])
        seen.update((recipe.name, m) for m in recipe.supported_modifiers
                    if request.get(m))
    families = {(biotsr.get_method(n).family, m) for n, m in seen}
    assert ("protein", "aa_grouping") in families or ("calpha", "aa_grouping") in \
        families or any(m == "aa_grouping" for _, m in families)
    for modifier in ("mirror", "size_filter", "sse", "aa_grouping"):
        assert any(m == modifier for _, m in seen), modifier


def test_previous_selection_membership_is_decided_by_residue_identity(tmp_path):
    """MSE carries MET's label 13, inside the old radix-20 window, so a label
    test cannot tell an MSE triangle from a MET one: the window check in
    ``_to_rc1_radix`` accepts it.  Membership by residue identity excludes it.

    1UBQ's MET 1 rewritten as a HETATM MSE: every triangle containing it is
    new, and the remaining ones are exactly those of the structure without it.
    """
    source = os.path.join(STRUCTURES, "1UBQ.pdb")
    with open(source) as fh:
        lines = fh.readlines()
    mse = tmp_path / "MSE1.pdb"
    without = tmp_path / "NO1.pdb"
    with open(str(mse), "w") as a, open(str(without), "w") as b:
        for line in lines:
            first = line.startswith("ATOM") and line[21] == "A" and line[22:26].strip() == "1"
            if first:
                a.write("HETATM" + line[6:17] + "MSE" + line[20:])
                continue
            a.write(line)
            b.write(line)
    triplets = biotsr.analyze(str(mse), method="protein", chain="A").triplets
    assert len(triplets) == 70300
    by_identity = [t for t in triplets if in_previous_selection(t)]
    by_label = [t for t in triplets if all(4 <= c <= 23 for c in t.labels)]
    assert len(triplets) - len(by_identity) == 2775          # C(75, 2), every MSE 1 triangle
    assert len(by_label) == 70300                             # the label test keeps them all
    mse_triangles = [t for t in by_label if not in_previous_selection(t)]
    _to_rc1_radix(mse_triangles[0].key, 20)                   # ... and the window check passes them
    reference = biotsr.analyze(str(without), method="protein", chain="A", keep_triplets=False)
    remapped = {}
    for t in by_identity:
        remapped[_to_rc1_radix(t.key, 20)] = remapped.get(_to_rc1_radix(t.key, 20), 0) + 1
    assert remapped == {_to_rc1_radix(k, 20): n for k, n in reference.frequencies.items()}
