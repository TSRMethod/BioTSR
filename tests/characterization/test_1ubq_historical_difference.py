"""1UBQ chain A: the standard result, and exactly how it differs from TSR.py's.

BioTSR does not reproduce historical programs.  It does state, and test, how its
standard protein C-alpha result differs from the historical program it replaces,
so the difference is a documented scientific decision rather than drift.

The comparison data is one column -- the integer key of each triangle, in
candidate order -- extracted from the output of ``TSR-Package/tsr_package/TSR.py``
(commit ``fd3c425a``) on the same 1UBQ file (see
``tests/data/reference/README.md``).  Nothing of that program is run here.

Two standardisations apply.  The label radix: ``TSR.py`` encodes its 20 residue
codes at radix 20, the standard encodes the residue labels -- which add the
supported modified residues -- at radix 30 (ADR-107).  1UBQ has no modified
residue, so the labels themselves are identical and the radix only renames each
key; the comparison is therefore made on the historical radix, by re-encoding
each standard triangle's labels and classes at 20.  What then remains is the
angle convention: ``TSR.py`` divides by the literal ``3.14`` where the standard
uses ``math.pi`` (ADR-002).  1UBQ's 76-residue C-alpha span never reaches
132.45 A, so the 35-class distance-scale correction (ADR-014) cannot change a
key here.
"""

from __future__ import annotations

import gzip
import hashlib
import math
import os

import pytest

import biotsr
from biotsr.binning.theta import THETA_SCHEME_29
from biotsr.encoding.key import KeyFormulaSpec, encode
from biotsr.io import native

from ..conftest import REFERENCE_DATA, _require

HISTORICAL_KEYS = os.path.join(
    REFERENCE_DATA, "historical", "1UBQ_A_TSR-Package_TSR_keys.txt.gz"
)
#: SHA-256 of the decompressed key column (one key per line).
HISTORICAL_KEYS_SHA256 = "4a33134c42dce5e38da8aff585ded1dec20f37157e8cec83ac111afc69ef676f"

#: The historical program's angle divisor. It appears in this test, which
#: describes the difference, and nowhere in the package.
HISTORICAL_DIVISOR = 3.14

#: The historical program's key formula: ``numOfLabels = 20``.
HISTORICAL_FORMULA = KeyFormulaSpec(n_labels=20, n_theta=29, n_distance=35)


def at_historical_radix(t):
    """A standard triangle's key, written at the historical radix 20."""
    return encode(t.labels, t.theta_class, t.distance_class, HISTORICAL_FORMULA)


@pytest.fixture(scope="module")
def standard(ubq_pdb):
    return biotsr.analyze(ubq_pdb, method="protein", chain="A")


@pytest.fixture(scope="module")
def historical_keys():
    with gzip.open(_require(HISTORICAL_KEYS), "rb") as fh:
        data = fh.read()
    assert hashlib.sha256(data).hexdigest() == HISTORICAL_KEYS_SHA256
    return [int(line) for line in data.decode("ascii").split()]


def test_the_standard_result_has_the_expected_shape(standard):
    assert standard.counts.n_points == 76
    assert standard.counts.n_triangles == math.comb(76, 3) == 70300
    assert standard.counts.n_unique_keys == 45349
    assert not standard.skipped
    assert not standard.warnings
    assert all(set(c.by_category) == {"regular"} for c in standard.keys.values())


def test_the_historical_column_is_complete(historical_keys):
    assert len(historical_keys) == 70300


def test_the_standard_keys_are_written_at_radix_30(standard):
    """The radix renames every key (ADR-107): the same labels and classes."""
    spec = biotsr.STANDARD_PROTEIN_CALPHA_V1.key_formula()
    assert spec.n_labels == 30
    for t in standard.triplets:
        assert t.key == encode(t.labels, t.theta_class, t.distance_class, spec)
        assert max(t.labels) <= 23, "1UBQ has only the 20 regular residues"


def test_the_difference_is_exactly_the_documented_one(standard, historical_keys):
    """At the historical radix, 1,144 of 70,300 keys change; 45,353 historical
    distinct keys become 45,349, with 529 appearing and 533 disappearing."""
    changed = sum(
        1 for t, h in zip(standard.triplets, historical_keys)
        if at_historical_radix(t) != h
    )
    assert changed == 1144
    historical = set(historical_keys)
    current = {at_historical_radix(t) for t in standard.triplets}
    assert len(set(standard.frequencies)) == len(current)
    assert len(historical) == 45353
    assert len(current) == 45349
    assert len(current - historical) == 529
    assert len(historical - current) == 533


def test_every_change_is_in_the_theta_class_alone(standard, historical_keys):
    """Labels and distance class never differ; where a key differs, only the
    theta class does."""
    for t, h in zip(standard.triplets, historical_keys):
        theta_class = h % 29 + 1
        distance_class = (h // 29) % 35 + 1
        if at_historical_radix(t) == h:
            assert (t.theta_class, t.distance_class) == (theta_class, distance_class)
            continue
        assert t.distance_class == distance_class
        assert t.theta_class != theta_class
        assert encode(t.labels, theta_class, distance_class, HISTORICAL_FORMULA) == h


def test_the_historical_theta_classes_follow_from_the_divisor_alone(standard,
                                                                    historical_keys):
    """Rescaling each standard pre-fold angle by pi/3.14, folding as the historical
    program did and binning on the unchanged 29-class scheme reproduces every
    historical theta class -- so the angle convention explains the whole
    difference."""
    mismatched = []
    for t, h in zip(standard.triplets, historical_keys):
        historical_raw = t.theta_raw * math.pi / HISTORICAL_DIVISOR
        folded = (historical_raw if historical_raw <= 90
                  else abs(180 - historical_raw))
        if THETA_SCHEME_29.classify(folded) != h % 29 + 1:
            mismatched.append((t.theta_raw, folded))
    # A last-bit rounding difference could only matter at a bin boundary; none
    # occurs on this structure.
    assert mismatched == []


def test_the_largest_angle_shift_matches_the_predicted_magnitude(standard):
    shift = max(
        abs(t.theta_raw * math.pi / HISTORICAL_DIVISOR - t.theta_raw)
        for t in standard.triplets
    )
    assert shift == pytest.approx(0.0909, abs=5e-4)


def test_the_distance_correction_cannot_be_exercised_by_this_structure(standard):
    assert max(t.max_distance for t in standard.triplets) < 132.45
    assert biotsr.STANDARD_PROTEIN_CALPHA_V1.distance_scheme.classify(132.45) == 35


def test_the_standard_file_names(standard):
    prefix = native.default_prefix(standard)
    assert prefix == "1UBQ_A_protein"
    assert native.keys_filename(prefix) == "1UBQ_A_protein.keys"
    assert native.triplets_filename(prefix) == "1UBQ_A_protein.triplets"
