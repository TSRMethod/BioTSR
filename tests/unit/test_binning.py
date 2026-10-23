"""Bin-scheme boundary behaviour for the standard theta and Ca distance schemes."""

from __future__ import annotations

import pytest

from biotsr.binning.base import LadderBinScheme
from biotsr.binning.distance import (
    CALPHA_BOUNDARIES,
    DISTANCE_SCHEME_CALPHA_35,
)
from biotsr.binning.theta import THETA_BOUNDARIES_29, THETA_SCHEME_29

EPS = 1e-9


# ---------------------------------------------------------------------------
# Theta -- 29 classes
# ---------------------------------------------------------------------------


def test_theta_scheme_shape():
    assert THETA_SCHEME_29.n_classes == 29
    assert len(THETA_BOUNDARIES_29) == 29
    assert THETA_SCHEME_29.max_emitted_class == 29


@pytest.mark.parametrize("index", range(len(THETA_BOUNDARIES_29)))
def test_theta_boundaries_are_lower_inclusive_upper_exclusive(index):
    cut = THETA_BOUNDARIES_29[index]
    assert THETA_SCHEME_29.classify(cut - EPS) == index + 1
    expected_at_cut = index + 2 if index + 1 < len(THETA_BOUNDARIES_29) else 29
    assert THETA_SCHEME_29.classify(cut) == expected_at_cut


def test_theta_zero_is_class_one():
    assert THETA_SCHEME_29.classify(0.0) == 1


def test_theta_at_exactly_90_is_a_total_function_returning_29():
    """The legacy if/elif ladder raises UnboundLocalError here; this one does not.

    29 matches the historical TSR.py generator form.
    """
    assert THETA_SCHEME_29.classify(90.0) == 29
    assert THETA_SCHEME_29.classify(89.999999) == 29
    assert THETA_SCHEME_29.classify(1e6) == 29


def test_the_anomalously_wide_bin_is_preserved_not_corrected():
    """72.69-79.2 is unusually wide compared with its neighbours, identically so in
    every legacy repo.  It is an open domain question, not a transcription bug to
    silently fix -- so the table is reproduced verbatim."""
    i = THETA_BOUNDARIES_29.index(72.69)
    assert THETA_BOUNDARIES_29[i + 1] == 79.2
    widths = [
        THETA_BOUNDARIES_29[n + 1] - THETA_BOUNDARIES_29[n]
        for n in range(len(THETA_BOUNDARIES_29) - 1)
    ]
    assert widths[i] == max(widths), "the 72.69-79.2 bin is the widest in the table"
    assert widths[i] > 1.2 * sorted(widths)[-2], "and by a clear margin"
    assert widths[i] > 2.9 * widths[i - 1], "against its immediate neighbour, ~3x"


# ---------------------------------------------------------------------------
# Ca distance -- the standard 35 classes
# ---------------------------------------------------------------------------


def test_the_standard_calpha_scheme_has_thirty_four_cut_points():
    assert DISTANCE_SCHEME_CALPHA_35.boundaries == CALPHA_BOUNDARIES
    assert len(CALPHA_BOUNDARIES) == 34


@pytest.mark.parametrize(
    "value,expected",
    [
        (0.0, 1),
        (3.82999, 1),
        (3.83, 2),
        (83.25999, 33),
        (83.26, 34),
        (132.44, 34),
        (132.449999, 34),
        (132.45, 35),
        (200.0, 35),
        (1e9, 35),
    ],
)
def test_standard_calpha_scheme(value, expected):
    assert DISTANCE_SCHEME_CALPHA_35.classify(value) == expected


def test_standard_scheme_has_no_unreachable_class():
    assert DISTANCE_SCHEME_CALPHA_35.n_classes == 35
    assert DISTANCE_SCHEME_CALPHA_35.max_emitted_class == 35
    assert not DISTANCE_SCHEME_CALPHA_35.allow_unreachable_classes


# ---------------------------------------------------------------------------
# Fingerprints
# ---------------------------------------------------------------------------


def test_identical_boundaries_but_different_overflow_fingerprint_differently():
    """The decisive hazard: a boundary-array-only hash would call these equal."""
    merged_top = LadderBinScheme(
        name="distance_calpha_35",
        boundaries=CALPHA_BOUNDARIES,
        n_classes=35,
        overflow_class=34,
        allow_unreachable_classes=True,
        unreachable_reason="test fixture: the top two classes merged",
    )
    assert merged_top.boundaries == DISTANCE_SCHEME_CALPHA_35.boundaries
    assert merged_top.classify(132.45) == 34
    assert DISTANCE_SCHEME_CALPHA_35.classify(132.45) == 35
    assert merged_top.content_fingerprint() != DISTANCE_SCHEME_CALPHA_35.content_fingerprint()


def test_fingerprint_ignores_the_name():
    renamed = LadderBinScheme(
        name="something_else_entirely",
        boundaries=CALPHA_BOUNDARIES,
        n_classes=35,
        overflow_class=35,
    )
    assert renamed.content_fingerprint() == DISTANCE_SCHEME_CALPHA_35.content_fingerprint()


def test_fingerprint_changes_when_a_boundary_moves():
    moved = list(CALPHA_BOUNDARIES)
    moved[0] = 3.84
    other = LadderBinScheme(
        name="distance_calpha_35",
        boundaries=tuple(moved),
        n_classes=35,
        overflow_class=35,
    )
    assert other.content_fingerprint() != DISTANCE_SCHEME_CALPHA_35.content_fingerprint()


# ---------------------------------------------------------------------------
# Construction validation
# ---------------------------------------------------------------------------


def test_unreachable_class_is_rejected_without_an_explicit_opt_out():
    with pytest.raises(ValueError, match="can only emit"):
        LadderBinScheme(
            name="oops", boundaries=(1.0, 2.0), n_classes=5, overflow_class=2
        )


def test_opt_out_requires_a_reason():
    with pytest.raises(ValueError, match="unreachable_reason"):
        LadderBinScheme(
            name="oops",
            boundaries=(1.0, 2.0),
            n_classes=5,
            overflow_class=2,
            allow_unreachable_classes=True,
        )


def test_overflow_class_cannot_exceed_the_declared_radix():
    with pytest.raises(ValueError, match="exceeds declared"):
        LadderBinScheme(
            name="oops", boundaries=(1.0, 2.0), n_classes=2, overflow_class=3
        )


def test_boundaries_must_be_strictly_increasing():
    with pytest.raises(ValueError, match="strictly increasing"):
        LadderBinScheme(
            name="oops", boundaries=(2.0, 1.0), n_classes=2, overflow_class=2
        )
