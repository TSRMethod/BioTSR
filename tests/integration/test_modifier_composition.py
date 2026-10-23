"""Modifier composition: all eight combinations, one engine.

The claim under test is architectural, not numeric: the package needs **no
combinatorial subclasses**.  Every combination is the same selector, the same
candidate strategy, the same canonical ordering, the same geometry, the same
binning, the same encoder and the same engine loop, reading different
configuration.

The eight combinations are exercised end to end, their identities checked for
correct separation, and -- in ``tests/characterization`` -- each is compared
byte-for-byte against a legacy reference produced by the historical
implementation with the matching flags.
"""

from __future__ import annotations

import itertools
import os

import pytest

import biotsr
from biotsr.domain.result import DEFAULT_CATEGORY
from biotsr.modifiers.mirror import MIRROR_CATEGORY

from ..conftest import FIXTURE_PDB, _require

CUTOFF = 25.0

#: (label, aa_grouping, mirror, size_filter)
COMBINATIONS = [
    ("regular", False, False, None),
    ("grouping", True, False, None),
    ("mirror", False, True, None),
    ("size_filter", False, False, CUTOFF),
    ("grouping+mirror", True, True, None),
    ("grouping+size_filter", True, False, CUTOFF),
    ("mirror+size_filter", False, True, CUTOFF),
    ("grouping+mirror+size_filter", True, True, CUTOFF),
]


@pytest.fixture(scope="module")
def ubq_path():
    path = os.path.join(FIXTURE_PDB, "1UBQ.pdb")
    _require(path)
    return path


@pytest.fixture(scope="module")
def results(ubq_path):
    out = {}
    for label, grouping, mirror, size_filter in COMBINATIONS:
        out[label] = biotsr.analyze(
            ubq_path,
            method="protein",
            chain="A",
            aa_grouping=grouping,
            mirror=mirror,
            size_filter=size_filter,
        )
    return out


def test_all_eight_combinations_are_covered():
    assert len(COMBINATIONS) == 8
    flags = {(g, m, s is not None) for _, g, m, s in COMBINATIONS}
    assert flags == set(itertools.product([False, True], repeat=3))


# ---------------------------------------------------------------------------
# Every combination runs, on the same pipeline
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label,grouping,mirror,size_filter", COMBINATIONS)
def test_each_combination_produces_a_coherent_result(
    results, label, grouping, mirror, size_filter
):
    r = results[label]
    assert r.counts.n_candidates == 70300
    assert r.counts.n_degenerate == 0
    assert r.counts.ledger_balances
    assert r.counts.n_triangles == sum(r.frequencies.values())
    assert r.counts.n_unique_keys == len(r.keys)
    assert sum(r.category_totals().values()) == r.counts.n_triangles


@pytest.mark.parametrize("label,grouping,mirror,size_filter", COMBINATIONS)
def test_every_combination_reports_the_same_provenance_shape(
    results, label, grouping, mirror, size_filter
):
    r = results[label]
    expected = []
    if grouping:
        expected.append("aa_grouping")
    if mirror:
        expected.append("mirror")
    if size_filter is not None:
        expected.append("size_filter")
    assert list(r.annotations["modifiers"]) == expected
    assert r.annotations["base_profile"] == "standard_protein_calpha_v1"
    assert r.method == "protein"
    assert r.selection.selector == "protein_calpha"


def test_no_combinatorial_classes_exist_in_the_package():
    """The point of the whole exercise, asserted mechanically.

    Checked against defined class and function *names* in the AST, not raw text:
    the modifier docstrings deliberately name the classes that do not exist
    ("no ``GroupedMirrorTSR``"), and that prose is the documentation of the
    decision rather than a violation of it.
    """
    import ast

    src = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "src",
        "biotsr",
    )
    defined = []
    for root, _dirs, files in os.walk(src):
        for name in sorted(files):
            if not name.endswith(".py"):
                continue
            tree = ast.parse(open(os.path.join(root, name)).read())
            for node in ast.walk(tree):
                if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
                    defined.append(node.name)

    lowered = [d.lower() for d in defined]
    # No name may mention two modifiers at once -- that is what a combinatorial
    # class looks like.
    tokens = ("group", "mirror", "filter")
    offenders = [
        d
        for d, low in zip(defined, lowered)
        if sum(1 for t in tokens if t in low) > 1
    ]
    assert offenders == [], "combinatorial modifier definitions found: %r" % offenders
    assert "MirrorTSR" not in defined
    assert "GroupedTSR" not in defined


def test_one_registered_method_serves_every_combination(results):
    """All eight combinations run under the single `protein` recipe.

    Other methods exist (atomic ones among them), but none of them is a
    per-combination variant of this one -- that is the claim being made.
    """
    assert {r.method for r in results.values()} == {"protein"}
    assert "protein" in biotsr.list_methods()


# ---------------------------------------------------------------------------
# What each modifier does, and does not, change
# ---------------------------------------------------------------------------


def test_size_filter_is_the_only_modifier_that_changes_the_triangle_count(results):
    unfiltered = {l for l, _g, _m, s in COMBINATIONS if s is None}
    filtered = {l for l, _g, _m, s in COMBINATIONS if s is not None}
    assert {results[l].counts.n_triangles for l in unfiltered} == {70300}
    assert {results[l].counts.n_triangles for l in filtered} == {60742}
    assert {results[l].counts.n_filtered_by_size for l in unfiltered} == {0}
    assert {results[l].counts.n_filtered_by_size for l in filtered} == {9558}


def test_filtered_triangles_are_never_reported_as_degenerate(results):
    for label, _g, _m, size_filter in COMBINATIONS:
        r = results[label]
        assert r.counts.n_degenerate == 0
        assert r.counts.n_valid_geometry == 70300
        if size_filter is not None:
            assert r.counts.n_filtered_by_size > 0
        reasons = {s.reason for s in r.skipped}
        assert "size_filtered" not in reasons


def test_mirror_is_the_only_modifier_that_produces_negative_keys(results):
    for label, _g, mirror, _s in COMBINATIONS:
        negatives = [k for k in results[label].keys if k < 0]
        if mirror:
            assert negatives, "%s should have negative keys" % label
        else:
            assert negatives == [], "%s must never produce a negative key" % label


def test_mirror_is_the_only_modifier_that_produces_a_mirror_category(results):
    for label, _g, mirror, _s in COMBINATIONS:
        categories = set(results[label].category_totals())
        if mirror:
            assert categories == {DEFAULT_CATEGORY, MIRROR_CATEGORY}
        else:
            assert categories == {DEFAULT_CATEGORY}


def test_grouping_is_the_only_modifier_that_changes_the_label_scheme(results):
    """Both residue schemes are radix 30 (ADR-107): grouping changes the scheme,
    never the radix, and no other modifier changes either."""
    for label, grouping, _m, _s in COMBINATIONS:
        scheme = results[label].key_space.label_scheme
        assert scheme.radix == 30, label
        assert scheme.name == ("amino_acid_group_23" if grouping else "amino_acid_30"), label


def test_grouping_reduces_the_distinct_key_count(results):
    """Fewer labels means fewer distinguishable triangles -- a sanity check that
    grouping actually did something."""
    assert results["grouping"].counts.n_unique_keys < results["regular"].counts.n_unique_keys
    assert (
        results["grouping+mirror"].counts.n_unique_keys
        < results["mirror"].counts.n_unique_keys
    )


def test_grouping_and_mirror_are_not_independent(results):
    """A real interaction, worth stating rather than assuming away.

    Grouping changes labels, labels drive canonical ordering, ordering decides
    which vertex is the apex, and the apex decides theta_raw -- which is what
    mirror keys its sign off.  So the mirrored fraction legitimately differs
    between grouped and ungrouped runs.  Composition is not the product of the
    parts, and the engine never assumes it is.
    """
    plain_mirrored = results["mirror"].category_totals()[MIRROR_CATEGORY]
    grouped_mirrored = results["grouping+mirror"].category_totals()[MIRROR_CATEGORY]
    assert plain_mirrored != grouped_mirrored


# ---------------------------------------------------------------------------
# Identity separation
# ---------------------------------------------------------------------------


def test_the_eight_combinations_occupy_eight_distinct_identity_pairs(results):
    pairs = {
        label: (r.key_space.fingerprint, r.sampling.fingerprint)
        for label, r in results.items()
    }
    assert len(set(pairs.values())) == 8


def test_key_changing_modifiers_change_the_key_space(results):
    base = results["regular"].key_space.fingerprint
    assert results["grouping"].key_space.fingerprint != base
    assert results["mirror"].key_space.fingerprint != base
    assert results["size_filter"].key_space.fingerprint == base


def test_sampling_changing_modifiers_change_the_sampling_profile(results):
    base = results["regular"].sampling.fingerprint
    assert results["size_filter"].sampling.fingerprint != base
    assert results["grouping"].sampling.fingerprint == base
    assert results["mirror"].sampling.fingerprint == base


def test_incompatible_results_refuse_to_be_compared(results):
    from biotsr.domain.keyspace import IncompatibleKeySpace, IncompatibleSampling

    with pytest.raises(IncompatibleKeySpace):
        results["regular"].key_space.require_compatible(results["grouping"].key_space)
    with pytest.raises(IncompatibleSampling):
        results["regular"].sampling.require_compatible(results["size_filter"].sampling)


# ---------------------------------------------------------------------------
# Determinism, under every combination
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label,grouping,mirror,size_filter", COMBINATIONS)
def test_repeated_runs_are_byte_identical(
    results, ubq_path, label, grouping, mirror, size_filter
):
    """Identical input and configuration must give identical bytes.

    No modifier may make the result depend on dict or set traversal order.
    """
    from biotsr.io import native

    def serialize(r):
        keys = "".join(
            "%d\t%s\t%d\n" % (k, c, n)
            for k in sorted(r.keys)
            for c, n in sorted(r.keys[k].by_category.items())
        )
        triplets = "".join(native.format_native_triplet(t) for t in r.triplets)
        return keys, triplets, tuple(r.key_order)

    fresh = biotsr.analyze(
        ubq_path,
        method="protein",
        chain="A",
        aa_grouping=grouping,
        mirror=mirror,
        size_filter=size_filter,
    )
    assert serialize(results[label]) == serialize(fresh)


def test_modifier_order_of_specification_does_not_matter(ubq_path):
    """Keyword order is not a pipeline order -- the hooks are fixed."""
    a = biotsr.analyze(
        ubq_path, method="protein", chain="A",
        aa_grouping=True, mirror=True, size_filter=CUTOFF,
    )
    b = biotsr.analyze(
        ubq_path, method="protein", chain="A",
        size_filter=CUTOFF, mirror=True, aa_grouping=True,
    )
    assert a.frequencies == b.frequencies
    assert a.key_order == b.key_order
    assert a.key_space.fingerprint == b.key_space.fingerprint


# ---------------------------------------------------------------------------
# Streaming still works under every combination
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label,grouping,mirror,size_filter", COMBINATIONS)
def test_streaming_matches_the_retained_run(
    ubq_path, label, grouping, mirror, size_filter
):
    """Size filtering must reject BEFORE anything downstream does work, so a
    filtered candidate never reaches the sink."""
    seen = []
    streamed = biotsr.analyze(
        ubq_path,
        method="protein",
        chain="A",
        aa_grouping=grouping,
        mirror=mirror,
        size_filter=size_filter,
        keep_triplets=False,
        triplet_sink=seen.append,
    )
    assert streamed.triplets is None
    assert len(seen) == streamed.counts.n_triangles
    if size_filter is not None:
        assert all(t.max_distance < size_filter for t in seen)
