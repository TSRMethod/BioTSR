"""Mixed-side cross candidates: exact counts, exact composition, exact order.

A test that only checks "some candidates were
produced" is not adequate here: the whole scientific content of cross candidate
generation is *which* triples exist and *how many*, so every assertion below is
an exact count or an exact set.
"""

from __future__ import annotations

from itertools import combinations
from math import comb

import pytest

from biotsr.cross.candidates import (
    COMPOSITION_RULE,
    count,
    count_from_sides,
    count_via_complement,
    mixed_side,
    side_index_lists,
)
from biotsr.domain.candidate import TriangleCandidate


def sides(a, b):
    """``a`` points of side A followed by ``b`` points of side B."""
    return ["A"] * a + ["B"] * b


# ===========================================================================
# The counting identity
# ===========================================================================


@pytest.mark.parametrize("a,b", [(1, 2), (2, 1), (2, 2), (3, 4), (4, 3), (0, 5),
                                 (5, 0), (1, 1), (0, 0), (10, 1), (7, 9)])
def test_the_two_spellings_of_the_mixed_side_count_agree(a, b):
    """``a*C(b,2) + b*C(a,2)  ==  C(a+b,3) - C(a,3) - C(b,3)``."""
    assert count(a, b) == count_via_complement(a, b)


@pytest.mark.parametrize("a,b", [(1, 2), (2, 1), (2, 2), (3, 4), (6, 5)])
def test_the_count_matches_the_closed_form_written_out(a, b):
    assert count(a, b) == a * comb(b, 2) + b * comb(a, 2)


@pytest.mark.parametrize("a,b", [(1, 2), (2, 1), (2, 2), (3, 4), (5, 6), (1, 1),
                                 (0, 4), (4, 0), (8, 3)])
def test_the_generator_emits_exactly_the_counted_number(a, b):
    tags = sides(a, b)
    assert len(list(mixed_side(tags))) == count(a, b)
    assert count_from_sides(tags) == count(a, b)


@pytest.mark.parametrize("a,b", [(1, 2), (2, 2), (3, 4), (4, 5)])
def test_the_generator_emits_exactly_the_non_monochromatic_triples(a, b):
    """Set equality against the brute-force definition, not just the count."""
    tags = sides(a, b)
    n = a + b
    expected = {
        (i, j, k)
        for i, j, k in combinations(range(n), 3)
        if len({tags[i], tags[j], tags[k]}) == 2
    }
    assert {c.as_tuple() for c in mixed_side(tags)} == expected


# ===========================================================================
# Composition: never AAA, never BBB
# ===========================================================================


@pytest.mark.parametrize("a,b", [(3, 3), (4, 2), (2, 4), (5, 5), (6, 1)])
def test_no_monochromatic_candidate_is_ever_emitted(a, b):
    tags = sides(a, b)
    for candidate in mixed_side(tags):
        composition = tuple(
            sorted(tags[i] for i in (candidate.i, candidate.j, candidate.k))
        )
        assert composition in (("A", "A", "B"), ("A", "B", "B")), composition


@pytest.mark.parametrize("a,b", [(3, 3), (4, 5)])
def test_both_allowed_compositions_are_actually_produced(a, b):
    tags = sides(a, b)
    seen = {
        tuple(sorted(tags[i] for i in c.as_tuple())) for c in mixed_side(tags)
    }
    assert seen == {("A", "A", "B"), ("A", "B", "B")}


@pytest.mark.parametrize("a,b", [(3, 3), (4, 2)])
def test_the_counts_of_the_two_compositions_are_exactly_aC2b_and_bC2a(a, b):
    tags = sides(a, b)
    aab = abb = 0
    for c in mixed_side(tags):
        if sorted(tags[i] for i in c.as_tuple()) == ["A", "A", "B"]:
            aab += 1
        else:
            abb += 1
    assert aab == b * comb(a, 2)
    assert abb == a * comb(b, 2)


def test_a_single_sided_pool_produces_nothing():
    """No side B means no cross triangle exists, not "all C(n,3) of side A"."""
    assert list(mixed_side(["A"] * 10)) == []
    assert count_from_sides(["A"] * 10) == 0


def test_three_sides_are_refused_rather_than_guessed():
    with pytest.raises(ValueError) as excinfo:
        side_index_lists(["A", "B", "C"])
    assert "exactly two sides" in str(excinfo.value)


# ===========================================================================
# Order and presentation
# ===========================================================================


@pytest.mark.parametrize("a,b", [(2, 3), (3, 4), (5, 2)])
def test_lexicographic_order_is_ascending_and_deterministic(a, b):
    tags = sides(a, b)
    produced = [c.as_tuple() for c in mixed_side(tags)]
    assert produced == sorted(produced)
    assert produced == [c.as_tuple() for c in mixed_side(tags)]


def test_the_lexicographic_generator_leaves_presentation_unset():
    """The standard generator never sets a presentation order."""
    for c in mixed_side(sides(2, 3)):
        assert c.presentation is None
        assert c.vertices() == c.as_tuple()


def test_presentation_must_be_a_permutation_of_the_indices():
    with pytest.raises(ValueError) as excinfo:
        TriangleCandidate(0, 1, 2, (0, 1, 3))
    assert "permutation" in str(excinfo.value)


def test_the_composition_rule_has_one_canonical_name():
    assert COMPOSITION_RULE == "mixed_side_AAB_ABB"


# ===========================================================================
# Interleaved layouts -- the generator must not assume A-block-then-B-block
# ===========================================================================


@pytest.mark.parametrize(
    "tags",
    [
        ["A", "B", "A", "B", "A"],
        ["B", "A", "A", "B"],
        ["A", "B", "B", "A", "B", "A"],
    ],
)
def test_interleaved_side_tags_still_give_the_exact_count_and_set(tags):
    """The generator must not assume an A-block-then-B-block layout: it is
    defined over side tags, whatever their positions."""
    a = tags.count("A")
    b = tags.count("B")
    expected = {
        (i, j, k)
        for i, j, k in combinations(range(len(tags)), 3)
        if len({tags[i], tags[j], tags[k]}) == 2
    }
    produced = [c.as_tuple() for c in mixed_side(tags)]
    assert set(produced) == expected
    assert len(produced) == count(a, b)
    assert produced == sorted(produced)
