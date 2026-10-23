"""The 18 categories: a frozen table, proved closed-form over its own domain.

The table was recovered from the historical ``SSE_TSR.py`` classifier; the tests
here check it independently -- a bijection onto eighteen named categories,
exhaustive over every partition three vertices can realise, order-independent
and determined by the signature alone.
"""

from __future__ import annotations

import itertools

import pytest

from biotsr.sse import (
    CATEGORY_BY_SIGNATURE,
    N_SSE_CATEGORIES,
    SSE_CATEGORIES,
    category_index,
    category_signature,
    classify,
    describe_category,
)
from biotsr.sse.states import (
    HELIX,
    NONE,
    SHEET,
    SSEAssignment,
    SSEElementId,
    UNASSIGNED,
)


def helix(serial):
    """A residue in the helix with this ``HELIX`` serial."""
    return SSEAssignment(HELIX, SSEElementId(serial=serial))


def strand(serial, sheet_id=""):
    """A residue in one individual β-strand, the standard ``(sheetID, serNum)``."""
    return SSEAssignment(SHEET, SSEElementId(serial=serial, sheet_id=sheet_id))


# Every partition three vertices can realise: coil, plus up to three distinct
# helices and three distinct β-strands. Three vertices cannot need a fourth id.
VOCABULARY = (
    [UNASSIGNED]
    + [helix(str(i)) for i in (1, 2, 3)]
    + [strand(str(i), "A") for i in (1, 2, 3)]
)


def _token(assignment):
    return assignment.token


# ---------------------------------------------------------------------------
# The table itself
# ---------------------------------------------------------------------------


def test_there_are_exactly_eighteen_categories_and_they_are_unique():
    assert len(SSE_CATEGORIES) == N_SSE_CATEGORIES == 18
    assert len(set(SSE_CATEGORIES)) == 18


def test_the_category_names_carry_their_own_numbering_in_order():
    """``SSE_CATEGORIES`` is the historical key-file column order, 1..18."""
    for i, name in enumerate(SSE_CATEGORIES):
        assert int(name.split("_")[0]) == i + 1
        assert category_index(name) == i


def test_the_signature_table_is_a_bijection_onto_the_eighteen():
    assert len(CATEGORY_BY_SIGNATURE) == 18
    assert sorted(CATEGORY_BY_SIGNATURE.values()) == sorted(SSE_CATEGORIES)


def test_the_signature_table_is_exhaustive_over_its_own_domain():
    """Every arithmetically possible signature has an entry -- no default branch.

    Enumerated from the definition rather than from the table: split 3 vertices
    across the three states, then count how many distinct helices and distinct
    β-strands they can occupy. If the table had a hole, some triangle would have no category.
    """
    expected = set()
    for n_helix in range(4):
        for n_strand in range(4 - n_helix):
            n_coil = 3 - n_helix - n_strand
            for d_helix in range(1, n_helix + 1) or [0]:
                for d_strand in range(1, n_strand + 1) or [0]:
                    expected.add((n_helix, n_strand, n_coil, d_helix, d_strand))
    assert set(CATEGORY_BY_SIGNATURE) == expected


def test_every_category_has_prose():
    for name in SSE_CATEGORIES:
        assert describe_category(name)


# ---------------------------------------------------------------------------
# Ordering / permutation semantics
# ---------------------------------------------------------------------------


def test_classification_does_not_depend_on_vertex_order():
    """The canonical triplet order cannot leak into the category.

    This is why the engine can hand the annotator its canonically ordered
    vertices without the ordering convention becoming part of the SSE science.
    """
    for combo in itertools.product(VOCABULARY, repeat=3):
        answers = {classify(p) for p in itertools.permutations(combo)}
        assert len(answers) == 1, (combo, answers)


def test_classification_is_deterministic():
    combo = (helix("1"), strand("3", "A"), UNASSIGNED)
    assert {classify(combo) for _ in range(50)} == {"18_1a1b1c"}


# ---------------------------------------------------------------------------
# Spot checks that read as science rather than as table lookups
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "assignment,expected",
    [
        (helix("1"), "1_3a1"),
        (strand("1", "A"), "7_3b1"),
        (strand("1", "BET"), "7_3b1"),
    ],
)
def test_three_vertices_in_one_element(assignment, expected):
    assert classify((assignment,) * 3) == expected


def test_all_coil_is_a_single_category():
    assert classify((UNASSIGNED, UNASSIGNED, UNASSIGNED)) == "17_3c"


def test_distinct_helix_counts_separate_the_first_three_categories():
    h1, h2, h3 = (helix(s) for s in ("1", "2", "3"))
    assert classify((h1, h1, h1)) == "1_3a1"
    assert classify((h1, h1, h2)) == "2_3a2"
    assert classify((h1, h2, h3)) == "3_3a3"


def test_one_of_each_state_is_the_eighteenth():
    assert classify((helix("2"), strand("5", "A"), UNASSIGNED)) == "18_1a1b1c"


def test_signature_is_the_only_thing_a_category_depends_on():
    """Two triangles with different group tokens but the same shape agree."""
    a = (helix("1"), helix("2"), UNASSIGNED)
    b = (helix("7"), helix("9"), UNASSIGNED)
    assert category_signature(a) == category_signature(b)
    assert classify(a) == classify(b) == "5_2a21c"


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


def test_a_coil_residue_cannot_carry_a_group():
    with pytest.raises(ValueError, match="coil residue has no group"):
        SSEAssignment(NONE, SSEElementId(serial="3"))


def test_an_unknown_state_is_refused():
    with pytest.raises(ValueError, match="unknown secondary-structure state"):
        SSEAssignment("TURN", SSEElementId(serial="1"))


def test_a_group_must_be_a_structured_element_id_not_a_string():
    """The guard that keeps ``sheetID + serNum`` from ever being pre-joined."""
    with pytest.raises(TypeError, match="SSEElementId"):
        SSEAssignment(SHEET, "A1")


def test_a_category_needs_exactly_three_vertices():
    with pytest.raises(ValueError, match="exactly 3 vertices"):
        classify((UNASSIGNED, UNASSIGNED))


# ---------------------------------------------------------------------------
# What makes two β-strands the same β-strand (ADR-063)
# ---------------------------------------------------------------------------
#
# The scientific object the categories distinguish is the INDIVIDUAL β-strand.
# A `SHEET` record names one with two fields, because `serNum` is the strand's
# index inside its sheet and restarts at 1 in the next one. These tests are
# the whole truth table of that pair, and they are what separates `7_3b1` from
# `9_3b3` on a real multi-sheet chain.


def test_same_sheet_same_serial_is_one_strand():
    a = strand("1", "A")
    assert category_signature((a, a, a)) == (0, 3, 0, 0, 1)
    assert classify((a, a, a)) == "7_3b1"


def test_same_sheet_different_serial_is_different_strands():
    """Strands 1, 2 and 4 of sheet ``A`` are three strands -- inside one sheet.

    This is the case ``sheetID`` alone would destroy, and it is why sheet id is
    a *disambiguator* here and never the group itself.
    """
    a, b, c = strand("1", "A"), strand("2", "A"), strand("4", "A")
    assert category_signature((a, b, c)) == (0, 3, 0, 0, 3)
    assert classify((a, b, c)) == "9_3b3"
    assert classify((a, a, b)) == "8_3b2"


def test_different_sheet_same_serial_is_different_strands_under_the_standard():
    """Sheet ``A`` strand 1 and sheet ``B`` strand 1 are two physical strands."""
    a, b = strand("1", "A"), strand("1", "B")
    assert a.group != b.group
    assert category_signature((a, a, b)) == (0, 3, 0, 0, 2)
    assert classify((a, a, b)) == "8_3b2"


def test_different_sheet_different_serial_is_different_strands():
    a, b = strand("1", "A"), strand("2", "B")
    assert classify((a, a, b)) == "8_3b2"
    assert classify((a, b, strand("3", "B"))) == "9_3b3"


def test_a_sheet_id_ending_in_a_digit_cannot_merge_with_the_serial():
    """No concatenation anywhere: identity is componentwise.

    ``("A1", "1")`` and ``("A", "11")`` would both spell ``A11`` if the pair were
    ever joined into a string. They are different strands, and they stay
    different.
    """
    a, b = strand("1", "A1"), strand("11", "A")
    assert a.group != b.group
    assert classify((a, a, b)) == "8_3b2"


def test_helix_identity_is_the_serial():
    """A HELIX serNum is unique file-wide, so a helix needs no disambiguator."""
    assert helix("1").group == SSEElementId(serial="1")
    assert helix("1").group.sheet_id == ""
    assert classify((helix("1"), helix("1"), helix("2"))) == "2_3a2"
