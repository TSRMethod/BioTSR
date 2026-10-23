"""The 18 SSE triplet categories.

**Recovered from source, not from memory.**  ``determine_type`` in
``TSR-Package/tsr_package/SSE_TSR.py`` is a chain of overlapping ``if`` blocks
with chained ``!=`` comparisons, so what it computes is not readable
off the page.  It was instead *executed* over its entire reachable input space --
all 343 triples drawn from ``{NONE_0} u {HELIX_1..3} u {SHEET_1..3}``, which is
every partition three vertices can realise -- and the result is:

* a **total** function: no input returns the unclassified sentinel ``0``;
* **permutation invariant**: all 6 orderings of every triple give one answer, so
  the canonical vertex order the engine produces does not enter the category;
* **exactly 18** reachable codes.

Executing it also collapses the chained-comparison mess into a closed form.  The
category depends on one signature and nothing else::

    (n_helix, n_sheet, n_coil, n_distinct_helices, n_distinct_strands)

and that signature space has exactly 18 inhabitants, each mapping to one code --
a bijection, with no default branch and nothing left over.  That is what
:data:`CATEGORY_BY_SIGNATURE` is, and
``tests/unit/test_sse_categories.py`` re-derives the whole table from the legacy
function itself, so the two can never drift apart silently.

The names are the legacy ones because they are the published vocabulary -- they
are the column headers of every historical ``.3Dkeys_*_SSE`` file.  Reading them:
``a`` = helix, ``b`` = **β-strand**, ``c`` = coil; ``3a1`` is "3 vertices, helix,
1 distinct helix"; ``2a21b`` is "2 vertices in 2 distinct helices, 1 in a
β-strand".

``b`` is a **strand**, not a sheet.  The object the categories distinguish is the
individual β-strand -- one β-sheet contains one or more of them -- so ``9_3b3``
says three *different strands*, which is satisfiable inside a single sheet and is
meant to be (ADR-063).  The legacy comments beside these codes say "sheet" where
they mean "strand"; that wording is historical and is not repeated here.  Which
PDB fields name an individual strand is
``biotsr.sse.assignment``'s convention, not this module's business: this module
counts how many *distinct* groups three vertices fall in and never inspects one.

**A category is not a key.**  It is the bucket a key's frequency is counted in
(``CategoryCounts``), which is why annotating a run with SSE cannot move a single
integer key -- see ADR-064 and ``tests/integration/test_sse_composition.py``.
"""

from __future__ import annotations

from typing import Mapping, Sequence, Tuple

from .states import HELIX, NONE, SHEET, SSEAssignment

__all__ = [
    "SSE_CATEGORIES",
    "CATEGORY_BY_SIGNATURE",
    "N_SSE_CATEGORIES",
    "category_signature",
    "classify",
    "category_index",
    "describe_category",
]

#: The 18 categories, in the legacy numeric order -- which is the column order of
#: the historical key file, so the tuple doubles as that file's schema.
SSE_CATEGORIES = (
    "1_3a1",
    "2_3a2",
    "3_3a3",
    "4_2a11c",
    "5_2a21c",
    "6_1a2c",
    "7_3b1",
    "8_3b2",
    "9_3b3",
    "10_2b11c",
    "11_2b21c",
    "12_1b2c",
    "13_2a21b",
    "14_2a11b",
    "15_2b21a",
    "16_2b11a",
    "17_3c",
    "18_1a1b1c",
)

N_SSE_CATEGORIES = 18

#: ``(n_helix, n_sheet, n_coil, n_distinct_helices, n_distinct_strands) -> code``.
#:
#: Exhaustive by construction: the ten ways to split 3 vertices across the three
#: states, refined by how many distinct helices and distinct β-strands they fall
#: in, give 3+3+1+2+2+2+2+1+1+1 = 18 signatures. There is no default and no
#: fallback,
#: because there is nothing for one to catch.
CATEGORY_BY_SIGNATURE: Mapping[Tuple[int, int, int, int, int], str] = {
    # three helix vertices
    (3, 0, 0, 1, 0): "1_3a1",
    (3, 0, 0, 2, 0): "2_3a2",
    (3, 0, 0, 3, 0): "3_3a3",
    # two helix + one coil
    (2, 0, 1, 1, 0): "4_2a11c",
    (2, 0, 1, 2, 0): "5_2a21c",
    # one helix + two coil
    (1, 0, 2, 1, 0): "6_1a2c",
    # three beta-strand vertices
    (0, 3, 0, 0, 1): "7_3b1",
    (0, 3, 0, 0, 2): "8_3b2",
    (0, 3, 0, 0, 3): "9_3b3",
    # two beta-strand + one coil
    (0, 2, 1, 0, 1): "10_2b11c",
    (0, 2, 1, 0, 2): "11_2b21c",
    # one beta-strand + two coil
    (0, 1, 2, 0, 1): "12_1b2c",
    # two helix + one beta-strand
    (2, 1, 0, 2, 1): "13_2a21b",
    (2, 1, 0, 1, 1): "14_2a11b",
    # one helix + two beta-strand
    (1, 2, 0, 1, 2): "15_2b21a",
    (1, 2, 0, 1, 1): "16_2b11a",
    # three coil
    (0, 0, 3, 0, 0): "17_3c",
    # one of each
    (1, 1, 1, 1, 1): "18_1a1b1c",
}

#: What each category says, in the corrected scientific vocabulary: ``b`` counts
#: **individual β-strands**, not β-sheets (ADR-063).  The legacy source's own
#: comments say "sheet"; that historical wording is preserved where it is quoted
#: as evidence and is not restated here as though it were the science.
_DESCRIPTIONS = {
    "1_3a1": "3 vertices from the same helix",
    "2_3a2": "3 vertices from 2 distinct helices",
    "3_3a3": "3 vertices from 3 distinct helices",
    "4_2a11c": "2 vertices from the same helix, 1 from coil",
    "5_2a21c": "2 vertices from 2 distinct helices, 1 from coil",
    "6_1a2c": "1 vertex from a helix, 2 from coil",
    "7_3b1": "3 vertices from the same β-strand",
    "8_3b2": "3 vertices from 2 distinct β-strands",
    "9_3b3": "3 vertices from 3 distinct β-strands",
    "10_2b11c": "2 vertices from the same β-strand, 1 from coil",
    "11_2b21c": "2 vertices from 2 distinct β-strands, 1 from coil",
    "12_1b2c": "1 vertex from a β-strand, 2 from coil",
    "13_2a21b": "2 vertices from 2 distinct helices, 1 from a β-strand",
    "14_2a11b": "2 vertices from the same helix, 1 from a β-strand",
    "15_2b21a": "2 vertices from 2 distinct β-strands, 1 from a helix",
    "16_2b11a": "2 vertices from the same β-strand, 1 from a helix",
    "17_3c": "all 3 vertices from coil",
    "18_1a1b1c": "1 vertex from a helix, 1 from a β-strand, 1 from coil",
}


def category_signature(
    assignments: Sequence[SSEAssignment],
) -> Tuple[int, int, int, int, int]:
    """The five-number signature of one triangle's three vertex assignments.

    The last two numbers are how many *distinct* helices and how many *distinct*
    individual β-strands the three vertices fall in.  Distinctness is decided by
    :class:`~biotsr.sse.states.SSEElementId` equality -- componentwise, so two
    strands that share a serial in different sheets are two strands under the
    standard convention and one under the historical one.
    """
    if len(assignments) != 3:
        raise ValueError(
            "an SSE category is defined for exactly 3 vertices; got %d"
            % len(assignments)
        )
    helices = set()
    strands = set()
    n_helix = n_sheet = n_coil = 0
    for a in assignments:
        if a.state == HELIX:
            n_helix += 1
            helices.add(a.group)
        elif a.state == SHEET:
            n_sheet += 1
            strands.add(a.group)
        elif a.state == NONE:
            n_coil += 1
        else:  # pragma: no cover - SSEAssignment validates this at construction
            raise ValueError("unknown state %r" % a.state)
    return (n_helix, n_sheet, n_coil, len(helices), len(strands))


def classify(assignments: Sequence[SSEAssignment]) -> str:
    """The SSE category of one triangle.  Total: every input has an answer."""
    signature = category_signature(assignments)
    try:
        return CATEGORY_BY_SIGNATURE[signature]
    except KeyError:  # pragma: no cover - the table is exhaustive over its domain
        raise AssertionError(
            "signature %r has no category, which cannot happen: the table covers "
            "every way 3 vertices can be distributed. This means the signature "
            "function and the table have diverged." % (signature,)
        )


def category_index(category: str) -> int:
    """0-based column index of a category in the legacy 18-column layout."""
    return SSE_CATEGORIES.index(category)


def describe_category(category: str) -> str:
    return _DESCRIPTIONS[category]
