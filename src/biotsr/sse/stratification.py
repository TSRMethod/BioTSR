"""Reading the SSE stratification off a result.  Reshaping, never analysis.

The analysis layer consumes ``key -> total frequency`` and is defined over
exactly that: every evidenced operation in the family -- feature matrices, the
four common-key operations, the five similarity measures, the clustering -- sums a
key's categories before it does anything, and not one of them distinguishes a
mirror bucket from a regular one or one SSE bucket from another.  That is a fact
about the sources, not an oversight, and ``AnalysisSample.from_result`` is
written to honour it.

So SSE support does **not** define an SSE feature representation.  No authoritative
source establishes one, and inventing a 18x-wider feature vector would be
inventing science, not migrating it (ADR-067).  What it does instead is
keep the stratification intact on the result and expose it in the shape a future,
source-backed analysis would want -- clearly labelled as data, in a module that
imports nothing from ``biotsr.analysis`` and offers no measure, no distance and no
clustering.

The two views never get confused because they are different shapes:

* ``result.frequencies``          -- ``{key: n}``, what analysis consumes;
* ``stratified_key_counts(...)``  -- ``{category: {key: n}}``, what SSE adds.

Summing the second over categories gives the first, and
:func:`check_conservation` asserts exactly that.

**The two views also have different combination rules, and that asymmetry is the
point.**  A category is a statement about β-strand identity, so two results
annotated under different strand-identity conventions do *not* agree about what
``7_3b1`` means: on 1TUP chain A the same triangle is ``7_3b1`` under
``legacy_serial_only`` and ``8_3b2`` under ``sheet_id_and_serial`` (ADR-063).
An integer key means the same thing either way, because the convention is in
neither identity -- so:

* **ordinary key-frequency analysis may mix them freely.**
  ``AnalysisSample.from_result`` sums the categories away and keeps only
  ``KeySpace`` and ``SamplingProfile``, so the convention never reaches a
  measure. That is structural, not a promise;
* **anything stratified by category may not.** There is no such analysis today
  (ADR-067 defines no SSE feature representation), which is exactly why the
  guard exists before the analysis does:
  :func:`require_common_strand_identity` is what one must call, and
  :class:`IncompatibleSSEStratification` is what it raises.
"""

from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence, Tuple

from ..domain.errors import IncompatibleRepresentation
from ..domain.result import TSRResult
from .categories import SSE_CATEGORIES

__all__ = [
    "IncompatibleSSEStratification",
    "stratified_key_counts",
    "category_totals_in_order",
    "check_conservation",
    "sse_annotation_summary",
    "strand_identity_of",
    "require_common_strand_identity",
]


class IncompatibleSSEStratification(IncompatibleRepresentation):
    """Raised when category breakdowns that do not mean the same thing would be
    combined.

    The sibling of ``IncompatibleKeySpace`` and ``IncompatibleSampling``, for the
    axis those two deliberately do not cover: SSE annotation provenance is in
    neither identity, so neither of those exceptions can ever fire on it.
    """


def stratified_key_counts(result: TSRResult) -> Mapping[str, Mapping[int, int]]:
    """``{category: {key: frequency}}``, one inner mapping per observed category.

    Sparse in both directions: a category no triangle reached is absent, and a
    key absent from a category's mapping has frequency zero there.  The union of
    the inner mappings' keys is ``result.keys``, and the sums reconcile exactly.
    """
    out: Dict[str, Dict[int, int]] = {}
    for key, counts in result.keys.items():
        for category, n in counts.by_category.items():
            if n:
                out.setdefault(category, {})[key] = n
    return out


def category_totals_in_order(result: TSRResult) -> Tuple[Tuple[str, int], ...]:
    """The 18 categories and their totals, in the canonical 1..18 order.

    Every category is present, including the ones with zero, because a dense,
    fixed-width vector is what makes two structures' distributions comparable
    without first agreeing on which categories each of them happened to reach.
    """
    totals = result.category_totals()
    return tuple((category, totals.get(category, 0)) for category in SSE_CATEGORIES)


def check_conservation(result: TSRResult) -> None:
    """Raise unless the categories partition the retained triangles exactly.

    Two things are asserted, and together they are the frequency-conservation
    property SSE support claims:

    1. every key's category vector sums to that key's own frequency;
    2. the category totals sum to ``n_triangles``.

    Both hold unconditionally, with no missing-SSE escape clause, because coil is
    a *state* rather than an absence: a residue in no ``HELIX``/``SHEET`` record
    is classified, not skipped (ADR-065). There is nothing for a triangle to fall
    through -- which is precisely what makes this checkable at all, and what the
    legacy key file fails.
    """
    for key, counts in result.keys.items():
        total = sum(counts.by_category.values())
        if total != counts.total:
            raise AssertionError(
                "key %d: category vector sums to %d but its frequency is %d"
                % (key, total, counts.total)
            )
    grand = sum(result.category_totals().values())
    if grand != result.n_triangles:
        raise AssertionError(
            "category totals sum to %d but %d triangles were retained"
            % (grand, result.n_triangles)
        )


def sse_annotation_summary(result: TSRResult) -> Mapping[str, object]:
    """Everything a reader needs to judge an SSE result, in one mapping.

    Includes ``has_records``: without it, "every triangle is ``17_3c``" reads as
    a finding when it may only mean the structure carried no annotation.
    """
    annotation = result.annotations.get("sse")
    return {
        "annotated": annotation is not None,
        "assignment": annotation,
        "coverage": result.annotations.get("sse_coverage"),
        "category_totals": dict(category_totals_in_order(result)),
        "n_triangles": result.n_triangles,
        "n_categories_observed": sum(
            1 for _, n in category_totals_in_order(result) if n
        ),
    }


def strand_identity_of(result: TSRResult) -> str:
    """Which β-strand identity convention produced this result's categories.

    Raises if the result carries no SSE annotation at all: an unannotated run has
    no categories, so asking which convention named them is a category error
    rather than a question with a default answer.
    """
    annotation = result.annotations.get("sse")
    if not annotation:
        raise IncompatibleSSEStratification(
            "this result carries no SSE annotation, so it has no strand-identity "
            "convention and nothing to stratify. Re-run with sse=True."
        )
    return annotation["strand_identity"]


def require_common_strand_identity(
    results: Sequence[TSRResult],
    *,
    context: Optional[str] = None,
) -> str:
    """Refuse to combine category breakdowns that do not mean the same thing.

    **Call this before any operation that reads categories across more than one
    result** -- comparing category vectors, stacking stratified key counts into a
    matrix, correlating category distributions.  ``KeySpace`` and
    ``SamplingProfile`` cannot catch this: the strand-identity convention is
    deliberately in neither identity, because it moves no key and draws no
    triangle (ADR-063), so two results that differ *only* here have identical
    fingerprints and identical integer keys -- and genuinely incomparable
    categories.

    Ordinary key-frequency analysis does not need this and must not use it:
    mixing the conventions there is correct, and ``AnalysisSample.from_result``
    already makes it unavoidable by summing the categories away.

    Returns the one convention they share, so a caller can record it.
    """
    if not results:
        raise ValueError(
            "no results: a strand-identity convention cannot be established"
        )
    identities = [(r, strand_identity_of(r)) for r in results]
    first_result, first = identities[0]
    for other_result, other in identities[1:]:
        if other != first:
            raise IncompatibleSSEStratification(
                "%sSSE categories from these results are not comparable: %s was "
                "annotated under strand identity %r and %s under %r. A category "
                "counts DISTINCT INDIVIDUAL beta-strands, and the two conventions "
                "disagree about when two sheet residues are in the same strand -- "
                "%r pairs the SHEET record's sheetID with its serNum, %r uses the "
                "serNum alone and so merges strands that share a serial across "
                "different sheets (ADR-063). Their integer keys ARE comparable "
                "and their KeySpace and SamplingProfile fingerprints match, "
                "because the convention is in neither identity; it is the "
                "per-category breakdown that is not. BioTSR annotates under "
                "%r only, so re-run the other result with this version, or "
                "compare result.frequencies instead of the categories."
                % (
                    "%s: " % context if context else "",
                    _name(first_result),
                    first,
                    _name(other_result),
                    other,
                    "sheet_id_and_serial",
                    "legacy_serial_only",
                    "sheet_id_and_serial",
                )
            )
    return first


def _name(result: TSRResult) -> str:
    """A readable identifier for one result, for use in a refusal message."""
    chain = getattr(result, "chain", None)
    return "%s/%s" % (result.structure_id.value, chain) if chain else (
        result.structure_id.value
    )
