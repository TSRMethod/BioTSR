"""Secondary-structure annotation for the protein C-alpha methods.

SSE is an **annotation on the shared result**, never a second engine.  A run with
SSE enabled uses the same selection, the same candidate enumeration, the same
canonical ordering, the same geometry, the same bin schemes, the same key formula
and the same serializers as the run without it -- and produces the *same integer
keys*.  What changes is only which bucket each key's frequency lands in:
``CategoryCounts`` carries one implicit ``"regular"`` bucket for a plain run, two
for mirror, and up to eighteen for SSE, with no schema change (ADR-036, ADR-064).

Three modules, one responsibility each:

* :mod:`biotsr.sse.states` -- the three states, the structured identity of one
  individual helix or β-strand, and one residue's assignment;
* :mod:`biotsr.sse.assignment` -- ``HELIX``/``SHEET`` records to residue states;
* :mod:`biotsr.sse.categories` -- three states to one of the 18 categories.

**What the categories distinguish is the individual β-strand**, not the whole
β-sheet -- a sheet contains one or more strands (ADR-063).  A strand is therefore
named by ``(sheetID, serNum)``, because a ``SHEET`` record's ``serNum`` is its
index *within* its sheet and restarts in the next one.  ``SSE_TSR.py``'s
``serNum``-only identity is kept as the named ``legacy_serial_only`` convention
so historical outputs stay byte-exact; ``sheetID`` alone is refused, because it
would merge the strands of a sheet.

Scope (ADR-036): the 18 categories are defined over
protein residues, so they apply to the C-alpha methods -- ``protein`` and
``aa_selection``, with grouping, mirror and size filtering in any combination --
and to nothing else.  An atom inside a residue has no secondary structure of its
own, and a ligand, a glycan and a nucleotide have none at all; extending the
concept to those families needs its own science, not a reused table.
"""

from __future__ import annotations

from .assignment import (
    HELIX_IDENTITY_BY_SERIAL,
    STRAND_IDENTITIES,
    STRAND_IDENTITY_SHEET_ID_AND_SERIAL,
    SSE_ASSIGNMENT_VERSION,
    SSEAnnotator,
    SSEAssignmentMap,
    build_assignment_map,
)
from .categories import (
    CATEGORY_BY_SIGNATURE,
    N_SSE_CATEGORIES,
    SSE_CATEGORIES,
    category_index,
    category_signature,
    classify,
    describe_category,
)
from .stratification import (
    IncompatibleSSEStratification,
    category_totals_in_order,
    check_conservation,
    require_common_strand_identity,
    sse_annotation_summary,
    strand_identity_of,
    stratified_key_counts,
)
from .states import (
    COIL,
    HELIX,
    NONE,
    SHEET,
    SSE_STATES,
    UNASSIGNED,
    SSEAssignment,
    SSEElementId,
)

__all__ = [
    "HELIX",
    "SHEET",
    "COIL",
    "NONE",
    "SSE_STATES",
    "SSEAssignment",
    "SSEElementId",
    "UNASSIGNED",
    "SSE_CATEGORIES",
    "N_SSE_CATEGORIES",
    "CATEGORY_BY_SIGNATURE",
    "category_signature",
    "category_index",
    "classify",
    "describe_category",
    "STRAND_IDENTITY_SHEET_ID_AND_SERIAL",
    "STRAND_IDENTITIES",
    "HELIX_IDENTITY_BY_SERIAL",
    "SSE_ASSIGNMENT_VERSION",
    "SSEAssignmentMap",
    "SSEAnnotator",
    "build_assignment_map",
    "stratified_key_counts",
    "category_totals_in_order",
    "check_conservation",
    "sse_annotation_summary",
    "strand_identity_of",
    "require_common_strand_identity",
    "IncompatibleSSEStratification",
]
