"""Where a legacy read of a structure stops.

The legacy family implements **two** stopping rules, and the difference is not
cosmetic -- it decides whether a ligand or a glycan that sits after the peptide's
``TER`` record is visible at all.

``LEGACY_SINGLE_PASS``
    Stop at the first of: an ``ENDMDL`` record; a ``TER`` record **whose chain id
    equals the requested chain**; a ``MODEL`` record numbered above 1.  Used by
    ``TSR-Package`` (protein Cα), ``TSR_AMINOACID_PACKAGE``,
    ``TSR_NUCLEOTIDE_PACKAGE`` and ``WebsiteTSR``'s two residue selectors.

``FIRST_MODEL_ONLY``
    Stop at the first ``ENDMDL`` or ``MODEL > 1`` only -- **never** at ``TER``.
    Used by ``WebsiteTSR``'s ligand selector and by ``CarbohydratesTSR``, whose
    own source says why: *"Keep scanning after TER because carbohydrate HETATM
    records can appear later in the same chain."*  Applying the peptide rule to a
    HETATM family would silently select nothing at all.

Because the first rule's second condition depends on which chain was requested,
this is a selection rule, not a parsing rule -- which is why it lives here and not
in the backend.  The parser records the events; these functions decide the cutoff.
"""

from __future__ import annotations

from typing import Optional

from ..domain.identity import ChainId
from ..structure.model import Structure

__all__ = [
    "LEGACY_SINGLE_PASS",
    "FIRST_MODEL_ONLY",
    "TRUNCATION_RULES",
    "legacy_truncation_line",
    "truncation_line",
]

#: Names recorded in provenance.
LEGACY_SINGLE_PASS = "legacy_single_pass_break"
FIRST_MODEL_ONLY = "first_model_only"

TRUNCATION_RULES = (LEGACY_SINGLE_PASS, FIRST_MODEL_ONLY)


def _model_cutoff(structure: Structure):
    for event in structure.truncation_events:
        if event.kind == "ENDMDL":
            yield event.file_order
        elif (
            event.kind == "MODEL"
            and event.model_number is not None
            and event.model_number > 1
        ):
            yield event.file_order


def legacy_truncation_line(structure: Structure, chain: ChainId) -> Optional[int]:
    """The ``file_order`` at which the ``LEGACY_SINGLE_PASS`` reader would stop.

    Records at or after the returned line are not read.  ``None`` means the whole
    file is read.
    """
    candidates = list(_model_cutoff(structure))
    for event in structure.truncation_events:
        if event.kind == "TER" and event.chain == chain:
            candidates.append(event.file_order)
    return min(candidates) if candidates else None


def first_model_truncation_line(structure: Structure) -> Optional[int]:
    """The cutoff under ``FIRST_MODEL_ONLY``: models only, never ``TER``."""
    candidates = list(_model_cutoff(structure))
    return min(candidates) if candidates else None


def truncation_line(
    structure: Structure, chain: ChainId, rule: str = LEGACY_SINGLE_PASS
) -> Optional[int]:
    """Dispatch to the named truncation rule.

    The rule is always named explicitly by the caller's profile; there is no
    default that silently picks the peptide rule for a HETATM family.
    """
    if rule == LEGACY_SINGLE_PASS:
        return legacy_truncation_line(structure, chain)
    if rule == FIRST_MODEL_ONLY:
        return first_model_truncation_line(structure)
    raise ValueError(
        "unsupported truncation rule %r; this package implements %s"
        % (rule, ", ".join(repr(r) for r in TRUNCATION_RULES))
    )
