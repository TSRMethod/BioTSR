"""Selection: which points participate, and how they group into one calculation."""

from .atomic import (
    STANDARD_AMINOACID_RESIDUE_ADMISSION,
    SCOPE_PER_ENTITY,
    SCOPE_POOLED,
    SELECTOR_ATOM_SELECTION,
    SELECTOR_RESIDUE_ATOMS,
    WEBSITE_TSR_ADMISSION,
    AtomAdmission,
    AtomSelectionRequest,
    parse_atom_selection,
    residue_targets,
    select_atom_selection,
    select_residue_atoms,
)
from .base import AltLocPolicy, SelectionResult
from .completeness import (
    EXACT_LEGACY_ATOM_COUNT,
    LEGACY_AMINOACID_ATOM_COUNTS,
    NO_COMPLETENESS_POLICY,
    Admission,
    CompletenessPolicy,
    ExactAtomCountPolicy,
)
from .protein import SCOPE, SELECTOR_NAME, select_calpha
from .truncation import LEGACY_SINGLE_PASS, legacy_truncation_line

__all__ = [
    "Admission",
    "AltLocPolicy",
    "AtomAdmission",
    "AtomSelectionRequest",
    "CompletenessPolicy",
    "EXACT_LEGACY_ATOM_COUNT",
    "ExactAtomCountPolicy",
    "STANDARD_AMINOACID_RESIDUE_ADMISSION",
    "LEGACY_AMINOACID_ATOM_COUNTS",
    "LEGACY_SINGLE_PASS",
    "NO_COMPLETENESS_POLICY",
    "SCOPE",
    "SCOPE_PER_ENTITY",
    "SCOPE_POOLED",
    "SELECTOR_ATOM_SELECTION",
    "SELECTOR_NAME",
    "SELECTOR_RESIDUE_ATOMS",
    "SelectionResult",
    "WEBSITE_TSR_ADMISSION",
    "legacy_truncation_line",
    "parse_atom_selection",
    "residue_targets",
    "select_atom_selection",
    "select_calpha",
    "select_residue_atoms",
]
