"""Cross side definitions -- which points each side of a cross contributes.

A cross side is a *selector plus its arguments*, nothing more.  Every selector
named here already exists and is used unchanged by the non-cross methods; this
module adds no selection logic of its own, it only says which existing selector a
given side of a given cross family uses and with which admission rule.  The drug
side of drug-protein cross is the same selection as regular drug TSR.

Every side is chain-scoped: both chains of a cross calculation are named by the
caller, and no other chain is scanned (ADR-049).

The residue vocabularies are literal transcriptions of the sets in the
historical contact programs (``Nucleotide-Protein/drug_protein_3A.py`` and
``CarbohydratesTSR``'s ``drug_protein_cutArmst.py``), kept here rather than in
a config file so the scope is visible at the point it is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from ..labels.schemes import LabelScheme
from ..selection import atomic as atomic_selection
from ..selection import protein as protein_selection
from ..selection.atomic import AtomAdmission
from ..selection.base import AltLocPolicy, SelectionResult
from ..selection.completeness import (
    MODIFIED_AMINO_ACIDS,
    NO_COMPLETENESS_POLICY,
    CompletenessPolicy,
)
from ..selection.truncation import FIRST_MODEL_ONLY, LEGACY_SINGLE_PASS
from ..structure.model import Structure

__all__ = [
    "PROTEIN_RESIDUE_NAMES",
    "NUCLEOTIDE_RESIDUE_NAMES",
    "STANDARD_CROSS_PROTEIN_ADMISSION",
    "STANDARD_CROSS_ENTITY_ADMISSION",
    "SIDE_SELECTOR_PROTEIN_CHAIN",
    "SIDE_SELECTOR_ENTITY",
    "SIDE_SELECTOR_ENTITY_CHAIN",
    "SIDE_SELECTOR_RESIDUE_TYPES",
    "SIDE_SELECTOR_SPECIFIC_RESIDUES",
    "SIDE_SELECTOR_ATOM_SELECTION",
    "SIDE_SELECTOR_CALPHA_CHAIN",
    "SIDE_SELECTOR_CALPHA_SELECTION",
    "CROSS_SIDE_SELECTORS",
    "CrossSideSpec",
    "select_side",
]

#: ``Nucleotide-Protein/drug_protein_3A.py`` and
#: ``CarbohydratesTSR/.../drug_protein_cutArmst.py``, verbatim: the twenty
#: standard residues plus the three phosphorylated ones.
PROTEIN_RESIDUE_NAMES = (
    "GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "SER",
    "THR", "MET", "CYS", "HIS", "LYS", "ARG", "ASP", "GLU", "ASN", "GLN",
    "TPO", "SEP", "PTR",
)

#: ``Nucleotide-Protein/drug_protein_3A.py``'s ``drug_residues``: RNA then DNA.
NUCLEOTIDE_RESIDUE_NAMES = ("G", "A", "C", "U", "DG", "DA", "DC", "DT")


#: The protein side of the carbohydrate- and nucleotide-protein cross:
#: ``ATOM`` records plus the recognised modified amino acids (TPO, SEP, PTR;
#: ADR-101), **non-hydrogen only** by the element field (columns 77-78), the
#: ordinary altLoc allow-list, read up to the protein chain's own ``TER`` like
#: every other protein side.  Until 1.0.0rc5 this side read to the end of model
#: 1, as the two historical contact generators did, so a residue written after
#: the chain's ``TER`` (a free phosphotyrosine, say) joined the protein; since
#: 1.0.0rc6 it does not (ADR-109).  The molecule side keeps reading the whole of
#: model 1, so a glycan or nucleotide written after the protein's ``TER`` is
#: still selectable.  Hydrogen exclusion is a settled standard rule for all
#: cross calculations -- contact detection, interface identification, pools and
#: vertices alike (ADR-050).
STANDARD_CROSS_PROTEIN_ADMISSION = AtomAdmission(
    name="standard_cross_protein",
    excluded_elements=("H", "D"),
    accepted_alt_locs=("", "A"),
    truncation=LEGACY_SINGLE_PASS,
    hetero_residue_names=MODIFIED_AMINO_ACIDS,
)

#: The standard entity side of carbohydrate and nucleotide cross.  ``ATOM`` **and**
#: ``HETATM``, because a nucleotide is an ``ATOM`` record and a glycan a
#: ``HETATM`` one; heavy atoms only by the element field; the whole of model 1,
#: because a glycan or a nucleic-acid strand may follow the protein's ``TER``.
STANDARD_CROSS_ENTITY_ADMISSION = AtomAdmission(
    name="standard_cross_entity",
    excluded_elements=("H", "D"),
    accepted_alt_locs=("", "A"),
    hetero_allowed=True,
    truncation=FIRST_MODEL_ONLY,
)


SIDE_SELECTOR_PROTEIN_CHAIN = "cross_protein_chain_atoms"
SIDE_SELECTOR_ENTITY = "cross_entity_atoms"
SIDE_SELECTOR_ENTITY_CHAIN = "cross_entity_chain_atoms"
SIDE_SELECTOR_RESIDUE_TYPES = "cross_residue_type_atoms"
SIDE_SELECTOR_SPECIFIC_RESIDUES = "cross_specific_residue_atoms"
SIDE_SELECTOR_ATOM_SELECTION = "cross_atom_selection"
SIDE_SELECTOR_CALPHA_CHAIN = "cross_calpha_chain"
SIDE_SELECTOR_CALPHA_SELECTION = "cross_calpha_selection"

#: Every selector a cross side may use, for validation and provenance.
CROSS_SIDE_SELECTORS = (
    SIDE_SELECTOR_PROTEIN_CHAIN,
    SIDE_SELECTOR_ENTITY,
    SIDE_SELECTOR_ENTITY_CHAIN,
    SIDE_SELECTOR_RESIDUE_TYPES,
    SIDE_SELECTOR_SPECIFIC_RESIDUES,
    SIDE_SELECTOR_ATOM_SELECTION,
    SIDE_SELECTOR_CALPHA_CHAIN,
    SIDE_SELECTOR_CALPHA_SELECTION,
)

#: Selectors whose points are **residues** (one Ca per residue), not generic
#: atoms.  They label in the residue label space at radix 30; every other
#: selector labels atoms by element (radix 112).  Mixing the two inside one
#: cross is refused rather than resolved -- see ``api.analyze_cross``.
CALPHA_SIDE_SELECTORS = (SIDE_SELECTOR_CALPHA_CHAIN, SIDE_SELECTOR_CALPHA_SELECTION)


@dataclass(frozen=True)
class CrossSideSpec:
    """One side of a cross: a selector, its chain, and its arguments.

    ``role`` is a human name for provenance ("protein", "drug", "carbohydrate",
    "nucleotide") and never affects selection.
    """

    selector: str
    chain: str
    role: str = ""
    entity_name: Optional[str] = None
    entity_id: Optional[str] = None
    entity_names: Tuple[str, ...] = ()
    residue_types: Tuple[str, ...] = ()
    residues: Tuple[str, ...] = ()
    atoms: Tuple[str, ...] = ()
    atom_selection: str = ""

    def __post_init__(self) -> None:
        if self.selector not in CROSS_SIDE_SELECTORS:
            raise ValueError(
                "unknown cross side selector %r; known selectors: %s"
                % (self.selector, ", ".join(CROSS_SIDE_SELECTORS))
            )

    @property
    def chain_label(self) -> str:
        """The chain as provenance records it: the one the caller named."""
        return self.chain

    @property
    def is_calpha(self) -> bool:
        return self.selector in CALPHA_SIDE_SELECTORS

    def describe(self) -> str:
        bits = [
            self.selector,
            "chain=%s" % self.chain,
        ]
        if self.role:
            bits.insert(0, self.role)
        if self.entity_name:
            # No identifier means every copy on the chain; said so, so the
            # side never reads as if one unnamed copy had been chosen.
            bits.append("entity=%s%s" % (
                self.entity_name,
                "_%s" % self.entity_id if self.entity_id else "(all copies)"))
        if self.entity_names:
            bits.append("entity_names=%s" % ";".join(self.entity_names))
        if self.residue_types:
            bits.append("residue_types=%s" % ";".join(self.residue_types))
        if self.residues:
            bits.append("residues=%s" % ";".join(self.residues))
        if self.atoms:
            bits.append("atoms=%s" % ";".join(self.atoms))
        if self.atom_selection:
            bits.append("atom_selection=%s" % self.atom_selection)
        return " ".join(bits)


def select_side(
    structure: Structure,
    spec: CrossSideSpec,
    label_scheme: LabelScheme,
    *,
    admission: AtomAdmission,
    completeness: CompletenessPolicy = NO_COMPLETENESS_POLICY,
) -> SelectionResult:
    """Run one cross side's selector.  Pure dispatch -- no selection logic here."""
    sel = spec.selector

    if sel == SIDE_SELECTOR_PROTEIN_CHAIN:
        # `entity_names=()` means "no residue-name filter", which is the
        # record-type-only rule `Drug_TSR.py`/`Cross_TSR.py` use; the two contact
        # generators pass the explicit 23-name vocabulary instead.
        return atomic_selection.select_entity_atoms(
            structure, spec.chain, label_scheme,
            admission=admission,
            entity_names=spec.entity_names or None,
            completeness=completeness,
            selector=SIDE_SELECTOR_PROTEIN_CHAIN,
        )

    if sel == SIDE_SELECTOR_ENTITY:
        if not spec.entity_name:
            raise ValueError(
                "cross side %r needs entity_name (and usually entity_id)" % sel)
        return atomic_selection.select_entity_atoms(
            structure, spec.chain, label_scheme,
            admission=admission,
            entity_names=(spec.entity_name,),
            entity_seqs=(str(spec.entity_id),) if spec.entity_id not in (None, "") else None,
            atom_names=tuple(spec.atoms) or None,
            completeness=completeness,
            selector=SIDE_SELECTOR_ENTITY,
        )

    if sel == SIDE_SELECTOR_ENTITY_CHAIN:
        return atomic_selection.select_entity_atoms(
            structure, spec.chain, label_scheme,
            admission=admission,
            entity_names=tuple(spec.entity_names) or None,
            completeness=completeness,
            selector=SIDE_SELECTOR_ENTITY_CHAIN,
        )

    if sel == SIDE_SELECTOR_RESIDUE_TYPES:
        if not spec.residue_types:
            raise ValueError("cross side %r needs residue_types" % sel)
        return atomic_selection.select_entity_atoms(
            structure, spec.chain, label_scheme,
            admission=admission,
            entity_names=tuple(spec.residue_types),
            completeness=completeness,
            selector=SIDE_SELECTOR_RESIDUE_TYPES,
        )

    if sel == SIDE_SELECTOR_SPECIFIC_RESIDUES:
        refs, _ = atomic_selection.parse_residue_refs(";".join(spec.residues))
        if not refs:
            raise ValueError("cross side %r needs residues, e.g. 'LYS_6;GLU18'" % sel)
        return atomic_selection.select_entity_atoms(
            structure, spec.chain, label_scheme,
            admission=admission,
            entity_refs=refs,
            completeness=completeness,
            selector=SIDE_SELECTOR_SPECIFIC_RESIDUES,
        )

    if sel == SIDE_SELECTOR_ATOM_SELECTION:
        requests, _ = atomic_selection.parse_atom_selection(spec.atom_selection)
        if not requests:
            raise ValueError(
                "cross side %r needs atom_selection, e.g. 'LYS_48_NZ;HIS_68_NE2'"
                % sel
            )
        return atomic_selection.select_atom_selection(
            structure, spec.chain, requests, label_scheme,
            admission=admission,
        )

    if sel == SIDE_SELECTOR_CALPHA_CHAIN:
        return protein_selection.select_calpha(
            structure, spec.chain, label_scheme,
            altloc_policy=AltLocPolicy.BLANK_OR_A,
            truncation=LEGACY_SINGLE_PASS,
        )

    if sel == SIDE_SELECTOR_CALPHA_SELECTION:
        labels, _ = protein_selection.parse_residue_selection(";".join(spec.residues))
        if not labels:
            raise ValueError(
                "cross side %r needs residues, e.g. '10-12;48;70-71'" % sel)
        return protein_selection.select_calpha_residues(
            structure, spec.chain, labels, label_scheme,
            altloc_policy=AltLocPolicy.BLANK_OR_A,
            truncation=LEGACY_SINGLE_PASS,
        )

    raise ValueError("unhandled cross side selector %r" % sel)
