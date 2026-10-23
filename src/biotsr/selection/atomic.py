"""Atomic selection -- the shared foundation for every non-Cα family.

    Structure -> chain/residue target -> atom admission -> completeness policy
              -> eligible atoms -> labelled atomic points

Nothing here is amino-acid-specific except the two recipes' *configuration*.  The
selector reads an admission policy, a label scheme and (optionally) a
completeness policy; the drug, carbohydrate and nucleotide families reuse this
module by supplying different ones, not by copying it.

Two scopes, because the legacy family really does implement two:

* :func:`select_residue_atoms` -- **per-entity**: one calculation per residue, all
  of that residue's admitted heavy atoms.  This is
  ``TSR_AMINOACID_PACKAGE``'s scope, and the only one with a completeness rule.
* :func:`select_atom_selection` -- **pooled**: the explicitly named atoms, however
  many residues they span, triangulated together as one set.  This is
  ``WebsiteTSR``'s ``residue_atom`` scope, and the one the historical dataset sheet's
  ``aa_atom_selection`` column drives.  It has no completeness rule, because
  asking for four atoms of a residue and then rejecting the residue for having
  only four would be incoherent.

Every point carries full physical identity (``AtomIdentity``), never a bare atom
name.  Within a single residue atom names happen to be unique, which is why the
legacy per-residue code gets away with a name-keyed dict; pooling across residues
is exactly where that assumption breaks, and it is the demonstrated cross-family
defect.  Here it cannot arise: identity is structural, not nominal.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import List, Optional, Sequence, Tuple

from ..domain.errors import InvalidRequest

from ..domain.identity import ChainId
from ..domain.point import LabeledPoint, PointSet
from ..domain.result import SkipRecord
from ..domain.selection_record import STRUCTURE_NOTES
from ..labels.schemes import LabelScheme, UnlabelledPointPolicy
from ..structure.model import Structure
from . import integrity
from .base import SelectionResult
from .completeness import (
    MODIFIED_AMINO_ACIDS,
    NO_COMPLETENESS_POLICY,
    Admission,
    CompletenessPolicy,
)
from .truncation import FIRST_MODEL_ONLY, LEGACY_SINGLE_PASS, truncation_line

__all__ = [
    "SELECTOR_RESIDUE_ATOMS",
    "SELECTOR_ATOM_SELECTION",
    "SELECTOR_RESIDUE_TYPE_ATOMS",
    "SELECTOR_SPECIFIC_RESIDUE_ATOMS",
    "SELECTOR_ENTITY_ATOMS",
    "SELECTOR_ENTITY_CHAIN_ATOMS",
    "SELECTOR_ENTITY_ATOM_SELECTION",
    "SCOPE_PER_ENTITY",
    "SCOPE_POOLED",
    "ATOM_SET_NOTE",
    "AtomAdmission",
    "STANDARD_AMINOACID_RESIDUE_ADMISSION",
    "WEBSITE_TSR_ADMISSION",
    "STANDARD_HETERO_ADMISSION",
    "AtomSelectionRequest",
    "ResidueRef",
    "parse_atom_selection",
    "parse_residue_types",
    "parse_residue_refs",
    "parse_residue_id",
    "parse_atom_names",
    "residue_targets",
    "hetero_amino_acid_residues",
    "entity_targets",
    "select_residue_atoms",
    "select_atom_selection",
    "select_entity_atoms",
]

SELECTOR_RESIDUE_ATOMS = "amino_acid_residue_atoms"
SELECTOR_ATOM_SELECTION = "amino_acid_atom_selection"
#: ``aa_ind`` / ``aa_ind_mix`` -- residues chosen by amino-acid *type*.
SELECTOR_RESIDUE_TYPE_ATOMS = "residue_type_atoms"
#: ``aa_ind_specific`` -- residues chosen by explicit ``NAME_SEQ`` identity.
SELECTOR_SPECIFIC_RESIDUE_ATOMS = "specific_residue_atoms"
#: ``drug_name``/``drug_id``, ``carb_name``/``carb_id``, ``nt_name``/``nt_id`` --
#: one named entity instance.
SELECTOR_ENTITY_ATOMS = "entity_atoms"
#: Polymer scope: every admitted entity of a listed set, pooled across the chain.
SELECTOR_ENTITY_CHAIN_ATOMS = "entity_chain_atoms"
#: ``drug_atom_selected`` / ``carb_atom_selected`` / ``nt_atom_selected``.
SELECTOR_ENTITY_ATOM_SELECTION = "entity_atom_selection"
SCOPE_PER_ENTITY = "per_entity"
SCOPE_POOLED = "pooled"

#: The warning that names the selected instances whose atoms differ from their
#: type's standard set (ADR-114).  It describes the structure file -- every
#: atom present was used -- so it never makes a result partial.
ATOM_SET_NOTE = STRUCTURE_NOTES[0] + (" (as the structure file has them; every "
                                      "atom present is used): %s")


# ---------------------------------------------------------------------------
# Atom admission -- which records are eligible at all
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtomAdmission:
    """Which ``ATOM`` records an implementation is willing to look at.

    The alternate-location rule is the same for every method: an atom whose
    altLoc indicator is blank or ``A`` is admitted, and ``B``, ``C`` and every
    later indicator is excluded, with a skip record naming the atom.  There is
    no fallback: a residue whose atoms exist only as ``B`` (or later)
    conformers contributes none of them, and a method with a completeness rule
    then excludes the residue with that reason recorded.  To study another
    conformation, prepare a structure file whose alternate-location indicators
    select it (ADR-095).

    Hydrogen and deuterium are recognised from the element field alone (PDB
    columns 77-78, never the charge columns), compared without regard to
    letter case (``h`` is hydrogen).  An atom whose element field is blank is
    admitted here, so it is still *found*; the selector then refuses the
    calculation, because such an atom cannot be told from a hydrogen
    (``integrity.refuse_blank_elements``, ADR-108).  Atom names are the full
    PDB field, columns 13-16, for every admission (ADR-108).
    """

    name: str
    excluded_elements: Tuple[str, ...] = ("H", "D")
    accepted_alt_locs: Tuple[str, ...] = ("", "A")
    hetero_allowed: bool = False
    #: When true only ``HETATM`` records are eligible -- the ligand/carbohydrate
    #: rule. ``hetero_allowed`` alone would also admit ``ATOM`` records, which is
    #: not what either lineage does.
    hetero_required: bool = False
    #: Which stopping rule this lineage's single pass obeys.
    truncation: str = LEGACY_SINGLE_PASS
    #: ``HETATM`` residue names read as if they were ``ATOM`` records, by an
    #: admission that otherwise reads ``ATOM`` records only: the recognised
    #: modified amino acids (``completeness.MODIFIED_AMINO_ACIDS``), which RCSB
    #: files write as ``HETATM``.  A name, never a pattern: every other
    #: ``HETATM`` record stays excluded (ADR-101).
    hetero_residue_names: Tuple[str, ...] = ()

    def accepts(self, atom) -> Tuple[bool, str]:
        """Return ``(accepted, reason_if_not)``."""
        if self.hetero_required and not atom.hetero:
            return False, "atom_record_excluded"
        if atom.hetero and not (self.hetero_allowed or self.hetero_required
                                or atom.identity.residue_name in self.hetero_residue_names):
            return False, "hetatm_excluded"

        if (atom.element or "").strip().upper() in self.excluded_elements:
            return False, "hydrogen_or_deuterium"

        if atom.identity.alt_loc not in self.accepted_alt_locs:
            return False, "altloc_rejected"
        return True, ""


#: ``amino_acid_atomic``: H/D from the element field, and the standard
#: alternate-location rule (blank or ``A``).  Until 1.0.0rc5 this admission
#: read its hydrogen test from columns 78-80, which include the charge field,
#: so a hydrogen written ``H0`` or ``H1+`` was counted as a heavy atom; since
#: 1.0.0rc6 it reads columns 77-78 like every other method (ADR-108).  Until
#: 1.0.0rc3 it excluded only ``B`` conformers (ADR-095).
STANDARD_AMINOACID_RESIDUE_ADMISSION = AtomAdmission(
    name="standard_amino_acid_residue",
    excluded_elements=("H", "D"),
    hetero_residue_names=MODIFIED_AMINO_ACIDS,
)

#: ``WebsiteTSR``'s two residue selectors: the element field, an explicit
#: altLoc allow-list, ``ATOM`` records plus the recognised modified amino acids
#: (ADR-101), the peptide ``TER`` rule.
WEBSITE_TSR_ADMISSION = AtomAdmission(
    name="website_tsr",
    excluded_elements=("H", "D"),
    hetero_residue_names=MODIFIED_AMINO_ACIDS,
)

#: The standard admission for every HETATM entity family (drug, carbohydrate):
#: the element field, the altLoc allow-list, and the model-only stopping rule.
STANDARD_HETERO_ADMISSION = AtomAdmission(
    name="standard_hetero_entity",
    excluded_elements=("H", "D"),
    hetero_required=True,
    truncation=FIRST_MODEL_ONLY,
)


# ---------------------------------------------------------------------------
# Atom-selection request parsing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtomSelectionRequest:
    """One ``RESIDUE_SEQNUM_ATOM`` entry, e.g. ``LYS_48_NZ``."""

    residue_name: str
    residue_seq: str
    atom_name: str

    def token(self) -> str:
        return "%s_%s_%s" % (self.residue_name, self.residue_seq, self.atom_name)


#: A PDB residue identifier: the sequence number (``resSeq``, possibly
#: negative) followed by at most one insertion-code letter (``iCode``).
_RESIDUE_ID = re.compile(r"^(-?\d+)([A-Za-z]?)$")
#: A residue name in the compact ``NAMEID`` spelling: letters only, so the
#: boundary with the number is unambiguous.
_COMPACT_REF = re.compile(r"^([A-Za-z]{1,3})(-?\d+[A-Za-z]?)$")
#: A residue name in the ``NAME_ID`` spelling: the PDB's three columns.
_RESIDUE_NAME = re.compile(r"^[A-Za-z0-9]{1,3}$")


def parse_residue_id(value, what: str = "residue id") -> str:
    """Normalise a PDB residue identifier: ``"404"``, ``"52A"``, ``"-5"``.

    An identifier is the residue's own number in the file plus its insertion
    code, never an ordinal position ("the second copy").  Leading zeros are
    dropped (``"007"`` is residue 7); an insertion code is kept as written and
    matched exactly, so ``"52"`` selects residue 52 and never 52A.  Anything
    else -- a space, a separator, two letters -- is refused rather than
    reinterpreted.
    """
    text = str(value).strip()
    match = _RESIDUE_ID.match(text)
    if not match:
        raise InvalidRequest(
            "%s %r is not a PDB residue identifier: write the residue number "
            "as it appears in the structure file, optionally followed by one "
            "insertion-code letter (e.g. 404, 52A, -5)" % (what, value)
        )
    return "%d%s" % (int(match.group(1)), match.group(2))


def parse_atom_selection(value) -> Tuple[Tuple[AtomSelectionRequest, ...], Tuple[str, ...]]:
    """Parse ``aa_atom_selection`` into requests, returning ``(requests, warnings)``.

    ``LYS_48_NZ;HIS_68_NE2`` -> two requests.  Each token is
    ``RESIDUE_ID_ATOM``: a residue name, its PDB residue identifier
    (:func:`parse_residue_id`, so ``LYS_52A_NZ`` names an insertion-code
    residue) and an atom name.  Duplicates are dropped, preserving first-seen
    order.  A malformed token is **refused** with its spelling: losing an atom
    would change the triangle set, and guessing what was meant could change it
    differently.  ``warnings`` is always empty; it is kept so callers that
    report parse notes need no change.
    """
    if value is None:
        return (), ()
    if isinstance(value, (list, tuple)):
        tokens = [str(v) for v in value]
    else:
        tokens = [t for t in str(value).replace("\n", ";").split(";")]

    requests: List[AtomSelectionRequest] = []
    seen = set()
    for raw in tokens:
        token = raw.strip()
        if not token:
            continue
        parts = token.split("_", 2)
        if (len(parts) != 3 or not _RESIDUE_NAME.match(parts[0].strip())
                or not parts[2].strip()):
            raise InvalidRequest(
                "atom selection token %r is not RESIDUE_ID_ATOM: write the "
                "residue name, its residue number and the atom name separated "
                "by underscores, e.g. LYS_48_NZ (or LYS_52A_NZ for an "
                "insertion code)" % token
            )
        residue, seq, atom = (p.strip() for p in parts)
        seq = parse_residue_id(seq, "residue id in atom selection token %r" % token)
        request = AtomSelectionRequest(residue.upper(), seq, atom.upper())
        if request in seen:
            continue
        seen.add(request)
        requests.append(request)
    return tuple(requests), ()


# ---------------------------------------------------------------------------
# Shared machinery
# ---------------------------------------------------------------------------




def _admitted_chain_atoms(structure: Structure, chain_id, admission: AtomAdmission):
    """Chain atoms in file order, up to this lineage's own truncation point."""
    cutoff = truncation_line(structure, chain_id, admission.truncation)
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.chain != chain_id:
            continue
        yield atom, cutoff


def residue_targets(
    structure: Structure,
    chain: str,
    admission: AtomAdmission = STANDARD_AMINOACID_RESIDUE_ADMISSION,
) -> Tuple[Tuple[str, str], ...]:
    """Every ``(residue_seq_label, residue_name)`` the chain offers, in file order.

    Mirrors the legacy first pass: a residue is *offered* whenever any of its
    admitted atoms appears.  Whether it is *used* is the completeness policy's
    decision, made later and recorded.  The legacy pass keys a dict by sequence
    label, so a later residue with the same label would overwrite an earlier one;
    order here is first appearance, which is what the legacy dict preserves too.
    """
    chain_id = ChainId(chain)
    seen = {}
    for atom, _cutoff in _admitted_chain_atoms(structure, chain_id, admission):
        ok, _reason = admission.accepts(atom)
        if not ok:
            continue
        label = atom.identity.residue.legacy_seq_label
        seen[label] = atom.identity.residue_name
    return tuple(seen.items())



def hetero_amino_acid_residues(
    structure: Structure,
    chain: str,
    truncation: str = LEGACY_SINGLE_PASS,
) -> Tuple[Tuple[str, str], ...]:
    """``(seq_label, residue_name)`` of the chain's ``HETATM`` amino-acid residues.

    A ``HETATM`` residue counts when it carries the backbone atoms ``N``, a
    carbon ``CA`` and ``C`` -- a modified amino acid such as ``TPO``, ``CME`` or
    ``MSE`` in an RCSB file, never a calcium ion or an ordinary ligand.  Read in
    file order, up to the named truncation rule's cutoff.

    This **names** the residues a selection cannot represent, so a result can
    say so; it never admits an atom.  Which of them a method reads is its
    admission's ``hetero_residue_names`` (ADR-101).
    """
    chain_id = ChainId(chain)
    cutoff = truncation_line(structure, chain_id, truncation)
    names = {}
    order = []
    for atom in structure.atoms:
        if cutoff is not None and atom.file_order >= cutoff:
            continue
        if atom.chain != chain_id or not atom.hetero:
            continue
        key = (atom.identity.residue.legacy_seq_label, atom.identity.residue_name)
        if key not in names:
            names[key] = set()
            order.append(key)
        if atom.atom_name == "CA" and atom.element.strip().upper() != "C":
            continue
        names[key].add(atom.atom_name)
    return tuple(key for key in order if {"N", "CA", "C"} <= names[key])

def _label_point(atom, label_scheme, unlabelled_policy, skipped, selector):
    identity = atom.identity
    # The scheme decides how a point is looked up; the atomic scheme reads the
    # atom's element field and never its name (ADR-113).
    code = label_scheme.resolve_atom(identity.atom_name, atom.element)
    if code is None:
        if unlabelled_policy == UnlabelledPointPolicy.ERROR:
            raise InvalidRequest(
                "atom %s (residue %s %s, atom name %r, element %r) cannot be "
                "labelled: the element labels %s have no label for that element "
                "(atoms are labelled by their element field, PDB columns 77-78; "
                "the atom name is not used). A calculation is not computed "
                "without an atom it selected (ADR-108): correct the atom's "
                "element field if it is wrong, or narrow the selection to leave "
                "the atom out."
                % (identity.describe(), identity.residue_name,
                   identity.residue.legacy_seq_label, identity.atom_name,
                   (atom.element or "").strip(), label_scheme.name)
            )
        skipped.append(
            SkipRecord(
                stage="labelling",
                reason="unlabelled_point",
                detail="%s has no code in %s"
                % (identity.describe(), label_scheme.name),
            )
        )
        return None
    return LabeledPoint(identity=identity, coord=atom.coord, label=code)


def _context(structure: Structure, chain: str) -> str:
    return "chain %s of %s" % (chain, structure.structure_id.value)


def _check_records(structure: Structure, chain: str, records) -> None:
    """The input checks every atomic selection applies to the records it uses."""
    context = _context(structure, chain)
    integrity.refuse_blank_elements(records, context)
    integrity.refuse_duplicate_atoms(records, context)


def _record_type_ok(admission: AtomAdmission, record) -> bool:
    """Whether ``admission``'s record-type rule would read ``record``."""
    if admission.hetero_required and not record.hetero:
        return False
    if record.hetero and not (admission.hetero_allowed or admission.hetero_required
                              or record.residue_name in admission.hetero_residue_names):
        return False
    return True


def refuse_unreadable(structure: Structure, chain: str, admission: AtomAdmission,
                      names=None) -> None:
    """Refuse when a record this admission would read on ``chain`` has an
    unreadable residue number; ``names`` restricts it to those residue names."""
    chain_id = ChainId(chain)
    wanted = {str(n).strip().upper() for n in names} if names else None
    integrity.refuse_unreadable_numbers(
        structure, chain_id, truncation_line(structure, chain_id, admission.truncation),
        lambda r: _record_type_ok(admission, r)
        and (wanted is None or r.residue_name.upper() in wanted),
        _context(structure, chain))


# ---------------------------------------------------------------------------
# Per-entity scope: one residue's atoms
# ---------------------------------------------------------------------------


def select_residue_atoms(
    structure: Structure,
    chain: str,
    residue_seq: str,
    residue_name: str,
    label_scheme: LabelScheme,
    *,
    admission: AtomAdmission = STANDARD_AMINOACID_RESIDUE_ADMISSION,
    completeness: CompletenessPolicy = NO_COMPLETENESS_POLICY,
    truncation: str = LEGACY_SINGLE_PASS,
    unlabelled_policy: str = UnlabelledPointPolicy.ERROR,
) -> SelectionResult:
    """Select one residue's admitted heavy atoms as labelled points.

    ``residue_seq`` is the legacy sequence label -- number plus insertion code,
    e.g. ``"15"`` or ``"52A"`` -- matched against
    ``ResidueIdentity.legacy_seq_label``.  The residue **name must match too**,
    exactly as the legacy predicate requires: a sequence number alone does not
    identify a residue if the file disagrees about what is there.

    When the completeness policy rejects the residue the result carries **zero
    points** and a `SkipRecord` explaining why, rather than raising -- reproducing
    the legacy outcome (no files at all) while making the reason inspectable
    instead of invisible.
    """
    if truncation != LEGACY_SINGLE_PASS:
        raise ValueError(
            "unsupported truncation rule %r; this selector implements only %r"
            % (truncation, LEGACY_SINGLE_PASS)
        )

    chain_id = ChainId(chain)
    points: List[LabeledPoint] = []
    skipped: List[SkipRecord] = []
    cutoff_seen = -1

    refuse_unreadable(structure, chain, admission, names=(residue_name,))
    admitted = []
    for atom, cutoff in _admitted_chain_atoms(structure, chain_id, admission):
        cutoff_seen = -1 if cutoff is None else cutoff
        if atom.identity.residue.legacy_seq_label != residue_seq:
            continue
        if atom.identity.residue_name != residue_name:
            continue
        ok, reason = admission.accepts(atom)
        if not ok:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason=reason,
                    detail="%s (policy=%s)" % (atom.identity.describe(), admission.name),
                )
            )
            continue
        admitted.append(atom)
    # Every admitted atom counts towards the completeness rule, so every one is
    # checked before any is labelled.
    _check_records(structure, chain, admitted)
    for atom in admitted:
        point = _label_point(
            atom, label_scheme, unlabelled_policy, skipped, SELECTOR_RESIDUE_ATOMS
        )
        if point is not None:
            points.append(point)

    verdict: Admission = completeness.admit(residue_name, len(points))
    if not verdict.admitted:
        skipped.append(
            SkipRecord(
                stage="completeness",
                reason=verdict.reason or "rejected",
                detail="%s_%s: %s" % (residue_seq, residue_name, verdict.detail),
            )
        )
        points = []

    return SelectionResult(
        points=PointSet(points),
        selector=SELECTOR_RESIDUE_ATOMS,
        chains_requested=(chain,),
        chains_used=(chain,) if points else (),
        residues_used=("%s_%s" % (residue_seq, residue_name),) if points else (),
        skipped=tuple(skipped),
        truncated_at_line=cutoff_seen,
    )


# ---------------------------------------------------------------------------
# Pooled scope: explicitly named atoms across residues
# ---------------------------------------------------------------------------


def select_atom_selection(
    structure: Structure,
    chain: str,
    requests: Sequence[AtomSelectionRequest],
    label_scheme: LabelScheme,
    *,
    admission: AtomAdmission = WEBSITE_TSR_ADMISSION,
    truncation: str = LEGACY_SINGLE_PASS,
    unlabelled_policy: str = UnlabelledPointPolicy.ERROR,
) -> SelectionResult:
    """Select exactly the named atoms, pooled into one calculation.

    Points come out in **structure order**, not request order -- the legacy
    implementation iterates the PDB records and keeps the ones it wants, and that
    order is what fixes the candidate enumeration and therefore the output.

    Two things are reported rather than silently absorbed, because each changes
    which triangles exist:

    * a request whose residue **name** disagrees with the structure is skipped
      (the sequence number alone is not an identity);
    * a request never matched at all is listed in ``atoms_missing``.

    A request matched by more than one record -- the same atom written twice,
    or written twice differently -- is refused, naming both records; no
    occurrence is preferred (ADR-108).
    """
    if truncation != LEGACY_SINGLE_PASS:
        raise ValueError(
            "unsupported truncation rule %r; this selector implements only %r"
            % (truncation, LEGACY_SINGLE_PASS)
        )

    chain_id = ChainId(chain)
    wanted = {(r.residue_seq, r.atom_name): r for r in requests}
    points: List[LabeledPoint] = []
    skipped: List[SkipRecord] = []
    matched = set()
    cutoff_seen = -1

    refuse_unreadable(structure, chain, admission,
                      names={r.residue_name for r in requests})
    chosen = []
    for atom, cutoff in _admitted_chain_atoms(structure, chain_id, admission):
        cutoff_seen = -1 if cutoff is None else cutoff
        key = (atom.identity.residue.legacy_seq_label, atom.atom_name)
        if key not in wanted:
            continue
        ok, reason = admission.accepts(atom)
        if not ok:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason=reason,
                    detail="%s (policy=%s)" % (atom.identity.describe(), admission.name),
                )
            )
            continue
        request = wanted[key]
        if atom.identity.residue_name != request.residue_name:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="residue_name_mismatch",
                    detail="%s was requested but the structure has %s at that "
                    "position" % (request.token(), atom.identity.residue_name),
                )
            )
            continue
        chosen.append((key, atom))
    _check_records(structure, chain, [atom for _key, atom in chosen])
    for key, atom in chosen:
        point = _label_point(
            atom, label_scheme, unlabelled_policy, skipped, SELECTOR_ATOM_SELECTION
        )
        if point is not None:
            matched.add(key)
            points.append(point)

    missing = tuple(
        sorted(r.token() for k, r in wanted.items() if k not in matched)
    )
    for token in missing:
        skipped.append(
            SkipRecord(
                stage="selection",
                reason="requested_atom_not_found",
                detail="%s was requested but never matched in chain %s" % (token, chain),
            )
        )

    residues = []
    for p in points:
        label = "%s_%s" % (p.identity.residue.legacy_seq_label, p.identity.residue_name)
        if label not in residues:
            residues.append(label)

    return SelectionResult(
        points=PointSet(points),
        selector=SELECTOR_ATOM_SELECTION,
        chains_requested=(chain,),
        chains_used=(chain,) if points else (),
        residues_used=tuple(residues),
        atoms_missing=missing,
        skipped=tuple(skipped),
        warnings=(),
        truncated_at_line=cutoff_seen,
    )


# ---------------------------------------------------------------------------
# Entity-scoped selection -- the shared foundation for every family
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResidueRef:
    """One explicitly named residue instance: ``LYS_6`` or ``GLU18``.

    The historical dataset sheet's ``aa_ind_specific`` column really does use both
    spellings in the same cell (``"LYS_6;GLU18"``), so both are parsed and
    neither is normalized away.
    """

    residue_name: str
    residue_seq: str

    def token(self) -> str:
        return "%s_%s" % (self.residue_name, self.residue_seq)


def parse_residue_types(value) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """Parse ``aa_ind`` / ``aa_ind_mix`` into residue-type names.

    ``"TYR"`` -> ``("TYR",)``;  ``"TYR;PHE"`` -> ``("TYR", "PHE")``.
    Duplicates are dropped preserving order; nothing is validated against a
    residue vocabulary here, because a modified residue such as ``TPO`` is a
    legitimate type and the completeness policy is what decides admission.
    """
    if value is None:
        return (), ()
    if isinstance(value, (list, tuple)):
        tokens = [str(v) for v in value]
    else:
        tokens = str(value).replace(",", ";").replace("\n", ";").split(";")
    types: List[str] = []
    warnings: List[str] = []
    for raw in tokens:
        token = raw.strip().upper()
        if not token:
            continue
        if token in types:
            continue
        types.append(token)
    return tuple(types), tuple(warnings)


def parse_residue_refs(value) -> Tuple[Tuple[ResidueRef, ...], Tuple[str, ...]]:
    """Parse ``aa_ind_specific`` into explicit residue references.

    Two spellings are supported, because the historical dataset sheets use both
    in one cell (``"LYS_6;GLU18"``)::

        LYS_6     ->  ResidueRef("LYS", "6")     NAME_ID: any 1-3 character name
        GLU18     ->  ResidueRef("GLU", "18")    NAMEID: letters-only name
        GLU18A    ->  ResidueRef("GLU", "18A")   insertion code kept

    The identifier is the residue's PDB number (:func:`parse_residue_id`), never
    a position.  A residue name containing a digit must use the ``NAME_ID``
    spelling, since ``NAMEID`` would not say where the name ends; any token
    matching neither spelling is **refused**, naming it, rather than skipped or
    reinterpreted.  The two spellings of one residue are one residue.
    ``warnings`` is always empty.
    """
    if value is None:
        return (), ()
    if isinstance(value, (list, tuple)):
        tokens = [str(v) for v in value]
    else:
        tokens = str(value).replace(",", ";").replace("\n", ";").split(";")

    refs: List[ResidueRef] = []
    for raw in tokens:
        token = raw.strip()
        if not token:
            continue
        if "_" in token:
            name, _, seq = token.partition("_")
            name = name.strip()
            if not _RESIDUE_NAME.match(name):
                name = ""
        else:
            match = _COMPACT_REF.match(token)
            name, seq = (match.group(1), match.group(2)) if match else ("", "")
        if not name or not seq.strip():
            raise InvalidRequest(
                "residue token %r is neither NAME_ID (e.g. LYS_6) nor NAMEID "
                "(e.g. GLU18): name one residue by its three-letter name and "
                "its residue number in the structure file; a name containing a "
                "digit needs the NAME_ID spelling" % token
            )
        seq = parse_residue_id(seq, "residue id in token %r" % token)
        ref = ResidueRef(name.upper(), seq)
        if ref in refs:
            continue
        refs.append(ref)
    return tuple(refs), ()


def parse_atom_names(value) -> Tuple[str, ...]:
    """Parse ``drug_atom_selected`` / ``carb_atom_selected`` / ``nt_atom_selected``.

    ``"N5;C6;O7"`` -> ``("N5", "C6", "O7")``.  Commas work as well as semicolons,
    matching ``WebsiteTSR``'s own ``split_on_semicolons``.  Case is preserved
    exactly: a PDB atom name is upper-case as written and is matched as written
    (``CA`` and ``Ca`` are different names).
    """
    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        tokens = [str(v) for v in value]
    else:
        tokens = str(value).replace(",", ";").replace("\n", ";").split(";")
    names: List[str] = []
    for raw in tokens:
        token = raw.strip()
        if token and token not in names:
            names.append(token)
    return tuple(names)


def entity_targets(
    structure: Structure,
    chain: str,
    admission: AtomAdmission,
    *,
    entity_names: Optional[Sequence[str]] = None,
) -> Tuple[Tuple[str, str], ...]:
    """Every ``(seq_label, entity_name)`` the chain offers under this admission.

    The generalization of :func:`residue_targets` to any family: a nucleotide
    residue, a ligand copy and a carbohydrate monomer are all "one entity
    instance", and which records are eligible is the admission policy's business,
    not this function's.
    """
    chain_id = ChainId(chain)
    wanted = {n.strip().upper() for n in entity_names} if entity_names else None
    seen = {}
    for atom, _cutoff in _admitted_chain_atoms(structure, chain_id, admission):
        ok, _reason = admission.accepts(atom)
        if not ok:
            continue
        if wanted is not None and atom.identity.residue_name.upper() not in wanted:
            continue
        label = atom.identity.residue.legacy_seq_label
        seen[label] = atom.identity.residue_name
    return tuple(seen.items())


def _entity_key(atom) -> Tuple[str, str, str]:
    """One entity *instance*. The chain is part of it.

    Without the chain, two residues with the same name and sequence number in
    different chains would share a group -- which is exactly the class of identity
    collapse this package exists to prevent, and it becomes reachable as soon as
    a selection scope spans chains.
    """
    return (
        atom.identity.residue_name,
        atom.identity.residue.legacy_seq_label,
        atom.identity.chain.value,
    )


def select_entity_atoms(
    structure: Structure,
    chain: str,
    label_scheme: LabelScheme,
    *,
    admission: AtomAdmission,
    entity_names: Optional[Sequence[str]] = None,
    entity_seqs: Optional[Sequence[str]] = None,
    entity_refs: Optional[Sequence[ResidueRef]] = None,
    atom_names: Optional[Sequence[str]] = None,
    completeness: CompletenessPolicy = NO_COMPLETENESS_POLICY,
    selector: str = SELECTOR_ENTITY_ATOMS,
    unlabelled_policy: str = UnlabelledPointPolicy.ERROR,
    require_entity_name: bool = True,
) -> SelectionResult:
    """Select the atoms of one or more named entity instances, pooled.

    This is the **one** selector every non-Cα family configures rather than
    copies.  It expresses each of the legacy scopes as a combination of filters:

    ===========================================  ==========================================
    ``entity_names={"LEU"}, entity_seqs={"15"}``  one residue (``TSR_AMINOACID_PACKAGE``)
    ``entity_names={"TYR"}``                      ``aa_ind`` -- every TYR of the chain
    ``entity_names={"TYR","TPO"}``                ``aa_ind_mix``
    ``entity_refs=[LYS_6, LYS_48]``             ``aa_ind_specific``
    ``entity_names={"ESW"}, entity_seqs={"404"}`` one drug copy (``WebsiteTSR`` ligand)
    ``… + atom_names={"N4","C5",…}``              ``drug_atom_selected``
    ``entity_names=CARBOHYDRATE_NAMES``           the carbohydrate polymer scope
    ===========================================  ==========================================

    **Completeness is applied per entity instance**, against that instance's own
    admitted, labelled atom count -- which is exactly what
    ``CarbohydratesTSR`` does (``_load_atoms_for_selection`` groups by
    ``(name, seq)`` and drops incomplete groups individually), and what the
    single-residue amino-acid and nucleotide analyzers do for their one target.
    A rejected instance contributes nothing and leaves a ``SkipRecord`` saying
    why; it never aborts the others.

    **Atom ordering** is grouped by entity instance in first-appearance order,
    then file order within the instance.  For a normal PDB, where an entity's
    records are contiguous, that *is* file order; stating it explicitly is what
    keeps the order deterministic when they are not.

    ``atom_names`` narrows an entity to a named subset.  Following
    ``WebsiteTSR``'s ligand selector, an **empty or absent** list means *every*
    admitted heavy atom -- the regular, whole-entity calculation -- and a name
    requested but not present is reported in ``atoms_missing``, never ignored.
    """
    chain_id = ChainId(chain)
    name_filter = (
        {n.strip().upper() for n in entity_names} if entity_names else None
    )
    seq_filter = {str(s).strip() for s in entity_seqs} if entity_seqs else None
    ref_filter = (
        {(r.residue_name.upper(), r.residue_seq) for r in entity_refs}
        if entity_refs
        else None
    )
    atom_filter = {a for a in atom_names} if atom_names else None

    groups = {}
    order = []
    skipped: List[SkipRecord] = []
    cutoff_seen = -1
    matched_atom_names = set()

    names_read = (name_filter if name_filter is not None
                  else {r[0] for r in ref_filter} if ref_filter is not None else None)
    refuse_unreadable(structure, chain, admission, names=names_read)
    for atom, cutoff in _admitted_chain_atoms(structure, chain_id, admission):
        cutoff_seen = -1 if cutoff is None else cutoff
        residue_name = atom.identity.residue_name
        seq_label = atom.identity.residue.legacy_seq_label
        if name_filter is not None and residue_name.upper() not in name_filter:
            continue
        if seq_filter is not None and seq_label not in seq_filter:
            continue
        if ref_filter is not None and (residue_name.upper(), seq_label) not in ref_filter:
            continue
        ok, reason = admission.accepts(atom)
        if not ok:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason=reason,
                    detail="%s (policy=%s)" % (atom.identity.describe(), admission.name),
                )
            )
            continue
        key = _entity_key(atom)
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(atom)

    # The records the result depends on: with a completeness rule, every
    # admitted atom of every instance (each one counts); without one, only the
    # atoms that become points. Those are checked before any is labelled, and
    # an atom outside them is neither checked nor labelled (ADR-108).
    counted = completeness is not NO_COMPLETENESS_POLICY
    if not counted and atom_filter is not None:
        for key in order:
            groups[key] = [a for a in groups[key] if a.atom_name in atom_filter]
    _check_records(structure, chain, [a for key in order for a in groups[key]])

    points: List[LabeledPoint] = []
    residues_used: List[str] = []
    #: Admitted instances whose atoms differ from their type's standard set,
    #: when the completeness policy describes them (ADR-114).
    atom_set_notes: List[str] = []
    describe = getattr(completeness, "differences", None)
    #: ``(NAME_ID, atoms requested but absent)`` per admitted instance, when an
    #: atom list narrows the instances -- the list applies to each of them.
    missing_by_instance: List[Tuple[str, Tuple[str, ...]]] = []
    for key in order:
        entity_name, seq_label = key[0], key[1]
        labelled: List[LabeledPoint] = []
        for atom in groups[key]:
            point = _label_point(
                atom, label_scheme, unlabelled_policy, skipped, selector
            )
            if point is None:
                continue
            labelled.append(point)

        verdict: Admission = completeness.admit(entity_name, len(labelled))
        if not verdict.admitted:
            skipped.append(
                SkipRecord(
                    stage="completeness",
                    reason=verdict.reason or "rejected",
                    detail="%s_%s: %s" % (seq_label, entity_name, verdict.detail),
                )
            )
            continue
        if describe is not None:
            lacking, extra = describe(entity_name, [p.identity.atom_name for p in labelled])
            if lacking or extra:
                atom_set_notes.append("%s_%s %s" % (entity_name, seq_label, "; ".join(
                    filter(None, ("lacks %s" % ", ".join(lacking) if lacking else "",
                                  "has %s in addition" % ", ".join(extra) if extra else "")))))

        # Atom narrowing happens AFTER completeness, never before: asking for
        # four atoms of a fourteen-atom monomer must not make it look incomplete.
        if atom_filter is not None:
            kept = []
            present = set()
            for point in labelled:
                if point.identity.atom_name in atom_filter:
                    matched_atom_names.add(point.identity.atom_name)
                    present.add(point.identity.atom_name)
                    kept.append(point)
            lacking = tuple(a for a in atom_names if a not in present)
            if lacking:
                missing_by_instance.append(("%s_%s" % (entity_name, seq_label), lacking))
            labelled = kept
            if not labelled:
                continue

        points.extend(labelled)
        residues_used.append("%s_%s" % (seq_label, entity_name))

    missing: Tuple[str, ...] = ()
    if atom_filter is not None and len(order) > 1 and missing_by_instance:
        # Several instances, one atom list: what is missing is said per
        # instance, so "C8 is absent from NAG_5" is not reported as "C8 is
        # absent" when NAG_1 and NAG_2 have it.
        missing = tuple(
            "%s in %s" % (name, instance)
            for instance, lacking in missing_by_instance for name in lacking
        )
        for token in missing:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="requested_atom_not_found",
                    detail="atom %s was requested but is not present (chain %s)"
                    % (token, chain),
                )
            )
    elif atom_filter is not None:
        missing = tuple(sorted(atom_filter - matched_atom_names))
        for name in missing:
            skipped.append(
                SkipRecord(
                    stage="selection",
                    reason="requested_atom_not_found",
                    detail="atom %r was requested but never matched in chain %s"
                    % (name, chain),
                )
            )

    unmatched_entities: List[str] = []
    if require_entity_name and ref_filter is not None:
        # `residues_used` entries are "<seq>_<NAME>"; compare on (NAME, seq).
        found = set()
        for entity in residues_used:
            seq_label, _, entity_name = entity.partition("_")
            found.add((entity_name.upper(), seq_label))
        for ref in entity_refs or ():
            if (ref.residue_name.upper(), ref.residue_seq) not in found:
                unmatched_entities.append(ref.token())
                skipped.append(
                    SkipRecord(
                        stage="selection",
                        reason="requested_entity_not_found",
                        detail="%s was requested but no admitted instance was found "
                        "in chain %s" % (ref.token(), chain),
                    )
                )

    chains_requested = (chain,)
    chains_used = (chain,) if points else ()

    return SelectionResult(
        points=PointSet(points),
        selector=selector,
        chains_requested=chains_requested,
        chains_used=chains_used,
        residues_used=tuple(residues_used),
        atoms_missing=missing + tuple(unmatched_entities),
        skipped=tuple(skipped),
        warnings=(ATOM_SET_NOTE % "; ".join(atom_set_notes),) if atom_set_notes else (),
        truncated_at_line=cutoff_seen,
        atoms_missing_by_instance=tuple(missing_by_instance),
    )
