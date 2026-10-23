"""Selection protocol and result types.

Selection decides *which points participate and how they are grouped into one
calculation*.  It is structural/identity work and may never consult triangle
geometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

from ..domain.point import PointSet
from ..domain.result import SkipRecord
from ..labels.schemes import MODIFIED_RESIDUE_CODES
from .completeness import MODIFIED_AMINO_ACIDS

__all__ = [
    "SelectionResult",
    "AltLocPolicy",
    "CALPHA_MODIFIED_RESIDUES",
    "ModifiedResiduePolicy",
    "InputRules",
    "SeqLabelPolicy",
    "processing_identity",
]


class AltLocPolicy:
    """The alternate-location policy, recorded in provenance.

    One rule for every method (ADR-095): keep atoms whose altLoc indicator is
    blank or ``A``; exclude ``B``, ``C`` and every later indicator.  There is no
    fallback to another conformer when ``A`` is absent.
    """

    #: ``line[16] == 'A' or line[16] == ' '``.
    BLANK_OR_A = "blank_or_a"


#: The modified residues the Ca methods read and label (ADR-107): every residue
#: of the residue label table other than the 20 regular amino acids.  RCSB
#: files write them as ``HETATM`` records, so the Ca selectors accept a
#: ``HETATM`` Ca of exactly these names when its element field is carbon or
#: blank, and an ``ATOM`` record of them as they accept any ``ATOM`` record.
#: Alphabetical.
CALPHA_MODIFIED_RESIDUES: Tuple[str, ...] = tuple(
    sorted(name for name, _code in MODIFIED_RESIDUE_CODES))


class ModifiedResiduePolicy:
    """Which modified amino acids written as ``HETATM`` a protein selection reads.

    Recorded in every result's provenance (``source.modified_residue_policy``)
    and in every completion record, because it decides which atoms exist for a
    method that reads protein atoms, and it changed without the key space
    changing (ADR-101, ADR-105).  It has one value per **label domain**, because
    the two kinds of method read different residues (ADR-107):

    * ``ATOMIC`` -- every method that labels atoms (the atomic protein
      selections and the protein sides of the atomic cross methods) reads TPO,
      SEP and PTR, the residues with defined atom counts;
    * ``CALPHA`` -- every method that labels residues (``protein``,
      ``aa_selection``, ``protein_protein_calpha_cross``) reads the supported
      modified residues of the residue label table.

    Every other ``HETATM`` residue stays outside a protein selection.  A result
    or record without the field was written before it existed -- under a rule
    that read no ``HETATM`` amino acid, or (briefly, before the rule was
    recorded) under the atomic one; the files cannot tell which.  It is
    reported as *not recorded*, never taken to be a current rule.

    A Ca result written before 1.0.0rc5 records ``ATOMIC``'s value, which was
    then the one rule for every method although the Ca methods read none of
    those residues; it differs from ``CALPHA``, so such a result is never
    reused for a Ca calculation now (and its key space differs too).
    """

    #: ``hetatm_amino_acids:PTR,SEP,TPO`` -- derived from
    #: ``completeness.MODIFIED_AMINO_ACIDS``, so the record cannot drift from
    #: the allow-list the atomic admissions apply.
    ATOMIC = "hetatm_amino_acids:" + ",".join(sorted(MODIFIED_AMINO_ACIDS))
    #: ``hetatm_amino_acids:CME,...,TRF;ca_element:C_or_blank`` -- the residue
    #: names derived from the residue label table (``CALPHA_MODIFIED_RESIDUES``),
    #: and the rule a ``HETATM`` Ca of one of them is read by: its own element
    #: field carbon or blank, nothing else of the residue required.  The suffix
    #: tells these results apart from those of an earlier 1.0.0rc5 build that
    #: also required the backbone atoms N and C, which recorded the names alone
    #: (ADR-107).
    CALPHA = ("hetatm_amino_acids:" + ",".join(CALPHA_MODIFIED_RESIDUES)
              + ";ca_element:C_or_blank")

    @classmethod
    def for_label_domain(cls, domain: str) -> str:
        """The rule a method whose label scheme is in ``domain`` applies."""
        from ..labels.schemes import LabelDomain

        return cls.CALPHA if domain == LabelDomain.PROTEIN_AMINO_ACID else cls.ATOMIC


class InputRules:
    """How structure records are read and checked before any point is chosen.

    Recorded in every result (``source.input_rules``) and completion record,
    because these rules changed in 1.0.0rc6 without the key space changing:
    the same request can select different atoms, or be refused, under them
    (ADR-108, ADR-109, ADR-110).  A result or record without the field was
    written by 1.0.0rc5 or earlier: *not recorded*, never read as the current
    rules, so it is not reused on resume and not mixed silently in analysis.
    """

    #: Atom names from PDB columns 13-16 (the Calpha test keeps columns 14-16);
    #: element from columns 77-78, case-insensitive; refusal of blank element
    #: fields, repeated, conflicting and coincident atoms, unreadable residue
    #: numbers, UNK and unlabelled atoms in a selection; protein sides stop at
    #: their chain's TER; the narrowed Calpha cross interface; SSE endpoints
    #: matched with their insertion codes.
    CURRENT = "input_rules_v2"


def processing_identity(label_domain: str, completeness_policy=None):
    """The processing rules a calculation is drawn under, for its label domain.

    The alternate-location rule (ADR-095) is the same for every method; the
    modified-residue rule depends on whether the method labels residues or
    atoms (ADR-101, ADR-107).  A completion record stores it, so
    ``--if-exists skip`` keeps an existing calculation only when it was computed
    under the same rules -- the package version alone cannot tell, since one
    version number has covered more than one rule set (ADR-105).
    ``biotsr.api.processing_identity`` takes a method name instead.

    A method with a completeness rule records it too (``completeness_policy``,
    since 1.0.0rc7): the rule decides which molecules take part with the atoms
    they have, and ``nucleotide_chain``'s changed in 1.0.0rc7 (ADR-114).  A
    method without one -- every Calpha method -- records exactly what it did
    before.
    """
    rules = {
        "altloc_policy": AltLocPolicy.BLANK_OR_A,
        "modified_residue_policy": ModifiedResiduePolicy.for_label_domain(label_domain),
        "input_rules": InputRules.CURRENT,
    }
    if completeness_policy is not None:
        rules["completeness_policy"] = completeness_policy
    return rules


def altloc_policy_of(admission) -> str:
    """The :class:`AltLocPolicy` an atom admission applies, for provenance.

    Derived from the admission's own fields, so the recorded policy is the one
    that selected the atoms; an admission matching no named policy is an error
    rather than a guess.
    """
    accepted = getattr(admission, "accepted_alt_locs", None)
    if accepted is not None and tuple(accepted) == ("", "A"):
        return AltLocPolicy.BLANK_OR_A
    raise ValueError("atom admission %r applies no named altLoc policy"
                     % getattr(admission, "name", admission))


class SeqLabelPolicy:
    """Which residues a Ca selector will admit, by the shape of their sequence label.

    A selection-side policy, so it enters ``SamplingProfile``: it changes *which*
    points take part, never what a key means.

    ``ACCEPT_ALL`` is the standard and what ``TSR.py`` does -- every residue is
    admitted whatever its sequence label looks like.

    (A historical SSE program also refused every residue with an insertion
    code; that rule is not reproduced.)
    """

    ACCEPT_ALL = "accept_all"


@dataclass(frozen=True)
class SelectionResult:
    """Admitted points plus a full account of what was not admitted."""

    points: PointSet
    selector: str
    chains_requested: Tuple[str, ...]
    chains_used: Tuple[str, ...]
    residues_used: Tuple[str, ...] = ()
    atoms_missing: Tuple[str, ...] = ()
    skipped: Tuple[SkipRecord, ...] = ()
    warnings: Tuple[str, ...] = ()
    truncated_at_line: int = -1
    #: ``(NAME_ID, atoms)`` per selected instance lacking requested atoms, when
    #: an atom list narrows one or more molecule instances.
    atoms_missing_by_instance: Tuple[Tuple[str, Tuple[str, ...]], ...] = ()
