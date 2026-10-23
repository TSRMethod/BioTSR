"""Entity-admission (completeness) policies.

A completeness policy answers one question: *given the atoms an entity actually
contributed, does this entity take part at all?*  It is a separate, swappable
object rather than a check embedded in the geometry engine, because the legacy
family implements three genuinely different behaviours -- always-on
(amino-acid/nucleotide atomic), mode-dependent (carbohydrate), and absent
(protein Cα, `WebsiteTSR`) -- and the architecture must be able to express all
three under one interface (ADR-004).

The legacy exact-atom-count behaviour is **intentional** (ADR-004), so for each
method that has one the standard policy is that method's own legacy rule,
reproduced exactly.  The one chain-scope method that had a count rule,
``nucleotide_chain``, now admits each nucleotide with the atoms it has and names
every nucleotide whose atoms differ from the standard set (ADR-114), as the
chain-only rule for terminal residues requires (ADR-111).  Nothing here
special-cases a terminal residue, repairs a missing atom or infers a
coordinate.  Rejection is always recorded with a reason; it is never a silent
drop.

Note the rule is **exact equality, not a minimum**.  A residue carrying one atom
too many is rejected exactly as firmly as one missing an atom -- which is how
1UBQ's C-terminal Gly76 (4 expected, 5 present because of `OXT`) comes to be
excluded.  Reading it as "at least N" would silently admit 76 residues where
history admits 75.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

__all__ = [
    "Admission",
    "CompletenessPolicy",
    "ExactAtomCountPolicy",
    "NO_COMPLETENESS_POLICY",
    "LEGACY_AMINOACID_ATOM_COUNTS",
    "MODIFIED_AMINO_ACIDS",
    "EXACT_LEGACY_ATOM_COUNT",
    "LEGACY_NUCLEOTIDE_ATOM_COUNTS",
    "EXACT_LEGACY_NUCLEOTIDE_ATOM_COUNT",
    "STANDARD_NUCLEOTIDE_ATOM_NAMES",
    "AvailableAtomsPolicy",
    "AVAILABLE_NUCLEOTIDE_ATOMS",
    "DNA_RESIDUE_NAMES",
    "RNA_RESIDUE_NAMES",
    "NUCLEOTIDE_RESIDUE_NAMES",
    "LEGACY_CARBOHYDRATE_ATOM_COUNTS",
    "EXACT_LEGACY_CARBOHYDRATE_ATOM_COUNT",
    "CARBOHYDRATE_RESIDUE_NAMES",
]


@dataclass(frozen=True)
class Admission:
    """Whether an entity participates, and -- when it does not -- why."""

    admitted: bool
    reason: Optional[str] = None
    detail: str = ""


class CompletenessPolicy:
    """The protocol: ``admit(entity_name, n_atoms) -> Admission``."""

    name = "abstract"

    def admit(self, entity_name: str, n_atoms: int) -> Admission:  # pragma: no cover
        raise NotImplementedError


@dataclass(frozen=True)
class _NoCompleteness(CompletenessPolicy):
    """Every entity participates. Used where the legacy method has no such rule."""

    name: str = "none"

    def admit(self, entity_name: str, n_atoms: int) -> Admission:
        return Admission(admitted=True)


#: Protein Cα and `WebsiteTSR`'s selectors never implemented a completeness
#: concept, so there is nothing to default to and none is invented.
NO_COMPLETENESS_POLICY = _NoCompleteness()


#: The exact atom counts from ``TSR_AMINOACID_PACKAGE/.../analyzer.py:118``,
#: transcribed verbatim. Twenty standard residues plus three phosphorylated ones.
LEGACY_AMINOACID_ATOM_COUNTS: Mapping[str, int] = {
    "GLY": 4, "ALA": 5, "VAL": 7, "LEU": 8,
    "ILE": 8, "PRO": 7, "PHE": 11, "TYR": 12,
    "TRP": 14, "SER": 6, "THR": 7, "MET": 8,
    "CYS": 6, "HIS": 10, "LYS": 9, "ARG": 11,
    "ASP": 8, "GLU": 9, "ASN": 8, "GLN": 9,
    "TPO": 11, "SEP": 10, "PTR": 16,
}

#: The modified amino acids BioTSR recognises: the three phosphorylated residues
#: this table gives exact heavy-atom counts for (phosphothreonine, phosphoserine,
#: phosphotyrosine), which the historical cross scripts also list as protein
#: residues (``cross/sides.py``).  Their atoms are labelled by element like every
#: other atom (``P``, ``O1P``..``O3P``, ``OG1``, ``OH`` ...).  RCSB files write them
#: as ``HETATM`` records, so the amino-acid admissions accept ``HETATM`` records
#: of exactly these names (ADR-101); no other modified residue is read as an
#: amino acid by the atomic methods, because none has a defined atom count.  The
#: Calpha methods label residues, not atoms, and read a longer list -- the
#: supported modified residues of the residue label table
#: (``selection.base.CALPHA_MODIFIED_RESIDUES``, ADR-107); that list is never
#: used to decide atomic completeness.
MODIFIED_AMINO_ACIDS: Tuple[str, ...] = ("TPO", "SEP", "PTR")


@dataclass(frozen=True)
class ExactAtomCountPolicy(CompletenessPolicy):
    """Admit an entity only when its atom count exactly matches the table.

    An entity whose name is absent from the table is rejected too, whatever its
    atom count -- which is how the legacy implementation treats every
    non-standard residue, and why a modified residue silently produces no output.
    """

    expected: Mapping[str, int]
    name: str = "exact_legacy_atom_count"

    def admit(self, entity_name: str, n_atoms: int) -> Admission:
        if entity_name not in self.expected:
            return Admission(
                admitted=False,
                reason="unlisted_entity_type",
                detail=(
                    "%s is not one of the %d residue types the completeness table "
                    "covers" % (entity_name, len(self.expected))
                ),
            )
        wanted = self.expected[entity_name]
        if n_atoms != wanted:
            return Admission(
                admitted=False,
                reason="atom_count_mismatch",
                detail="%s: expected exactly %d atoms, got %d"
                % (entity_name, wanted, n_atoms),
            )
        return Admission(admitted=True)


#: The amino-acid atomic family's own legacy rule, which is also its standard
#: (ADR-004: the behaviour is intentional).
EXACT_LEGACY_ATOM_COUNT = ExactAtomCountPolicy(expected=LEGACY_AMINOACID_ATOM_COUNTS)


# ---------------------------------------------------------------------------
# Nucleotide
# ---------------------------------------------------------------------------

#: Verbatim from ``TSR_NUCLEOTIDE_PACKAGE/nucleotide_tsr_package/analyzer.py:114``.
#: The same shape as the amino-acid table and the same exact-equality rule --
#: read from that source, not inferred from the sibling package.
LEGACY_NUCLEOTIDE_ATOM_COUNTS: Mapping[str, int] = {
    "G": 23, "A": 22, "C": 20, "U": 20,
    "DG": 22, "DA": 21, "DC": 19, "DT": 20,
}

#: The four deoxyribonucleotides the legacy table covers.
DNA_RESIDUE_NAMES = ("DA", "DC", "DG", "DT")
#: The four ribonucleotides the legacy table covers.
RNA_RESIDUE_NAMES = ("A", "C", "G", "U")
#: Both, which is the whole vocabulary the legacy analyzer can admit.
NUCLEOTIDE_RESIDUE_NAMES = tuple(sorted(DNA_RESIDUE_NAMES + RNA_RESIDUE_NAMES))

#: The nucleotide family's own legacy rule, which is also its standard (ADR-004)
#: for the residue-specific methods (``nucleotide``,
#: ``nucleotide_atom_selection``).
EXACT_LEGACY_NUCLEOTIDE_ATOM_COUNT = ExactAtomCountPolicy(
    expected=LEGACY_NUCLEOTIDE_ATOM_COUNTS,
    name="exact_legacy_nucleotide_atom_count",
)

_PHOSPHATE_SUGAR = ("P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
                    "C2'", "C1'")

#: The heavy-atom names of the eight standard nucleotides, as the PDB's
#: standard residue definitions spell them, without the 5'-terminal ``OP3``.
#: Each list has exactly the legacy count above (DA 21, DC 19, DG 22, DT 20,
#: A 22, C 20, G 23, U 20), which an import-time check confirms.  Used only to
#: *describe* a nucleotide's atoms in a chain calculation (ADR-114) -- never to
#: admit, reject, rename or repair one.
STANDARD_NUCLEOTIDE_ATOM_NAMES: Mapping[str, Tuple[str, ...]] = {
    "DA": _PHOSPHATE_SUGAR + ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"),
    "DC": _PHOSPHATE_SUGAR + ("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
    "DG": _PHOSPHATE_SUGAR + ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2",
                              "N3", "C4"),
    "DT": _PHOSPHATE_SUGAR + ("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6"),
    "A": _PHOSPHATE_SUGAR + ("O2'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2",
                             "N3", "C4"),
    "C": _PHOSPHATE_SUGAR + ("O2'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"),
    "G": _PHOSPHATE_SUGAR + ("O2'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2",
                             "N2", "N3", "C4"),
    "U": _PHOSPHATE_SUGAR + ("O2'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"),
}
assert {k: len(set(v)) for k, v in STANDARD_NUCLEOTIDE_ATOM_NAMES.items()} == dict(
    LEGACY_NUCLEOTIDE_ATOM_COUNTS)


@dataclass(frozen=True)
class AvailableAtomsPolicy(CompletenessPolicy):
    """Admit every instance of a listed type with the atoms it has.

    No atom count is required: an instance takes part with every admitted
    heavy atom it carries, at a chain end or inside the chain.  A type not in
    ``expected_atoms`` is still rejected (``unlisted_entity_type``), as under
    the exact-count rule, so the set of supported residue types is unchanged.

    :meth:`differences` names the atoms an instance lacks or carries beyond
    its type's standard set.  That is a description of the structure file --
    the calculation uses every atom present -- so it is reported as
    information and never makes a result partial.  No coordinate is invented,
    no atom is reconstructed and no chain end is detected.
    """

    expected_atoms: Mapping[str, Tuple[str, ...]]
    name: str = "available_atoms"

    def admit(self, entity_name: str, n_atoms: int) -> Admission:
        if entity_name not in self.expected_atoms:
            return Admission(
                admitted=False,
                reason="unlisted_entity_type",
                detail=(
                    "%s is not one of the %d residue types this method reads"
                    % (entity_name, len(self.expected_atoms))
                ),
            )
        return Admission(admitted=True)

    def differences(self, entity_name: str, atom_names) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
        """``(missing, extra)`` atom names of one instance against its standard set.

        ``missing`` follows the standard order; ``extra`` the order given.
        """
        expected = self.expected_atoms.get(entity_name, ())
        present = list(atom_names)
        missing = tuple(a for a in expected if a not in present)
        extra = tuple(a for a in present if a not in expected)
        return missing, extra


#: ``nucleotide_chain``: every available heavy atom of each standard nucleotide
#: of the chain, with every nucleotide whose atom set differs from the
#: standard named (ADR-114).
AVAILABLE_NUCLEOTIDE_ATOMS = AvailableAtomsPolicy(
    expected_atoms=STANDARD_NUCLEOTIDE_ATOM_NAMES,
    name="available_nucleotide_atoms",
)


# ---------------------------------------------------------------------------
# Carbohydrate
# ---------------------------------------------------------------------------

#: Verbatim from ``CarbohydratesTSR/key_generator/key_generation_using_id_chain_res_seq/
#: expected_atom_counts.json``'s ``"ligand"`` section.  This file does double duty
#: in the legacy script: it is both the **completeness table** for monomer scope
#: and the **vocabulary of what counts as a carbohydrate** for polymer scope
#: (``process_target`` passes ``sorted(ligand_atom_counts)`` as the chain-scope
#: name filter).  Both uses are preserved and named separately here so neither is
#: an accident of the other.
LEGACY_CARBOHYDRATE_ATOM_COUNTS: Mapping[str, int] = {
    "BGC": 11, "BMA": 11, "FUC": 10, "FUL": 10,
    "GAL": 11, "MAN": 11, "NAG": 14, "NDG": 14,
}

#: The carbohydrate vocabulary, which is the key set of the table above.
CARBOHYDRATE_RESIDUE_NAMES = tuple(sorted(LEGACY_CARBOHYDRATE_ATOM_COUNTS))

#: Monomer scope only.  The legacy polymer (``id_chain``) path passes
#: ``apply_ligand_atom_count_filter=False`` and admits every instance of a
#: configured name whatever its atom count -- so the polymer profile carries
#: ``NO_COMPLETENESS_POLICY``, which is a characterized fact, not an omission.
EXACT_LEGACY_CARBOHYDRATE_ATOM_COUNT = ExactAtomCountPolicy(
    expected=LEGACY_CARBOHYDRATE_ATOM_COUNTS,
    name="exact_legacy_carbohydrate_atom_count",
)
