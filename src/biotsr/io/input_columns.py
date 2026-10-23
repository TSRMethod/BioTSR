"""The historical dataset-sheet vocabulary: column -> method, and header aliases.

The TSR group's historical dataset sheets describe one PDB entry per row, with a
column per calculation, and come with a dictionary defining each column.  This module records what each of those columns means
and which method serves it, so a caller working from such a sheet never has to
guess and never silently renames data.  It is not a reader itself: the command line reads such a sheet
method-aware with ``biotsr batch --sheet tsr`` (``biotsr.cli.sheet``, whose
built-in layout follows this table), and reads one-run-per-row manifests
(``biotsr.cli.manifest``), whose aliases accept these column names.

**This is deliberately not a CSV driver.**  It contains no file reading, no row
iteration and no defaults for blank cells: an undocumented or blank field has no
invented meaning, and a column with no implementation is listed as such rather
than quietly dropped.  What it provides is the mapping a driver would need, in
one auditable place.

**The known header mismatch is preserved as an explicit alias, not a rename.**
The column dictionary names the column ``carb_atom_selected``; the sheets'
actual header is ``car_atom_selected``.  Both are accepted, the canonical name is
recorded, and the data is never rewritten.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

__all__ = [
    "ColumnSpec",
    "DATASET_COLUMNS",
    "HEADER_ALIASES",
    "canonical_column",
    "method_for_column",
]


@dataclass(frozen=True)
class ColumnSpec:
    """What one dataset-sheet column means, and which method (if any) serves it."""

    column: str
    meaning: str
    #: The registered method name, or ``None`` when nothing implements it yet.
    method: Optional[str]
    #: Which column supplies the chain.
    chain_column: Optional[str] = None
    #: Other columns this one needs to be interpretable.
    companion_columns: Tuple[str, ...] = ()
    #: The ``analyze()`` keyword this column's value is passed as.
    parameter: Optional[str] = None
    #: Set for columns that exist but have no implementation, with the reason.
    unimplemented_reason: str = ""


#: Each column's meaning, in this package's words, with the method that serves
#: it.  The format is the historical one; the sheets themselves are not
#: distributed (``tests/data/manifests/dataset_sheet.csv`` is a synthetic example).
DATASET_COLUMNS: Mapping[str, ColumnSpec] = {
    spec.column: spec
    for spec in (
        ColumnSpec(
            column="protein",
            meaning=(
                "PDB ID, for every type of key and triplet calculation. No length "
                "is assumed anywhere in this package: identifiers longer than "
                "four characters are accepted unchanged."
            ),
            method=None,
            parameter="structure_id",
        ),
        ColumnSpec(
            column="aa_chain",
            meaning="the protein chain for every aa_* calculation",
            method=None,
            parameter="chain",
        ),
        ColumnSpec(
            column="drug_chain", meaning="the chain for drug calculations",
            method=None, parameter="chain",
        ),
        ColumnSpec(
            column="carb_chain",
            meaning="the chain for carbohydrate monomer-residue and polymer-chain "
                    "calculations",
            method=None, parameter="chain",
        ),
        ColumnSpec(
            column="nt_chain",
            meaning="the chain for nucleotide (residue), DNA (chain) and RNA "
                    "(chain) calculations",
            method=None, parameter="chain",
        ),
        ColumnSpec(
            column="aa_chain_2",
            meaning="the interacting protein's chain, for the *_2 columns",
            method=None, parameter="chain",
        ),
        ColumnSpec(
            column="aa_group",
            meaning=(
                "annotation for the proteins in aa_chain (a family or target name). "
                "NOT the aa_grouping modifier -- ADR-023."
            ),
            method=None,
            unimplemented_reason=(
                "a dataset annotation label; it feeds no capability and must never "
                "be bound to the aa_grouping modifier"
            ),
        ),
        ColumnSpec(
            column="aa_ind",
            meaning="one specified amino-acid type in aa_chain, all non-H/D atoms",
            method="aa_ind", chain_column="aa_chain", parameter="residue_types",
        ),
        ColumnSpec(
            column="aa_ind_mix",
            meaning="two or more specified amino-acid types in aa_chain",
            method="aa_ind_mix", chain_column="aa_chain", parameter="residue_types",
        ),
        ColumnSpec(
            column="aa_ind_specific",
            meaning="specified amino-acid types with specific residue numbers",
            method="aa_ind_specific", chain_column="aa_chain", parameter="residues",
        ),
        ColumnSpec(
            column="aa_ind_2", meaning="aa_ind against aa_chain_2",
            method="aa_ind", chain_column="aa_chain_2", parameter="residue_types",
        ),
        ColumnSpec(
            column="aa_ind_mix_2", meaning="aa_ind_mix against aa_chain_2",
            method="aa_ind_mix", chain_column="aa_chain_2", parameter="residue_types",
        ),
        ColumnSpec(
            column="aa_ind_specific_2", meaning="aa_ind_specific against aa_chain_2",
            method="aa_ind_specific", chain_column="aa_chain_2", parameter="residues",
        ),
        ColumnSpec(
            column="drug_name",
            meaning="the drug in drug_chain, with drug_id; all non-H/D atoms",
            method="drug", chain_column="drug_chain",
            companion_columns=("drug_id",), parameter="entity_name",
        ),
        ColumnSpec(
            column="drug_id", meaning="the drug's residue sequence number",
            method="drug", chain_column="drug_chain", parameter="entity_id",
        ),
        ColumnSpec(
            column="carb_name",
            meaning="the carbohydrate in carb_chain, with carb_id",
            method="carbohydrate", chain_column="carb_chain",
            companion_columns=("carb_id",), parameter="entity_name",
        ),
        ColumnSpec(
            column="carb_id", meaning="the carbohydrate's residue sequence number",
            method="carbohydrate", chain_column="carb_chain", parameter="entity_id",
        ),
        ColumnSpec(
            column="nt_name", meaning="the nucleotide in nt_chain, with nt_id",
            method="nucleotide", chain_column="nt_chain",
            companion_columns=("nt_id",), parameter="entity_name",
        ),
        ColumnSpec(
            column="nt_id", meaning="the nucleotide's residue sequence number",
            method="nucleotide", chain_column="nt_chain", parameter="entity_id",
        ),
        ColumnSpec(
            column="aa_selection",
            meaning="the listed amino acids in aa_chain, using only Ca atoms",
            method="aa_selection", chain_column="aa_chain", parameter="residues",
        ),
        ColumnSpec(
            column="aa_selection_2", meaning="aa_selection against aa_chain_2",
            method="aa_selection", chain_column="aa_chain_2", parameter="residues",
        ),
        ColumnSpec(
            column="aa_atom_selection",
            meaning="the specific amino acids and atoms in aa_chain",
            method="amino_acid_atom_selection", chain_column="aa_chain",
            parameter="atoms",
        ),
        ColumnSpec(
            column="aa_atom_selection_2",
            meaning="aa_atom_selection against aa_chain_2",
            method="amino_acid_atom_selection", chain_column="aa_chain_2",
            parameter="atoms",
        ),
        ColumnSpec(
            column="drug_atom_selected",
            meaning="the specific drug atoms in drug_chain/drug_name/drug_id",
            method="drug_atom_selection", chain_column="drug_chain",
            companion_columns=("drug_name", "drug_id"), parameter="atoms",
        ),
        ColumnSpec(
            column="carb_atom_selected",
            meaning="the specific carbohydrate atoms in carb_chain/carb_name/carb_id",
            method="carbohydrate_atom_selection", chain_column="carb_chain",
            companion_columns=("carb_name", "carb_id"), parameter="atoms",
        ),
        ColumnSpec(
            column="nt_atom_selected",
            meaning="the specific nucleotide atoms in nt_chain/nt_name/nt_id",
            method="nucleotide_atom_selection", chain_column="nt_chain",
            companion_columns=("nt_name", "nt_id"), parameter="atoms",
        ),
    )
}

#: Header spellings that are **not** the documented name.  Preserved as explicit
#: aliases rather than renamed in the data.
#:
#: ``car_atom_selected`` is the header historical dataset sheets actually carry;
#: their column dictionary names the same field ``carb_atom_selected``.  Both
#: resolve; the data is never edited.
HEADER_ALIASES: Mapping[str, str] = {
    "car_atom_selected": "carb_atom_selected",
    # Historical per-script variants, from WebsiteTSR's own CLI overrides.
    "atom_selected": "drug_atom_selected",
    "aa": "aa_selection",
    "chain": "aa_chain",
}


def canonical_column(header: str) -> str:
    """Resolve a header spelling to its documented name.

    Raises rather than guessing: a header this package has never seen needs a
    documented meaning, not a pattern match.
    """
    name = str(header).strip()
    if name in DATASET_COLUMNS:
        return name
    if name in HEADER_ALIASES:
        return HEADER_ALIASES[name]
    raise KeyError(
        "unknown dataset column %r; documented columns are %s (aliases: %s)"
        % (header, ", ".join(sorted(DATASET_COLUMNS)),
           ", ".join(sorted(HEADER_ALIASES)))
    )


def method_for_column(header: str) -> Optional[str]:
    """The registered method a column drives, or ``None`` if none does.

    ``None`` means "this column has no implementation", which is a statement, not
    a silent omission: :data:`DATASET_COLUMNS` records why.
    """
    return DATASET_COLUMNS[canonical_column(header)].method
