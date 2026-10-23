# Synthetic manifests

`dataset_sheet.csv` is a small, synthetic example of the historical TSR
dataset-sheet format: one PDB entry per row, one column per calculation
(`biotsr.io.input_columns`). It was written for the test suite. Every value in
it (chains, ligands, residue and atom selections, labels and comments) was chosen
from the public PDB entries in `../structures/`, and none was copied from any
dataset. It keeps the format's real quirks on purpose: the `car_atom_selected`
header spelling, `comment_*` columns, and both residue spellings (`LYS_207`,
`SER266`) in one cell.

`wide_sheet.csv` is a second synthetic sheet in the same format, written for the
wide-sheet workflow (`biotsr batch --sheet tsr --method ...`). Its rows share a
structure and a protein chain across several ligand or glycan rows, name some
molecules without an identifier (every copy on the chain), fill the second
protein chain's columns (`aa_chain_2`, `aa_ind_2`, ...) in one row, and mix both
residue spellings. Every value was chosen from the public structures in
`../structures/`; none was copied from any dataset.
