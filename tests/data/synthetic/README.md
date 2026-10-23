# Synthetic structures

`SYNPHOS.pdb` is **not a real structure**. `make_synthetic.py` in this folder
writes it, deterministically (no randomness, no network); running the script
again reproduces the file byte for byte. Its geometry is arbitrary: residues
are laid out along two parallel paths and each atom sits at a fixed offset from
its residue's centre. What the file is for is its records, which no small
public entry combines:

| chain | residue | record | what it exercises |
|---|---|---|---|
| A | TPO 3, SEP 5, PTR 7 | `HETATM` | the phosphorylated amino acids every protein method reads (ADR-101, ADR-107) |
| A | CME 8 | `HETATM` | a modified residue the Cα methods read (ADR-107) and the atomic methods do not: named, never dropped |
| A | GLY 10A, ALA 10B | `ATOM` | insertion codes in identifiers and ranges (ADR-102) |
| A | GLY 1 ... VAL 12 | `ATOM` | ordinary residues around them |
| B | TPO 3 among five others | `ATOM`/`HETATM` | a second chain for the cross methods |
| A | LIG 101, CA 102 (after `TER`) | `HETATM` | a ligand, and a calcium ion whose atom is also named `CA` |

The tests that read it are `tests/integration/test_modified_residue_identifier_row_rules.py`,
`tests/integration/test_processing_identity.py` and
`tests/integration/test_modified_residue_calpha.py`; the tutorial
`docs/tutorials/05_batch_manifest.md` uses it too.
