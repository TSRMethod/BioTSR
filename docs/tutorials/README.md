# BioTSR tutorials

Short, runnable walk-throughs of the established BioTSR workflows. For an
example of every method and of the dataset workflow, see the
[user guide](../USER_GUIDE.md); for complete references, the
[command-line guide](../CLI_GUIDE.md) and the
[Python API guide](../API_GUIDE.md).

## Before you start

1. **Install BioTSR** from a clone of the repository:
   `python -m venv .venv`, `source .venv/bin/activate`, `pip install .`
   Check it with `biotsr --version`.
2. **Get the example structures.** The tutorials use unmodified RCSB PDB
   entries tracked in the repository under `tests/data/structures/`:
   `1UBQ.pdb` (ubiquitin), `1BRS.pdb` (barnase–barstar), `3PTB.pdb`
   (trypsin with benzamidine), `1BNA.pdb` (a B-DNA dodecamer), `7LF9.pdb` (a
   glycoprotein; tutorials 02, 04 and 05) and `6G9M.pdb` and `2AC0.pdb` (a
   kinase with a ligand and a protein–DNA complex; the wide sheet of tutorial
   05). Tutorial 05 also uses `SYNPHOS.pdb`, a **synthetic** structure under
   `tests/data/synthetic/` (not a real entry; its README says how it was
   written), whose chain A holds the modified amino acids TPO, SEP, PTR and
   CME as `HETATM` records and the insertion-code residues 10A and 10B. Copy
   them into an empty working directory, and also into a `pdb/`
   subdirectory, which the manifest examples (05, 07, 08) read from:

   ```
   mkdir -p tutorial/pdb
   cp tests/data/structures/{1UBQ,1BRS,3PTB,1BNA,7LF9,6G9M,2AC0}.pdb tutorial/
   cp tests/data/synthetic/SYNPHOS.pdb tutorial/
   cp tutorial/*.pdb tutorial/pdb/
   cd tutorial
   ```

   Without a clone, `biotsr fetch 1UBQ 1BRS 3PTB 1BNA 7LF9 6G9M 2AC0 --pdb-dir .`
   downloads the same entries from the RCSB PDB (the tutorials themselves
   never use the network); `SYNPHOS.pdb` exists only in the repository.
3. **Run the commands from that directory.** Each page starts from a fresh
   directory, so run the pages in separate directories or delete the outputs
   between them: a calculation whose files are already in `--output` is
   refused rather than overwritten (`--if-exists` changes that).

To use your own structures, replace the file name, chain and target
(`--entity-name`, `--residue`, …) with yours. `biotsr methods NAME` lists the
options each method needs.

## The tutorials

| tutorial | what it covers |
|---|---|
| [01 — Protein Cα TSR](01_protein_calpha.md) | one structure, one chain, the three Cα modifiers, the output files |
| [02 — Atomic TSR](02_atomic.md) | one residue, one ligand, one copy or every copy of a molecule, identifiers |
| [03 — SSE-TSR](03_sse.md) | secondary-structure annotation, and why it changes no key |
| [04 — Cross TSR](04_cross.md) | a drug–protein interface; one glycan copy or all copies against a protein chain |
| [05 — A dataset manifest](05_batch_manifest.md) | many structures in one command; a wide dataset sheet with methods chosen by `--method`; participants and cross sides; modified amino acids, duplicate rows and residue ranges |
| [06 — Load a saved result](06_load_result.md) | reading native output back into Python |
| [07 — Common keys and similarity](07_analysis.md) | the analysis workflow, and what has to match for samples to be compared |
| [08 — Similarity, clustering and the clustered heat map](08_similarity_clustering.md) | saved results to Jaccard matrices, the tree, the ordered heat map and its provenance; an optional flat cut by a stated rule (Python; the figure needs `biotsr[plot]`) |

**Every command and every code block on these pages is executed by the test
suite** (`tests/integration/test_tutorials.py`), and every output block is
checked against what the command actually printed. If a page stops being true,
the tests fail.

Conventions used by the harness, so you can read the pages as plain Markdown:

* a `bash` block is run line by line from a fresh working directory that holds
  the tutorial's PDB files;
* a `text` block directly after a `bash` block is output that must appear;
* a `python` block is executed, in order, in one namespace per page;
* a block announced with `<!-- file: NAME -->` is written to `NAME` first;
* a block announced with `<!-- requires: plot -->` draws a figure and runs only
  where matplotlib is installed (`pip install "biotsr[plot]"`).
