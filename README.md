# BioTSR

**BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology.** It is developed by the **TSR Research Group**
(<https://tsrmethod.github.io/>).

**New to BioTSR? Start with the [user guide](https://github.com/TSRMethod/BioTSR/blob/main/docs/USER_GUIDE.md)**:
installation, all twenty methods with an example each, and the dataset workflow
from download to analysis. `biotsr --help` prints the same link.

> **Version 1.0.0.** Licensed under the Apache License 2.0 (see
> [License](#license)). How to cite it: [Citation](#citation).

## What TSR is

TSR describes a 3D structure by the triangles its points form. Every set of
three points — Cα atoms, or individual heavy atoms — is one triangle. Its
shape (one angle class and one distance class) and the labels of its three
vertices are encoded together as a single integer, the **TSR key**. A structure
is then represented by which keys occur and how often, and structures are
compared through the keys they share.

## What BioTSR is

One Python package, one command line and one shared scientific core — the
**TSR Engine** — replacing a family of independent historical TSR programs.
BioTSR:

* runs every TSR method family through the same shared geometry, binning and
  encoding framework;
* runs each method under one **standard** profile — the accepted science, as
  set by the TSR Research Group — so the same request always means the same
  calculation;
* records, with every result, a fingerprint of what its keys *mean* (the key
  space) and of how its triangles were *drawn* (the sampling profile), and
  refuses to compare results whose keys are not comparable;
* names every calculation's files by one readable prefix in one shared output
  folder, and records each calculation's completion last, so an interrupted
  run never leaves a result that looks complete;
* computes keys with the standard library alone; its one third-party runtime
  dependency, joblib, runs batch rows in parallel when you ask for it.

## Method families

| family | methods | vertices |
|---|---|---|
| Protein Cα | `protein`, `aa_selection` | one Cα per residue |
| Atomic — protein | `amino_acid_atomic`, `amino_acid_atom_selection`, `aa_ind`, `aa_ind_mix`, `aa_ind_specific` | protein heavy atoms |
| Atomic — ligand, carbohydrate, nucleotide | `drug`, `drug_atom_selection`, `carbohydrate`, `carbohydrate_chain`, `carbohydrate_atom_selection`, `nucleotide`, `nucleotide_chain`, `nucleotide_atom_selection` | heavy atoms |
| Cross TSR | `drug_protein_cross`, `carbohydrate_protein_cross`, `nucleotide_protein_cross`, `protein_protein_cross`, `protein_protein_calpha_cross` | two sides, every triangle spans both |

A Cα point is labelled by its residue: the 20 regular amino acids and eleven
modified residues RCSB files write as `HETATM` — `TPO`, `SEP`, `PTR`, `NIY`,
`M3L`, `TRF`, `CSO`, `LLP`, `CME` and `CSD` with codes of their own, and `MSE`
with methionine's — encoded at radix 30 so that no two residue combinations
share a key. Cα results from pre-release builds that used radix 20 and no
modified residues are in another key space; compute them again (the
[changelog](https://github.com/TSRMethod/BioTSR/blob/main/CHANGELOG.md) says
which results this concerns).

An atom in the atomic methods is labelled by its element, read from the PDB
element field (columns 77–78) — carbon 106, nitrogen 108, oxygen 109,
phosphorus 110, sulfur 111, chlorine 93, … at radix 112 — never by its atom
name. Atomic results from pre-release builds that labelled atoms by name are
in another key space and are not compared with current results, so compute
them again.

BioTSR reads PDB files strictly: the full atom-name field and the element
field, with no guessing. A calculation whose selection contains a residue named
`UNK`, an atom without an element symbol, an atom written twice, or a residue
number it cannot read is refused with the records named, rather than computed
without them; records outside the selection never stop a calculation. The
[known limitations](https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md)
list every such rule.

Modifiers compose on one call: amino-acid grouping and SSE (secondary-structure)
annotation on the Cα methods, mirror signing and size filtering almost
everywhere. Grouping puts residues of one group under one code (ASP with GLU,
SER with THR, MSE with MET; each other modified residue on its own). The [capability matrix](https://github.com/TSRMethod/BioTSR/blob/main/docs/CAPABILITY_MATRIX.md) — generated from
the method registry — says exactly which method supports what.

## Installation

BioTSR supports CPython 3.9–3.13 on Linux and macOS. CI runs the full default
test suite on each of 3.9–3.13 on Linux and on macOS (arm64). Windows is not
supported (it is neither tested nor validated). An ordinary
installation brings its one runtime dependency, joblib, which runs batch rows in
parallel; the scientific core itself uses only the standard library.

<!-- not-executed -->
```bash
pip install "biotsr @ git+https://github.com/TSRMethod/BioTSR@v1.0.0"   # the 1.0.0 release, from GitHub
pip install .                                                         # from a clone of https://github.com/TSRMethod/BioTSR
```

The wheel and source distribution are also attached to the
[1.0.0 release](https://github.com/TSRMethod/BioTSR/releases/tag/v1.0.0) on
GitHub. **From PyPI**, `pip install biotsr` installs the same release; it is
available when <https://pypi.org/project/biotsr/> lists version 1.0.0.

Batches run one row at a time unless you ask for more with `--jobs N` or
`--jobs auto`; nothing extra needs installing for that. One optional extra,
`pip install "biotsr[plot]"` (matplotlib), draws the clustered heat-map figure,
whose numbers are all available without it.

## First run

`1UBQ.pdb` is any local PDB-format file; `biotsr fetch 1UBQ --pdb-dir .`
downloads it from the RCSB (see [Getting structures](#getting-structures)).

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

```text
method     protein
profile    standard_protein_calpha_v1
triangles  70300
keys       45349 distinct
output     results/1UBQ_A_protein (4 files, completion record 1UBQ_A_protein.run.json)
```

`biotsr methods` lists all twenty methods; `biotsr methods protein` describes
one: its target arguments, its modifiers and the standard profile it runs. `python -m biotsr` is the same program, for a cluster
where the console script is not on `PATH`.

## Many structures: a manifest

A manifest is a CSV in which **one row is one run**, with the `biotsr run` flags
as columns. An existing dataset table such as `sample_details.csv` works as a
manifest if it has one row per calculation: historical column names are
accepted, other columns can be mapped (`--column label=group`), carried along as
metadata (`--keep-column group1`) or dropped (`--ignore-column notes`), and
`--chain` fills the rows that leave the chain blank, and `--method`
(repeatable) chooses the methods for every row of a manifest without a method
column. A wide dataset sheet
with one column set per method is read with `--sheet tsr` (or `--sheet custom`),
the methods named on the command line (`--method`, repeatable); where a row
could narrow a side of a cross method in several ways, `--cross-side` says
which. The
[CLI guide](https://github.com/TSRMethod/BioTSR/blob/main/docs/CLI_GUIDE.md), §7, has the details, including `--dry-run`,
verified resuming with `--if-exists skip`, and parallel rows with `--jobs`.

<!-- file: samples.csv -->
```csv
sample_id,pdb_id,method,chain,label
ubiquitin,1UBQ,protein,A,ubiquitin
barnase,1BRS,protein,A,barnase
barstar,1BRS,protein,D,barstar
```

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --no-triplets
biotsr analyze common-keys --from-results batch/ --output analysis_out/
```

```text
common keys 9523 over 3 samples
```

## Getting structures

BioTSR reads PDB-format files you supply (mmCIF is not supported), and
downloads them only when asked:

<!-- not-executed -->
```bash
biotsr fetch 1UBQ 1BRS 3PTB --pdb-dir pdb/
biotsr fetch --manifest sample_details.csv --pdb-dir pdb/
biotsr batch --manifest samples.csv --pdb-dir pdb/ --download-missing --output batch/
```

Every download is checked — complete transfer, gzip integrity, a complete PDB
file for the right entry — before it is installed, a valid cached file is
reused, and every id is reported in `fetch_report.csv`. `biotsr fetch` exits 0
only when every id was downloaded, refreshed or already cached.

## Output

Every calculation writes into the `--output` folder directly, its files sharing
one prefix built from the request (`structure_chain[_chain2]_method[_target][_modifiers]`,
for example `1UBQ_A_protein` or `7LF9_D_carbohydrate_NAG_1`):

| file | contents |
|---|---|
| `PREFIX.keys` | every key and its frequency |
| `PREFIX.triplets` | one row per triangle (omit with `--no-triplets`) |
| `PREFIX.meta.json` | the key space, the sampling profile and the selection, with their fingerprints |
| `PREFIX.run.json` | the completion record: the exact command, and every file with its size; written last |

The bin definitions are in the metadata, not the file name. `biotsr describe
results/` lists the calculations in a folder, and `biotsr describe
results/1UBQ_A_protein` shows how one was produced. A batch also writes
`batch_summary.csv`, one row per request. Output folders in the per-sample
directory layout of early pre-release builds are not read; recompute them.

## Python API

```python
import biotsr

result = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", mirror=True)
result.n_triangles              # 70300
result.frequencies              # {key: frequency}
result.key_space.fingerprint    # what the integers mean

biotsr.io.native.save(result, "results/")       # results/1UBQ_A_protein_mirror.keys, ...
loaded = biotsr.load_result("results/1UBQ_A_protein_mirror.keys")
assert loaded.frequencies == result.frequencies
```

Most command-line flags are the keyword of the same name (`--size-filter 25` ⇄
`size_filter=25.0`). The exceptions: `--pdb` is the first positional argument;
`--no-triplets` is `keep_triplets=False`; and `--output` belongs to saving
(`biotsr.io`), not to the calculation. An argument the method does not read is
refused by name, in Python as on the command line. Cross methods use
`biotsr.analyze_cross(...)`; the analysis layer — feature matrices, common
keys, similarity, clustering — is `biotsr.analysis`.

## One standard per method

Each method runs under its own **standard** profile, the accepted TSR science;
`biotsr methods NAME` names it, and every result records it. There is no
option to choose another. Where the standard differs from a historical TSR
program (for example, some historical programs divided angles by the literal
3.14 where BioTSR uses `math.pi`), the method descriptions and the
architecture decisions say so.

Results from different key spaces are never silently combined: comparing a
mirror-signed result with an unsigned one, or a Cα result with an atomic one,
raises an error naming the difference. Results processed with different
settings (a size filter on one and not the other, say) are refused by every
analysis unless you opt in with `--allow-mixed-processing`. A partial result —
one whose selection did not get everything it asked for — exits with status 4
and is refused by every analysis unless you acknowledge it, separately, with
`--allow-partial`. What was selected — structures, chains, residue types,
molecules, residues, atoms — never blocks a comparison.

## Documentation

The [documentation map](https://github.com/TSRMethod/BioTSR/blob/main/docs/README.md) lists everything. The main pages:

| | |
|---|---|
| [Tutorials](https://github.com/TSRMethod/BioTSR/blob/main/docs/tutorials/README.md) | runnable walk-throughs of every established workflow |
| [User guide](https://github.com/TSRMethod/BioTSR/blob/main/docs/USER_GUIDE.md) | all twenty methods, each with an example, and the dataset workflow from download to analysis |
| [Command-line guide](https://github.com/TSRMethod/BioTSR/blob/main/docs/CLI_GUIDE.md) | every command, the manifest format, exit codes, HPC use |
| [Python API guide](https://github.com/TSRMethod/BioTSR/blob/main/docs/API_GUIDE.md) | run, inspect, save, load, analyse |
| [Capability matrix](https://github.com/TSRMethod/BioTSR/blob/main/docs/CAPABILITY_MATRIX.md) | what each method supports, generated from the registry |
| [Reproducibility note](https://github.com/TSRMethod/BioTSR/blob/main/docs/REPRODUCIBILITY_NOTE.md) | which command computes which published method |
| [Known limitations](https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md) | what 1.0 does not do, stated precisely |
| [Methods](https://github.com/TSRMethod/BioTSR/blob/main/docs/methods/README.md) | cross TSR, SSE-TSR and analysis: the definitions the standard profiles implement |
| [Architecture](https://github.com/TSRMethod/BioTSR/blob/main/docs/ARCHITECTURE.md), [decisions](https://github.com/TSRMethod/BioTSR/blob/main/docs/architecture/README.md) | how the package is built, and why |
| [Releasing](https://github.com/TSRMethod/BioTSR/blob/main/docs/RELEASING.md) | versioning, Python support, CI, what ships in the package |
| [Licensing and provenance](https://github.com/TSRMethod/BioTSR/blob/main/docs/LICENSING.md) | where each part of BioTSR comes from, and under which terms |

The reference outputs and validation tools that check the standard results
are part of the source repository (under `validation/`), not of the installed
package.

## Citation

If you use BioTSR, please cite it using the metadata in
[`CITATION.cff`](https://github.com/TSRMethod/BioTSR/blob/main/CITATION.cff) (GitHub's "Cite this repository" button reads
it):

> Khajouie, P., Rauniyar, K., Raghavan, V. V., & Xu, W. *BioTSR* (version
> 1.0.0) [Computer software]. https://github.com/TSRMethod/BioTSR

## License

BioTSR is licensed under the **Apache License 2.0** — see [`LICENSE`](https://github.com/TSRMethod/BioTSR/blob/main/LICENSE).
Copyright 2026 TSR Research Group.

BioTSR includes or derives from third-party material that keeps its own
licence: tables from `TSR-Package` and the element codes of the `CarbohydratesTSR`
atom lexicon (MIT),
and a clustering implementation that follows SciPy's (BSD 3-Clause). Their
notices, and attribution for the historical TSR repositories BioTSR draws on,
are in [`NOTICE`](https://github.com/TSRMethod/BioTSR/blob/main/NOTICE).

The PDB structures under `tests/data/structures/` are unmodified entries from
the RCSB Protein Data Bank. `tests/data/synthetic/SYNPHOS.pdb` is not a real
entry: it is a synthetic structure written for the tests (its README says how).
