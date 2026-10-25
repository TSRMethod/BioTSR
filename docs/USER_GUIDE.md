# BioTSR user guide

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology. This guide shows how to use it, method by
method: which command or Python call answers which question, what goes in and
what comes out. The [command-line guide](CLI_GUIDE.md) and the
[Python API guide](API_GUIDE.md) are the complete references; this page links to
them wherever it summarises.

**Every command and code block on this page is executed by the test suite**
(`tests/integration/test_tutorials.py`), from a directory that holds the
example structures (§2.1), except the few marked as not executed because they
need the network, a cluster or a clone. Output shown after a command is checked
against what the command printed, with runs of spaces collapsed.

## Contents

- [The methods at a glance](#the-methods-at-a-glance)
- [Scientific methods and workflow features](#scientific-methods-and-workflow-features)
1. [Installation](#1-installation)
2. [Structures: your files, or downloaded](#2-structures-your-files-or-downloaded)
3. [Protein Cα TSR](#3-protein-cα-tsr)
4. [Selected residues and amino-acid grouping](#4-selected-residues-and-amino-acid-grouping)
5. [Amino-acid atomic TSR](#5-amino-acid-atomic-tsr)
6. [SSE-TSR](#6-sse-tsr)
7. [Drug (ligand) TSR](#7-drug-ligand-tsr)
8. [Carbohydrate TSR](#8-carbohydrate-tsr)
9. [Nucleotide TSR](#9-nucleotide-tsr)
10. [Cross TSR](#10-cross-tsr)
11. [Modifiers](#11-modifiers)
12. [Arguments, entities and selections](#12-arguments-entities-and-selections)
13. [Datasets: manifests, sheets and batches](#13-datasets-manifests-sheets-and-batches)
14. [Memory and disk](#14-memory-and-disk)
15. [Output files](#15-output-files)
16. [Loading previous results](#16-loading-previous-results)
17. [Analysis](#17-analysis)
18. [Common workflows](#18-common-workflows)
19. [Troubleshooting](#19-troubleshooting)
20. [Where to go next](#20-where-to-go-next)

---

## The methods at a glance

BioTSR 1.0 has **twenty methods** in three families. `biotsr methods`
lists them and `biotsr methods NAME` describes one: the target arguments it
takes, its modifiers and the standard profile it runs. The [capability matrix](CAPABILITY_MATRIX.md),
generated from the method registry, is the authoritative table of what each
method supports.

| method | computes TSR over | target arguments | example |
|---|---|---|---|
| `protein` | every Cα of one chain | — | [§3](#3-protein-cα-tsr) |
| `aa_selection` | the Cα of listed residues | `--residues '1-10;20'` | [§4](#4-selected-residues-and-amino-acid-grouping) |
| `amino_acid_atomic` | the heavy atoms of one residue, or of each residue in turn | `--residue` with `--residue-name`; `--residue-name` alone for each residue of that name; neither for every residue | [§5](#5-amino-acid-atomic-tsr) |
| `amino_acid_atom_selection` | named atoms across residues | `--atoms 'LYS_11_NZ;LEU_15_CB'` | [§5](#5-amino-acid-atomic-tsr) |
| `aa_ind` | every heavy atom of every residue of one type | `--residue-types LEU` | [§5](#5-amino-acid-atomic-tsr) |
| `aa_ind_mix` | every heavy atom of residues of two or more types | `--residue-types 'LEU;ILE'` | [§5](#5-amino-acid-atomic-tsr) |
| `aa_ind_specific` | every heavy atom of named residues | `--residues 'LEU_15;LEU50'` | [§5](#5-amino-acid-atomic-tsr) |
| `drug` | one ligand copy, or every copy of that ligand on the chain | `--entity-name`, optional `--entity-id` (omitted: every copy) | [§7](#7-drug-ligand-tsr) |
| `drug_atom_selection` | named atoms of one ligand copy, or of every copy | `--entity-name`, optional `--entity-id`, `--atoms` | [§7](#7-drug-ligand-tsr) |
| `carbohydrate` | one carbohydrate monomer, or every copy of that sugar on the chain | `--entity-name`, optional `--entity-id` | [§8](#8-carbohydrate-tsr) |
| `carbohydrate_chain` | a whole carbohydrate chain, pooled | optional `--entity-names` | [§8](#8-carbohydrate-tsr) |
| `carbohydrate_atom_selection` | named atoms of one monomer, or of every copy | `--entity-name`, optional `--entity-id`, `--atoms` | [§8](#8-carbohydrate-tsr) |
| `nucleotide` | one nucleotide residue, or every copy of that nucleotide on the chain | `--entity-name`, optional `--entity-id` | [§9](#9-nucleotide-tsr) |
| `nucleotide_chain` | a whole DNA or RNA strand, pooled | optional `--entity-names` | [§9](#9-nucleotide-tsr) |
| `nucleotide_atom_selection` | named atoms of one nucleotide, or of every copy | `--entity-name`, optional `--entity-id`, `--atoms` | [§9](#9-nucleotide-tsr) |
| `drug_protein_cross` | ligand × protein interface | `--chain-2`, `--entity-name`, optional `--entity-id` | [§10](#10-cross-tsr) |
| `carbohydrate_protein_cross` | carbohydrate × protein interface | `--chain-2`, `--entity-name`, optional `--entity-id` | [§10](#10-cross-tsr) |
| `nucleotide_protein_cross` | nucleotide × protein interface | `--chain-2`, `--entity-name`, optional `--entity-id` | [§10](#10-cross-tsr) |
| `protein_protein_cross` | protein × protein interface, all heavy atoms | `--chain-2` | [§10](#10-cross-tsr) |
| `protein_protein_calpha_cross` | protein × protein interface, Cα of interface residues | `--chain-2` | [§10](#10-cross-tsr) |

Every method also needs `--pdb`, `--chain` and `--output`. An argument a
method does not read is refused by name (§12.2).

## Scientific methods and workflow features

The twenty methods are the TSR calculations the research group's historical
programs and dataset definitions describe, implemented on one engine. Each runs
under one **standard profile**: the accepted science, as set by the TSR
Research Group. There is no option to choose another (§12.1). Where the
standard differs from a historical program — some divided angles by the literal
3.14 where BioTSR uses `math.pi`, for example — the
[methods records](methods/README.md) say so.

The rest of BioTSR 1.0 is **workflow**, and changes no key: downloading and
caching structures (§2.4), configurable manifest columns and explicit chain
defaults (§13.2–13.3), dry runs (§13.4), one shared output folder with
completion records and duplicate-output checks (§13.5, §15), parallel batches
(§13.6), verified resume and interruption handling (§13.8), wide dataset sheets
read method by method (§13.9), and streamed triplet output (§14).

BioTSR reads PDB-format files only, on Linux and macOS; Windows is not
supported.

---

## 1. Installation

BioTSR needs Python 3.9 or later and supports CPython 3.9–3.13 on Linux and
macOS. CI runs the full default test suite on each of 3.9–3.13 on Linux and on
macOS (arm64). Windows is not supported.

<!-- not-executed -->
```bash
git clone https://github.com/TSRMethod/BioTSR
cd BioTSR
python -m venv .venv
source .venv/bin/activate
pip install .
```

The 1.0.0 release also installs without a clone, from GitHub —
`pip install "biotsr @ git+https://github.com/TSRMethod/BioTSR@v1.0.0"`, or the
wheel attached to the
[release](https://github.com/TSRMethod/BioTSR/releases/tag/v1.0.0) — and from
PyPI, `pip install biotsr`, once <https://pypi.org/project/biotsr/> lists
version 1.0.0. Every way, the installation brings BioTSR's
one runtime dependency, **joblib**, which runs batch rows in parallel
(`--jobs`, §13.6); the calculation itself uses only the standard library.
Nothing else is needed to run every method in parallel. The optional
clustered heat-map figure needs matplotlib: `pip install ".[plot]"` from a
checkout, or `pip install "biotsr[plot]"` (§17).

Check the installation:

```bash
biotsr --version        # prints: biotsr 1.0.0
biotsr methods
```

```text
METHOD FAMILY MODIFIERS
protein calpha aa_grouping, mirror, size_filter, sse
drug_protein_cross cross size_filter
```

If the `biotsr` command is not found (common on a compute node),
`python -m biotsr` is the same program.

---

## 2. Structures: your files, or downloaded

BioTSR reads **PDB-format** files (`.pdb` or `.ent`). mmCIF is not supported;
convert it to PDB first. How a file is read — which model, which alternate
locations, where a chain ends — is in §12.5.

### 2.1 Example structures

The repository tracks unmodified RCSB PDB entries under
`tests/data/structures/`. This page uses 1UBQ (ubiquitin), 1BRS
(barnase–barstar), 3PTB (trypsin with benzamidine), 1BNA (B-DNA), 5WT9 and
7LF9 (glycoproteins) and 2AC0 (a protein–DNA complex). The examples of
modified amino acids and insertion codes use `SYNPHOS.pdb`, a synthetic
structure under `tests/data/synthetic/` (not a real entry; its README says how
it was written): chain A holds the `HETATM` residues TPO 3, SEP 5, PTR 7 and
CME 8, and the residues 10A and 10B between 10 and 11. To follow along, copy
them into a working directory, and into a `pdb/` subdirectory for the manifest
examples:

<!-- not-executed -->
```bash
mkdir -p tsr-work/pdb
cp tests/data/structures/*.pdb tsr-work/
cp tests/data/synthetic/SYNPHOS.pdb tsr-work/
cp tsr-work/*.pdb tsr-work/pdb/
cd tsr-work
```

To use your own structures, replace the file name, chain and target with yours.

### 2.2 One structure

One structure is one `biotsr run` command, or one `biotsr.analyze(...)` call:

<!-- not-executed -->
```bash
biotsr run --method protein --pdb path/to/structure.pdb --chain A --output results/
```

Every method needs a chain; BioTSR never guesses one.

### 2.3 Several structures

`biotsr run` takes **one** structure file. For several there are two routes:

* **the command line**: describe the runs in a manifest and run `biotsr batch`
  (§13);
* **Python**: loop over your structures and call `biotsr.analyze` for each.

```python
import biotsr

structures = {
    "ubiquitin": ("1UBQ.pdb", "A"),
    "barnase": ("1BRS.pdb", "A"),
    "barstar": ("1BRS.pdb", "D"),
}
results = {}
for name, (path, chain) in structures.items():
    results[name] = biotsr.analyze(path, method="protein", chain=chain,
                                   keep_triplets=False)
    biotsr.io.native.save(results[name], "loop_results/", prefix=name,
                          triplets=False)     # loop_results/ubiquitin.keys, ...

assert results["ubiquitin"].n_triangles == 70300
assert results["barstar"].n_triangles == 105995
```

### 2.4 Downloading structures

Nothing uses the network unless you ask. `biotsr fetch` downloads entries from
the RCSB PDB into a directory, which is also the cache it checks first:

<!-- not-executed -->
```bash
biotsr fetch 1UBQ 4HHB --pdb-dir pdb/
biotsr fetch --id-file ids.txt --pdb-dir pdb/
```

`--manifest` reads the ids from **one column** of any CSV (`pdb_id`, an alias
such as `protein`, or the column named by `--column pdb_id=HEADER`) and checks
no other column, so an existing dataset sheet works as it is:

<!-- file: ids.csv -->
```csv
protein,notes
1UBQ,ubiquitin
1BRS,barnase-barstar
1ubq,listed twice
```

```bash
biotsr fetch --manifest ids.csv --pdb-dir pdb/ --dry-run
biotsr fetch --manifest ids.csv --pdb-dir pdb/
```

```text
ids 3 requested, 2 distinct
result 2 cached
dry run nothing was downloaded or written
report pdb/fetch_report.csv
```

Here both entries are already in `pdb/`, so nothing is downloaded. In general:

* each download is verified (complete transfer, gzip integrity, a complete PDB
  file for the right entry) before it is installed under its final name;
* a valid file already in `--pdb-dir` is reused and never overwritten;
  `--refresh` downloads again and replaces it only after the new file passes
  every check;
* `--download-workers N` (1 to 16, default 4), `--timeout SECONDS` (default
  30) and `--retries N` (default 3) control the downloads, independently of
  `--jobs`. An invalid value is refused, never replaced by a default;
* `--timeout` limits each wait — for a connection, or for the next piece of a
  response — **not** a whole download. A server that keeps sending, however
  slowly, is never timed out, so there is no upper bound on how long one
  download can take. A failed attempt is retried after an exponential backoff
  of up to 60 s;
* Ctrl-C (or SIGTERM) stops downloads promptly: queued ids never start, no
  retry follows, a backoff wait ends, a connection being made is abandoned and
  a response being received is stopped at once. Only a host-name lookup and a
  TLS handshake step in progress are waited for. Every id is reported
  (`interrupted` or `not_started` if it did not finish), and no partial file is
  ever installed;
* every id gets a row in `fetch_report.csv`. The command exits `0` only when
  every id was downloaded, refreshed or already cached, and `3` otherwise —
  including `refresh_failed`, where `--refresh` failed and the old file was
  kept;
* an obsolete entry (one with an `OBSLTE` record) is downloaded, with a note in
  the report; an entry RCSB serves only as mmCIF is reported as
  `pdb_format_unavailable`;
* Ctrl-C or SIGTERM stops the command (exit `130`) with the report written; a
  download already in flight ends at its next network read, and no new attempt
  starts.

`biotsr batch ... --download-missing` fetches whatever a manifest needs before
computing (§13). From Python: `biotsr.fetch_structures(["1UBQ"], "pdb/")`. The
[CLI guide](CLI_GUIDE.md), §7a, lists every report status.

---

## 3. Protein Cα TSR

One point per residue (its Cα atom), labelled by amino acid; all C(n,3)
triangles of one chain.

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

```text
method protein
profile standard_protein_calpha_v1
triangles 70300
keys 45349 distinct
output results/1UBQ_A_protein (4 files, completion record 1UBQ_A_protein.run.json)
```

The four files share the prefix `1UBQ_A_protein`: the keys, the metadata
sidecar, the triplets and the completion record (§15).

The same from Python:

```python
import biotsr

result = biotsr.analyze("1UBQ.pdb", method="protein", chain="A")
assert result.n_triangles == 70300          # C(76, 3)
assert result.n_unique_keys == 45349
key, frequency = next(iter(result.frequencies.items()))   # key -> frequency
assert result.key_space.label_scheme.name == "amino_acid_30"
```

**Residue labels.** A Cα point is labelled by its residue: the 20 regular amino
acids have codes 4–23, and eleven modified residues have codes of their own —
phosphothreonine `TPO` 24, phosphoserine `SEP` 25, phosphotyrosine `PTR` 26,
3-nitrotyrosine `NIY` 27, trimethyllysine `M3L` 28, formyl-tryptophan `TRF` 29,
S-hydroxycysteine `CSO` 30, lysine-pyridoxal phosphate `LLP` 31,
S,S-(2-hydroxyethyl)thiocysteine `CME` 32 and 3-sulfinoalanine `CSD` 33 —
while selenomethionine `MSE` shares methionine's 13. The 30 codes are encoded at
radix 30, so no two residue combinations can share a key (ADR-107). RCSB files
write these residues as `HETATM` records; a Cα method reads one when its CA
atom's element field is carbon or blank — no other atom of the residue is
needed, as for an `ATOM` record — and reads no other `HETATM` record: a
ligand, or a calcium ion whose atom is also named `CA`, is never a point. Any
other modified residue is named as `not represented:` (§5), and so is a
supported one present without a readable Cα, with the reason (no CA atom; a CA
only as alternate location B or later; a `HETATM` CA declaring another
element). A residue `aa_selection` names that the chain contains is never
reported as "not found".

```python
synthetic = biotsr.analyze("SYNPHOS.pdb", method="protein", chain="A", keep_triplets=False)
assert synthetic.selection.residues_used[:9] == ("1", "2", "3", "4", "5", "6", "7", "8", "9")
assert synthetic.warnings == ()          # TPO 3, SEP 5, PTR 7 and CME 8 are points
```

**Cα results from pre-release builds at radix 20 must be computed again.**
Those builds encoded the 20 regular residues at radix 20 and did not read
modified residues, so every Cα key had another value. Those results load, but
their key space differs: an analysis refuses to combine them with current ones, and
`--if-exists skip` refuses to reuse them, leaving their files untouched;
compute them again with `--if-exists overwrite` or into a new `--output`
(§13.8).

---

## 4. Selected residues and amino-acid grouping

**`aa_selection`** runs Cα TSR over a listed set of residues. Ranges use `-`
and items are separated by `;`. Quote the list, because `;` is special to the
shell:

```bash
biotsr run --method aa_selection --pdb 1UBQ.pdb --chain A --residues '1-10;20;30-35' --output results/
```

```text
method aa_selection
triangles 680
```

A residue is written as its PDB number with an optional one-letter insertion
code (`48`, `52A`, `-5`; leading zeros are dropped). A range is two plain
numbers, ascending, both ends included (`10-12`, `-5--3`). A range with an
insertion code at either end (`27A-29`), a descending range (`12-10`) and any
other malformed item are refused, naming the item. So is a range that spans a
residue with an insertion code, because whether `9-11` should include `10A` has
no single obvious answer: list such a residue to include it, or list the
residues one by one to leave it out. On SYNPHOS chain A (§2.1):

```bash
biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --residues '9-11;10A;10B' --output ranges/
biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --residues '9;10;11' --output ranges/
```

```text
triangles 10
triangles 1
```

```python
listed = biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A",
                        residues="10B;11;10A;9-10", keep_triplets=False)
assert listed.selection.residues_used == ("9", "10", "10A", "10B", "11")  # file order
try:
    biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A", residues="9-11")
except biotsr.InvalidRequest as error:
    assert "span residue(s) 10A, 10B" in str(error)
else:
    raise AssertionError("a range over insertion-code residues must be refused")
```

The order of a list never matters: the points follow the structure file. The
command line, a manifest and Python share this grammar.

**Amino-acid grouping** is not a separate method. It is the `--aa-grouping`
modifier on either Cα method. It groups the 20 regular amino acids into 13
groups (for example ASP with GLU, SER with THR), puts MSE with MET, and keeps
each other modified residue in a group of its own: 23 grouped codes, also at
radix 30 ([§11](#11-modifiers)). Residues of one group share a code, so their
triangles share keys on purpose; grouped and ungrouped keys are different key
spaces and are never compared.

---

## 5. Amino-acid atomic TSR

Atomic methods use individual heavy atoms as vertices, each labelled by its
element — read from the element field, PDB columns 77–78, never from the atom
name — with 18 distance classes. The element labels ship with the package
(`element_labels_v1.csv`: carbon 106, nitrogen 108, oxygen 109, phosphorus
110, sulfur 111, chlorine 93, calcium 99, …, radix 112); you never supply them.
An atom whose element field is blank, or names no element in the table
(`X`), stops a calculation that selects it (§12.5).

The atom name is the atom's identity — in `--atoms`, in the triplet file, in
every record — but it never decides the label: an ATP's `PA` is phosphorus
(110), a ligand carbon called `CL1` is carbon (106) and a chlorine called
`CL35` is chlorine (93), whatever the name might suggest. BioTSR takes the
element field as written and corrects nothing; a file with a wrong element
field needs correcting before use.

**Atomic results from pre-release builds that labelled atoms by name must be
computed again.** Those builds labelled atoms through their names, so some
atoms had other labels and every atomic key space had another fingerprint. Those results load, but an
analysis refuses to combine them with current ones and `--if-exists skip`
refuses to reuse them, leaving their files untouched (§13.8). The labels that
did not change — the keys of most calculations are identical — are no reason
to mix them: the key space records what a label means.

| method | what it pools | target options |
|---|---|---|
| `amino_acid_atomic` | the atoms of **one residue**, one calculation per residue | `--residue` with `--residue-name`; `--residue-name` alone for each residue of that name; neither for every residue |
| `amino_acid_atom_selection` | exactly the named atoms, across residues | `--atoms 'RES_SEQ_ATOM;…'` |
| `aa_ind` | every atom of every residue of **one** type | `--residue-types LEU` |
| `aa_ind_mix` | every atom of every residue of **two or more** types | `--residue-types 'LEU;ILE'` |
| `aa_ind_specific` | every atom of named residues | `--residues 'LEU_15;LEU50'` |

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --residue 15 --residue-name LEU --output results/
```

```text
method amino_acid_atomic
triangles 56
keys 46 distinct
```

`--residue` is the residue's identifier in the file: its number, with an
insertion code if it has one (`52A`). Give `--residue-name` alone to run one
calculation per residue of that name, and neither to run every residue of the
chain. Each admitted residue gets its own prefix in the same output folder
(`1UBQ_A_amino_acid_atomic_perresidue-LEU_LEU_15`, …), and the sweep has one
completion record. A residue whose atom count does not match its type is
declined, and the reason is recorded:

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --residue-name LEU --output results/
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --output results/
```

```text
entities 9
output results/1UBQ_A_amino_acid_atomic_perresidue-LEU (28 files, completion record 1UBQ_A_amino_acid_atomic_perresidue-LEU.run.json)
entities 75
```

The selection methods:

```bash
biotsr run --method amino_acid_atom_selection --pdb 1UBQ.pdb --chain A --atoms 'LYS_11_NZ;LEU_15_CB;ILE_13_CG1;LYS_6_NZ' --sample-id ubq_atom_selection --output results/
biotsr run --method aa_ind --pdb 1UBQ.pdb --chain A --residue-types LEU --output results/
biotsr run --method aa_ind_mix --pdb 1UBQ.pdb --chain A --residue-types 'LEU;ILE' --output results/ --no-triplets
biotsr run --method aa_ind_specific --pdb 1UBQ.pdb --chain A --residues 'LEU_15;LEU50' --output results/
```

`--sample-id` sets the output prefix; it is useful when the default prefix,
built from the selection, would be long. It must be a plain file-name stem
(§13.5).

A residue is named by its three-letter name and its identifier in the file:
`LEU_15` and `LEU15` are the same residue (a name containing a digit needs the
`NAME_ID` spelling), and an atom by `RES_ID_ATOM` (`LYS_11_NZ`, or
`LYS_52A_NZ` with an insertion code). A token that is neither is refused, not
skipped.

A selection item that is missing is named in a warning, and a selection that
matches nothing is refused (§12.4).

**Modified amino acids.** RCSB files write modified amino acids as `HETATM`
records. The atomic methods read exactly three of them — phosphothreonine
`TPO`, phosphoserine `SEP` and phosphotyrosine `PTR`, which have exact
heavy-atom counts (11, 10 and 16), labelled by element like every atom — in the amino-acid
atomic methods above and on the protein side of the atomic cross methods
(§10). The alternate-location, identifier and completeness rules apply to them
unchanged, and none is converted to its unmodified residue. Any other modified
residue stored as `HETATM` (`CME`, `MSE`, …) has no defined atom count and is
**not supported** by the atomic methods: it is named in a warning beginning
`not supported:`, never dropped silently, and a request of unsupported residues
only is refused with that reason. The Cα methods (§3, §4) label residues, not
atoms, and read a longer list: the eleven modified residues of the residue
labels (§3), CME and MSE included. A modified residue outside that list is
named as `not represented:`. The two lists are separate on purpose: a residue
label says nothing about how many atoms a complete residue has.

<!-- not-executed -->
```bash
biotsr run --method aa_ind_mix --pdb SYNPHOS.pdb --chain A --residue-types 'TYR;TPO;CME' --output modified/ --no-triplets
biotsr run --method protein --pdb SYNPHOS.pdb --chain A --output modified/ --no-triplets
```

The first command writes its result, prints its warning
(`biotsr: warning: not supported: residue type CME …`) and exits `4`, because
the result is partial (§12.4); the second reads all four modified residues of
chain A as Cα points and exits `0`. From Python:

```python
mixed = biotsr.analyze("SYNPHOS.pdb", method="aa_ind_mix", chain="A",
                       residue_types="TYR;TPO;CME", keep_triplets=False)
assert mixed.warnings == (
    "not supported: residue type CME (a modified residue stored as HETATM; of the "
    "modified amino acids BioTSR reads only TPO, SEP, PTR, which have defined atom "
    "counts)",)
calpha = biotsr.analyze("SYNPHOS.pdb", method="protein", chain="A", keep_triplets=False)
assert calpha.warnings == ()
assert {"3", "5", "7", "8"} <= set(calpha.selection.residues_used)   # TPO SEP PTR CME
```

The selection record (§12.3) lists what the request named that the structure
does not contain under `not_found`, and what it contains but the method cannot
read under `unsupported`:

```python
named = biotsr.analyze("SYNPHOS.pdb", method="aa_ind_specific", chain="A",
                       residues="TPO_3;CME_8;LYS_99", keep_triplets=False)
record = named.annotations["selection_request"]["sides"][0]
assert record["instances"] == ["TPO_3"]
assert record["not_found"] == ["LYS_99"] and record["unsupported"] == ["CME_8"]
```

A residue you named that the completeness rule declines is reported as
`requested but excluded by the method's completeness rule`.

---

## 6. SSE-TSR

`--sse` annotates each Cα triangle with one of **18 secondary-structure
categories**, read from the structure's own `HELIX` and `SHEET` records. A
file without those records has no secondary structure to report. It is
offered by `protein` and `aa_selection`, and refused elsewhere.

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --sse --output sse/
```

```text
sse on, strand identity sheet_id_and_serial
triangles 70300
keys 45349 distinct
```

SSE is an **annotation**: the keys and both identity fingerprints are exactly
those of the run without it. What changes is the category each key's frequency
is counted in. The keys file gains a real `category` column, and the triplets
file gains per-vertex secondary-structure states.

```python
annotated = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", sse=True,
                           keep_triplets=False)
totals = annotated.category_totals()      # {"18_1a1b1c": 14256, ...}
assert sum(totals.values()) == annotated.n_triangles
assert annotated.frequencies == result.frequencies
```

The categories are defined in the [SSE-TSR methods record](methods/SSE_TSR.md).

---

## 7. Drug (ligand) TSR

A ligand copy is identified by its chain, residue name (`--entity-name`) and
residue identifier (`--entity-id`):

```bash
biotsr run --method drug --pdb 3PTB.pdb --chain A --entity-name BEN --entity-id 1 --output results/
```

```text
method drug
triangles 84
keys 38 distinct
```

To see which ligands a chain holds:

```python
targets = biotsr.list_entity_targets("3PTB.pdb", chain="A", method="drug")
assert ("1", "BEN") in targets             # (residue number, residue name)
```

Without `--entity-id`, every copy of the name on the chain is used together, in
one calculation (prefix `…_BEN_allcopies`); 3PTB chain A holds one BEN, so here
that is the same atoms, recorded as an all-copies selection (§12.3).

`drug_atom_selection` uses only the named atoms of the ligand:

```bash
biotsr run --method drug_atom_selection --pdb 3PTB.pdb --chain A --entity-name BEN --entity-id 1 --atoms 'C1;C2;C3;C4;N1;N2' --output results/
```

```text
triangles 20
```

---

## 8. Carbohydrate TSR

One monomer, identified like a ligand. 7LF9 chain D carries one glycan,
NAG 1 – NAG 2 – BMA 3 – MAN 4 – NAG 5 – GAL 6 – MAN 7. With `--entity-id`
the calculation is that one NAG; without it, every NAG on chain D (1, 2 and 5)
in one combined calculation, whose triangles may span the copies:

```bash
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id 1 --output results/
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --output results/
```

```text
method carbohydrate
triangles 364
output results/7LF9_D_carbohydrate_NAG_1
triangles 11480
output results/7LF9_D_carbohydrate_NAG_allcopies
```

Every copy of one sugar is not the whole glycan: `--entity-name MAN` is MAN 4
and MAN 7 only, while `carbohydrate_chain` pools every residue of the chain,
different sugars included.

The whole carbohydrate chain, pooled across residues and sugar types. Add
`--entity-names` to keep only the named residue types:

```bash
biotsr run --method carbohydrate_chain --pdb 5WT9.pdb --chain A --output results/
biotsr run --method carbohydrate_chain --pdb 5WT9.pdb --chain A --entity-names NAG --output nag_only/
```

```text
triangles 8436
triangles 3276
```

Selected atoms of one monomer, or of each copy when the identifier is left out:

```bash
biotsr run --method carbohydrate_atom_selection --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id 1 --atoms 'C1;C2;N2;O5;O6' --output results/
biotsr run --method carbohydrate_atom_selection --pdb 7LF9.pdb --chain D --entity-name NAG --atoms 'C1;C2;N2' --output results/
```

```text
output results/7LF9_D_carbohydrate_atom_selection_NAG_allcopies_atoms-C1+C2+N2
```

---

## 9. Nucleotide TSR

One nucleotide residue:

```bash
biotsr run --method nucleotide --pdb 1BNA.pdb --chain A --entity-name DC --entity-id 3 --output results/
```

```text
method nucleotide
triangles 969
```

A whole DNA or RNA strand is `nucleotide_chain`. Every standard nucleotide of
the chain takes part with the heavy atoms it has — a 5′ nucleotide without its
phosphate, or an incomplete nucleotide inside the chain, included — and each
one whose atoms differ from the standard set is named in a warning; nothing is
reconstructed (ADR-114). A full strand of a 12-mer is already about 2.4 million
triangles, so this example keeps only its cytosines with `--entity-names`:

```bash
biotsr run --method nucleotide_chain --pdb 1BNA.pdb --chain A --entity-names DC --output results/
```

```text
triangles 62196
biotsr: warning: atom sets that differ from the standard nucleotide (as the structure file has them; every atom present is used): DC_1 lacks P, OP1, OP2
```

DC 1, the 5′ end, has 16 of its 19 atoms; they are all used. The warning
describes the file, not the calculation: the result is complete, and its
status is `ok`.

Selected atoms of one residue:

```bash
biotsr run --method nucleotide_atom_selection --pdb 2AC0.pdb --chain E --entity-name DT --entity-id 7 --atoms "P;C1';N1;C2;O4;C7;C6" --output results/
```

---

## 10. Cross TSR

A cross method draws triangles across **two sides**, and every triangle has
vertices on both. An atom joins a side's pool only if an atom of the other
side lies within the contact cutoff. The definitions are in the
[cross-TSR methods record](methods/CROSS_TSR.md).

**Both sides are always named.**

* `--chain` is the first side's chain: the ligand, carbohydrate or nucleotide
  chain, or the first protein chain.
* `--chain-2` is the second side's chain: the protein. It is required by every
  cross method.
* `--entity-name` / `--entity-id` identify the non-protein side by the same
  rule as the single-entity methods (§12.3): with an identifier that one copy,
  without it every copy of the name on `--chain`.

| method | first side | second side | example |
|---|---|---|---|
| `drug_protein_cross` | ligand | protein | 3PTB, BEN 1 in chain A, protein chain A |
| `carbohydrate_protein_cross` | carbohydrate | protein | 7LF9, NAG 1 in chain C, protein chain A |
| `nucleotide_protein_cross` | nucleotide | protein | 2AC0, DG 8 in chain E, protein chain B |
| `protein_protein_cross` | protein (all heavy atoms) | protein | 1BRS, barnase A and barstar D |
| `protein_protein_calpha_cross` | protein (Cα of interface residues) | protein | 1BRS, barnase A and barstar D |

```bash
biotsr run --method drug_protein_cross --pdb 3PTB.pdb --chain A --chain-2 A --entity-name BEN --entity-id 1 --output cross/
biotsr run --method carbohydrate_protein_cross --pdb 7LF9.pdb --chain C --chain-2 A --entity-name NAG --entity-id 1 --output cross/
biotsr run --method nucleotide_protein_cross --pdb 2AC0.pdb --chain E --chain-2 B --entity-name DG --entity-id 8 --output cross/
biotsr run --method protein_protein_cross --pdb 1BRS.pdb --chain A --chain-2 D --contact-cutoff 4.0 --output cross/ --no-triplets
biotsr run --method protein_protein_calpha_cross --pdb 1BRS.pdb --chain A --chain-2 D --contact-cutoff 8.0 --output cross/
```

```text
output cross/3PTB_A_A_drug_protein_cross_BEN_1
output cross/7LF9_C_A_carbohydrate_protein_cross_NAG_1
output cross/2AC0_E_B_nucleotide_protein_cross_DG_8
output cross/1BRS_A_D_protein_protein_cross_cutoff4
output cross/1BRS_A_D_protein_protein_calpha_cross_cutoff8
```

A prefix names both chains in argument order, then the method, the entity, and
any setting you gave (`cutoff4`).

A protein side of the four atomic cross methods reads TPO, SEP and PTR and
names any other modified residue as not supported, on that side
(`not supported: CME_8 (…) (protein side, chain A)`);
`protein_protein_calpha_cross` represents no modified residue, and its
interface rule reads `ATOM` records only (§5, *Modified amino acids*).

**Contact settings.**

* `--contact-cutoff` (Å) decides **which atoms enter the pools**. It is a
  sampling choice: it changes which triangles are drawn, never what a key means.
  Each family has its own default. Results with different cutoffs are compared
  only with `--allow-mixed-processing` (§17).
* It does **not** limit triangle edges. Use `--size-filter` for that.
* `--contact-filter` (`both_sides`, `side_a_only`, `side_b_only`) chooses which
  side the cutoff narrows. `--no-contact-filter` skips the contact scan
  (refused by `protein_protein_calpha_cross`, whose interface is defined by the
  cutoff).

**Narrowing a side.** Before the contact scan, each side can be narrowed:

| argument | narrows | value |
|---|---|---|
| `--atoms` | the entity side of a drug, carbohydrate or nucleotide cross | atom names, `'C1;C2;N1'` |
| `--residue-types` | a protein side to residues of given types (all heavy atoms) | `'ARG;LYS'` |
| `--residues` | a protein side to named residues (all heavy atoms) | `'ASP_189;SER_190'` |
| `--atom-selection` | a protein side to named atoms | `'ARG_59_NH1;LYS_27_NZ'` |
| `--calpha-residues` | a protein side to the Cα of listed residues (`protein_protein_calpha_cross`) | `'50-110'` |
| `--residue-types-2`, `--residues-2`, `--atom-selection-2`, `--calpha-residues-2` | the **second** protein side of a protein–protein cross | as above |

A protein side takes at most one narrowing argument. Which arguments each
cross method takes:

| method | arguments it takes |
|---|---|
| `drug_protein_cross`, `carbohydrate_protein_cross`, `nucleotide_protein_cross` | `--chain-2`, `--contact-cutoff`, `--contact-filter`, `--no-contact-filter`, `--entity-name`, `--entity-id`, `--entity-names`, `--atoms`, `--residue-types`, `--residues`, `--atom-selection` |
| `protein_protein_cross` | `--chain-2`, `--contact-cutoff`, `--contact-filter`, `--no-contact-filter`, `--residue-types`, `--residues`, `--atom-selection` and their `_2` forms |
| `protein_protein_calpha_cross` | `--chain-2`, `--contact-cutoff`, `--calpha-residues`, `--calpha-residues-2` |

Anything else is refused by name: the methods with one protein side refuse the
`_2` forms, and `protein_protein_calpha_cross` takes only the Cα forms.
`--entity-names` (several entity types, pooled) cannot be combined with
`--entity-name`, `--entity-id` or `--atoms`.

```bash
biotsr run --method drug_protein_cross --pdb 3PTB.pdb --chain A --chain-2 A --entity-name BEN --entity-id 1 --residues 'ASP_189;SER_190;GLY_219' --output narrowed/
biotsr run --method protein_protein_cross --pdb 1BRS.pdb --chain A --chain-2 D --residue-types 'ARG;LYS' --residue-types-2 'ASP;GLU' --contact-cutoff 5.0 --output narrowed/ --no-triplets
biotsr run --method protein_protein_calpha_cross --pdb 1BRS.pdb --chain A --chain-2 D --calpha-residues '50-110' --calpha-residues-2 '29-44' --contact-cutoff 8.0 --output narrowed/
```

```text
triangles 63
triangles 5950
triangles 10560
```

In Python, cross methods use `biotsr.analyze_cross`, with the same arguments
spelled with underscores (`chain_2=`, `residue_types_2=`, …):

```python
cross = biotsr.analyze_cross("3PTB.pdb", method="drug_protein_cross",
                             chain="A", chain_2="A",
                             entity_name="BEN", entity_id="1")
assert cross.n_triangles == 84
```

---

## 11. Modifiers

Modifiers are options on the same command, and they combine freely:

| modifier | option | methods | effect |
|---|---|---|---|
| amino-acid grouping | `--aa-grouping` | `protein`, `aa_selection` | 23 grouped residue labels (13 groups of the regular amino acids; MSE with MET; each other modified residue on its own) instead of the individual ones, both at radix 30 (a different key space) |
| mirror | `--mirror` | every non-cross method | negate a key whose angle before folding exceeds 90° (a different key space) |
| size filter | `--size-filter Å` | every method | keep only triangles whose longest edge is below the cutoff (same key space, different processing) |
| SSE annotation | `--sse` | `protein`, `aa_selection` | count keys in 18 secondary-structure categories (§6) |

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --aa-grouping --mirror --size-filter 25.0 --output modified/
```

```text
modifiers aa_grouping=True, mirror=True, size_filter=25.0
```

A modifier a method does not support is refused, with a message naming what it
does support:

```python
try:
    biotsr.analyze("3PTB.pdb", method="drug", chain="A",
                   entity_name="BEN", entity_id="1", sse=True)
except ValueError as error:
    assert "does not support 'sse'" in str(error)
else:
    raise AssertionError("an unsupported modifier must be refused")
```

The contact cutoff of the cross methods (§10) is a processing setting, not a
modifier. Modifiers and settings you give are part of the output prefix
(`1UBQ_A_protein_aagroup_mirror_size25`).

---

## 12. Arguments, entities and selections

### 12.1 One standard profile per method

Every method runs under its own **standard profile**, the accepted BioTSR
science. `biotsr methods NAME` names it, the completion summary prints it, and
the completion record and the metadata sidecar record it:

```bash
biotsr methods drug
```

```text
profile standard_drug_v1
```

There is no `--profile` option. In Python, `profile=` accepts only that
method's own standard profile (`biotsr.get_method(NAME).default_profile`) or
`None`; anything else, including a profile name, is refused:

```python
try:
    biotsr.analyze("1UBQ.pdb", method="protein", chain="A",
                   profile="standard_protein_calpha_v1")
except biotsr.InvalidRequest as error:
    assert "Omit profile=" in str(error)
else:
    raise AssertionError("a profile name must be refused")
```

Native output is the only format, and every native result can be loaded back
(§16).

### 12.2 Arguments are checked per method

An argument the method does not read is refused by name, never silently
ignored, on the command line and in Python alike. `biotsr methods NAME` lists
the target arguments a method takes (its `targets` line):

```python
try:
    biotsr.analyze("1UBQ.pdb", method="protein", chain="A", residues="1-10")
except biotsr.InvalidRequest as error:
    assert "does not use residues" in str(error)
else:
    raise AssertionError("an argument the method does not read must be refused")
```

The cross methods' arguments are listed in §10.

### 12.3 Molecules: one copy, or every copy on the chain

`drug`, `carbohydrate`, `nucleotide`, their `*_atom_selection` variants and the
entity side of the three entity cross methods select molecules by chain, name
and identifier:

* the chain is always required, and only that chain is searched;
* `--entity-name` with `--entity-id` is **that one copy** (scope
  `one_instance`, target and prefix `NAG_1`);
* `--entity-name` alone is **every copy of that name on the chain, in one
  combined calculation** (scope `all_copies`, target and prefix
  `NAG_allcopies`), even when the chain holds one copy. Triangles may span
  copies; that is the intended calculation. The completeness rule applies to
  each copy, and a declined copy is recorded and warned about;
* `--atoms` applies to each selected copy, and an atom a copy lacks is named
  per copy (`requested but not found: XX9 in NAG_1, …`);
* every copy of one sugar is not the whole glycan: `carbohydrate_chain` and
  `--entity-names` pool every residue of a vocabulary, different sugars
  included.

```python
nag = biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D",
                     entity_name="NAG", keep_triplets=False)
assert nag.annotations["target"] == "NAG_allcopies"
side = nag.annotations["selection_request"]["sides"][0]
assert side["scope"] == "all_copies"
assert side["instances"] == ["NAG_1", "NAG_2", "NAG_5"]   # the copies used
assert nag.n_triangles == 11480

ben = biotsr.analyze("3PTB.pdb", method="drug", chain="A",
                     entity_name="BEN", entity_id="1", keep_triplets=False)
assert ben.annotations["target"] == "BEN_1"
```

An identifier is a PDB residue identifier: the residue number as written in
the file, which may be negative, with at most one insertion-code letter
(`404`, `52A`, `-5`). Leading zeros are dropped (`005` is `5`), and `52` never
matches `52A`. A position ("the second copy"), any other spelling, and an
identifier whose records appear in two separate places on the chain are
refused, as is an id or a name the chain does not hold, naming what is there:

```python
try:
    biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D",
                   entity_name="NAG", entity_id="second")
except biotsr.InvalidRequest as error:
    assert "is not a PDB residue identifier" in str(error)
else:
    raise AssertionError("a position is not an identifier")
```

The selection is recorded in the sidecar (`provenance.selection_request`: per
side the chain, scope, arguments, the instances used and those excluded) and
as the provenance `target`, in the sidecar and the completion record. The
record is provenance only: it never decides whether results can be compared
(§17).

For amino acids the same principle holds: `amino_acid_atomic` with `--residue`
and `--residue-name` is that residue, with `--residue-name` alone one
calculation per residue of that name (§5); `aa_ind_specific` names residues as
`NAME_ID` or `NAMEID`; `aa_ind` pools residue types over the chain.

### 12.4 Empty and partial selections

* An explicit selection that matches **nothing** is refused, with a message
  naming the entity or the items.
* Requested items that are **missing** are named in a warning,
  `requested but not found: …`, printed by the command line and recorded in
  the completion record (each result's `warnings`) and on `result.warnings`.
  Items the structure holds but the method cannot read are named the same
  way, as `not supported: …` or `not represented: …` (§5, *Modified amino
  acids*), and a named residue the completeness rule declines as
  `requested but excluded by the method's completeness rule: …`.
* A selection with **fewer than three points** forms no triangle. The result
  is written, with zero triangles and a `no triangles: …` warning.

<!-- not-executed -->
```bash
biotsr run --method aa_selection --pdb 1UBQ.pdb --chain A --residues '1-3;200' --output partial/
```

This prints `biotsr: warning: requested but not found: 200` and
`triangles 1`, writes the result, and exits `4`, because the result is
partial (below).

```bash
biotsr run --method aa_selection --pdb 1UBQ.pdb --chain A --residues '1-2' --output partial/
```

```text
biotsr: warning: no triangles: the selection has 2 points and a triangle needs three
triangles 0
```

```python
three = biotsr.analyze("1UBQ.pdb", method="aa_selection", chain="A",
                       residues="1-3;200", keep_triplets=False)
assert three.n_triangles == 1
assert three.warnings == ("requested but not found: 200",)
two = biotsr.analyze("1UBQ.pdb", method="aa_selection", chain="A",
                     residues="1-2", keep_triplets=False)
assert two.n_triangles == 0
assert two.warnings[0].startswith("no triangles: ")
```

In a batch, a request with no triangles has the status `no_triangles` in
`batch_summary.csv`. It is not a failure and does not change the exit code
(§13.7), unless its selection also lacked something it named: such a result
is partial as well, the summary line counts it
(`2 with no triangles (1 of them partial)`), and the command exits `4`. The completeness rules of the methods that have one (an exact atom
count per residue type, §5) are unchanged; a residue or entity they exclude is
reported, and if that leaves fewer than three points the result carries the
`no triangles` warning.

A result whose selection did not get everything it asked for is **partial**:
a requested residue, atom, molecule or residue type is not found, a modified
residue is not supported or not represented, a residue or molecule copy you
named was declined by the completeness rule, or a selected copy lacks atoms
you named for it. The completeness rule declining a residue among *all*
residues of a type, or one copy among *all* copies of a molecule, does not
make a result partial; that is reported in the skip records. Each side of the
selection record (`provenance.selection_request` in `PREFIX.meta.json`) lists
what was left out as `not_found`, `unsupported` and `declined` (each only when
non-empty), and each result in the completion record has an `omissions` entry,
`{}` when the result is complete.

A partial result is written and usable, but it is never reported as the
complete selection that was asked for. In a batch its status is
**`partial`**: it is written like `ok`, its `message` names what is missing,
and it is counted separately in the summary line. `biotsr run` and
`biotsr batch` exit `4` when a result is partial and nothing failed (§13.7),
and every analysis refuses partial samples unless you pass `--allow-partial`
(§17):

<!-- file: named.csv -->
```csv
pdb_id,chain,method,residues,sample_id
SYNPHOS,A,aa_ind_specific,TPO_3;TYR_6,complete
SYNPHOS,A,aa_ind_specific,TPO_3;CME_8,partly
```

<!-- not-executed -->
```bash
biotsr batch --manifest named.csv --pdb-dir pdb/ --output named/ --no-triplets
```

The command prints `requests 1 ok, 0 skipped, 0 failed, 1 partial` and exits
`4`. The same, from a script that checks the exit status:

```python
import contextlib
import io

from biotsr.cli import EXIT_INCOMPLETE, main

log = io.StringIO()
with contextlib.redirect_stderr(log):
    code = main(["batch", "--manifest", "named.csv", "--pdb-dir", "pdb/",
                 "--output", "named/", "--no-triplets"])
assert code == EXIT_INCOMPLETE == 4
assert "requests   1 ok, 0 skipped, 0 failed, 1 partial" in log.getvalue()
```

```python
import csv

with open("named/batch_summary.csv") as fh:
    status = {r["sample_id"]: (r["status"], r["message"]) for r in csv.DictReader(fh)}
assert status["complete"] == ("ok", "")
assert status["partly"][0] == "partial"
assert status["partly"][1].startswith("not supported: CME_8")

import json

with open("named/partly.run.json") as fh:
    omissions = json.load(fh)["results"][0]["omissions"]
assert omissions["primary"] == {"unsupported": ["CME_8"]}
with open("named/complete.run.json") as fh:
    assert json.load(fh)["results"][0]["omissions"] == {}
```

### 12.5 How a structure is read

* **Model 1 only.** Every method stops reading at the first `ENDMDL` record or
  a `MODEL` record numbered above 1.
* **`TER` ends a chain for protein and nucleotide reading.** The protein Cα,
  amino-acid atomic and nucleotide methods, and the protein side of every
  cross method, stop at the `TER` record of the requested chain. The `HETATM`
  families — `drug` and the carbohydrate methods — and the molecule side of
  every entity cross method read past `TER` to the end of model 1, because
  ligand and glycan records usually follow it. (Pre-release builds before
  ADR-109 read past `TER` on the protein side of `carbohydrate_protein_cross`
  and `nucleotide_protein_cross` too.) A `TER` record ends the chain its chain identifier names; one
  without a chain identifier ends none — check files written by tools that
  write such records.
* **Atom names and elements.** An atom's name is the full PDB atom-name field,
  columns 13–16, and its element the element field, columns 77–78 (never the
  charge in 79–80), read without regard to letter case. Only an atom named
  exactly `CA` (in columns 14–15) is a Cα; calcium and `CA1` never are.
* **`ATOM` records for amino acids, plus three modified ones.** The amino-acid
  and protein reading takes `ATOM` records; the amino-acid atomic methods and
  the protein side of the atomic cross methods also take the `HETATM` records
  of TPO, SEP and PTR, and no other `HETATM` record (§5, *Modified amino
  acids*).
* **Alternate locations, one rule for every method.** An atom whose altLoc is
  blank or `A` is kept; `B`, `C` and every later indicator are excluded. There
  is no fallback to another conformer when `A` is absent: a residue whose side
  chain exists only as `B` or `C` loses those atoms, and a method with a
  completeness rule then excludes the residue, with the reason recorded. Every
  result records `source.altloc_policy = "blank_or_a"`, and excluded atoms
  leave a result warning (printed with `--verbose`), for example on 1BRS chain
  D, whose SER 14 has two conformers:
  `alternate locations: kept blank and A; excluded 2 atom record(s) with indicator B`.
  To study another conformation, prepare a PDB file whose altLoc indicators
  select it: keep only the wanted conformer, relabelled `A` or blank.
* **Hydrogen and deuterium** are never vertices, and never contact atoms.
* **Metal ions** (Mn, Zn, …) are labelled as their elements (manganese 100,
  zinc 105), whether written as `HETATM` or, in a non-conformant file, as
  `ATOM` records.
* **What is refused rather than computed** (ADR-108, ADR-109). A calculation is
  refused, with the records named and what to do, when what it selects
  contains: a residue named `UNK` (no calculation is done for UNK — narrow the
  selection to leave it out); an atom with a blank element field (except a
  `HETATM` Cα); an atom written twice, or two atoms at exactly the same
  coordinates; a record whose residue number cannot be read (hybrid-36); or an
  atom whose element the element table does not list (`X`). Records outside the selection are never
  checked, so a problem in another chain does not stop the calculation. A
  protein–protein cross of a chain with itself, and a cross whose two sides
  share an atom, are refused too; a molecule side may share the protein's
  chain identifier.
* **Terminal residues** (ADR-111). Cα methods keep terminal Cα points. The
  per-residue and residue-type atomic methods keep their exact atom count, so a
  C-terminal residue with `OXT` (1UBQ GLY 76) or a 5′ nucleotide without its
  phosphate is excluded and named. Cross methods use terminal atoms, `OXT`
  included, like any other. `nucleotide_chain` has no count: every nucleotide
  takes part with the atoms it has, and one whose atoms differ from the
  standard set is named in a warning (ADR-114).

---

## 13. Datasets: manifests, sheets and batches

### 13.1 The manifest

A **manifest** is a CSV file in which **one row is one request**. Its columns
are the `biotsr run` options without the leading dashes, and the file can have
any name. (A wide dataset sheet, one row per structure, is read with `--sheet`,
§13.9.)

<!-- file: samples.csv -->
```csv
sample_id,pdb_id,method,chain,label
ubiquitin,1UBQ,protein,A,ubiquitin
barnase,1BRS,protein,A,barnase
barstar,1BRS,protein,D,barstar
```

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --no-triplets
```

```text
requests 3 ok, 0 skipped, 0 failed
summary batch/batch_summary.csv
```

Every calculation is written directly into `batch/`, under its prefix — here
the sample ids: `ubiquitin.keys`, `ubiquitin.meta.json`, `ubiquitin.run.json`,
and the same for `barnase` and `barstar` (§15).

| column | required? | meaning |
|---|---|---|
| `pdb_path` or `pdb_id` | one of them | `pdb_path` is a path used as given; `pdb_id` is looked up in `--pdb-dir` as `<id>.pdb` or `<id>.ent`, in any case |
| `chain` | **yes** (or `--chain`, §13.3) | the chain the run is scoped to |
| `method` | yes, unless `--method` is given (never both) | one of the names `biotsr methods` lists |
| `sample_id` | no | the output prefix (§13.5); a deterministic one is made if blank |
| `label` | no | free text: a class or group name, used by `--group-by label` |
| `chain_2` | cross methods | the second chain |
| method options | as the method needs | `residue`, `residues`, `residue_types`, `atoms`, `entity_name`, `entity_id`, `entity_names`, `mirror`, `aa_grouping`, `size_filter`, `sse`, `contact_cutoff`, the cross narrowing columns, … |

A blank cell means "not set", never a default. There is no `profile` column:
every method runs its standard profile, and a sheet that has one is refused
until you drop it or pass `--ignore-column profile`. Each row carries its own
method, so one manifest can mix methods and families:

<!-- file: mixed.csv -->
```csv
sample_id,pdb_id,method,chain,chain_2,entity_name,entity_id,mirror,label
ubq_calpha,1UBQ,protein,A,,,,,ubiquitin
ubq_mirror,1UBQ,protein,A,,,,true,ubiquitin
trypsin_ben,3PTB,drug,A,,BEN,1,,ligand
trypsin_ben_cross,3PTB,drug_protein_cross,A,A,BEN,1,,interface
```

```bash
biotsr batch --manifest mixed.csv --pdb-dir pdb/ --output mixed/ --no-triplets
```

```text
requests 4 ok, 0 skipped, 0 failed
```

The methods may instead be chosen on the command line: `--method NAME`, which
is repeatable, runs each listed method on every row. A manifest whose rows
name methods together with `--method` is refused, naming the remedies (drop
`--method`, or `--ignore-column method`). The [CLI guide](CLI_GUIDE.md), §7,
lists every column.

**How the CSV is read.**

* UTF-8, with or without a byte-order mark (as spreadsheet programs often
  write).
* A non-blank cell beyond the header, or under an empty header, is refused,
  naming the line: it would otherwise be dropped silently. A list inside one
  cell separates its items with `;` (`1-10;20`), or the cell is quoted.
* Two columns with the same header name are refused.

### 13.2 Columns with other names, kept and ignored columns

An unknown column is an error until you give it a role; columns whose names
start with `comment` are ignored. Column names used in the historical TSR
dataset sheets are accepted as aliases (`protein` for `pdb_id`, `aa_chain` for
`chain`, `aa_group` for `label`, `drug_name` for `entity_name`, …; the full
list is in `biotsr.io.input_columns` and in the error message). For any other
column:

| option | meaning |
|---|---|
| `--column FIELD=HEADER` | read manifest field `FIELD` from the column `HEADER` |
| `--keep-column HEADER` | carry the column into each completion record and `batch_summary.csv` as **metadata**; it never reaches a calculation |
| `--ignore-column HEADER` | drop the column |

A dataset sheet whose `group` column is a class label and `group1` its numeric
code:

<!-- file: sheet.csv -->
```csv
protein,chain,group,group1
1UBQ,A,alpha_beta,0
1BRS,A,alpha_beta,0
1BRS,D,alpha_beta,0
```

```bash
biotsr batch --manifest sheet.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --output sheet_out/ --no-triplets
```

```text
requests 3 ok, 0 skipped, 0 failed
```

`--column label=group1 --keep-column group` works equally well. A label is an
annotation: it never turns on `--aa-grouping` or any other scientific option.

### 13.3 Command-line values are row defaults; chains are explicit

Options given on the command line (`--mirror`, `--size-filter`, `--chain A`)
apply to rows that leave that column blank; anything a row states wins.
`--method` is different: it chooses the methods for every row (§13.1). **Every request needs its chain** — and a cross method its second chain —
from the row, a mapped column (`--column chain=HEADER`, `--map`) or `--chain`
for rows that leave it blank. BioTSR never assumes chain A or the first chain
of a structure. A request whose chain is still blank is not run: it is listed
in `batch_summary.csv` with the status `chain_not_specified` and the reason,
the other requests run, and the batch exits 4 (ADR-112). A sheet with no chain
column needs `--chain`:

<!-- file: ids_only.csv -->
```csv
protein,notes
1UBQ,ubiquitin
1BRS,barnase
```

```bash
biotsr batch --manifest ids_only.csv --method protein --chain A --ignore-column notes --pdb-dir pdb/ --output chain_default/ --no-triplets
```

```text
requests 2 ok, 0 skipped, 0 failed
```

`batch_summary.csv` records where each chain came from (`chain_source` is
`row` or `default`):

```python
import csv

with open("chain_default/batch_summary.csv") as fh:
    rows = list(csv.DictReader(fh))
assert [r["chain_source"] for r in rows] == ["default", "default"]
```

### 13.4 Dry run

`--dry-run` plans a batch and changes nothing: no download, no output
folder, no file. It reports where chains come from, which columns are read,
kept or unused, which structures are present or would be downloaded, and
every request with its source rows, method, participant, chains, selection,
output prefix and outcome (would run, would be skipped, inapplicable, or would
fail with the reason):

```bash
biotsr batch --manifest sheet.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --download-missing --output planned/ --no-triplets --dry-run
```

```text
dry run nothing will be computed, downloaded or written
chains 3 from the row, 0 from --chain
structures 3 present, 0 distinct would be downloaded
outcome 3 would run, 0 would be downloaded then run, 0 would be skipped (verified), 0 inapplicable, 0 chain not specified, 0 would fail
1 2 protein - A - 1UBQ_A_protein [would run]
```

```python
import os

assert not os.path.exists("planned")
```

It exits `0` when every request would run and `3` when some would fail.
Checks that need the structure itself (that the chain exists in it, that a
selection matches its atoms) are listed as deferred.

### 13.5 Output names and duplicate outputs

Every calculation writes its files directly into `--output` under one
**prefix**: `PREFIX.keys`, `PREFIX.meta.json`, `PREFIX.triplets` (unless
`--no-triplets`) and the completion record `PREFIX.run.json` (§15). The prefix
is the row's `sample_id`, or else a deterministic name built from the request:
structure, chain(s), method, target and the modifiers and settings that change
the numbers — `1UBQ_A_protein`, `1UBQ_A_protein_mirror`,
`7LF9_D_carbohydrate_NAG_1`, `7LF9_D_carbohydrate_NAG_allcopies`,
`1BRS_A_D_protein_protein_cross_A-ARG+LYS_D-ASP+GLU_cutoff5`. When the
readable form cannot represent the request exactly (an unusual character, a
lower-case chain) or is longer than 96 characters, a request fingerprint
`_h<10 hex digits>` is appended. The [CLI guide](CLI_GUIDE.md), §7, *Output
names*, has the full grammar.

A `sample_id` must be **one plain file-name stem**:

* no `/` or `\`, no absolute path, not `.` or `..`, and no leading `.`
  (reserved for hidden and staging names);
* not empty, no leading or trailing spaces, no control characters, and at
  most 200 bytes;
* not a name the output folder uses for itself, in any letter case:
  `analysis`, `batch_summary.csv`, `batch_failures.csv`, `fetch_report.csv`,
  `.biotsr-batch.lock`.

An invalid id fails its request (exit `2` for `biotsr run`). The check is made
while planning, so `--dry-run` reports it without writing anything.

**Duplicate samples and duplicate computations.** A sample is one entry of the
output folder: its prefix and its label. Two rows are merged only when they
ask for **the same calculation of the same sample** — method, structure,
chains, selection, modifiers and sample id all equal, the sample id being the
row's `sample_id` or, when that is blank, the automatic prefix. The merged
request is computed once, the summary lists every source row, and its message
says `one sample from 2 source requests (line …)`. Rows with **distinct
`sample_id` values** are distinct samples even when their calculation is
identical: each is computed and written under its own prefix, and an analysis
reads both. Rows that would merge but carry **different labels** are refused
before anything runs (one output carries one label; give them distinct
`sample_id` values). A `--keep-column` whose values differ between merged rows
records each row's value (`line 2=x; line 3=y`).

<!-- file: samples_twice.csv -->
```csv
pdb_id,chain,method,residue_types,sample_id,label
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_a,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_b,ubiquitin
```

```bash
biotsr batch --manifest samples_twice.csv --pdb-dir pdb/ --output samples_twice/ --no-triplets
```

```text
requests 3 ok, 0 skipped, 0 failed
```

```python
with open("samples_twice/batch_summary.csv") as fh:
    twice = list(csv.DictReader(fh))
assert [(r["rows"], r["sample_id"]) for r in twice] == [
    ("2;3", "1UBQ_A_aa_ind_LYS"), ("4", "ubq_tyr_a"), ("5", "ubq_tyr_b")]
assert "one sample from 2 source requests (line 2, 3)" in twice[0]["message"]
```

Two **different** requests that would write the same
prefix are refused before anything runs. Prefixes that differ only in letter
case count as the same, because they are the same files on macOS's default
file system:

<!-- file: dup.csv -->
```csv
sample_id,pdb_id,method,chain
barnase,1BRS,protein,A
Barnase,1BRS,protein,D
```

```
$ biotsr batch --manifest dup.csv --pdb-dir pdb/ --output dup_out/
biotsr: error: 1 output prefix(es) in dup.csv would each be written by more than one different request: barnase (rows 2, 3)
```

From Python the same check raises `batch.BatchError`:

```python
from biotsr import batch

try:
    batch.run("dup.csv", "dup_out/", pdb_dir="pdb/", triplets=False)
except batch.BatchError as error:
    assert "more than one different request" in str(error)
else:
    raise AssertionError("two requests writing one prefix must be refused")
assert not os.path.exists("dup_out")
```

### 13.6 One request at a time, or several

A batch computes **one request at a time** by default (`--jobs 1`). `--jobs N`
computes N requests at a time in worker processes, and `--jobs auto` uses the
CPUs this process may use: the scheduler's allocation (Slurm, PBS, LSF, Grid
Engine), its CPU affinity, its cgroup quota or the machine's CPU count,
whichever is smallest. More workers than requests are never started. joblib,
installed with BioTSR, runs the workers; nothing else is needed.

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch_parallel/ --no-triplets --jobs 2
```

```text
workers 2 processes
requests 3 ok, 0 skipped, 0 failed
```

The results are identical to `--jobs 1`, file for file, and the summary keeps
manifest order:

```python
for prefix in ("ubiquitin", "barnase", "barstar"):
    for suffix in (".keys", ".meta.json"):
        with open("batch/" + prefix + suffix, "rb") as serial, \
                open("batch_parallel/" + prefix + suffix, "rb") as parallel:
            assert serial.read() == parallel.read()
```

Each worker holds one request's calculation, so peak memory is about N times
that of one request (§14). `--download-workers` sets download concurrency
separately.

### 13.7 Errors and failed requests

Mistakes in the request as a whole (an unknown column, a row without a chain,
duplicate outputs, a method chosen twice, an invalid download setting) print
one line starting `biotsr: error:` and exit `2` before anything runs. Within a
batch, every request is attempted by default: a failing request is recorded in
`batch_failures.csv` against its manifest line(s), it leaves no files in the
output folder, and the command exits `3`. `--on-error stop` halts at the first
failure instead, writes the summary, and exits `3` (or `2` when a request could
not even be planned, before anything was computed). A `partial` or
`no_triangles` request (§12.4) and an `inapplicable` one (§13.9) are not
failures; a `partial` one makes the command exit `4` when nothing failed.

| exit code | meaning |
|---|---|
| `0` | success |
| `1` | an internal error: a defect in BioTSR (`--debug` prints the traceback) |
| `2` | a user or input error |
| `3` | a batch in which some request failed, was stopped by a worker-pool failure or never started; or `biotsr fetch` with an id not downloaded, refreshed or cached |
| `4` | completed and nothing failed, but at least one result is partial (§12.4): `biotsr run` or `biotsr batch` wrote a partial result (a `partial` request, or a `no_triangles` one whose selection also lacked something it named) or, with `--if-exists skip`, kept one; or `biotsr analyze` analysed partial samples under `--allow-partial`. The results are written and usable |
| `130` | interrupted (Ctrl-C or SIGTERM) |

`3` takes precedence over `4`: a batch with a failed request exits `3` even if
another request is partial. An analysis that refuses a set exits `2`.
`--dry-run` is unaffected. Until exit code `4` was introduced a partial result
exited `0`, so a script that stops on any non-zero status now stops on a
partial result; accept `4` where a partial result is acceptable.

### 13.8 Resume, interruption and concurrent runs

An existing calculation is never overwritten silently. `--if-exists skip`
resumes a batch, keeping an existing calculation only when its completion
record is complete and records exactly the same request (structure content,
method, chains, selections, modifiers, contact settings) computed under the
same processing rules, every file it lists is present at the size it was
written, and it holds what is asked for now:

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --no-triplets --if-exists skip
```

```text
requests 0 ok, 3 skipped, 0 failed
```

A calculation that fails the check is reported as a failed request with the
difference named, never silently kept or replaced. For example, outputs written
with `--no-triplets` do not satisfy a request for triplets:

```python
report = batch.run("samples.csv", "batch/", pdb_dir="pdb/", if_exists="skip")
assert report.n_failed == 3
assert "written with --no-triplets" in report.failures[0]["message"]
```

`--if-exists overwrite` replaces the calculation's files, and removes files of
the older version that the new one no longer produces.

**Processing rules.** The same request can produce different numbers when the
rules that admit atoms change. The amino-acid atomic methods, and the protein
sides of the atomic cross methods, now read the modified amino acids TPO, SEP
and PTR written as `HETATM`; outputs written by earlier builds of the same
version left them out. Every completion record therefore holds the rules it
was computed under, in `request.processing`
(`{"altloc_policy": "blank_or_a", "modified_residue_policy": "hetatm_amino_acids:PTR,SEP,TPO"}`
for an atomic method; a Cα method records
`hetatm_amino_acids:CME,CSD,CSO,LLP,M3L,MSE,NIY,PTR,SEP,TPO,TRF;ca_element:C_or_blank`, the modified
residues it reads), and `--if-exists skip` keeps a calculation only when they equal the current
rules; the package version is not used for this. A record that does not state
them is never taken to be current. The request fails with
`output … exists but cannot be reused: its completion record does not state
the processing rules it was computed under …`, and a record that states other
rules with `… it was computed under different processing rules (…)`.
Recompute such outputs with `--if-exists overwrite` (or into a new
`--output`); they can still be analysed deliberately with
`--allow-mixed-processing` (§17).

**The batch reports describe the latest run.** `batch_summary.csv` (and
`batch_failures.csv`, when a request failed) belong to the output folder, not
to a calculation: every `biotsr batch` into the folder rewrites them for its
own requests, and removes a `batch_failures.csv` left by an earlier run when
nothing failed this time. `--if-exists skip` never changes a calculation's
own files — keys, triplets, sidecar and completion record — whether it keeps
them or refuses them. A `--dry-run` writes nothing at all. Copy the reports
first if you need the earlier run's.

`skip` also requires each result's key space to be the one this version
computes for the request. A Cα output of a pre-release build with residue
labels at radix 20, and an atomic output of a pre-release build that labelled
atoms by name, therefore fails with `… its keys were computed in a different
key space …`, alone or followed by the processing message above; it cannot be
analysed with current results either, so compute it again (§3, §5). An output
whose key space and recorded processing rules are this version's is kept,
whichever build wrote it.

**A kept partial result stays partial.** When `--if-exists skip` keeps an
output that is itself partial (§12.4), the request is `skipped`, its message
says what the kept result lacks, the summary line counts it
(`2 skipped (1 of them partial)`), and the command exits `4`:

```python
report = batch.run("named.csv", "named/", pdb_dir="pdb/", triplets=False,
                   if_exists="skip")
assert (report.n_skipped, report.n_skipped_partial) == (2, 1)
kept = {r["sample_id"]: r["message"] for r in report.records}
assert kept["complete"] == "verified existing output of the same request"
assert kept["partly"].startswith("verified existing output of the same request; "
                                 "kept partial result: not supported: CME_8")
```

**What is never touched.** BioTSR writes, replaces and removes only files that a
completion record of that prefix lists. Your own files in the output folder —
and a file whose name differs from a planned one only in letter case — are
left alone, even with `--if-exists overwrite`: a calculation that would have to
replace one is refused, naming it. Symbolic links are never written through.

**How a calculation is installed.** It is computed into a staging directory
private to the invocation (`.biotsr-staging-…` inside `--output`), then moved
into place under a per-prefix lock: first an intent record
(`"complete": false`), then each file, and last the completion record with
`"complete": true` and every file's size. A calculation is complete exactly
when its record says so.

**Interruption.** Ctrl-C, or SIGTERM from a scheduler's time limit, stops the
batch at any point after planning — while downloading, while checking existing
outputs or while computing. It removes only its own unfinished staging
directories, writes the summary with `interrupted` and `not_started` requests
(and `fetch_report.csv` if it was downloading), and exits `130`.

**A worker that dies.** If a worker process is killed (for example by the
out-of-memory killer), the summary is still written: finished requests keep
their results, requests that were running are `failed` with a message naming
the pool failure, requests never started are `not_started`, and the batch exits
`3`.

Either way, re-run with `--if-exists skip` to finish the rest:

<!-- not-executed -->
```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --jobs auto --if-exists skip
```

**What cannot be caught.** SIGKILL and power loss stop the process before it
can write anything more. What it leaves is a staging directory, which a later
run on the same machine removes once that process is gone (ones left by
another machine are kept and reported), or an incomplete completion record,
which every reader refuses. `--if-exists skip` or `overwrite` computes such a
calculation again; the default `fail` refuses it, naming it.

**Concurrent runs.** Several `biotsr run` invocations may share an `--output`;
if two finish the same prefix, the first stays and the second reports the
clash. Only one `biotsr batch` at a time may write into an output folder (a
second is refused before computing); give concurrent batches their own
`--output`. See the [CLI guide](CLI_GUIDE.md), §7, *Resuming*, and §12 for
SLURM examples.

### 13.9 A wide dataset sheet

A sheet with **one row per structure and many calculations per row** —
several chain columns (`aa_chain`, `drug_chain`, `carb_chain`, `nt_chain`,
`aa_chain_2`) and several alternative selections, as in the historical TSR
dataset sheets — is read with `--sheet`. Which calculations to run is not in
the sheet: you choose them with `--method`, which is required and repeatable,
and each method reads **only its own columns**:

<!-- file: wide.csv -->
```csv
protein,aa_chain,drug_chain,carb_chain,aa_group,aa_ind,drug_name,drug_id,carb_name,carb_id,comment_1
3PTB,A,A,,protease,,BEN,1,,,a ligand row
7LF9,B,,D,glycoprotein,TYR,,,NAG,1,one NAG of the glycan
7LF9,B,,D,glycoprotein,TYR,,,NAG,,every NAG on chain D
1BRS,A,,,barnase,LYS,,,,,no glycan
```

```bash
biotsr batch --manifest wide.csv --sheet tsr --method carbohydrate --method aa_ind --inapplicable skip --pdb-dir pdb/ --output from_wide/ --no-triplets
```

```text
manifest wide.csv (4 source rows, 7 requests)
unused drug_chain, drug_id, drug_name (read by none of the selected methods)
requests 4 ok, 0 skipped, 0 failed, 3 inapplicable
```

```python
with open("from_wide/batch_summary.csv") as fh:
    planned = [(r["rows"], r["method"], r["status"], r["sample_id"])
               for r in csv.DictReader(fh)]
assert planned == [
    ("2", "carbohydrate", "inapplicable", ""),
    ("2", "aa_ind", "inapplicable", ""),
    ("3", "carbohydrate", "ok", "7LF9_D_carbohydrate_NAG_1"),
    ("3;4", "aa_ind", "ok", "7LF9_B_aa_ind_TYR"),
    ("4", "carbohydrate", "ok", "7LF9_D_carbohydrate_NAG_allcopies"),
    ("5", "carbohydrate", "inapplicable", ""),
    ("5", "aa_ind", "ok", "1BRS_A_aa_ind_LYS"),
]
```

Here `carbohydrate` read `carb_chain`, `carb_name` and `carb_id`, and `aa_ind`
read `aa_chain` and `aa_ind`; the ligand columns were read by neither and are
reported as unused. Row 4 names NAG without an identifier, so it is every NAG
on chain D (§12.3). Rows 3 and 4 ask for the same `aa_ind` calculation, which
is computed once and lists both rows. A row that does not describe a method's
calculation — the 3PTB row for `carbohydrate`, a row with an empty `aa_ind`
for `aa_ind` — fails by default; `--inapplicable skip` records it as
`inapplicable` instead. It is never silently left out. A row that describes a
calculation but leaves its chain blank fails, naming the column.

`--sheet tsr` reads the TSR dataset-sheet columns:

| method | chain column | other columns |
|---|---|---|
| `protein`, `amino_acid_atomic` | `aa_chain` | — |
| `aa_selection` / `amino_acid_atom_selection` | `aa_chain` | `aa_selection` / `aa_atom_selection` |
| `aa_ind` / `aa_ind_mix` / `aa_ind_specific` | `aa_chain` | `aa_ind` / `aa_ind_mix` / `aa_ind_specific` |
| `drug`, `drug_atom_selection` | `drug_chain` | `drug_name`, `drug_id` (+ `drug_atom_selected`) |
| `carbohydrate`, `carbohydrate_atom_selection` | `carb_chain` | `carb_name`, `carb_id` (+ `car_atom_selected`) |
| `nucleotide`, `nucleotide_atom_selection` | `nt_chain` | `nt_name`, `nt_id` (+ `nt_atom_selected`) |
| `carbohydrate_chain` / `nucleotide_chain` | `carb_chain` / `nt_chain` | — |
| the three entity cross methods | the entity's chain column, and `aa_chain` | the entity's name and id columns; each side's narrowing column as chosen with `--cross-side` (below) |
| `protein_protein_cross`, `protein_protein_calpha_cross` | `aa_chain` and `aa_chain_2` | each side's narrowing column as chosen with `--cross-side` |

`protein` is the structure and `aa_group` the label (an annotation only; it
never turns on `--aa-grouping` and never chooses a method). Three options do
three different things: `--method` chooses the calculations, `--sheet tsr`
chooses this known column layout, and `--map`, `--map-file` and
`--sheet custom` adapt other headers. A `method` column and `--chain` are
refused with `--sheet`.

**Your own headers.** `--map METHOD[.PARTICIPANT].FIELD=HEADER` (repeatable)
and `--map-file` add or override entries of the layout, `--column FIELD=HEADER`
maps the fields every request shares (`pdb_id`, `label`, …), and
`--sheet custom` starts from no layout and reads only what you map:

<!-- file: glycans.csv -->
```csv
entry,glycan_chain,sugar,sugar_id,group
7LF9,D,NAG,1,glycoprotein
7LF9,D,MAN,,glycoprotein
```

```bash
biotsr batch --manifest glycans.csv --sheet custom --method carbohydrate --column pdb_id=entry --column label=group --map carbohydrate.chain=glycan_chain --map carbohydrate.entity_name=sugar --map carbohydrate.entity_id=sugar_id --pdb-dir pdb/ --output custom/ --no-triplets
```

```text
requests 2 ok, 0 skipped, 0 failed
```

**Participants.** The amino-acid methods have a second participant, the `_2`
columns. `--participant 1` (the default) reads `aa_chain` and the
participant-1 columns (`aa_ind`, `aa_selection`, …), `--participant 2` reads
`aa_chain_2` and the `_2` columns (`aa_ind_2`, `aa_selection_2`, …), and
`--participant all` reads both, as separate requests. A row that fills a `_2`
selection but leaves `aa_chain_2` blank fails (`aa_ind (participant 2) needs
chain, and the row leaves aa_chain_2 blank; BioTSR never infers a chain`). A
cross method has one participant, which covers both of its sides.

**Cross methods: choose each side.** A cross method reads its chains and its
molecule's name and identifier from the row. Each side can also be narrowed by
several columns of the sheet — the protein side of `drug_protein_cross` by
`aa_ind`, `aa_ind_mix`, `aa_ind_specific` or `aa_atom_selection`, the drug by
`drug_atom_selected` — and the sheet does not say which one a cross
calculation means. BioTSR neither picks one nor ignores them: while a row
populates any of a side's columns, the run is refused until you choose.

| method | sides | candidate columns |
|---|---|---|
| `drug_protein_cross`, `carbohydrate_protein_cross`, `nucleotide_protein_cross` | `drug` / `carbohydrate` / `nucleotide`, and `protein` | `drug_atom_selected` / `car_atom_selected` / `nt_atom_selected`; `aa_ind`, `aa_ind_mix`, `aa_ind_specific`, `aa_atom_selection` |
| `protein_protein_cross` | `1` and `2` | the `aa_*` columns above; their `_2` forms |
| `protein_protein_calpha_cross` | `1` and `2` | `aa_selection`; `aa_selection_2` |

`--cross-side METHOD.SIDE=COLUMN` narrows a side by that column, and
`--cross-side METHOD.SIDE=whole` uses the whole side; a `--map` of the side's
field also counts as the choice. This sheet's 3PTB row fills two protein-side
columns:

<!-- file: wide_cross.csv -->
```csv
protein,aa_chain,drug_chain,aa_chain_2,aa_group,aa_ind,aa_ind_specific,aa_ind_2,drug_name,drug_id
3PTB,A,A,,protease,TYR,ASP_189;SER_190;GLY_219,,BEN,1
1BRS,A,,D,barnase-barstar,LYS,,ASP,,
```

```python
try:
    batch.run("wide_cross.csv", "wide_cross/", pdb_dir="pdb/", sheet="tsr",
              methods=["drug_protein_cross"], dry_run=True)
except batch.BatchError as error:
    assert ("the protein side of drug_protein_cross: aa_ind on line 2; "
            "aa_ind_specific on line 2") in str(error)
else:
    raise AssertionError("an unchosen cross side must be refused")
```

Choose it, and preview: under each request line, one indented line per side
gives its resolved selection, and the populated columns the choice leaves
unused are named:

```bash
biotsr batch --manifest wide_cross.csv --sheet tsr --method drug_protein_cross --cross-side drug_protein_cross.protein=aa_ind_specific --inapplicable skip --pdb-dir pdb/ --output wide_cross/ --no-triplets --dry-run
```

```text
drug, chain A: BEN copy 1, all atoms
protein, chain A: residues=ASP_189;SER_190;GLY_219
not applied: aa_ind (the protein side uses aa_ind_specific)
```

Computed, the summary's `selection` holds both sides, `side text | side text`,
and its `message` the columns not applied:

```bash
biotsr batch --manifest wide_cross.csv --sheet tsr --method drug_protein_cross --cross-side drug_protein_cross.protein=aa_ind_specific --inapplicable skip --pdb-dir pdb/ --output wide_cross/ --no-triplets
```

```python
with open("wide_cross/batch_summary.csv") as fh:
    (cross_row,) = [r for r in csv.DictReader(fh) if r["status"] == "ok"]
assert cross_row["selection"] == ("drug, chain A: BEN copy 1, all atoms | "
                                  "protein, chain A: residues=ASP_189;SER_190;GLY_219")
assert cross_row["message"] == "not applied: aa_ind (the protein side uses aa_ind_specific)"
assert cross_row["n_triangles"] == "63"
```

The drug side needed no choice: no row fills `drug_atom_selected`. A row whose
chosen column is blank does not describe the narrowed calculation: it is
inapplicable, never computed on the whole side. In Python the choice is
`cross_sides={"drug_protein_cross.protein": "aa_ind_specific"}` of
`batch.run`. The [CLI guide](CLI_GUIDE.md), §7, *Wide dataset sheets* and
*Cross methods in a wide sheet*, has the details and the mapping rules.

---

## 14. Memory and disk

What dominates a run's memory is its per-triangle (triplet) records, and
whether they are held depends on how it runs:

* **streamed** (each line written as its triangle is computed): `biotsr run`,
  `biotsr batch` and `batch.run(...)` from Python, for every method except the
  `amino_acid_atomic` sweep without `--residue`; cross methods included;
* **never built**: `--no-triplets`, and `biotsr analyze common-keys
  --manifest`, which writes keys only;
* **held until the run finishes**: the `amino_acid_atomic` sweep over every
  residue, and `biotsr.analyze(...)`, `biotsr.analyze_cross(...)` and
  `biotsr.analyze_entities(...)` with the default `keep_triplets=True` (pass
  `keep_triplets=False` and a `triplet_sink` to stream; API guide,
  *Very large results*).

What is always in memory is small by comparison: one run's key counts (one
entry per distinct key, not per triangle), and in a batch one summary line per
request. The key file, its `.meta.json` sidecar, the completion record,
`batch_summary.csv` and `batch_failures.csv` are written whole when their
calculation or batch finishes.
Analysis is different: `biotsr analyze common-keys` and `biotsr.analysis` hold
every sample's key counts at once, so their memory grows with the number of
samples times the distinct keys per sample.

Streaming saves memory, **not disk**: a triplet file is exactly as large
either way, and it is by far the largest output, so for a large dataset plan
disk for triplet files or use `--no-triplets`. Triangle counts grow with the
cube of the number of points, so sizes projected from a few structures are a
rough guide only. Measured peaks for one large chain, and the file sizes they
wrote, are in the [CLI guide](CLI_GUIDE.md), §7, *Memory and disk*; the
remaining limits are in [known limitations](KNOWN_LIMITATIONS.md), §5.

---

## 15. Output files

`--output` is one shared folder. Every calculation writes its files directly
into it, under one prefix (§13.5). For `results/1UBQ_A_protein`:

| file | contents |
|---|---|
| `1UBQ_A_protein.keys` | the keys: `key`, `category`, `frequency`, one row per key and category, sorted by key |
| `1UBQ_A_protein.triplets` | one row per triangle: each vertex's identity and label, the angle and distance classes, raw geometry, coordinates and key (omit with `--no-triplets`; columns below) |
| `1UBQ_A_protein.meta.json` | the metadata sidecar: key space (with both bin schemes), sampling profile, both fingerprints, the selection record, source-file SHA-256, counts, the target used |
| `1UBQ_A_protein.run.json` | the **completion record**: the exact command, BioTSR version, profile, parameters, source rows of a batch request, the normalised request, every file with its size, and the warnings |

`.keys` and `.triplets` are tab-separated UTF-8 text with one header line; the
suffix names the content, not a different format. The file names carry no bin
counts: the sidecar records the full angle and distance bin definitions.

A batch adds `batch_summary.csv` (and `batch_failures.csv`, and
`fetch_report.csv` with `--download-missing`) to the same folder. Nothing else
in the folder is touched: files a completion record does not list are never
replaced or removed (§13.8).

**The completion record** is written last. A calculation is complete only
when its record says `"complete": true`; an interrupted run leaves either no
record or an incomplete one, which `describe` reports and `analyze
--from-results` refuses.

**Fingerprints.** The sidecar records two fingerprints. The `KeySpace`
fingerprint identifies what the integer keys mean; the `SamplingProfile`
fingerprint identifies how the triangles were drawn, from what. Results can be
compared only if their key spaces match (§17).

For a per-residue sweep (`amino_acid_atomic` without `--residue`), each
residue's files carry the residue after the sweep's prefix
(`1UBQ_A_amino_acid_atomic_perresidue_LEU_15.keys`), and the sweep has one
completion record (`1UBQ_A_amino_acid_atomic_perresidue.run.json`).

**Older outputs.** Results written by early pre-release builds — one directory per
sample with `run.json` inside, and names such as `1UBQ_A.keys_29_35` — are not
read by this version: `describe`, `analyze --from-results` and
`biotsr.load_result` refuse them with a message saying so. Recompute them.

**The triplet file** is tab-separated UTF-8 text, with a header line and one
line per triangle. Its 47 columns are the same for every method:

| columns | contents |
|---|---|
| `key` | the TSR key of the triangle |
| `label1`, `label_name1`, `structure1`, `side1`, `residue1`, `chain1`, `seq1`, `atom1` | vertex 1: the integer label code the key was built from, the label scheme's name for that code, the structure id, the side (`primary` for a one-sided method, `A` or `B` in a cross method), the residue name actually present, the chain, the residue number with any insertion code, and the atom name |
| the same eight with `2` and `3` | vertices 2 and 3, in the canonical vertex order the key uses |
| `theta_class`, `theta`, `theta_raw` | the angle class (1–29), the folded angle in degrees (0–90) and the angle before folding, which decides the mirror sign |
| `distance_class`, `max_distance` | the distance class and the longest edge, in ångström |
| `x1` … `z3` | the three vertices' coordinates |
| `category` | the frequency bucket the triangle is counted in: `regular`, `mirror` with `--mirror`, or the SSE category with `--sse` |
| `sse_state1` … `sse_group3`, `sse_category` | with `--sse`, each vertex's secondary-structure state and group and the triangle's category; `-` in all seven when SSE was not requested |

Floating-point values are written so that they read back exactly. Under
`--aa-grouping`, `label_name` is the **group's representative** residue name
(an ILE vertex grouped with LEU shows `LEU`), while `residue` is the residue
actually present (`ILE`). In an atomic method `label_name` is the element
symbol the label stands for (`C`, `N`, `P`, `Cl`, `Zn`, …), while `atom` is
the atom's own name (`PA`, `CL1`, …).

---

## 16. Loading previous results

`biotsr describe` reads an output folder (one line per calculation: prefix,
method, triangles, complete or incomplete), a prefix or its `PREFIX.run.json`
(how the calculation was produced, including the source rows of a batch
request), or a `.meta.json` sidecar (what the numbers mean):

```bash
biotsr describe results/
biotsr describe results/1UBQ_A_protein
biotsr describe results/1UBQ_A_protein.meta.json
```

```text
1UBQ_A_protein protein 70300 triangles complete
command biotsr run
prefix 1UBQ_A_protein
complete True
key space da39d73510c54a4c07529e17f34449789584fad4ea79b020b1badd3b5550d30b
```

To get the result back as an object, load its keys file (or its sidecar). Both
identities and their fingerprints are restored:

```python
loaded = biotsr.load_result("results/1UBQ_A_protein.keys")
assert loaded.frequencies == result.frequencies
assert loaded.key_space.fingerprint == result.key_space.fingerprint
assert loaded.sampling.fingerprint == result.sampling.fingerprint
assert loaded.triplets is None        # per-triangle records are not restored
```

Every output this version writes is native and can be loaded this way.
Results in the per-sample layout of early pre-release builds are refused (§15).

---

## 17. Analysis

The analysis layer works on results from any method, provided the results
share a **key space**.

**Common keys**, from the command line, over saved results:

```bash
biotsr analyze common-keys --from-results batch/ --output analysis_out/
```

```text
loaded 3 saved result(s) from 3 calculation(s) in batch/
common keys 9523 over 3 samples
```

The report is written to `analysis_out/analysis/common_keys/` as
`all_samples.csv` (per sample: `sample_id`, `method`, `total_keys`,
`common_keys`, `sum_common_keys_freq`, and the key-space, processing and
sampling fingerprints) and `all_samples.json` (the keys, the same rows, each
sample's method, selection and processing settings, and a `comparison`
summary).

* **With `--from-results`** nothing is computed. Completion records
  (`PREFIX.run.json`) are found recursively under the folder (hidden
  directories, such as staging ones, and `analysis/` folders are never read),
  and each calculation is reported under its prefix. An incomplete record, a
  file that does not match the size recorded for it, and a `.keys` file no
  complete record lists are refused, as are results in the per-sample layout
  of early pre-release builds. Flags that would compute — `--method`, `--chain`, modifiers,
  `--pdb-dir`, `--column`/`--keep-column`/`--ignore-column`, the sheet options,
  download options, `--jobs`, `--if-exists`, `--on-error` — are refused.
* **With `--manifest samples.csv --pdb-dir pdb/`** the samples are computed
  first, as `biotsr batch` would, writing **keys only** (no triplet files).
* **`--group-by label`** reports one intersection per label: the manifest's
  `label` column (`aa_group` in a TSR sheet), or with `--from-results` the label
  recorded in each completion record (CLI guide, §8).

**The comparison rule.** Every analysis — common keys, occurrence tables,
feature matrices, similarity and clustering — applies one rule, to computed
and loaded results alike:

1. **Key definitions must match, always.** The same integer must mean the same
   triangle: label and point scheme, angle and distance bins, encoding,
   ordering and mirror convention. This is decided from the key space, never
   from method names, so different methods that share one compare (every atomic
   method labels atoms by element in 29 × 18 bins), while the same method with
   `--mirror` or `--aa-grouping` does not. Nothing relaxes this.
2. **Processing settings must match, unless you acknowledge the difference.**
   Candidate strategy, completeness and composition rules, size filter,
   degenerate and sequence-label policies, contact cutoff, filter, rule and
   projection, pool identity and order, the alternate-location rule, and the
   modified-residue rule (`modified_residue_policy`: which modified amino
   acids written as `HETATM` were read). A result written before the
   modified-residue rule was recorded shows it as `not recorded`, never
   assumed to be the current rule, so comparing it with a current result
   needs `--allow-mixed-processing`; its key definitions are unchanged.
3. **Partial samples must be acknowledged.** A sample whose selection did not
   get everything it asked for (§12.4) is refused, naming each partial sample
   and what it lacks, unless you pass `--allow-partial`
   (`allow_partial=True` in Python). This is separate from
   `--allow-mixed-processing`: neither implies the other, and neither relaxes
   step 1.
4. **What was selected never blocks a comparison**: structures, chains,
   amino-acid types (LYS-only against TYR-only), ligands, carbohydrates,
   nucleotides, residue numbers, atom subsets, one copy against all copies, and
   the method name. Each sample's selection is recorded in the report. A
   complete LYS-only sample and a complete TYR-only sample compare without any
   option.
5. **Samples need not share any key**: an empty intersection is a valid
   result.

The full rule is in the [analysis methods record](methods/ANALYSIS.md), §1.1.

A set whose processing differs — one sample with a size filter and one
without, say — is refused, because its frequencies were drawn differently,
unless you pass `--allow-mixed-processing`:

<!-- file: mixed_sampling.csv -->
```csv
sample_id,pdb_id,method,chain,size_filter
plain,1UBQ,protein,A,
filtered,1BRS,protein,A,25
```

```bash
biotsr analyze common-keys --manifest mixed_sampling.csv --pdb-dir pdb/ --output mixed_ck/ --allow-mixed-processing
```

```text
common keys 17445 over 2 samples
biotsr: warning: the samples were not all processed with the same settings (2 processing policies; differing in: size_filter).
```

Without `--allow-mixed-processing` the same command stops with an error naming
the difference (`size_filter (None vs 25.0)`). With it, the report is written
and records the warning and, per sample, the fingerprints, the processing
settings, the selection and the warnings of the sample's own result. In
Python, every analysis function takes `allow_mixed_processing=True`. (In
early pre-release builds this option was `--allow-mixed-sampling` /
`allow_mixed_sampling=`, and selections needed it too; the old flag now stops
with an error naming the new one.)

A set with a **partial** sample is refused in the same way until you
acknowledge it. Chain B of `SYNPHOS.pdb` has no SER, so a TYR, LYS and SER
selection there is partial, while the same selection on chain A is complete.
With both computed into `out/` (a batch of the two rows, which exits `4`):

<!-- not-executed -->
```bash
biotsr analyze common-keys --from-results out/ --output partial_ck/ --allow-partial
```

Without `--allow-partial` the command stops (exit `2`) with
`1 of the 2 samples are partial -- their selections did not get everything
they asked for: chainB (not found: SER). …`. With it, nothing is dropped: the
report's `comparison` records `partial_allowed` and `partial_samples` (each
partial sample's omissions), each sample's entry records `partial` and
`omissions`, a warning names them, and the command exits `4`. In Python:

```python
from biotsr.analysis import AnalysisSample, PartialSamples, common_key_intersection

chains = [AnalysisSample.from_result(
              biotsr.analyze("SYNPHOS.pdb", method="aa_ind_mix", chain=chain,
                             residue_types="TYR;LYS;SER", keep_triplets=False),
              "chain" + chain)
          for chain in ("A", "B")]
assert [s.is_partial for s in chains] == [False, True]
try:
    common_key_intersection(chains)
except PartialSamples as refused:
    assert list(refused.partial) == ["chainB"]
else:
    raise AssertionError("a partial sample must be acknowledged")

allowed = common_key_intersection(chains, allow_partial=True)
assert allowed.n_common_keys == 258
assert allowed.comparison["partial_samples"]["chainB"]["primary"] == {"not_found": ["SER"]}
assert "chainB (not found: SER)" in allowed.warnings[0]
```

**Feature matrices, similarity and clustering** are in Python:

```python
from biotsr.analysis import (AnalysisSample, clustermap_spec,
                             common_key_intersection, cut_hierarchy,
                             distance_matrix, feature_matrix, hierarchy)

names = {"ubiquitin": "batch/ubiquitin.keys",
         "barnase": "batch/barnase.keys",
         "barstar": "batch/barstar.keys"}
samples = [AnalysisSample.from_result(biotsr.load_result(path), name, label=name)
           for name, path in names.items()]

assert common_key_intersection(samples).n_common_keys == 9523
fm = feature_matrix(samples)          # one column per key, in numeric order
dm = distance_matrix(fm)              # generalized Jaccard distance by default
assert round(dm.value(1, 2), 4) == 0.8674     # barnase vs barstar
tree = hierarchy(dm)                  # average linkage and leaf order; no cluster count
spec = clustermap_spec(dm, [s.label for s in samples])   # dendrograms + ordered heat map
assert spec.assignment is None        # no partition unless you ask for one
two = cut_hierarchy(tree, n_clusters=2)   # an optional flat cut, by an explicit rule
```

**The whole workflow, executed:** [tutorial 08](tutorials/08_similarity_clustering.md)
loads a batch's saved results, matches them to class labels, builds the feature
matrix, computes generalized and binary Jaccard similarities and distances,
exports them, builds the tree, draws the clustered heat map and records
provenance. In brief:

* **Python only.** The command line computes results and common keys;
  feature matrices, similarity, clustering and the figure are in
  `biotsr.analysis`.
* **Generalized versus binary Jaccard.** Generalized Jaccard weighs each key by
  its frequency (`Σ min / Σ max`); binary Jaccard compares only which keys are
  present (`|A ∩ B| / |A ∪ B|`). Three more measures are available (`wu`,
  `sarika`, `cosine`).
* **Similarity versus distance.** Every measure is a similarity in [0, 1];
  `.to_distance()` gives `D = 1 − J`, and nothing converts implicitly. Values are
  fractions unless you ask for percent (`.scaled(100)`), which every label and
  export then states.
* **Loading saved results and labels.** Read each request's prefix from
  `batch_summary.csv` (its `sample_id` column; keep the `ok` rows, and a
  `partial` row only once its `message` says the result is what you want, and
  then pass `allow_partial=True`), load
  `PREFIX.keys` with `biotsr.load_result`, and take class labels from your own
  table matched by sample id — never from a file or sample name. Rows the
  completeness rule excluded appear there as failed rows naming the rule.
* **The tree needs no cluster count.** `hierarchy` (average linkage) and
  `clustermap_spec(matrix)` give the dendrograms and the ordered heat map, with
  an "observed class" stripe if you pass labels. A flat partition is optional
  and needs an explicit rule: `cut_hierarchy(tree, n_clusters=t)` (`maxclust`)
  or `cut_hierarchy(tree, height=h)`, or the same `n_clusters=` / `height=` on
  `clustermap_spec`, which adds a separate "cluster" stripe and, when every
  sample is labelled, the adjusted Rand index. Never take the count from the
  known labels: an uncut tree does not say how many groups there are.
* **Files to inspect.** The square matrices (`<measure>_<similarity|distance>.csv`,
  sample ids on both axes), the tree (linkage and leaf order), the cluster
  assignments when you cut it (`sample_id, cluster, leaf_position,
  class_label`), and the provenance JSON.
* **The figure** needs matplotlib: `pip install "biotsr[plot]"`. Without it,
  every number is still computed and drawing raises an error naming that
  command.
* **Size.** Similarity and clustering are quadratic in the number of samples and
  run in one process; `--jobs` speeds up the structure computations of a
  batch, not the analysis. Measured timings are in the
  [analysis methods record](methods/ANALYSIS.md), §5.8.

The [Python API guide](API_GUIDE.md), §6–7, lists the functions; the
definitions are in the [analysis methods record](methods/ANALYSIS.md).

**Comparability.** Integers from different key spaces mean different things,
so combining them is refused, not silently computed as zero overlap:

```python
mirrored = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", mirror=True,
                          keep_triplets=False)
try:
    feature_matrix(samples + [AnalysisSample.from_result(mirrored, "mirrored")])
except biotsr.IncompatibleKeySpace as error:
    assert "mirror_convention" in str(error)
else:
    raise AssertionError("mixing key spaces must be refused")
```

The same guard refuses a Cα result next to an atomic one. A processing
difference only, such as a different contact cutoff or size filter, is refused
by every analysis function unless `allow_mixed_processing=True`, as above.
What was selected never is:

```python
lys = biotsr.analyze("1BRS.pdb", method="aa_ind", chain="A",
                     residue_types="LYS", keep_triplets=False)
tyr = biotsr.analyze("3PTB.pdb", method="aa_ind", chain="A",
                     residue_types="TYR", keep_triplets=False)
pair = common_key_intersection([AnalysisSample.from_result(lys, "lys"),
                                AnalysisSample.from_result(tyr, "tyr")])
assert pair.n_samples == 2            # LYS-only against TYR-only: no option needed
```

---

## 18. Common workflows

**I have one PDB and one chain.**
`biotsr run --method protein --pdb my.pdb --chain A --output results/` (§3).

**I have a dataset sheet with hundreds of PDB ids.** If it has one row per
calculation, use it as a manifest. Check the plan first, then run it,
downloading what is missing (§2.4, §13):

<!-- not-executed -->
```bash
biotsr batch --manifest my_sheet.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --download-missing --output results/ --no-triplets --dry-run
biotsr batch --manifest my_sheet.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --download-missing --output results/ --no-triplets --jobs auto --if-exists skip
```

Add `--chain A` if the sheet has no chain column. If it is a wide,
one-row-per-protein TSR sheet, read it with `--sheet tsr` and choose the
methods with `--method` (§13.9), previewing the requests with `--dry-run`
first:

<!-- not-executed -->
```bash
biotsr batch --manifest my_sheet.csv --sheet tsr --method carbohydrate --method aa_ind --inapplicable skip --pdb-dir pdb/ --download-missing --output results/ --dry-run
```

On a cluster, use
`python -m biotsr` and see the SLURM examples in the
[CLI guide](CLI_GUIDE.md), §12.

**I need SSE-TSR.** Add `--sse` to `protein` or `aa_selection`. The PDB file
must carry `HELIX`/`SHEET` records (§6).

**I need protein–ligand cross keys.**
`biotsr run --method drug_protein_cross --pdb my.pdb --chain L --chain-2 P --entity-name LIG --entity-id 1 --output results/`,
with the ligand's chain as `--chain` and the protein's as `--chain-2` (§10).
Leave out `--entity-id` to use every copy of the ligand on chain L (§12.3).

**I need to compare several generated results.**
`biotsr analyze common-keys --from-results results/ --output analysis_out/`,
or load them in Python and use `biotsr.analysis` (§17,
[tutorial 08](tutorials/08_similarity_clustering.md)). All results must share a
key space (same label scheme and bins, same modifiers) and their processing
settings, unless you pass `--allow-mixed-processing`
(`allow_mixed_processing=True` in Python), and a partial sample needs
`--allow-partial` (`allow_partial=True`); what was selected may differ freely
(§17).

**I need the method behind a published TSR result.** The
[reproducibility note](REPRODUCIBILITY_NOTE.md) says which command computes
each published method under the standard.

---

## 19. Troubleshooting

Ordinary mistakes print one line starting `biotsr: error:` and exit with status
2; `--debug` adds a traceback when you suspect a defect.

| symptom | cause and fix |
|---|---|
| `chain 'Z' is not in 1UBQ, so the selection produced no atoms; chains in this structure: A` | the chain is not in the file; use one of the chains listed |
| `nothing to compute: … the method's completeness rule excluded the selection (1_DC: DC: expected exactly 19 atoms, got 16)` | the target exists but lacks atoms the method requires (here the 5′-terminal phosphate); leave it out, or use the `*_atom_selection` method |
| `nothing to compute: … no atom pair across the two sides lies within 3.5 A` | both sides selected atoms, but none is in contact; check the chains, or widen `--contact-cutoff` |
| `no chain given`, or a batch status `chain_not_specified` | every run needs its chain (and a cross run `--chain-2`), from the row, `--column chain=HEADER`, `--map`, or `--chain` for rows that leave it blank; none is guessed. A batch skips such a request, lists it in `batch_summary.csv` and exits 4 |
| `… includes residue(s) named UNK …` | no calculation is done for UNK; narrow the selection to leave it out (`aa_selection` residues, named residues, residue types), or give the residue its real name |
| `… blank element field …`, `… repeated records …`, `… conflicting records …`, `… exactly the same coordinates …`, `… residue number BioTSR cannot read …` | the selection includes records BioTSR cannot read reliably; correct the structure file (element symbols in columns 77–78, each atom written once, residue numbers below 10000), or leave those residues or atoms out |
| `… cannot be labelled …` | the element labels have no label for that atom's element field (blank, `X`, or not an element); correct columns 77–78, or leave the atom out |
| `… name the same chain` / `… select N of the same atom(s) …` | a cross needs two disjoint sides; name two different protein chains, or change one side's selection |
| `SSE annotation of chain … is refused: …` | a `HELIX`/`SHEET` record of that chain names a residue that is not in it (or cannot be read); correct the record, or run without `--sse` |
| `method 'drug' needs --entity-name` | the method needs a target; `biotsr methods NAME` lists its target options |
| `method 'drug' does not support 'aa-grouping'` followed by `it supports: --mirror, --size-filter` | that modifier is not defined for that method; check the [capability matrix](CAPABILITY_MATRIX.md) |
| `method 'protein' does not use --residues, so it would be silently ignored; it is refused instead` | the method does not read that argument; the message lists what it does take (§12.2) |
| `method 'drug_protein_cross' does not use --residue-types-2 …` | narrow the protein side of a drug, carbohydrate or nucleotide cross with the unsuffixed argument (§10) |
| `cross method 'drug_protein_cross' needs --chain-2` | every cross method names both chains (§10) |
| `entity BEN 7 is not in chain A of 3PTB: copies of BEN there have ids 1` | the id does not exist; use one listed, or leave out `--entity-id` for every copy (§12.3) |
| `--entity-id 'second' is not a PDB residue identifier` | give the residue number as written in the file, with an insertion code if it has one (`404`, `52A`, `-5`) (§12.3) |
| `residue token 'LEU' is neither NAME_ID (e.g. LYS_6) nor NAMEID (e.g. GLU18)` | name each residue with its identifier (§5) |
| `method 'amino_acid_atomic' was given --residue without --residue-name` | give both, `--residue-name` alone for each residue of that name, or neither for every residue (§5) |
| `none of the requested residues …`, `no residue of type TRP …` | the selection matches nothing in that chain; check the chain and the spelling (§12.4) |
| `biotsr: warning: requested but not found: …` | some selection items are not in the chain; the result was computed from the rest (§12.4) |
| `biotsr: warning: no triangles: …`, or a batch row with status `no_triangles` | the selection has fewer than three points (§12.4) |
| `not supported: …` or `not represented: …`, or a batch row with status `partial` | the chain holds a modified residue the method cannot read, or something you named is missing or excluded; the result was computed without it, and the command exits `4` (§5, §12.4) |
| `the residue range(s) 9-11 of chain A span residue(s) 10A, 10B, which have an insertion code` | list those residues to include them (`9-11;10A;10B`), or list the residues one by one to leave them out (§4) |
| `residue range '27A-29' is not supported`, `residue range '12-10' is descending` | a range is two plain residue numbers, ascending; list an insertion-code residue on its own (§4) |
| `the sheet populates columns that could narrow a cross calculation, and BioTSR does not choose among them` | choose each named side with `--cross-side METHOD.SIDE=COLUMN` or `=whole` (§13.9) |
| `… are requested by rows with different labels` | give those rows distinct `sample_id` values, or correct the label (§13.5) |
| `unrecognized arguments: --profile …` or `--format …` | there are no profile or format options; every method runs its standard profile and writes native output (§12.1) |
| `manifest … has a 'profile' column, which this version of BioTSR does not read` | drop the column, or pass `--ignore-column profile` |
| `… value '…' is in column 4, which has no header` | a list in one cell was separated with commas; use `;` or quote the cell (§13.1) |
| `manifest … has more than one column named '…'` | rename or remove the duplicate header |
| `sample id '…' cannot be used as an output prefix` | use a plain file-name stem (§13.5) |
| `results do not share a key space; their integer keys are not comparable. Differences: …` | you combined results with different methods or modifiers; compare like with like |
| `… share key definitions but were processed with different settings: …` | the samples differ in a processing setting (size filter, contact cutoff, completeness rule, …); analyse them separately, or pass `--allow-mixed-processing` (§17) |
| `… processed with different settings: modified_residue_policy ('hetatm_amino_acids:PTR,SEP,TPO' vs not recorded)` | one result was written before the modified-residue rule was recorded and may have left out TPO, SEP and PTR; recompute it with `--if-exists overwrite`, or pass `--allow-mixed-processing` to compare it as it is (§13.8, §17) |
| `… its keys were computed in a different key space …`, or an analysis refusing a result's key space | an atomic result of a pre-release build that labelled atoms by name (label scheme `atom_lexicon_standard_112`), or a Cα result of a pre-release build with residue labels at radix 20; compute it again (`--if-exists overwrite`) (§3, §5). No option compares them with current results |
| `N of the M samples are partial -- their selections did not get everything they asked for: …` | some samples lack something they asked for (the message names each and what it lacks); fix the requests or the structures, or pass `--allow-partial` to analyse them as computed (§17) |
| `--allow-mixed-sampling was renamed --allow-mixed-processing` | use the new name; selections never need it (§17) |
| `… in the per-sample directory layout with the bin counts in the file name …` | that output's layout is not read by this version; recompute it into a new `--output` (§15) |
| `… is not part of a completed calculation: no completion record lists it`, or `the calculation … is not complete` | an interrupted run left it; recompute it (`--if-exists skip` or `overwrite`), or move it out of the folder (§13.8) |
| `--from-results reads results that were already computed, so --method would have no effect` | drop the computing flags with `--from-results`, or use `--manifest` (§17) |
| `--download-workers must be from 1 to 16`, `--timeout must be greater than 0 seconds` | invalid download settings are refused, never replaced (§2.4) |
| `cannot read structure missing.pdb: No such file or directory` | the `--pdb` path is wrong |
| `no structure file for '9XYZ' in pdb/` | the manifest's `pdb_id` has no `9XYZ.pdb` or `9XYZ.ent` in `--pdb-dir`; add `--download-missing` or run `biotsr fetch` |
| `the manifest identifies structures by pdb_id but no directory was given` | pass `--pdb-dir`, or use a `pdb_path` column |
| `manifest … has unrecognised column(s): …` | a column name is misspelled, or it is one of your own columns: read it with `--column FIELD=HEADER`, keep it with `--keep-column`, or drop it with `--ignore-column`; a wide sheet needs `--sheet` (§13.9) |
| `… names a method in N row(s) … and --method also chooses methods` | choose the methods in one place: drop `--method`, or pass `--ignore-column method` (§13.1) |
| `the row has no carb_name / carb_id for carbohydrate` | that row does not describe the method's calculation; pass `--inapplicable skip` to record such rows as `inapplicable` (§13.9) |
| `… would each be written by more than one different request` | give those rows distinct `sample_id` values; prefixes differing only in letter case collide too (§13.5) |
| `… exists but cannot be reused: its completion record does not state the processing rules it was computed under …` | the output was written before the processing rules were recorded and may have left out TPO, SEP and PTR; recompute it with `--if-exists overwrite`, or into a new `--output` (§13.8) |
| `… exists but cannot be reused: it was computed under different processing rules (…)` | recompute it with `--if-exists overwrite`, or into a new `--output` (§13.8) |
| `… exists but cannot be reused: …` | `--if-exists skip` found an output of a different request; the message names the difference |
| `output … already exists in …` | a run never overwrites; pass `--if-exists skip` or `--if-exists overwrite`, or choose another `--output` |
| `… already exists in … and is not part of an earlier version of this calculation` | a file BioTSR did not write as this calculation has that name; it is never replaced: remove it, or give the request its own `sample_id` |
| `another batch is writing into …` | a batch is already running into that `--output`; wait for it, or use another output directory |
| a batch exits with status 3 | some requests failed, or a worker process died and requests were left `failed` or `not_started`; each failure is listed in `batch_failures.csv` with its manifest line(s); re-run with `--if-exists skip` |
| a batch exits with status 130 | it was interrupted; re-run with `--if-exists skip` |
| `biotsr run`, `batch` or `analyze` exits with status 4 | nothing failed, but at least one result is partial (or a resumed run kept a partial one, or an analysis used partial samples under `--allow-partial`); the results are written; `batch_summary.csv` and the warnings say what each lacks (§12.4, §13.7) |
| `biotsr fetch` exits with status 3 | an id was not downloaded, refreshed or cached (`refresh_failed` included); see `fetch_report.csv` |
| `analysis/analysis/common_keys` when `--output analysis/` | `analyze common-keys` always writes into `OUTPUT/analysis/common_keys/`; choose an output name such as `analysis_out/` |

---

## 20. Where to go next

| | |
|---|---|
| [Tutorials](tutorials/README.md) | eight short, runnable walk-throughs |
| [Command-line guide](CLI_GUIDE.md) | every command and option, the full manifest format, exit codes, HPC use |
| [Python API guide](API_GUIDE.md) | run, inspect, save, load, analyse from Python |
| [Capability matrix](CAPABILITY_MATRIX.md) | which method supports what, generated from the registry |
| [Reproducibility note](REPRODUCIBILITY_NOTE.md) | which command computes which published method |
| [Methods](methods/README.md) | the scientific definitions of cross TSR, SSE-TSR and analysis |
| [Known limitations](KNOWN_LIMITATIONS.md) | what 1.0 does not do |
| [Architecture](ARCHITECTURE.md) | how the package is built |
