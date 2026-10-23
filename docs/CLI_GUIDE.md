# The `biotsr` command line

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology.

A guide for someone who wants TSR keys out of a PDB file and does not intend to
read the package's source.  Every command here is parsed by the validation
suite, and the workflows are executed by the test suite; the examples are not
illustrative, they are tested.

---

## 1. Install

```bash
pip install .                                       # from a clone of https://github.com/TSRMethod/BioTSR
pip install dist/biotsr-1.0.0rc7-py3-none-any.whl   # or from a built wheel
```

`pip install biotsr` will be the installation once 1.0 is published on PyPI;
it is not yet. That is the whole installation, parallel execution included.  The package has
one third-party dependency, joblib, which runs batch rows in worker processes
when `--jobs` asks for more than one (ADR-086); the calculation itself uses only
the standard library — no pandas, no numpy, no scipy, no BioPython — so it
installs in seconds and runs on Linux and macOS wherever Python 3.9+ runs,
including a compute node with no scientific stack.  Windows is not supported.
(On a node with no network,
install from a wheelhouse: `pip wheel . -w wheels/` in a clone on a machine
that has a network, then `pip install --no-index --find-links wheels/ biotsr`.)

Two equivalent entry points:

```bash
biotsr --help           # the console script
python -m biotsr --help # the same program, for an environment whose bin/ is not on PATH
```

Use `python -m biotsr` on a cluster if `biotsr` is not found: a console script only
works when the environment's `bin` directory is on `PATH`, which under a job
scheduler it often is not.

For development, install in place:

```bash
pip install -e ".[dev]"   # from the repository root
```

---

## 2. First output in one command

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

```
method     protein
profile    standard_protein_calpha_v1
triangles  70300
keys       45349 distinct
output     results/1UBQ_A_protein (4 files, completion record 1UBQ_A_protein.run.json)
```

The summary goes to **stderr**; the path of the completion record goes to
**stdout**, so this works:

```bash
RECORD=$(biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/ --quiet)
```

`--output` is one shared folder.  Every calculation writes its files directly
into it, under one **prefix** (here `1UBQ_A_protein`):

| file | what it is |
|---|---|
| `1UBQ_A_protein.keys` | `key`, `category`, `frequency` — one row per (key, category), sorted numerically by key |
| `1UBQ_A_protein.meta.json` | the self-describing sidecar: key space (with both bin schemes), sampling profile, both fingerprints, selection record, source SHA-256, the candidate ledger |
| `1UBQ_A_protein.triplets` | one row per triangle, with full per-vertex identity and geometry (every column is explained in the [user guide](USER_GUIDE.md), §15); left out with `--no-triplets` |
| `1UBQ_A_protein.run.json` | the **completion record**: the command, the arguments, the profile, the normalised request, every file with its size, and any warnings |

`.keys` and `.triplets` are tab-separated UTF-8 text with one header line; the
suffix names the content, not a different format.  The file names carry no bin
counts: the sidecar records the full bin definitions.  Nothing else in the
folder is touched, so a second calculation simply adds its own files next to
these.

Outputs written by 1.0.0rc3 or earlier — one directory per sample, with
`run.json` inside and names such as `1UBQ_A.keys_29_35` — are not read by this
version: `describe`, `analyze --from-results` and `biotsr.load_result` refuse
them with a message saying so.  Recompute them.

---

## 2a. What you supply, and what BioTSR supplies

The line between the two is not a convenience: it is what makes a BioTSR result
reproducible by someone who has only the package and your inputs.

### For a single run you supply

| what | how |
|---|---|
| a **structure file** | `--pdb path/to/1UBQ.pdb` |
| a **method** | `--method protein` (`biotsr methods` lists all twenty) |
| the **chain / entity / selection** the method requires | `--chain A`, plus `--chain-2`, `--entity-name`, `--residues`, `--atoms` … as the method demands |
| optional **scientifically supported modifiers** | `--mirror`, `--aa-grouping`, `--size-filter`, `--sse` |
| cross methods' **sampling parameter** | `--contact-cutoff` (which atoms enter the pool; not a modifier, see §5) |

### For a batch run you supply

| what | how |
|---|---|
| a **dataset manifest** or **wide sheet** | `--manifest sample_details.csv` (add `--sheet tsr` or `--sheet custom` for a wide sheet, §7) |
| the **methods**, when the file does not name them | `--method NAME`, repeatable |
| for a cross method from a wide sheet, **each side's selection** when the row could narrow it | `--cross-side METHOD.SIDE=COLUMN` or `=whole` (§7) |
| the **structure files** it names | `pdb_path` per row, or `pdb_id` plus `--pdb-dir` |

`sample_details.csv` is only an example filename — the manifest is identified by
`--manifest`, and nothing in the package depends on what it is called.

### You do **not** supply — ever — any of this

These are part of the scientific definition of the software. They are versioned
with BioTSR, they ship inside the installed package, and a normal workflow has
no argument for them:

* the **element labels** (`element_labels_v1.csv`, element → label; every
  atomic method labels an atom by its PDB element field, never by its name) —
  loaded through `importlib.resources`, never from the working directory;
* the **key formula** and the canonical-ordering rules;
* the **theta-bin** and **distance-bin** definitions;
* the **standard profile** of every method;
* the **amino-acid label tables**, SSE categories and every other standard
  scientific lookup table.

If a TSR implementation requires its user to bring a label table, the user is being
made responsible for a scientific constant. BioTSR does not do that: three
historical programs did, and had to be `cd`-ed into before they would run at
all.

## 3. Finding out what the package can do

```bash
biotsr methods
```

```
METHOD                         FAMILY   MODIFIERS
aa_ind                         atomic   mirror, size_filter
aa_ind_mix                     atomic   mirror, size_filter
aa_ind_specific                atomic   mirror, size_filter
aa_selection                   calpha   aa_grouping, mirror, size_filter, sse
amino_acid_atom_selection      atomic   mirror, size_filter
amino_acid_atomic              atomic   mirror, size_filter
carbohydrate                   atomic   mirror, size_filter
carbohydrate_atom_selection    atomic   mirror, size_filter
carbohydrate_chain             atomic   mirror, size_filter
carbohydrate_protein_cross     cross    size_filter
drug                           atomic   mirror, size_filter
drug_atom_selection            atomic   mirror, size_filter
drug_protein_cross             cross    size_filter
nucleotide                     atomic   mirror, size_filter
nucleotide_atom_selection      atomic   mirror, size_filter
nucleotide_chain               atomic   mirror, size_filter
nucleotide_protein_cross       cross    size_filter
protein                        calpha   aa_grouping, mirror, size_filter, sse
protein_protein_calpha_cross   cross    size_filter
protein_protein_cross          cross    size_filter
```

One method in full, including which target arguments it takes (its `targets`
line) and the standard profile it runs:

```bash
biotsr methods drug
```

```
drug
  family     atomic
  scope      pooled
  ...
  modifiers  --mirror, --size-filter
  targets    --entity-id, --entity-name
  profile    standard_drug_v1
```

And for a script — each record has the keys `name`, `family`, `scope`,
`description`, `modifiers`, `target_arguments`, `profile` and
`profile_description`:

```bash
biotsr methods --json | jq -r '.[] | select(.family=="cross") | .name'
```

This listing is **generated from the method registry**.  It is not a second
hand-maintained list, so it cannot say the package supports something it does
not, or omit something it does.

> Candidate-strategy extensions (for example CA-kNN) are **not part of this
> release** and are absent from the CLI on purpose.  There is no flag that turns
> them on.

---

## 4. One structure

### Protein Cα, plain

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

### The three Cα modifiers, composed

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/ \
    --aa-grouping --mirror --size-filter 25.0
```

* `--aa-grouping` — the grouped residue labels instead of the individual ones:
  the 13 groups of the 20 regular amino acids, MSE with MET, and each other
  supported modified residue in a group of its own (23 codes), at radix 30
  like the individual labels and in a different key space.
* `--mirror` — negate a key when its pre-fold angle exceeds 90°.  Off unless
  asked for; with it off **no key is ever negative**.
* `--size-filter A` — keep only triangles whose longest edge is strictly below
  `A` ångström.

All eight combinations work.  There is no per-combination method name and there
is not meant to be.

### Secondary structure

```bash
biotsr run --method protein --pdb 1TUP.pdb --chain A --output results/ --sse
```

```
method     protein
profile    standard_protein_calpha_v1
sse        on, strand identity sheet_id_and_serial
triangles  1235780
keys       431473 distinct
```

`--sse` is an **annotation, not a modifier**.  It does not change a single key:
the keys, their signs, the candidate ledger and both identity fingerprints are
identical to the same run without it.  What it changes is which of the 18
secondary-structure categories each key's frequency is counted in.  It composes
freely with all three modifiers, and it is offered only by `protein` and
`aa_selection` — the 18 categories are defined over protein residues, so asking
for it anywhere else is refused by name.

### A ligand

```bash
biotsr run --method drug --pdb 3PTB.pdb --chain A --output results/ \
    --entity-name BEN --entity-id 1
```

`--entity-id 1` names that one copy (prefix `3PTB_A_drug_BEN_1`).  Without
`--entity-id`, every copy of the name on the chain is used together; see
*Entities* below.

### Every residue of a chain, atomically

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --output results/
```

Omitting `--residue` and `--residue-name` runs one calculation per residue of
the chain — 75 of 1UBQ's 76, since the terminal glycine carries an extra atom
and the completeness policy declines it (the reason is recorded, not inferred
from an absence).  Each admitted residue gets its own prefix,
`1UBQ_A_amino_acid_atomic_perresidue_LEU_15` and so on, and the sweep has one
completion record, `1UBQ_A_amino_acid_atomic_perresidue.run.json`.
`--residue-name` alone sweeps the residues of that name (prefixes
`1UBQ_A_amino_acid_atomic_perresidue-LEU_LEU_15`, …); both together run one
residue:

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --output results/ \
    --residue-name LEU
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --output results/ \
    --residue 15 --residue-name LEU
```

`--residue` is the residue's identifier in the file: its number, with an
insertion code if it has one (`52A`).  `--residue` without `--residue-name` is
refused.

### Residue and atom selections

```bash
# Cα of listed residues, with ranges
biotsr run --method aa_selection --pdb 1UBQ.pdb --chain A --output results/ \
    --residues '10-12;48;70-71'

# named atoms of named residues, pooled
biotsr run --method amino_acid_atom_selection --pdb 1UBQ.pdb --chain A \
    --output results/ --atoms 'LYS_11_NZ;LEU_15_CB;ILE_13_CG1'

# every heavy atom of every residue of a type
biotsr run --method aa_ind_mix --pdb 1UBQ.pdb --chain A --output results/ \
    --residue-types 'LEU;ILE'

# every heavy atom of named residues: NAME_ID or NAMEID
biotsr run --method aa_ind_specific --pdb 1UBQ.pdb --chain A --output results/ \
    --residues 'LEU_15;LEU50'
```

Quote selections: `;` is a shell metacharacter.  A residue is named by its
three-letter name and its identifier in the file: `LEU_15` and `LEU15` are the
same residue (a name containing a digit needs the `NAME_ID` spelling), and an
atom by `RES_ID_ATOM` (`LYS_11_NZ`, or `LYS_52A_NZ` with an insertion code).
A token that is neither — `LEU` alone, say — is refused, not skipped.

### Residue identifiers and ranges

A residue identifier is the residue's PDB number with an optional one-letter
insertion code — `48`, `52A`, `-5` — wherever one is written: `--residue`,
`--entity-id`, `NAME_ID` references (`LYS_52A`, `LYS52A`), atom tokens
(`LYS_52A_NZ`) and the items of an `aa_selection` list.  Leading zeros are
dropped (`048` is `48`) and the insertion code is matched exactly (`52` never
matches `52A`).  An `aa_selection` item may also be a range of plain numbers,
ascending, both ends included: `10-12`, `-5--3`.  Refused, naming the item:

* a range with an insertion code at either end (`27A-29`);
* a descending range (`12-10`);
* any other malformed item (`9-`, `ten`);
* a range that spans a residue of the structure with an insertion code —
  `9-11` on a chain that holds `10A` — unless that residue is listed too.
  Whether a numeric range should include it has no single obvious answer, so
  you say: list it to include it (`9-11;10A`), or list the residues one by one
  to leave it out (`9;10;11`).

`SYNPHOS.pdb`, a synthetic structure in `tests/data/synthetic/`, has residues
`10A` and `10B` between `10` and `11` on chain A:

```bash
# 9, 10, 10A, 10B and 11
biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --output results/ \
    --residues '9-11;10A;10B'
# 9, 10 and 11 only
biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --output results/ \
    --residues '9;10;11'
```

```
$ biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --residues '9-11' --output results/
biotsr: error: the residue range(s) 9-11 of chain A span residue(s) 10A, 10B, which have an insertion code: a numeric range names plain residue numbers only. List them to include them (e.g. 9-11;10A;10B), or list the residues one by one to leave them out
```

List order never matters: the points follow the structure file.  The command
line, a manifest's `residues` column and the Python API share one parser, so
the same list is accepted or refused the same way everywhere.

### Arguments are checked per method

Every argument you give must be one the method reads.  Anything else is
refused by name rather than silently ignored, and the message lists what the
method does take:

```
$ biotsr run --method protein --pdb 1UBQ.pdb --chain A --residues '1-5' --output results/
biotsr: error: method 'protein' does not use --residues, so it would be silently ignored; it is refused instead. It takes no target argument beyond the chain.
  `biotsr methods protein` describes it
```

### Modified amino acids

RCSB files write modified amino acids as `HETATM` records.  The atomic methods
read exactly three of them — phosphothreonine `TPO`, phosphoserine `SEP` and
phosphotyrosine `PTR`, which have exact heavy-atom counts (11, 10 and 16) and
whose atoms are labelled by element like any other — in the amino-acid atomic methods (`amino_acid_atomic`,
`amino_acid_atom_selection`, `aa_ind`, `aa_ind_mix`, `aa_ind_specific`) and on
the protein side of the atomic cross methods (`drug_protein_cross`,
`carbohydrate_protein_cross`, `nucleotide_protein_cross`,
`protein_protein_cross`).  The alternate-location, chain, identifier and
completeness rules apply to them unchanged, and none is converted to its
unmodified residue.

* **Any other modified residue** stored as `HETATM` — `CME` (6G9M chain A,
  residue 161), `MSE`, … — has no defined atom count and is **not
  supported** by the atomic methods.  It is named in a warning beginning `not supported:` and listed
  under `unsupported` in the side's selection record; it is never dropped
  silently.  A request of unsupported residues only is refused with that
  reason.
* **The Cα methods** (`protein`, `aa_selection`,
  `protein_protein_calpha_cross`) label residues, not atoms, and read the
  eleven modified residues of the residue labels — `TPO` 24, `SEP` 25, `PTR`
  26, `NIY` 27, `M3L` 28, `TRF` 29, `CSO` 30, `LLP` 31, `CME` 32, `CSD` 33,
  and `MSE` with methionine's 13 — at radix 30 (ADR-107).  A `HETATM` Cα is a
  point only when its residue is one of these and the CA atom's element field
  is carbon or blank; nothing else of the residue is required, as for an
  `ATOM` record.  A ligand or an ion with an atom named `CA` never is.  The
  interface rule of `protein_protein_calpha_cross` reads the same residues'
  heavy atoms.  Any other modified residue is named in a warning beginning
  `not represented:`, and so is a supported one without a readable Cα, with
  the reason (`no CA atom`, a CA only as altLoc B, or a `HETATM` CA whose
  element is not carbon).
* A named residue that the completeness rule declines is reported as
  `requested but excluded by the method's completeness rule`, not as not
  found.

`SYNPHOS.pdb` (`tests/data/synthetic/`) holds `HETATM` TPO 3, SEP 5, PTR 7 and
CME 8 on chain A:

```
$ biotsr run --method aa_ind_mix --pdb SYNPHOS.pdb --chain A --residue-types 'TYR;TPO;CME' --output results/
biotsr: warning: not supported: residue type CME (a modified residue stored as HETATM; of the modified amino acids BioTSR reads only TPO, SEP, PTR, which have defined atom counts)

$ biotsr run --method protein --pdb SYNPHOS.pdb --chain A --output results/

$ biotsr run --method aa_ind --pdb SYNPHOS.pdb --chain A --residue-types CME --output results/
biotsr: error: nothing to compute in chain A of SYNPHOS: not supported: residue type CME (...)
```

The first run pools TYR 6 and TPO 3 and names CME; the second is the chain's
Cα, all four modified residues included, with no warning; the third has
nothing it can read.

### Entities

`drug`, `carbohydrate`, `nucleotide`, their `*_atom_selection` variants and the
entity side of the three entity cross methods select molecules by **chain,
name and identifier**:

* `--chain` is always required, and only that chain is searched;
* `--entity-name` with `--entity-id` is **that one copy** (scope
  `one_instance`, prefix `…_NAG_1`);
* `--entity-name` without `--entity-id` is **every copy of that name on the
  chain, in one combined calculation** (scope `all_copies`, prefix
  `…_NAG_allcopies`), even when there is only one copy.  Triangles may span
  copies; that is the intended calculation.

7LF9 chain D carries one glycan, NAG 1 – NAG 2 – BMA 3 – MAN 4 – NAG 5 – GAL 6 –
MAN 7:

```bash
# one copy: NAG 1, 364 triangles
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG \
    --entity-id 1 --output results/
# every NAG on chain D (1, 2 and 5) together: 11480 triangles
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG \
    --output results/
# the whole glycan, every sugar: a different method
biotsr run --method carbohydrate_chain --pdb 7LF9.pdb --chain D --output results/
```

Every copy of one sugar is not the whole glycan: `--entity-name MAN` on chain D
is MAN 4 and MAN 7 only, while `carbohydrate_chain` (and `--entity-names`) pool
every residue of the vocabulary, different sugars included.

The copies used are recorded in the sidecar's selection record
(`provenance.selection_request`: scope and `instances`, here `NAG_1`, `NAG_2`,
`NAG_5`) and as the provenance `target` (`NAG_1` or `NAG_allcopies`).  The
completeness rule applies to each copy; a declined copy is recorded and
warned about.  `--atoms` applies to each selected copy, and an atom a copy lacks
is named per copy (`requested but not found: XX9 in NAG_1, …`).

An identifier is the residue number as written in the file, which may be
negative, with at most one insertion-code letter: `404`, `52A`, `-5`.  Leading
zeros are dropped (`005` is `5`), and `52` never matches `52A`.  Anything else —
a position such as `second`, a malformed id, an identifier whose records appear
in two separate places on the chain — is refused, as is an id or a name that is
not in the chain, naming what is there:

```
$ biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id second --output results/
biotsr: error: --entity-id 'second' is not a PDB residue identifier: write the residue number as it appears in the structure file, optionally followed by one insertion-code letter (e.g. 404, 52A, -5)
```

### Empty and partial selections

* A selection that matches **nothing** is refused, naming the entity or the
  items, and the command exits `2`.
* Items that are **missing** from an otherwise valid selection are named in a
  warning, and the result is computed from the rest:
  `biotsr: warning: requested but not found: 200`.  The warning is also
  recorded in the completion record (each result's `warnings`).  Items the
  structure holds but the method cannot read are named the same way, as `not
  supported: …` or `not represented: …` (*Modified amino acids*, above).
* A selection with **fewer than three points** forms no triangle: the result is
  written with zero triangles and a `no triangles: …` warning.  In a batch,
  that request's status is `no_triangles`; it is not a failure and leaves the
  exit code unchanged — unless the selection also lacked something it named
  (`requested but not found: 200; no triangles: …`), when the result is also
  partial and the command exits 4 (below).
* The completeness rules of the methods that have one (an exact atom count per
  residue type) are unchanged, and what they exclude is reported; a residue
  you named that they decline is `requested but excluded by the method's
  completeness rule`.  `nucleotide_chain` has none: each nucleotide takes part
  with the atoms it has, and one whose atoms differ from the standard set is
  named (`atom sets that differ from the standard nucleotide …: DC_1 lacks P,
  OP1, OP2`).  That note describes the structure file and never makes a result
  partial (ADR-114).
* A result whose selection did not get everything it asked for — a named
  residue, atom, molecule or residue type not found, a modified residue the
  method cannot read (`not supported`, `not represented`), a residue or
  molecule copy you named that the completeness rule declined, or an atom
  named for a molecule that a selected copy lacks — is **partial**.  The
  completeness rule declining a residue among *all* residues of a type, or a
  copy among *all* copies of a molecule, does not make a result partial; that
  is reported in the skip records and the selection record's `excluded`.  The
  selection record of each side (`provenance.selection_request` in
  `PREFIX.meta.json`) lists what was left out as `not_found`, `unsupported`
  and `declined`, each only when non-empty (and atoms a copy lacks under
  `atom_subset.missing`), and each result in the completion record has an
  `omissions` entry (`{}` when the result is complete).  For a result written
  before these fields existed, the warnings above are read instead.
* A partial result is written and usable, but it is never reported as the
  complete selection that was asked for.  In a batch its status is
  **`partial`**: it is written like `ok`, counted separately in the summary
  line (`… 1 partial`), and its `message` names what is missing.  `biotsr run`
  and `biotsr batch` exit **4** when a result is partial and nothing failed
  (§10), and every analysis refuses partial samples unless you pass
  `--allow-partial` (§8).  `biotsr run` prints these warnings without
  `--verbose`.

### How a structure is read

Model 1 only: reading stops at the first `ENDMDL`, or at a `MODEL` record
numbered above 1.  Protein and nucleotide reading also stops at the `TER`
record of the requested chain; the `HETATM` families (`drug`, the carbohydrate
methods) and the entity side of the cross methods read past `TER`, because
ligand and glycan records usually follow it.  Amino-acid and protein reading
takes `ATOM` records, and in the atomic methods also the `HETATM` records of
TPO, SEP and PTR (*Modified amino acids*, above).  Hydrogen and deuterium are
never vertices.  Metal ions read from `HETATM` records are labelled as their elements
(manganese 100, zinc 105).

**Alternate locations: one rule for every method.**  An atom whose altLoc is
blank or `A` is kept; `B`, `C` and every later indicator are excluded, and
there is no fallback to another conformer when `A` is absent, so a residue
whose side chain exists only as `B` or `C` loses those atoms (and a method with
a completeness rule then excludes it, with the reason recorded).  Every result
records `source.altloc_policy = "blank_or_a"`, and excluded atoms are counted in
a result warning, printed with `--verbose`:

```
$ biotsr run --method amino_acid_atomic --pdb 1BRS.pdb --chain D --residue 14 --residue-name SER --output results/ --verbose
  warning: alternate locations: kept blank and A; excluded 2 atom record(s) with indicator B
```

To study another conformation, prepare a PDB file whose altLoc indicators
select it — keep only the wanted conformer, relabelled `A` or blank.  The
[user guide](USER_GUIDE.md), §12.5, has the details.

---

## 5. Two sides: cross TSR

```bash
# protein-protein, all atom
biotsr run --method protein_protein_cross --pdb 1BRS.pdb \
    --chain A --chain-2 D --contact-cutoff 4.0 --output results/

# protein-protein, Ca interface
biotsr run --method protein_protein_calpha_cross --pdb 1BRS.pdb \
    --chain A --chain-2 D --contact-cutoff 8.0 --output results/

# drug-protein: --chain is the ligand's chain, --chain-2 the protein's
biotsr run --method drug_protein_cross --pdb 3PTB.pdb \
    --chain A --chain-2 A --entity-name BEN --entity-id 1 \
    --contact-cutoff 4.0 --output results/
```

**`--contact-cutoff` and `--size-filter` are different things, and the CLI keeps
them apart:**

* `--contact-cutoff` decides **which atoms enter the pool**.  It is a
  *processing* setting: changing it leaves the `KeySpace` fingerprint untouched
  and changes the `SamplingProfile` and processing fingerprints, so results with
  different cutoffs are compared only with `--allow-mixed-processing` (§8).  Any
  positive value is accepted, and a value you set is part of the output prefix
  (`…_cutoff4`).
* `--size-filter` decides **which triangles survive**, by longest edge, after the
  geometry.  It is a modifier and is off by default.

An atom admitted because it lies 3 Å from the other side can still form a
triangle with a 20 Å edge, and every implementation keeps that triangle.  If you
want short triangles, ask for short triangles.

`--contact-filter` chooses which side the cutoff narrows (`both_sides`,
`side_a_only`, `side_b_only`; the standard is `both_sides`).
`--no-contact-filter` skips the contact scan entirely.  `--chain-2` is required
by every cross method.  The entity side of a drug, carbohydrate or nucleotide
cross follows the molecule rule of §4: with `--entity-id` one copy, without it
every copy of the name on `--chain` (the side is then described as
`entity=NAG(all copies)`).  A protein side of an atomic cross reads TPO, SEP
and PTR and names any other modified residue as not supported;
`protein_protein_calpha_cross` represents none (§4, *Modified amino acids*).
`--calpha-residues` follows the range rules of §4, *Residue identifiers and
ranges*.

Each side can also be narrowed before the contact scan.  `--atoms` narrows the
entity side of a drug, carbohydrate or nucleotide cross.  A protein side takes
at most one of `--residue-types 'ARG;LYS'` (residues of those types, all heavy
atoms), `--residues 'ASP_189;SER_190'` (named residues, all heavy atoms),
`--atom-selection 'ARG_59_NH1;LYS_27_NZ'` (named atoms) or
`--calpha-residues '50-110'` (the Cα of listed residues, for
`protein_protein_calpha_cross`, which takes only this form).  In a
protein–protein cross the unsuffixed arguments narrow the first side and
`--residue-types-2`, `--residues-2`, `--atom-selection-2` and
`--calpha-residues-2` the second; the methods with one protein side refuse the
`_2` forms.  In the output prefix each side's narrowing is preceded by its chain
(`1BRS_A_D_protein_protein_cross_A-ARG+LYS_D-ASP+GLU_cutoff5`); the two chains
of a protein–protein cross always differ (ADR-109).  The arguments each cross
method takes:

| method | takes |
|---|---|
| `drug_protein_cross`, `carbohydrate_protein_cross`, `nucleotide_protein_cross` | `--chain-2`, `--contact-cutoff`, `--contact-filter`, `--no-contact-filter`, `--entity-name`, `--entity-id`, `--entity-names`, `--atoms`, `--residue-types`, `--residues`, `--atom-selection` |
| `protein_protein_cross` | `--chain-2`, `--contact-cutoff`, `--contact-filter`, `--no-contact-filter`, `--residue-types`, `--residues`, `--atom-selection`, `--residue-types-2`, `--residues-2`, `--atom-selection-2` |
| `protein_protein_calpha_cross` | `--chain-2`, `--contact-cutoff`, `--calpha-residues`, `--calpha-residues-2` |

`--entity-names` cannot be combined with `--entity-name`, `--entity-id` or
`--atoms`.

```bash
biotsr run --method protein_protein_cross --pdb 1BRS.pdb \
    --chain A --chain-2 D --residue-types 'ARG;LYS' --residue-types-2 'ASP;GLU' \
    --contact-cutoff 5.0 --output results/
```

---

## 6. The standard profile

Every method runs under its own **standard profile**, the accepted TSR science
as set by the TSR Research Group.  There is nothing to choose: there is no
`--profile` or `--format` option, and a manifest `profile` column is refused.
`biotsr methods NAME` names the profile and describes it, the completion
summary prints it, and the completion record and the metadata sidecar record
it, so a folder of results never leaves you guessing which convention produced
it.

Output is always the native format: the key file, its `.meta.json` sidecar and
(unless `--no-triplets`) the triplet file, all of which `biotsr describe` and
`biotsr.load_result` read back.  Where a standard convention differs from a
historical TSR program — `math.pi` where some programs divided angles by the
literal 3.14, for instance — the method records under `docs/methods/` say so.

---

## 7. A whole dataset

### The manifest

**One row is one request.**  The columns are the flags without their dashes.
(A wide dataset sheet, one row per structure, is read with `--sheet`; see
*Wide dataset sheets* below.)

```csv
sample_id,pdb_id,method,chain,label,mirror,size_filter
ubq_plain,1UBQ,protein,A,ubiquitin,,
ubq_mirror,1UBQ,protein,A,ubiquitin,true,25.0
tup_sse,1TUP,protein,A,p53,,
ben,3PTB,drug,A,ligand,,
```

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/
```

Columns, by role:

**Structure** — one of these is required:

| column | meaning |
|---|---|
| `pdb_path` | path to the file, used as given |
| `pdb_id` | identifier resolved against `--pdb-dir` (`<id>.pdb` or `<id>.ent`, any case) |
| `structure_id` | override the identifier taken from the filename |

Giving both `pdb_path` and `pdb_id` is legitimate and useful: the path finds the
file, the id names the structure, so a dataset of `download_4f2q.pdb` files still
reports `4F2Q`.

**Universal**

| column | meaning |
|---|---|
| `method` | the row's method; leave the column out and choose methods with `--method` instead (never both) |
| `chain` | **required** — no method infers a chain |
| `chain_2` | the second chain, for cross methods |
| `sample_id` | the output prefix (see *Output names*); defaults to a deterministic one |
| `label` | free text, carried into the completion record and used by `--group-by label` |

**Modifiers** — `aa_grouping`, `mirror`, `sse` (truth values:
`true/yes/y/1/on` or `false/no/n/0/off`), `size_filter` (a number).

**Target selection**, whichever the method takes — `residue`, `residue_name`,
`atoms`, `residues`, `residue_types`, `entity_name`, `entity_id`, `entity_names`.

**Cross** — `contact_cutoff`, `contact_filter`, `no_contact_filter`,
`atom_selection`, `calpha_residues`, and the `_2` forms of the second side.

There is no `profile` column: a sheet that has one is refused, with the
remedy (drop it, or `--ignore-column profile`).

Rules:

1. **A blank cell is "not set", never a default.**  Defaults belong to the
   package, or to a row default you state on the command line (below).
2. **An unknown column is an error**, naming the remedies.  A column you
   believed was driving the calculation is never silently dropped.  Columns
   whose names start with `comment` are ignored, so annotated sheets can be used
   directly.
3. **A row is one request per method.**  A *wide* dataset sheet — one row per
   PDB entry carrying a dozen calculations with four different chain columns,
   as in the historical TSR dataset sheets — is not a manifest: read it with
   `--sheet` (*Wide dataset sheets*, below).
4. **Every cell belongs to a column.**  A non-blank cell beyond the header, or
   under an empty header, is refused with its line number: it would otherwise
   be dropped silently.  Usually a list was separated with commas inside an
   unquoted cell; separate list items with `;` (`1-10;20`) or quote the cell.
5. **Header names are unique.**  A duplicate header is refused.
6. **UTF-8**, with or without a byte-order mark.

Aliases accept the column names used in the historical TSR dataset sheets, so
a one-request-per-row manifest need not rename them: `protein`→`pdb_id`,
`aa_chain`/`drug_chain`/`carb_chain`/`nt_chain`→`chain`, `aa_chain_2`→`chain_2`,
`aa_group`→`label`, `aa_selection`→`residues`, `aa_ind`→`residue_types`,
`drug_name`/`carb_name`/`nt_name`→`entity_name`,
`car_atom_selected`/`carb_atom_selected`/`aa_atom_selection`→`atoms`, and others.
Two aliases of the *same* column populated in one row is an error naming both.
In a manifest every chain alias means `chain`, so a sheet with several chain
columns needs `--sheet`, which reads each method's own chain column.

### Columns with other names, and columns to carry along

A dataset sheet usually has columns BioTSR does not know.  Give each one a role:

| option | meaning |
|---|---|
| `--column FIELD=HEADER` | read manifest field `FIELD` from the column `HEADER` (repeatable) |
| `--keep-column HEADER` | carry the column through as **metadata**: it is written to each completion record and added to `batch_summary.csv`, and it never reaches a calculation (repeatable) |
| `--ignore-column HEADER` | drop the column (repeatable) |

For example, a sheet with the columns `protein,chain,group,group1`, where
`group` is a class label (`alpha`/`beta`) and `group1` its 0/1 encoding:

```bash
biotsr batch --manifest sample_details.csv --method protein \
    --column label=group --keep-column group1 \
    --pdb-dir pdb/ --output results/ --no-triplets
```

`--column label=group1 --keep-column group` works equally well, with the numeric
code as the label.  A label is an annotation: it is used by `--group-by label`,
and it never turns on `--aa-grouping` or any other scientific option.

Anything ambiguous is refused before a single row is read, with the remedy in
the message: a header that does not exist, two columns mapped to one field, one
column mapped to two fields, a mapped field that the sheet also has under its
own or an alias name (drop the mapping, or `--ignore-column` the other one), or
a column given two roles.  From Python the same options are `columns=`,
`keep_columns=` and `ignore_columns=` of `batch.run` (`from biotsr import batch`).

### Command-line values are row defaults

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/ \
    --mirror --size-filter 25.0
```

Anything the row states wins; anything it leaves blank takes the command-line
value.  That includes the chain:

```bash
biotsr batch --manifest ids.csv --method protein --chain A \
    --pdb-dir pdb/ --output results/
```

`--chain A` fills only the rows with no chain of their own, and
`batch_summary.csv` records which applied in its `chain_source` column (`row` or
`default`).  **Every row must end up with a chain.**  If some do not, the whole
batch is refused before anything runs, naming the lines: BioTSR never assumes
chain A or the first chain of a structure.  A chain the structure does not have
fails that request, and the message lists the chains it does have.

### Choosing methods on the command line

`--method NAME` is repeatable, and each source row becomes one request per
method.  In a one-request-per-row manifest every column a row fills must be one
that each listed method reads (a `residue_types` column is refused by
`protein`), so several methods are most useful on a wide sheet (below).

The methods come from **one** place.  A manifest whose rows name a method
together with `--method` is refused at planning, naming the lines (drop
`--method` to run each row's own method, or `--ignore-column method` to run the
`--method` list on every row); a manifest with neither is refused too.  With
`--method`, a row that sets none of the arguments a method needs (a row without
`residue_types` for `aa_ind`, say) does not describe that calculation: it
fails, or with `--inapplicable skip` is recorded as `inapplicable`.

### Wide dataset sheets

A wide sheet has one row per structure and the columns of many calculations
at once: several chain columns (`aa_chain`, `drug_chain`, `carb_chain`,
`nt_chain`, `aa_chain_2`) and several alternative selections (`aa_ind`,
`aa_ind_mix`, `carb_name`, …).  Which calculations to run is not in the sheet;
you choose them with `--method`, which is required:

```bash
biotsr batch --manifest sheet.csv --sheet tsr --method carbohydrate --method aa_ind \
    --inapplicable skip --pdb-dir pdb/ --download-missing --output results/ --dry-run
```

Three options do three different things, and none of them stands in for
another:

* `--method` chooses **which calculations** run (repeatable);
* `--sheet tsr` chooses **how the columns are read**: the known layout of the
  TSR dataset sheets, below;
* `--map`, `--map-file` and `--sheet custom` **adapt other headers** to that
  reading (*Your own columns*, below).

`--sheet tsr` reads the TSR dataset-sheet columns **method by method**: each
method reads only its own columns, so one sheet serves several methods whose
chain columns differ.

| method | chain | other fields |
|---|---|---|
| `protein`, `amino_acid_atomic` | `aa_chain` | — |
| `aa_selection` | `aa_chain` | `residues` ← `aa_selection` |
| `amino_acid_atom_selection` | `aa_chain` | `atoms` ← `aa_atom_selection` |
| `aa_ind`, `aa_ind_mix` | `aa_chain` | `residue_types` ← `aa_ind`, `aa_ind_mix` |
| `aa_ind_specific` | `aa_chain` | `residues` ← `aa_ind_specific` |
| `drug` (`drug_atom_selection`) | `drug_chain` | `entity_name` ← `drug_name`, `entity_id` ← `drug_id` (`atoms` ← `drug_atom_selected`) |
| `carbohydrate` (`carbohydrate_atom_selection`) | `carb_chain` | `carb_name`, `carb_id` (`atoms` ← `car_atom_selected`) |
| `carbohydrate_chain` | `carb_chain` | — |
| `nucleotide` (`nucleotide_atom_selection`) | `nt_chain` | `nt_name`, `nt_id` (`atoms` ← `nt_atom_selected`) |
| `nucleotide_chain` | `nt_chain` | — |
| `drug_protein_cross`, `carbohydrate_protein_cross`, `nucleotide_protein_cross` | the entity's chain column, `chain_2` ← `aa_chain` | the entity's name and id columns; each side's narrowing column as chosen with `--cross-side` (*Cross methods in a wide sheet*, below) |
| `protein_protein_cross`, `protein_protein_calpha_cross` | `aa_chain`, `chain_2` ← `aa_chain_2` | each side's narrowing column as chosen with `--cross-side` |

The structure comes from `protein` (`pdb_id`) and the label from `aa_group`,
which is an annotation only and never turns on `--aa-grouping` or chooses a
method.  The atom
column of the carbohydrate is spelled `car_atom_selected` in the TSR sheets and
`carb_atom_selected` in their column dictionary; either is read, and a sheet
with both is refused.

* **Participants.**  The amino-acid methods have a second participant, the
  `_2` columns (`aa_chain_2`, `aa_ind_2`, `aa_selection_2`, …).  `--participant
  1` (the default) reads `aa_chain` and the participant-1 columns (`aa_ind`,
  `aa_selection`, …), `2` reads `aa_chain_2` and the `_2` columns, and `all`
  reads both, as separate requests.  A row that fills a `_2` selection but
  leaves `aa_chain_2` blank is not run; it gets the status
  `chain_not_specified`, naming the column:

  ```
  [chain not specified: aa_ind (participant 2) needs chain, and the row leaves aa_chain_2 blank; BioTSR never assumes a chain: give it in the row, map a column to it (--column chain=HEADER or --map), or pass --chain for rows that leave it blank]
  ```

  A cross method has one participant covering both of its sides; its two
  chains and its sides' columns are described under *Cross methods in a wide
  sheet*, below.
* **Your own columns.**  `--map METHOD[.PARTICIPANT].FIELD=HEADER` (repeatable)
  adds an entry or overrides one, and `--map-file CSV` gives the same entries
  as a file with the header `method,participant,field,column`.  `--column
  FIELD=HEADER` maps the fields every request shares (`pdb_id`, `pdb_path`,
  `structure_id`, `sample_id`, `label`).  `--sheet custom` starts from no
  layout at all and reads only what you map:

  ```bash
  biotsr batch --manifest mine.csv --sheet custom --method carbohydrate \
      --column pdb_id=entry --map carbohydrate.chain=glycan_chain \
      --map carbohydrate.entity_name=sugar --map carbohydrate.entity_id=sugar_id \
      --pdb-dir pdb/ --output results/
  ```

  An unknown method, field or column, a field mapped twice, a participant other
  than 1 or 2, participant 2 of a cross method and a malformed `--map` are each
  refused with the remedy.
* **Rows that do not describe a calculation.**  A row with none of a method's
  defining columns filled — a ligand row for `carbohydrate`, say; for an
  `*_atom_selection` method the atom column — fails by default, or with
  `--inapplicable skip` is recorded with the status `inapplicable`.  It is never
  silently left out.  A row that describes a calculation but leaves its chain
  blank (after any `--chain` default) is not run: its status is
  `chain_not_specified`, naming the column, and the batch exits 4.  BioTSR
  never infers a chain.
* **Shared calculations.**  Rows asking for the same calculation of the same
  sample — the same protein chain repeated on every ligand row — are computed
  once, and the summary lists every source row (`4;5;6;7`).  Rows with
  distinct sample ids stay distinct samples (*Two rows, one calculation*,
  below).
* **Unused columns.**  Columns no selected method reads are listed (`unused`)
  in the summary printed at the start; `comment*` columns are ignored.
* **What a sheet refuses.**  A `method` column (methods come from `--method`;
  pass `--ignore-column method` to drop it) and `--chain` (each method reads
  its own chain column; map one with `--map METHOD.chain=HEADER`).
  `--cross-side` is refused without `--sheet`.

A `sample_id` column names the request's prefix.  When one source row
becomes several requests, each prefix is `<sample_id>_<method>`, with `_p2`
added for participant 2.

### Cross methods in a wide sheet: choose each side

`--sheet tsr` reads a cross method's chains and its molecule's name and
identifier from the row.  Each side can also be narrowed by several columns of
the sheet, which were written for one-sided methods, and the sheet does not
say which one a cross calculation means.  BioTSR neither picks one nor ignores
them: while a row populates any candidate column of a side, the run is refused
until you choose that side.

| method | side | candidate columns |
|---|---|---|
| `drug_protein_cross` | `drug` | `drug_atom_selected` |
| `carbohydrate_protein_cross` | `carbohydrate` | `car_atom_selected` (`carb_atom_selected`) |
| `nucleotide_protein_cross` | `nucleotide` | `nt_atom_selected` |
| the three above | `protein` | `aa_ind`, `aa_ind_mix`, `aa_ind_specific`, `aa_atom_selection` |
| `protein_protein_cross` | `1` / `2` | `aa_ind`, `aa_ind_mix`, `aa_ind_specific`, `aa_atom_selection` / their `_2` forms |
| `protein_protein_calpha_cross` | `1` / `2` | `aa_selection` / `aa_selection_2` |

(`biotsr.cli.sheet.CROSS_SIDE_COLUMNS` holds the same table.)  On
`tests/data/manifests/wide_sheet.csv`, whose 6G9M row fills `aa_ind`,
`aa_ind_mix`, `aa_ind_specific`, `aa_atom_selection` and `drug_atom_selected`:

```
$ biotsr batch --manifest wide_sheet.csv --sheet tsr --method drug_protein_cross --inapplicable skip --pdb-dir pdb/ --output cross/ --dry-run
biotsr: error: the sheet populates columns that could narrow a cross calculation, and BioTSR does not choose among them: the drug side of drug_protein_cross: drug_atom_selected on line 2 | the protein side of drug_protein_cross: aa_atom_selection on line 2; aa_ind on line 2,3; aa_ind_mix on line 2; aa_ind_specific on line 2
```

`--cross-side METHOD.SIDE=COLUMN` (repeatable) narrows that side by that
column, and `--cross-side METHOD.SIDE=whole` uses the whole side.  A `--map` of
the side's field (`--map protein_protein_cross.residue_types=aa_ind`) counts
as the choice too; choosing one side twice, by two `--cross-side` options or by
`--cross-side` and `--map`, is refused, as are an unknown side, a column that
does not narrow that side, and a method not among the `--method` choices.

```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method drug_protein_cross \
    --cross-side drug_protein_cross.protein=aa_ind \
    --cross-side drug_protein_cross.drug=whole \
    --inapplicable skip --pdb-dir pdb/ --output cross/ --dry-run
```

The preview prints, under each request line, one indented line per side with
its resolved selection, and the populated columns the choice leaves unused:

```
  1    2        drug_protein_cross           1   A+A     drug, chain A: ESW copy 404, al... 6G9M_A_A_drug_protein_cross_ESW_404_s1-TYR  [would run]
         drug, chain A: ESW copy 404, all atoms
         protein, chain A: residue_types=TYR
         not applied: drug_atom_selected (the drug side is used whole); aa_ind_mix (the protein side uses aa_ind); aa_ind_specific (the protein side uses aa_ind); aa_atom_selection (the protein side uses aa_ind)
```

A protein–protein cross names its sides `side 1` and `side 2`
(`side 1, chain A: residues=LYS_19;TYR24`).  In `batch_summary.csv` the
`selection` of a cross request is both sides, `side text | side text`
(`drug, chain A: ESW copy 404, all atoms | protein, chain A: whole chain`),
and its `message` repeats the `not applied: …` list.  A row whose chosen
column is blank does not describe the narrowed calculation: it is
inapplicable, never computed on the whole side.  From Python, the same choice
is `cross_sides=` of `batch.run` (API guide, §8).

### Output names

Every calculation writes its files directly into `--output` under one
**prefix**: `PREFIX.keys`, `PREFIX.meta.json`, `PREFIX.triplets` (unless
`--no-triplets`) and the completion record `PREFIX.run.json`.  The prefix is
the row's `sample_id` (or `biotsr run --sample-id`), or else a deterministic
name built from the request:

```
STRUCTURE_CHAIN[_CHAIN2]_METHOD[_TARGET][_MODIFIERS]
```

| part | examples |
|---|---|
| target: one copy, every copy | `7LF9_D_carbohydrate_NAG_1`, `7LF9_D_carbohydrate_NAG_allcopies` |
| with an atom subset | `…_NAG_1_atoms-C1+C2+N2` |
| residue types, named residues | `1BRS_A_aa_ind_LYS`, `…_aa_ind_mix_TYR+ARG`, `…_aa_ind_specific_LEU15+LEU50` |
| residue list, named atoms | `…_aa_selection_res-1-10+20`, `…_amino_acid_atom_selection_atoms-LYS11.NZ+LEU15.CB` |
| one residue, a sweep | `…_amino_acid_atomic_LEU_15`, `…_perresidue` or `…_perresidue-LEU` (members `…_perresidue-LEU_LEU_15`) |
| entity vocabulary | `…_carbohydrate_chain_names-NAG+FUC` |
| cross narrowing | `1BRS_A_D_protein_protein_cross_A-ARG+LYS_D-ASP+GLU` |
| modifiers and settings | `aagroup`, `mirror`, `size25`, `sse`, `cutoff4.5`, `filter-side_a_only`, `nocontactfilter` |

So a plain and a mirrored run of one chain get different prefixes
(`1UBQ_A_protein`, `1UBQ_A_protein_mirror`).  When the readable form cannot
represent the request exactly — a character outside `A-Z a-z 0-9 . _ + -`, a
prime (written `p`), a lower-case chain identifier, a residue name containing a
digit in a compact token — or is longer than 96 characters, a request
fingerprint `_h<10 hex digits>` is appended, so two different requests never
share a name.

A sample id must be **one file-name stem**: no `/` or `\`, no absolute path,
not `.` or `..`, no leading `.` (reserved for hidden and staging names), not
empty, no leading or trailing spaces, no control characters, at most 200 bytes,
and not one of the names the folder uses for itself (`analysis`,
`batch_summary.csv`, `batch_failures.csv`, `fetch_report.csv`,
`.biotsr-batch.lock`, in any letter case).  An invalid id fails its request in a
batch and exits `2` from `biotsr run`; it is checked while planning, so
`--dry-run` reports it without writing.

A batch in which two **different** requests would write the same prefix is
refused before anything runs, naming the rows; give those rows their own
`sample_id` values.  Prefixes that differ only in letter case count as the
same, because they are the same files on macOS's default file system, and so
does a prefix that extends a per-residue sweep's prefix.  Identical requests
of one sample are not a clash: they are computed once (below).

### Two rows, one calculation: samples and computations

A **sample** is one entry of the output folder: its prefix and its label.  A
batch merges two requests only when they are **the same calculation of the
same sample** — method, structure, chains, selection, modifiers and sample id
all equal, the sample id being the row's `sample_id` or, when that is blank,
the automatic prefix, which names exactly the calculation.  The merged request
is computed once; its `rows` lists every source row (and `participant` every
participant), and its `message` says
`one sample from 2 source requests (line 2, 3)`.

* Rows with **distinct `sample_id` values** are distinct samples, even when
  the calculation is identical: each is computed and written under its own
  prefix, and an analysis reads both.
* Rows that would merge but carry **different labels** are refused before
  anything runs, because one output carries one label: give them distinct
  `sample_id` values, or correct the label.
* A `--keep-column` whose values differ between merged rows records each
  row's value: `line 2=first entry; line 3=second entry`.

```csv
pdb_id,chain,method,residue_types,sample_id,label
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_a,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_b,ubiquitin
```

Lines 2 and 3 are one sample, `1UBQ_A_aa_ind_LYS`, computed once with `rows`
`2;3`; lines 4 and 5 are two samples, `ubq_tyr_a` and `ubq_tyr_b`, each
computed and written.  Had lines 2 and 3 different labels:

```
biotsr: error: 1 calculation(s) are requested by rows with different labels: 1UBQ_A_aa_ind_LYS: label 'group1' on line 2; label 'group2' on line 3
  one output is one sample and carries one label; give those rows distinct sample_id values to keep them as separate samples, or correct the label
```

### Checking a batch before running it

```bash
biotsr batch --manifest sample_details.csv --method protein \
    --column label=group --keep-column group1 \
    --pdb-dir pdb/ --output results/ --dry-run
```

`--dry-run` plans the batch and changes nothing: no download, no output
folder, no file.  It prints how many source rows and requests there are, where
the chains come from, the columns read under other names, kept as metadata or
unused, which structures are present or would be downloaded, and the outcome
of every request, one line each:

```
  req  rows     method                       pt  chains  selection                          prefix / outcome
  5    4        carbohydrate                 1   D       entity_name=NAG; entity_id=1       7LF9_D_carbohydrate_NAG_1  [would run]
  6    4;5;6;7  aa_ind                       1   B       residue_types=TYR                  7LF9_B_aa_ind_TYR  [would run]
  9    7        carbohydrate                 1   D       entity_name=MAN                    7LF9_D_carbohydrate_MAN_allcopies  [would run]
  10   8        carbohydrate                 1   -       -                                  -  [inapplicable: the row has no carb_name / carb_id for carbohydrate]
```

(the request number, its source row(s), method, participant, chains, selection,
output prefix and outcome: would run, would be downloaded then run, would be
skipped, inapplicable, or would fail with the reason).  A cross request is
followed by one indented line per side with its resolved selection
(*Cross methods in a wide sheet*).  What can only be
checked once a structure is read (that the chain exists in it, that a selection
matches its atoms) is listed as deferred.  It exits `0` when every request
would run, and `3` when some would fail.

### What a batch writes

```
results/
    batch_summary.csv               one row per request, in order, with its status
    batch_failures.csv              only when something failed
    fetch_report.csv                only with --download-missing
    7LF9_D_carbohydrate_NAG_1.keys  every calculation's files, side by side
    7LF9_D_carbohydrate_NAG_1.meta.json
    7LF9_D_carbohydrate_NAG_1.triplets
    7LF9_D_carbohydrate_NAG_1.run.json
    ...
```

`batch_summary.csv` columns: `request`, `rows`, `participant`, `sample_id` (the
output prefix), `status`, `method`, `profile`, `structure_id`, `chain`,
`chain_2`, `selection`, `n_results`, `n_triangles`, `n_unique_keys`, `output`
(folder and prefix), `message`, `chain_source`, then any `--keep-column`
columns; `batch_failures.csv` has the same columns.  `rows` holds the line
numbers in your manifest (`4;5;6;7` for a calculation several rows share), so
every result and every failure points back at the lines that caused it.  For a
cross request `selection` holds both sides, `side text | side text`.
`status` is one of:

| status | meaning |
|---|---|
| `ok` | computed and written |
| `partial` | computed and written like `ok`, but something the request asked for — a named item, or a residue of the chain or side it covers — is not in the result: not found, not supported, not represented, or excluded by the completeness rule (§4, *Empty and partial selections*); `message` names each.  Counted separately (`… 1 partial`).  Not a failure, but the batch exits 4 when nothing failed (§10) |
| `no_triangles` | computed and written, but the selection had fewer than three points; `message` says so.  Not a failure.  When the selection also lacked something it named, the result is partial as well: the summary line counts it (`2 with no triangles (1 of them partial)`) and the batch exits 4 when nothing failed |
| `skipped` | an existing output was verified to be this exact request, computed under the current processing rules (`--if-exists skip`).  When the kept output is itself partial, `message` reads `verified existing output of the same request; kept partial result: …`, naming what it lacks, and the summary line counts it (`2 skipped (1 of them partial)`) |
| `inapplicable` | the source row does not describe this method's calculation, and `--inapplicable skip` was given; `message` says why.  Not a failure |
| `chain_not_specified` | the request needs a chain (a cross method also `chain_2`) that neither its row nor `--chain` gives; nothing was computed and no chain was assumed; `message` names the column.  Not a failure, but the batch exits 4 when nothing failed (§10) (ADR-112) |
| `failed` | not computed, or its worker process died; `message` says why |
| `interrupted` | it was being computed when the run was stopped |
| `not_started` | the run stopped (or the worker pool failed) before it began |

The batch keeps only this one line per request in memory: each result is
written and released before the next, so a dataset of thousands of structures
does not grow the process.

### When one request fails

The default is **continue**: every request is attempted, every failure is
recorded against its own source row(s) in `batch_failures.csv`, and the process
exits **3**.

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/ \
    --on-error stop     # halt at the first failure instead
```

With `stop`, a request that can be seen to fail before anything is computed (an
unknown method, a missing structure) stops the batch before any computation,
and the command exits **2**; otherwise the batch stops after the first request
that fails while computing, and exits **3**.  Either way the summary is written
and lists every request, and nothing is omitted silently.  A failed request
leaves **no** files in the output folder — not an empty file, not a partial one.

### Resuming, overwriting and what is never touched

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/ \
    --if-exists skip
```

`--if-exists` is `fail` (the default — refuse and change nothing), `overwrite`
(replace the calculation's files) or `skip`.

`skip` keeps an existing calculation only when it is **verified** to be the
result of exactly this request.  Its completion record `PREFIX.run.json` holds
the full normalised request — the structure's SHA-256 (its content, not its
path), the method, both chains, every selection, every modifier, every
contact setting and the processing rules that admit atoms — and `skip` requires the record to be complete, to match,
every file it lists to be present at the size recorded when it was written (so
a file cut short or appended to is caught), and the existing output to cover
what is asked for now: a calculation written with `--no-triplets` does not
satisfy a request for triplets.  A calculation that fails the check is reported
as a failed request, with the difference named; it is never silently kept or
replaced.  Its name, or matching fingerprints, are never enough on their own.

**Processing rules.**  The same request can produce different numbers when
the rules that admit atoms change: the amino-acid atomic methods, and the
atomic cross methods' protein sides, now read the modified amino acids TPO,
SEP and PTR written as `HETATM`, and outputs written by earlier builds of
the same version left them out (ADR-101).  So the request in every completion
record holds the rules it was computed under (`request.processing`):

```json
"processing": {"altloc_policy": "blank_or_a", "completeness_policy": "exact_legacy_atom_count",
               "input_rules": "input_rules_v2", "modified_residue_policy": "hetatm_amino_acids:PTR,SEP,TPO"}
```

That is `aa_ind_specific`'s.  The input rules (`input_rules_v2`) are recorded
since 1.0.0rc6 (ADR-108), and a method's completeness rule since 1.0.0rc7
(ADR-114) — a method without one, such as every Cα method, records none.  The
Cα methods read more modified residues (ADR-107), so their records hold
`"modified_residue_policy": "hetatm_amino_acids:CME,CSD,CSO,LLP,M3L,MSE,NIY,PTR,SEP,TPO,TRF;ca_element:C_or_blank"`.

`skip` keeps a calculation only when these equal the current rules; the
package version is not used for this, because one version number has covered
more than one rule set.  A record that does not state them is never taken to
be current: in a batch the request is `failed` (exit 3) with this message,
and `biotsr run` stops with it (exit 2).  `skip` also requires each result's key space
to be the one this version computes for the request, so a Cα output written
before 1.0.0rc5 — radix 20, no modified residues — and every atomic output
written before 1.0.0rc7 — atoms labelled by name (ADR-113) — is refused with
`its keys were computed in a different key space (…)` and recomputed only with
`--if-exists overwrite`; when the processing rules differ too, both reasons
are given, key space first.  A Cα output of 1.0.0rc5 or rc6 is kept.  Here
`old/` holds an atomic output written by an earlier build, before the
processing rules were recorded:

```
$ biotsr run --method aa_ind_specific --pdb SYNPHOS.pdb --chain A --residues 'SER_4;TPO_3;TYR_6' \
      --sample-id SYNPHOS_A_specific --no-triplets --output old/ --if-exists skip
biotsr: error: output SYNPHOS_A_specific exists but cannot be reused: its keys were computed in a different key space (fingerprint 1187735d9085; this version computes 1ecbef1b8334 for this request), so they do not mean what this version's keys mean; and its completion record does not state the processing rules it was computed under (it was written before they were recorded, and outputs that old may have excluded the modified amino acids TPO, SEP and PTR written as HETATM), so it is not taken to match the current rules (altloc_policy=blank_or_a, completeness_policy=exact_legacy_atom_count, input_rules=input_rules_v2, modified_residue_policy=hetatm_amino_acids:PTR,SEP,TPO)
  --if-exists skip only keeps a calculation produced by exactly this request; pass --if-exists overwrite to recompute it, choose another --output, or give the row its own sample_id
```

A record that states different rules is refused the same way, as `it was
computed under different processing rules (…)`.  Recompute such an output with
`--if-exists overwrite`, or into a new `--output`.  When only its processing
rules differ it can still be analysed as it is, deliberately, with
`--allow-mixed-processing` (§8); when its key space differs, nothing compares
it with current results.

**A kept result that is partial stays partial.**  When `skip` keeps an output
that is itself partial (§4, *Empty and partial selections*), its row is
`skipped`, its message says so, and the command exits 4, as it did when the
result was computed:

```
biotsr: warning: m.csv line 3 (aa_ind_mix): verified existing output of the same request; kept partial result: requested but not found: residue type SER

requests   0 ok, 2 skipped (1 of them partial), 0 failed
```

`biotsr run --if-exists skip` does the same for its one calculation.

**How a calculation is installed.**  It is computed into a staging directory
private to the invocation (`.biotsr-staging-<pid>-<host>-<token>/` inside
`--output`), then installed under a per-prefix lock: an intent record is
written (`"complete": false`, listing the files to come), each file is moved
into place in one step, files of an older version of the same calculation that
the new one no longer produces are removed, and the record is rewritten with
`"complete": true` and every file's size.  The record is written last, so a
calculation is complete exactly when its record says so.

**What is never touched.**  BioTSR writes, replaces and removes only files that
a completion record of that prefix lists.  Anything else in the folder — your
own notes, another program's output, a file whose name differs from a planned
one only in letter case — is left alone, even with `--if-exists overwrite`: a
calculation that would have to replace such a file is refused, naming it, and a
name that collides case-insensitively with an existing file is refused before
anything is computed.  Symbolic links are never written through, replaced or
deleted.

**Interruption.**  Ctrl-C or SIGTERM (a scheduler's time limit) stops the batch
at any point after planning — while downloading, while checking existing
outputs, or while computing.  It removes only its own staging directories,
writes the summary with `interrupted` and `not_started` rows (and
`fetch_report.csv`, if it was downloading), and exits `130`.  Downloads stop
promptly: queued ids never start, no retry follows, and a connection being made
or a response being received is stopped at once; only a host-name lookup or a
TLS handshake step in progress is waited for (`--timeout` limits each wait,
not a whole download).

**A worker that dies.**  If a worker process is killed — by the out-of-memory
killer, for example — the summary is still written: finished requests keep
their results, requests that were running are `failed` with a message naming
the pool failure, requests never started are `not_started`, and the batch
exits `3`.

Either way, re-running with `--if-exists skip` finishes the rest.

**What cannot be caught.**  SIGKILL and power loss end the process before it
can write a summary.  What it leaves is either a staging directory, which a
later run on the same machine removes once that process is gone (ones left by
another machine are kept and reported), or an incomplete completion record,
which every reader refuses.  An incomplete calculation is computed again by
`--if-exists skip` or `overwrite`; the default `fail` refuses it, naming it.

Several invocations may write into one `--output` at the same time, provided
their prefixes differ.  If two finish the *same* prefix, the first one's files
stay and the second reports the clash; only `--if-exists overwrite` replaces a
finished calculation.  Batch reports are different: there is one
`batch_summary.csv` per output folder, so a second `biotsr batch` into a folder
where another batch is still running is refused before it computes anything (a
`.biotsr-batch.lock` file names the running batch; one left by a process that
has exited on the same machine is taken over).  Give concurrent batches their
own `--output`, as the array job below does.  `--dry-run` takes no lock.

The folder's reports describe its **latest** batch: each `biotsr batch` into
it rewrites `batch_summary.csv` for its own requests, writes
`batch_failures.csv` when a request failed and removes one left by an earlier
run when none did.  A calculation's own files are never changed by
`--if-exists skip`, whether it keeps them or refuses them; a refused request
says so (`… exists but cannot be reused: … (it is left untouched; pass
--if-exists overwrite to recompute it, or choose a new --output)`).
`--dry-run` writes nothing.

### Several requests at a time

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/ --jobs 4
```

`--jobs N` computes N requests at a time in worker processes; `--jobs auto` uses the
CPUs this process may use — the scheduler's allocation (Slurm, PBS, LSF, Grid
Engine), its CPU affinity, its cgroup quota or the machine's CPU count,
whichever is smallest — rather than every CPU in the machine.  More workers
than requests are never started.  The default is `--jobs 1`: one request at a
time, in this process.  More than one worker uses joblib, which is installed with
BioTSR; in an environment where it has gone missing, `--jobs 4` stops with the
command that reinstalls it rather than quietly running one request at a time.
The results are identical either way: the same files, byte for byte, and the
same summary in the same (manifest) order; only the command recorded in the
completion records differs.

Each worker holds one request's calculation at a time, so peak memory is about
N times that of one request (see *Memory and disk* below).

### Memory and disk

What dominates a calculation's memory is the per-triangle (triplet) records, and
whether they are held depends on how it runs:

| run | triplet records |
|---|---|
| `biotsr run`, `biotsr batch` and `batch.run(...)` with triplets, every method except the per-residue sweep below, cross methods included | **streamed**: each line is written to the triplet file as its triangle is computed; only the key counts stay in memory |
| `--no-triplets`, and the keys-only samples of `analyze common-keys` | never built; only the key counts are kept |
| `amino_acid_atomic` without `--residue` (the per-residue sweep, one result per residue) | **held** for every residue until the sweep finishes |
| Python: `biotsr.analyze(...)`, `analyze_cross(...)` or `analyze_entities(...)` with the default `keep_triplets=True` | **held** on the result; pass `keep_triplets=False, triplet_sink=...` to stream (API guide, *Very large results*) |

Always in memory, whatever the mode: the run's key counts (one entry per
distinct key, not per triangle) and, in a batch, one summary line per request.
The key file, the `.meta.json` sidecar, the completion record, `batch_summary.csv` and
`batch_failures.csv` are written whole at the end of their run or batch.
`biotsr analyze common-keys` holds every sample's key counts at once, so its
memory grows with samples × distinct keys per sample.

Measurements, to give the scale (not a guarantee for other inputs): 1B8X
chain A, 2.9 million triangles, method `protein`, one `biotsr run` process,
macOS arm64, Python 3.13.  Streamed, the peak resident memory was 0.54 GB;
before triplets were streamed, the same run held its triplet list and peaked at
3.05 GB, writing a byte-identical file.  Other structures and methods will
differ; a held run grows with the number of triangles.

Streaming saves memory, not disk: the triplet file is exactly as large either
way, and it is by far the largest output.  In the run above the triplet file
was 705 MB and the key file 14 MB.  For a large dataset, plan disk for the
triplet files or leave them out with `--no-triplets`; estimates projected from a
handful of structures are only a rough guide, because triangle counts grow with
the cube of the number of points selected.

### Downloading structures

BioTSR never uses the network unless asked.  To fetch what a manifest needs:

```bash
biotsr batch --manifest sample_details.csv --method protein \
    --column label=group --keep-column group1 \
    --pdb-dir pdb/ --download-missing --output results/ --no-triplets
```

`--download-missing` downloads every structure named by `pdb_id` that is not in
`--pdb-dir` yet, before computing, and writes `fetch_report.csv`; a structure
that cannot be downloaded fails its rows and is reported, and the others go
ahead.  `--download-workers`, `--timeout` and `--retries` control the
downloads, independently of `--jobs`, with the same rules as `biotsr fetch`
(§7a).  A cached file
that looks truncated (no `END` record, or no longer matching the checksum
recorded when it was downloaded) fails its rows rather than being used; `biotsr
fetch --refresh` replaces it.  To download without computing, see §7a.

---

## 7a. Getting structures

```bash
biotsr fetch 1UBQ 4HHB --pdb-dir pdb/
biotsr fetch --id-file ids.txt --pdb-dir pdb/
biotsr fetch --manifest sample_details.csv --pdb-dir pdb/ --download-workers 8
biotsr fetch --manifest other.csv --column pdb_id=PDB --pdb-dir pdb/
biotsr fetch --manifest sample_details.csv --pdb-dir pdb/ --dry-run
```

`--manifest` reads **one column** — `pdb_id`, an alias such as `protein`, or the
column named by `--column pdb_id=HEADER` — and nothing else, so any dataset sheet
works as it is, whatever other columns it has.  Ids are case-insensitive and each
entry is downloaded once however often it appears.

Each entry is downloaded from RCSB over HTTPS as `.pdb.gz` and installed as
`<ID>.pdb` only after every check passes: the transfer is complete, the gzip
stream ends properly and its integrity check matches, and the decompressed file
has atom records, ends with an `END` record, names the right entry and is read
by the same reader every calculation uses.  The file is written under a
temporary name and moved into place in one step, so an interrupted download
never leaves a partial file behind.  Transient failures (timeouts, dropped
connections, server errors, rate limits) are retried with exponential backoff.

**Download settings**, validated identically by `fetch`, `batch`, `analyze`
and `biotsr.fetch_structures` — an invalid value is refused, never replaced by
a default or clamped:

| option | values | meaning |
|---|---|---|
| `--download-workers N` | whole number, 1–16 (default 4) | concurrent downloads, independent of `--jobs` |
| `--timeout SECONDS` | finite, greater than 0 (default 30) | how long to wait for a connection, or for the next piece of a response.  It limits each connection and each read, **not** a whole download |
| `--retries N` | whole number, 0 or more (default 3) | further attempts after a transient error |

A failed attempt is retried after an exponential backoff of up to 60 s.
Because `--timeout` limits each wait rather than the whole transfer, a server
that keeps sending, however slowly, is never timed out: no total time limit
applies to one download.  Ctrl-C or SIGTERM stops it promptly (see
*Interruption* above).

A valid file already in `--pdb-dir` is reused and never overwritten.
`--refresh` downloads again and replaces the cached file only after the new one
passes every check; if the refresh fails, the previous file is kept.  Each
download's source URL, UTC time and SHA-256 are recorded in
`--pdb-dir/.biotsr-fetch/`, and every row's outcome in `fetch_report.csv`:

| status | meaning |
|---|---|
| `downloaded` / `refreshed` | a new, verified file is installed |
| `cached` | a valid file was already there |
| `would_download` | `--dry-run` only: not in the cache, would be downloaded |
| `invalid_cached` | a file is there but looks incomplete or changed; it is not used (use `--refresh`) |
| `not_found` | no such entry |
| `pdb_format_unavailable` | RCSB has the entry only as mmCIF, which BioTSR does not read |
| `invalid_id` | not a PDB id |
| `failed` | the download did not succeed; `message` says why |
| `refresh_failed` | `--refresh` did not succeed; the previous file is kept, but the id counts as a failure |
| `interrupted` / `not_started` | the command was stopped while, or before, fetching it |

An obsolete entry (one with an `OBSLTE` record) is downloaded like any other,
with a note in `message`.

`biotsr fetch` exits `0` only when every id was `downloaded`, `refreshed` or
`cached`, and `3` otherwise — `refresh_failed` included.  Ctrl-C or SIGTERM
stops it with the report written and exits `130`.  From Python:
`biotsr.fetch_structures(["1UBQ", "4HHB"], "pdb/", workers=4)`.

---

## 8. Analysis: the keys a set of samples shares

```bash
biotsr analyze common-keys --manifest samples.csv --pdb-dir pdb/ --output results/
```

```
common keys 9523 over 3 samples
analysis   results/analysis/common_keys
```

Writes `results/analysis/common_keys/all_samples.json` (the key list, the
per-sample report — method, selection, processing settings, fingerprints and
warnings — and a `comparison` summary) and `all_samples.csv` (`sample_id`,
`method`, `total_keys`, `common_keys`, `sum_common_keys_freq`,
`key_space_fingerprint`, `processing_fingerprint`, `sampling_fingerprint`).
With `--manifest` the samples are computed exactly as `biotsr batch` computes
them — the same columns, `--sheet` and `--method` options, defaults, `--jobs`
and `--download-missing` apply — into the same shared folder, but **keys
only**: no triplet files are written, so there is no `--no-triplets` option.
They are then read back from the outputs written, so a calculation kept by a
verified `--if-exists skip` counts like a freshly computed one.  The analysis
holds every sample's key counts at once, so its memory grows with samples ×
distinct keys per sample.

Per group, using the label (a manifest's `label` column, `aa_group` in a TSR
sheet, or the label recorded in each completion record with `--from-results`):

```bash
biotsr analyze common-keys --manifest samples.csv --pdb-dir pdb/ --output results/ \
    --group-by label
```

`--on-error` defaults to **stop** here, not continue.  An intersection silently
taken over a reduced sample set is a wrong answer, not a partial one.

### What has to match, and what never blocks

Every analysis — common keys here, and in Python feature matrices, occurrence
tables, similarity and clustering — applies one comparison rule, to computed
and loaded results alike:

1. **Key definitions must match, always.**  The samples must share one key
   space: label and point scheme, theta and distance bins, encoding, ordering,
   mirror convention.  This is decided from the key-space definition, never
   from method names, so methods that share a key space compare (every atomic
   method — `aa_ind`, `drug`, `carbohydrate`, … — labels atoms by element in 29 ×
   18 bins), and the same method with `--mirror` or `--aa-grouping` does not.
   A Cα result written before 1.0.0rc5 (residue labels at radix 20), or an
   atomic result written before 1.0.0rc7 (atoms labelled by name), has a
   different key space from a current one, so the two are never compared:
   compute it again (ADR-107, ADR-113).  No option relaxes this — neither
   `--allow-mixed-processing` nor `--allow-partial`.
2. **Processing settings must match unless you say otherwise.**  How triangles
   are drawn from what was selected — candidate strategy, completeness and
   composition rules, size filter, degenerate policy, sequence-label policy,
   contact cutoff, contact filter and rule, contact projection, pool identity
   and order, the alternate-location rule and the modified-residue rule
   (`modified_residue_policy`: which modified amino acids written as `HETATM`
   were read) — must agree, or the comparison is refused.
   `--allow-mixed-processing` compares such a set anyway; the report then
   carries a warning naming each differing setting and records every sample's
   settings.  A result written before the modified-residue rule was recorded
   shows it as `not recorded`, which is never taken to be the current rule:
   comparing it with a current result needs `--allow-mixed-processing`, and
   the message names `modified_residue_policy ('hetatm_amino_acids:PTR,SEP,TPO'
   vs not recorded)`.  Its key definitions are unchanged, so it still shares
   the key space.
3. **Partial samples must be acknowledged.**  A sample whose selection did not
   get everything it asked for (§4, *Empty and partial selections*) is not the
   selection that was requested, so a set containing one is refused, naming
   each partial sample and what it lacks, unless you pass `--allow-partial`.
   With it nothing is dropped: the report's `comparison` records
   `partial_allowed` and `partial_samples` (each sample's omissions), each
   sample's entry records `partial` and `omissions`, a warning names them, and
   the command exits 4.  This is separate from `--allow-mixed-processing`:
   neither implies the other, and neither relaxes step 1.
4. **What was selected never blocks a comparison.**  Structures, chains,
   amino-acid types (LYS-only against TYR-only), ligands, carbohydrates,
   nucleotides, residue numbers, atom subsets, one copy against all copies, and
   the method name are what a dataset is made of; each sample's selection is
   recorded in the report.  A *complete* LYS-only sample and a *complete*
   TYR-only sample compare without any option; only a selection that fell
   short of its own request needs one (step 3).
5. **Samples need not share any key.**  An empty intersection is a valid
   result.

```
$ biotsr analyze common-keys --manifest mixed.csv --pdb-dir pdb/ --output results/
biotsr: error: samples 'plain' and 'filtered' share key definitions but were processed with different settings: size_filter (None vs 25.0). ...
  pass --allow-mixed-processing to compare them anyway; the report then records each sample's settings and a warning
```

Two `aa_ind` samples of different residue types, by contrast, compare without
any option:

```bash
biotsr run --method aa_ind --pdb 1BRS.pdb --chain A --residue-types LYS --output types/ --no-triplets
biotsr run --method aa_ind --pdb 3PTB.pdb --chain A --residue-types TYR --output types/ --no-triplets
biotsr analyze common-keys --from-results types/ --output types_ck/
```

A partial sample is refused until it is acknowledged.  Chain B of the
synthetic structure `SYNPHOS.pdb` has no SER, so this batch writes `chainA` as
`ok` and `chainB` as `partial` and exits 4:

```csv
pdb_id,chain,method,residue_types,sample_id
SYNPHOS,A,aa_ind_mix,TYR;LYS;SER,chainA
SYNPHOS,B,aa_ind_mix,TYR;LYS;SER,chainB
```

```bash
biotsr batch --manifest m.csv --pdb-dir . --output out/ --no-triplets
biotsr analyze common-keys --from-results out/ --output a/
```

```
biotsr: warning: m.csv line 3 (aa_ind_mix): requested but not found: residue type SER

requests   1 ok, 0 skipped, 0 failed, 1 partial
...
biotsr: error: 1 of the 2 samples are partial -- their selections did not get everything they asked for: chainB (not found: SER). ...
  pass --allow-partial to analyse them as computed; the report then records what each one lacks. Or fix the requests (or the structures) so every selection is complete. Partial results stay in out/ either way
```

The analysis exits 2.  With `--allow-partial` it runs, warns, and exits 4:

```bash
biotsr analyze common-keys --from-results out/ --output a/ --allow-partial
```

```
common keys 258 over 2 samples
biotsr: warning: 1 of the 2 samples are partial and were explicitly allowed: chainB (not found: SER). They are compared as computed, not as the complete selections that were requested.
analysis   a/analysis/common_keys
```

In Python every analysis function takes `allow_mixed_processing=True` and
`allow_partial=True`.  `--allow-mixed-processing` was called
`--allow-mixed-sampling` (`allow_mixed_sampling=` in Python) in 1.0.0rc3,
where selections needed it too; the old flag now stops with an error naming
the new one.

### Analysing results you already computed

`--from-results` reads previously saved results instead of recomputing them:

```bash
biotsr run -m protein --pdb 1UBQ.pdb -c A -o results/ --no-triplets
biotsr run -m protein --pdb 1BRS.pdb -c A -o results/ --no-triplets
biotsr analyze common-keys --from-results results/ --output analysis_out/
```

```
loaded     2 saved result(s) from 2 calculation(s) in results/
common keys 20511 over 2 samples
analysis   analysis_out/analysis/common_keys
```

The folder is searched recursively for **completion records**
(`PREFIX.run.json`), so point it at the `--output` folder a `biotsr run` or
`biotsr batch` wrote, or at a parent of several.  Hidden directories (staging
among them) and `analysis/` folders are never read.  An incomplete record, a
listed file that is missing or no longer has its recorded size, and a `.keys`
file that no complete record lists are each refused by name, so an incomplete
output is never analysed.  Results in the layout of 1.0.0rc3 or earlier
(per-sample directories, names like `1UBQ_A.keys_29_35`) are refused with a
message saying they must be recomputed.

**Choosing calculations.** A shared folder often holds several key spaces (Cα
and atomic methods side by side), which cannot be compared with each other.
`--method NAME` (repeatable) keeps only the calculations of those methods:

```bash
biotsr analyze common-keys --from-results results/ --method carbohydrate --output carb_out/
```

The calculations not selected are still checked for completeness.

**Nothing is computed**, so no structure is re-read, no `--pdb-dir` is needed,
and the flags that would compute are refused rather than ignored: `--chain`,
the modifiers, `--pdb-dir`, `--column`, `--keep-column`,
`--ignore-column`, `--sheet`, `--map`, `--map-file`, `--participant`,
`--cross-side`, `--inapplicable`, the download options, `--jobs`, `--if-exists` and
`--on-error`.

**Each sample is reported under its prefix** — for a per-residue sweep, under
each residue's own prefix.  If two calculations in different folders of the
tree have one prefix, the command refuses and names both files, because a
report with two rows under one id cannot be matched back to what produced them.

**Compute once, analyse many times** is the point.  The reports are
byte-identical to the recomputing path, which the validation suite checks.

Each sample keeps the key space and settings it was computed under, so the same
checks apply as for a freshly computed set — loading launders nothing.
`--manifest` and `--from-results` are mutually exclusive: one computes the
dataset, the other reads it.

### What is still Python-only

Feature matrices, similarity measures, hierarchical clustering and the clustered
heat map are implemented in `biotsr.analysis` and used from Python — over
computed *and* loaded results alike.  The default hierarchical analysis needs
no cluster count: average linkage, dendrograms and the ordered heat map; a flat
partition is optional and needs an explicit rule (a number of clusters or a cut
height).  BioTSR 1.0 exposes common keys on the command line; for the rest,
[tutorial 08](tutorials/08_similarity_clustering.md) goes from a batch's saved
results to exported matrices, the tree and the figure, and
[`API_GUIDE.md`](API_GUIDE.md) §7 lists the functions.  The figure needs
`pip install "biotsr[plot]"` (matplotlib).  `--jobs` parallelises the structure
computations of a batch, not these analyses.

---

## 9. Reading back what you produced

A folder lists its calculations, one line each — prefix, method, triangles and
whether the calculation is complete (and keys only):

```bash
biotsr describe results/
```

```
1UBQ_A_amino_acid_atomic_LEU_15                          amino_acid_atomic                  56 triangles  complete
1UBQ_A_amino_acid_atomic_perresidue-LEU                  amino_acid_atomic                 504 triangles  complete  (keys only)
1UBQ_A_protein                                           protein                         70300 triangles  complete
```

A prefix, or its `PREFIX.run.json`, shows the completion record — how the
calculation was produced, including the source rows of a batch request:

```bash
biotsr describe results/1UBQ_A_protein
```

```
command          biotsr run
biotsr version   1.0.0rc7
method           protein
prefix           1UBQ_A_protein
complete         True
profile          standard_protein_calpha_v1
parameters
  chain            A
  method           protein
  pdb              1UBQ.pdb
results          1
  1UBQ_A_protein           70300 triangles, 45349 keys  keyspace da39d73510c5  sampling 866885436e9a
```

Point it at the `.meta.json` sidecar (or the `.keys` file) instead for the full
key space and sampling profile — what the numbers mean:

```bash
biotsr describe results/1UBQ_A_protein.meta.json
```

`--json` prints the file itself, formatted, for piping into `jq` (for a folder,
the list of calculations).

`describe` prints; it does not load.  To get the result itself back as an object,
use the Python API, with the `.keys` or the `.meta.json` path:

```python
import biotsr
loaded = biotsr.load_result("results/1UBQ_A_protein.keys")
```

Every result this version writes can be loaded this way.  It restores every
key, both identity fingerprints, the profile that produced the result and the
SSE annotation provenance.  Triplet records and first-seen key order are
deliberately not restored, and the object says so.  Results written by
1.0.0rc3 or earlier are refused, with a message saying they must be
recomputed.  See [`API_GUIDE.md`](API_GUIDE.md) §4.

---

## 10. Errors and exit codes

| code | meaning |
|---|---|
| `0` | success (a batch request with status `no_triangles` whose selection lacked nothing it named, or `inapplicable` under `--inapplicable skip`, is a success) |
| `1` | unexpected internal error — a defect in `biotsr`, not in the command |
| `2` | user, input or configuration error, including invalid download settings; also `--on-error stop` when a request could not even be planned |
| `3` | a batch with requests that failed, were left `failed` or `not_started` by a worker-pool failure, or never started (see `batch_summary.csv` and `batch_failures.csv`), including `--on-error stop` after a computed request fails; or `biotsr fetch` with an id that was not downloaded, refreshed or cached, `refresh_failed` included (see `fetch_report.csv`) |
| `4` | completed and nothing failed, but at least one result is partial (§4, *Empty and partial selections*) or a batch request had no chain: `biotsr run` with a partial result, `biotsr batch` with a `partial` or `chain_not_specified` request — or a `no_triangles` one whose selection also lacked something it named — and none failed, either of them when `--if-exists skip` kept a partial result, and `biotsr analyze` when it analysed partial samples (only with `--allow-partial`).  The results are written and usable; the code lets a script tell them from complete success.  `3` takes precedence: a batch with a failed request exits `3` even if another is partial |
| `130` | interrupted by Ctrl-C or SIGTERM, for every command; a batch interrupted at any point after planning writes its summary first, and `biotsr fetch` its report |

An analysis that refuses a set — partial samples without `--allow-partial`,
or differing processing without `--allow-mixed-processing` — exits `2`.
`--dry-run` computes nothing, and its exit status is unaffected by partial
results.  Before exit code `4` existed, a partial result exited `0`; a script
that treats every non-zero status as a failure (`set -e`, a scheduler's job
state) now stops on a partial result, so accept `4` explicitly where a partial
result is acceptable.

Ordinary mistakes print one line and no traceback:

```
$ biotsr run --method protein --pdb 1UBQ.pdb --chain Z --output results/
biotsr: error: chain 'Z' is not in 1UBQ, so the selection produced no atoms; chains in this structure: A. ...

$ biotsr run --method drug --pdb 3PTB.pdb --chain A --output results/
biotsr: error: method 'drug' needs --entity-name
  e.g. --entity-name ESW --entity-id 404 (omit --entity-id for every copy on the chain)

$ biotsr run --method drug --pdb 3PTB.pdb --chain A --entity-name BEN \
      --entity-id 1 --sse --output results/
biotsr: error: method 'drug' does not support 'sse'
  it supports: --mirror, --size-filter

$ biotsr run --method drug --pdb 3PTB.pdb --chain A --entity-name BEN \
      --entity-id 7 --output results/
biotsr: error: entity BEN 7 is not in chain A of 3PTB: copies of BEN there have ids 1
```

Use `--debug` to get a traceback when you suspect a defect rather than a typo.

---

## 11. Verbosity

Default: one summary per command; in a batch, one line per *failed* request
and none per successful one.  Warnings about the request itself —
`requested but not found: …`, `requested but excluded by the method's
completeness rule: …`, `not supported: …`, `not represented: …` and `no
triangles: …` — and the note naming nucleotides whose atoms differ from the
standard set (`atom sets that differ from the standard nucleotide …`) are
printed by default, by `biotsr run` and `biotsr batch`
alike (in a batch, with the manifest line), and recorded in the completion
record.

* `--quiet` / `-q` — errors only (and, on stdout, the completion record or
  output folder).
* `--verbose` / `-v` — one line per calculation, plus the package's own
  warnings (missing atoms, alternate locations excluded,
  residues the completeness rule declined, degenerate triangles).

Results go to stdout; everything else goes to stderr.

---

## 12. Scripts and HPC

The CLI is built for non-interactive use: it never prompts, never reads a config
file, writes results only under `--output` (downloads go to `--pdb-dir`, and
`biotsr fetch --report` writes where it is told), has no GUI dependency, and
returns a deterministic exit status.

**One job, several CPUs.**  The simplest way to run a whole dataset on a
cluster is one job that asks for several CPUs and lets `--jobs auto` use exactly
those:

```bash
#!/bin/bash
#SBATCH --job-name=biotsr
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/tsr-%j.out

set -euo pipefail

python -m biotsr batch \
    --manifest samples.csv \
    --pdb-dir  /scratch/$USER/pdb \
    --output   /scratch/$USER/tsr-results \
    --jobs auto \
    --if-exists skip \
    --quiet
```

`--jobs auto` reads the allocation (here `SLURM_CPUS_PER_TASK=16`), not the
node's total.  If the job hits its time limit, the scheduler's SIGTERM stops the
batch cleanly with its summary written (exit 130), whatever phase it is in; if a
worker is killed for exceeding `--mem`, the summary is still written and the
batch exits 3.  Either way, submitting the same script again resumes it,
keeping every calculation `--if-exists skip` verifies and recomputing any it
finds incomplete.  Only SIGKILL (for example when a scheduler's grace period
runs out) cannot be caught; its leftover staging directories are removed by the
next run on the same node, or reported.  Download structures first
with `biotsr fetch` on a node that has network access; compute nodes often do
not.  A batch with a `partial` request and no failure exits 4, which a
scheduler records as a failed job; read `batch_summary.csv` to see which
results are partial and what each lacks.

**A job array, one task per row.**  When the queue favours many small jobs:

```bash
#!/bin/bash
#SBATCH --job-name=biotsr
#SBATCH --array=2-501          # manifest line numbers; line 1 is the header
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/tsr-%A_%a.out

set -euo pipefail

MANIFEST=samples.csv
OUTPUT=/scratch/$USER/tsr-results/task_${SLURM_ARRAY_TASK_ID}

# One row of the manifest, with the header, as its own one-row manifest.
ROW=$(mktemp)
trap 'rm -f "$ROW"' EXIT
head -1 "$MANIFEST"                       >  "$ROW"
sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST" >> "$ROW"

python -m biotsr batch \
    --manifest "$ROW" \
    --pdb-dir  /scratch/$USER/pdb \
    --output   "$OUTPUT" \
    --if-exists skip \
    --quiet

# 0 = done, 3 = this row failed and said why in batch_failures.csv,
# 4 = done, but the result is partial (its batch_summary.csv message says what it lacks)
exit $?
```

Each task writes into its own output folder, so the tasks' summaries do not
overwrite each other; `biotsr analyze common-keys --from-results` searches the
parent folder recursively for completion records.  `--if-exists skip` makes the array re-runnable: a
task whose verified output is already there costs nothing.  `python -m biotsr`
rather than `biotsr` because a compute node's `PATH` may not include the
environment's `bin`.

### Overhead

Measured for 1.0 (macOS arm64, Python 3.13):

| | |
|---|---|
| `import biotsr` in a fresh process | ~77 ms, paid once per process |
| orchestration around one calculation | 3–358 ms, 0.6–4.3 % of a run with real computation in it |

The CLI computes nothing, so the only thing that scales with dataset size is the
TSR computation itself.

---

## 13. What the CLI is, precisely

A thin orchestration layer over `biotsr.api` and `biotsr.io`.  It contains no
geometry, no ordering, no binning, no key formula, no contact selection, no SSE
category logic and **no scientific default of its own** — `--mirror` defaults to
*unset*, not to `False`, so what you did not ask for is never passed and the
package's own default applies.

Two gates hold that in place:

* the validation suite runs every representative workflow
  through the CLI and through the Python API and requires the files to be
  byte-identical and the fingerprints, counts and stratification to match;
* structural tests refuse a numeric constant, a serializer filename or an import
  of a scientific internal anywhere in `biotsr/cli/`.

If you need something the CLI cannot express, the Python API has it, and the
mapping is one-to-one: every long flag is the `analyze()` keyword with
underscores written as dashes.
