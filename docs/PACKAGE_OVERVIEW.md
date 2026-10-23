# BioTSR — package overview

**BioTSR is the unified software implementation of the Triangular Spatial Relationship (TSR) methodology.**  It is developed by the **TSR Research Group** (<https://tsrmethod.github.io/>).

This page describes BioTSR **1.0.0rc5**, a release candidate of 1.0.0.

Its shared scientific core is the **TSR Engine**: one geometry, encoding, labelling and candidate layer that every method runs on.

**What the package contains.** The shared scientific core; ordinary protein Cα
TSR with its three composable modifiers; the shared **atomic** foundation and
the non-cross atomic families on it — amino-acid atomic, the five amino-acid
selection modes of the research group's dataset sheets, drug, carbohydrate and
nucleotide/DNA/RNA; and, on the same core, **five cross methods**:
drug–protein, carbohydrate–protein, nucleotide/DNA/RNA–protein, protein–protein
(all-atom) and protein–protein (Cα). Twenty methods, one TSR Engine — plus one
analysis layer, `biotsr.analysis`, that serves all of them, and
**secondary-structure annotation** in `biotsr.sse`.

Each method runs under **one standard profile**, which fixes every scientific
convention it uses; `biotsr methods NAME` shows it. There is no profile to
choose and no alternative output format (ADR-087).

The **`biotsr` command line** is part of the package — a thin orchestration layer
over the same Python API, adding no science, so a dataset can be processed
without writing Python. See **[CLI_GUIDE.md](CLI_GUIDE.md)**.

A result can be **saved and loaded back**. `biotsr.io.save` writes the native
format; `biotsr.load_result` reads it, restoring both scientific identities with
their exact fingerprints, the profile that produced the result and the SSE
annotation provenance — so `biotsr.analysis` runs over previously computed results
without re-reading a structure file. What a loaded result deliberately does not
carry, triplet records and first-seen key order, is stated on the object
(ADR-075). For ordinary use, start at **[API_GUIDE.md](API_GUIDE.md)**.

Candidate-strategy extensions (for example CA-kNN) are not part of this release,
and are deliberately **absent rather than stubbed**.

SSE adds a residue-level annotation and an 18-category classifier over individual
helices and β-strands. It adds **no** geometry, ordering, binning, encoding, key
arithmetic or PDB parsing of its own, and it changes no key: `analyze(..., sse=True)`
and the same call without it produce identical keys, identical counts and
identical `KeySpace` *and* `SamplingProfile` fingerprints. Structural tests
enforce each of those absences.

Cross adds side definitions, contact-pool construction and a mixed-side candidate
strategy. It likewise adds **no** geometry, ordering, binning, encoding or
serialization of its own, and architectural tests enforce that.

The core imports nothing outside the Python standard library, and neither does
the analysis layer. `biotsr.analysis.render` is the single exception — it draws
the clustered heat map with matplotlib (the `plot` extra), and nothing else in
the package imports it.

## Install

```bash
pip install .                                     # from a clone of https://github.com/TSRMethod/BioTSR
pip install dist/biotsr-1.0.0rc5-py3-none-any.whl  # or from a built wheel
pip install biotsr                                # from PyPI, once published
```

One third-party dependency, joblib, for parallel batches (ADR-086); the
calculation itself needs no pandas, numpy or Biopython, so it installs in seconds
and runs on a compute node with no scientific stack. Add `-e` for an editable
development install. Linux and macOS are supported; Windows is not (it is
neither tested nor validated). Input is PDB format only.

## First output, from the command line

```bash
biotsr run --method protein --pdb tests/data/structures/1UBQ.pdb --chain A --output results/
```

```
method     protein
profile    standard_protein_calpha_v1
triangles  70300
keys       45349 distinct
output     results/1UBQ_A_protein (4 files, completion record 1UBQ_A_protein.run.json)
```

`results/` now holds four files sharing the prefix `1UBQ_A_protein`: the key
frequencies (`.keys`), the per-triangle file (`.triplets`), a self-describing
`.meta.json` sidecar carrying both identity fingerprints, and the completion
record `.run.json`, which states the exact command and every file it wrote.
Every later calculation with the same `--output` writes its own prefix into the
same folder.

```bash
biotsr methods                  # the 20 methods and their modifiers
biotsr methods protein          # one in full: arguments, modifiers, standard profile
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/
```

`python -m biotsr ...` is the same program, for an environment whose `bin`
directory is not on `PATH`.

**[CLI_GUIDE.md](CLI_GUIDE.md)** is the full command-line guide: every
workflow, the dataset manifest, the output contract, exit codes and an HPC
example.

## First output, from Python

The examples below run from the repository root, against the structures under
`tests/data/structures/`.

```python
import biotsr

path = "tests/data/structures/1UBQ.pdb"
result = biotsr.analyze(path, method="protein", chain="A")

result.n_triangles          # 70300  = C(76,3)
result.n_unique_keys        # 45349
result.frequencies[5059832] # 1, the frequency of one key
result.key_space.fingerprint
result.counts.n_degenerate  # 0
```

Every command-line flag is the `analyze()` keyword with underscores written as
dashes (`size_filter=25.0` ⇄ `--size-filter 25.0`), so the two surfaces are one
vocabulary. An argument the method does not read is refused by name, in Python
and on the command line; `biotsr methods NAME` lists the ones it does read.

### Modifiers

Three composable keyword modifiers. All eight combinations run on the same
engine — there is no `GroupedMirrorTSR` class, and there is not meant to be:

```python
biotsr.analyze(path, method="protein", chain="A", aa_grouping=True)   # 23 grouped codes, radix 30
biotsr.analyze(path, method="protein", chain="A", mirror=True)        # negate when theta_raw > 90
biotsr.analyze(path, method="protein", chain="A", size_filter=25.0)   # keep max_distance < 25 Å
biotsr.analyze(path, method="protein", chain="A",
               aa_grouping=True, mirror=True, size_filter=25.0)       # compose freely
```

| Modifier | Hook | Changes key meaning? | Changes which triangles? |
|---|---|---|---|
| `aa_grouping` | labels, before ordering | **yes** (scheme + radix) | no |
| `size_filter` | after geometry, before binning | no | **yes** |
| `mirror` | after encoding | **yes** (sign) | no |
| `sse` | annotation, **after** the signed key | **no** | **no** |

`aa_grouping` is **Cα only** and is rejected — naming the required label
domain — anywhere else. `mirror` and `size_filter` read only `theta_raw` and
`max_distance` respectively, so every family reuses them unchanged.

### Secondary structure

`sse=True` annotates each triangle with one of 18 helix/strand/coil categories,
read from the structure's own `HELIX`/`SHEET` records. It is the only entry in
the table above on **neither** identity — it changes no key:

```python
plain     = biotsr.analyze(path, method="protein", chain="A")
annotated = biotsr.analyze(path, method="protein", chain="A", sse=True)

annotated.frequencies           == plain.frequencies            # True
annotated.key_space.fingerprint == plain.key_space.fingerprint  # True
annotated.sampling.fingerprint  == plain.sampling.fingerprint   # True

annotated.category_totals()     # {'7_3b1': 196, '18_1a1b1c': 14256, ...}
sum(annotated.category_totals().values()) == annotated.n_triangles   # always
annotated.triplets[0].sse.category                                   # '7_3b1'
annotated.triplets[0].sse.groups                # ('BET/2', 'BET/2', 'BET/2')
annotated.annotations["sse"]["strand_identity"] # 'sheet_id_and_serial'
```

The categories distinguish **individual β-strands** (a β-sheet holds one or
more), so a strand is named by its `SHEET` record's `(sheetID, serNum)` pair.
`serNum` alone restarts in every sheet; the historical `SSE_TSR.py` used it
alone and so merged strands of different sheets that share a serial number
(ADR-063). The convention is recorded in `annotations["sse"]["strand_identity"]`
and in the metadata sidecar, beside the category totals it qualifies.

SSE composes with all three modifiers, on `protein` and `aa_selection`, and is
refused by name everywhere the categories are undefined — atomic, ligand,
carbohydrate, nucleotide and cross. Under grouping the label space collapses from
20 to 13 while each residue keeps its *own* secondary structure; under mirror the
sign still follows `theta_raw` and never the category.

If a structure carries no `HELIX`/`SHEET` records, every triangle is `17_3c` —
and the result says so, with a warning and `annotations["sse"]["has_records"]`,
rather than letting a missing input read as a finding.

### Size filtering

**`size_filter` is optional and OFF by default for every method.** When enabled,
a triangle is kept iff `max_distance < cutoff`. That rule was established by
executing the historical Cα program, the only historical TSR program with a size
filter; for every other family it is a stated standard policy (ADR-024,
ADR-035).

The candidate ledger always balances, and a filtered triangle is never reported
as a degenerate one:

```python
r = biotsr.analyze(path, method="protein", chain="A", size_filter=25.0)
r.counts.n_candidates        # 70300
r.counts.n_valid_geometry    # 70300
r.counts.n_degenerate        #     0   no definable angle
r.counts.n_filtered_by_size  #  9558   valid, but out of range
r.counts.n_triangles         # 60742   encoded
r.counts.ledger_balances     # True
r.category_totals()          # {'regular': 60742}
```

### Amino-acid atomic TSR

A different family on the same engine — different labels, different bins, a
different scope and a completeness rule, all as configuration:

```python
one = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                     residue="15", residue_name="LEU")
one.counts.n_triangles                 # 56  = C(8,3) over LEU's heavy atoms
one.key_space.distance_scheme.radix    # 18, never 35
one.key_space.label_scheme.radix       # 112 (a radix; the alphabet is 109 codes)

from biotsr.api import analyze_entities
every = analyze_entities(path, method="amino_acid_atomic", chain="A")
len(every)                             # 75 of 1UBQ's 76 residues
```

Gly76 is absent because it carries an `OXT` — five atoms against the completeness
table's four. The rule is **exact equality, not a minimum**, and the rejection is
recorded rather than silent:

```python
every = analyze_entities(path, method="amino_acid_atomic", chain="A",
                         include_rejected=True)
[s.detail for s in every["76_GLY"].skipped if s.stage == "completeness"]
# ['76_GLY: GLY: expected exactly 4 atoms, got 5']
```

Atom selection pools the named atoms across residues into one calculation:

```python
sel = biotsr.analyze(path, method="amino_acid_atom_selection", chain="A",
                     atoms="LEU_8_CD1;LEU_43_CD2;ILE_44_CD1;LYS_48_NZ;HIS_68_NE2;VAL_70_CG1")
sel.n_triangles                        # 20  = C(6,3)
```

`mirror` and `size_filter` work here unchanged. `aa_grouping` does not — it is
Cα-only, and is refused by label domain rather than by a per-family list.

### The amino-acid selection modes

```python
biotsr.analyze(path, method="aa_ind",          chain="A", residue_types="TYR")
biotsr.analyze(path, method="aa_ind_mix",      chain="A", residue_types="TYR;PHE")
biotsr.analyze(path, method="aa_ind_specific", chain="A", residues="LYS_6;GLU18")
biotsr.analyze(path, method="aa_selection",    chain="A", residues="10-12;48;70-71")
```

The first three pool every eligible non-hydrogen atom of every matching residue
into one calculation, in the **atomic** 112-label space. `aa_selection` is a **Cα**
method: one Cα per residue, labelled by residue identity at **radix 30**. The
label space follows the family, never the selection strategy.

An explicit selection that matches nothing is refused, naming what was asked
for. When only some requested items exist, the result carries a warning that
begins `requested but not found:`; a selection of fewer than three points gives
a `no triangles:` warning (ADR-089). In a batch, a calculation that lacks
something it asked for has the status `partial`, written like `ok` but never
reported as complete: the command exits 4, and an analysis uses it only when
acknowledged with `--allow-partial` (ADR-106).

Residue identifiers are PDB numbers with an optional insertion code (`52A`),
and `aa_selection` ranges are ascending plain numbers; a range that would span
an insertion-code residue the list does not name is refused (ADR-102). Of the
modified amino acids stored as `HETATM`, the amino-acid atomic methods and the
protein sides of the atomic cross methods read exactly `TPO`, `SEP` and `PTR`,
which have defined atom counts; any other is named as `not supported:`
(ADR-101). The Cα methods label residues: the 20 regular amino acids and eleven
modified residues (TPO, SEP, PTR, NIY, M3L, TRF, CSO, LLP, CME, CSD with codes
of their own, MSE with methionine's), 30 codes at radix 30, read from `HETATM`
records when the CA atom's element is carbon or blank; any other modified
residue, and a supported one without a readable Cα, is named as
`not represented:` with the reason (ADR-107). Cα results computed before 1.0.0rc5
were encoded at radix 20 and must be computed again.

### Drug, carbohydrate and nucleotide TSR

```python
S = "tests/data/structures/"

# One ligand copy, identified by (chain, name, id).
biotsr.analyze(S + "6G9M.pdb", method="drug", chain="A", entity_name="ESW", entity_id="404")
biotsr.analyze(S + "6G9M.pdb", method="drug_atom_selection", chain="A",
               entity_name="ESW", entity_id="404", atoms="C3;C5;O6;N8;C9;C18;CL1")

# Carbohydrate: one monomer, or the whole chain pooled inter-residue.
biotsr.analyze(S + "5WT9.pdb", method="carbohydrate", chain="A", entity_name="NAG", entity_id="1")
biotsr.analyze(S + "5WT9.pdb", method="carbohydrate_chain", chain="A")

# Nucleotide: one residue, or a DNA/RNA strand.
biotsr.analyze(S + "1BNA.pdb", method="nucleotide", chain="A", entity_name="DC", entity_id="3")
biotsr.analyze(S + "1BNA.pdb", method="nucleotide_chain", chain="A")

# What is actually in this chain?
biotsr.list_entity_targets(S + "7LF9.pdb", chain="D", method="carbohydrate_chain")
# (('1', 'NAG'), ('2', 'NAG'), ('3', 'BMA'), ('4', 'MAN'), ('5', 'NAG'), ('6', 'GAL'), ('7', 'MAN'))
```

The chain is always required and only that chain is searched. A name **with**
an `entity_id` selects that one copy; a name **without** one selects every copy
of that name on the chain, in one combined calculation — even when there is only
one copy — and triangles may then span copies (ADR-096). The scope is recorded,
and so is every copy used:

```python
biotsr.analyze(S + "6G9M.pdb", method="drug", chain="A",
               entity_name="ESW").annotations["target"]       # 'ESW_allcopies'
r = biotsr.analyze(S + "5WT9.pdb", method="carbohydrate", chain="A", entity_name="NAG")
r.annotations["selection_request"]["sides"][0]["instances"]   # ['NAG_1', 'NAG_2']
```

An `entity_id` is a PDB residue identifier: the residue number (it may be
negative) and an optional insertion-code letter — `404`, `52A`, `-5`. Leading
zeros are dropped, `52` never matches `52A`, and an ordinal such as `second` is
refused. A named atom subset (`atoms=`) applies to each selected copy, and an
atom a copy lacks is reported per copy. Every copy of one sugar is not the
whole glycan: `carbohydrate` with `entity_name="MAN"` on 7LF9 chain D uses MAN 4
and MAN 7, while `carbohydrate_chain` and `nucleotide_chain` pool every residue
of the chain. All three families
are the same engine with a different admission policy and a different
completeness policy. All use **29 theta classes, 18 distance classes and radix
112**, and atoms labelled by element. `nucleotide_chain` takes each nucleotide
with the atoms it has and names any whose atoms differ from the standard set
(ADR-114). `mirror` and `size_filter` compose here
exactly as everywhere else.

### Cross TSR — two sides, one framework

**Both chains are always named.** A cross is between two user-specified sides, and
`chain_2` is required — for a glycan or a nucleotide the entity rarely sits in the
protein's chain, so there is nothing sensible to default to (ADR-049).

```python
# drug-protein: chain is the drug's chain, chain_2 the protein chain
biotsr.analyze_cross(S + "6G9M.pdb", method="drug_protein_cross", chain="A",
                     chain_2="A", entity_name="ESW", entity_id="404")

# carbohydrate-protein and nucleotide/DNA/RNA-protein
biotsr.analyze_cross(S + "7LF9.pdb", method="carbohydrate_protein_cross", chain="C",
                     chain_2="A", entity_name="NAG", entity_id="1")
biotsr.analyze_cross(S + "1KX5.pdb", method="nucleotide_protein_cross", chain="J",
                     chain_2="D", entity_name="DT", entity_id="50")

# protein-protein, all atom; the _2 arguments narrow the second side
biotsr.analyze_cross(S + "1BRS.pdb", method="protein_protein_cross", chain="A",
                     chain_2="D", residues="LYS_27;ARG_59", residues_2="ASP_35;ASP_39")

# protein-protein, Ca: interface residues by ANY heavy-atom pair within the
# cutoff, each then represented by its Ca -- which itself need not be within it
biotsr.analyze_cross(S + "1BRS.pdb", method="protein_protein_calpha_cross",
                     chain="A", chain_2="D")
```

An atom of either side enters the pool **iff** it has at least one non-hydrogen
atom of the other side within the cutoff. Symmetric, and atom-level: an atom is
not admitted because a different atom of its residue or ligand touched. Hydrogen
is excluded from every stage.

A cross triangle always contains points from **both** sides — `A+A+B` or `A+B+B`,
never `A+A+A` or `B+B+B` — and there are exactly `a·C(b,2) + b·C(a,2)` of them.
They are enumerated directly, not filtered out of `C(a+b,3)`.

The **contact cutoff** decides which atoms enter the pool. It is *not* a
constraint on the resulting triangle's edges: an atom admitted because it sits
3 Å from one atom of the other side may form a triangle with a 20 Å edge. To
constrain triangle size, use `size_filter`, a separate modifier applied after
geometry and off by default:

```python
pp = dict(method="protein_protein_cross", chain="A", chain_2="D")
biotsr.analyze_cross(S + "1BRS.pdb", contact_cutoff=4.5, **pp)  # a SamplingProfile fact:
                                                               # the KeySpace is unchanged
biotsr.analyze_cross(S + "1BRS.pdb", size_filter=12.0, **pp)    # a triangle constraint
```

Pool membership is keyed by the **full physical `AtomIdentity`**, so `LYS 24 NZ`
and `LYS 28 NZ` are two atoms. Two historical cross programs keyed their pools by
`(side, atom_name)` and so collapsed distinct atoms with the same name (ADR-041);
the pool records that nothing was lost:

```python
r = biotsr.analyze_cross(S + "1KX5.pdb", method="nucleotide_protein_cross",
                         chain="J", chain_2="D", entity_name="DT", entity_id="50")
r.annotations["cross_pools"]["pool_identity"]              # 'full_atom_identity'
r.annotations["cross_pools"]["n_atoms_lost_to_identity"]   # 0
```

### Saving and streaming

Write the native output — `PREFIX.keys`, `PREFIX.triplets` and a
`PREFIX.meta.json` sidecar carrying the full `KeySpace` and `SamplingProfile`.
The prefix defaults to `native.default_prefix(result)`, the same grammar the
command line uses (`1UBQ_A_protein` here); pass `prefix=` to choose one:

```python
from biotsr.io import native
native.save(result, "out/")                       # out/1UBQ_A_protein.keys, ...
loaded = biotsr.load_result("out/1UBQ_A_protein.keys")
loaded.key_space.fingerprint == result.key_space.fingerprint   # True
```

For runs producing millions of triangles, stream the triplets instead of
retaining them:

```python
with native.NativeTripletStreamWriter("out/1UBQ_A_protein.triplets") as w:
    streamed = biotsr.analyze(path, method="protein", chain="A",
                              keep_triplets=False, triplet_sink=w.write)
```

## What is implemented

| Layer | Module | Content |
|---|---|---|
| domain | `domain/` | `StructureId`, `ChainId`, `ResidueIdentity`, `AtomIdentity`, `Point`, `PointSet`, `TriangleCandidate`, `KeySpace`, `SamplingProfile`, `TSRResult`, the error hierarchy |
| structure | `structure/` | format-neutral model + the one PDB reader, `pdb_legacy` |
| selection | `selection/` | protein Cα selectors (whole chain and residue selection); one shared entity selector serving every atomic scope; atom-admission, truncation and completeness policies |
| labels | `labels/` | the residue labels (30 codes) and the grouped residue labels (23 codes), both at radix 30, the window rule that keeps every scheme collision-free, and the 112-radix element labels (`element_labels_v1.csv`, a packaged resource): atoms are labelled by their element field, ADR-113 |
| modifiers | `modifiers/` | AA grouping, mirror, size filter — three named hooks, no per-combination classes; plus the `ANNOTATION` hook SSE attaches to |
| sse | `sse/` | the three states, the element identity that names an individual helix or β-strand, the `HELIX`/`SHEET` → residue assignment, the 18 categories as a signature→code bijection, and a reshaping helper that is explicitly not an analysis |
| candidates | `candidates/` | exhaustive `C(n,3)` |
| cross | `cross/` | side specs, contact-pool construction under full `AtomIdentity`, the mixed-side candidate strategy, and the runner that hands them to the shared core |
| geometry | `geometry/` | canonical ordering (four cases, one convention), triangle geometry, the midpoint median, `GeometryResult` |
| binning | `binning/` | 29-class theta; the 35-class Cα distance scheme; the 18-class atomic scheme |
| encoding | `encoding/` | the shared mixed-radix key formula |
| methods | `methods/` | registry + twenty declarative recipes + the per-method argument table + the shared core loop |
| profiles | `profiles/` | the standard profiles — Cα, atomic and cross — all the same `Profile` type, no subclasses |
| io | `io/` | the native writer and reader; output prefixes (`naming`); the dataset-sheet column map (`input_columns`) |
| retrieval | `retrieval.py` | opt-in download of PDB-format entries from RCSB, with transfer, gzip and content verification and a validated cache — the only module that uses the network |

## What is deliberately not implemented

Candidate-strategy extensions (for example CA-kNN), mmCIF input, Windows, and a
structure reader that uses the PDB specification's atom-name columns
throughout. (Downloading PDB-format entries is available on request —
`biotsr fetch`, `biotsr.fetch_structures` — and a batch can use several worker
processes with `--jobs`; see the command-line guide.)

BioTSR does not reproduce the output files of the historical TSR programs. Where
a standard convention differs from a historical program's — the angle divisor,
the top Cα distance class, the grouped radix, the cysteine sulfur label, β-strand
identity, cross pool identity and contact filtering — the difference and its
reason are recorded in the [architecture decisions](architecture/DECISIONS.md)
and the [methods pages](methods/README.md).

SSE is implemented for the two protein Cα methods and **deliberately undefined**
for every other family: an atom inside a residue has no secondary structure of
its own, and a ligand, a glycan and a nucleotide have none at all.

Mirror and SSE are not offered for cross methods: no historical cross program
has either, and neither is extended by analogy. Mixing a Cα side with an atomic
side is refused: no established key formula combines the residue label space
with the 112-code atom space.

Some capabilities have **no historical counterpart**: size filtering for every
atomic family, the drug mirror, the DNA/RNA chain scope, carbohydrate and
nucleotide atom narrowing, and the amino-acid type-selection modes. Their rules
are stated as standard policy, not presented as reproduced history.

Asking for something a method does not support raises an error naming what *is*
supported rather than silently doing nothing:

```python
biotsr.analyze(path, method="protein", chain="A", annotate="sse")
# ValueError: method 'protein' does not support 'annotate'; it supports
#             aa_grouping, mirror, size_filter, sse

biotsr.analyze(S + "3PTB.pdb", method="drug", chain="A", entity_name="BEN",
               entity_id="1", sse=True)
# ValueError: method 'drug' does not support 'sse'; it supports
#             mirror, size_filter
```

## Two identities

`KeySpace` answers "can these integers be compared at all?"; `SamplingProfile`
answers "were these counts produced the same way?" They are separate objects with
separate content fingerprints, and they are not merged into one version string.

The reason is a real hazard, not a hypothetical: the historical `TSR-Package`
Cα distance binner and the standard 35-class scheme share **every** cut point
and both declare 35 classes, yet produce different keys for any triangle at or
beyond 132.45 Å, because only their top class differs. A human-assigned version
string would have been identical for both, so fingerprints hash boundaries,
class count, interval semantics and overflow behaviour together.

## How the standard differs from the historical programs

| | standard (BioTSR) | historical programs |
|---|---|---|
| angle divisor | `math.pi` | the literal `3.14` in most programs |
| Cα distance classes | 35 (`d ≥ 132.45 → 35`) | 34 effective in `TSR-Package` (top two merged) |
| grouped radix | 13 | 20 in `TSR-Package` |
| cysteine `SG` label | 111 (sulfur) | 72 (Seaborgium) |
| β-strand identity | `(sheetID, serNum)` | `serNum` only in `SSE_TSR.py` |
| cross pool identity | full physical atom identity | atom name only in two cross programs |
| cross contact filter | symmetric, atom-level | protein side only in `Drug_TSR.py` |
| grouped triplet residue column | label name *and* true residue | group representative only |
| per-key category counts | one counter per key | shared counter in `SSE_TSR.py` (a defect) |

On 1UBQ chain A the angle divisor is the only difference that applies: 1,144 of
70,300 triangle keys differ from the historical `TSR.py` output (ADR-002).

## Tests

```bash
pytest                     # default suite
pytest -m slow             # the large reference structures only
pytest -m ""               # everything
```

The command-line and API gate scripts and the release checks are part of the
source repository, not of the package (see [RELEASING.md](RELEASING.md)):

```bash
python validation/scripts/gates_cli.py            # command-line gates
python validation/scripts/gates_api_hardening.py  # save/load and API gates
```

The suite covers all four canonical-ordering cases including both tie branches,
the degenerate-geometry cases, bin boundaries, encoding against hand-calculated
keys, fingerprint sensitivity *and* insensitivity, deterministic serialization,
real PDB parsing, every method and all eight Cα modifier combinations, and
standard reference results for 42 calculations over all 20 methods
(`tests/data/reference/standard_baseline.json`), compared by keys-file SHA-256,
counts and both fingerprints.

It also covers the **round trip**: every representative family saved and loaded
back, both fingerprints unchanged, `serialize(deserialize(serialize(r)))`
byte-identical to `serialize(r)`, and every refusal the reader must make — an
unversioned or future schema, a truncated or reordered keys file, a sidecar
whose stated fingerprint does not recompute from its own contents, and a keys
file paired with a different result's sidecar.

## Analysis

`biotsr.analysis` consumes a `TSRResult` from any method family — or any
`key -> count` mapping carrying the two identities — and never asks which engine
produced it.

```python
from biotsr.analysis import (
    AnalysisSample, feature_matrix, distance_matrix, clustermap_spec,
)

# One nucleotide TSR result per residue of the B-DNA dodecamer 1BNA, chain A.
# The 5'-terminal DC 1 lacks a phosphate, so the completeness rule leaves it
# with no points; the other eleven residues are the samples.
results = {}
for entity_id, name in biotsr.list_entity_targets(S + "1BNA.pdb", chain="A",
                                                  method="nucleotide"):
    r = biotsr.analyze(S + "1BNA.pdb", method="nucleotide", chain="A",
                       entity_name=name, entity_id=entity_id)
    if r.n_triangles:
        results["1BNA_A_%s" % entity_id] = (r, name)

# The class of each sample is the base the structure itself names, kept from
# list_entity_targets -- never parsed back out of a sample id.
samples = [
    AnalysisSample.from_result(r, sample_id, label=base)
    for sample_id, (r, base) in results.items()
]

fm = feature_matrix(samples)            # columns ordered by integer key value
dm = distance_matrix(fm)                # 1 - sum(min)/sum(max)

# The default result: the average-linkage tree and its leaf order. No number
# of clusters is needed, and none is reported.
spec = clustermap_spec(dm, [s.label for s in samples])
spec.hierarchy.ordered_sample_ids()     # dendrogram leaf order
spec.assignment                         # None: no partition, so no ARI

# Optional: an explicit cut. Four clusters, because four bases are present, so
# the ARI measures agreement at that granularity; it does not discover that
# there are four groups.
cut = clustermap_spec(dm, [s.label for s in samples], n_clusters=4)
cut.assignment.labels                   # flat clusters (average linkage, maxclust)
cut.adjusted_rand                       # 1.0: the four clusters are the four bases
```

A shared `KeySpace` is a **precondition**, decided from the key definitions and
never from method names: grouped and ungrouped Cα results, or mirror-signed and
unsigned ones, are refused with the differing field named, while `protein` and
`aa_selection` (one Cα key space) compare. A difference in **processing** — a
contact cutoff, size filter, completeness rule or alternate-location rule — is
refused unless `allow_mixed_processing=True`, which records a warning and every
sample's settings. What was selected never blocks a comparison: structures,
chains, residue types (LYS-only against TYR-only), molecules, residue numbers,
atoms, one copy or all copies; an empty key intersection is a valid result
(ADR-094).

Four separately named common-key operations, five similarity measures, average
linkage, the dendrogram leaf order and the optional `maxclust` or cut-height
partition are implemented on the
standard library; how far they agree with SciPy is stated precisely in
[methods/ANALYSIS.md](methods/ANALYSIS.md), §5.7 (ADR-060). Drawing the figure
is optional:

```python
from biotsr.analysis.render import render_clustermap   # needs biotsr[plot]
render_clustermap(spec, "clustermap.png")
```

Without matplotlib, everything above still works — the numbers live on `spec`,
not in the figure. The executed, end-to-end version is
[tutorial 08](tutorials/08_similarity_clustering.md).
