# BioTSR architecture

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology. This page describes how the package is built.
It covers what exists and why the parts are separated the way they are. The
reasoning behind each decision is in the
[architecture decisions](architecture/DECISIONS.md), and the scientific
definitions are in the [methods pages](methods/README.md).

## 1. One pipeline, many methods

Every TSR method is the same pipeline with different settings. A structure goes
in, and a `TSRResult` comes out:

```text
Structure I/O          structure/        read the PDB file into a domain model
  ↓
Selection              selection/, cross/  which atoms or residues take part, admission
  ↓
Candidate strategy     candidates/, cross/ which index triples are proposed: all C(n,3)
  ↓
Labels                 labels/           point → integer code        ← aa_grouping
  ↓
Canonical ordering     geometry/ordering.py  which vertex is which, from labels and coordinates
  ↓
Geometry               geometry/         edge lengths, angle, longest edge
  ↓                                                                  ← size_filter
Binning                binning/          angle class and distance class
  ↓
Encoding               encoding/         (labels, classes) → one integer key
  ↓                                                                  ← mirror
Annotation             sse/              which category a key's count lands in ← sse
  ↓
TSRResult              domain/result.py  keys, per-category counts, identities, provenance
  ↓
Serialization / analysis   io/, analysis/
```

The modifiers attach at fixed, named points, and data dependencies force that
order:

* amino-acid grouping swaps the label scheme, so it acts before ordering;
* size filtering needs the longest edge, so it acts after geometry;
* mirror needs the pre-fold angle, so it signs a key that is already complete;
* SSE annotation comes after the key and cannot change it.

## 2. The TSR Engine and the method registry

**The TSR Engine** (`methods/engine.py`) is the one loop that wires the stages
together. It contains no science of its own: no ordering rule, bin table, key
arithmetic or modifier branch. Each scientific decision lives in the component
named for it. Modifiers are resolved into the profile before the loop starts,
so all combinations of grouping, mirror and size filtering run one code path.

**The method registry** (`methods/registry.py`, `methods/builtin/`) is a small,
closed set of named, declarative **recipes**. A recipe picks the family's
selector, label scheme, bin schemes, candidate strategy and supported
modifiers, and contains no arithmetic. `biotsr methods` and the
[capability matrix](CAPABILITY_MATRIX.md) are generated from it, so neither can
claim a capability the package lacks.

**Cross TSR** (`cross/`) adds only what is genuinely two-sided: side
definitions, contact pools keyed on full atom identity, and the mixed-side
candidate strategy. Geometry, ordering, binning and encoding are the shared
components.

## 3. One standard profile per method

A **profile** fixes every scientific convention a method runs under
(`profiles/`): the angle, square-root and median conventions, the label and bin
schemes, the admission and completeness policies, and the modifiers in effect.

Each of the 20 methods runs under exactly one profile, its **standard** profile
(`biotsr methods NAME` shows it). Modifiers derive a profile from it before the
engine runs. A caller cannot choose another profile: `profile=` in Python
accepts only `None` or the method's own default profile, and the command line
has no profile option (ADR-087). Arguments are checked per method
(`methods/arguments.py`): one the method does not read is refused by name.

## 4. Two identities: `KeySpace` and `SamplingProfile`

Every result carries two fingerprinted identities (`domain/keyspace.py`):

| identity | answers | contains, for example |
|---|---|---|
| **`KeySpace`** | can these integers be compared at all? | label scheme, angle and distance schemes, angle and square-root conventions, key formula, mirror convention, numerical precision |
| **`SamplingProfile`** | were these counts produced the same way? | selection, scope, completeness policy, candidate strategy, size filter, contact cutoff and filter |

Each fingerprint hashes every scheme's name, content fingerprint and radix, and
every convention's name. Those names are not themselves semantic: the choice is
conservative, so a renamed but otherwise identical scheme would be refused,
never wrongly accepted. Two bin schemes that share a name but differ in one
overflow rule therefore get different fingerprints — and renaming a scheme or a
convention also changes the fingerprint, so those identifiers are never renamed
casually. The two identities are
separate objects and are never merged into one version string.

**Structure-interpretation provenance** (`domain/interpretation.py`) records
which reading of the structure file produced a result. Two readings of the same
file can name an atom differently while both identities stay equal. It is
recorded outside both identities and checked when results are combined
(ADR-079).

## 5. Comparability guards

The analysis layer reduces every result to one `AnalysisSample`: an id, a
key → count mapping and the two identities (`analysis/samples.py`).

* A shared **`KeySpace` is a precondition** of every operation that touches key
  coordinates. Mixing, for example, a mirror-signed and an unsigned result, or
  a Cα and an atomic result, raises `IncompatibleKeySpace` naming the
  difference.
* **Key definitions must match, always.** The decision is taken from the
  `KeySpace` itself, never from method names: different methods that share a
  key space (`aa_ind`, `drug` and `carbohydrate` all label atoms by element at
  radix 112 in 29 × 18 bins; `protein` and `aa_selection` share the Cα
  key space) pass this check, and the same method under `mirror` or
  `aa_grouping` does not. No option relaxes it.
* **Processing settings are checked separately**
  (`analysis/samples.py::check_comparable`, ADR-094). The settings that change
  how triangles are drawn from whatever was selected (`PROCESSING_FIELDS`:
  candidate strategy, completeness and composition rules, size filter,
  degenerate-triangle and sequence-label policies, contact cutoff, filter, rule
  and projection, pool identity and order) and the alternate-location and
  modified-residue rules the result recorded (`source.altloc_policy`,
  `source.modified_residue_policy`; ADR-105) form a processing fingerprint. A
  modified-residue rule that a result written before it existed does not
  record is reported as *not recorded* and never equals the current rule. A
  difference is refused with
  `IncompatibleProcessing` (a subclass of `IncompatibleSampling`) unless the
  caller passes `allow_mixed_processing=True` (`--allow-mixed-processing`);
  then a warning names each differing setting and every sample's settings are
  recorded. Methods that share a key space can still differ here: the
  amino-acid and carbohydrate methods apply their own atom-count completeness
  rules and `drug` applies none, so `aa_ind` against `drug` needs the opt-in.
* **Partial samples are acknowledged separately** (ADR-106). A sample whose
  selection did not get everything it asked for — read from the selection
  record's `not_found`, `unsupported`, `declined` and `atom_subset.missing`
  (`domain/selection_record.py::selection_omissions`), or, for a result written
  before those fields, from its warnings — is refused with `PartialSamples` (a
  `TSRError`) unless the caller passes `allow_partial=True`
  (`--allow-partial`); then nothing is dropped, the report records each
  sample's omissions, and a warning names them. It is independent of
  `allow_mixed_processing`, and neither relaxes the key-space check. A
  complete sample of a different selection needs no acknowledgement.
* **What was selected never blocks a comparison**: structures, chains,
  amino-acid types, ligands, carbohydrates, nucleotides, residue numbers, atom
  subsets, one copy or all copies, and the method name. These are recorded per
  sample — the `TARGET_FIELDS` of the `SamplingProfile` and the selection
  record (`domain/selection_record.py`), which is provenance only — and an
  empty key intersection is a valid result. The full `SamplingProfile`
  fingerprint is still recorded, but it hashes chain names and free-text side
  descriptions, so it is not used to decide comparability.
* Category-stratified SSE analysis has its own guard,
  `IncompatibleSSEStratification`, because the β-strand identity convention is
  recorded in neither fingerprint.

## 6. Package resources

The scientific constants ship inside the package and are versioned with it:

* the element label table (`labels/resources/element_labels_v1.csv`,
  `element,label`, validated when loaded, with `importlib.resources`, never
  from the working directory; every atomic method labels an atom by its element
  field through it, ADR-113);
* the amino-acid label tables, bin schemes, key formula and SSE categories;
* the built-in profiles.

No normal workflow takes any of these as an argument.

## 7. The serialization boundary

`TSRResult` is the source of truth; text files are an export format.

* `io/native.py` writes the native format: the files of one calculation share
  one prefix — `PREFIX.keys`, an optional `PREFIX.triplets`, and a
  `PREFIX.meta.json` sidecar carrying both identities, the selection record,
  the source SHA-256 and the candidate ledger. The bin counts are not in the
  file name; the sidecar records the full bin definitions. The prefix grammar
  (`structure_chain[_chain2]_method[_target][_modifiers]`) is in `io/naming.py`
  and is shared by the command line and `native.default_prefix(result)`. The
  sidecar declares a `result_schema_version` (2) that is independent of the
  package version.
* `io/native_read.py` reads native output back (`biotsr.load_result`, from
  `PREFIX.keys` or `PREFIX.meta.json`). It recomputes the fingerprints and
  counts and checks that the files belong together. It defaults nothing, and it
  refuses a schema version it does not know, including the files written by
  1.0.0rc3 and earlier release candidates. Per-triangle records and first-seen
  key order are deliberately not restored (ADR-075).
* The native format is the only output format; every file BioTSR writes can be
  read back.

## 8. The analysis layer

`biotsr.analysis` provides feature matrices, the four common-key operations,
five similarity measures, and average-linkage clustering that follows SciPy's
algorithms (the scope of the agreement: `methods/ANALYSIS.md`, §5.7). The
default clustering result is the tree and its leaf order (`hierarchy`), which
needs no number of clusters; a flat partition is an optional step with an
explicit rule (`cut_hierarchy`, ADR-099). It
consumes results only; the core never imports it. The numerical code uses only
the standard library. Only the optional figure (`analysis/render.py`) imports a
third-party library, matplotlib, through the `plot` extra.

## 9. The command line

`biotsr.cli` is a thin orchestration layer over `biotsr.api` and `biotsr.io`.
It holds no geometry, binning, key arithmetic or scientific default. A flag you
do not pass is not passed, so the package's own default applies. An
equivalence test requires CLI and Python runs to produce byte-identical files.
The CLI does own the operational contract around a calculation: manifests and
wide dataset sheets read method by method (ADR-070, ADR-093, ADR-097), the one
shared output folder with a prefix and a `PREFIX.run.json` completion record
per calculation (ADR-098), batch reporting, interruption and worker failures
(ADR-072, ADR-092), and verified downloads (`retrieval.py`). The request a
completion record holds includes the processing rules that admit atoms
(`selection/base.py::processing_identity`: the alternate-location and
modified-residue rules, the latter per label domain), so `--if-exists skip`
keeps a calculation only when it was computed under the current rules and in
the key space this version computes for it, never on the package version alone
(ADR-105, ADR-107). A partial result is written and usable but exits `4`
(`cli/errors.py::EXIT_INCOMPLETE`), also when it has no triangles and when a
resumed run keeps one (`BatchReport.n_incomplete` counts all three); a batch
with a failure exits `3`, which takes precedence (ADR-106).

## 10. Package layout

| module | responsibility |
|---|---|
| `api.py` | the public entry points: `analyze`, `analyze_cross`, `analyze_entities`, `load_result`, `list_*` |
| `domain/` | value objects and identities; the dependency sink, importing nothing else in the package |
| `structure/` | the structure model and the one PDB reader, `pdb_legacy` |
| `selection/` | selectors, completeness policies, model/`TER` truncation rules |
| `labels/` | label schemes: the residue labels and the packaged element label table |
| `candidates/` | candidate strategies (exhaustive C(n,3)) |
| `geometry/` | canonical ordering, triangle geometry, conventions |
| `binning/` | angle and distance bin schemes |
| `encoding/` | the key encoder |
| `modifiers/` | amino-acid grouping, mirror, size filter, and the hook vocabulary |
| `sse/` | secondary-structure assignment and the 18 categories |
| `cross/` | two-sided selection, contact pools, mixed-side candidates |
| `methods/` | the TSR Engine, the registry, the built-in recipes and the per-method argument table |
| `profiles/` | the standard profiles |
| `io/` | native write and read, output prefixes (`naming.py`), and the dataset-sheet column vocabulary |
| `analysis/` | features, common keys, similarity, clustering, the optional figure |
| `cli/` | the `biotsr` command, manifests and wide sheets, the output folder and completion records, the batch engine and the parallel backend |
| `retrieval.py` | opt-in, verified download of PDB-format structures (`fetch_structures`) — the only network access |

The dependency direction points one way: `cli` → `api` → `methods` → the
scientific components → `domain`, and `analysis` → `io`/`domain`. Structural
tests enforce the boundaries that matter most: `domain` imports nothing above
it, `analysis` and `cli` import no scientific internals, and SSE adds no
geometry, binning or encoding.

## 11. Where to read more

* [Architecture decisions](architecture/DECISIONS.md): ADR-001 to ADR-100,
  the reason for each choice above.
* [Methods](methods/README.md): the scientific definitions the standard
  profiles implement.
* [Reproducibility note](REPRODUCIBILITY_NOTE.md): which command runs which
  method, and what a saved result records.
* [Platform validation](PLATFORM_VALIDATION.md): what is identical on every
  platform, and what may differ in the final printed digits.
