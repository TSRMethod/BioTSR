# BioTSR — Python API guide

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology.

Ordinary use of the package, start to finish, without reading an architecture
document. Seven things: run a method, look at what came back, save it, load it
back, check its provenance, build a feature representation, run an analysis.

Everything on this page is executed by
`tests/integration/test_api_guide.py`, which runs each numbered
example and asserts what it claims. A line that stops being true stops the test
suite.

For the command line, see `CLI_GUIDE.md`. For which command computes which
published method, see `REPRODUCIBILITY_NOTE.md`. For why anything is the way it
is, see `ARCHITECTURE.md` and `architecture/DECISIONS.md` — but not for this.

---

## Install

```bash
pip install .                                       # from a clone of https://github.com/TSRMethod/BioTSR
pip install dist/biotsr-1.0.0rc7-py3-none-any.whl   # or from a built wheel
```

`pip install biotsr` will work once 1.0 is published on PyPI; it is not yet.
Python 3.9 or later, on Linux or macOS. The installation brings one third-party dependency,
joblib, used only to run batch rows in parallel; the calculation itself uses the
standard library alone.

---

## 1. Run a method from Python

```python
import biotsr

result = biotsr.analyze("1UBQ.pdb", method="protein", chain="A")
```

That is the whole call. `method` names one of twenty; `chain` is required and is
never guessed.

Modifiers are keyword arguments on the same call and compose freely — there is no
per-combination class and no per-combination method name:

```python
result = biotsr.analyze(
    "1UBQ.pdb", method="protein", chain="A",
    aa_grouping=True,      # the grouped residue labels (23 codes), radix 30
    mirror=True,           # negate a key whose pre-fold angle exceeds 90 degrees
    size_filter=25.0,      # keep triangles whose longest edge is under 25 A
)
```

Secondary structure is an **annotation**, not a second engine. `sse=True`
produces the *same integer keys*; what changes is which of the eighteen
categories each key's frequency is counted in:

```python
annotated = biotsr.analyze("1TUP.pdb", method="protein", chain="A", sse=True)
plain = biotsr.analyze("1TUP.pdb", method="protein", chain="A")
assert annotated.frequencies == plain.frequencies          # identical keys
assert annotated.key_space.fingerprint == plain.key_space.fingerprint
```

Methods that take a target take it instead of a whole chain:

```python
leu15 = biotsr.analyze("1UBQ.pdb", method="amino_acid_atomic", chain="A",
                    residue="15", residue_name="LEU")
ligand = biotsr.analyze("3PTB.pdb", method="drug", chain="A",
                     entity_name="BEN", entity_id="1")
```

Cross TSR between two sides goes through its own entry point, because its
arguments are routed by role:

```python
cross = biotsr.analyze_cross(
    "1BRS.pdb", method="protein_protein_cross",
    chain="A", chain_2="D", contact_cutoff=5.0,
)
```

`biotsr.list_methods()` returns all twenty; `biotsr.list_cross_methods()` the five
cross families. An argument a method does not take raises `InvalidRequest`,
naming what it does take — never a silent no-op. Every cross method needs
`chain_2`.

**One profile per method.** Each method runs under its own standard profile,
`biotsr.get_method(NAME).default_profile`. `profile=` exists so a script can
state it explicitly: it accepts that profile or `None`, and refuses anything
else — a profile name, another method's profile, a derived profile — with
`InvalidRequest`. `analyze_cross` has no `mode=` argument.

**Molecule instances.** The single-entity methods (`drug`, `carbohydrate`,
`nucleotide`, their `*_atom_selection` variants, and the entity side of
`drug_protein_cross`, `carbohydrate_protein_cross` and
`nucleotide_protein_cross`) always need `chain`, and search only that chain.
`entity_name` with `entity_id` is that one copy (`target` `'BEN_1'`);
`entity_name` alone is **every** copy of that name on the chain, in one
combined calculation, even when there is one copy (`target` `'NAG_allcopies'`):

```python
nag = biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D", entity_name="NAG")
nag.annotations["target"]                                      # 'NAG_allcopies'
nag.annotations["selection_request"]["sides"][0]["instances"]  # ['NAG_1', 'NAG_2', 'NAG_5']
```

`entity_id` is a PDB residue identifier — a residue number, possibly negative,
and an optional insertion-code letter: `"404"`, `"52A"`, `"-5"`. `"52"` never
matches `52A`; an ordinal or malformed id raises `InvalidRequest`. Named atoms
(`atoms=`) apply to each selected copy, and an atom a copy lacks is reported
per copy. `biotsr.list_entity_targets(path, chain=..., method=...)` lists what
a chain holds.

**Residues.** `amino_acid_atomic` with `residue` and `residue_name` is that
residue; `analyze_entities(..., method="amino_acid_atomic", residue_name="LEU")`
runs one calculation per LEU residue. `aa_ind_specific` names residues as
`NAME_ID` or `NAMEID` (`"LYS_6;GLU18"`); a token matching neither raises
`InvalidRequest`.

**Residue identifiers and ranges.** A residue identifier is its PDB number with
an optional insertion code (`"48"`, `"52A"`, `"-5"`; leading zeros dropped) in
`residue`, `entity_id`, `NAME_ID` tokens (`"LYS_52A"`), atom tokens
(`"LYS_52A_NZ"`) and `residues=` of `aa_selection`, whose items may also be
ascending ranges of plain numbers (`"10-12"`, `"-5--3"`). A range with an
insertion code at either end, a descending range, a malformed item, and a range
that spans an insertion-code residue of the structure unless that residue is
also listed, each raise `InvalidRequest` naming the item. The order of a list
never matters. The command line and a manifest use the same parser:

```python
biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A",
               residues="9-11;10A;10B").selection.residues_used
# ('9', '10', '10A', '10B', '11')
biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A", residues="9-11")
# InvalidRequest: ... span residue(s) 10A, 10B, which have an insertion code ...
```

(`SYNPHOS.pdb` is a synthetic structure in `tests/data/synthetic/`.)

**Modified amino acids.** Of the modified amino acids RCSB files write as
`HETATM`, the amino-acid atomic methods and the protein sides of the atomic
cross methods read exactly `TPO`, `SEP` and `PTR` (11, 10 and 16 heavy atoms),
under the same alternate-location and completeness rules; any other modified
residue (`CME`, `MSE`, …) is not supported by them. The Cα methods label
residues and read eleven modified residues — `TPO`, `SEP`, `PTR`, `NIY`, `M3L`,
`TRF`, `CSO`, `LLP`, `CME`, `CSD` with codes 24–33 and `MSE` with methionine's
13, at radix 30 (ADR-107) — when the CA atom's element is carbon or blank; any
other modified residue, and a supported one without a readable Cα, is not
represented, with the reason. Neither is dropped silently: the
result carries a warning beginning `"not supported: "` or
`"not represented: "`, and a request of unsupported residues only raises
`InvalidRequest` with that reason.

```python
calpha = biotsr.analyze("SYNPHOS.pdb", method="protein", chain="A")
calpha.key_space.label_scheme.radix               # 30
calpha.source.modified_residue_policy
# 'hetatm_amino_acids:CME,CSD,CSO,LLP,M3L,MSE,NIY,PTR,SEP,TPO,TRF;ca_element:C_or_blank'
biotsr.api.processing_identity("protein")["modified_residue_policy"]   # the same
```

Cα results computed before 1.0.0rc5 (radix 20, no modified residues) have
another key space, so no analysis combines them with current ones: compute
them again.

**Atomic labels.** Every atomic method labels a heavy atom by its element
field (columns 77–78, any case), never by its atom name (ADR-113). One scheme
object, `biotsr.ELEMENT_LABELS_112`, serves all of them; it is built from the
packaged table `element_labels_v1.csv` (`element,label`: 110 elements, 109
labels, radix 112):

```python
from biotsr import ELEMENT_LABELS_112 as labels

labels.resolve_atom("PA", "P")       # 110: ATP's alpha phosphorus is phosphorus
labels.resolve_atom("CL1", "C")      # 106: a carbon named CL1 is carbon
labels.resolve_atom("CL35", "CL")    # 93: chlorine, whatever its name
labels.resolve_atom("CA", "CA")      # 99: a calcium ion
labels.resolve_atom("X1", "X")       # None: an unknown element has no label
labels.label_source, labels.resolver # ('element', 'element_field_v1')
```

Both enter the scheme's content fingerprint, so atomic results computed before
1.0.0rc7, when atoms were labelled by name (`atom_lexicon_standard_112`), are
in another key space and are never combined with current ones — compute them
again. `biotsr.labels.elements.parse_element_labels(text)` validates a table in
the same format (header, two fields, chemical-case symbols, positive labels,
no duplicate symbol, no shared label but H/D) and raises `ElementTableError`.

**`nucleotide_chain`** takes every standard nucleotide of the chain with the
heavy atoms it has, and names in one warning each nucleotide whose atoms
differ from the standard set (ADR-114):

```python
dna = biotsr.analyze("1BNA.pdb", method="nucleotide_chain", chain="A", keep_triplets=False)
dna.counts.n_points        # 243: DC 1, which lacks its 5' phosphate, takes part
dna.warnings[-1]
# "atom sets that differ from the standard nucleotide (as the structure file has
#  them; every atom present is used): DC_1 lacks P, OP1, OP2"
biotsr.api.processing_identity("nucleotide_chain")["completeness_policy"]
# 'available_nucleotide_atoms'
```

The warning does not make the result partial. `nucleotide` and
`nucleotide_atom_selection` keep the exact atom count.

**Input that is refused rather than computed** (ADR-108, ADR-109). Every one of
these raises `InvalidRequest` naming the records and the remedy, and applies
only to what the calculation selects — a problem in another chain never stops
it:

* a residue named `UNK` in the selection (a whole chain containing one, a named
  UNK residue, type or atom, a cross side): no calculation is done for UNK;
* an atom with a blank element field (columns 77–78) in an atomic calculation
  or among a cross's contact atoms; a `HETATM` Cα may still have one;
* an atom written twice (repeated or conflicting records), or two atoms at
  exactly the same coordinates;
* a record whose residue number cannot be read (hybrid-36);
* an atom whose element field names no element in the label table (`X`), so it
  cannot be labelled (this used to raise `KeyError`);
* a protein–protein cross of a chain with itself, or a cross whose two sides
  share an atom;
* with `sse=True`, a `HELIX`/`SHEET` record of the chain whose endpoints cannot
  be found by number and insertion code.

Atom names are read from the full PDB field (columns 13–16) and elements from
columns 77–78, case-insensitively; only an atom named exactly `CA` in columns
14–15 is a Cα. These input rules are recorded with every result, and are part
of the processing identity that resume and analysis compare:

```python
calpha.source.input_rules                                  # 'input_rules_v2'
biotsr.api.processing_identity("protein")["input_rules"]   # the same
```

A result written by 1.0.0rc5 or earlier records none (`None`, *not recorded*):
it is not reused by `--if-exists skip`, and comparing it with a current result
needs `allow_mixed_processing=True`.

**Selections.** An explicit selection that matches nothing raises
`InvalidRequest`. Requested items that are missing leave a warning on
`result.warnings` starting `"requested but not found: "` (or, for a named
residue the completeness rule declines, `"requested but excluded by the
method's completeness rule: "`), and a selection of fewer than three points
gives a result with no triangles and a warning starting `"no triangles: "`.

Each side's selection record,
`result.annotations["selection_request"]["sides"][i]`, says what was asked for
and what was used: `chain`, `scope`, `arguments`; for the scopes that name
residues or molecules, `instances` (what contributed atoms) and `excluded`
(requested instances that contributed nothing); and, only when non-empty, `not_found` (what the request named that
the structure does not contain), `unsupported` (what it contains but the
method cannot read) and `declined` (a residue or molecule copy you named that
the completeness rule rejected):

```python
named = biotsr.analyze("SYNPHOS.pdb", method="aa_ind_specific", chain="A",
                       residues="TPO_3;CME_8;LYS_99")
side = named.annotations["selection_request"]["sides"][0]
side["instances"], side["not_found"], side["unsupported"]
# (['TPO_3'], ['LYS_99'], ['CME_8'])
```

A result whose record lists any of these, or atoms a selected copy lacks
(`atom_subset.missing`), is **partial**: it did not get everything it asked
for. A residue the completeness rule rejects among *all* residues of a type,
or a copy among *all* copies of a molecule, does not make it partial.
What was left out is recorded per result in the completion record
(`omissions`, `{}` when nothing was) and is read back as
`AnalysisSample.omissions`; for a result written before these fields existed,
it is read from the warnings above, `not supported:` and `not represented:`.
Every analysis refuses a partial sample unless it is acknowledged (§6).

### Very large results

`keep_triplets=False` keeps the keys and the counts and discards the per-triangle
records, which is what dominates memory. Pass `triplet_sink` alongside it to
stream the records straight to a file instead:

```python
from biotsr.io.native import NativeTripletStreamWriter

with NativeTripletStreamWriter("1UBQ_A_protein.triplets") as sink:
    result = biotsr.analyze("1UBQ.pdb", method="protein", chain="A",
                         keep_triplets=False, triplet_sink=sink.write)
```

The result is otherwise identical — same keys, same counts, same provenance.

---

## 2. Inspect a `TSRResult`

```python
result.frequencies          # {key: total frequency}  -- what analysis consumes
result.keys                 # {key: CategoryCounts}   -- the per-category breakdown
result.n_triangles          # triangles that produced a key
result.n_unique_keys
result.category_totals()    # {"regular": 70300}, or the 18 SSE buckets
result.triplets             # per-triangle records, or None
result.warnings             # e.g. "requested but not found: ...", "no triangles: ..."
result.skipped              # every candidate that did not contribute, and why
```

The candidate ledger balances exactly, and nothing is dropped silently:

```python
c = result.counts
assert c.n_candidates == c.n_degenerate + c.n_filtered_by_size + c.n_triangles
assert c.ledger_balances
```

**The two identities.** `KeySpace` answers *can these integers be compared at
all?*; `SamplingProfile` answers *were these counts produced the same way?* They
are separate objects with separate fingerprints and are never merged:

```python
result.key_space.fingerprint     # changes if the keys mean something different
result.sampling.fingerprint      # changes if triangles were drawn differently
```

A size filter changes the second and not the first — the keys still mean what
they meant, they were just sampled differently. Grouping changes the first, twice
over: a different label scheme and a different radix. The sampling fingerprint
also hashes the chain names and side descriptions, so analysis does not use it
to decide comparability (§6).

---

## 3. Save a result

```python
from biotsr.io import native, save

native.default_prefix(result)    # '1UBQ_A_protein'
written = save(result, "results/")
# {'keys': 'results/1UBQ_A_protein.keys',
#  'metadata': 'results/1UBQ_A_protein.meta.json',
#  'triplets': 'results/1UBQ_A_protein.triplets'}
save(result, "results/", prefix="ubiquitin", triplets=False)   # your own prefix
```

Three files, or two with `triplets=False`, side by side in the folder you name
and sharing one prefix. The default prefix is the one the command line gives the
same request — `structure_chain[_chain2]_method[_target][_modifiers]`, for
example `7LF9_D_carbohydrate_NAG_1` or `1UBQ_A_protein_aagroup_mirror_size25`.
`save` replaces files of the same names; it writes no completion record, which
only the command line adds. The keys file is
`key<TAB>category<TAB>frequency`, sorted numerically by key then alphabetically
by category — deterministic, so two runs of the same request produce
byte-identical files. The `.meta.json` sidecar carries everything else: both
identities, the profile, the selection, the source, the ledger, the SSE
provenance.

The bin counts are not in the file name: the sidecar records the full key
space, both bin schemes included, so nothing has to be recovered from a name.

---

## 4. Load the saved result

```python
loaded = biotsr.load_result("results/1UBQ_A_protein.keys")
```

Either file names the pair — the keys file or its sidecar. What comes back is a
real `TSRResult`:

```python
assert loaded.frequencies == result.frequencies
assert loaded.key_space.fingerprint == result.key_space.fingerprint
assert loaded.sampling.fingerprint == result.sampling.fingerprint
assert loaded.counts == result.counts
```

**Two things a loaded result deliberately does not carry**, and it says so rather
than pretending:

```python
loaded.triplets        # None  -- alt_loc is per-atom and is in no column
loaded.key_order       # ()    -- first-seen order is not serialized
loaded.annotations["loaded_from"]
# {'keys': '/abs/path/.../1UBQ_A_protein.keys', 'metadata': '...', 'result_schema_version': 2,
#  'triplets_loaded': False, 'key_order_loaded': False}
```

Nothing is fabricated for either. Re-run the calculation if you need the
per-triangle records.

**The reader refuses rather than guesses.** A sidecar with no
`result_schema_version` (a sidecar from before the native format was versioned), a version from
the future, a file written by 1.0.0rc3 or an earlier release candidate (schema 1,
`*.keys_29_35` names: recompute it), a truncated or reordered keys file, a sidecar whose stated
fingerprint does not recompute from its own contents, or a keys file paired with
a *different* result's sidecar — each raises with the mismatch named:

```python
try:
    biotsr.load_result("results/sample_A.keys")
except biotsr.MalformedResult as exc:
    ...        # the files do not belong together, or one is corrupt
except biotsr.UnsupportedFormat as exc:
    ...        # this build cannot read that schema version
```

---

## 5. Verify provenance

Everything needed to interpret a result rides on the object. Nothing is recovered
by parsing a filename:

```python
loaded.method                                       # 'protein'
loaded.annotations["profile"]                       # 'standard_protein_calpha_v1'
loaded.annotations["base_profile"]
loaded.annotations["target"]                        # 'A'; e.g. 'BEN_1' (one copy), 'NAG_allcopies'
loaded.annotations["modifiers"]                     # ('mirror',) if mirror was on
loaded.source.path
loaded.source.sha256                # of the structure file, as it was read
loaded.source.parser_backend
loaded.selection.chains_used
loaded.key_space.numerical_precision               # 'float64'
```

For an SSE run, the annotation provenance too:

```python
from biotsr.sse import (
    category_totals_in_order, check_conservation,
    require_common_strand_identity, strand_identity_of,
)

strand_identity_of(loaded)              # 'sheet_id_and_serial'
loaded.annotations["sse"]["has_records"]  # did the file assert ANY structure?
dict(category_totals_in_order(loaded))    # all 18, including the zeros
check_conservation(loaded)                # categories partition the triangles exactly
```

`has_records` matters: without it, "every triangle is `17_3c`" reads as a finding
when it may only mean the structure carried no `HELIX`/`SHEET` record.

**Comparing SSE categories across results requires compatible conventions**, and
`KeySpace`/`SamplingProfile` cannot catch it — the strand-identity convention is
in neither identity, so two results that differ only there have identical
fingerprints and genuinely incomparable categories:

```python
require_common_strand_identity([a, b])   # call before ANY category comparison
```

Ordinary key-frequency analysis does not need this and must not use it: mixing
the conventions there is correct, because the categories are summed away.

---

## 6. Build a feature representation

```python
from biotsr.analysis import AnalysisSample, feature_matrix

samples = [AnalysisSample.from_result(r) for r in results]
fm = feature_matrix(samples)            # columns ordered by integer key value

fm.sample_ids
fm.columns                              # the key coordinates, deterministic
fm.index_of("1UBQ_A")                   # row index for a sample id
fm.row(0)                               # the dense vector, materialised on demand
fm.value(0, 5059832)                    # one cell; a key not present is zero
fm.comparison                           # key space, processing policies and who shares each
fm.samples                              # per sample: fingerprints, settings, selection, warnings
fm.to_header_csv()                      # sample_id, then one column per key
```

`AnalysisSample` takes a `TSRResult` from **any** method family, computed or
loaded — the analysis layer never asks which engine produced it. It reduces to an
id, `{key: count}`, and the two identities.

**One comparison rule for every analysis** (`biotsr.analysis.check_comparable`):

1. **Key definitions, always.** A shared `KeySpace` is a precondition and is
   refused rather than coerced: grouped and ungrouped Cα results, or
   mirror-signed and unsigned ones, raise `IncompatibleKeySpace` with the
   differing field named, whatever else is passed. It is decided from the key
   definitions, never from method names: `protein` and `aa_selection` share one
   key space, and every atomic method shares another.
2. **Processing settings, unless acknowledged.** A difference in a setting that
   changes how triangles are drawn — size filter, contact cutoff, filter or
   rule, completeness rule, candidate strategy, the alternate-location rule,
   the modified-residue rule (`modified_residue_policy`)
   (`biotsr.analysis.PROCESSING_FIELDS` and `processing_policy(sample)`) —
   raises `biotsr.analysis.IncompatibleProcessing` (an `IncompatibleSampling`)
   naming it, unless you pass `allow_mixed_processing=True`; the matrix then
   carries a warning naming each differing setting (`fm.warnings`), and every
   sample's settings are recorded. A result written before the
   modified-residue rule was recorded has `None` there, reported as
   `not recorded` and never taken to be the current rule.
3. **Partial samples, unless acknowledged.** A sample whose selection did not
   get everything it asked for (`sample.is_partial`; `sample.omissions` says
   what it lacks, `biotsr.analysis.describe_omissions(...)` in one line)
   raises `biotsr.analysis.PartialSamples`, naming each partial sample, unless
   you pass `allow_partial=True`. Its `.partial` maps each sample id to its
   omissions. With the acknowledgement nothing is dropped: `fm.comparison`
   records `partial_allowed` and `partial_samples`, each entry of
   `fm.samples` records `partial` and `omissions`, and a warning names them.
   This is independent of `allow_mixed_processing`: neither implies the
   other, and neither relaxes step 1.

```python
from biotsr.analysis import PartialSamples

chains = [AnalysisSample.from_result(
              biotsr.analyze("SYNPHOS.pdb", method="aa_ind_mix", chain=c,
                             residue_types="TYR;LYS;SER", keep_triplets=False),
              "chain" + c)
          for c in ("A", "B")]              # chain B has no SER
try:
    feature_matrix(chains)
except PartialSamples as refused:
    refused.partial                     # {'chainB': {'primary': {'not_found': ['SER']}, 'messages': [...]}}
fm = feature_matrix(chains, allow_partial=True)
fm.comparison["partial_samples"]        # the same omissions, recorded in the report
```

What was selected never blocks a comparison: structures, chains, amino-acid
types (LYS-only against TYR-only), ligands, carbohydrates, nucleotides, residue
numbers, atom subsets, one copy or all copies, and the method name are recorded
per sample and compare freely: a complete LYS-only sample and a complete
TYR-only sample need no acknowledgement. The rule is stated in full in the
[analysis methods record](methods/ANALYSIS.md), §1.1.

Two column orderings exist (`numeric_ascending`, the standard, and
`first_encounter`). They are the same object under a permutation and change no
number any measure produces.

---

## 7. Run an analysis

All of this works identically on computed results and on loaded ones.

**Common keys** — four separately named operations, because source gives the
phrase four meanings:

```python
from biotsr.analysis import (
    common_key_intersection, common_key_intersection_by_group,
    key_occurrence_percent, key_occurrence_percent_one_vs_all,
)

intersection = common_key_intersection(samples)
intersection.keys              # the keys present in EVERY sample, ascending
intersection.n_common_keys
intersection.rows              # per sample: total keys, common keys, their frequency
intersection.samples           # per sample: key-space, processing and sampling fingerprints,
                               # processing settings, selection (chain, scope, instances)
intersection.comparison        # key space, processing policies, methods, selection scopes
intersection.warnings
```

The four common-key operations apply the same comparison rule as the feature
matrix (above), and take `allow_mixed_processing=True` and
`allow_partial=True` the same way. An empty
intersection is a valid result (`intersection.is_empty`), not an error. Similarity
and clustering work on a feature matrix, so they inherit its check, and its
warnings and per-sample records travel with their results.

**Similarity and distance** — five measures, `generalized_jaccard` by default:

```python
from biotsr.analysis import (
    MEASURES, distance_matrix, generalized_jaccard, similarity_matrix,
)

generalized_jaccard(samples[0].counts, samples[1].counts)   # sum(min)/sum(max)
sm = similarity_matrix(fm)                                  # or ..., "wu"
dm = distance_matrix(fm)                                    # 1 - similarity
sorted(MEASURES)   # binary_jaccard, cosine, generalized_jaccard, sarika, wu
```

Every measure is a **similarity** in [0, 1]; a distance is always the explicit
`D = 1 − J`. Counts must be finite and non-negative (anything else raises
`ValueError`); an explicit zero is an absent key; two empty samples have
similarity 1, an empty and a non-empty one 0, and a matrix over a set with an
empty sample carries a warning.

Both are stored **condensed** — the upper triangle, row-major — because the
diagonal is a definition rather than a computed value:

```python
sm.values           # the condensed tuple, length n*(n-1)/2
sm.value(0, 1)      # one pair, by row index
sm.square()         # the full n x n form, diagonal included
sm.kind             # 'similarity' or 'distance'
sm.scale, sm.units  # 1.0 and 'fraction' -- unless you asked for .scaled(100)
sm.value_label      # 'generalized Jaccard similarity (0-1)'
sm.to_csv()         # square, sample ids on both axes, kind in the first header cell
sm.to_long_csv()    # one row per pair, with measure, kind and scale
sm.provenance()     # formula, conversion, scale, fingerprints, per-sample records
```

`dm.scaled(100)` is the percent form: its scale is 100, its diagonal and
conversions follow (`D = 100 − J`), and its exports and labels say percent. A
value outside a matrix's scale is refused at construction.

**Clustering** — average linkage and the dendrogram leaf order. The tree is the
result, and it needs no number of clusters:

```python
from biotsr.analysis import clustermap_spec, cut_hierarchy, hierarchy

tree = hierarchy(dm)            # average linkage + leaf order; no partition
tree.ordered_sample_ids()       # the dendrogram's leaves, top to bottom
tree.to_linkage_csv()           # step, left, right, height, size
tree.to_order_csv()             # sample_id, leaf_position

spec = clustermap_spec(dm, [s.label for s in samples])   # labels are optional
spec.ordered_sample_ids()  # the figure's rows, top to bottom
spec.assignment            # None: no partition was asked for, so no ARI
spec.to_row_names_csv()    # the plotted order
spec.provenance()          # linkage, leaf order, flat_partition None, ...
```

A flat partition is optional and needs one explicit cutting rule, never
inferred:

```python
cut_hierarchy(tree, n_clusters=3)  # maxclust: the lowest cut leaving at most 3
cut_hierarchy(tree, height=0.5)    # cut every merge above 0.5 (the matrix's units)

spec = clustermap_spec(dm, [s.label for s in samples], n_clusters=3)  # or height=0.5
spec.assignment.labels     # flat clusters, in sample order
spec.adjusted_rand         # ARI against your labels (a partition, every sample labelled)
spec.to_assignments_csv()  # sample_id, cluster, leaf_position[, class_label]
```

`cluster(dm, 3)` is shorthand for `cut_hierarchy(hierarchy(dm), n_clusters=3)`.
`n_clusters` and `height` cannot both be given, and without either,
`to_assignments_csv()` and `to_grouped_csv()` raise. Never choose the cut from
the known classes: setting `n_clusters` to the number of classes measures
agreement at that granularity, it is not an independent discovery of how many
classes there are. When merges tie at the cut, fewer clusters can form, and the
assignment carries a warning. The matrix is clustered as given (fractions, or
percent if you scaled it).

Class labels are optional and always supplied by you, in the matrix's sample
order — typically from your own table matched by sample id, as
[tutorial 08](tutorials/08_similarity_clustering.md) does. Nothing derives a
class from a name. If your sample ids encode the class by a convention of your
dataset, `label_from_name(name, token)` applies that convention explicitly;
the token index has no default.

---

## 8. Run a dataset, and fetch structures

A whole manifest runs from Python with the engine `biotsr batch` uses, under the
same rules (see the command-line guide, §7):

```python
from biotsr import batch

report = batch.run(
    "samples.csv", "results/",
    pdb_dir="pdb/",
    methods=["protein"],              # each row runs once per method named here
    columns={"label": "group"},       # read the label from the column "group"
    keep_columns=["group1"],          # carried into each completion record, never scientific
    triplets=False,
    jobs=1,                           # or N, or "auto": rows in worker processes
)
report.n_ok, report.n_failed, report.summary_path
for failure in report.failures:
    print(failure["rows"], failure["message"])
```

A manifest may instead name each row's method in a `method` column, but not
both. A wide dataset sheet is read with `sheet="tsr"` (or `"custom"`) and
`methods=[...]`; `maps`, `participant` and `inapplicable` are the `--map`,
`--participant` and `--inapplicable` options of the command-line guide, §7.
`cross_sides` is `--cross-side`: when a row could narrow a side of a cross
method in several ways, it chooses the column, or `"whole"` for the whole side,
as a mapping or as strings in the command-line spelling:

```python
report = batch.run(
    "wide_sheet.csv", "cross/", pdb_dir="pdb/", sheet="tsr",
    methods=["drug_protein_cross"], inapplicable="skip",
    cross_sides={"drug_protein_cross.protein": "whole",
                 "drug_protein_cross.drug": "whole"},
    # or cross_sides=["drug_protein_cross.protein=whole",
    #                 "drug_protein_cross.drug=whole"]; a column name, such as
    # "drug_protein_cross.protein": "aa_ind", narrows the side by that column
)
report.n_ok, report.n_partial          # (1, 1): 6G9M chain A holds CME 161
```

(`wide_sheet.csv` is the synthetic sheet `tests/data/manifests/wide_sheet.csv`
of [tutorial 05](tutorials/05_batch_manifest.md).)

The sides are `drug`, `carbohydrate` or `nucleotide`, and `protein`, for the
molecule-protein methods, and `1` and `2` for `protein_protein_cross` and
`protein_protein_calpha_cross`; `biotsr.cli.sheet.CROSS_SIDE_COLUMNS` lists
each side's candidate columns. A side that a row could narrow and that is not
chosen raises `batch.BatchError` naming the columns and lines. Each record's
`selection` for a cross request is both sides, `"side text | side text"`.

`batch.read_manifest(...)` resolves the rows without running them, and
`batch.run(..., dry_run=True)` plans without computing, downloading or writing.
A problem with the request as a whole — an ambiguous column, a row without a
chain, a method both in the sheet and in `methods`, two requests writing the
same prefix, an invalid download setting — raises `batch.BatchError`. Per-row
problems, such as an invalid `sample_id` (§7 of the command-line guide,
*Output names*), are reported as failed rows. A row
with fewer than three points has the status `no_triangles` and is not a
failure. A calculation that lacks something it asked for — not found, not
supported, not represented, or excluded by the completeness rule — has the
status `partial` (`batch.STATUS_PARTIAL`): it is written like `ok`, counted by
`report.n_partial`, its `message` names what is missing, and it is not a
failure (the command line exits `4` for it). With `if_exists="skip"`, an
existing output is kept only when its completion record states the current
processing rules (`request.processing`); an output kept that is itself
partial stays `skipped`, its `message` ends `kept partial result: …`, and
`report.n_skipped_partial` counts it. `report.n_incomplete` counts every
request whose result does not hold the whole requested selection: `partial`
rows, kept partial outputs, and `no_triangles` rows whose selection also
lacked something it named (`report.n_no_triangles_partial`); the command line
exits `4` when it is non-zero and nothing failed. Two rows are computed once only when they are the same calculation of
the same sample; rows with distinct `sample_id` values are separate samples,
and rows that would merge with different labels raise `batch.BatchError`.

Structures are downloaded only when asked:

```python
import biotsr

report = biotsr.fetch_structures(["1UBQ", "4HHB"], "pdb/", workers=4)
for record in report.records:
    print(record.pdb_id, record.status, record.path)
```

Each file is verified before it is installed and a valid cached file is reused;
`refresh=True` downloads again and replaces the cached file only after the new
one passes every check.  `report.ok` is true only when every id was
downloaded, refreshed or cached; a failed refresh (`refresh_failed`, old file
kept) is a failure.  `workers` (1–16), `timeout` (seconds, finite and greater
than 0, per connection and per read — not per download) and `retries` (a whole
number, 0 or more) are validated the same way as on the command line: an
invalid value raises, it is never replaced.  `batch.run(...,
download_missing=True)` fetches what a manifest needs before computing.

---

## Compute once, analyse many times

The point of §3–§7 together:

```python
import glob
import biotsr
from biotsr.analysis import AnalysisSample, common_key_intersection
from biotsr.io import save

# once, wherever the structures are
for path in glob.glob("pdb/*.pdb"):
    result = biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
    save(result, "results/", triplets=False)      # results/<PREFIX>.keys, .meta.json

# thereafter, as often as you like -- no structure is re-read
samples = [
    AnalysisSample.from_result(biotsr.load_result(p))
    for p in sorted(glob.glob("results/*.keys"))
]
common_key_intersection(samples)
```

Every sample keeps the key space, sampling profile and selection record it was
computed under, so mixing incompatible ones is refused exactly as it would be for
a freshly computed set. Loading launders nothing.

The command line does the same thing:

```bash
biotsr run -m protein --pdb 1UBQ.pdb -c A -o runs/ --no-triplets
biotsr analyze common-keys --from-results runs/ -o analysis/
```

and produces byte-identical reports to the recomputing path, which a gate checks.
There, each sample is reported under its **prefix** (`1UBQ_A_protein`; a member
prefix for each residue of a per-residue sweep). `--from-results` reads the
completion records (`PREFIX.run.json`) under the folder, and refuses an
incomplete record, a listed file that is missing or changed in size, and a
`.keys` file that no complete record lists — so a file written by `save()`
alone is not picked up there; load it with `biotsr.load_result`.

---

## Errors

Six distinguishable categories, one shallow tree:

```python
biotsr.TSRError                        # the base: "your request was wrong"
├── biotsr.InvalidRequest              # bad input or configuration
├── biotsr.IncompatibleRepresentation  # two valid things that cannot be combined
│   ├── biotsr.IncompatibleKeySpace            # the keys mean different things
│   ├── biotsr.IncompatibleSampling            # the counts were drawn differently
│   │   └── biotsr.analysis.IncompatibleProcessing  # processing settings differ
│   ├── biotsr.IncompatibleModifier            # the modifier cannot apply here
│   └── biotsr.IncompatibleSSEStratification   # the categories are not comparable
├── biotsr.analysis.PartialSamples     # samples whose selection fell short, not acknowledged
├── biotsr.MalformedResult             # a saved result is corrupt or mismatched
└── biotsr.UnsupportedFormat           # a schema version this build cannot read
```

Every one of them is also a `ValueError`, so `except ValueError` keeps working.

Anything **outside** this tree is a defect in the package, not in your call. That
is the distinction, and it is why the CLI can map the tree to exit code 2 and
everything else to exit code 1.

---

## What is public

| Import path | Contents |
|---|---|
| `biotsr` | `analyze`, `analyze_cross`, `analyze_entities`, `load_structure`, `load_result`, `fetch_structures`, `list_methods`, `list_cross_methods`, `get_method`, `list_entity_targets`, the result and identity model, the standard profiles, the error hierarchy |
| `biotsr.batch` | `run`, `read_manifest`, `BatchReport`, `BatchError`, the status names (`STATUS_OK`, `STATUS_PARTIAL`, `STATUS_NO_TRIANGLES`, `STATUS_SKIPPED`, `STATUS_INAPPLICABLE`, `STATUS_FAILED`, `STATUS_INTERRUPTED`, `STATUS_NOT_STARTED`) — a dataset manifest from Python; import it explicitly (`from biotsr import batch`): `import biotsr` alone does not load it |
| `biotsr.analysis` | the comparison rule (`check_comparable`), feature matrices, common keys, similarity, the clustering hierarchy and optional flat cuts |
| `biotsr.sse` | states, categories, annotation provenance, the stratification guard |
| `biotsr.io` | `save`, `load_result`, `native` (the one output format) |
| `biotsr.profiles` | `Profile`, `STANDARD_PROTEIN_CALPHA_V1` |
| `biotsr.cli` | `main`, `CLIError`, the exit codes |

Everything else — `binning`, `candidates`, `cross`, `encoding`, `geometry`,
`labels`, `methods`, `modifiers`, `selection`, `structure` — is internal. Their
user-facing values are re-exported above; import from there and the package can
be reorganised under you without breaking your script.

`biotsr.analysis.render` is the one module that imports a third-party plotting
library, matplotlib (`pip install "biotsr[plot]"`). It is never imported for
you, and importing `biotsr.analysis` works without it.
`render_clustermap(spec, "clustermap.png", dpi=150)` draws the heat map, both
dendrograms, the **observed class** stripe (when labels were given) and the
**cluster** stripe (only when the spec was cut by an explicit rule), each with
its own legend, and a labelled colour bar on a standalone matplotlib `Figure`,
which it returns; it never touches
pyplot, the backend or your other figures. Without matplotlib it raises
`MissingRenderDependency` naming the install command.

---

## What is deliberately not here

Candidate-strategy extensions (for example CA-kNN): every method enumerates all
triangles of its selected points. SSE for any family but the two protein Cα
methods — an atom inside a residue has no secondary structure of its own, and a
ligand, a glycan and a nucleotide have none at all. Mirror and SSE for cross. An
SSE-stratified feature representation, because no authoritative source
establishes one. Sub-selection of atoms within a single residue. Reproducing
the outputs of the historical TSR programs: every method runs its standard
profile only.

Asking for one of these raises an error naming what *is* supported, rather than
silently doing nothing.
