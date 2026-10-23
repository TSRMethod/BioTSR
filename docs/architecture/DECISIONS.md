# Architecture decision records

These records state the architecture and scientific decisions behind BioTSR 1.0
and the reasoning for each. The numbering is stable: source comments and tests
cite records as `ADR-0xx`, so every number continues to resolve, including
records that were superseded and records that described functionality the 1.0
package does not provide.

Status values:

* **Accepted** — in force as written.
* **Accepted, amended for 1.0** — in force; restated for the standard-only 1.0
  package, with what no longer applies removed.
* **Superseded by ADR-xxx** — replaced by the record named; kept for
  traceability.
* **Historical: compatibility functionality removed in 1.0** — the record
  described a convention of the historical TSR programs, or a mechanism for
  reproducing one. BioTSR 1.0 implements standard conventions only (ADR-087).
  The entry states what the historical convention was and, where relevant, why
  the standard differs.

Several records name the historical TSR programs (`TSR-Package`,
`TSR_AMINOACID_PACKAGE`, `WebsiteTSR`, `CarbohydratesTSR`,
`TSR_NUCLEOTIDE_PACKAGE` and the cross-key programs). They do so as scientific
rationale — to say what those programs computed and how the standard relates to
it. BioTSR does not reproduce their output.

---

## ADR-001 — Composable pipeline with declarative method recipes

**Status:** Accepted

**Decision.** BioTSR is one composable pipeline — selection → candidates →
labels → canonical ordering → geometry → binning → encoding → result — and each
named scientific method is a declarative recipe that fixes that method's
selector, label scheme, bin schemes, candidate strategy and supported
modifiers. A recipe contains no arithmetic.

**Rationale.** Four independent historical implementations of the
ordering + geometry + key-encoding core (`TSR-Package/TSR.py`,
`TSR_AMINOACID_PACKAGE/analyzer.py`, `CarbohydratesTSR`'s key-triplet generator
and `WebsiteTSR/key_triplets_protein.py`) are mathematically identical: the
same descending label sort, the same four ordering cases, the same `>=`
tie-breaks, the same geometry and the same mixed-radix key formula. What varies
between methods is at the edges — which points are selected, how they are
grouped and labelled, which distance-bin table is used and which triples are
candidates. Method-specific classes with shared helpers are what produced four
copies of the same ordering logic; a generic pipeline alone hides which method
is running. The recipe layer keeps the core single and the method legible.

---

## ADR-002 — `math.pi` is the angle convention

**Status:** Accepted, amended for 1.0

**Decision.** Angles are converted with `math.pi` (`PI_CONVENTION`). The
angle convention is a value passed once into the geometry stage and enters the
`KeySpace` fingerprint (ADR-003). No code path in the package contains the
literal `3.14`, and no flag, environment variable or profile produces angles
with it.

**Rationale.** Eight of the ten angle-computing historical implementations
divide by the literal `3.14`. That changes every raw angle by about 0.05% and,
on nine measured structures, moves 0.00%–3.57% of theta bins and integer keys
per structure — distinct-key membership genuinely changes. On 1UBQ chain A,
1,144 of 70,300 triangle keys move between the `3.14` and `math.pi`
conventions, and the 45,353 distinct historical keys become 45,349 (a pinned
test states this difference against the historical key column). Because
`180/3.14 = 57.3248…`, the historical pre-fold angle also ranges over
`[0, 180.0913…]` rather than `[0, 180]`; under `math.pi` that overshoot cannot
occur. A boolean `legacy=` flag threaded through the arithmetic was rejected:
it scatters a scientific choice across the math layer.

---

## ADR-003 — Key-space identity is a computed fingerprint over scientific configuration only

**Status:** Accepted

**Decision.** Two value objects, because there are two questions:

* `KeySpace` — *can these integers be compared at all?* It carries the angle
  convention, the theta, distance and label schemes (each by name, content
  fingerprint and radix), the key-formula version and the mirror convention;
  later records added the square-root convention (ADR-026), the
  canonical-ordering version (ADR-038), the median convention (ADR-039) and the
  numerical precision (ADR-053).
* `SamplingProfile` — *were these counts produced the same way?* It carries
  the selection, scope, completeness policy, candidate strategy, composition
  rule, size filter, degenerate-triangle policy and, for cross methods, the
  contact settings (ADR-043).

Each exposes a stable `fingerprint` over the canonical serialization of every
field that defines it. Renaming a scheme or a convention changes the
fingerprint, so those identifiers are part of the scientific definition. Both
identities are attached to every result and every native file. Operations that
compare key coordinates refuse a `KeySpace` mismatch with no override
(ADR-054); how sampling differences are treated is ADR-055 and ADR-090.
Nothing that cannot change a key or a count — worker count, output directory,
verbosity, cache location, input paths — enters either identity.

**Rationale.** Three verified hazards show that a key integer is meaningless
without its configuration: the angle convention (ADR-002); a label-scheme swap
with an unchanged radix (`TSR-Package/TSR.py` swaps its 20 residue labels for
13 group labels under `aa_grouping` while keeping `numOfLabels = 20`, so
numerically identical keys mean different things; ADR-015); and two distance
binners with identical declared size and cut points that differ at the top
class (ADR-014). The last two defeat any identity built from
`(dLen, dTheta, numOfLabels)` or a human-assigned version string. Content
hashing catches them. The composition rule, selection, scope and size filter
change which triangles exist, not what a key means, so they belong to
`SamplingProfile`.

---

## ADR-004 — Entity completeness is an explicit admission policy; each method's lineage rule is its standard

**Status:** Accepted (which methods apply it, and to terminal residues, is stated
by ADR-111)

**Decision.** Completeness is a predicate evaluated after an entity's atoms are
collected and before the entity contributes any point:
`CompletenessPolicy(entity_type, entity_name, collected_atoms) -> Admit | Reject(reason)`.
A rejected entity contributes no points; whether that ends the calculation
follows from scope. Every rejection is recorded in the result's `skipped` list
with its reason. For each method, the standard policy is the exact-atom-count
rule of that method's historical lineage, as decided by the TSR Research Group:
not "all available atoms", not "required atoms only", no terminal-residue or
`OXT` exception and no automatic repair, unless that lineage itself does so.

**Rationale.** The historical programs apply completeness in three different
ways — always on (`TSR_AMINOACID_PACKAGE`, where a rejected residue produces no
output at all, which is why 1UBQ's C-terminal Gly76 with its `OXT` yields
nothing), mode-dependent (`CarbohydratesTSR`, which bypasses it for ligands) and
absent (`WebsiteTSR`). A swappable policy object represents all three without
restructuring, and recording every rejection makes the rule visible rather than
silent.

---

## ADR-005 — Selection and labelling are separate stages; unlabelled points are a policy, not a silent drop

**Status:** Accepted

**Decision.** Selection produces structurally identified points; a
`LabelScheme` then assigns integer codes as a distinct stage. Unlabelled points
are handled by an explicit policy (reject the entity, drop the point and record
it, or fail), and every dropped point is recorded in the result.

**Rationale.** In `WebsiteTSR` an atom whose name is absent from the lexicon is
dropped during selection with only a printed warning, so labelling acts as a
second, invisible filter; another historical program writes `"N/A"` and fails
downstream. Separating the stages makes label coverage visible and testable.

---

## ADR-006 — One canonical-ordering component, no per-method variants

**Status:** Accepted

**Decision.** `canonical_order(labels, coords)` exists exactly once, is pure
and stateless, and takes no method flag. It is covered by branch-complete tests,
including the two exact geometric-tie branches in the "two labels equal" cases,
where `>=` makes the lower original index win.

**Rationale.** The historical implementations agree on every ordering decision
— sort direction, the four case predicates, the three geometric tie-breaks and
the comparison operators. One subtlety must be preserved: the label→index map is
keyed by label, so duplicate labels overwrite, which is safe only because each
branch reads a label that is unique within the triple. A rewrite that "fixes"
the map would read wrong indices.

---

## ADR-007 — Scope (grouping) is orthogonal to membership (selection)

**Status:** Accepted

**Decision.** Scope — per-entity or pooled — is a selection-level choice
independent of which entities qualify.

**Rationale.** `WebsiteTSR`'s `residue_ca` pools the Cα of every listed residue
into one point cloud, and `CarbohydratesTSR`'s `id_chain` pools every qualifying
carbohydrate residue's atoms (verified: 38 atoms give exactly `C(38,3)` = 8,436
triangles). They are the same operation with different membership rules. The
same separation is what makes `nucleotide_chain` the nucleotide selector under
pooled scope rather than a new method.

---

## ADR-008 — Core package carries no heavy dependencies; analysis, ML and HPC stay outside it

**Status:** Accepted (amended by ADR-086)

**Decision.** The key-generation core uses the Python standard library only.
Analysis (frequencies, common keys, similarity, clustering) is a layer that
consumes results; the core never imports it. Machine-learning pipelines and HPC
implementations are outside the package. The analysis layer is itself
standard-library-only, with one rendering carve-out (ADR-059). The package's
single runtime dependency, joblib, serves parallel batch execution only
(ADR-086).

**Rationale.** The historical programs variously depend on pandas, scipy,
seaborn, matplotlib, TensorFlow, gradient-boosting libraries, h5py, Biopython
and an MPI/HDF5 toolchain. Key generation needs none of them, and a dependency-free
core installs and runs on compute nodes with no scientific stack.

---

## ADR-009 — Structure parsing sits behind a domain model

**Status:** Accepted, amended for 1.0 (atom-name and element columns amended by
ADR-108)

**Decision.** The geometry engine never sees a PDB line: a structure backend
produces the domain model. BioTSR 1.0 ships one reader, `pdb_legacy`, which
reads fixed-column PDB records. Since 1.0.0rc6 the reader gives every atom the
PDB specification's atom-name field (columns 13–16) and element field
(columns 77–78), for every method; only the Cα test reads columns 14–16
(ADR-108). The reading that produced a result is recorded and guarded
(ADR-079). The core takes no Biopython dependency. mmCIF input is not
supported.

**Rationale.** The historical family contains several parsers, including a
fixed-column reader with off-by-one column slices and a `WebsiteTSR` reader that
corrected those slices to the specification. Keeping parsing behind a domain
model confines such choices to one layer and makes the reading a recorded fact
rather than an accident of which script ran.

---

## ADR-010 — SSE is an annotation of the frequency vector, not a distinct key space

**Status:** Accepted

**Decision.** Secondary structure is a per-triangle classifier whose output
partitions a key's frequency. It does not enter the `KeySpace` fingerprint
(implemented as ADR-064).

**Rationale.** The historical `SSE_TSR.py` uses the identical key formula and
constants (`dTheta=29, dLen=35, numOfLabels=20`) as `TSR.py`, and on 1UBQ both
produce the same 45,353 distinct keys, with the SSE program adding a per-key
category breakdown. SSE results therefore remain directly comparable to
unannotated results of the same key space.

---

## ADR-011 — Carbohydrate–protein cross uses the shared cross composition rule

**Status:** Accepted

**Decision.** Carbohydrate–protein cross TSR is a pairing under the shared
cross engine (ADR-012) on the same terms as the other cross families: the
mixed-side composition rule (`2 carbohydrate + 1 protein` and
`1 carbohydrate + 2 protein`, i.e. `a·C(b,2) + b·C(a,2)` triangles), the
atomic key space, and full physical atom identity (ADR-041).

**Rationale.** The historical carbohydrate–protein cross-key program is
algorithmically identical to the nucleotide–protein one (same composition rule,
same `dLen=18, dTheta=29, numOfLabels=112` key space, same four-case ordering),
differing only in cosmetic lines. Its composition rule was verified numerically
against its own output on 67 carbohydrate targets over 7 structures. The rule
was therefore taken from source evidence, not invented.

---

## ADR-012 — Cross-entity TSR is one engine with a pluggable pool builder

**Status:** Accepted

**Decision.** One mixed-side candidate strategy serves every cross method.
Per-pair variation is expressed as a pool specification (side selectors,
contact rule, atom identity). The merged-pool form is used: same mathematics as
the historical dual-anchor loops, with a one-line invariant instead of a
duplicated body.

**Rationale.** The historical `Cross_TSR.py` and `Drug_TSR.py` (one atom from A
with an unordered pair from B, then the mirror pass) and the cross-key scripts
(all triples of a merged pool, skipping monochromatic ones) produce the same
set of triangles, each exactly once:
`C(a+b,3) − C(a,3) − C(b,3) = a·C(b,2) + b·C(a,2)`. On 1BRS A×D this is
1,230,714 triangles, matching the historical row count with no duplicates; on
3PTB, 441. What genuinely differs is pool construction — the cutoff value, which
sides it narrows, and how an atom is identified — which is exactly what the pool
specification makes explicit (ADR-041, ADR-051).

---

## ADR-013 — Contact cutoffs build pools; they impose no per-triangle constraint

**Status:** Accepted

**Decision.** The contact cutoff is a selection-time pool filter, never a
triangle predicate. It is recorded in `SamplingProfile` (ADR-043).

**Rationale.** In every historical cross program the cutoff decides pool
membership before enumeration, and a pooled atom may form a triangle with any
other pool member. The historical 1BRS output, generated at a 5 Å cutoff,
contains triangles with a longest edge of 28.08 Å. Adding a per-triangle
distance constraint would look like a fix and would change every cross result;
if one is wanted it is a separate, named operation (size filtering, ADR-044).

---

## ADR-014 — The standard Cα distance scheme has 35 classes

**Status:** Accepted, amended for 1.0

**Decision.** `DISTANCE_SCHEME_CALPHA_35` — `83.26 <= d < 132.45 → 34`,
`d >= 132.45 → 35` — is the distance scheme of every Cα method (regular,
grouped, mirrored, size-filtered, residue-selected, SSE-annotated and Cα
protein–protein cross).

**Rationale.** `TSR-Package/TSR.py`'s binner shares these cut points but merges
the top two classes (it can only emit 1–34 although it declares `dLen = 35`),
while `WebsiteTSR`'s Cα binner distinguishes them. The TSR Research Group
designated the 35-class behaviour as the standard. Keys therefore differ from
`TSR.py`'s for any triangle with a longest edge of at least 132.45 Å; on 1UBQ,
whose Cα span never reaches that, no key changes for this reason. Atomic
distance schemes are a separate, family-level rule (ADR-021).

---

## ADR-015 — The grouped amino-acid radix is 13

**Status:** Accepted, amended for 1.0; amended by ADR-107: the grouped scheme
is now `AMINO_ACID_GROUP_23`, 23 codes over 4..33 at radix 30, and the radix
follows the span of the codes, not their count

**Decision.** `AMINO_ACID_GROUP_13` declares radix (`n_labels`) 13, matching its
13 codes.

**Rationale.** `TSR-Package/TSR.py` encodes its 13 group codes in a 20-slot
radix under `aa_grouping`, so a grouped key and an ungrouped key share one
numeric space with different label meanings. The TSR Research Group designated
as standard the behaviour in which the radix matches the label count. The group
table itself is `TSR.py`'s; only the radix differs. A standard grouped
`KeySpace` differs from an ungrouped one both in label-scheme content and in
radix.

---

## ADR-016 — The cross atom-identity defect is demonstrated, not latent

**Status:** Accepted

**Decision.** Cross pools identify atoms by full physical identity (ADR-012,
ADR-041). This record holds the evidence that the requirement corrects a real
defect.

**Rationale.** The historical nucleotide–protein cross-key program keys its pool
by `(side, atom_name)`. On PDB 1KX5 (139 real contact rows) 12 contact groups
contain a collision between two distinct physical atoms (for example `LYS 24 NZ`
and `LYS 28 NZ`), under-counting triangles by 66 (44%) across those groups when
regenerated under name-only versus physical identity with everything else held
fixed. The sparser 1TUP interface does not expose the collision, which is why a
dense interface (1KX5) is part of the test structures.

---

## ADR-017 — The atom-identity defect recurs in the carbohydrate–protein lineage

**Status:** Accepted

**Decision.** No change to ADR-012 or ADR-016. This record extends the evidence
to a second entity family.

**Rationale.** The historical carbohydrate–protein cross-key program carries the
identical `(side, atom_name)` pool mechanism. On a historical carbohydrate
contact table of 719 contact rows over 7 structures and 67 targets, 12 targets
(18%) show a real collision; across the table, 12,957 triangles under name-only
identity become 18,406 under physical identity. One further key difference in
one target was traced to canonical-ordering Case 1 (all three labels equal),
where swapping the base endpoints replaces the raw angle by its supplement;
under `math.pi` that folds to a bit-identical theta, and it changed a key only
because the historical `3.14` overshoot crossed the 79.2° bin edge. The defect is
inherent to name-only pooling, not to one dataset.

---

## ADR-018 — Mirror and amino-acid grouping are user-facing composable modifiers routed to distinct internal mechanisms

**Status:** Accepted

**Decision.** The public vocabulary has three layers: selection (which points),
candidate strategy (which triangles) and modifiers (mirror, amino-acid
grouping, size filtering, SSE — keyword arguments on one shared recipe). There
are no per-combination classes. Internally mirror is a `MirrorConvention` in
encoding, grouping is a label-scheme swap, size filtering is a candidate
predicate and SSE is a result annotator. Whether a modifier applies is checked
against what the points carry (ADR-023), not against a per-selector list.

**Rationale.** The TSR Research Group asked for modifier composability to be the
stated model. `register_method`'s `supported_modifiers` expresses per-method
compatibility without subclassing. A modifier proposed later is placed in
`KeySpace` or `SamplingProfile` by ADR-003's test: does it change what a key
means, or only which triangles were drawn?

---

## ADR-019 — Degenerate triangles are represented explicitly (`GeometryResult`)

**Status:** Accepted (the consumption policy is ADR-033)

**Decision.** Triangle geometry returns a `GeometryResult` (`valid`,
`theta_raw`, `theta`, `max_distance`, `degeneracy_reason`) rather than always an
angle. Three guards populate it: clamping of the `acos` ratio within a small
tolerance (numerical safety only), a zero-distance guard checked before any
division, and a pairwise-coincidence check. The encoder only receives valid
results.

**Rationale.** Eight of nine historical geometry implementations have no
protection against collinear or coincident triples: a `ZeroDivisionError` or an
`acos` domain error would end the run. A synthetic sweep of 200,000
near-collinear triples at PDB coordinate magnitudes found the `acos` domain
failure in 27.7% of near-exactly-collinear trials. A coincident point outside the
selected base pair completes silently with `theta = 0` in every historical
implementation. No occurrence was found in about 86,000 real triplets across
six datasets, but the risk is demonstrated.

---

## ADR-020 — Cross atom-identity safety is classified per implementation

**Status:** Accepted

**Decision.** Statements about whether a cross implementation preserves atom
identity are made per file and per pool-construction step, with evidence, never
per repository or per "cross code" as a category.

**Rationale.** Safety does not follow repository or naming: `Cross_TSR.py` and
`Drug_TSR.py` index atoms by position and are safe; the contact generators
write full composite identities and are safe in themselves; the two cross-key
programs collapse those identities back to `(side, atom_name)` and are not. A
full-identity intermediate file does not prove that the program reading it is
safe. BioTSR's own pools use full identity by construction (ADR-041).

---

## ADR-021 — Distance-bin count is a family property: 35 classes for Cα, 18 for every atomic family

**Status:** Accepted

**Decision.** Every family uses 29 theta classes. The distance scheme follows
the family:

| Family | Theta | Distance |
|---|---:|---:|
| Protein Cα (all modifiers, residue selection, Cα protein–protein cross) | 29 | **35** |
| Amino-acid atomic and the amino-acid selection modes | 29 | **18** |
| Drug / ligand | 29 | **18** |
| Carbohydrate | 29 | **18** |
| Nucleotide | 29 | **18** |
| Atomic cross methods | 29 | **18** |

The 35-class scheme is never applied to an atomic family and the 18-class
scheme never to a Cα calculation; because each scheme's content fingerprint
enters `KeySpace`, a substitution produces a different key space and the
comparison guard fires.

**Rationale.** Decided by the TSR Research Group. Cα–Cα distances span a far
wider range than intra-residue or intra-ligand atom distances, which is the
rationale `WebsiteTSR`'s own Cα binner states. The historical
`Cross_TSR.py`/`Drug_TSR.py` 17-class table is a coarser historical variant, not
the standard (ADR-040).

---

## ADR-022 — Modifiers attach at named pipeline hooks, resolved into the profile

**Status:** Accepted

**Decision.** Three hooks, in pipeline order:

```text
Selection -> Candidate strategy
  -> LABELS            aa_grouping   (swap the label scheme)
  -> Canonical ordering
  -> Geometry
  -> CANDIDATE_FILTER  size_filter   (predicate on max_distance)
  -> Binning -> Encoding
  -> KEY_SIGN          mirror        (sign from theta_raw)
  -> Features/annotations -> Result
```

Modifiers are resolved into a derived profile before the engine runs, so the
engine has no modifier branches and every combination is one code path. SSE
attaches at a fourth, annotation hook downstream of the key (ADR-064).

**Rationale.** Each position is forced by a data dependency: grouping precedes
ordering because ordering compares labels; size filtering needs `max_distance`;
mirror needs `theta_raw`. Grouping and mirror are therefore not independent —
grouping changes labels, labels decide the apex, and the apex decides
`theta_raw`. Measured on 1UBQ chain A: 33,747 mirrored triangles ungrouped
against 30,770 grouped.

---

## ADR-023 — Amino-acid grouping is restricted to the protein amino-acid label domain

**Status:** Accepted

**Decision.**

1. Applicability is checked against the label scheme's declared `domain`.
   `AMINO_ACID_20` and `AMINO_ACID_GROUP_13` are `protein_amino_acid`; the
   atom lexicon is `atom_lexicon`. `aa_grouping` requires `protein_amino_acid`
   and raises `IncompatibleModifier` otherwise, naming both domains.
2. The dataset-sheet column `aa_group` is an annotation (a family or target
   name), never bound to the `aa_grouping` modifier. In a manifest it is an
   alias of `label`.

**Rationale.** The TSR Research Group clarified that grouping applies to protein
Cα TSR only, even though atomic points carry residue metadata. Checking the
domain rather than the selector keeps the rule correct as selectors are added
and rejects every atomic family with no list to maintain. The similarly named
sheet column invites exactly the confusion rule 2 prevents.

---

## ADR-024 — Size filtering is a strict `<` predicate on `max_distance`, and a filtered candidate is not a degenerate one

**Status:** Accepted, amended for 1.0

**Decision.** A triangle is kept iff `max_distance < cutoff`, tested after
geometry and before binning. Rejected candidates are counted:
`ResultCounts.n_filtered_by_size` is separate from `n_degenerate`, and the
ledger balances (`n_candidates == n_degenerate + n_filtered_by_size + n_triangles`).

**Rationale.** The historical `size_filter` (in `TSR-Package`'s Cα and SSE
programs only) was characterized by executing it with the cutoff placed exactly
on an observed `max_distance` of 1UBQ:

| Question | Answer |
|---|---|
| Criterion | `max_distance`, the longest of the three edges, in ångströms |
| Default | `10000`, a sentinel meaning "off" |
| Comparison | **strict `<`**: at a cutoff of exactly 15.303005783178676, 14,174 triangles kept and all 26 with that `max_distance` excluded; one ULP higher admits exactly those 26 |
| Placement | gates the whole body, so rejected candidates are never binned, encoded or signed |
| Rejected candidates counted? | no — 56,126 of 70,300 vanish without trace |
| Interaction with mirror or grouping | none |

BioTSR keeps the predicate and the placement and departs in one respect:
rejected candidates are reported. A degenerate triple has no definable angle; a
size-filtered triple is well defined and merely outside the requested range, so
reporting one as the other would misdescribe the structure. Default-off and
family applicability are ADR-035.

---

## ADR-025 — Triplet records carry both the label name and the residue present

**Status:** Accepted, amended for 1.0

**Decision.** `TripletRecord` carries the integer `labels` the key was built
from, the scheme's `label_names` for those codes, and the vertices' physical
`identities`. The native triplet file writes `label_name` and `residue` as
separate columns. Under amino-acid grouping, `label_name` is the group's
representative residue name and `residue` is the residue actually present;
without grouping they coincide. Keys are unaffected either way.

**Rationale.** The historical triplets writer printed each vertex's residue name
by reverse-looking-up its label code in the active table. Under grouping that
lookup is many-to-one and returns the first residue listed for the group, so
historical grouped triplet files report a grouped ILE as `LEU` and a grouped GLN
as `ASN`. Writing both columns means no one has to reconstruct physical identity
from a collapsed label.

---

## ADR-026 — The square-root spelling is a key-affecting convention

**Status:** Accepted

**Decision.** The square-root spelling is a named `SqrtConvention` on the
profile and enters `KeySpace`. `SQRT_POW_HALF` (`x ** 0.5`) is the Cα family's;
`SQRT_MATH_SQRT` is every atomic family's. Canonical ordering uses the same
convention as the geometry stage.

**Rationale.** The Cα lineage (`TSR-Package`) writes `** 0.5` while both atomic
lineages write `math.sqrt`. `** 0.5` is evaluated by C `pow`, which IEEE-754 does
not require to be correctly rounded; `math.sqrt` is correctly rounded. Over
200,000 random magnitudes they disagree by one ULP on 233 (0.12%), and one ULP
in a distance can move a triangle across a bin boundary. The difference was
found by byte comparison: eight of 75 amino-acid residues differed only in the
last digit of theta. A numeric-tolerance comparison would have hidden it. A Cα
and an atomic result therefore differ in `KeySpace` for three independent
reasons: distance scheme, label scheme and square-root convention.

---

## ADR-027 — Atomic selection is a configured selector, not a second TSR engine

**Status:** Accepted, amended for 1.0

**Decision.** Every difference between the atomic families and protein Cα is
configuration of shared components:

| Difference | Expressed as |
|---|---|
| 18 vs 35 distance classes | a different `LadderBinScheme` (ADR-021) |
| atom lexicon vs amino-acid labels | a different `LabelScheme`, loaded as a package resource |
| `math.sqrt` vs `** 0.5` | a `SqrtConvention` (ADR-026) |
| per-residue vs whole-chain scope | a different selector and `Method.scope` |
| exact-atom-count admission | a `CompletenessPolicy` (ADR-004) |
| altLoc and element rules | an `AtomAdmission` policy |

Ordering, geometry, binning, encoding, the engine, the modifiers and the result
model are reused unchanged; a test asserts that no family-specific geometry,
ordering or encoder class exists. `selection/atomic.py` is family-agnostic: the
drug, carbohydrate and nucleotide families supply different admission, label
and completeness policies.

**Rationale.** The amino-acid atomic family differs from Cα in seven measured
ways, two of which would have produced wrong output if inherited by assumption.
None of them needs a separate engine.

---

## ADR-028 — A capability with no historical counterpart is a stated standard policy

**Status:** Accepted, amended for 1.0

**Decision.** Where the TSR Research Group requires a capability that no
historical program implements, BioTSR implements it with an explicitly stated
policy and says that the policy is a standard choice, not reproduced history. It
is tested for self-consistency (for size filtering: the candidate ledger and the
boundary behaviour) rather than against a historical reference, and no
historical reference is fabricated for it.

**Rationale.** Size filtering on the amino-acid atomic family is one example:
no amino-acid atomic program has a `size_filter`, and in-residue atom
selection does not exist in `TSR_AMINOACID_PACKAGE`. Running the Cα program and
relabelling its output, or silently assuming the Cα boundary rule is "the"
history of a family that has none, were both rejected. Atomic size filtering
uses strict `<` on `max_distance` because that is the only boundary rule
characterized anywhere in the family (ADR-024), stated as a choice.

---

## ADR-029 — The atom lexicon's radix, its unreachable slots and its loading order

**Status:** Accepted, amended for 1.0 (the `SG` treatment is ADR-030; amended by
ADR-107: no override is needed for the mismatch). The lexicon itself was
replaced by the element label table in 1.0.0rc7 (ADR-113); the radix 112 and
the 109 codes over 4..112 are unchanged.

**Decision.** The atom lexicon is a 112-radix label scheme. The
radix/alphabet mismatch was declared explicitly through
`LabelScheme.allow_radix_mismatch` with a written reason; since ADR-107 a
scheme is valid when its codes fit in one radix window, which the lexicon's
109 codes over 4..112 do, and the override is gone. The lexicon is loaded
in a fixed order (the generic column, then the specific column, which wins) and
ships as a package resource read through `importlib.resources`, never from the
working directory.

**Rationale.** Measured on the historical lexicon: `numOfLabels = 112` is an
encoder radix; the codes resolve to 109 distinct values, consecutive over 4–112,
so three radix slots are unreachable and the encoding stays injective. The
2,959 rows of the historical CSV expand to 3,034 atom-name spellings because
both columns are loaded, the second overwriting the first; reversing the passes
would change which code wins. The historical lexicon also bound the atom name
`SG` to Seaborgium (72), so every cysteine sulfur was labelled 72; that is a
scientific question about what a cysteine sulfur's label is, and it is settled
by ADR-030. Loading from the package also closes the historical defect that
required changing into the package directory before the analyzer would run.

---

## ADR-030 — The standard lexicon labels `SG` as sulfur

**Status:** Superseded by ADR-113 (1.0.0rc7). Atoms are labelled by element, so
cysteine's `SG` (element `S`) is 111 and a Seaborgium record (element `SG`) is
72 without any row for the name. Kept as the record of how the lexicon was
corrected.

**Decision.** The packaged resource `atom_lexicon_standard_v1.csv` is the
standard atom lexicon. It carries the correction decided by the TSR Research
Group:

```text
Sg,72,Sg      <- Seaborgium keeps its own chemical case
S,111,SG      <- sulfur gamma, the atom name
```

The historical lexicon files are not part of the package.

**Rationale.** The lexicon already disambiguates thirteen collisions of the same
shape by carrying a second row whose generic column is the upper-cased spelling
and whose code is the organic-element code (`Ca,99,CA` for calcium, then
`CA,106,CA` for carbon alpha). `SG` is the one case in the file where that
second row is missing. Case sensitivity alone could not fix it, because the
historical row's third column already spells `SG` in upper case; the data had
to change, and the resolver had to know which namespace it is asking (ADR-031).
The standard resource derives from the newer historical lexicon (which adds nine
atom names to the older one and changes no code) and differs from it in exactly
one mapping, `SG: 72 → 111`. Measured impact: 60 of 24,993 amino-acid atomic
triplets on 6G9M chain A, all in its six cysteine residues, and none on 1UBQ,
which has no cysteine. Because the scheme's content fingerprint enters
`KeySpace`, corrected keys cannot be silently compared with historical ones.

---

## ADR-031 — Standard lexical resolution is namespace-aware and driven by the physical atom

**Status:** Superseded by ADR-113 (1.0.0rc7): the atom-name namespace is no
longer consulted. Kept as the record of the resolver it replaced.

**Decision.** The lexicon is loaded into two namespaces (`LabelNamespaces`): the
chemical-element column and the PDB atom-name column. The standard scheme
resolves a physical atom from both the atom name and the element field:

| record | atom name | element | namespace | label |
|---|---|---|---|---:|
| cysteine sulfur gamma | `SG` | `S` | atom name | **111** |
| Seaborgium | `SG` | `SG` | element | **72** |
| carbon alpha | `CA` | `C` | atom name | **106** |
| a calcium ion | `CA` | `CA` | element | **99** |
| a ligand chlorine | `CL1` | `CL` | atom name | **93** |

1. **Elemental identity.** When the atom name and the element field name the
   same token, the record is that element (a monatomic ion, a metal); the
   element column is consulted in chemical case, which keeps `Sg` apart from
   `SG`. If it has no chemical-case entry, resolution falls through.
2. **Biomolecular atom identity.** Otherwise, the PDB atom-name namespace,
   matched exactly as the file spells it.
3. **Unknown name** → the element the structure declares, in chemical case.
4. **No element field** → the name read as a chemical symbol.

`namespace_aware` and the element table enter the scheme's content fingerprint.
No source line special-cases `CYS`/`SG`; a test asserts that.

**Rationale.** The historical loader merged both namespaces into one flat
dictionary, under which `Sg` and `SG` are the same identifier and a calcium ion
named `CA` would be labelled 106 (carbon). Step 1 fires only on positive
structural evidence, so a biological atom is never captured by an element that
merely upper-cases to the same letters. Across 642,456 atom occurrences in the
amino-acid atomic, drug, carbohydrate and nucleotide reference structures the
elemental branch fires 37,405 times (backbone `C`, `N`, `O`, nucleotide `P`) and
never produces a different label from a name-only rule: it makes the resolver
correct for elements it has not yet met and moves no existing label. The
thirteen pre-disambiguated pairs are documented, not edited.

---

## ADR-032 — The label space is a property of the family: Cα 20, grouped Cα 13, every non-Cα atomic family 112

**Status:** Accepted, amended for 1.0; amended by ADR-107: the Cα label space is
the residue labels at radix 30, grouped or not

**Decision.**

| Family | Labels | Radix | Distance classes |
|---|---|---:|---:|
| Every Cα method — regular, residue-selected (`aa_selection`), SSE-annotated, mirrored, size-filtered, Cα protein–protein cross | residue identity | **20** | 35 |
| Cα with amino-acid grouping | grouped residue identity | **13** | 35 |
| Every non-Cα atomic method — amino-acid atomic, `aa_ind`/`aa_ind_mix`/`aa_ind_specific`, atom selections, drug, carbohydrate, nucleotide and the atomic cross methods | atom lexicon | **112** | 18 |

112 is the encoding radix, not the number of occupied codes (ADR-029).

**Rationale.** Decided by the TSR Research Group. Two readings of "a protein
Cα point" existed historically: `TSR-Package` labels it by residue identity;
`WebsiteTSR`'s `residue_ca` selector looks up the literal atom name `CA` in the
atom lexicon, so every point carries code 106 and the key's label axis carries
no information. The standard interpretation of a Cα point is residue identity
in the 20-label space. Every atomic family shares one label-scheme fingerprint,
so atomic keys are mutually comparable across families, while Cα and atomic
keys never are.

---

## ADR-033 — Defined linear geometry is valid and encoded; only undefined geometry is skipped

**Status:** Accepted

**Decision.** When the denominator is valid and the geometry yields
`θ_raw = 0°` or `180°`, the triangle is encoded like any other; folding maps both
to theta class 1, and mirror still keys off `θ_raw` (0 regular, 180 mirrored).
When θ cannot be evaluated — the apex exactly on the base midpoint, a zero-length
base edge, any coincident pair or any other zero denominator — no key is
generated, the triangle is skipped and a `DegeneracyReason` is recorded. θ is
never forced to zero: `theta` and `theta_raw` are `None`. The `acos` argument is
clamped only after the denominator is known valid and only for excursions below
`DEFAULT_ACOS_TOLERANCE` (1e-9); a larger excursion is reported as
`acos_domain_error`.

**Rationale.** Decided by the TSR Research Group, refining ADR-019. A zero-area
triangle can still have a well-defined θ; conflating "very flat but real" with
"undefined" is not defensible. Tests pin each case, including an end-to-end run
of four collinear Cα points (all four triangles encoded) and the same points with
one moved onto a midpoint (one triple skipped as `apex_at_base_midpoint`).

---

## ADR-034 — The five amino-acid selection modes are five recipes over shared selectors

**Status:** Accepted; the admission sentence is superseded by ADR-095 (1.0.0rc4)

**Decision.**

| Column | Method | Points | Label space |
|---|---|---|---|
| `aa_ind` | `aa_ind` | every eligible non-hydrogen atom of every admitted residue of **one** type | atomic, 112 |
| `aa_ind_mix` | `aa_ind_mix` | the same, for **two or more** types | atomic, 112 |
| `aa_ind_specific` | `aa_ind_specific` | the same, for named residue instances (`LYS_6;GLU18`) | atomic, 112 |
| `aa_selection` | `aa_selection` | **Cα only**, for residue numbers with ranges (`10-12;48;70-71`) | **Cα, 20** |
| `aa_atom_selection` | `amino_acid_atom_selection` | exactly the atoms named as `RESIDUE_SEQNUM_ATOM`, pooled | atomic, 112 |

The first three pool every matching residue's atoms into one candidate set.
Arity is validated: `aa_ind` takes exactly one type and `aa_ind_mix` at least
two. `aa_ind_specific` accepts both spellings found in dataset sheets (`LYS_6`
and `GLU18`) and matches the residue name as well as the number; a mismatch
with the structure is reported. The `_2` dataset columns are the same methods
applied to the second chain of a protein–protein cross (ADR-045). The
type-selection recipes use the specification-correct admission (`ATOM` records,
element field, altLoc blank or `A`); `amino_acid_atomic` keeps its lineage's
admission rule (ADR-004). The two rules select identically on every test
structure.

**Rationale.** The column definitions of the research group's dataset sheets
and its clarifications settle the five meanings. Pooled scope is forced by
`aa_ind_mix`: mixing two residue types is meaningful only if their atoms meet in
the same triangles. No historical program selects residues by type, so these
recipes have no lineage admission to inherit and take the
specification-correct one rather than a documented parsing defect.

---

## ADR-035 — Size filtering is an optional, off-by-default, family-agnostic modifier

**Status:** Accepted, amended for 1.0

**Decision.** Size filtering is optional and off by default for every
registered method. When enabled it is available to every method family,
resolved into the profile like the other modifiers (ADR-022). It does not enter
`KeySpace` — it changes which triangles were drawn, not what a key means — and
it does enter `SamplingProfile`.

**Rationale.** A historical size filter exists only in the protein Cα lineage
(`TSR-Package`'s Cα and SSE programs), where its semantics were established by
execution (ADR-024). The amino-acid atomic, drug, carbohydrate and nucleotide
lineages have no size-filter parameter or code, so for those families the same
strict `<` rule is a stated standard policy (ADR-028). The same distinction
applies to mirror, where the family without a historical implementation is drug.
Tests assert, for every registered method, that `size_filter` is off by default,
is a supported modifier and leaves the `KeySpace` fingerprint unchanged.

---

## ADR-036 — SSE is a residue-level annotation and applies to the protein Cα methods only

**Status:** Accepted

**Decision.** SSE applies to `protein` and `aa_selection`, with any combination
of grouping, mirror and size filtering. The 18 SSE triplet categories are not
applied to intra-residue atomic TSR, drug, carbohydrate, nucleotide or cross
methods; requesting `sse` there is refused by name, listing what the method
supports. SSE changes no integer key: it is a category on the key's frequency
counter (`CategoryCounts`).

**Rationale.** Secondary structure is a property of a residue. An atom inside a
residue has none of its own, and a ligand, glycan or nucleotide has none at all.
Extending the categories to those families would need its own scientific
definition, not a reused table.

---

## ADR-037 — "Mixed-side cross candidate" is the canonical term

**Status:** Accepted

**Decision.** A mixed-side cross candidate is a triangle with points from both
sides: `A+A+B` or `A+B+B`, never `A+A+A` or `B+B+B`. "Non-monochromatic cross
candidate" is an accepted gloss. One-sided methods carry no composition rule;
cross methods declare theirs in `SamplingProfile.composition_rule`.

**Rationale.** "Cross" alone is ambiguous between "there are two pools" and
"each triangle spans both pools"; the standard needs one term.

---

## ADR-038 — The historical `TSR-Package` cross ordering convention

**Status:** Historical: compatibility functionality removed in 1.0

**What the convention was.** `TSR-Package`'s `Drug_TSR.py` and `Cross_TSR.py`
ordered triangle vertices differently from the rest of the family in three
ways: lexicon labels compared as strings (so label 93 sorted above 108), an
all-distinct test written as the chained `l1 != l2 != l3` (true when the first
and third labels are equal), and strict `>` instead of `>=` on distance ties.

**In 1.0.** There is one ordering function and one convention,
`STANDARD_ORDERING` (numeric collation, three pairwise distinctness tests,
`>=` ties), used by every method including every cross method. Its version
string, `legacy_four_case_v1`, enters `KeySpace.canonical_ordering_version`, so
a change to the ordering semantics would change every fingerprint rather than
silently change keys. String collation is not a defensible way to order integer
labels, which is why it was never a standard.

---

## ADR-039 — The median `s3` is a key-affecting convention

**Status:** Accepted, amended for 1.0

**Decision.** The apex-to-base-midpoint median `s3` is computed from the base
midpoint's coordinates (`MEDIAN_FROM_MIDPOINT`). The convention enters
`KeySpace.median_convention`, like `sqrt_convention`.

**Rationale.** Algebraically there is one median; numerically there are two.
Most historical programs build the midpoint from coordinates. `Drug_TSR.py` and
`Cross_TSR.py` instead used Apollonius' formula `½·√(2b² + 2c² − a²)` over the
three already-rounded edge lengths, which differs in the last bits and, because
`s3` is a denominator, can move a triangle across a theta-bin boundary. That
formula's radicand can also go negative for a near-degenerate triple, where the
historical program failed. BioTSR implements only the midpoint spelling, so
every method shares one median convention.

---

## ADR-040 — The historical 17-class cross distance scale

**Status:** Historical: compatibility functionality removed in 1.0

**What the convention was.** `Drug_TSR.py` and `Cross_TSR.py` declare
`dLen = 17`, with a ladder of sixteen 1 Å cut points plus an overflow class. Read
from source, it is intentional and internally consistent — a coarser scale than
the family's 18-class one, not a typo.

**In 1.0.** Every atomic family, including every atomic cross method, uses the
18-class scheme (ADR-021), so standard cross keys and non-cross atomic keys are
built on the same scale. The 17-class scheme is not part of the package.

---

## ADR-041 — Cross contact pools key on full physical `AtomIdentity`

**Status:** Accepted, amended for 1.0

**Decision.** Pool membership is keyed by full physical identity — structure,
model, chain, residue name, residue number, insertion code, atom name and altLoc
(`PoolIdentity.FULL_ATOM_IDENTITY`), the only identity implemented. Repeated
references to the same physical atom collapse to one pool member, and two
different atoms can never share a slot, because the slot key is the identity.

**Rationale.** The historical nucleotide– and carbohydrate–protein cross-key
programs keyed their pool by side and atom name, so two atoms sharing a name
collapsed and the later contact row overwrote the earlier coordinates. Measured
against those programs: on 1KX5 target `50_1KX5_J_DT`, `LYS 24 NZ` and
`LYS 28 NZ` share one slot, giving 4 + 4 pool atoms, 48 triangles and 43
distinct keys where physical identity gives 4 + 5, 70 and 62; over the historical
1TUP + 1KX5 contact table, 12 affected targets and 254 triangles against 320;
over a historical carbohydrate contact table (67 targets, 7 structures), 12
affected targets and 12,957 triangles against 18,406. `Drug_TSR.py` and
`Cross_TSR.py` index by position and were never affected (ADR-020).

---

## ADR-042 — Cross pools are built in memory from `AtomIdentity`, with no intermediate file

**Status:** Accepted, amended for 1.0

**Decision.** The cross flow is structure parsing → side-specific selection →
contact-pool construction over `AtomIdentity` → mixed-side candidates → shared
ordering → geometry → binning → encoding → result, with no file in between. No
contact CSV is written or read.

**Rationale.** The historical nucleotide and carbohydrate pipelines were two
programs joined by a contact CSV, and that file could not carry full identity:
it has no altLoc column. On 2AC0 target `12_2AC0_F_DG`, both conformers of
THR 123 `OG1` reach the table as two rows with the same token, so even full
identity reconstructed from the table collapses them (1 + 2 pool atoms, one
candidate), while reading the structure does not (1 + 3, three candidates).
Counts derived from a historical contact CSV are therefore a lower bound on
physical identity.

---

## ADR-043 — Contact cutoff, contact filter, pool identity and side selectors are `SamplingProfile`, not `KeySpace`

**Status:** Accepted

**Decision.** `SamplingProfile` carries the cross method, both side selectors,
the contact cutoff, the contact-filter rule, the contact projection, the pool
identity and the pool order. `KeySpace` carries none of them.

**Rationale.** A 3.0 Å pool and a 3.5 Å pool over the same structure produce
keys that mean the same thing; the wider run simply drew more triangles. Putting
the cutoff in `KeySpace` would declare those runs incommensurable, which is
false. A test asserts that changing the cutoff moves the sampling fingerprint
and leaves the key-space fingerprint unchanged.

---

## ADR-044 — The contact cutoff is pool membership and a different concept from size filtering

**Status:** Accepted (an earlier addendum on the one-sided drug filter is superseded by ADR-051)

**Decision.**

| | contact cutoff | size filter |
|---|---|---|
| acts on | an atom, before any triangle exists | a triangle, after its geometry is known |
| decides | pool membership | whether a formed triangle is kept |
| comparison | `<=` | strict `<` on `max_distance` (ADR-024) |
| default | the method's reference value | off |

An atom enters the pool because it is within the cutoff of one atom of the other
side; the resulting triangle may have edges far longer than the cutoff and is
kept. Size filtering is available on every cross method.

**Rationale.** Both are distances in ångströms and both reduce output, which
makes them easy to conflate. Every historical cross program applies the cutoff
with `<=` to pool membership and never re-checks triangle edges (3.5 Å for drug,
5 Å chain–chain, 3 Å nucleotide, configurable carbohydrate). A test builds
exactly the case of one A–B pair inside the cutoff and another edge of the
resulting triangle outside it, and asserts the triangle survives.

---

## ADR-045 — Protein–protein cross is fully atomic; a Cα side is never mixed with an atomic side

**Status:** Accepted (the Cα–Cα part is ADR-052)

**Decision.** `protein_protein_cross` is an atomic method: radix 112, 18
distance classes, both sides labelled from the atom lexicon. The two chains name
the sides, and the `_2` arguments narrow the second side with the existing
amino-acid selectors (ADR-034). A cross that would put a Cα residue side
against an atomic side is refused, naming both label spaces.

**Rationale.** `Cross_TSR.py`, the only historical protein–protein
implementation, takes every atom of each chain and labels from the atom lexicon.
No source and no decision exists for a key formula combining a 20-code residue
space with a 112-code atom space, and relabelling a Cα point as a carbon atom
would change what the key means (ADR-032).

---

## ADR-046 — Cross contributes selection, pooling and candidates; it contributes no geometry, ordering, binning, encoding or serialization

**Status:** Accepted

**Decision.** `biotsr/cross/` contains side definitions (`sides.py`),
contact-pool construction (`pools.py`), the mixed-side candidate strategy
(`candidates.py`) and the runner that hands them to the shared core
(`runner.py`). Tests assert that no cross-specific geometry, ordering or encoder
class exists, that no cross module calls `acos`, `atan2`, `degrees` or
`radians` or contains a bin ladder or the key radix, and that the runner calls
the shared core. The mixed-side strategy is a generator: for a same-side pair
the third index is drawn only from the opposite side, so a monochromatic triple
is never proposed and `count(a, b) = a·C(b,2) + b·C(a,2)` is exact.

**Rationale.** Several cross families arriving together is when a codebase grows
several near-identical engines. ADR-027 settled the same question for the atomic
families. A new cross family needs a recipe and a pool policy, not an engine.

---

## ADR-047 — Mirror and SSE are not extended to cross

**Status:** Accepted, amended for 1.0 (the coordinate-precision part is superseded by ADR-053)

**Decision.** No cross method supports `mirror` or `sse`; requesting either is
refused by name. The 18 protein SSE categories are not extended to atomic cross
methods (ADR-036).

**Rationale.** No historical cross program has a mirror flag. That the modifier
reads only `theta_raw` is a reason it *could* be offered, not evidence that it
*should*. This record originally also kept coordinate precision out of
`KeySpace`; ADR-053 reverses that.

---

## ADR-048 — Each cross family kept its lineage's chain scope

**Status:** Superseded by ADR-049

**Summary.** The historical cross programs disagree on which chains the protein
side is drawn from: `Drug_TSR.py` and `Cross_TSR.py` scope each side to one
chain, while both contact generators feeding the cross-key programs scan the
whole structure. On 1KX5 the chain-J DNA contacts protein chain D; on 2AC0 the
chain-E nucleotides contact chains B and C. This record kept each lineage's
scope and left the question open. ADR-049 answers it.

---

## ADR-049 — Both participating chains are user-specified for every cross method

*Amended by ADR-109: the two protein chains of a protein–protein cross must
differ, and no cross may select one atom on both sides.*

**Status:** Accepted, amended for 1.0 (supersedes ADR-048)

**Decision.** Every cross method — drug–protein, carbohydrate–protein,
nucleotide–protein, protein–protein all-atom and protein–protein Cα — takes both
chains from the caller and draws points from those two chains only. `chain_2`
is required and never defaulted. No whole-structure scan is implemented.

**Rationale.** Decided by the TSR Research Group: both participating chains are
specified by the user, and unrelated chains are not scanned. A default
`chain_2` would usually be wrong for a glycan or a nucleotide, which rarely
share the protein's chain, and would be wrong silently. Measured over the
historical reference targets: naming the chain changed none of the 11
carbohydrate–protein rows (the historical scan reached exactly one protein chain
on each) and 4 of the 56 nucleotide–protein rows, from two causes — on 1KX5 the
historical scan merged two protein chains of the nucleosome into one key set,
which one named chain pair cannot do, and on 2AC0 the altLoc rule of ADR-050.

---

## ADR-050 — Hydrogen is excluded from every cross calculation, and records are read under the specification's element field

**Status:** Accepted, amended for 1.0

**Decision.** Hydrogen and deuterium are excluded from contact detection,
interface identification, contact pools and vertex generation. Every cross
admission reads the specification's element field (columns 77–78, compared
without regard to case; a blank one refuses the calculation, ADR-108) and keeps
altLoc blank or `A`. Every protein side stops at its chain's `TER`; the molecule
side of the carbohydrate– and nucleotide–protein cross
(`STANDARD_CROSS_ENTITY_ADMISSION`) reads the whole of model 1. *(Amended by
ADR-109: until 1.0.0rc5 the protein side of those two methods also read to the
end of model 1.)*

**Rationale.** Decided by the TSR Research Group. Every historical cross lineage
drops hydrogen, but the contact generators did so through an off-by-one element
slice and applied no altLoc test, and `Drug_TSR.py`/`Cross_TSR.py` kept
deuterium. On the historical reference targets the only effect was on
`12_2AC0_F_DG`, where the historical read admits both conformers of THR 123
`OG1` and the standard admits conformer A only — the altLoc rule every other
standard selector already applies.

---

## ADR-051 — Contact filtering is symmetric and acts at the atom level

**Status:** Accepted (supersedes the ADR-044 addendum)

**Decision.** For every atomic cross method:

```
A_i is retained  iff  ∃ B_j : distance(A_i, B_j) ≤ cutoff
B_j is retained  iff  ∃ A_i : distance(A_i, B_j) ≤ cutoff
```

over non-hydrogen atoms, at the atom level — never residue-level or
ligand-level. `SamplingProfile.contact_filter_rule` records it as
`symmetric_atom_contact_filter`. An explicit entity atom selection is applied
first and the selected atoms are then contact-filtered. A one-sided filter
(`contact_filter="side_a_only"` or `"side_b_only"`) or no filter
(`no_contact_filter=True`) is available only when the caller asks for it, and is
recorded under its own rule name.

**Rationale.** Decided by the TSR Research Group: an entity or residue is not
admitted whole because one of its atoms makes contact. `Drug_TSR.py` narrowed
only the protein and admitted every heavy atom of the ligand. On 6G9M chain A,
ESW 404 at 3.5 Å, that one-sided rule gives pools of 10 protein and 37 ligand
atoms and 8,325 triangles; the symmetric rule gives 10 and 10, and 900. The 27
dropped ligand atoms are exactly those with no protein atom within 3.5 Å.

---

## ADR-052 — Cα protein–protein cross: interface residues by heavy-atom contact, projected onto Cα, in the radix-20 key space

**Status:** Accepted (resolves the open part of ADR-045; amended by ADR-107: the
Cα key space is radix 30, and supported modified residues take part; amended by
ADR-109: a narrowed side narrows the interface)

**Decision.** `protein_protein_calpha_cross` uses a two-step pool rule:

1. **membership** — residue R on chain A and residue S on chain B are interface
   residues if at least one pair of their non-hydrogen atoms (side chains
   included) is within the cutoff;
2. **representation** — each interface residue contributes its Cα, whatever the
   Cα's own distance.

`SamplingProfile.contact_projection` records
`PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA`, against `ATOM_LEVEL` for every
other cross method. Points are labelled in the 20-code amino-acid space on the
35-class Cα scale, so the key space is identical to ordinary protein Cα TSR. An
interface residue without a usable Cα contributes no point and is recorded as
`interface_residue_without_calpha`; no other atom is substituted. `aa_grouping`
is not offered.

**Rationale.** Decided by the TSR Research Group. Testing Cα–Cα distance instead
would drop exactly the long-side-chain contacts that define most protein–protein
interfaces. The pinning test is built so the difference is unmissable: LYS 1
reaches GLU 1 through NZ at 2 Å while the two Cα atoms are 12 Å apart, at a 3 Å
cutoff; the standard rule pools both residues and a Cα–Cα rule pools none.

---

## ADR-053 — Numerical precision is a `KeySpace` field

**Status:** Accepted, amended for 1.0 (supersedes the coordinate-precision part of ADR-047)

**Decision.** `KeySpace.numerical_precision` records the floating-point
precision of the coordinates and arithmetic. BioTSR computes in double precision
(`"float64"`) everywhere, for every method.

**Rationale.** `Bio.PDB` stores coordinates as float32, so the two historical
`TSR-Package` cross programs worked from single-precision coordinates. The
perturbation is about 1e-6 Å and occasionally decisive: on 1BRS A/D, 2 of
1,230,714 triangles fall on opposite sides of a theta-bin boundary under the two
precisions (37.339998 against the 37.34 cut; 52.430001 against 52.43). Two runs
over the same file under the same settings can therefore disagree on a key,
which is exactly what `KeySpace` exists to rule out — the same argument that put
`sqrt_convention` and `median_convention` there. Compact storage, if ever
wanted, is a serialization concern and not this field.

---

## ADR-054 — A shared `KeySpace` is a precondition of every analysis that touches key coordinates

**Status:** Accepted

**Decision.** `require_common_key_space()` is defined once, in
`analysis/samples.py`, and is called through `check_comparable()` by
`feature_matrix`, `common_key_intersection`, `common_key_intersection_by_group`,
`key_occurrence_percent` and `key_occurrence_percent_one_vs_all`; similarity,
clustering and the clustermap inherit it through the feature matrix. A mismatch
raises `IncompatibleKeySpace` naming every differing field. There is no
override.

**Rationale.** The analysis layer's columns are key coordinates. Combining two
key spaces would put two meanings in one column and then sum across it. The
historical feature builders union whatever key files a directory holds, so a
directory with grouped and ungrouped files would silently produce a matrix whose
columns mean two things.

---

## ADR-055 — A `SamplingProfile` difference is preserved and exposed by feature matrices and similarity

**Status:** Superseded by ADR-090 (1.0.0rc3)

**Decision.** `feature_matrix` permits sampling variation and attaches a
`sampling_variation()` summary — the number of distinct profiles, whether they
are uniform, which fields vary and which samples carry which fingerprint. The
summary travels onto similarity matrices and clustering results. An analysis
whose specification demands identical sampling passes
`require_same_sampling=True` and receives `IncompatibleSampling`. The common-key
operations instead require one sampling rule unless mixed sampling is explicitly
allowed (ADR-090).

**Rationale.** Comparing the same interface at two contact cutoffs, or filtered
against unfiltered sampling, may be exactly the question. A blanket prohibition
would be wrong and blanket silence careless; whether a difference is admissible
belongs to the analysis specification.


*Superseded (1.0.0rc3).* Feature matrices, similarity and clustering now apply
the same sampling-rule check as the common-key operations (ADR-090); the
`require_same_sampling` parameter is gone. A study that compares two cutoffs
on purpose passes `allow_mixed_sampling=True` and the result says so.
---

## ADR-056 — Feature columns are ordered by integer key value

**Status:** Accepted, amended for 1.0

**Decision.** The standard column order is increasing integer key value
(`numeric_ascending`). `first_encounter` (insertion order, as one historical
feature-selection workflow used) is also offered. Hash order is not offered.

**Rationale.** The historical analysis programs used three orderings:
`list(set(...))` hash order, which changes between runs under a different
`PYTHONHASHSEED`; `sorted()` over decimal strings, which puts `"10"` before
`"9"` and every negative mirror key before every positive one; and first
encounter. Every measure the analysis layer implements is a sum over all columns
with the same permutation applied to both operands, so column order changes no
number; the standard order is chosen for reproducibility and readability. A
non-reproducible order cannot be offered as a choice.

---

## ADR-057 — Every measure is defined as a similarity; distance is an explicit conversion

**Status:** Accepted

**Decision.** All five measures are similarities in [0, 1]. `SimilarityMatrix`
carries its `kind` (`"similarity"` or `"distance"`), and its diagonal follows by
definition (1.0 or 0.0). Conversion is explicit (`to_distance()`,
`to_similarity()`); `average_linkage` refuses a similarity matrix rather than
inverting it. When the denominator is zero, the similarity is 1.0.

**Rationale.** The historical HPC program writes `1 − J` unless asked
otherwise, and a historical workflow wrote variables named as distances to files
named `*_similarity*.csv`: the historical "similarity" files contain distances.
A layer that inherited that naming would inherit the ambiguity. The zero
denominator follows the HPC program; the Python programs instead printed an
error and silently shortened the row, which is a defect.

---

## ADR-058 — The feature matrix is stored sparsely and is dense-identical

**Status:** Accepted

**Decision.** Rows are `key -> count` mappings; `row(i)` and `dense()` expand
them with explicit zeros in column order. A key absent from a row and a key with
count zero are the same fact; an explicit zero is dropped at construction.

**Rationale.** The mathematical object is unchanged — an *n × |K|* matrix over
the union of keys. The measures are computed over the sparse intersection
because `min(a_k, 0) = 0` and `Σ max = Σa + Σb − Σ min`, an identity rather than
an approximation. Canonicalising zeros makes every presence-based measure agree
that presence means `value > 0`. Measured on one large acceptance case (5.6%
non-zero cells): 2.6 MiB stored against 5.9 MiB dense. On a small, dense case
(31.9% non-zero, 1BNA) sparse storage is larger, which is recorded rather than
hidden.

---

## ADR-059 — The numerical analysis layer is standard-library-only; `analysis/render.py` is the single carve-out

**Status:** Accepted, amended for 1.0.0rc3

**Decision.** Feature matrices, common keys, the five measures, average
linkage, `maxclust`, leaf ordering and the adjusted Rand index are implemented
on the standard library (ADR-060). `analysis/render.py` alone may import a
third-party library, and only matplotlib — the `plot` extra, which declares
exactly that; nothing else in the package imports it. Without it a figure raises
`MissingRenderDependency` naming `pip install "biotsr[plot]"`, and the numerical
result is complete without it.

**Amendment (1.0.0rc3).** The renderer used `seaborn.clustermap` and imported
pandas, while the `plot` extra declared only matplotlib and seaborn; and
`seaborn.clustermap` needs scipy at run time, which seaborn does not install, so
a clean `pip install "biotsr[plot]"` could not draw. It also switched the
caller's matplotlib backend and closed every open figure. It now draws with
matplotlib alone, on a standalone `Figure` (no pyplot), from BioTSR's own
linkage and leaf order.

**Rationale.** Key generation must run where scipy, scikit-learn and seaborn do
not. Tests hold the line: the import-purity test skips exactly that one path, a
second test bounds `render.py` to matplotlib and compares it with the `plot`
extra, and a third asserts no other module imports it.

---

## ADR-060 — The standard-library clustering follows scipy's algorithms; `maxclust` follows its definition

**Status:** Accepted, amended for 1.0.0rc3; amended by ADR-099 (1.0.0rc4)

**Decision.** Average linkage follows scipy's algorithms rather than textbook
ones: the `nn_chain` search with its "prefer the previous chain element" rule,
the stable sort of the linkage matrix by height, union-find relabelling, and the
Lance–Williams update in scipy's operand order
`(n_x·d_xi + n_y·d_yi)/(n_x + n_y)`. The leaf order follows `dendrogram`'s
default traversal. `maxclust` follows its definition: the **lowest** threshold
— below every merge, or one of the merge heights — that leaves at most *t*
clusters, so with distinct merge heights exactly `min(t, n)` clusters form, and
fewer only when merges tie at the cut (the result then carries a warning).
Flat clusters are numbered in scipy's `cluster_monocrit` traversal order. The
number of clusters is always the caller's choice; nothing infers it.

**Amendment (1.0.0rc3).** Until rc2 `maxclust` reproduced scipy 1.11.1's binary
search, which never evaluated the lowest candidate cut: it could return fewer
clusters than requested even with distinct heights (612 of 9,413 cuts in the
evidence below), and two clusters for `n == 2, t == 1` (an out-of-bounds read).
`docs/methods/ANALYSIS.md` already stated the definition, so the code and the
documentation disagreed. The cut now follows the definition, as later scipy
releases do; linkage and leaf order are unchanged.

**Evidence and its scope** (`validation/scripts/analysis_scipy_equivalence.py`,
records `validation/evidence/analysis/SCIPY_EQUIVALENCE_scipy-<version>.json`):
400 random matrices up to n = 45, a third built from a three-value alphabet so
that ties are common, compared with two scipy releases.

* **scipy 1.11.1** (scikit-learn 1.6.1, Python 3.9): 400/400 linkage matrices
  bit-identical and 400/400 leaf orders identical; on the same tree the
  `maxclust` partitions are identical except the 624 cuts where scipy 1.11.1
  misses the cut (612 with fewer clusters than the definition gives, 12 with
  two clusters for `n = 2, t = 1`); 400/400 ARI values identical.
* **scipy 1.18.1** (scikit-learn 1.9.1, Python 3.13): its linkage differs in 177
  of 400 matrices — 175 contain tied distances; in the two without any tie the
  merges are identical and one merge height differs by 1.1e-16 — and its leaf
  order in 12; **on the same
  tree** the `maxclust` partition is identical in all 9,413 cuts and the leaf
  order in all 400 (scipy numbers the clusters differently in 380 cuts); 400/400
  ARI values identical.

So agreement is claimed precisely: the linkage with scipy 1.11.1; the
`maxclust` partition and the leaf order, given the same tree, with scipy 1.18.1;
nothing for scipy releases not checked. Where distances tie, two releases of
scipy can return different trees, and BioTSR's is one valid, deterministic
choice among them.

---

## ADR-061 — The historical clustering input was a quantized distance matrix

**Status:** Historical: compatibility functionality removed in 1.0

**What the convention was.** The historical pipeline stored the distance matrix
at single precision, clamped it to [0, 1] and printed it with three decimals;
the clustering program read that text file. Its clustering input was therefore
quantized to three decimals. On the historical program's key sets for 1UBQ's 75
amino-acid residues, that changes the flat clusters, the leaf order and the
adjusted Rand index (0.4132 quantized against 0.4188 at full precision); on
1BNA's eleven nucleotides it changes nothing.

**In 1.0.** Clustering operates on whatever matrix it is given, at full double
precision; the package offers no quantization step. A study comparing against a
historical figure should state that the historical input was quantized.

---

## ADR-062 — Four common-key operations keep four names

**Status:** Accepted

**Decision.** Four separately named functions — `common_key_intersection`,
`common_key_intersection_by_group`, `key_occurrence_percent` and
`key_occurrence_percent_one_vs_all` — with no `mode=` parameter. All four are
presence-based; frequency appears only as a reported sum over the intersection.
The empty intersection is a valid result. Percentages are stored unrounded;
rounding is an export concern.

**Rationale.** Four historical programs share the phrase "common keys" and
compute four different things: the global intersection, the group-wise
intersection, the per-key occurrence percentage and the one-versus-all table.
The first two return a subset of keys and the last two a table over every key;
collapsing them would let a caller who wanted a subset receive a table of a
different shape. No historical program computes a pairwise intersection or a
frequency-threshold subset, so neither is implemented.

---

## ADR-063 — An individual β-strand is `(sheetID, serNum)`

**Status:** Accepted, amended for 1.0

**Decision.** The identity of an individual β-strand is the pair
`(sheetID, serNum)` from the `SHEET` record, carried as
`biotsr.sse.states.SSEElementId`, a two-field value object compared
componentwise (so `("A1", "1")` and `("A", "11")` stay distinct). The convention
is named `sheet_id_and_serial`
(`biotsr.sse.assignment.STRAND_IDENTITY_SHEET_ID_AND_SERIAL`) and is the only
one implemented. `sheetID` is used only to disambiguate strands; strands of the
same sheet remain distinct elements. The convention is recorded in the SSE
annotation provenance (`annotations["sse"]["strand_identity"]`) and in neither
identity: it moves no key and draws no triangle.

**Rationale.** The TSR Research Group clarified that the SSE categories
distinguish individual β-strands, and that one β-sheet contains one or more
strands. Per the PDB specification `serNum` is a strand's index within its
sheet and restarts at 1 in every sheet, so the historical `SSE_TSR.py`, which
used `serNum` alone, merged strands of different sheets that share a serial.
`sheetID` alone was rejected: it would collapse every strand of a sheet into one
group. Measured over eleven historical reference runs, `serNum` alone versus
`(sheetID, serNum)`:

| | triangles | changed | % |
|---|---|---|---|
| 1UBQ (8 modifier combinations) | 524,168 | 0 | 0.0000 |
| 1BRS A (regular, mirror) | 408,312 | 0 | 0.0000 |
| 1TUP A | 1,235,780 | 27,228 | 2.2033 |
| **total** | **2,168,260** | **27,228** | **1.2558** |

1UBQ and 1BRS chain A each carry one β-sheet, so nothing can change; 1TUP chain
A carries two whose serials 1–4 collide. Every change is toward more distinct
strands, and a triangle changes category iff it contains two sheet residues
sharing a serial under different sheet ids. Helix identity is unchanged: a
`HELIX` record's `serNum` is unique across the file, and substituting `helixID`
changes no category on any of these runs.

---

## ADR-064 — SSE is an annotation on the shared result, in neither identity

**Status:** Accepted

**Decision.**

| | changes what a key means? | changes which triangles? | identity |
|---|---|---|---|
| `aa_grouping` | yes | no | `KeySpace` |
| `mirror` | yes | no | `KeySpace` |
| `size_filter` | no | yes | `SamplingProfile` |
| **`sse`** | **no** | **no** | **neither** |

For identical labels and geometry, an annotated and an unannotated run produce
the same integer keys and share both fingerprints. Three structural choices make
this true by construction: the annotation hook (`ModifierHook.ANNOTATION`) runs
after the final signed key exists; the annotator receives the three ordered
`AtomIdentity` values and nothing else — no coordinate, label, bin or key (a test
reads the call site and asserts the argument list); and `CategoryCounts` carried
a per-key category breakdown from the start, so no schema changed. `sse` is
listed in `supported_modifiers` only so that requesting it where the categories
are undefined is refused by name. With mirror and SSE together, the category
takes the bucket and the mirror state is carried by the key's sign.

**Rationale.** The historical `SSE_TSR.py` was a 452-line copy of `TSR.py` that
then drifted. There is no second engine, no `SSETSR` class and no duplicated
ordering, geometry, binning, key arithmetic or parsing; structural tests enforce
each absence.

---

## ADR-065 — Coil is a state, not an absence; a missing annotation is reported rather than reclassified

**Status:** Accepted

**Decision.** A residue no `HELIX` or `SHEET` record covers is classified as
coil, so the classification is total and the category totals equal the retained
triangle count exactly. Whether any record covered the chain is recorded beside
the states: `SSEAssignmentMap.has_records`, an explicit warning on the result
when none did, `-` rather than `NONE` in the native triplet file for an
unannotated run, and the helix/sheet/coil coverage in the sidecar.

**Rationale.** The historical program classifies an uncovered residue as coil,
and it is the only implementation, so the rule comes from source. The
three-state collapse cannot by itself distinguish "in no record" from "this file
carries no secondary-structure records at all": a predicted model or a stripped
file yields every triangle as `17_3c`, which reads like a finding and is a
missing input. That is closed beside the states, without a fourth state.

---

## ADR-066 — The historical SSE aggregation defect

**Status:** Historical: compatibility functionality removed in 1.0

**What the defect was.** `SSE_TSR.py` allocated its 18-element category counter
once, outside the triangle loop, then both mutated it through an alias and
rebound it, so each key's vector ended up aliasing a list created for another
key. On 1UBQ chain A the historical per-key category file sums to 111,017
against a true 70,300 triangles, and 18,119 of 45,353 rows disagree with their
own key's frequency; the same inflation is present in every historical SSE run
measured. The per-triangle category column is unaffected and correct.

**In 1.0.** `CategoryCounts` gives every key its own counter, so there is no
shared mutable state for an alias to form in and the defect cannot arise. The
package does not write the historical per-key category file. Analyses derived
from historical per-key category files should take the defect into account;
analyses derived from historical per-triangle records are unaffected.

---

## ADR-067 — No SSE feature representation is defined; the stratification is preserved and exposed

**Status:** Accepted, amended for 1.0

**Decision.** No SSE feature representation (an 18-times-wider vector, an
SSE-stratified Jaccard, per-category clustering) is defined. Instead:

* the stratification is preserved on the result — `CategoryCounts` per key,
  `TripletSSE` per triangle, the assignment summary in provenance;
* `biotsr.sse.stratification` exposes it as `{category: {key: n}}`, the dense
  18-vector in canonical order, and a conservation check, and offers no measure,
  distance or clustering;
* `AnalysisSample.from_result` sums the categories, so an SSE result and its
  unannotated twin reduce to the same sample in every existing analysis;
* `biotsr.sse.require_common_strand_identity` refuses to combine
  category-stratified results annotated under different strand-identity
  conventions, raising `IncompatibleSSEStratification` and naming what is still
  comparable.

**Rationale.** Every evidenced analysis — feature matrices, common keys, the
five measures, clustering — sums a key's categories first; none reads a category
breakdown. Building one would be inventing science. The strand-identity
convention is deliberately in neither identity (ADR-063), so `KeySpace` and
`SamplingProfile` cannot catch a mismatch; the dedicated guard can, and costs
nothing to have before any category-level analysis exists.

---

## ADR-068 — One canonical executable, `biotsr`, and `python -m biotsr` beside it

**Status:** Accepted

**Decision.** One console script, `biotsr` (`biotsr = "biotsr.cli:main"`), with
subcommands generated from the method registry. `python -m biotsr` is an exactly
equivalent second entry point, for environments whose `bin` directory is not on
`PATH`.

**Rationale.** The historical family has a Python function with import-time
side effects, a 19-flag `argparse` program, and per-family `sys.argv` scripts,
with the same concept spelled three ways. Per-method scripts would undo the
registry; a single `--entity-type` switch conflates method, selector and label
space; free-form `--dtheta`/`--dlen`/`--num_labels` integers would let a user
invent an unrecorded key space. Registering a method makes it available on the
command line with no new script or flag, and `biotsr methods` cannot describe a
method the package will not run.

---

## ADR-069 — The CLI is a thin wrapper and holds no scientific default

**Status:** Accepted, amended for 1.0

**Decision.** `biotsr/cli/` orchestrates and computes nothing. It calls
`biotsr.api`, `biotsr.io`, `biotsr.profiles` and the method registry, and
contains no geometry, ordering, binning, key formula, contact selection, SSE
category logic or serializer, and no default for any scientific argument: a flag
the user does not pass is not passed, so the package's own default applies. The
CLI hands the API a structure path, not a pre-parsed structure. This is
enforced by an equivalence test that runs representative workflows through the
CLI and the Python API and requires byte-identical files and matching
fingerprints, and by structural tests over every CLI module (no float literals,
no serializer filename fragments, no imports of scientific internals, no
CLI-side scientific defaults).

**Rationale.** The standing risk of a command line over a scientific library is
a second opinion accumulating — a default here, a special case there — until two
paths disagree. The historical family shows the end state: a `10000` sentinel
for "off", per-program hard-coded cutoffs, and four copies of the ordering
logic. A scientific change lands in one place and the CLI follows.

---

## ADR-070 — The canonical dataset manifest is long, not wide

**Status:** Accepted; extended by ADR-097 (1.0.0rc4), which reads wide sheets method-aware

**Decision.** One manifest row is exactly one `biotsr run`, with the column
names the flags use. A blank cell means "not set", never a default; an unknown
column is an error naming the accepted ones; two aliases of one canonical column
populated in the same row is an error naming both. Aliases for the dataset
sheets' own spellings are accepted (`protein`→`pdb_id`, the chain columns→`chain`,
`aa_group`→`label`, `drug_name`/`carb_name`/`nt_name`→`entity_name`,
`car_atom_selected`→`atoms`, and others). What a manifest refuses outright is
ADR-093.

**Rationale.** The research group's dataset sheets are wide: one row is one PDB
entry carrying up to a dozen calculations, each with its own chain column, so a
row does not describe a run and "chain" is ambiguous within it. Driving the wide
sheet directly would mean guessing which chain belongs to which calculation;
per-method formats or `--*_column` rename flags would bring back the
fragmentation ADR-068 removes.

---

## ADR-071 — A run owns a directory; nothing is overwritten without being asked

**Status:** Superseded by ADR-098 (1.0.0rc4)

**Decision.** `--output` is required and is the only place the CLI writes. Each
sample owns a directory under it; the scientific files inside are written by
`biotsr.io` under its own names. Beside them the CLI writes `run.json`: the
command, the package version, the resolved profile and the SSE strand-identity
convention, each result's warnings, and only the arguments the user supplied.
Default sample names are built from structure, chain, second chain, method and
target. `--if-exists` is `fail` by default; `overwrite` replaces a sample
directory wholesale; `skip` makes an interrupted batch re-runnable. A sample is
built in a sibling staging directory (`.<sample>.partial-<pid>-<host>-<token>`)
and moved into place with one `os.replace`. The collision policy is applied
before computing. Sample names are checked as directory names inside the output
root (ADR-091).

**Rationale.** The historical programs write into results directories by
convention, a second run lands on top of the first, and no output records which
convention produced it. A self-describing, staged results tree is safe to
re-enter, and an interrupted run leaves something obviously incomplete rather
than something that looks finished.

---

## ADR-072 — A batch continues and attributes every failure; an analysis stops

**Status:** Accepted, amended for 1.0

**Decision.** `biotsr batch` defaults to continuing: every row is attempted and
appears in `batch_summary.csv` keyed by its manifest line number; every failure
also appears in `batch_failures.csv` with its reason; a failed sample leaves no
output directory. `--on-error stop` writes both reports before halting.
`biotsr analyze common-keys` stops on the first failure. Exit codes: 0 success;
1 internal error; 2 user or input error; 3 a batch with rows failed, interrupted
by a pool failure or never started, or a fetch with failures; 130 interrupted.
`--on-error stop` exits 3 after a computed row fails, and 2 if a row could not
even be planned.

**Rationale.** A research dataset has bad rows, and a 5,000-sample run should not
be lost to row 12. A key intersection silently taken over a reduced sample set
is a wrong answer that looks like a right one, so the analysis command does not
continue. Exit code 3 lets a scheduler tell partial success from success and
from failure; interruption and pool failures are ADR-092.

---

## ADR-073 — Only common keys are exposed on the command line in 1.0

**Status:** Accepted

**Decision.** The command line exposes one analysis, `biotsr analyze
common-keys`, driven by a manifest or by previously computed native results
(`--from-results`). The intersection requires a shared key space (ADR-054) and
one sampling rule (ADR-090, superseded by ADR-094); `--group-by label` gives the grouped form.
Similarity measures, feature matrices, clustering and the clustermap are
available from Python only.

**Rationale.** Similarity and feature matrices need a user-facing vocabulary for
measure and ordering choices, clustering needs that plus a cluster-count
vocabulary, and the clustermap needs the rendering dependencies (ADR-059).
Common keys is the most-used operation, is a pure function with its precondition
built in, and wraps thinly. What the CLI offers, it offers completely.

---

## ADR-074 — The native result format carries an explicit schema version, independent of the package version

**Status:** Accepted (the current schema is 2, ADR-098)

**Decision.** The sidecar declares `result_schema_version`, a top-level integer,
`1` when this was decided and `2` since 1.0.0rc4 (ADR-098). Schema 1 carries a `provenance` block (`profile`, `base_profile`,
`legacy_compatibility_profile` — always false in 1.0 and kept so the schema is
unchanged — `modifiers`, `target`, `truncated_at_line`, and `cross` for cross
results) and `selection.residues_used`. The version is an integer so a reader
can say "newer than anything I know"; it is not tied to `biotsr.__version__`; it
is on the sidecar only, because the columnar files declare their shape in their
header line; and an unversioned sidecar is refused, not assumed to be version 1.

**Rationale.** A reader that cannot tell one format from another must either
refuse everything or guess, and guessing can produce a `KeySpace` that parses
and means something else. Tying the schema to the package version would force a
format migration on every release or make the field useless.

---

## ADR-075 — What a loaded result carries, and what it deliberately does not

**Status:** Accepted, amended for 1.0

**Decision.** `biotsr.load_result` returns a real `TSRResult`. Every key and
per-category frequency, `KeySpace` and `SamplingProfile` with bit-identical
fingerprints, the profile provenance, the SSE annotation provenance, the
selection summary, the source provenance, the candidate ledger, the skip records
and the warnings round-trip exactly, and re-serializing a loaded result
reproduces the keys file and the sidecar byte for byte. Two things are outside
the boundary: `loaded.triplets is None` and `loaded.key_order == ()`. The result
records that it was loaded in `annotations["loaded_from"]`, with the two files
and the schema version.

**Rationale.** A second, reduced result type would fork every consumer. The
native triplet file writes six of `AtomIdentity`'s eight fields per vertex;
`alt_loc` is per-atom (one structure contributes both `""` and `"A"` atoms) and is
in no column, so reconstructing triplets would be a guess. First-seen key order
is used by no BioTSR output, and serializing it would cost about 408 KB against
a sidecar of under 4 KB on 1UBQ chain A — more than 100 times larger. Every analysis is
defined over `AnalysisSample` — a frequency vector plus the two identities — so a
loaded result is fully usable for analysis.

---

## ADR-076 — One shallow package error hierarchy, every member a `ValueError`

**Status:** Accepted

**Decision.** `biotsr/domain/errors.py`:

```
TSRError                        (ValueError)
├── InvalidRequest                              bad user configuration or input
├── IncompatibleRepresentation                  two results cannot be combined
│   ├── IncompatibleKeySpace
│   ├── IncompatibleSampling
│   ├── IncompatibleModifier
│   ├── IncompatibleSSEStratification
│   └── IncompatibleStructureInterpretation     (ADR-079)
├── MalformedResult                             corrupt or mismatched saved result
└── UnsupportedFormat                           unknown schema version
```

`TSRError` subclasses `ValueError`, so existing `except ValueError` code keeps
working. Programming errors (`AttributeError`, a `TypeError` from a wrongly typed
argument, an internal `AssertionError`) stay outside the tree, which is what lets
exit code 2 mean "your request" and exit code 1 mean "our defect" without a list
in the CLI. The base classes live in `domain`, the package's dependency sink.

**Rationale.** Before this, bad chains, missing files, malformed manifests and
out-of-range cutoffs were undifferentiated `ValueError`s, and the CLI kept its
own taxonomy of user errors, free to drift from the package.

---

## ADR-077 — The public API boundary is declared

**Status:** Accepted, amended for 1.0

**Decision.** The supported surface, stated in the README and the
[API guide](../API_GUIDE.md):

* `biotsr` — the façade (`analyze`, `analyze_cross`, `analyze_entities`,
  `load_result`, the `list_*` functions, `fetch_structures`), the result and
  identity model, the error hierarchy, `ATOM_LEXICON_STANDARD_112`;
* `biotsr.analysis` — the analysis layer, via its own `__all__`;
* `biotsr.sse` — annotation, categories and the stratification guard;
* `biotsr.io` — `save`, `load_result` and the native format;
* `biotsr.profiles` — `Profile` and `STANDARD_PROTEIN_CALPHA_V1`;
* `biotsr.cli` — `main`, `CLIError` and the exit codes.

Everything else (`binning`, `candidates`, `cross`, `encoding`, `geometry`,
`labels`, `methods`, `modifiers`, `selection`, `structure`) is internal.
`profile=` accepts only `None` or the method's own default profile (ADR-087).

**Rationale.** A boundary nobody states is whatever `dir()` happens to show. A
hundred-name top-level namespace is not a boundary either; five import paths are
what a user needs.

---

## ADR-078 — The native serialized-format compatibility policy

**Status:** Accepted, amended for 1.0

**Decision.** Anchored on `result_schema_version` (ADR-074):

* a reader accepts the versions it declares and refuses every other, newer or
  older (`UnsupportedFormat`); it never reads a file partially and never infers a
  schema from contents;
* within a version the format is append-only: no existing field changes name,
  path, type or meaning, and no new field is required; anything else is a version
  bump;
* a version bump obliges the new reader to read the old version for as long as
  it is listed in `SUPPORTED_RESULT_SCHEMA_VERSIONS`; dropping one is a
  deliberate, documented act.

`run.json` is not part of the format: it records the command, not the result,
and `biotsr.io.native` neither writes nor reads it. Python API compatibility is
ADR-077.

**Rationale.** "Can this build read this file?" should have a mechanical answer.
A permissive reader that ignores unknown fields cannot tell whether an unknown
field in an identity block is part of the identity.

---

## ADR-079 — Structure-interpretation provenance is recorded outside both identities and guarded at the analysis boundary

**Status:** Accepted, amended for 1.0

**Decision.**

1. The semantic that needs provenance is the atom-name interpretation,
   `structure_interpretation`. BioTSR 1.0 declares one, `legacy_columns_v1`; a
   second, `spec_columns_v1`, is named but not shipped, so results computed now
   are distinguishable from it should it ever exist.
2. It is derived from the recorded `SourceInfo.parser_backend` (with any
   `+qualifier` stripped), so every saved result already carries it and no schema
   changes.
3. It is in neither identity. It does not change what a key means, and adding a
   `KeySpace` field would change every fingerprint and make every saved result
   unreadable.
4. It is guarded where combining would be wrong: `AnalysisSample.from_result`
   carries it, and `require_common_structure_interpretation(...)` raises
   `IncompatibleStructureInterpretation` on a mixture. The module is
   `biotsr/domain/interpretation.py`, because the analysis layer may import
   `domain` and not `structure`. An undeclared backend is refused rather than
   assumed to match.

**Rationale.** Two readings of the same file can disagree about what an atom is
called while both fingerprints stay equal: on 1KX5 chain I with the manganese
ions named explicitly, reading atom names from columns `[12:16]` or `[13:16]`
turns `MN` into `N`, and 35 of 35 triangles carry different keys under identical
`KeySpace` and `SamplingProfile` fingerprints. With one interpretation shipped
the guard refuses nothing today, which is the only time it can be added without
retrofitting provenance onto results that never recorded it.

---

## ADR-080 — The package project is the repository root

**Status:** Accepted

**Decision.** The repository is an ordinary `src`-layout Python project:
`pyproject.toml`, `src/biotsr/` and `tests/` at the root, user documentation in
`docs/`, validation material in `validation/`. It installs with `pip install .`.

**Rationale.** A public package should not ask users to install from a
subdirectory whose name carries development history. No fingerprint includes a
path, so the layout is scientifically neutral.

---

## ADR-081 — Committed evidence is deterministic: timing and workstation paths are not evidence

**Status:** Accepted

**Decision.** Committed validation evidence carries scientific measurements and
no wall-clock timing; runtime observations are printed or written to untracked
files. Committed sidecars record input paths relative to the repository; an
input's identity is its `source.sha256`. Re-running a check on identical inputs
leaves the tree unchanged, which CI asserts. The package's own writer is
unchanged: a user's result records the path they gave.

**Rationale.** Timing values and absolute paths make every re-run "change" the
evidence on another machine, and an absolute path exposes a workstation account
name.

---

## ADR-082 — The historical narrow-slice metal-ion defect

**Status:** Historical: compatibility functionality removed in 1.0

**What the defect was.** A monatomic `MN` or `ZN` ion read through the narrow
historical atom-name slice (PDB columns 14–16) is spelled `N` and labelled
nitrogen (108) instead of manganese (100) or zinc (105), so every key of a
triangle containing it is chemically wrong. The `CarbohydratesTSR` lineage read
`HETATM` records that way.

**In 1.0.** Every atom is read with the specification's atom-name field,
columns 13–16 (ADR-108), and the namespace-aware resolver labels the ions as
their elements (ADR-031): manganese 100, zinc 105. Until 1.0.0rc5 the amino-acid
admissions still read `ATOM` records through the narrow slice, so a
non-conformant file writing an ion as an `ATOM` record (`ZN` from column 13)
had it read as `N`; since 1.0.0rc6 the defect is not reachable through any
method. The Cα test keeps columns 14–16, so calcium is never a Cα.

---

## ADR-083 — No licence until the provenance review was complete

**Status:** Superseded by ADR-085

**Summary.** A licensing and provenance review found a lexicon byte-identical to
files in repositories that carried no licence, a GPL-3.0 ancestor of the bin
tables inherited through `TSR-Package`, and clustering code that follows scipy's
implementation closely. This record withheld an open-source licence until those
questions were answered, and introduced the test that requires Apache-2.0 to be
claimed if and only if the canonical Apache text is present. ADR-085 records the
licence; the test remains.

---

## ADR-084 — Historical raw-float triplet text across platforms

**Status:** Historical: compatibility functionality removed in 1.0

**What it covered.** Historical triplet files print raw floating-point geometry
through `str(float)`, computed with the host C library's `acos` and `pow`. macOS
arm64 does not round these correctly on a measurable fraction of inputs while
glibc does, so the final digits of printed angles and lengths differed between
macOS and Linux, while keys, bins, mirror signs, categories, counts and both
fingerprints were identical. This record set a contract for comparing frozen
historical text: byte-exact on macOS arm64, and on Linux a measured tolerance on
the raw geometry columns only.

**In 1.0.** The package writes no historical text. What remains is a property of
the native format: keys files and sidecars carry integers and identities and are
identical on every platform, while the optional native triplet file prints raw
geometry whose final digits can differ by platform; the test suite compares
native triplet files byte for byte only against references produced on macOS
arm64. See [platform validation](../PLATFORM_VALIDATION.md).

---

## ADR-085 — BioTSR is licensed Apache-2.0; third-party material keeps its own terms

**Status:** Accepted (supersedes ADR-083)

**Decision.** BioTSR is licensed under the Apache License 2.0: the canonical
`LICENSE`, `license = { text = "Apache-2.0" }` with the OSI Apache classifier,
and `license: Apache-2.0` in `CITATION.cff`. `NOTICE` names the TSR Research
Group and reproduces the third-party notices the review identified — the MIT
notices of `TSR-Package` and `CarbohydratesTSR` and scipy's BSD-3-Clause notice —
with attribution for the historical repositories whose lexicon and constants
BioTSR ships. Third-party material is not relicensed. `NOTICE` carries the
notices so that setuptools' default licence-file patterns ship it in both the
wheel and the sdist. Details: [licensing](../LICENSING.md).

**Rationale.** The licence follows the ownership and provenance review recorded
in ADR-083. The test introduced there now holds on its "present" branch, and the
release-artifact checks require `LICENSE` and `NOTICE` in both artifacts.

---

## ADR-086 — joblib is a runtime dependency, used only for parallel batch execution

**Status:** Accepted (amends ADR-008)

**Decision.** `joblib>=1.3` is declared in `dependencies`. Batch execution is
serial by default (`--jobs 1`); parallelism is chosen explicitly. joblib is
imported only by `biotsr.cli.parallel`, inside a function, so the scientific core
imports nothing outside the standard library and `import biotsr` and `--jobs 1`
work where joblib is missing; asking for more than one worker there is an error
naming the installation command, never a silent fall-back. `biotsr[parallel]`
remains as an alias that adds nothing. Download concurrency
(`--download-workers`) uses standard-library threads.

**Rationale.** Parallel execution is a core capability for the datasets BioTSR
is used on, and an optional extra that users must know to request is a needless
obstacle. The import-hygiene tests allow joblib in that one module only, and the
release-artifact checks require `Requires-Dist: joblib>=1.3` and nothing else.

---

## ADR-087 — The public package implements standard conventions only

**Status:** Accepted

**Decision.**

* Each of the 20 methods runs under exactly one profile, its standard profile;
  `biotsr methods NAME` shows it. In Python, `profile=` (on `analyze`,
  `analyze_entities` and `analyze_cross`) accepts only `None` or that method's
  own default profile; a profile name, another method's profile or a derived
  profile raises `InvalidRequest`. The command line has no `--profile` and no
  `--format` option.
* Native output is the only output format, and it can be loaded back
  (`biotsr.load_result`, ADR-074, ADR-075).
* Arguments are checked per method. An argument the method does not read is
  refused by name, on the command line and in Python; `biotsr methods NAME` lists
  the method's target arguments. `--chain-2` is required for every cross method.
* The provenance field `legacy_compatibility_profile` is always false, and is
  kept so that the result schema is unchanged (ADR-078).

**Rationale.** With one convention per method, a result's meaning follows from
the method and the two recorded identities, and no reader has to ask which
profile produced it. Refusing an unread argument by name follows the package's
standing rule that an option which quietly does nothing is worse than one that
refuses. Where a standard convention differs from a historical program's, the
difference is stated as rationale in the records marked Historical and in the
methods pages.

---

## ADR-088 — One molecule per single-entity method

**Status:** Superseded by ADR-096 (1.0.0rc4)

**Decision.** The single-entity methods — `drug`, `carbohydrate`, `nucleotide`,
their `*_atom_selection` variants, and the entity side of
`drug_protein_cross`, `carbohydrate_protein_cross` and
`nucleotide_protein_cross` — compute over exactly one copy of the named entity.
`entity_id` may be omitted only when exactly one copy of the entity is in the
chain. Several copies are refused, listing their ids; an entity absent from the
chain is refused, naming the entities present; an explicit id that does not exist
is refused. The copy used is recorded as the result's `target` (for example
`BEN_1`), also when the id was omitted. The pooled methods (`carbohydrate_chain`,
`nucleotide_chain`, and `--entity-names` on cross methods) pool by design.

**Rationale.** Copies of a ligand, monosaccharide or nucleotide are separate
molecules. Pooling them silently would be a different calculation — triangles
spanning two molecules — under the name of a single-molecule one. Requiring the
id when it is ambiguous, and recording the copy used, keeps every result
attributable to one molecule.

---

## ADR-089 — Empty and partially matched selections

**Status:** Accepted (amended by ADR-101: a calculation that lacks something it
named is reported `partial` in a batch; and by ADR-106: such a calculation
exits 4 and every analysis asks for it to be acknowledged)

**Decision.**

* An explicit selection that matches nothing — atoms, residue types, specific
  residues, Cα residue numbers, entity atoms, an `entity_names` override, an
  empty cross side — is refused with a message naming the entity or the items.
* When some requested items are missing, the result carries a warning beginning
  `requested but not found: `, naming them; the CLI prints it and records it
  under the result's `"warnings"` in `run.json`.
* A selection with fewer than three points produces a warning beginning
  `no triangles: `, as does a cross calculation in which no mixed-side triangle
  can be formed. In a batch that row's status is `no_triangles`: it is not a
  failure and does not affect the exit code.
* Completeness rules (exact atom counts, ADR-004) are unchanged, and their
  exclusions are reported.

**Rationale.** An empty result from a mistyped selection is indistinguishable
from a genuine finding unless the software says which it is. Refusing what
matches nothing, and naming what is missing from what partly matches, turns a
silent difference in the triangle set into a visible one. A selection too small
to form a triangle is a legitimate, reportable outcome, not an error.

---

## ADR-090 — Every analysis requires one key space and one sampling rule, compared from structured records

**Status:** Superseded by ADR-094 (1.0.0rc4); previously amended for 1.0.0rc3 (supersedes ADR-055)

**Decision.** Every operation that compares samples — `feature_matrix` (and so
similarity, clustering and the clustermap), `common_key_intersection`,
`common_key_intersection_by_group`, `key_occurrence_percent`,
`key_occurrence_percent_one_vs_all` and `biotsr analyze common-keys` — applies
one function, `analysis/samples.py::check_comparable`:

* one key space, always (ADR-054); no option relaxes it;
* one sampling **rule**, unless the caller passes `allow_mixed_sampling=True`
  (`--allow-mixed-sampling`), in which case every result carries a warning
  naming the differing fields.

The rule (`comparison_rule`) is built from structured fields only:

* the method;
* every `SamplingProfile` field except `side_a_chain`, `side_b_chain` and the
  two free-text side descriptions;
* per side, from the selection record each result carries
  (`provenance.selection_request`, `domain/selection_record.py`): the role, the
  selector, the type-level vocabularies (`residue_types`, `entity_names`, as
  upper-case sets) and which kinds of explicit enumeration narrowed the side.

Targets — structure, chain, entity name and copy of a single-entity method,
residue numbers, atom tokens — are recorded per sample in every report and never
decide comparability. A result without a selection record (written before
1.0.0rc3) is compared only with the opt-in, unless its selector takes no
selection argument. Each report records, per sample, both fingerprints, the rule
fingerprint, the sampling settings, the selection rule and targets, and the
warnings of the sample's own result.

**Amendment (1.0.0rc3).** The rc2 rule reduced a cross side to "role and
selector kind" by deleting every `name=value` token from its free-text
description, which also deleted a ligand's atom subset, a protein side's
residue types and a pooled entity vocabulary; and one-sided methods recorded no
selection argument in their sampling profile at all, so `aa_ind` over LYS and
over TYR, or `carbohydrate_chain` with two different vocabularies, compared as
one rule. Feature matrices, similarity and clustering applied no sampling rule
(ADR-055). All three are corrected; `SamplingProfile.rule_fingerprint` and
`rule_dict` are removed, and the `SamplingProfile` itself, with its
fingerprint, is unchanged.

**Rationale.** A key-frequency comparison is a statement about how often a
triangle type occurs, and a common key that it occurs at all. If the samples
were drawn by different rules — a size filter, two contact cutoffs, two
residue-type filters, the whole ligand against a chosen part of it — the
difference can be an artefact of the sampling rather than of the structures.
Comparing the same rule across different structures, chains and molecules is
what a dataset is for, so those are targets. Free text is never parsed for the
decision because a description is not a contract: the rc2 parse shows how
easily it drops what matters. Mixing rules remains possible when a study
intends it, and is then stated on every result.

---

## ADR-091 — Sample ids are names inside the output root, checked against the file system

**Status:** Superseded by ADR-098 (1.0.0rc4): a sample id is now an output prefix, checked by the same rules

**Decision.** A sample id is one directory name directly inside `--output`: no
path separators, no absolute path, not `.` or `..`, no leading dot, no control
characters, not one of the names the output root uses for itself (`analysis`,
`batch_summary.csv`, `batch_failures.csv`, `fetch_report.csv`), and at most 200
bytes. A symbolic link at a sample's path is never written through, replaced or
deleted. The checks run during planning, so `--dry-run` reports an invalid id
without writing anything. An invalid id fails its row in a batch and exits 2 for
a single run.

**Rationale.** A sample id is joined onto the output root, and `--if-exists
overwrite` replaces — that is, deletes — the directory it names. An id that
climbs out of the root, names an absolute path, nests, collides with the root's
own reports or resolves through a symbolic link would let a run delete something
that is not its own. String checks alone cannot see a link, so the final path is
checked against the file system as well.

---

## ADR-092 — Interruption and worker-failure reporting

**Status:** Accepted; staging names amended by ADR-098 (1.0.0rc4)

**Decision.**

* SIGTERM is handled like Ctrl-C by every command (exit 130).
* A batch interrupted at any point after planning — while downloading, while
  checking existing outputs or while computing — writes its summary, with rows
  `interrupted` or `not_started`; when interrupted during downloads it also writes
  `fetch_report.csv`. `biotsr fetch` writes its report when interrupted.
* If a worker process dies (for example, out of memory), finished rows are kept,
  rows that were running are `failed` with a message naming the pool failure,
  rows never started are `not_started`, and the batch exits 3; rerunning with
  `--if-exists skip` resumes.
* Limits: SIGKILL and power loss cannot be caught. Their staging directories
  (`.<sample>.partial-<pid>-<host>-<token>`) are removed by a later run on the
  same machine once the owning process is gone; ones from another machine are
  left and reported. After Ctrl-C no new download attempt starts; a download in
  flight is stopped at once (amended in 1.0.0rc3, below).

*Amendment (1.0.0rc3).* The original text said a download in flight "ends at
its next network read (bounded by `--timeout`)". That was not a bound: a
response that keeps arriving in small pieces never waits `--timeout`, and the
body was read in one blocking call. Now a connection being made is polled for
cancellation, every connected socket is registered and shut down on
cancellation, and the body is read in pieces with a cancellation check between
them; only a host-name lookup and a TLS handshake step in progress are waited
for. `--timeout` limits each wait and is documented as such; there is no total
time limit on a download.

**Rationale.** A scheduler's time limit sends SIGTERM, and a long batch is often
stopped this way. A summary written at any interruption point means the record of
what was done is never lost with the run, and a resumable state is always
reachable. Staging directories carry the process id and host so that cleanup
never removes another live run's work.

---

## ADR-093 — CSV manifests refuse what they cannot read

**Status:** Accepted

**Decision.** A manifest is read as UTF-8, with or without a byte-order mark.
It is refused, with the line and the remedy, when it has non-blank cells beyond
the header or under an empty header (separate list items with `;` or quote the
cell); when header names are duplicated; and when it has a `profile` column
(drop it, or pass `--ignore-column profile`; every method runs under its
standard profile, ADR-087). Aliases, `--column`/`--keep-column`/`--ignore-column`
and the chain rules are unchanged (ADR-070); a chain is never inferred.

**Rationale.** An unquoted list such as `TYR,PHE` in a comma-separated file
spills into cells that have no column, and a duplicated header makes one value
shadow another; reading either would bind values to the wrong parameters without
a word. A `profile` column states a choice the package does not offer, so
accepting it silently would misdescribe what ran.

---

## ADR-094 — Comparability is one key definition; processing differences are acknowledged; selections never block

**Status:** Accepted (1.0.0rc4; supersedes ADR-090; amended by ADR-105: the
processing record includes the modified-residue rule; and by ADR-106: partial
samples need their own acknowledgement)

**Decision.** Every operation that compares samples — feature matrices (and so
similarity, clustering and the clustermap), the four common-key operations and
`biotsr analyze common-keys`, on live or loaded results — applies
`analysis/samples.py::check_comparable`:

1. **Key definitions, always.** The samples must share one `KeySpace`: the same
   number means the same triangle — point and label conventions, theta and
   distance bins (by content fingerprint), canonical ordering, key formula,
   mirror convention and numeric conventions. Comparability is decided from
   the key-space definition, never from method names: methods that share a key
   space compare (every atomic method labels from the 112-code lexicon in 29 ×
   18 bins), and the same method with `mirror` or `aa_grouping` does not. No
   option relaxes this (ADR-054).
2. **Processing, separately.** The `SamplingProfile` fields that say how
   triangles are drawn from whatever was selected — candidate strategy,
   completeness and composition rules, size filter, degenerate policy, sequence
   label policy, contact cutoff, contact filter and rule, contact projection,
   pool identity and order — and the alternate-location rule the result
   records must agree. When they differ the comparison is refused
   (`IncompatibleProcessing`, a subclass of `IncompatibleSampling`) unless the
   caller acknowledges it (`allow_mixed_processing=True`,
   `--allow-mixed-processing`); every result then carries a warning naming each
   differing setting.
3. **Selections never block.** Structures, chains, amino-acid types (LYS-only
   against TYR-only), ligands, carbohydrates, nucleotides, residue numbers, atom
   subsets, one copy against all copies, and the method name are what a dataset
   is made of. They are recorded per sample — the selection record
   (`provenance.selection_request`), the target and the `SamplingProfile`
   target fields (selector, scope, cross method, side chains and side
   descriptions) — and never decide comparability. Samples need not share a
   key; an empty intersection is a result.

Each report records, per sample, the method, the structure, the key-space,
full sampling and processing fingerprints, the processing settings, the
selection per side, and the sample's own warnings; and for the set, the
processing fingerprints and which samples share each, the differing settings,
the methods and the selection scopes.

**Fingerprints.** Neither identity fingerprint changed; the 42 frozen reference
calculations pin both. `KeySpace` hashes each scheme's name as well as its
content fingerprint, so two schemes identical in content but differently named
would be refused — a conservative error, never a wrong acceptance; no two
schemes in this release differ that way. The full `SamplingProfile`
fingerprint of a cross result includes the side chains and the side
descriptions (which name the ligand or residues), so for cross methods it
differs whenever the targets do, and it differs between selectors and scopes;
it remains an exact record of the run and is reported, and comparability is
decided by the processing fingerprint instead.

**Rationale.** A comparison needs the same key
meaning, and a processing difference changes what the frequencies measure, so
it must be stated; what was selected is the subject of a study, not a reason to
refuse it. rc3 treated residue-type sets, entity vocabularies, atom-subset
narrowing and the method as part of a "sampling rule" and refused, for
example, LYS-only against TYR-only without an opt-in.

---

## ADR-095 — One alternate-location rule: blank or A, no fallback

**Status:** Accepted (1.0.0rc4; supersedes the `amino_acid_atomic` admission
part of ADR-004 and ADR-034)

**Decision.** Every method admits an atom whose altLoc indicator is blank or
`A`, and excludes `B`, `C` and every later indicator. There is no fallback: a
residue or molecule whose atoms exist only as `B` (or later) conformers
contributes none of them, and a method with a completeness rule then excludes
it with the reason recorded. Excluded atoms leave skip records
(`altloc_rejected`) and the result a warning counting them by indicator; every
result records `source.altloc_policy = "blank_or_a"`, and the processing
comparison of ADR-094 includes it. To study another conformation, prepare a
structure file whose indicators select it. `amino_acid_atomic` keeps its
lineage's element token and completeness table (ADR-004).

**Change.** `amino_acid_atomic` previously excluded only `B` (its lineage's
`line[16:20]` test) and admitted `C` and later conformers: a residue with three
side-chain conformers had extra atoms and was excluded as incomplete, and a
residue whose side chain existed only as `B` and `C` was built from `C`. On a
synthetic leucine pair this is the before/after: rc3 excluded LEU 2 (blank +
A + C = 12 atoms) and computed LEU 3 from conformer C; rc4 computes LEU 2 from
blank + A (8 atoms) and excludes LEU 3 (4 atoms), naming the rule. No test
structure carries an altLoc beyond `B`, so no frozen reference changed.

**Rationale.** One rule across methods makes results from different methods
describe the same conformation, which the comparison rule now relies on; a
silent fallback to another conformer would mix conformations without a trace.

---

## ADR-096 — Molecule instances: an identifier names one copy; a name alone is every copy on the chain

**Status:** Accepted (1.0.0rc4; supersedes ADR-088)

**Decision.** For `drug`, `carbohydrate`, `nucleotide`, their `*_atom_selection`
variants and the entity side of the three entity cross methods:

* the chain is required and only that chain is searched; a chain is never
  inferred and chains are never combined;
* chain + name + identifier selects that one copy (scope `one_instance`);
* chain + name without an identifier selects every copy of that name on that
  chain for **one combined calculation** (scope `all_copies`), even when there
  is one copy; triangles may span copies, by design. All selected copies are
  recorded (`instances`, `NAME_ID` tokens), the target reads `NAME_allcopies`,
  a cross side reads `entity=NAME(all copies)`, and the CLI prefix ends in
  `NAME_allcopies`. The completeness rule applies per copy; a declined copy is
  recorded and warned about;
* an atom list applies to each selected copy; requested atoms a copy lacks are
  reported per copy.

Identifiers are PDB residue identifiers: the residue number as written (it may
be negative) with at most one insertion-code letter (`404`, `52A`, `-5`).
Leading zeros are dropped; `52` never matches `52A`; positions ("the second
copy") are refused, as is any other spelling; an identifier whose records
appear in two separate places on the chain is refused. All copies of one sugar
name is not the whole glycan: `carbohydrate_chain` (and `entity_names`) pool
every residue of a vocabulary, different sugars included.

For amino acids, an occurrence needs its identifier: `amino_acid_atomic` with
`residue` and `residue_name` is that residue; with `residue_name` alone it runs
one calculation per residue of that name (its per-residue scope is kept);
`aa_ind_specific` names residues as `NAME_ID` or `NAMEID` (letters-only names),
both normalised, and refuses any other token; `aa_ind` pools residue types as
before. Whole-chain, per-residue, pooled and cross definitions are not merged.

**Rationale.** rc3 refused a name
without an identifier when a chain held several copies, which made "every NAG
of chain D" unreachable; recording the scope and the instances keeps every
result attributable, and the `all_copies` scope is visible in the result, the
record and the file name, so it is never mistaken for a single-copy
calculation.

---

## ADR-097 — Methods are chosen when BioTSR is invoked; wide sheets are read method-aware

**Status:** Accepted (1.0.0rc4; extends ADR-070; cross narrowing amended by
ADR-104, merging of identical requests by ADR-103)

**Decision.** `biotsr batch --method NAME` is repeatable: every source row
becomes one request per method, and each request keeps its source row(s) and
method in its completion record, the summary and the analysis. A
one-request-per-row manifest may still name methods in a `method` column, but
not together with `--method` (refused at planning, naming both remedies).

A wide sheet (`--sheet tsr` for the TSR dataset-sheet columns, `--sheet custom`
for explicit mappings only) is read **method-aware**: per method and
participant a layout says which column supplies each argument, so each method
reads its own participant columns — `carbohydrate` reads `carb_chain`,
`carb_name`, `carb_id`; `aa_ind` reads `aa_chain` and `aa_ind`, or with
`--participant 2` `aa_chain_2` and `aa_ind_2`; a cross method's single
participant covers both sides (`carbohydrate_protein_cross`: `carb_chain`
with `aa_chain`). Cross narrowing columns are read only when chosen
explicitly (ADR-104).
`--map METHOD[.PARTICIPANT].FIELD=HEADER` and `--map-file` add or override
entries; unknown methods, fields or columns, a field mapped twice, and both
spellings of one column (`car_atom_selected`, `carb_atom_selected`) are
refused. A populated column never means "run another method". Columns no
selected method reads are reported as unused. A row that does not describe a
method's calculation fails, or with `--inapplicable skip` is recorded as
`inapplicable`; a row that describes it with a blank chain fails, naming the
column. Rows requesting the same calculation of the same sample share one
(ADR-103). `--dry-run` lists every
request with its source rows, method, participant, chains, selection, output
prefix and outcome. Canonical argument names (`chain`, `chain_2`, …) are
unchanged; `aa_group` stays an annotation.

**Rationale.** A wide sheet describes structures, not runs; which methods to
run is the analyst's choice for this invocation. Collapsing the chain columns
into one alias (ADR-070) made a wide row ambiguous; a per-method layout removes
the ambiguity without renaming the API to one dataset's headers.

---

## ADR-098 — One shared output folder; one prefix per calculation; completion recorded last

**Status:** Accepted (1.0.0rc4; supersedes ADR-071 and ADR-091; amended by
ADR-105: a completion record also holds the processing rules, and
`--if-exists skip` compares them)

**Decision.** Every calculation writes `PREFIX.keys`, `PREFIX.meta.json`,
`PREFIX.triplets` (unless keys-only) and its completion record
`PREFIX.run.json` directly into `--output`. The prefix names the structure,
the chain(s), the method, the target (`NAG_1`, `NAG_allcopies`, residue types,
named residues or atoms, each cross side's narrowing) and the modifiers and
settings that change the numbers (`biotsr/io/naming.py`); when the readable
form is inexact or long, a deterministic request fingerprint is appended.
`sample_id` overrides it, checked by the former sample-id rules. The bin
counts left the file names (`.keys_29_35`); the metadata keeps the full bin
definitions. File contents are unchanged: tab-separated UTF-8 text with one
header line.

A calculation is computed into a staging directory private to the invocation,
then installed under a per-prefix lock: an intent record (`"complete":
false`), each file moved into place, files of an older version of the same
calculation no longer produced removed, and the record rewritten complete with
every file's size. Planning refuses two different requests with one prefix
(case-insensitively, and a prefix extending a per-residue sweep's); a file at a
target name that no record of that prefix lists is never replaced or deleted,
even with `--if-exists overwrite`; links are never written through. `--if-exists
skip` keeps a calculation only when its record is complete, records exactly the
request, covers the requested files and every file has its recorded size, and
recomputes one an interrupted run left incomplete. Readers (`analyze
--from-results`, `describe`) read completion records, refuse incomplete ones
and any key file no record lists. The result schema is 2 and the record
schema 4; outputs of 1.0.0rc3 and earlier (per-sample directories, names with
bin counts, schema 1) are refused with a message saying so, and no legacy
reader is kept.

**Rationale.** One folder per study is the layout analysts use; a directory per
calculation added nothing but nesting. Without a directory to rename
atomically, completion is what the record says, and the record is written
last, so nothing partial is ever accepted.

---

## ADR-099 — The default hierarchical analysis needs no cluster count

**Status:** Accepted (1.0.0rc4; amends ADR-060)

**Decision.** The default workflow is distance matrix → average-linkage tree
→ dendrogram and ordered heat map (`hierarchy`, `clustermap_spec(matrix)`,
`render_clustermap`); it requires neither `n_clusters` nor `maxclust`, and an
uncut tree is not reported as a number of groups. A flat partition is
optional and needs an explicit rule: `cut_hierarchy(tree, n_clusters=t)`
(`maxclust`, ADR-060) or `cut_hierarchy(tree, height=h)` (every merge above
*h* cut); `clustermap_spec(..., n_clusters=t)` or `(..., height=h)` adds the
cluster stripe. The adjusted Rand index is computed only for a partition and
fully labelled samples. Observed class labels and algorithm clusters are
separate stripes, labelled "observed class" and "cluster". Nothing infers a
cluster count, least of all from the known labels.

**Rationale.** Requiring a count to draw the tree forced a choice the analysis
does not need, and invited choosing it from the labels, which makes agreement
with them a construction.

---

## ADR-100 — Redistribution of the historical tables is confirmed by the project owner

**Status:** Accepted (1.0.0rc4)

**Decision.** The project owner has confirmed proceeding with redistribution of
the tables BioTSR carries from the group's historical repositories — the
amino-acid and nucleotide exact-atom-count tables, the nucleotide residue
vocabulary, the bin and label tables, and the standard atom lexicon — as part
of BioTSR under Apache-2.0 (ADR-085). Their attribution, licence notices and
provenance records are preserved unchanged (`NOTICE`, `docs/LICENSING.md`);
the confirmation is the owner's and is recorded as such, and no third-party
permission document is claimed.

**Rationale.** The decision belongs to the project owner and is stated where
the licensing record is kept.

---

## ADR-101 — Modified amino acids: read where defined, named where not

**Status:** Accepted (1.0.0rc4; amended by ADR-105: the rule is recorded with
every result, by ADR-106: a partial calculation exits 4, and by ADR-107: the
Cα methods read and label the supported modified residues)

**Decision.** BioTSR recognises exactly three modified amino acids:
phosphothreonine `TPO`, phosphoserine `SEP` and phosphotyrosine `PTR`. They are
the residues the historical amino-acid package gives exact heavy-atom counts
for (11, 10 and 16, in the completeness table BioTSR already carries), the
historical cross scripts list as protein residues, and the shared atom lexicon
labels (`P`, `O1P`–`O3P`, `OG1`, `OH`, …). RCSB files write them as `HETATM`
records, which the historical readers never reached (they read `ATOM` lines
only). The amino-acid atomic admissions (`amino_acid_atomic`, `aa_ind`,
`aa_ind_mix`, `aa_ind_specific`, `amino_acid_atom_selection`) and the protein
sides of the atomic cross methods now admit `HETATM` records **of exactly these
names** (`AtomAdmission.hetero_residue_names`); every other `HETATM` record
stays excluded, the alternate-location rule (ADR-095), chain, identifiers,
truncation and completeness rules apply to them unchanged, and nothing is
converted to its unmodified residue.

Every other modified residue stored as `HETATM` — a residue carrying the
backbone atoms `N`, a carbon `CA` and `C`, such as `CME` or `MSE` — has no
atom count and is **not supported**. The Calpha methods label residues in the
20-code amino-acid space, which has no code for any modified residue, so they
represent none of them (`not represented`), TPO, SEP and PTR included; the
interface rule of `protein_protein_calpha_cross` keeps reading `ATOM` records
only. Such residues are never dropped silently: a request naming only them is
refused with the reason; otherwise the result carries a warning beginning
`not supported:` or `not represented:` that names them, and the side's
selection record lists them under `unsupported` (and what the structure does not
contain under `not_found`). A named residue the completeness rule declines is
reported as `requested but excluded by the method's completeness rule`, no
longer as not found. In a batch, a calculation that lacks anything it named —
not found, not supported, not represented, or excluded — has the status
`partial` (written like `ok`, counted separately; ADR-106 later gave it exit
code 4).

**Unresolved scientific choice.** Representing a modified residue in a Calpha
method would need a label for it in the 20-code space: either the unmodified
residue's label (TPO as THR), which silently converts it, or a new code. The
20 codes (4–23) already fill radix 20, one per remainder modulo 20, so a 21st
code either aliases two label triples onto one key or needs a larger radix,
which changes the value of every Calpha key. Neither is defined; the choice
belongs to the project owner and nothing is assumed.

**Fingerprints.** Unchanged. Admission names, the altLoc rule and the
completeness policy are those of rc4; the 42 frozen reference calculations
contain no modified residue and are unchanged. Results for structures with
`HETATM` TPO, SEP or PTR change: those residues are now part of the atomic
selections that include them.

---

## ADR-102 — One grammar for residue identifiers and ranges

**Status:** Accepted (1.0.0rc4)

**Decision.** A residue identifier is the residue's PDB number with an
optional one-letter insertion code (`48`, `52A`, `-5`; leading zeros dropped,
the insertion code matched exactly), in every place one is written: `residue`,
`entity_id`, `NAME_ID` residue references (`LYS_52A`, `LYS52A`), atom tokens
(`LYS_52A_NZ`) and `aa_selection` items. An `aa_selection` item may also be a
range of plain numbers, ascending, both ends included (`10-12`, `-5--3`),
expanding to plain numbers only. Refused with a message naming the item,
instead of being kept verbatim or skipped: a range with an insertion code at
either end (`27A-29`), a descending range (`12-10`), any other malformed item,
and a range that spans a residue with an insertion code in the structure
(`9-11` over `10A`) unless that residue is also listed (`9-11;10A`); to leave
it out, list the residues one by one. List order never matters: points follow
the structure file. The command line, a manifest and the Python API share the
parser. The selection record's `instances` for a Calpha list are the residues
actually used, and the requested ones it lacks are `excluded`.

**Rationale.** A range kept verbatim matched nothing and a descending range
expanded to nothing; both changed the triangle set silently. Whether a numeric
range includes insertion-code residues has no single obvious answer, so it is
asked rather than assumed.

---

## ADR-103 — A duplicate calculation is computed once; a sample is never merged away

**Status:** Accepted (1.0.0rc4; amends ADR-097)

**Decision.** A batch merges two requests only when they are the same
calculation **of the same sample**: method, structure, chains, selection,
modifiers and sample id (an explicit `sample_id`, or the automatic prefix,
which names exactly the calculation) all equal. The merged request lists every
source row and participant, and the summary says it serves several source
requests. Rows with distinct sample ids are distinct samples: each is computed
and written under its own prefix, and an analysis reads both. Requests that
would merge but carry different labels are refused at planning (one output
carries one label); a kept metadata column whose values differ records each
row's value (`line 11=x; line 12=y`). The output-prefix collision checks
(ADR-098) are unchanged.

**Rationale.** In a wide sheet one calculation often appears on several rows
(the protein of several ligand rows); computing it once is correct. A
different sample id or label, however, says the rows are different samples,
and silently choosing one row's label would change the analysis.

---

## ADR-104 — Cross sides from a wide sheet are chosen explicitly

**Status:** Accepted (1.0.0rc4; amends ADR-097)

**Decision.** The `tsr` layout reads a cross method's chains and molecule
name and identifier from the row. Each side of a cross method also has
candidate narrowing columns from the column dictionary
(`biotsr.cli.sheet.CROSS_SIDE_COLUMNS`): the protein side of the
molecule-protein methods `aa_ind`, `aa_ind_mix`, `aa_ind_specific`,
`aa_atom_selection`; the molecule side its atom-subset column
(`drug_atom_selected`, `car_atom_selected`, `nt_atom_selected`); side 1 and
side 2 of `protein_protein_cross` the `aa_*` and `aa_*_2` columns; and of
`protein_protein_calpha_cross` `aa_selection` and `aa_selection_2`. When a
row populates any candidate of a side, the side must be chosen —
`--cross-side METHOD.SIDE=COLUMN` (or `cross_sides=` in Python) narrows it by
that column, `=whole` uses the whole side, and a `--map` of the side's field
counts as the choice — otherwise the run is refused, naming the columns and
rows. A row whose chosen column is blank is inapplicable, never computed on
the whole side. The populated columns a choice leaves unused are named per
request (`not applied: …`) in the preview and the summary, and the preview and
the summary's `selection` show both sides' resolved selections.

**Rationale.** The sheet's narrowing columns were written for one-sided
methods; for a cross calculation they imply different selections and the sheet
does not say which one is meant. Reading none of them silently computed whole
sides while the sheet looked narrowed; reading one by precedence would invent
a rule. An explicit choice per side is the only reading that is both safe and
reproducible.

---

## ADR-105 — Processing rules are part of a calculation's identity

**Status:** Accepted (1.0.0rc4; amends ADR-094 and ADR-098; amended by ADR-107:
the modified-residue rule is recorded per label domain, and resume also
compares the key space)

**Decision.** Two identities are kept apart. The **key space** says what a key
number means (labels, bins, ordering, encoding); it is unchanged by this
decision, and two results whose key definitions agree still share it. The
**processing rules** say which atoms were admitted and how they were selected;
they are recorded with every calculation and compared wherever a saved result
is reused or compared:

* `selection.base.processing_identity()` is
  `{"altloc_policy": "blank_or_a", "modified_residue_policy":
  "hetatm_amino_acids:PTR,SEP,TPO"}` — the alternate-location rule (ADR-095)
  and the modified-residue rule (ADR-101), the latter derived from
  `completeness.MODIFIED_AMINO_ACIDS`. One value for every method, like the
  alternate-location rule.
* Every completion record's request holds it (`request.processing`), and
  `--if-exists skip` keeps an existing calculation only when the recorded
  rules equal the current ones. A record without them is never taken to be
  current; it is refused with the reason and recomputed only with
  `--if-exists overwrite`. A record of different rules is refused the same
  way.
* Every result's provenance records the modified-residue rule
  (`source.modified_residue_policy`), and the processing policy every
  analysis compares (`analysis.samples.processing_policy`) includes it beside
  the alternate-location rule. A result written before the field existed
  reads as *not recorded*; comparing it with a current result needs
  `allow_mixed_processing=True`, and the report names the difference.

The package version is not used to decide any of this: one version number
(1.0.0rc4) covered both the rule that excluded `HETATM` TPO, SEP and PTR and
the rule that reads them.

**Evidence.** `tests/data/legacy_outputs/before_adr101_hetatm_modified_residues/` is an unedited output of
an earlier development build that also reported `1.0.0rc4`: `aa_ind_specific`
of `SER_4;TPO_3;TYR_6` on the synthetic `SYNPHOS.pdb`, 18 atoms and 816
triangles, TPO 3 excluded. Its request differs from the one the current
version builds only in `processing`, so before this decision `--if-exists skip`
kept it; the current rule computes 29 atoms and 3654 triangles for the same
request. `tests/integration/test_processing_identity.py` checks that skip
refuses it and leaves it untouched, that overwrite recomputes it, and that
analysis keeps its key space but names the processing difference.

**Fingerprints.** Key-space fingerprints are unchanged. `SamplingProfile`
fingerprints are unchanged: the rule is recorded with the source, as the
alternate-location rule is, not in the profile. The analysis *processing*
fingerprint changes for every sample, because the compared policy has a new
field; that is the point of the change, not a side effect to hide.

**Rationale.** A resumed run that keeps an output computed under other
admission rules mixes two rule sets in one dataset without any sign of it, and
the numbers differ for exactly the structures the rule change concerns. Only a
record of the rules themselves can tell them apart.

---

## ADR-106 — Partial results are explicit through every stage

**Status:** Accepted (1.0.0rc4; amends ADR-089, ADR-094 and ADR-101)

**Decision.** A result is **partial** when its selection did not get everything
it named: a requested residue, atom, residue type or molecule the structure
does not contain (`not_found`); a modified residue the method cannot read
(`unsupported`, including every modified residue in a Calpha method); a named
residue or molecule copy the completeness rule rejected (`declined`); or
requested atoms a selected copy lacks (`atom_subset.missing`). A residue
rejected among all residues of a type, or a copy rejected among all copies, is
recorded (skip records, `excluded`) but does not make the result partial:
nothing that was named is missing.

* **Structured, not worded.** The selection record lists each kind per side;
  `selection_record.selection_omissions` derives the omissions from it, and
  from the result's warnings only for a result written before those fields
  existed. Every result entry of a completion record carries `omissions`
  (`{}` when complete), and the batch status is decided from the record.
* **Batch and resume.** The status is `partial` (ADR-101). A request whose
  verified existing output is kept by `--if-exists skip` stays `skipped`, but
  when that output is partial its message says so (`kept partial result: …`)
  and the summary line counts it.
* **Exit codes.** A new code, **4** (`EXIT_INCOMPLETE`): everything finished
  and nothing failed, but at least one result is partial (or a resumed run kept
  a partial one, or a row with no triangles also lacked something it named). `biotsr run`, `biotsr batch` and `biotsr analyze` return it;
  3 (a failure) takes precedence. Until now a partial result exited 0, so
  automation could not tell it from complete success.
* **Analysis.** `check_comparable` refuses a set containing partial samples
  (`PartialSamples`, naming each and what it lacks) unless
  `allow_partial=True` (`--allow-partial`). That acknowledgement is separate
  from `allow_mixed_processing`; neither implies the other, and neither relaxes
  the key-space rule. With it, nothing is dropped: the comparison records
  `partial_allowed` and every partial sample's omissions, each sample's
  provenance carries `partial` and `omissions`, and a warning names them —
  in feature matrices, similarity and clustering exports, and common-key
  reports alike. `biotsr analyze` then exits 4.

What was selected still never blocks a comparison: a complete LYS-only sample
and a complete TYR-only sample compare without any acknowledgement (ADR-094).

**Rationale.** A partial result is a legitimate, usable outcome, but it is not
the selection that was asked for. Reported only in a batch summary, it reached
an analysis looking exactly like a complete one. Making the acknowledgement
explicit, recorded and separate from the processing opt-in keeps both
decisions visible in the report.

---

## ADR-107 — Residue labels for modified residues, at radix 30

**Status:** Accepted (1.0.0rc5; amends ADR-015, ADR-032, ADR-052, ADR-101 and
ADR-105)

**Decision.** The protein Cα methods (`protein`, `aa_selection`,
`protein_protein_calpha_cross`) label residues from the TSR Research Group's
residue tables, which add eleven modified residues to the 20 regular amino
acids, and encode both tables at **radix 30**. The method names are unchanged;
there is no second mode.

| scheme | residues | codes | radix |
|---|---|---|---|
| `amino_acid_30` (default) | 20 regular (4..23); TPO 24, SEP 25, PTR 26, NIY 27, M3L 28, TRF 29, CSO 30, LLP 31, CME 32, CSD 33; MSE = MET's 13 | 30, over 4..33 | 30 |
| `amino_acid_group_23` (`aa_grouping`) | the 13 historical groups of the 20 (4..16); MSE in MET's group (14); each other modified residue on its own, with the same code as above | 23, over 4..33 | 30 |

The two schemes share a radix and remain two key spaces: the scheme name and
content fingerprint enter `KeySpace`. Grouping keeps its meaning — residues of
one group share one code, so their triangles share keys on purpose.

* **Why 30.** The encoder writes each label as the digit `code − 1` in base
  `n_labels`. With the 30 codes at radix 20, codes above 23 carry into the
  neighbouring digit: `(SEP, PHE, PHE) = (25, 4, 4)` and
  `(TPO, TPO, PHE) = (24, 24, 4)` both give label component 9663, so in the
  same angle and distance classes two different residue combinations got one
  key. Of the 4,960 ordered code triples, 255 label components were shared at
  radix 20; of the grouped table's 2,300, 494 at radix 13, and still 98 at
  radix 23, its number of codes, because its codes skip 17..23. At radix 30
  every code lies in one window of 30 consecutive integers, and all 4,960 and
  2,300 triples have distinct label components (checked exhaustively in
  `tests/unit/test_encoding.py`).
* **The window rule replaces the count rule (amends ADR-015).** `LabelScheme`
  used to require the number of distinct codes to equal the radix, with an
  override for the atom lexicon. That test is neither necessary (109 atom codes
  at radix 112 are collision-free) nor sufficient (23 grouped codes at radix 23
  collide). A scheme is now refused unless every code is an integer ≥ 1 and all
  codes fit in one window of `n_labels` consecutive integers — the condition
  under which no label digit can carry. The override
  (`allow_radix_mismatch`) is gone. `labels.schemes.label_triple_collisions`
  enumerates every canonically ordered triple through the encoder for anyone
  checking a table.
* **The encoder checks labels.** `KeyFormulaSpec.label_range` carries the
  scheme's lowest and highest code (every profile supplies it), and `encode`
  refuses a label outside it instead of producing a key that aliases another
  triple. It is not part of `KeySpace`: it adds a check and changes no key.
* **Selection (amends ADR-101 for the Cα methods).** The Cα selectors read the
  supported modified residues — an `ATOM` record of them like any other, and a
  `HETATM` Cα when the residue name is one of the eleven and the CA atom's
  element field is carbon or blank. No other atom of the residue is required
  — the calculation uses the Cα alone, and an `ATOM` record never needed one
  either. (An earlier 1.0.0rc5 build also required `N` and `C`; that
  requirement was withdrawn, because it dropped residues without reporting
  them and added no protection: a free amino acid has `N` and `C` too.) No
  other `HETATM` record is a Cα point: not a ligand, not an ion with an atom
  named `CA`, not a modified residue the table does not cover (still reported
  *not represented*, and the result is partial, as before). A supported
  residue present without a readable Cα — no `CA` atom, a `CA` only as altLoc
  `B` or later, or a `HETATM` `CA` declaring another element — is named with
  that reason and makes the result partial; a residue `aa_selection` names
  that the chain contains is never reported as not found. `ATOM` records are
  not element-checked, for any residue; that asymmetry is recorded, not
  changed. The chain,
  residue-instance, truncation and blank-or-`A` alternate-location rules are
  unchanged. `protein_protein_calpha_cross` decides interface residues with the
  heavy atoms of the same residues, so a modified residue that has a Cα point
  takes part in the interface like any other; one without a readable Cα is an
  interface residue without a Cα, recorded as such (ADR-052) and named.
* **The atomic methods are unchanged.** They label atoms from the atom lexicon
  and read TPO, SEP and PTR only, the residues with defined atom counts; the
  residue table is never used to decide atomic completeness. Their key spaces,
  bins and results are unchanged.
* **The modified-residue rule is recorded per label domain (amends ADR-105).**
  `ModifiedResiduePolicy.ATOMIC` stays `hetatm_amino_acids:PTR,SEP,TPO`; the Cα
  methods record `ModifiedResiduePolicy.CALPHA`,
  `hetatm_amino_acids:CME,CSD,CSO,LLP,M3L,MSE,NIY,PTR,SEP,TPO,TRF;ca_element:C_or_blank`
  (the suffix tells these results from the earlier build that required `N`
  and `C`, which recorded the names alone; such outputs are not reused).
  `biotsr.api.processing_identity(method)` takes the method. An atomic output
  of an earlier build is therefore still reusable; a Cα one is not.
* **Resume checks the key space.** `--if-exists skip` also requires every
  result entry's `KeySpace` fingerprint to equal the one this version computes
  for the request. Analyses already refuse results of different key spaces
  (ADR-094), so a Cα result written before this release cannot be compared with
  a current one; it has to be computed again.

**What changes for existing results.** Every Cα key is renamed: a triangle of
regular residues keeps its labels and classes, and its key is the same tuple
written at radix 30 instead of 20 (or 13 under grouping). A chain containing a
supported modified residue gains that residue's point and every triangle that
contains it; on 6G9M chain A, CME 161 adds 57,970 triangles (C(341, 2)) and the
result is no longer partial.

**Evidence.** `tests/data/reference/standard_baseline.json` is unchanged. The
C-alpha cases' current values are in `tests/data/reference/calpha_base30.json`,
derived without the changed code: the previous build reproduced the rc1 files,
every key was decoded at its old radix and re-encoded at 30, and for 6G9M the
triangles containing CME 161 were added from an independent calculation. The
changed code reproduces each file byte for byte, and
`tests/integration/test_standard_references.py` checks both the values and the
remap back to the rc1 files. Sampling fingerprints are unchanged; every Cα
key-space fingerprint changes.

**Not changed.** The atom lexicon (`atom_lexicon_standard_v1.csv`) and its
fingerprint, every atomic radix, and every angle and distance bin.

**Rationale.** A key is only useful if it can be read back as one residue
combination in one shape class. Codes past the radix make that impossible
silently; a radix that spans the codes restores it, and a check in the scheme
and in the encoder keeps it so.

---

## ADR-108 — Reading atom records: the specification's fields, and refusal of input that cannot be read reliably

**Status:** Accepted (1.0.0rc6; amends ADR-009, ADR-005 and ADR-082)

**Decision.**

* **Atom names** are the full PDB atom-name field, columns 13–16, for every
  atomic method and every atomic cross side (contact atoms and vertices). The
  Cα test alone reads columns 14–16: an atom is a Cα only when it is named
  exactly `CA` there, so calcium (written from column 13) and `CA1` never are.
* **Elements** are read from columns 77–78 only, never the charge columns
  79–80, and compared without regard to letter case for the hydrogen test and
  for labelling. Nothing infers an element from an atom name.
* A calculation is **refused**, naming the records, when the records it would
  use include: an atom with a blank element field (atomic calculations and the
  contact atoms of every cross; not a Cα point, whose `HETATM` CA may be blank,
  ADR-107); two records of one atom (model, chain, residue number, insertion
  code, atom name), whether repeated or conflicting; two different atoms at
  exactly the same coordinates; a record whose residue number cannot be read
  (hybrid-36 is not decoded); or an atom the label table cannot label (since
  1.0.0rc7, an element the element table does not list, ADR-113). Nothing
  is merged, chosen, repaired or dropped, and no distance tolerance is used.
  Records the calculation would not use are never checked.
* Results record `source.input_rules = "input_rules_v2"`, and the processing
  identity (ADR-105) carries it, so a result of 1.0.0rc5 or earlier (*not
  recorded*) is not reused on resume and not mixed silently in analysis.

**Rationale.** Decided by the TSR Research Group. The narrow slice respelled
two-letter elements; the amino-acid analyzer's columns 78–80 read into the
charge field, so a hydrogen written `H0` or `H1+` counted as a heavy atom and
failed the residue's exact count. A blank element cannot be told from a
hydrogen, and a label then comes from the name (`HG` reads as mercury). Repeated
and coincident records multiply triangle counts silently, and a skipped
unreadable residue reports an incomplete selection as complete. Each of these
gave a well-formed result that was wrong; refusing it is the only answer that
cannot be mistaken for a finding. A partial-output option for unlabelled
atoms was considered and not added.

**Effect.** No standard reference result changes: standard atoms start in
column 14, every reference file has element fields, and none has repeated,
coincident or unreadable records.

---

## ADR-109 — Selections: UNK, cross sides and where protein sides stop

**Status:** Accepted (1.0.0rc6; amends ADR-049, ADR-050 and ADR-052)

**Decision.**

* **UNK.** No calculation is done for a residue named `UNK`. A calculation is
  refused when UNK is part of what it selects: a whole chain that contains one,
  a named UNK residue, residue type or atom, a whole cross side, or a narrowed
  Cα-cross contact selection. It is never dropped silently with the rest
  computed, and never makes contact without a point. A chain the calculation
  does not select may contain UNK. Every other residue-name rule is unchanged
  and method-specific (KNOWN_LIMITATIONS §2).
* **Cross sides.** The two chains of `protein_protein_cross` and
  `protein_protein_calpha_cross` must differ. A molecule side (drug,
  carbohydrate, nucleotide) may share the protein's chain identifier, since the
  sides are told apart by what they select. Every cross refuses two sides that
  share a physical atom, naming the atoms; nothing is subtracted.
* **The narrowed Cα cross.** With `calpha_residues` / `calpha_residues_2`, the
  selections are applied on both sides first; interface residues are those
  whose heavy atoms are within the cutoff of the other side's *selected*
  residues' heavy atoms (ADR-052's rule), each represented by its supported Cα.
  A partner residue touching only unselected residues is not part of it.
* **Where protein sides stop.** Every protein side stops at its chain's `TER`,
  including the protein side of `carbohydrate_protein_cross` and
  `nucleotide_protein_cross`; the molecule side still reads the whole of model 1.
  A `TER` without a chain identifier ends no chain (no inference is made).
* **Residues outside a side's vocabulary are named.** The whole-chain protein
  side of the carbohydrate- and nucleotide-protein cross reads the 23 residue
  names of their historical contact programs; any other residue of the chain is
  reported as `not supported`, and the result is partial (ADR-106).

**Rationale.** Decided by the TSR Research Group. A UNK residue has no label
and no defined atom set; computing around it reports a different selection
than the one requested. A residue crossed against its own chain, or an atom on
both sides, counts each triangle several times (1BRS D/D SER 14: every
triangle six times). The narrowed Cα cross kept partner residues that touch
only unselected residues (1BRS, A 27 against D: 3 triangles now, 7 partner
residues and many more triangles before). The two protein sides that read past
`TER` admitted free residues written after the chain.

**Effect.** Whole-chain Cα and atomic cross results without UNK, overlap or
records after `TER` are unchanged (1BRS A/D: 41 points, 8,151 triangles, as
before). Narrowed Cα cross results change; 1BRS A {27, 83, 87, 102} against D:
798 → 336 triangles.

---

## ADR-110 — SSE endpoints are matched with their insertion codes, in chain order

**Status:** Accepted (1.0.0rc6)

**Decision.** A `HELIX`/`SHEET` record covers the residues of its chain from
its start residue to its end residue in the order the chain's records are
written (up to the chain's `TER`), each endpoint matched by chain, residue
number **and insertion code**. Residues with insertion codes take part in SSE
runs and are annotated when a record covers them. If an endpoint cannot be
resolved — absent from the chain, written in two separate places, with an
unreadable (for example hybrid-36) number, on another chain, or after the
other endpoint — SSE annotation of that chain is refused, naming the record;
there is no numeric fallback and no residue becomes coil by default. Insertion
codes are kept apart from alternate-location indicators and from residue
numbering. `sse.version` is `sse_assignment_v3`.

**Rationale.** Decided by the TSR Research Group. The numeric range of the
historical `SSE_TSR.py` never covered an insertion-coded residue and could
cover one outside the record (1AE8 chain L: helix 14C–14I was read as residue
14 helix and 14C–14I coil). File order, not number, is the residue order:
insertion codes can precede their number (3PTB 184A before 184).

**Effect.** Keys, points and both fingerprints are unchanged; only SSE
categories of residues in or next to insertion-coded ranges change.

---

## ADR-111 — Terminal residues follow each method's scope

**Status:** Accepted (1.0.0rc6; states how ADR-004 applies)

**Decision.**

* **Cα methods** keep the first and last residues' Cα points; they have no
  completeness or backbone requirement.
* **Residue-specific and residue-type atomic methods** (`amino_acid_atomic` and
  its per-residue sweep, `aa_ind`, `aa_ind_mix`, `aa_ind_specific`, and the
  single-molecule `nucleotide` and `carbohydrate` methods) keep their exact
  atom-count rule (ADR-004): a terminal variant with a different atom count —
  a C-terminal residue carrying `OXT`, a 5′ nucleotide without its phosphate —
  is excluded and named.
* **Chain-only atomic methods** include the first and last residues with the
  supported heavy atoms they have, `OXT` included; hydrogen is excluded and no
  atom is invented. `carbohydrate_chain` has no completeness rule and does so.
  `nucleotide_chain` applied the per-nucleotide atom count in 1.0.0rc6; since
  1.0.0rc7 it has none and names every nucleotide whose atoms differ from the
  standard set (ADR-114).
* **Cross methods** apply no completeness rule to either side; a terminal
  residue's atoms take part like any other.
* There is no first/last-residue detector and no carbohydrate-endpoint rule.

**Rationale.** Decided by the TSR Research Group: a chain sample describes the
whole chain, while a residue or residue-type sample compares chemically
identical instances, which a terminal variant is not.

---

## ADR-112 — Every calculation names its chains; a batch row without one is skipped

**Status:** Accepted (1.0.0rc6; amends ADR-070's batch handling)

**Decision.** Every calculation needs its chain, and a cross calculation its
second chain, from the command, the row, or an explicit `--chain` default; the
package never assumes chain A or the first chain. A single run without its
chain is refused with an actionable message. In a batch, a request whose
required chain is still blank after the explicit defaults is not run: its
status is `chain_not_specified`, it is listed in `batch_summary.csv` with the
reason, the other requests run, and the batch exits 4 (ADR-106).

**Rationale.** A batch of many rows should not be lost because one row is
incomplete, and a missing chain is not a failure of the calculation — none was
attempted — but it must never be filled in by guessing.

---

## ADR-113 — Atomic labels come from the element field

**Status:** Accepted (1.0.0rc7; supersedes ADR-030 and ADR-031, amends ADR-029,
ADR-032 and ADR-108)

**Decision.**

* Every atomic method — amino-acid atomic, residue-type, drug, carbohydrate,
  nucleotide and the atomic cross methods — labels a heavy atom by its element
  field: PDB columns 77–78, never the charge columns 79–80, read without regard
  to letter case. The atom name is not consulted for the label; it stays the
  atom's identity in selections, records and outputs.
* The labels come from one packaged table, `labels/resources/element_labels_v1.csv`,
  with two columns, `element,label`: 110 element symbols, 109 distinct labels,
  4..112. Every element keeps the code the historical atom lexicon gave it
  (carbon 106, nitrogen 108, oxygen 109, phosphorus 110, sulfur 111, chlorine
  93, calcium 99, holmium 37, …). Hydrogen and deuterium share 107, the one
  documented exception to one label per element; both are excluded before any
  atom is labelled, as before.
* The radix stays 112 (ADR-032). It is not reduced to the number of labels:
  the labels fit one window of 112 consecutive integers, which keeps the
  encoding collision-free (ADR-107).
* A blank element field, or one naming an element the table does not list
  (`X`, the superheavy elements the lexicon never had), cannot be labelled, and
  a calculation that would use such an atom is refused, naming it (ADR-108).
  Nothing is inferred from the name, repaired or checked against bond lengths.
* The table is validated when it is loaded: a wrong header, a row with the
  wrong number of fields, a symbol not in chemical case, a label that is not a
  positive integer, a symbol listed twice and two elements sharing a label
  (other than hydrogen and deuterium) are each refused. There is no second
  copy of the table.
* The scheme is `element_labels_112`. Its content fingerprint covers the table
  and the resolver (`label_source = element`, `resolver = element_field_v1`),
  so every atomic key space has a new fingerprint. The Cα residue labels and
  their key spaces are unchanged.

**Rationale.** Decided by the TSR Research Group. The historical lexicon's codes
are element classes — each of its 109 codes belongs to exactly one element,
and every standard amino-acid and nucleotide atom name is listed on its
element's code — but it reached the code through the atom name, so the label
depended on naming conventions rather than chemistry. In public RCSB entries
an ATP's `PA` and `PB` were protactinium (55) and lead (56) (1OL6), `N1A` of
NAD and `CB1` of a ligand were hydrogen (107) (1EVZ, 1UFJ), a ligand chlorine
`CL2` was carbon (4X6J), the calcium of photosystem II's Mn₄CaO₅ cluster
(`CA1`) was carbon (2AXT), and a chlorine named `CL35` could not be labelled,
which refused its calculation (3FSF); these records are in
`tests/unit/test_element_labels.py`. Atom names are not a reliable source of
the element: `CL1` is a chlorine in some ligands and a carbon in others (1BXQ). The element field is the
specification's statement of the element and is already required (ADR-108).

**Effect.** Every atomic key space has a new fingerprint, so an atomic result
computed before 1.0.0rc7 is not reused by `--if-exists skip` and is never
compared with a current one; no option relaxes that (ADR-094). The numbers
change only where an atom's name and element disagreed or it could not be
labelled: all 25 atomic standard reference results keep their keys, counts
and sampling fingerprints (`nucleotide_chain` apart, ADR-114); their triplet
files differ only in how a label name is spelled (`Cl` for `CL`). The record
of the reference expectations is `tests/data/reference/atomic_element_labels.json`.

---

## ADR-114 — `nucleotide_chain` takes every nucleotide with the atoms it has

**Status:** Accepted (1.0.0rc7; amends ADR-111)

**Decision.**

* `nucleotide_chain` has no atom-count rule. Each of the eight standard
  nucleotides of the chain (DA, DC, DG, DT, A, C, G, U) takes part with every
  admitted heavy atom it has, at a chain end or inside the chain.
* Every nucleotide whose atom names differ from its type's standard set — the
  PDB's heavy-atom names without the 5′-terminal `OP3` — is named in one
  warning, with the atoms it lacks and the ones it has in addition. The
  warning describes the structure file: every atom present was used, so the
  result is not partial (ADR-106) and the batch status stays `ok`.
* No coordinate is invented, no atom reconstructed, no chain end detected.
  A nucleotide type outside the eight is still not read, whatever vocabulary
  is given. The model, alternate-location, `TER`, `UNK` and input-integrity
  rules (ADR-095, ADR-108, ADR-109) are unchanged.
* The residue-specific methods `nucleotide` and `nucleotide_atom_selection`
  keep the exact count (ADR-004, ADR-111); the cross methods still apply no
  completeness rule.
* The method's completeness rule is `available_nucleotide_atoms`. It is in the
  sampling profile and, since 1.0.0rc7, in the processing identity of every
  method that has a completeness rule (ADR-105), so a `nucleotide_chain`
  output computed under the count rule is not reused.

**Rationale.** Decided by the TSR Research Group, applying the chain-only rule
for terminal residues (ADR-111): a chain calculation includes all available
supported atoms, including those of incomplete nucleotides within the chain,
and reports what is missing or extra; the residue-specific calculations keep
their completeness checks. Without a terminal-residue detector the count rule
cannot be kept for internal nucleotides and dropped for terminal ones, and a
detector was not adopted.

**Effect.** 1BNA chain A: 227 points before (DC 1, the 5′ nucleotide, has 16
of its 19 atoms and was dropped), 243 now — the 16 added points are DC 1's —
with the warning `DC_1 lacks P, OP1, OP2`. The triangles without a DC 1 atom
are the previous result's, key for key.

