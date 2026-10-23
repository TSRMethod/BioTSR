# Secondary-structure-annotated TSR (SSE-TSR)

This page describes the secondary-structure extension of the Triangular Spatial
Relationship descriptor as implemented in BioTSR 1.0. It states what is
computed and how. It reports no results and makes no biological or performance
claims.

Companion pages: [CROSS_TSR.md](CROSS_TSR.md) (cross-entity TSR) and
[ANALYSIS.md](ANALYSIS.md) (feature representation, common keys, similarity,
clustering). For the commands that run SSE-annotated TSR, see
[../REPRODUCIBILITY_NOTE.md](../REPRODUCIBILITY_NOTE.md); the decision records
are ADR-063 to ADR-067 in
[../architecture/DECISIONS.md](../architecture/DECISIONS.md).

---

## 1. Scope

SSE annotation is defined for the **protein Cα representation**: one point per
residue, placed at that residue's Cα atom, labelled by residue identity. It
applies to whole-chain Cα TSR and to residue-selected Cα TSR, with amino-acid
grouping, mirror-image signing and size filtering in any combination.

It is **not** defined for intra-residue atomic TSR, ligand, carbohydrate,
nucleotide or cross-entity TSR. An atom within a residue has no secondary
structure of its own, and a ligand, a glycan and a nucleotide have none at all;
extending the scheme to those representations would require a separately defined
categorisation rather than a reuse of the one described here.

---

## 2. Secondary-structure assignment

Secondary structure is taken from the deposited structure's own `HELIX` and
`SHEET` records. No secondary-structure assignment is computed from coordinates;
no external assignment program (DSSP, STRIDE) is invoked. The annotation
therefore reflects the depositor's assertion about the model.

Each record declares a residue range on one chain, from its start residue to
its end residue. Both endpoints are found in the chain by chain identifier,
residue number **and insertion code** (`initSeqNum`/`initICode`,
`endSeqNum`/`endICode`), and the record covers every residue of the chain from
the start residue to the end residue in the order the chain's records are
written, up to the chain's `TER` (ADR-110). Residues with insertion codes are
therefore covered when they lie inside a range (`14`–`15` covers `14A` and
`14B`), and a range from `14C` to `14I` covers those residues and not residue
14. File order, not number, decides: a file may place `184A` before `184`.
Records on other chains are ignored.

When a record of the analysed chain cannot be read exactly — an endpoint that
is not in the chain, that appears in two separate places, whose number cannot
be read (a hybrid-36 number), that lies on another chain, or a range that ends
before it starts — the SSE annotation of that chain is **refused**, naming the
record. There is no numeric fallback, and no residue silently becomes coil.
The same calculation without SSE is unaffected. (Until 1.0.0rc5 a record was
read as the numeric interval `[initSeqNum, endSeqNum]`, ignoring insertion
codes, so an insertion-coded residue was never covered.)

Three precedence rules complete the assignment:

1. a residue covered by both a `HELIX` and a `SHEET` record is assigned **helix**;
2. a residue covered by two records of the same kind is assigned the **first** in
   file order;
3. a residue covered by no record is assigned **coil**.

### Element identity

Two residues assigned the same state are further distinguished by which
*individual element* they belong to, since the categorisation separates "two
vertices in the same helix" from "two vertices in different helices", and
likewise for β-strands. The categories distinguish **individual β-strands**; a
β-sheet comprises one or more strands, and no category refers to a whole sheet.

An individual helix is identified by its `HELIX` record's `serNum`, which is
unique within a structure.

Individual β-strands were distinguished using the combination of the PDB `SHEET`
record's sheet identifier (`sheetID`) and the strand serial number (`serNum`),
preventing strand-number collisions between distinct β-sheets. This is necessary
because `serNum` is the strand's index *within its sheet* and restarts at 1 in
each sheet, so it does not identify a strand uniquely within a structure; the
sheet identifier is used solely to disambiguate strands, and strands of the same
sheet remain distinct elements. Two residues are therefore in the same β-strand
only when both fields agree.

Categories 7–12, 15 and 16 below are statements about **individual β-strands**.
A triangle whose three vertices lie in three different strands of a single
β-sheet is category 9 (`9_3b3`).

---

## 3. State representation

Three states are used: **helix**, **β-strand** and **coil**. Coil is the
absence of a `HELIX`/`SHEET` record, and there is no separate "unknown" state:
the assignment is total, so every residue in the calculation carries exactly one
state and every triangle is classifiable.

One limitation follows directly and should be controlled for. A structure that
carries no `HELIX` or `SHEET` records at all — for example a predicted model or a
stripped coordinate file — yields all-coil assignments and therefore places every
triangle in the all-coil category. That output is not distinguishable from a
genuine all-coil finding by the categories alone, so the presence of
secondary-structure records should be verified before interpreting an
SSE-annotated result. The implementation reports whether any record covered the
chain and emits an explicit warning when none did.

---

## 4. The 18 triplet categories

Each triangle carries three vertices, hence three states and, where applicable,
three element identities. The triangle is assigned to one of 18 categories
determined entirely by the signature

    (n_helix, n_strand, n_coil, n_distinct_helices, n_distinct_strands)

where the first three count vertices by state (summing to 3) and the last two
count the distinct elements the helix and β-strand vertices occupy. This signature
space has exactly 18 inhabitants, and the mapping to categories is a bijection:
the classification is total and has no default class.

| # | name | signature | definition |
|---|---|---|---|
| 1 | `1_3a1` | (3,0,0,1,0) | 3 vertices in the same helix |
| 2 | `2_3a2` | (3,0,0,2,0) | 3 vertices in 2 distinct helices |
| 3 | `3_3a3` | (3,0,0,3,0) | 3 vertices in 3 distinct helices |
| 4 | `4_2a11c` | (2,0,1,1,0) | 2 vertices in the same helix, 1 in coil |
| 5 | `5_2a21c` | (2,0,1,2,0) | 2 vertices in 2 distinct helices, 1 in coil |
| 6 | `6_1a2c` | (1,0,2,1,0) | 1 vertex in a helix, 2 in coil |
| 7 | `7_3b1` | (0,3,0,0,1) | 3 vertices in the same β-strand |
| 8 | `8_3b2` | (0,3,0,0,2) | 3 vertices in 2 distinct β-strands |
| 9 | `9_3b3` | (0,3,0,0,3) | 3 vertices in 3 distinct β-strands |
| 10 | `10_2b11c` | (0,2,1,0,1) | 2 vertices in the same β-strand, 1 in coil |
| 11 | `11_2b21c` | (0,2,1,0,2) | 2 vertices in 2 distinct β-strands, 1 in coil |
| 12 | `12_1b2c` | (0,1,2,0,1) | 1 vertex in a β-strand, 2 in coil |
| 13 | `13_2a21b` | (2,1,0,2,1) | 2 vertices in 2 distinct helices, 1 in a β-strand |
| 14 | `14_2a11b` | (2,1,0,1,1) | 2 vertices in the same helix, 1 in a β-strand |
| 15 | `15_2b21a` | (1,2,0,1,2) | 2 vertices in 2 distinct β-strands, 1 in a helix |
| 16 | `16_2b11a` | (1,2,0,1,1) | 2 vertices in the same β-strand, 1 in a helix |
| 17 | `17_3c` | (0,0,3,0,0) | all 3 vertices in coil |
| 18 | `18_1a1b1c` | (1,1,1,1,1) | 1 vertex in a helix, 1 in a β-strand, 1 in coil |

In the category names, `a` denotes helix, `b` an individual β-strand and `c`
coil. The numeric prefixes 1–18 are part of the names and fix the column order used when
the distribution is tabulated.

**The classification is invariant to vertex order.** All six permutations of a
triangle's three vertices give the same category, so the canonical vertex
ordering used to construct the key does not enter the categorisation.

---

## 5. Relationship between the annotation and the TSR key

The annotation is a **label on a triangle, not a component of its key**. For a
given structure, chain, residue labelling and geometry, the integer TSR key of a
triangle is identical whether or not secondary structure is annotated. The
annotation introduces no new integer key space, changes no key, changes no key's
sign, and changes neither the set of triangles enumerated nor the set retained.

What the annotation provides is a **stratification of the key frequencies**. An
ordinary TSR result reports, for each key, the number of triangles that produced
it. An SSE-annotated result reports, for each key, that same number broken down
across the 18 categories. Summing a key's categories recovers its ordinary
frequency exactly, and summing over all keys recovers the total number of
retained triangles:

    Σ_categories f(key, category) = f(key)          for every key
    Σ_keys Σ_categories f(key, category) = N        the retained triangle count

Both identities hold unconditionally, because the classification is total (§3):
no triangle is unclassified and none is discarded by the annotation.

Consequently an SSE-annotated result and the corresponding unannotated result are
directly comparable, and any downstream analysis defined over ordinary key
frequencies — feature matrices, common keys, similarity, clustering — consumes
either interchangeably, by summing the categories. No separate SSE feature
representation is defined here.

---

## 6. Interaction with the Cα representation's other options

**Residue selection.** When the calculation is restricted to an explicitly listed
set of residue positions, those residues retain their own secondary-structure
states. The representation remains the Cα one — one point per residue, labelled
by residue identity over the residue labels (the 20 regular amino acids and the
eleven supported modified residues) — and is not reinterpreted as an atomic
representation. A supported modified residue inside a `HELIX` or `SHEET` range
takes that range's state like any other residue.

**Amino-acid grouping.** Grouping replaces the residue labels with the grouped
labels (13 groups of the regular amino acids, MSE with MET, each other modified
residue on its own) and encodes keys at the same radix 30 in a different key
space. It changes what a key means. It does **not** change any residue's secondary structure: a residue's
state is read from its position in the chain, never from its label. A grouped and
an ungrouped run over the same chain therefore produce different keys and the
identical category distribution.

**Mirror-image signing.** When mirror handling is enabled, a triangle's key is
negated if and only if its *pre-fold* angle θ_raw exceeds 90°; the folded angle
θ ∈ [0°, 90°] determines the angle bin, as without mirroring. The secondary-
structure category takes no part in this decision, and the sign takes no part in
the category. Positive and negative keys are each stratified across the 18
categories independently.

**Size filtering.** Size filtering is optional and off by default. When enabled, a
triangle is retained if and only if its maximum inter-vertex distance is
**strictly less than** the cutoff. The decision is made from geometry alone,
before the key is formed and before any annotation is consulted, so annotation
cannot affect which triangles are retained. The order of operations is:
candidate → geometry → size-filter decision → retained triangle → key →
category.

---

## 7. Conventions, and how they differ from the original program

BioTSR computes SSE-TSR under one set of conventions: π taken as `math.pi`; the
35-class Cα distance scheme; the residue labels (30 codes), or the 23 grouped
codes, both at radix 30; IEEE double-precision arithmetic throughout; β-strands named by
`(sheetID, serNum)`; and deterministic aggregation in which every key owns its
own category counter. BioTSR does not reproduce the original SSE program's
output files. The original program differs in ways that are individually
attributable, and published results derived from it should be read with them in
mind.

Four differences are shared with the original non-SSE Cα program:

* π was taken as the literal `3.14`, which shifts computed angles and moves some
  triangles across angle-bin boundaries;
* the distance scheme's highest class was not reachable, giving 34 effective
  classes rather than 35;
* under grouping, the 13 grouped labels were encoded in the unchanged 20-symbol
  space;
* it read `ATOM` records only and labelled the 20 regular amino acids at radix
  20, so a modified residue was never a point; BioTSR reads eleven of them and
  encodes at radix 30 (ADR-107).

Three are specific to the original SSE program:

* **β-strand identity.** It identified the β-strand a residue belongs to by the
  `SHEET` record's `serNum` alone. Because `serNum` restarts at 1 in each sheet,
  strands of different sheets that share a serial number were treated as one
  element. On structures whose analysed chain carries a single β-sheet the two
  conventions agree exactly; on a chain with several sheets they differ wherever
  a triangle spans strands with a shared serial number, always in the direction
  of resolving strands the original merged. Over eleven reference runs used to
  characterise this difference, it affects 27,228 of 2,168,260 triangles
  (1.26%), all within the one multi-sheet chain in that set (1TUP chain A:
  27,228 of 1,235,780, 2.20%). No TSR key changes; only the category a triangle
  is counted under.
* **Point set.** It admitted a residue only when its sequence-number field
  parsed as a plain integer. Because the field it read spans the insertion code
  as well as the number, residues carrying an insertion code were excluded from
  the calculation entirely. BioTSR retains them, as the original non-SSE program
  did, so the two agree exactly on structures without insertion codes.
* **Aggregate category file.** Its per-key category file is affected by an
  aggregation defect: a single category vector was shared and re-bound across
  keys, so most keys report another key's counts. On 1UBQ chain A the reported
  vectors total 111,017 against a true 70,300 triangles, and 18,119 of 45,353
  rows disagree with their own key's frequency. **The per-triangle category
  assignments are unaffected** and correct. Analyses derived from the original
  per-key category file should be treated with this in mind; analyses derived
  from its per-triangle records are unaffected.

None of the three arises in BioTSR: β-strands are named by sheet and serial,
insertion-code residues are retained (and, since 1.0.0rc6, annotated when a
record's range includes them), and each key owns a private counter. The
decision records are ADR-063 to ADR-067 in
[../architecture/DECISIONS.md](../architecture/DECISIONS.md).
