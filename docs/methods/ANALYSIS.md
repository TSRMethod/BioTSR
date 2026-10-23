# TSR analysis

This page describes how TSR key sets are compared, summarised and clustered in
BioTSR 1.0. The measures and the clustering procedure are those of the TSR
group's earlier analysis programs; how far the clustering agrees with SciPy is
stated precisely in §5.7. A worked, executed example of the whole workflow —
saved results to feature matrix, similarity, clustering and figure — is
[tutorial 08](../tutorials/08_similarity_clustering.md).
It complements [CROSS_TSR.md](CROSS_TSR.md) and [SSE_TSR.md](SSE_TSR.md), which
describe how keys are produced. For the commands that run these analyses, see
[../REPRODUCIBILITY_NOTE.md](../REPRODUCIBILITY_NOTE.md); the decision records
are ADR-054 to ADR-062, ADR-094 (the comparison rule) and ADR-099 (the
hierarchy without a cluster count) in
[../architecture/DECISIONS.md](../architecture/DECISIONS.md).

This page states the scientific procedure; where a choice had alternatives, it
says which was taken and what distinguishes it. Implementation details are
described in [../ARCHITECTURE.md](../ARCHITECTURE.md).

It makes **no claim about biological performance or scientific superiority.**
One quantitative result is reported below — the recovery of nucleotide identity
on a single B-DNA dodecamer — and it is reported as what was computed on that
structure, not as evidence that the method outperforms anything.

---

## 1. What the analysis layer operates on

A TSR calculation produces, for one sample, a set of integer **keys** with an
integer **frequency** for each. A sample is whatever the calculation was asked
about: a structure, a chain, a residue, a ligand, a glycan, a nucleotide, or an
interface between two of these.

All analysis below is defined on that object and on nothing else. It does not
read coordinates, does not re-derive geometry, and does not know which TSR method
produced the keys. Consequently the same procedures apply unchanged to protein
Cα TSR, amino-acid grouped TSR, mirror-signed TSR, size-filtered TSR, amino-acid
atomic TSR, drug, carbohydrate and nucleotide TSR, and every cross-TSR family.

### 1.1 When two key sets may be compared

A key is an integer that encodes a triangle's labels, its angle class and its
distance class under a specific set of conventions. Two calculations that used
different conventions produce integers that are *not* comparable even when they
are numerically equal: under amino-acid grouping the same integer denotes a
different triple of residues, and under mirror signing a negative key denotes a
chirality that an unsigned calculation does not distinguish at all.

Every analysis on this page — the feature matrix, common keys, occurrence
tables, similarity and clustering, on computed or loaded results — applies one
rule, in this order:

1. **The key definitions must match. Always; no option relaxes it.** The inputs
   must share a **key space**: the angle convention, the square-root and median
   spellings, the angle and distance binning schemes with their content, the
   label scheme and its radix, the canonical vertex ordering, the key formula
   and the mirror convention, together with the floating-point precision the
   coordinates were held at. This is decided from the key-space definition,
   never from method names. Different methods that share a key space are
   comparable on this count — `aa_ind`, `drug`, `carbohydrate` and every other
   atomic method, cross methods included, label atoms by element at radix 112
   in 29 × 18 classes;
   `protein`, `aa_selection` and `protein_protein_calpha_cross` share the Cα key
   space — while the same method with `mirror` or `aa_grouping` has a different
   key space and is refused.
2. **The processing settings must match, unless the difference is
   acknowledged.** A difference in how triangles were drawn from whatever was
   selected leaves every key's meaning unchanged, but the frequencies — and
   which keys are present at all — answer a different question. The settings
   compared are the candidate strategy, the completeness and composition rules,
   the size filter, the degenerate-triangle and sequence-label policies, the
   contact cutoff, contact filter, contact rule and projection, the pool
   identity and order, and the alternate-location and modified-residue rules
   the result recorded (the modified-residue rule says which modified amino
   acids written as `HETATM` were read; a result written before it was
   recorded is reported as *not recorded* and is never taken to follow the
   current rule). A difference is refused unless the study explicitly allows it
   (`allow_mixed_processing=True` in Python, `--allow-mixed-processing` on the
   command line); every result produced under the opt-in carries a warning
   naming each differing setting, and records every sample's settings. The
   opt-in never relaxes step 1. Methods that share a key space can differ here:
   the amino-acid and carbohydrate methods apply their own atom-count
   completeness rules and `drug` applies none, so `aa_ind` against `drug` needs
   the opt-in.
3. **A partial sample must be acknowledged, separately.** A sample whose
   selection did not get everything it asked for — a requested residue, atom,
   molecule or residue type absent from the structure, a modified residue the
   method cannot read, a residue or molecule copy named in the request that
   the completeness rule rejected, or requested atoms a selected copy lacks —
   is not the selection that was requested. A set containing one is refused,
   naming each partial sample and what it lacks, unless the study explicitly
   allows it (`allow_partial=True` in Python, `--allow-partial` on the command
   line). Nothing is then dropped: the report records each partial sample's
   omissions and carries a warning naming them. This acknowledgement is
   independent of the one in step 2 — neither implies the other — and neither
   relaxes step 1. A residue the completeness rule rejects among all residues
   of a requested type, or a copy among all copies of a molecule, does not make
   a sample partial.

What a calculation **selected** never blocks a comparison: the structure, the
chains, the amino-acid types (a LYS-only set against a TYR-only set needs no
opt-in), the ligand, carbohydrate or nucleotide, the residue numbers, the atom
subset, one copy of a molecule or all copies, and the method name. These are
what a dataset is made of, and a drug dataset compares BEN in one entry with ESW
in another by design. A complete LYS-only sample and a complete TYR-only sample
compare without any acknowledgement; only a selection that fell short of its
own request needs one (step 3). Two samples need not share any key: an **empty
intersection is a valid result**, not an error.

Every report records, per sample: its id and method, the structure, the
key-space fingerprint, the processing fingerprint and settings, the full
sampling fingerprint and the sampling settings that are set, the selection it
actually used (per side: chain, scope, arguments and every molecule or residue
instance), its target and chains, its number of keys and total frequency, the
warnings the sample's own result carried, and whether it is partial and what
it lacks. The summary of the set records
the key-space fingerprint, whether processing was uniform and whether mixing was
allowed, which samples share each processing fingerprint, the differing
settings, the methods and selection scopes in the set, whether partial samples
were allowed, and the omissions of each partial sample.

Two of the recorded identifiers are not purely semantic. The key-space
fingerprint hashes each scheme's name beside its content fingerprint, which is
conservative: a renamed but otherwise identical scheme would be refused, never
wrongly accepted (no such case exists among the twenty methods). The full
sampling fingerprint hashes each side's chain and free-text description, so it
can differ between samples whose differences never block a comparison; it is
recorded for provenance and is not used to decide comparability — the
processing fingerprint is.

---

## 2. Feature representation

Given *n* samples, let *K* be the union of all keys occurring in any of them. The
**feature matrix** is the *n × |K|* matrix whose entry *(i, k)* is the frequency
of key *k* in sample *i*, or zero if that key does not occur in that sample.

Frequencies enter as raw integer counts. They are not normalised, not scaled, not
binarised and not log-transformed.

Columns are ordered by **increasing integer key value**, so that negative
mirror keys precede positive ones in numerical order. The ordering is a
presentational choice with no effect on any quantity defined below: every measure
in §4 is a sum over all columns evaluated with the same ordering on both
operands, and is therefore invariant under any permutation of the columns. The
ordering is fixed so that a matrix written to disk is reproducible.

The matrix is stored sparsely. A key absent from a sample and a key present with
frequency zero are the same fact, and no procedure below distinguishes them.

---

## 3. Common-key analysis

Four distinct quantities are defined. They answer different questions and are
never interchangeable.

**3.1 Global common keys.** The set of keys present in *every* sample of the set:
the intersection of the *n* key sets. Membership depends on presence only;
frequency plays no part. The empty intersection is a valid result.

For each sample, three quantities are reported alongside: the number of distinct
keys it contains, the size of the intersection, and the sum of that sample's
frequencies over the intersection.

**3.2 Group-wise common keys.** §3.1 applied independently within each group of a
partition of the samples — by residue type, by chain, or by any other grouping the
study defines.

**3.3 Key occurrence percentage.** For each key, the number of samples in which
it occurs at all, expressed as a percentage of the sample count. This is a
property of every key, not a subset of keys: a key occurring in two of three
samples has an occurrence of 66.7% and is not a common key in the sense of §3.1.

**3.4 One-versus-all occurrence percentage.** §3.3 computed for every group of a
partition and presented as one table, one row per key and one column per group.
A key absent from a group has occurrence 0% in that group.

---

## 4. Similarity and distance

All measures are defined below as **similarities** on the interval [0, 1], with 1
denoting identity. The corresponding distance is `1 − J` in every case, and the
conversion is always stated explicitly; a similarity matrix and a distance matrix
are never interchanged silently.

Let *a* and *b* be the frequency vectors of two samples over the shared key set,
with *a_k* the frequency of key *k*.

**4.1 Generalized Jaccard** (the primary measure):

```
J(a, b) = Σ_k min(a_k, b_k) / Σ_k max(a_k, b_k)
```

This is frequency-weighted: a key occurring five times in one sample and once in
the other contributes 1 to the numerator and 5 to the denominator. Sums are
computed exactly over integers; only the final division is in floating point.

**4.2 Binary Jaccard**, on presence sets *A* = {k : a_k > 0}:

```
J(a, b) = |A ∩ B| / |A ∪ B|
```

Frequencies are discarded.

**4.3 Wu**:

```
J(a, b) = Σ_k min(a_k, b_k) / min( Σ_k max(a_k, b_k), max(Σ a, Σ b) )
```

For non-negative counts `Σ max ≥ max(Σa, Σb)` always holds, so this reduces to
`Σ min / max(Σa, Σb)`.

**4.4 Sarika**:

```
J(a, b) = Σ_{k : a_k > 0 and b_k > 0} (a_k + b_k) / (Σ a + Σ b)
```

the fraction of total counted triangles that lie on keys the two samples share.

**4.5 Cosine**:

```
J(a, b) = (a · b) / (‖a‖ ‖b‖)
```

**4.6 Conventions common to all five.** Every measure is symmetric, gives 1 for a
sample against itself, and lies in [0, 1]. Frequencies must be finite and
non-negative; a negative or non-finite count is refused, and an explicit zero is
the same as an absent key. **Empty samples:** two samples with no keys have
similarity 1 (distance 0) — they are treated as identical — and an empty sample
against a non-empty one has similarity 0 (distance 1); a matrix over a set that
includes an empty sample carries a warning saying so. The diagonal of a
similarity matrix is 1 by definition and of a distance matrix 0 by definition;
neither is computed.

**4.7 Fractions and percent.** Every value is a fraction in [0, 1]. A matrix can
be expressed in percent explicitly (`scaled(100)`); it then records its scale,
its diagonal and conversions follow it (`D = 100 − J`), and every export and
figure label says "percent". A value outside a matrix's scale is refused, so a
percent value cannot be mistaken for a fraction.

**4.8 Which measure is used.** The generalized Jaccard (§4.1) is the primary
measure and the one used for the results reported in §8. §§4.2–4.5 are
implemented because they are defined in the historical source; no result in §8
uses them.

---

## 5. Clustering

**5.1 Input.** A generalized-Jaccard **distance** matrix, `D = 1 − J`, over the
samples. The conversion from similarity to distance is explicit.

**5.2 Algorithm.** Agglomerative hierarchical clustering with **average linkage**
(UPGMA) on the precomputed distance matrix. At each step the two clusters with
the smallest average inter-cluster distance are merged, and the distance from the
merged cluster *xy* to any remaining cluster *i* is

```
d(xy, i) = ( n_x · d(x, i) + n_y · d(y, i) ) / ( n_x + n_y )
```

where *n_x* and *n_y* are the cluster sizes. No other linkage method and no other
metric is used. There is no feature scaling, centering or transformation: the
matrix is clustered as given, and merge heights are in its units (fractions,
unless the caller passed a percent matrix).

**5.3 Determinism.** The procedure is fully deterministic. It contains no random
initialisation, no sampling and no seed: given the same distance matrix it
produces the same dendrogram, the same leaf order and, for a given cutting
rule, the same partition on any platform.

**The default result is the tree.** The hierarchy — the merges with their
heights, and the dendrogram leaf order (§5.5) — is the complete result of the
clustering and needs **no number of clusters**. An uncut tree describes nested
groupings at every height; it does not determine one final number of groups,
and BioTSR reports none. The merges can be exported (`step, left, right,
height, size`) together with the leaf order.

**5.4 Flat clusters, optional.** A flat partition is a separate step, taken
only when the analyst states an explicit cutting rule:

* **`maxclust` *t*** — cut the dendrogram at the lowest height that leaves at
  most *t* clusters. With distinct merge heights exactly *t* clusters form (for
  *t* up to the number of samples). Because the cut is by height, merges that
  tie at the cutting height are taken together, so with ties fewer than *t* may
  form; the result then says so.
* **height *h*** — cut every merge above *h*, in the matrix's units, so that no
  cluster contains two samples joined above *h*.

One rule, never both, and never a default.

**Choosing the rule.** The rule is the analyst's choice, and BioTSR never
infers it — in particular, it never takes a cluster count from the known class
labels. A study that sets *t* to the number of known classes is asking how well
the dendrogram, cut at that granularity, agrees with those classes (§5.6); it is
not an independent discovery of how many classes there are, and the agreement
is partly a consequence of the choice. To explore the number of groups, examine
the merge heights or compare several cuts.

**5.5 Leaf order.** The dendrogram leaf order is obtained by depth-first traversal
of the merge tree, visiting each node's first-recorded child before its second,
without sorting by cluster size or by merge height.

**5.6 Agreement with known classes.** Where a flat partition was made (§5.4) and
every sample carries a known class label, the agreement between the partition
and those labels is quantified by the **adjusted Rand index** (ARI), which is 1
for a perfect correspondence up to relabelling, approximately 0 for a partition
no better than chance, and may be negative. Without a partition there is no
ARI: the tree alone is not compared with the classes by any single number.
Class labels come from the analyst's own table — a manifest's `label` column,
or a file matched by sample id — and are never derived from file or sample
names.

**5.7 Agreement with SciPy.** The clustering is implemented on the Python
standard library, following SciPy's algorithms, and was compared with two SciPy
releases over 400 random matrices (a third of them full of ties), ADR-060:

* with **SciPy 1.11.1**, the linkage matrices and leaf orders are bit-identical
  in all 400; the flat partitions agree except where SciPy 1.11.1 does not find
  the `maxclust` cut (it returns fewer clusters than the definition in §5.4
  gives, or two clusters for two samples and *t* = 1);
* with **SciPy 1.18.1**, the linkage differs in 177 matrices — tie-breaking,
  apart from two last-bit differences in a merge height — and, **given the same
  tree**, the flat partitions and leaf orders agree in every case.

No agreement is claimed for SciPy releases that were not checked. With tied
distances, different SciPy releases build different but equally valid trees;
BioTSR's is deterministic.

**5.8 Scale.** Pairwise similarity and average linkage each take time
quadratic in the number of samples (*n(n − 1)/2* pair values, each pair costing
time linear in the two samples' distinct keys) and run in one process: `--jobs`
parallelises the *structure computations* of a batch, not this analysis. Every
sample's key counts are held in memory. Measured with BioTSR 1.0.0rc3 on one
Apple M2 Pro core (synthetic samples of 20,000 distinct keys each); the
similarity and linkage computations are unchanged through 1.0.0rc5:

| samples | generalized Jaccard matrix | peak memory |
|---|---|---|
| 100 | 13 s | 0.25 GB |
| 200 | 55 s | 0.55 GB |
| 400 | 204 s | 0.9 GB |

— about 1.3 ms per pair at that size, so 1,000 such samples (about 500,000
pairs) would take of the order of ten minutes, and the time grows with the
square of the number of samples and with the keys per sample. Average linkage
on a precomputed matrix took 0.4 s for 500 samples, 2.2 s for 1,000 and 8.5 s
for 2,000. Nothing larger was measured: for thousands of samples, time a subset
first. A clustered heat map of more than about a hundred rows is hard to read
even when it can be drawn; use the exported matrices and assignments.

---

## 6. Visualization

The clustered heat map ("clustermap") displays the distance matrix as given
(fractions in [0, 1], or percent if the caller scaled it, as its colour-bar
label states), with rows and columns ordered by the same dendrogram — one
linkage serves both axes — and the dendrogram drawn on both. By default that is
the whole figure, and it needs no number of clusters.

Two annotation stripes are optional, and are kept distinct in the figure and
its legend:

* **observed class** — the class label the analyst supplied for each sample,
  drawn only when labels were given (an unlabelled sample is left blank);
* **cluster** — the flat partition, drawn only when the tree was cut by an
  explicit rule (§5.4), and labelled with that rule (`maxclust 4`,
  `height 0.5`).

The figure is a display of quantities computed beforehand. The linkage, the leaf
order and, when a cut was asked for, the flat partition and the adjusted Rand
index are all determined without reference to any rendering, and every one of
them is reported independently of whether a figure is drawn.

---

## 7. Numerical precision and reproducibility

Frequencies are integers and the sums in §4 are exact; only the final divisions
are in floating point, at double precision. The clustering in §5 is likewise
computed at double precision.

BioTSR clusters the full-precision distance matrix it is given; it applies no
rounding. One point deserves care when comparing with published results. The
historical pipeline wrote its distance matrix to a text file with **three
decimal places**, after storing it at single precision, and its clustering read
that file, so the clustering input in those results was the distance matrix
quantized to three decimals. On small or well-separated datasets this makes no
difference; on larger ones it can change the partition, the leaf order and the
resulting ARI (ADR-061). A comparison with a historical figure should take this
into account.

---

## 8. Illustrative result

One quantitative result is stated here, as an illustration of the procedure on
real data and not as a performance claim.

The eleven nucleotide residues of the B-DNA dodecamer **1BNA**, chain A, that
pass the completeness rule (the 5'-terminal DC 1 carries no phosphate) were
described by nucleotide TSR at θ = 29, d = 18 classes, one key set per residue.
Pairwise generalized-Jaccard distances were computed (§4.1), average linkage
applied (§5.2), and a flat partition into four clusters taken by an explicit
`maxclust` cut (§5.4): *t* was set to the number of distinct bases present, so
the result measures agreement at that granularity (§5.4), not the discovery of
four groups. Class labels were
the base identities DA, DC, DG, DT.

The four clusters correspond exactly to the four bases: **adjusted Rand index
= 1.0**.

On the 75 amino-acid residues of **1UBQ** chain A that pass the completeness
rule, described by amino-acid atomic TSR at the same settings, the same
procedure with *t* set to the number of distinct residue types (18) gives an
adjusted Rand index of **0.521**. The 75 residue key sets share **no key in
common** in the sense of §3.1.

Quantizing either distance matrix to three decimals, as the historical pipeline
did (§7), changes neither result.

These are single-structure observations. No claim is made that they generalise,
and no comparison against any other descriptor or method was performed.
