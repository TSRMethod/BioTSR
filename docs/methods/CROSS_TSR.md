# Cross-TSR

This page describes the cross-entity TSR calculation as implemented in BioTSR
1.0: what is computed and how. It states the scientific procedure; where a
choice had alternatives, it says which alternative was taken and what
distinguishes it. Implementation details that do not change the result — module
layout, data structures, serialization formats — are deliberately absent; see
[../ARCHITECTURE.md](../ARCHITECTURE.md) for those.

It makes **no claim about biological performance or scientific superiority.**
Nothing here has been demonstrated to classify, cluster or predict anything
better than anything else, and nothing here says so. What is stated is what is
computed and how.

For the commands that run each cross method, see
[../REPRODUCIBILITY_NOTE.md](../REPRODUCIBILITY_NOTE.md).

---

## 1. The cross-TSR concept

Ordinary TSR describes one molecular object by the triangles formed among its own
points. **Cross-TSR describes an interaction** by the triangles formed *between*
two objects.

Two molecular sides are specified by the user — always two, always explicitly:

| cross type | side A | side B |
|---|---|---|
| protein–protein | one protein chain | a second protein chain |
| drug–protein | one ligand (a named residue instance) | one protein chain |
| carbohydrate–protein | one glycan residue, or a listed vocabulary of them | one protein chain |
| nucleotide / DNA / RNA–protein | one nucleotide residue, or a strand | one protein chain |

Both participating chains are named by the user. The calculation draws points
from those two sides only; no other chain of the structure takes part. The two
chains of a protein–protein cross must differ; a ligand, glycan or nucleotide
side may share a chain identifier with the protein. In every mode the two sides
must be disjoint: a calculation whose sides would share an atom is refused,
naming the atoms, and nothing is removed from either side to make it run
(ADR-109).

A cross triangle always spans the interface: it has vertices on both sides, never
on one. This is what makes the descriptor a description of the *interaction*
rather than of either partner.

---

## 2. Atomic contact definition

For the atomic cross modes, the interface is defined atom by atom over
**non-hydrogen atoms only**.

Let *A* be the set of non-hydrogen atoms of side A and *B* that of side B, and let
*c* be the contact cutoff in ångströms. Then

> an atom **a ∈ A** belongs to the contact pool **iff** there exists **b ∈ B**
> with **d(a, b) ≤ c**;
> an atom **b ∈ B** belongs to the contact pool **iff** there exists **a ∈ A**
> with **d(a, b) ≤ c**.

Three properties of this rule are worth stating explicitly, because each has a
plausible alternative that was **not** taken:

1. **It is symmetric.** Both sides are narrowed by the cutoff. Neither side enters
   whole.
2. **It acts at the atom level.** An atom is not admitted because a *different*
   atom of its residue, ligand or nucleotide makes contact. If a glycan residue
   has atoms E1, E2, E3 and only E1 has a qualifying protein contact, only E1
   becomes a vertex.
3. **One qualifying partner is enough.** An atom needs a single opposite-side
   atom within the cutoff, not a majority of them and not the nearest one.

Distances are Euclidean, computed from the deposited coordinates in double
precision, and compared with `≤`.

---

## 3. Protein–protein cross, all-atom

Two user-specified protein chains. Non-hydrogen atoms of each. Section 2's rule
applied directly: an atom of either chain is a vertex iff it has at least one
non-hydrogen atom of the other chain within the cutoff.

Interface residues are *not* admitted whole. If atoms A1 and A2 of one residue
have qualifying contacts and A3 does not, A1 and A2 are vertices and A3 is not.

Points are labelled by element (the element field, columns 77–78; radix 112,
ADR-113).

A protein side reads the chain's `ATOM` records and, of the modified amino
acids that structure files write as `HETATM` records, exactly phosphothreonine
(TPO), phosphoserine (SEP) and phosphotyrosine (PTR): the modified residues
with a defined heavy-atom count, labelled by element like every atom. They enter under the
same hydrogen, alternate-location and contact rules as every other atom, and
none is converted to its unmodified residue. Any other modified residue stored
as `HETATM` has no defined atom count and is not read; it is named with the
result (`not supported`) rather than dropped silently (ADR-101). The protein
side of the three entity cross methods (section 5) follows the same rule.

Every protein side reads its chain up to the chain's `TER` record. A residue
named `UNK` in a selected side — anywhere in a whole chain, or among named
residues — refuses the calculation: no calculation is done for UNK (ADR-109).
No completeness rule applies to a cross side, so a terminal residue's atoms,
including `OXT`, take part like any other (ADR-111). An atom with a blank
element field among a side's atoms refuses the calculation, because it cannot
be told from a hydrogen (ADR-108).

---

## 4. Protein–protein cross, Cα

The Cα mode uses a **different selection projection**, and the difference is the
scientifically important part of this section. It is two steps, and they ask two
different questions.

**Step 1 — which residues are at the interface.** Residue *R* on chain A and
residue *S* on chain B are interface residues if **at least one pair of their
non-hydrogen atoms** satisfies *d ≤ c*. Membership is decided over **all** heavy
atoms, side chains included — not over Cα atoms.

**Step 2 — how an interface residue is represented.** Each interface residue
contributes **its Cα** as the TSR point.

> **The Cα itself does not have to be within the cutoff**, and for a long side
> chain it usually is not.

```
        heavy-atom contact  (any non-hydrogen pair, d ≤ c)
                 ↓
          interface residue
                 ↓
         represented by its Cα
```

Why the distinction matters: a lysine can reach the opposite chain through its NZ
while its own Cα sits 8–12 Å away. Applying the cutoff to Cα atoms instead would
silently discard exactly the long-side-chain contacts that characterize most
protein–protein interfaces. The two rules are not approximations of each other;
on a constructed test case with a 3 Å cutoff, a 2 Å NZ–OE1 contact and a 12 Å
Cα–Cα separation, the standard rule selects both residues and a Cα–Cα rule
selects neither.

Cα points are **residues**, so they are labelled with the residue labels (the
20 regular amino acids and eleven modified residues, radix 30) on the 35-class
Cα distance scale — the same key space as ordinary protein Cα TSR, so Cα cross
keys and ordinary Cα keys are directly comparable. They are not encoded with
the atomic element labels.

If an interface residue carries no usable Cα, it contributes no point and the
omission is recorded; no other atom is substituted for it.

The residue labels give eleven modified residues codes of their own (TPO,
SEP, PTR, NIY, M3L, TRF, CSO, LLP, CME, CSD; MSE shares methionine's), and the
Cα mode reads them (ADR-107): a `HETATM` residue of one of these names whose CA
atom's element is carbon or blank has a Cα point, and step 1 counts its heavy
atoms, so it can be an interface residue like any other. One without a readable
Cα is still an interface residue if it makes contact; like any interface
residue without a Cα it contributes no point, and it is named with its reason. Any other modified
residue stored as `HETATM` contributes no point and no contact and is named
with the result (`not represented`).

When a side is narrowed to listed residues, the list follows the residue
grammar shared by every method: PDB numbers with an optional insertion code,
and ascending ranges of plain numbers that may not span an insertion-code
residue the list does not name (ADR-102). **Narrowing applies before the
interface is found** (ADR-109): step 1 compares the heavy atoms of the selected
residues of one side with those of the selected residues (or the whole chain)
of the other, so a residue is never at the interface only because it touches
residues the caller left out. On 1BRS, barnase residue 27 against barstar
chain D gives three interface residues of D and three triangles; before
1.0.0rc6 the partner interface was taken from the whole of chain A.

---

## 5. Drug–protein, carbohydrate–protein and nucleotide–protein cross

All three use the atomic rule of section 2, symmetrically, on the two
user-specified sides:

The entity side is named by its chain, which is always required and is the only
chain searched, and by the entity's name. With a residue identifier (the PDB
residue number and an optional insertion code, e.g. `404` or `52A`) it is that
**one copy**. Without one it is **every copy of that name on that chain**, taken
together as one side — even when there is a single copy — so the contact pool
and the cross triangles may draw on several copies (ADR-096). Every copy used is
recorded with the result, and the side is described as, for example,
`entity=NAG(all copies)`. The polymer scopes (a listed vocabulary of glycan
names, or a whole strand) pool every listed residue of the chain by design.

* **drug–protein** — one named ligand copy, or every copy of that name on its
  chain, against one protein chain. A
  ligand atom is retained only if it has a qualifying protein atom within the
  cutoff, and a protein atom only if it has a qualifying ligand atom. Where an
  explicit ligand atom selection is supplied, that selection is applied **first**
  and the selected atoms are then contact-filtered.
* **carbohydrate–protein** — one glycan residue, every copy of one
  monosaccharide name on the chain, or a listed vocabulary of monosaccharide
  names for a polymer scope, against one protein chain.
* **nucleotide / DNA / RNA–protein** — one nucleotide residue, every copy of one
  nucleotide name on the chain, or a whole strand, against one protein chain.

The protein sides are not identical in which residues they read. The protein
side of drug–protein (and of protein–protein) cross reads every `ATOM` residue
of its chain. The whole-chain protein side of carbohydrate–protein and
nucleotide–protein cross reads the 20 regular amino acids and TPO, SEP and PTR,
the vocabulary of those methods' historical contact programs; any other residue
of the chain is named with the result as `not supported`, and the result is
partial. Narrowing a protein side with residue types, residues or atoms reads
exactly what is named, in every mode. Each protein side stops at its chain's
`TER`, and the molecule side reads the whole of model 1, so a glycan or
nucleotide written after the protein's `TER` is still selectable.

For the carbohydrate and nucleotide modes in particular: **only the contacting
atoms of the entity residue are retained.** The residue is not expanded to all of
its atoms because one of them touches the protein.

Points on every side are labelled by element (radix 112, ADR-113).

---

## 6. Cross-triplet generation

Once the two contact pools are fixed, every mixed-side triple of pooled points is
enumerated. With *a* points in pool A and *b* in pool B, the candidates are

```
    A + A + B   and   A + B + B
```

and their number is exactly

```
    a·C(b,2) + b·C(a,2)  =  C(a+b,3) − C(a,3) − C(b,3)
```

Triples with all three vertices on one side (A+A+A, B+B+B) are **not** generated:
they describe one partner, not the interaction.

If no atom pair across the two sides lies within the cutoff, both pools are
empty and no triangle can be formed. That is a valid outcome, not an error: the
result has no keys and carries a warning beginning `no triangles:` that says
why.

**The contact cutoff is not re-imposed on the resulting triangles.** A vertex
entered the pool because it had *one* qualifying partner; the triangle built from
it may have edges far longer than the cutoff, and such triangles are kept. The
cutoff is a rule about pool membership and nothing else.

Optional **size filtering** — rejecting a triangle whose longest edge reaches a
threshold — is a separate operation applied after the geometry is known. It is
**off by default** and is never conflated with the contact cutoff: the cutoff
decides which atoms exist as vertices, the size filter decides which formed
triangles are kept.

---

## 7. Cutoff values

The contact cutoff is a **user-configurable parameter**, not a constant of the
method. The values below are the reference values used by the standard profiles
and are reported as defaults, not as properties of the science:

| cross type | reference cutoff |
|---|---|
| protein–protein (all-atom and Cα) | 5 Å |
| drug–protein | 3.5 Å |
| carbohydrate–protein | 3.5 Å |
| nucleotide / DNA / RNA–protein | 3.5 Å |

Every calculation records the cutoff actually used alongside its results, so a
reported key frequency is never separable from the sampling that produced it.

---

## 8. Hydrogen atoms

Hydrogen atoms are **excluded from every stage** of every cross calculation:
contact detection, interface-residue identification, contact pools and TSR vertex
generation alike. Deuterium is excluded on the same rule. Element assignment is
read from the PDB element field (columns 77–78), without regard to letter case;
an atom with a blank element field cannot be classified and refuses the
calculation (ADR-108).

---

## 9. TSR calculation

From the pooled, mixed-side triples onward, cross-TSR uses the **same** procedure
as ordinary TSR, with no cross-specific geometry, ordering, binning or encoding:

* the three vertices are put in a canonical order by the shared four-case rule
  over their labels and edge lengths, so that a triangle has one representation
  regardless of the order its vertices were enumerated in;
* an angle θ and a length are derived from the canonical triangle;
* θ is binned into 29 classes and the length into the scheme belonging to the
  point type — 18 classes for atomic points, 35 for Cα points;
* the label triple and the two bin classes are combined into a single integer key
  by the mixed-radix formula, at radix 112 for atomic points and radix 30 for Cα
  points.

These conventions are shared with non-cross TSR and are defined once for all
methods (see [../ARCHITECTURE.md](../ARCHITECTURE.md) and
[../architecture/DECISIONS.md](../architecture/DECISIONS.md)); they are
referenced rather than restated here so that the two cannot drift apart. All
computation is in **double precision (float64)**.

The output of a cross calculation is a frequency distribution over integer keys,
optionally accompanied by the per-triangle records that produced it. Each triangle
record carries the full physical identity of each vertex — structure, model,
chain, residue name, residue number, insertion code, atom name, alternate
location — and the side it came from, so any key can be traced back to the atoms
that generated it.

---

## 10. Relationship to the historical implementations

The cross calculations described above were previously implemented as four
separate programs, one per molecular pair. BioTSR does not reproduce their
output. The standard described on this page differs from them in ways that are
individually documented and individually measured:

| aspect | historical | standard |
|---|---|---|
| chain scope | two lineages scoped each side to a chain; two scanned the whole structure for the protein side | both chains always user-specified |
| contact filtering | one lineage narrowed only the protein side, admitting the whole ligand | symmetric on both sides |
| atom identity in the pool | two lineages keyed pool membership on the atom *name*, so two identically named atoms of different residues collapsed into one | full physical identity; distinct atoms cannot collapse |
| alternate locations | the contact generators applied no altLoc test, so both conformers of an atom could enter | altLoc blank or `A` only |
| distance classes | one lineage used 17 | 18 for all atomic modes |
| label ordering | one lineage compared labels as strings | numeric |
| median length | one lineage derived it from rounded edge lengths (Apollonius) | from the base midpoint's coordinates |
| angle divisor | most lineages used the literal 3.14 | π |
| numerical precision | one lineage worked from single-precision coordinates | double precision |

Each of these corrects a sampling, identity or numerical issue rather than
changing the definition: the descriptor being computed is the same descriptor.
Two of them are measured examples worth stating: on 6G9M chain A with ligand
ESW 404 at 3.5 Å, the one-sided historical rule pools 37 ligand atoms where the
symmetric rule pools 10 (8,325 triangles against 900); and on 1KX5, name-only
pool identity collapses two distinct lysine NZ atoms into one, under-counting
the triangles of the affected contact groups by 44%.

Every result records the conventions that produced it in its two fingerprinted
identities, so results computed under different key conventions cannot be
combined. The decision records for these choices, with their derivations, are
ADR-038 to ADR-053 and ADR-096 (molecule instances) in
[../architecture/DECISIONS.md](../architecture/DECISIONS.md).
