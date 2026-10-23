# 02 — Atomic TSR

Atomic methods use individual heavy atoms as vertices, each labelled by its
element — the PDB element field, columns 77–78, never the atom name — from
BioTSR's packaged element labels. You never supply them; they ship with the
package.

## One amino-acid residue

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --residue 15 --residue-name LEU --output results/
```

```text
method     amino_acid_atomic
profile    standard_aminoacid_atomic_v1
triangles  56
keys       46 distinct
output     results/1UBQ_A_amino_acid_atomic_LEU_15
```

Leucine has 8 heavy atoms, so C(8,3) = 56 triangles. `--residue` with
`--residue-name` is that one residue. With `--residue-name` alone the method
runs **once per residue of that name** — here one calculation for each of the
nine leucines of 1UBQ, written as members of one sweep — and with neither, once
per residue of the chain:

```bash
biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --residue-name LEU --output results/ --no-triplets
```

```text
entities   9
triangles  504
output     results/1UBQ_A_amino_acid_atomic_perresidue-LEU
```

## One ligand

A ligand copy is named by its chain, residue name and residue identifier:

```bash
biotsr run --method drug --pdb 3PTB.pdb --chain A --entity-name BEN --entity-id 1 --output results/
```

```text
method     drug
profile    standard_drug_v1
triangles  84
keys       38 distinct
output     results/3PTB_A_drug_BEN_1
```

## One copy, or every copy on the chain

A molecule — ligand, sugar or nucleotide — is selected by chain, name and
identifier, and the chain is always required: BioTSR never infers it and never
searches another chain.

* `--entity-name` with `--entity-id` is **that one copy** (prefix `…_NAG_1`);
* `--entity-name` alone is **every copy of that name on the chain, in one
  combined calculation** (prefix `…_NAG_allcopies`), even when the chain holds
  only one copy. Triangles may span copies; that is the intended calculation.

7LF9 chain D carries one glycan, NAG 1 – NAG 2 – BMA 3 – MAN 4 – NAG 5 – GAL 6
– MAN 7. Each NAG on its own has 14 heavy atoms, C(14,3) = 364 triangles;
without an identifier the three NAGs (1, 2 and 5) are pooled, 42 atoms and
11,480 triangles:

```bash
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id 1 --output results/
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id 2 --output results/
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --entity-id 5 --output results/
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name NAG --output results/
```

```text
output     results/7LF9_D_carbohydrate_NAG_1
output     results/7LF9_D_carbohydrate_NAG_2
output     results/7LF9_D_carbohydrate_NAG_5
triangles  11480
output     results/7LF9_D_carbohydrate_NAG_allcopies
```

Every copy of one sugar is **not** the whole glycan. `--entity-name MAN` is
MAN 4 and MAN 7 only; `carbohydrate_chain` is a different method that pools
every residue of the chain, different sugars included:

```bash
biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --entity-name MAN --output results/
biotsr run --method carbohydrate_chain --pdb 7LF9.pdb --chain D --output results/
```

```text
triangles  1540
output     results/7LF9_D_carbohydrate_MAN_allcopies
triangles  102340
output     results/7LF9_D_carbohydrate_chain
```

All of these calculations sit side by side in the one `results/` folder, one
line each:

```bash
biotsr describe results/
```

```text
7LF9_D_carbohydrate_NAG_1 carbohydrate 364 triangles complete
7LF9_D_carbohydrate_NAG_allcopies carbohydrate 11480 triangles complete
7LF9_D_carbohydrate_MAN_allcopies carbohydrate 1540 triangles complete
7LF9_D_carbohydrate_chain carbohydrate_chain 102340 triangles complete
```

## From Python

```python
import biotsr

leu15 = biotsr.analyze("1UBQ.pdb", method="amino_acid_atomic", chain="A",
                       residue="15", residue_name="LEU")
assert leu15.n_triangles == 56
assert leu15.key_space.label_scheme.radix == 112
assert leu15.key_space.distance_scheme.radix == 18

ben = biotsr.analyze("3PTB.pdb", method="drug", chain="A",
                     entity_name="BEN", entity_id="1")
assert ben.n_triangles == 84
assert ben.annotations["target"] == "BEN_1"

# Which molecules are in this chain? (identifier, name)
targets = biotsr.list_entity_targets("7LF9.pdb", chain="D", method="carbohydrate")
assert [t for t in targets if t[1] == "NAG"] == [("1", "NAG"), ("2", "NAG"), ("5", "NAG")]

# No identifier: every copy on the chain, in one calculation, and the copies
# actually used are recorded.
nag = biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D",
                     entity_name="NAG", keep_triplets=False)
side = nag.annotations["selection_request"]["sides"][0]
assert nag.annotations["target"] == "NAG_allcopies"
assert side["scope"] == "all_copies"
assert side["instances"] == ["NAG_1", "NAG_2", "NAG_5"]
assert nag.n_triangles == 11480                 # not 3 x 364: triangles span copies

# The same rule with a single copy: 3PTB chain A holds one BEN, and leaving
# the identifier out is still an all-copies selection, recorded as such.
all_ben = biotsr.analyze("3PTB.pdb", method="drug", chain="A",
                         entity_name="BEN", keep_triplets=False)
assert all_ben.annotations["target"] == "BEN_allcopies"
assert all_ben.frequencies == ben.frequencies   # the same atoms

# An argument the method does not read is refused, never ignored.
try:
    biotsr.analyze("3PTB.pdb", method="drug", chain="A", entity_name="BEN",
                   residues="1-10")
except biotsr.InvalidRequest as error:
    assert "does not use residues" in str(error)
else:
    raise AssertionError("an unused argument must be refused")
```

## Identifiers

An identifier is a **PDB residue identifier**: the residue number as written in
the file, which may be negative, with at most one insertion-code letter (`404`,
`52A`, `-5`). Leading zeros are dropped, so `005` is `5`, and `52` never
matches `52A`. A position ("the second NAG") or any other spelling is refused,
and so is an identifier the chain does not hold — naming the ones it does:

```python
assert biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D",
                      entity_name="NAG", entity_id="005",
                      keep_triplets=False).annotations["target"] == "NAG_5"

for entity_id, message in (("second", "is not a PDB residue identifier"),
                           ("3", "copies of NAG there have ids 1, 2, 5")):
    try:
        biotsr.analyze("7LF9.pdb", method="carbohydrate", chain="D",
                       entity_name="NAG", entity_id=entity_id)
    except biotsr.InvalidRequest as error:
        assert message in str(error), error
    else:
        raise AssertionError("%r must be refused" % entity_id)
```

The selection is recorded in the result's sidecar (per side: the chain, the
scope, the arguments, the copies used and any excluded). It is provenance only:
one copy against all copies never blocks a comparison
([tutorial 07](07_analysis.md)).
