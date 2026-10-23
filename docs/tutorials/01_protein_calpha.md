# 01 — Protein Cα TSR

The protein Cα method draws every triangle over the Cα atoms of one chain,
labels each vertex by its residue — one of the 20 regular amino acids or one of
eleven supported modified residues, such as phosphoserine `SEP` — and encodes each triangle's shape — an
angle class and a distance class — together with its three labels as one
integer **TSR key**.

## From the command line

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

```text
method     protein
profile    standard_protein_calpha_v1
triangles  70300
keys       45349 distinct
output     results/1UBQ_A_protein (4 files, completion record 1UBQ_A_protein.run.json)
```

76 residues give C(76,3) = 70,300 triangles, which produce 45,349 distinct keys.
The summary goes to stderr; stdout is the path of the completion record,
`results/1UBQ_A_protein.run.json`.

`--output` is one shared folder. The calculation writes its files directly
into it, under one **prefix** built from the request (structure, chain,
method):

| file | contents |
|---|---|
| `1UBQ_A_protein.keys` | each key and its frequency |
| `1UBQ_A_protein.triplets` | one row per triangle: vertices, geometry, classes, key |
| `1UBQ_A_protein.meta.json` | the key space (with both bin definitions) and sampling profile, with their fingerprints |
| `1UBQ_A_protein.run.json` | the completion record: the exact command, the request and every file with its size |

A second calculation into the same folder adds its own files beside these,
under its own prefix.

## The three Cα modifiers

The modifiers are flags on the same command, and they combine freely:

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output modified/ --aa-grouping --mirror --size-filter 25.0
```

* `--aa-grouping`: the grouped residue labels (13 groups of the regular amino
  acids, MSE with MET, each other modified residue on its own) instead of the
  individual ones, at the same radix 30;
* `--mirror`: negate a key whose angle before folding exceeds 90°;
* `--size-filter 25.0`: keep only triangles whose longest edge is under 25 Å.

The modifiers are part of the prefix, so a modified run never shares a name
with the plain one, even in the same folder:

```text
output     modified/1UBQ_A_protein_aagroup_mirror_size25 (4 files, completion record 1UBQ_A_protein_aagroup_mirror_size25.run.json)
```

## From Python

```python
import biotsr

result = biotsr.analyze("1UBQ.pdb", method="protein", chain="A")
assert result.n_triangles == 70300
assert result.n_unique_keys == 45349

grouped = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", aa_grouping=True)
assert result.key_space.label_scheme.radix == grouped.key_space.label_scheme.radix == 30

# Different label schemes mean the integers are not comparable, even at the
# same radix -- and the fingerprints say so.
assert grouped.key_space.fingerprint != result.key_space.fingerprint
```

Every command-line flag is the `analyze()` keyword with dashes for underscores
(`--size-filter` ⇄ `size_filter`).
