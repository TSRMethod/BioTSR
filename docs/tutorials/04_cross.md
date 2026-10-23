# 04 — Cross TSR

A cross method draws triangles across **two sides**, such as a ligand and the
protein around it. Every triangle has vertices on both sides. An atom joins the
pool only if at least one heavy atom of the other side lies within the contact
cutoff.

Both chains are always named: `--chain` is the ligand's chain and `--chain-2`
the protein's.

```bash
biotsr run --method drug_protein_cross --pdb 3PTB.pdb --chain A --chain-2 A --entity-name BEN --entity-id 1 --output results/
```

```text
method     drug_protein_cross
profile    standard_atomic_cross_v1
triangles  84
keys       72 distinct
output     results/3PTB_A_A_drug_protein_cross_BEN_1
```

The prefix names both chains in argument order, then the method and the
ligand copy.

## One copy, or every copy on the chain

The ligand, carbohydrate or nucleotide side of a cross method selects its
molecule by the same rule as the single-entity methods
([tutorial 02](02_atomic.md)): `--entity-name` with `--entity-id` is that one
copy; `--entity-name` alone is every copy of that name on `--chain`, pooled
into one side. With the glycan on 7LF9 chain D against protein chain B:

```bash
biotsr run --method carbohydrate_protein_cross --pdb 7LF9.pdb --chain D --chain-2 B --entity-name NAG --entity-id 1 --output results/
biotsr run --method carbohydrate_protein_cross --pdb 7LF9.pdb --chain D --chain-2 B --entity-name NAG --output results/
```

```text
triangles  160
output     results/7LF9_D_B_carbohydrate_protein_cross_NAG_1
triangles  4185
output     results/7LF9_D_B_carbohydrate_protein_cross_NAG_allcopies
```

The first side of the second run is NAG 1, 2 and 5 of chain D, together;
the protein side is chain B, and no other chain is searched for either.

## From Python

```python
import biotsr

cross = biotsr.analyze_cross("3PTB.pdb", method="drug_protein_cross",
                             chain="A", chain_2="A",
                             entity_name="BEN", entity_id="1")
assert cross.n_triangles == 84

# Every NAG on chain D as the carbohydrate side; the copies used are recorded.
glycan = biotsr.analyze_cross("7LF9.pdb", method="carbohydrate_protein_cross",
                              chain="D", chain_2="B", entity_name="NAG",
                              keep_triplets=False)
side_a, side_b = glycan.annotations["selection_request"]["sides"]
assert side_a["scope"] == "all_copies"
assert side_a["instances"] == ["NAG_1", "NAG_2", "NAG_5"]
assert side_b["chain"] == "B" and side_b["scope"] == "whole_chain"
assert "entity=NAG(all copies)" in glycan.sampling.side_a_selector

# The contact cutoff is a processing setting: it changes which atoms are
# drawn, never what a key means.
wider = biotsr.analyze_cross("3PTB.pdb", method="drug_protein_cross",
                             chain="A", chain_2="A",
                             entity_name="BEN", entity_id="1",
                             contact_cutoff=4.5)
assert wider.key_space.fingerprint == cross.key_space.fingerprint
assert wider.sampling.fingerprint != cross.sampling.fingerprint
```

The contact cutoff selects atoms. It does not limit triangle edges; use
`--size-filter` for that. Because it changes which triangles are drawn,
results with different cutoffs are compared only when you acknowledge the
difference (`--allow-mixed-processing`, [tutorial 07](07_analysis.md)); results
that differ only in what they selected — one NAG or all of them, chain B or
another chain — compare freely.
