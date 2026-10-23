# 03 — SSE-TSR

`--sse` annotates each Cα triangle with one of 18 secondary-structure
categories, read from the structure's own `HELIX` and `SHEET` records.

It is an **annotation**: the integer keys, their counts and both identity
fingerprints are exactly those of the same run without it. What changes is
which category each key's frequency is counted in.

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --sse --output results/
```

```text
sse        on, strand identity sheet_id_and_serial
triangles  70300
keys       45349 distinct
```

## From Python

```python
import biotsr

plain = biotsr.analyze("1UBQ.pdb", method="protein", chain="A")
annotated = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", sse=True)

# The same keys and the same identities ...
assert annotated.frequencies == plain.frequencies
assert annotated.key_space.fingerprint == plain.key_space.fingerprint
assert annotated.sampling.fingerprint == plain.sampling.fingerprint

# ... counted in secondary-structure categories instead of one bucket.
totals = annotated.category_totals()
assert sum(totals.values()) == annotated.n_triangles
assert set(plain.category_totals()) == {"regular"}
```

SSE is offered by the two protein Cα methods, `protein` and `aa_selection`,
and refused by name everywhere else: an atom, a ligand or a nucleotide has no
secondary structure of its own.
