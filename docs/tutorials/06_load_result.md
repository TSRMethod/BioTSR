# 06 — Load a saved result

Native output can be read back into a `TSRResult`, with both identity
fingerprints restored, so you can analyse results without recomputing them.

```bash
biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
```

```text
results/1UBQ_A_protein.run.json
```

The command prints the path of the completion record; the result's files sit
beside it under the same prefix. Load a result from its `.keys` file (or its
`.meta.json` sidecar):

```python
import biotsr

loaded = biotsr.load_result("results/1UBQ_A_protein.keys")
fresh = biotsr.analyze("1UBQ.pdb", method="protein", chain="A")

assert loaded.frequencies == fresh.frequencies
assert loaded.key_space.fingerprint == fresh.key_space.fingerprint
assert loaded.sampling.fingerprint == fresh.sampling.fingerprint

again = biotsr.load_result("results/1UBQ_A_protein.meta.json")
assert again.frequencies == loaded.frequencies
```

A loaded result does not carry the per-triangle records or the order in which
keys were first seen, and it says so: operations that need them refuse a loaded
result instead of guessing.

```python
assert loaded.triplets is None
```

To see what produced a calculation without loading it, describe its prefix
(or its `.run.json`):

```bash
biotsr describe results/1UBQ_A_protein
```

```text
command          biotsr run
method           protein
prefix           1UBQ_A_protein
complete         True
```

and to list every calculation in a folder, describe the folder:

```bash
biotsr describe results/
```

```text
1UBQ_A_protein protein 70300 triangles complete
```

Every result this version writes can be loaded this way. The loader refuses a
key file paired with another result's sidecar, or a file that was cut short,
rather than guessing. Outputs written by a release candidate before 1.0.0rc4
— one directory per sample, with names such as `1UBQ_A.keys_29_35` — are
refused with a message saying they must be recomputed.
