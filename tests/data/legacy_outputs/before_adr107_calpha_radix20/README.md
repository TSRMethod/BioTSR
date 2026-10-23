# A protein Cα output written before the residue labels moved to radix 30

These three files are a real calculation written by the **previous
implementation** of BioTSR, an earlier development build
(`1.0.0rc4`). It labelled Cα points with the 20 regular amino-acid codes at
radix 20 and read no modified residue (ADR-107 changed both). The files are
kept unedited, so the tests can check that the current version never reuses
them and never compares them with current Cα results.

They were produced from the synthetic structure
`tests/data/synthetic/SYNPHOS.pdb` (SHA-256 `71e2eaf8…`, the same file the
tests use). The command was run in a folder holding a copy of that file, with
the older source tree on `PYTHONPATH`:

```
biotsr run --method protein --pdb SYNPHOS.pdb --chain A \
    --sample-id SYNPHOS_A_protein --output out --no-triplets
```

The older implementation read 10 Cα points and 120 triangles, warned
`not represented: TPO_3, SEP_5, PTR_7, CME_8 (…)` and exited 4 (partial). The
current version reads 14 points, TPO 3, SEP 5, PTR 7 and CME 8 included, and
finds 364 triangles for the same request, with no warning.

Its completion record states the processing rules it was computed under, with
the modified-residue rule `hetatm_amino_acids:PTR,SEP,TPO` (then one rule for
every method), and its result's key space is the radix-20 one
(`amino_acid_20`).

| file | SHA-256 |
|---|---|
| `SYNPHOS_A_protein.keys` | `5fee880c45e88d9e04f4dd35896a22d44ab61403d0de762bd766981d0d433a0b` |
| `SYNPHOS_A_protein.meta.json` | `832368ae4d9014f6ba06613474f5b28e667ce0321789870d93226f9fb2018cc6` |
| `SYNPHOS_A_protein.run.json` | `3c5b64e012cfbf15a7f51ab2a6cbcac0f0ecba03a8220cea23171bd925fbbd2a` |

The tests that read these files are in
`tests/integration/test_modified_residue_calpha.py`.
