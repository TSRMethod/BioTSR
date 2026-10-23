# Outputs written before atomic labels came from the element field

These nine files are three real calculations written by the **previous
implementation** of BioTSR, an earlier development build (`1.0.0rc6`), which
labelled atoms by atom name through the historical atom lexicon and gave each
nucleotide of a `nucleotide_chain` calculation an exact atom count. The current
version labels atoms by element (ADR-113) and includes every nucleotide of a
chain with the atoms it has (ADR-114). The files are kept unedited, so the
tests can check which old outputs the current version reuses, which it
refuses, and that refusing leaves them untouched.

The commands were run in a folder holding a copy of
`tests/data/synthetic/SYNPHOS.pdb` (SHA-256 `71e2eaf8…`) and `1BNA_A_1-3.pdb`,
which is chain A residues 1–3 of `tests/data/structures/1BNA.pdb` followed by
its own `TER` line (SHA-256 `f0888b7f…`; the tests write the same bytes), with
the older source tree on `PYTHONPATH`:

```
biotsr run --method drug --pdb SYNPHOS.pdb --chain A --entity-name LIG \
    --entity-id 101 --sample-id SYNPHOS_A_LIG101 --output out --no-triplets
biotsr run --method nucleotide_chain --pdb 1BNA_A_1-3.pdb --chain A \
    --sample-id 1BNA_A_1-3_chain --output out --no-triplets
biotsr run --method protein --pdb SYNPHOS.pdb --chain A \
    --sample-id SYNPHOS_A_protein --output out --no-triplets
```

| calculation | then | now |
|---|---|---|
| `drug`, LIG 101 | 8 points, 56 triangles, atomic key space `1187735d…` | the same keys, in the element key space: not reused, never compared with current results |
| `nucleotide_chain`, 1BNA A 1–3 | 41 points (DC 1 excluded by its count), 10,660 triangles | 57 points, 29,260 triangles: not reused |
| `protein` (Cα) | 14 points, 364 triangles, key space `da39d735…` | unchanged: reused by `--if-exists skip` |

| file | SHA-256 |
|---|---|
| `SYNPHOS_A_LIG101.keys` | `ee12d9567a6ed5f2bbe485b19e305534b3c740bc56b59b6b69d9d0086fa962c6` |
| `SYNPHOS_A_LIG101.meta.json` | `8b4f829926a6d46a4dcb48aa19dac3f73ea5b96151f54e3f5b629ea1a0886164` |
| `SYNPHOS_A_LIG101.run.json` | `0bd4b5291b1e02a8c6948cc189c5759265c99c129799130e8c310787ff101d4d` |
| `1BNA_A_1-3_chain.keys` | `f1128f301a09878bcef95df7770d68a72a55ede28a2d828da9f5372c21ebb119` |
| `1BNA_A_1-3_chain.meta.json` | `54142133b0617d3181c1e6e78276fd32bc3e1baaedb720aaa8d815e4c487b54c` |
| `1BNA_A_1-3_chain.run.json` | `169d2b6f1b3d1ee9c437ae7232a408686ef35d78b0377f22c47d6dab3740cbf0` |
| `SYNPHOS_A_protein.keys` | `2139564d8e13f73778fcf2b42f6b689a3f7f63bf5c817ec62cee7a7113bc789c` |
| `SYNPHOS_A_protein.meta.json` | `d9558665aef17dba58a2e595ecf1882e645ce68cccbaeac46104bdede9477752` |
| `SYNPHOS_A_protein.run.json` | `1f223955dc4eee604ea81e12d23f8133d29e901d17cde80f181df8b2ebb3ca5d` |

The tests that read these files are in
`tests/integration/test_previous_outputs_element_labels.py`.
