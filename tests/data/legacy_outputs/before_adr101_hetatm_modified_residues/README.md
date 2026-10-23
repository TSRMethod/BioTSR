# An output written under the older modified-residue rule

These three files are a real calculation written by an **older implementation**
of BioTSR: an earlier development build that also reported version
`1.0.0rc4`. That implementation did not read the modified amino acids TPO, SEP
and PTR when a file writes them as `HETATM` records (ADR-101 changed this).
The files are kept unedited, so the tests can check that the current version
never reuses them as if they were current (ADR-105).

They were produced from the synthetic structure
`tests/data/synthetic/SYNPHOS.pdb` (SHA-256 `71e2eaf8…`, the same file the
tests use). The command was run in a folder holding a copy of that file, with
the older source tree on `PYTHONPATH`:

```
biotsr run --method aa_ind_specific --pdb SYNPHOS.pdb --chain A \
    --residues 'SER_4;TPO_3;TYR_6' --sample-id SYNPHOS_A_specific \
    --output out --no-triplets
```

The selection names two ordinary residues and one modified residue. The older
rule read 18 atoms (SER 4 and TYR 6) and 816 triangles, and warned
`requested but not found: TPO_3`. The current rule reads TPO 3 as well and
finds 3654 triangles for the same request.

The completion record has no `processing` entry, and the `.meta.json` source
record has no `modified_residue_policy`. Neither existed then. The request is
otherwise identical to the one the current version builds, down to the
structure's SHA-256 and the version number. So a check that compares only the
request, or only the version, would keep this output.

| file | SHA-256 |
|---|---|
| `SYNPHOS_A_specific.keys` | `a0c8a6e53cad47f598e9a29bb4c850f3ba2372e41822de79f33684989600d310` |
| `SYNPHOS_A_specific.meta.json` | `9ddee469713437abc6853aa1b7066627fcf15bb6df9fcc24b4dd3f23cf035e5c` |
| `SYNPHOS_A_specific.run.json` | `61210c5455ee51e48c5b09f2e022d31398ef9d64fdd54737deb0c7f73c34a2be` |

The tests that read these files are in `tests/integration/test_processing_identity.py`.
