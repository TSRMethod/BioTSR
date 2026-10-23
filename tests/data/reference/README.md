# Fixed reference data for the test suite

## `standard_baseline.json`

Standard results for 42 calculations over the structures in `../structures/`,
covering all 20 methods and the supported modifiers (AA grouping, mirror, size
filter, SSE annotation). They were computed with BioTSR 1.0.0rc1, the last build
before the historical compatibility code was removed, on macOS arm64 with
CPython 3.13.2; the four largest protein C-alpha cases (2AC0, 3EJ1, 6G9M, 7LF9
chain A) carry the keys-file hashes the release checks recorded for the same
candidate.

For each calculation the file records the SHA-256 of the native keys file, the
point, triangle and distinct-key counts, both identity fingerprints, the
residues used and the category totals, and, for calculations that kept their
triangle records, the SHA-256 of the native triplet file on that platform.
`tests/integration/test_standard_references.py` recomputes every case and requires
the same values. The file is never regenerated to match changed code.

The atom and residue selections in it were chosen for the tests from the public
structures; none was copied from any dataset.

## `calpha_base30.json`

The current expected values of the 17 protein Cα cases of
`standard_baseline.json`, whose own values stay frozen at rc1. Since 1.0.0rc5
the Cα residue labels are encoded at radix 30 and the Cα methods read the
supported modified residues (ADR-107), so every Cα key has another value. The
values here were derived without the changed code: the previous build
recomputed each case and reproduced the frozen rc1 files, every key was decoded
at its old radix (20, or 13 under grouping) and re-encoded at radix 30, and for
6G9M chain A — whose CME 161 is now a point — the 57,970 triangles containing
it were added from an independent calculation. The current package reproduces
every keys and triplet file byte for byte. `rc1_radix` is the radix the frozen
value was written at; `test_standard_references.py` also maps the current keys
of the old selection's triangles — those whose residues are all regular amino
acids, by identity — back to it and requires the rc1 file, and counts the
added triangles separately (`added_triangles`).

## `atomic_element_labels.json`

The current expected values of the 25 atomic cases of `standard_baseline.json`,
whose own values stay frozen at rc1. Since 1.0.0rc7 every atomic method labels
an atom by its element field (ADR-113), which gives every atomic key space a
new fingerprint, and `nucleotide_chain` includes each nucleotide with the atoms
it has (ADR-114). Each value was compared with rc1 before it was written, and
`test_standard_references.py` repeats the comparison: apart from the key-space
fingerprint, every value is rc1's, except that the three 6G9M ESW triplet files
spell chlorine's label name `Cl` (rc1: `CL`; `rc1_label_names` spells it back,
which gives rc1's file) and `nucleotide_chain_1BNA_A` includes DC 1 (243
points, 227 in rc1). That case's keys were derived without the changed code —
the previous build (1.0.0rc6) computed the same chain with its count rule
switched off — and its triangles without a DC 1 atom reproduce the rc1 keys
file. Triplet hashes are given only where the test compares them (120,000
triangles or fewer).

## `historical/1UBQ_A_TSR-Package_TSR_keys.txt.gz`

One column: the integer key of each of the 70,300 triangles of 1UBQ chain A, in
candidate order, as written by `TSR-Package/tsr_package/TSR.py` (commit
`fd3c425a`; TSR-Package is MIT-licensed, © 2024 Poorya Khajouie) when run on
`../structures/1UBQ.pdb`. It was extracted from the first column of that run's
`1UBQ.triplets_29_35` output (SHA-256 `c07056c1654e…6f`). The decompressed column
has SHA-256 `4a33134c42dce5e38da8aff585ded1dec20f37157e8cec83ac111afc69ef676f`.

`tests/characterization/test_1ubq_historical_difference.py` uses it to show that
the standard result differs from the historical program's in exactly the
documented way (the angle divisor) and in no other.
