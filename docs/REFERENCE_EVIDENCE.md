# Reference evidence

What the test suite compares BioTSR's results against, case by case, and what
kind of evidence each comparison is. Two kinds are kept apart throughout:

* **Regression evidence** shows that a result has not changed since a known
  build. It cannot show that the result is scientifically right; it shows that
  nothing moved.
* **Independent evidence** compares BioTSR with something it does not share code
  with — a calculation done by hand, or the output of another program — and so
  says something about correctness.

## 1. The standard reference calculations (regression evidence)

`tests/data/reference/standard_baseline.json` holds **42 calculations** over the
pinned public structures in `tests/data/structures/`, covering all 20 methods
and the four modifiers (AA grouping, mirror, size filter, SSE annotation).
`tests/integration/test_standard_references.py` recomputes each one and requires
the stored values.

**Where the values came from.** They were computed by BioTSR **1.0.0rc1** — the
last build before the historical compatibility code was removed — on macOS
arm64 with CPython 3.13.2, and frozen before any of that code was removed. The
keys-file hashes of the four largest protein C-alpha cases (2AC0, 3EJ1, 6G9M and
7LF9, chain A) are the ones the release checks had already recorded for that
candidate. The file is never regenerated to match changed code: a difference is
a finding to explain, not a value to update. Because the values come from an
earlier BioTSR build, they are regression evidence only.

**The protein Cα cases since 1.0.0rc5.** ADR-107 moved the Cα residue labels to
radix 30 and made the Cα methods read the supported modified residues, which
changes every one of the 17 Cα cases' keys. `standard_baseline.json` was not
edited; their current values are in `tests/data/reference/calpha_base30.json`.
Three kinds of evidence stand behind them, and the test checks each:

1. **Triangles of the old selection — through the remap.** The previous build
   recomputed every case and reproduced the frozen rc1 files byte for byte;
   every key was decoded at its old radix (20, or 13 under grouping) and
   re-encoded at 30. In the test, the current triangles whose vertices are all
   regular amino acids — decided by residue identity, not by label, because
   MSE shares MET's code — are mapped back to the old radix and must reproduce
   the rc1 keys file byte for byte.
2. **Triangles with a newly read modified residue — separately.** Only 6G9M
   chain A has one (CME 161). Its 57,970 triangles (C(341, 2)) were computed by
   an independent calculation that shares no code with BioTSR — its own
   geometry, ordering, binning and encoding, with the same residue tables at
   radix 30 — and added to the prediction. They are never mapped back to the
   old radix, where codes above 23 would collide. The test counts them. That
   independent calculation applied the same rule for *which* residue is read,
   so it confirms CME 161's triangles and keys, not the admission rule itself;
   the rule is covered by its own fixtures
   (`tests/integration/test_modified_residue_calpha.py`).
3. **The whole updated output — against the combined expectation.** The current
   keys file (and, on macOS arm64, the triplet file) of every case must equal
   the stored hash of the remapped old triangles plus the added ones.

Sampling fingerprints are unchanged; key-space fingerprints change.

**The atomic cases since 1.0.0rc7.** ADR-113 labels every atom by its element
field instead of its name, which gives every atomic key space a new
fingerprint, and ADR-114 lets `nucleotide_chain` include each nucleotide with
the atoms it has. `standard_baseline.json` was not edited; the 25 atomic
cases' current values are in `tests/data/reference/atomic_element_labels.json`.
The evidence for them, checked by the test:

1. **Everything but the key-space fingerprint is rc1's** for 24 of the 25:
   keys files, counts, residues, sampling fingerprints. The labels of the
   atoms in these cases did not change (their names and element fields agree),
   so the keys could not. A check reads the reference file against rc1 and
   refuses any other difference.
2. **Triplet files** are rc1's, or differ only in how a label name is spelled:
   the three 6G9M ESW cases name chlorine `Cl` where rc1 wrote `CL`. Spelled
   back, each is rc1's file byte for byte (macOS arm64).
3. **`nucleotide_chain_1BNA_A`** gains DC 1's 16 atoms (243 points; 227 in
   rc1). Its keys were derived without the changed code: the previous build
   (1.0.0rc6) computed the same chain with its count rule switched off, and the
   current keys file equals that build's. In the test, the triangles without a
   DC 1 atom reproduce the rc1 keys file byte for byte, and the added ones
   number exactly C(243, 3) − C(227, 3) = 438,216. This is regression evidence
   for the selection change, not independent evidence for the keys.

**Counting.** The test file holds 48 test items: the 42 calculations, each a
parameter of `test_the_standard_result_is_unchanged`, and six checks
(`test_the_rc1_reference_file_is_unchanged`,
`test_every_calpha_case_has_a_base30_expectation_and_nothing_else_does`,
`test_every_atomic_case_has_a_current_expectation_that_differs_from_rc1_only_as_approved`,
`test_every_method_is_covered`, `test_every_supported_modifier_is_covered_somewhere`,
`test_previous_selection_membership_is_decided_by_residue_identity`).
A calculation with more than 250,000 triangles is marked slow:

| run | calculations | checks | test items |
|---|---|---|---|
| default (`pytest`) | 34 | 6 | **40** |
| slow (`pytest -m slow`, part of release validation) | 8 | 0 | **8** |
| total | **42** | 6 | 48 |

So "40 + 8" counts test items and "42" counts calculations; no calculation runs
in both sets, and every one runs in exactly one.

**What is compared.** For every one-structure calculation:

* the SHA-256 of the native keys file — every key, its category (regular,
  mirror or SSE category) and its frequency, byte for byte;
* the number of points, of triangles and of distinct keys;
* both identity fingerprints — the key space and the sampling profile;
* the residues used, in order, and the per-category totals;
* on macOS arm64 only, and only for calculations of at most 120,000 triangles,
  the SHA-256 of the native triplet file (raw floating-point geometry, which is
  platform-specific in its last digits; see [PLATFORM_VALIDATION.md](PLATFORM_VALIDATION.md)).
  32 calculations carry this check.

For the 17 protein Cα calculations the values compared are those of
`calpha_base30.json`, and the mapping of the old-selection triangles back to
the rc1 file, and the count of the added ones (above), are checked as well.
For the 25 atomic calculations they are those of `atomic_element_labels.json`,
with the checks against rc1 above.

The one per-residue sweep (`amino_acid_atomic_1UBQ_A_sweep`) compares the list
of residues, each residue's keys-file SHA-256, both fingerprints and the total
triangle count.

| # | case id | method | structure, chain | modifiers / target | triangles | run in | compared |
|---|---|---|---|---|---|---|---|
| 1 | `protein_1UBQ_A` | `protein` | 1UBQ A | — | 70,300 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 2 | `protein_1UBQ_A_grouping` | `protein` | 1UBQ A | aa_grouping=True | 70,300 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 3 | `protein_1UBQ_A_mirror` | `protein` | 1UBQ A | mirror=True | 70,300 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 4 | `protein_1UBQ_A_size25` | `protein` | 1UBQ A | size_filter=25.0 | 60,742 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 5 | `protein_1UBQ_A_all_modifiers` | `protein` | 1UBQ A | aa_grouping=True, mirror=True, size_filter=25.0 | 60,742 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 6 | `protein_1UBQ_A_sse` | `protein` | 1UBQ A | sse=True | 70,300 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 7 | `protein_1UBQ_A_sse_grouping_mirror` | `protein` | 1UBQ A | aa_grouping=True, mirror=True, sse=True | 70,300 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 8 | `protein_1BRS_D` | `protein` | 1BRS D | — | 105,995 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 9 | `aa_selection_1UBQ_A` | `aa_selection` | 1UBQ A | residues=10-12;48;70-71 | 20 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 10 | `aa_selection_1UBQ_A_mirror` | `aa_selection` | 1UBQ A | mirror=True, residues=1-30 | 4,060 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 11 | `aa_selection_1UBQ_A_sse` | `aa_selection` | 1UBQ A | residues=1-30, sse=True | 4,060 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 12 | `amino_acid_atomic_1UBQ_A_LEU15` | `amino_acid_atomic` | 1UBQ A | residue=15, residue_name=LEU | 56 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 13 | `amino_acid_atomic_1UBQ_A_LEU15_mirror_size` | `amino_acid_atomic` | 1UBQ A | mirror=True, residue=15, residue_name=LEU, size_filter=4.0 | 28 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 14 | `amino_acid_atomic_1UBQ_A_sweep` | `amino_acid_atomic` | 1UBQ A | — | 4,860 | default | per-residue keys SHA-256 (75 residues), residue list, both fingerprints, total triangles |
| 15 | `amino_acid_atom_selection_1UBQ_A` | `amino_acid_atom_selection` | 1UBQ A | atoms=LEU_8_CD1;LEU_43_CD2;ILE_44_CD1;LYS_48_NZ;HIS_68_NE2;VAL_70_CG1 | 20 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 16 | `aa_ind_1UBQ_A_LEU` | `aa_ind` | 1UBQ A | residue_types=LEU | 59,640 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 17 | `aa_ind_mix_1UBQ_A_LYS_ARG` | `aa_ind_mix` | 1UBQ A | residue_types=LYS;ARG | 198,485 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 18 | `aa_ind_specific_1UBQ_A` | `aa_ind_specific` | 1UBQ A | residues=LYS_6;LYS_11;LEU_15 | 2,600 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 19 | `drug_3PTB_A_BEN1` | `drug` | 3PTB A | entity_id=1, entity_name=BEN | 84 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 20 | `drug_6G9M_A_ESW404` | `drug` | 6G9M A | entity_id=404, entity_name=ESW | 7,770 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 21 | `drug_6G9M_A_ESW404_mirror` | `drug` | 6G9M A | entity_id=404, entity_name=ESW, mirror=True | 7,770 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 22 | `drug_atom_selection_3PTB_A_BEN1` | `drug_atom_selection` | 3PTB A | atoms=C1;C2;C3;C4;N1, entity_id=1, entity_name=BEN | 10 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 23 | `carbohydrate_5WT9_A_NAG1` | `carbohydrate` | 5WT9 A | entity_id=1, entity_name=NAG | 364 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 24 | `carbohydrate_chain_5WT9_A` | `carbohydrate_chain` | 5WT9 A | — | 8,436 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 25 | `carbohydrate_atom_selection_7LF9_D_NAG1` | `carbohydrate_atom_selection` | 7LF9 D | atoms=C1;C2;N2;O5;O6, entity_id=1, entity_name=NAG | 10 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 26 | `nucleotide_1BNA_A_DC1` | `nucleotide` | 1BNA A | entity_id=1, entity_name=DC | 0 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 27 | `nucleotide_2AC0_E_DT7` | `nucleotide` | 2AC0 E | entity_id=7, entity_name=DT | 1,140 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 28 | `nucleotide_chain_1BNA_A` | `nucleotide_chain` | 1BNA A | — | 1,923,825 (rc1); 2,362,041 with DC 1 since 1.0.0rc7 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 29 | `nucleotide_atom_selection_2AC0_E_DT7` | `nucleotide_atom_selection` | 2AC0 E | atoms=P;C1';N1;C2;O4;C7;C6, entity_id=7, entity_name=DT | 35 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 30 | `drug_protein_cross_6G9M_A_A` | `drug_protein_cross` | 6G9M A×A | chain_2=A, entity_id=404, entity_name=ESW | 900 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 31 | `drug_protein_cross_3PTB_A_A` | `drug_protein_cross` | 3PTB A×A | chain_2=A, contact_cutoff=4.5, entity_id=1, entity_name=BEN | 8,073 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 32 | `carbohydrate_protein_cross_7LF9_C_A` | `carbohydrate_protein_cross` | 7LF9 C×A | chain_2=A, entity_id=1, entity_name=NAG | 288 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 33 | `nucleotide_protein_cross_2AC0_E_B` | `nucleotide_protein_cross` | 2AC0 E×B | chain_2=B, entity_id=8, entity_name=DG | 126 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 34 | `protein_protein_cross_1BRS_A_D` | `protein_protein_cross` | 1BRS A×D | chain_2=D | 1,230,714 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 35 | `protein_protein_cross_1BRS_A_D_size20` | `protein_protein_cross` | 1BRS A×D | chain_2=D, size_filter=20.0 | 1,016,887 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 36 | `protein_protein_cross_1BRS_A_D_filtered` | `protein_protein_cross` | 1BRS A×D | chain_2=D, contact_filter=side_a_only, residue_types=LYS;ARG;GLU;ASP | 10,936,233 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 37 | `protein_protein_calpha_cross_1BRS_A_D` | `protein_protein_calpha_cross` | 1BRS A×D | chain_2=D | 8,151 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 38 | `protein_protein_calpha_cross_1BRS_A_D_cut8` | `protein_protein_calpha_cross` | 1BRS A×D | chain_2=D, contact_cutoff=8.0 | 49,140 | default | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals; triplet SHA-256 (macOS arm64) |
| 39 | `protein_2AC0_A_acceptance` | `protein` | 2AC0 A | — | 1,274,196 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 40 | `protein_3EJ1_A_acceptance` | `protein` | 3EJ1 A | — | 4,235,315 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 41 | `protein_6G9M_A_acceptance` | `protein` | 6G9M A | — | 6,550,610 (rc1); 6,608,580 with CME 161 since 1.0.0rc5 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |
| 42 | `protein_7LF9_A_acceptance` | `protein` | 7LF9 A | — | 1,456,935 | `-m slow` | keys SHA-256, points/triangles/keys, both fingerprints, residues used, category totals |

Case 26 (`nucleotide_1BNA_A_DC1`) is a calculation with no triangles on
purpose: the 5′-terminal nucleotide has no phosphate, so the method's
completeness rule excludes it, and the reference pins that outcome.

## 2. Independent evidence

| evidence | compared with | what it shows |
|---|---|---|
| `tests/characterization/test_1ubq_historical_difference.py` | the key of each of the 70,300 triangles of 1UBQ chain A written by the earlier program `TSR-Package/tsr_package/TSR.py` (commit `fd3c425a`), extracted into `tests/data/reference/historical/` | written at the earlier program's radix 20 (BioTSR encodes the residue labels at radix 30, ADR-107; 1UBQ has no modified residue, so the labels are the same), exactly 1,144 keys differ, all in the angle class, and every one follows from the documented angle-divisor change (`3.14` → `math.pi`, ADR-012) alone; no other difference |
| `tests/unit/test_geometry.py`, `test_binning.py`, `test_encoding.py`, `test_atomic.py` | geometry, bin boundaries and key arithmetic worked out by hand | the triangle geometry, folding, classes and key formula on cases whose answers are known independently of the code |
| `tests/unit/test_analysis_math.py` | Jaccard values, distances and cluster partitions worked out by hand | the analysis measures, their conversions and the clustering on hand-calculable inputs |
| `validation/scripts/analysis_scipy_equivalence.py` | SciPy 1.11.1 and 1.18.1, scikit-learn 1.6.1 and 1.9.1 (run by hand; records in `validation/evidence/analysis/`) | the scope stated in [methods/ANALYSIS.md](methods/ANALYSIS.md), §5.7 |
| the counts in the reference cases | combinatorics: C(*n*, 3) triangles for *n* points with no filter | e.g. 1UBQ chain A, 76 C-alpha atoms, 70,300 triangles |

None of this is a validation of TSR as a scientific method; that is the subject
of the TSR group's publications. It shows that BioTSR computes what its
documentation says it computes.
