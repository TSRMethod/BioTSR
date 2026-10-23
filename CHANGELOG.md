# Changelog

All notable user-visible changes to BioTSR. Versions follow
[PEP 440](https://peps.python.org/pep-0440/); the release process is in
[docs/RELEASING.md](https://github.com/TSRMethod/BioTSR/blob/main/docs/RELEASING.md).

## [1.0.0rc7] — release candidate (unreleased)

The seventh release candidate. It changes how atoms are labelled in every
atomic method and how `nucleotide_chain` admits nucleotides; everything else is
as described under 1.0.0rc6 below. **Every atomic result must be recomputed**:
atomic results written by 1.0.0rc6 or earlier are in another key space, so
`--if-exists skip` does not reuse them and no analysis compares them with
rc7 results, whatever opt-in is given. Cα results are unaffected: their key
spaces and processing rules are unchanged, and `--if-exists skip` keeps them.

### Changed — atomic labels come from the element field (ADR-113)

- **Every atomic method labels a heavy atom by its element field** (PDB
  columns 77–78, any letter case), no longer by its atom name. This covers
  `amino_acid_atomic`, `amino_acid_atom_selection`, `aa_ind`, `aa_ind_mix`,
  `aa_ind_specific`, the drug, carbohydrate and nucleotide methods, and the
  atomic cross methods. Atom names still identify atoms in selections, records
  and outputs.
- The labels are those the historical atom lexicon gave each element —
  carbon 106, nitrogen 108, oxygen 109, phosphorus 110, sulfur 111, chlorine
  93, calcium 99, holmium 37 … — so **no label was renumbered**; the radix stays
  112. What changes is the label of an atom whose name used to choose another
  element: ATP `PA`/`PB` (were 55 and 56) and every other phosphorus are 110;
  names such as `N1A`, `CB1`, `CB2` (were 107, the hydrogen label) take their
  element's label; a chlorine is 93 whatever it is called (`CL2` of 4X6J was
  106; `CL35` of 3FSF could not be labelled and refused its calculation); a
  carbon named `CL1` is 106 (was 93); calcium `CA1` is 99 (was 106).
- **One packaged table**, `labels/resources/element_labels_v1.csv`
  (`element,label`; 110 elements, 109 labels; hydrogen and deuterium share
  107 and are excluded before labelling, as before). It replaces
  `atom_lexicon_standard_v1.csv`. The table is validated when it is loaded —
  a malformed row, a symbol not in chemical case, a duplicate symbol or an
  unexpected shared label is refused.
- An atom whose element field is blank or names an element the table does not
  list (`X`) still refuses a calculation that selects it (ADR-108). Nothing is
  inferred from the atom name.
- **Key spaces.** The scheme is `element_labels_112`; its fingerprint records
  the table and the rule (`label_source = element`, `resolver =
  element_field_v1`). Every atomic key space has a new fingerprint (`1ecbef1b…`
  for the plain atomic methods; `1187735d…` before). The 25 atomic standard
  reference results keep their keys and counts (`nucleotide_chain` apart);
  triplet files spell a label name as its element symbol (`Cl` where rc6
  wrote `CL`).
- **API.** `biotsr.ELEMENT_LABELS_112` replaces `biotsr.ATOM_LEXICON_STANDARD_112`;
  the module `biotsr.labels.lexicon` (with `audit_lexicon` and the namespace
  loader) is removed, and `biotsr.labels.elements` loads the table.
  `LabelDomain.ATOM_ELEMENT` (`"atom_element"`) replaces
  `LabelDomain.ATOM_LEXICON`.

### Changed — `nucleotide_chain` takes every nucleotide with the atoms it has (ADR-114)

- `nucleotide_chain` no longer requires each nucleotide's exact atom count:
  each standard nucleotide takes part with every heavy atom it has, at a chain
  end or inside the chain. 1BNA chain A: 243 points (227 before; the 16 added
  are DC 1's, which lacks its 5′ phosphate).
- Every nucleotide whose atoms differ from its standard set is named in one
  warning — `atom sets that differ from the standard nucleotide (as the
  structure file has them; every atom present is used): DC_1 lacks P, OP1,
  OP2` — shown by `biotsr run` and in the batch summary. It does not make the
  result partial: the status stays `ok` and the exit code 0.
- `nucleotide` and `nucleotide_atom_selection` keep the exact count; a
  nucleotide type outside DA, DC, DG, DT, A, C, G, U is still not read.
- The completeness rule (`available_nucleotide_atoms`) is part of the
  processing identity of every method that has one, so rc6 `nucleotide_chain`
  outputs are not reused.

### Changed — resume messages

- When an output is in another key space and was computed under other
  processing rules, `--if-exists skip` names both reasons, key space first.

## [1.0.0rc6] — release candidate (unreleased)

The sixth release candidate. It changes how structure records are read and
which selections are refused; keys, label tables and every other rule are as
described under 1.0.0rc5 below. Results written by 1.0.0rc5 or earlier record
no `input_rules`, so they are not reused by `--if-exists skip` and are not
mixed silently with rc6 results in an analysis (ADR-108).

### Changed — reading atom records (ADR-108)

- **Atom names are the full PDB field, columns 13–16,** in every atomic method
  and every atomic cross side. A name starting in column 13 — a two-letter
  element such as an `ATOM`-record zinc `ZN`, or a ligand `CL1` — is read as
  written instead of losing its first letter. Standard protein and nucleic-acid
  atoms are unaffected. The Cα test is unchanged: only an atom named exactly
  `CA` in columns 14–15 is a Cα (calcium, written from column 13, never is,
  nor is `CA1`).
- **The element is read from columns 77–78 only, without regard to letter
  case.** `amino_acid_atomic` used to read columns 78–80, which include the
  charge, so hydrogens written `H0`/`H1+` counted as heavy atoms and failed the
  residue's exact atom count; they are now excluded like any hydrogen. `h` is
  hydrogen.
- **Refused, naming the records, instead of computed:**
  - a selected atom with a blank element field (atomic calculations and the
    contact atoms of every cross; a `HETATM` Cα may still have a blank element,
    ADR-107);
  - two selected records of one atom (repeated or conflicting), and two
    different selected atoms at exactly the same coordinates — nothing is
    merged, preferred or repaired, and there is no distance tolerance;
    `aa_selection` and `amino_acid_atom_selection` no longer keep the first of
    two matching records;
  - a record the selection would include whose residue number cannot be read
    (hybrid-36 numbers are not decoded); such records used to be skipped with a
    warning and the result reported complete;
  - a selected atom the atom table cannot label, now as an `InvalidRequest`
    (exit 2) naming the atom, its residue and the remedy, instead of a
    `KeyError`. An atom outside the selection — an unrequested atom of a ligand
    narrowed with `atoms=` — is no longer labelled at all.

### Changed — selections (ADR-109)

- **UNK:** no calculation is done for a selection that includes a residue named
  `UNK` — a whole chain that contains one, a named UNK residue, residue type or
  atom, or a cross side or Cα-cross contact selection that contains one. It is
  refused with the residues named; the rest is not computed without it. An
  unselected chain containing UNK does not affect other calculations.
- **Protein–protein cross needs two different chains** (`protein_protein_cross`,
  `protein_protein_calpha_cross`). A molecule side (drug, carbohydrate,
  nucleotide) may share a chain identifier with the protein.
- **The two sides of any cross may not share an atom.** An overlap is refused
  and named; nothing is subtracted from either selection.
- **The narrowed Cα cross narrows its interface.** With `calpha_residues` or
  `calpha_residues_2`, only the heavy atoms of the selected residues decide
  interface residues, on both sides; a partner residue that touches only
  unselected residues is no longer part of it. Whole-chain runs are unchanged.
- **The protein side of the carbohydrate- and nucleotide-protein cross stops
  at the protein chain's `TER`,** like every other protein side. A residue
  written after it (a free phosphotyrosine, say) no longer joins the protein;
  a glycan or nucleotide after it is still selectable as the molecule side.
- **Residues outside a side's vocabulary are named.** The whole-chain protein
  side of the carbohydrate- and nucleotide-protein cross reads the 20 regular
  amino acids and TPO, SEP, PTR; any other residue of the chain is now listed
  as `not supported:` (the result is partial), instead of being left out
  silently.

### Changed — secondary structure and insertion codes (ADR-110)

- A `HELIX`/`SHEET` record's endpoints are matched by chain, residue number
  and insertion code, and the record covers the residues between them in the
  chain's own order. Residues with insertion codes inside a range are now
  annotated (they were always coil), and a record no longer covers a residue
  outside its endpoints. When an endpoint cannot be resolved — absent from the
  chain, written in two places, unreadable, on another chain, or after the end
  — SSE annotation of that chain is refused, naming the record. Keys, points
  and fingerprints are unchanged; `sse.version` is `sse_assignment_v3`.

### Changed — chains and batches (ADR-112)

- **A batch request without its chain is skipped, not failed and not
  guessed.** A manifest or sheet row that leaves a required chain (or, for a
  cross method, the second chain) blank, with no `--chain` default, gets status
  `chain_not_specified` in `batch_summary.csv`; the other rows run, and the
  batch exits 4. Before, a manifest with such a row was refused as a whole and
  a sheet row failed.

### Not changed

- Terminal residues (ADR-111): Cα methods keep terminal Cα points; per-residue
  and residue-type atomic methods keep the exact-atom-count rule, so a terminal
  residue with `OXT` (or a 5′ nucleotide without its phosphate) is excluded and
  named; cross methods apply no completeness rule. `nucleotide_chain` keeps
  its per-nucleotide atom-count rule for now — see the known limitations.
- Bare or chain-less `TER` records still end nothing (documented limitation).
- The packaged atom table is unchanged.

## [1.0.0rc5] — release candidate (unreleased)

The fifth release candidate. It changes the protein Cα label space; everything
else is as described under 1.0.0rc4 below.

### Changed — protein Cα residue labels, at radix 30 (ADR-107)

- **The Cα methods label modified residues.** `protein`, `aa_selection` and
  `protein_protein_calpha_cross` label residues from the TSR Research Group's
  residue tables: the 20 regular amino acids keep codes 4–23, and eleven
  modified residues are added — phosphothreonine `TPO` 24, phosphoserine `SEP`
  25, phosphotyrosine `PTR` 26, 3-nitrotyrosine `NIY` 27, trimethyllysine `M3L`
  28, formyl-tryptophan `TRF` 29, S-hydroxycysteine `CSO` 30, lysine-pyridoxal
  phosphate `LLP` 31, `CME` 32 and 3-sulfinoalanine `CSD` 33, with
  selenomethionine `MSE` sharing methionine's 13. The method names are
  unchanged.
- **Radix 30.** The 30 codes are encoded at radix 30 (`amino_acid_30`). At the
  old radix 20 codes above 23 carried into the neighbouring label digit, and
  different residue combinations (for example SEP–PHE–PHE and TPO–TPO–PHE)
  could share a key. All 4,960 possible label triples now have keys of their
  own.
- **Grouping** (`--aa-grouping`) uses the grouped table (`amino_acid_group_23`):
  the 13 historical groups of the regular amino acids, MSE with MET, and each
  other modified residue in a group of its own — 23 codes spanning 4–33, also
  at radix 30. Grouped and ungrouped keys share the radix but remain different
  key spaces. All 2,300 grouped label triples have keys of their own.
- **Modified residues written as `HETATM` take part.** A Cα method reads a
  `HETATM` Cα when the residue is one of the eleven and the CA atom's element
  field is carbon or blank — no other atom of the residue is required; an
  `ATOM` record of them is read like any other (until now it stopped the
  calculation). Ligands, ions with an atom named `CA` and modified residues
  outside the list are never points; the last are still named as
  `not represented:`. A supported residue present in the chain without a
  readable Cα (no CA atom, a CA only as altLoc B or later, or a `HETATM` CA
  declaring another element) is named with that reason, and a residue that
  `aa_selection` names and the chain contains is never reported as "not
  found". The chain, residue-instance, truncation
  and blank-or-`A` alternate-location rules are unchanged. The Cα cross
  decides interface residues with the same residues' heavy atoms.
- **Collision-free encoding is enforced.** A label scheme is refused unless all
  its codes fit in one window of `n_labels` consecutive integers ≥ 1; the old
  rule (number of codes equals the radix) and its override
  (`allow_radix_mismatch`) are removed. The encoder refuses a label outside the
  scheme's codes. `biotsr.labels.label_triple_collisions` checks a table
  exhaustively.
- **Provenance and resume.** The Cα methods record
  `modified_residue_policy = hetatm_amino_acids:CME,CSD,CSO,LLP,M3L,MSE,NIY,PTR,SEP,TPO,TRF;ca_element:C_or_blank`;
  the atomic methods keep `hetatm_amino_acids:PTR,SEP,TPO`.
  `biotsr.api.processing_identity` now takes the method name.
  `--if-exists skip` also requires each result's key space to be the one this
  version computes for the request.
- **Public names.** `AMINO_ACID_20` and `AMINO_ACID_GROUP_13` are replaced by
  `AMINO_ACID_30` and `AMINO_ACID_GROUP_23`; `biotsr.labels` also exports
  `STANDARD_AMINO_ACIDS`, `MODIFIED_RESIDUE_CODES`, `ordered_label_triples` and
  `label_triple_collisions`.

### Fixed

- **The Cα cross no longer stops on an atom the atom lexicon cannot label.**
  `protein_protein_calpha_cross` reads interface heavy atoms only for contact;
  an atom whose element is outside the lexicon used to raise `KeyError`
  although its label is never used. It now decides contact like any other
  atom.
- **`analyze()` refuses a cross method by name.** Passing a cross method to
  `biotsr.analyze` fell through to the Cα selector and failed with an
  unrelated `KeyError` on the first residue; it now raises `InvalidRequest`
  naming `analyze_cross`. The command line was not affected.
- **Clearer refusals.** An atomic calculation with an atom the lexicon cannot
  label, and a Cα calculation with a residue name the labels do not cover,
  still stop, now with a message naming the atom or residue and the remedy
  (the old message pointed to an option the API does not offer).
- **Documented readings of unusual input** (KNOWN_LIMITATIONS §2): non-integer
  residue numbers, exact-case altLoc and element fields, duplicate atom names,
  insertion codes in SSE ranges, `nucleotide_chain`'s 5′ nucleotide, and the
  cross-method protein sides; ADR-082 now says where the narrow-slice metal
  defect is still reachable.

### Regenerate older Cα results

Every protein Cα key changes value: a triangle of regular residues keeps its
labels and classes, but its key is written at radix 30 instead of 20 (13 under
grouping). A chain with a supported modified residue also gains that residue's
point and every triangle containing it — on 6G9M chain A, CME 161 adds 57,970
triangles. Cα results computed by earlier versions load, but their key space
differs: no analysis combines them with current results, and
`--if-exists skip` refuses to reuse them and leaves them untouched. Compute them
again with `--if-exists overwrite`, or into a new `--output`. Atomic results —
including the atomic cross methods — are unchanged and remain reusable.

### Unchanged

The atom lexicon and its fingerprint, every atomic radix, and every angle and
distance bin. The frozen reference file `tests/data/reference/standard_baseline.json`
is unchanged; the Cα cases' current values are in
`tests/data/reference/calpha_base30.json`, derived from it by an explicit
re-encoding.

## [1.0.0rc4] — release candidate (unreleased)

The fourth release candidate of BioTSR 1.0: the unified software
implementation of the Triangular Spatial Relationship (TSR) methodology,
replacing a family of independent historical TSR programs with one package.
BioTSR 1.0 is a **standard-only** package: every method runs under one
standard profile, the accepted TSR science as set by the TSR Research Group,
and writes one self-describing output format. 1.0.0 has not been published.

### What BioTSR 1.0 provides

- **One TSR Engine.** Every method runs on the same geometry, canonical vertex
  ordering, angle and distance binning, labelling and integer key encoding.
- **Twenty methods in three families** — protein Cα, atomic (amino-acid, ligand,
  carbohydrate, nucleotide/DNA/RNA, and residue/atom selections), and five cross
  (two-sided) methods: drug–protein, carbohydrate–protein,
  nucleotide–protein, protein–protein all-atom and protein–protein Cα. See the
  [capability matrix](https://github.com/TSRMethod/BioTSR/blob/main/docs/CAPABILITY_MATRIX.md).
- **One standard profile per method.** `biotsr methods NAME` names it and every
  result records it. Where the standard differs from a historical program (for
  example `math.pi` where some historical programs divided by the literal
  3.14), the method records and architecture decisions say so.
- **Composable modifiers** on one call: amino-acid grouping (Cα), mirror
  signing, and size filtering; no per-combination method names.
- **SSE-TSR.** Secondary-structure annotation of Cα triangles into 18 categories
  from the structure's own `HELIX`/`SHEET` records, distinguishing individual
  β-strands. An annotation: it changes no key.
- **Two identities on every result.** A `KeySpace` fingerprint (what the integer
  keys mean) and a `SamplingProfile` fingerprint (how the triangles were drawn).
  Results whose keys are defined differently are refused, never silently
  combined; processing differences are compared only when acknowledged; what
  was selected never prevents a comparison.
- **Analysis.** Feature matrices, four common-key operations, five similarity
  measures, the average-linkage tree with its dendrogram and ordered heat map
  (no cluster count needed), and optional flat clusters by an explicit rule —
  over results from any method family, on the standard library alone (the
  figure itself needs the optional `plot` extra, matplotlib). An executed
  tutorial goes from saved results to exported matrices, the tree, the figure
  and its provenance.
- **Native serialization and deserialization.** A self-describing key file,
  triplet file and metadata sidecar sharing one prefix (`PREFIX.keys`,
  `PREFIX.triplets`, `PREFIX.meta.json`); `biotsr.load_result` reads a result
  back with both fingerprints restored, so analyses run over saved results
  without recomputation.
- **Command line.** `biotsr run`, `biotsr batch` (dataset manifests and wide
  dataset sheets, methods chosen with `--method`),
  `biotsr fetch`, `biotsr analyze common-keys`, `biotsr methods`,
  `biotsr describe`, and `python -m biotsr` for clusters.
- **Datasets.** One-request-per-row manifests and wide dataset sheets, whose
  methods are chosen when BioTSR is invoked and whose columns are resolved per
  method (`--sheet`, `--map`); manifest columns with other names are mapped
  (`--column FIELD=HEADER`), kept as metadata (`--keep-column`) or dropped
  (`--ignore-column`); no chain is ever guessed; `--dry-run` previews every
  request; every calculation writes into one shared output folder under a
  descriptive prefix with a completion record written last; `--if-exists
  skip` resumes only over calculations whose record shows exactly the same
  request with complete files; two requests that would write the same prefix
  are refused before anything runs; `batch_summary.csv` lists every request
  with its source rows and status. `batch.run` (`from biotsr import batch`)
  does the same from Python.
- **Parallel batches.** `--jobs N` or `--jobs auto` computes requests in worker
  processes through joblib, with results identical to one process; triplet
  files are written as they are computed.
- **Structure retrieval, opt-in.** `biotsr fetch`, `biotsr batch
  --download-missing` and `biotsr.fetch_structures` download PDB-format entries
  from RCSB over HTTPS, verify every file before installing it, reuse valid
  cached files and report every id.
- **One runtime dependency.** Python 3.9–3.13 on Linux and macOS. joblib,
  installed with BioTSR, runs batch rows in parallel; the calculation uses the
  standard library only. `biotsr[parallel]` remains as an alias that adds
  nothing.
- **Licensed under the Apache License 2.0.** `NOTICE` carries the notices of
  third-party material, which keeps its own terms: `TSR-Package` and
  `CarbohydratesTSR` (MIT) and SciPy (BSD 3-Clause).

### Changed since 1.0.0rc3

1. **Comparison rule: key definitions always, processing on acknowledgement,
   selections never** (ADR-094). Every analysis requires one key space, decided
   from the key-space definition rather than method names, so different
   methods sharing a key space compare and the same method with another key
   convention does not. Processing settings (size filter, contact cutoff and
   rule, completeness rule, candidate strategy, pool rules, alternate-location
   rule) are checked separately and need `allow_mixed_processing=True` /
   `--allow-mixed-processing` when they differ. What was selected — structures,
   chains, amino-acid types (LYS-only against TYR-only), ligands,
   carbohydrates, residue numbers, atom subsets, copies, the method — never
   blocks a comparison and is recorded per sample; an empty intersection is a
   result. `allow_mixed_sampling` / `--allow-mixed-sampling` are renamed (the
   old flag is refused with the new name); reports carry a `comparison`
   summary and per-sample processing fingerprints and settings. Neither
   identity fingerprint changed.
2. **One alternate-location rule** (ADR-095). Every method keeps blank and `A`
   and excludes `B`, `C` and later, with no fallback to another conformer.
   *Scientific change for `amino_acid_atomic`*, which used to keep `C` and
   later conformers: a residue with three conformers is now computed from
   blank + `A` instead of being excluded as incomplete, and a residue whose
   side chain exists only as `B`/`C` is now excluded instead of being built
   from `C`. Excluded atoms are counted in a result warning. No frozen
   reference result changed (no test structure has an altLoc beyond `B`).
3. **Molecule instances** (ADR-096). A molecule named with its identifier is
   that copy; named without one it is **every copy of that name on the chain,
   in one combined calculation** (`NAME_allcopies`; triangles may span copies),
   where rc3 refused several copies. Identifiers are PDB residue numbers with
   an optional insertion code, never positions; malformed and ambiguous ones
   are refused. Atom lists apply to each copy, with missing atoms reported per
   copy. `amino_acid_atomic --residue-name NAME` alone runs each residue of that
   name. Malformed residue and atom tokens (`aa_ind_specific`, atom selections)
   are now refused instead of skipped with a warning. The selection record
   states each side's scope and selected instances.
4. **Methods chosen outside the CSV; wide sheets read method-aware** (ADR-097).
   `--method` is repeatable and expands each row into one request per method;
   `--sheet tsr|custom`, `--map`, `--map-file`, `--participant` and
   `--inapplicable` read a wide sheet whose rows carry several participant
   chains and alternative selections, each method reading its own columns.
   Unused columns are reported, identical requests from several rows are run
   once, and `--dry-run` previews every request with its source rows, chains,
   selection, output prefix and outcome.
5. **One shared output folder** (ADR-098). Every calculation writes
   `PREFIX.keys`, `PREFIX.meta.json`, `PREFIX.triplets` and its completion
   record `PREFIX.run.json` into `--output`, e.g.
   `7LF9_D_carbohydrate_NAG_1.keys`; the `_29_35` suffix is gone (the bin
   definitions are in the metadata); file contents are unchanged. Staging,
   locking, collision checks (including letter case), interruption, overwrite
   and resume are rebuilt for the shared folder; files BioTSR did not write are
   never replaced. The workers of a parallel batch share one staging
   directory that only the batch removes (a worker removing it while empty
   raced another worker's staging; macOS reported `OSError: [Errno 22]`).
   Result schema 2, record schema 4: outputs of rc3 and earlier
   are not read and must be recomputed. `biotsr analyze common-keys
   --from-results DIR --method NAME` analyses one method's calculations of a
   shared folder that holds several key spaces.
6. **Hierarchical analysis without a cluster count** (ADR-099). The tree, the
   dendrogram and the ordered heat map need no `n_clusters`;
   `hierarchy()`, `cut_hierarchy(n_clusters=|height=)` and an optional flat
   cut in `clustermap_spec` replace rc3's required count. Observed classes and
   algorithm clusters are separate, labelled stripes; no ARI without a
   partition.
7. **Redistribution confirmed** (ADR-100). The project owner has confirmed
   redistribution of the historical tables; attribution and notices are
   unchanged.
8. **Modified amino acids** (ADR-101). Phosphothreonine `TPO`, phosphoserine
   `SEP` and phosphotyrosine `PTR` -- the modified residues with a defined
   atom count (11, 10, 16) and lexicon labels -- are read from `HETATM`
   records by the amino-acid atomic methods and the protein sides of the
   atomic cross methods; until now RCSB files' `HETATM` copies were never
   seen. **Scientific change:** results that include those residues change
   (for example `aa_ind_mix TYR;TPO` on a phosphorylated chain now includes
   TPO); the 42 frozen references contain none and are unchanged. Any other
   modified residue is named as `not supported`, and the Calpha methods name
   every one as `not represented`, in a warning and the selection record
   (`unsupported`, `not_found`); a request of only such residues is refused.
   A batch calculation that lacks anything it named is `partial`, counted
   separately; `biotsr run` prints those warnings. A named residue the
   completeness rule declines is reported as excluded, no longer as not found.
9. **One residue grammar** (ADR-102). Identifiers are PDB numbers with an
   optional insertion code everywhere; `aa_selection` ranges are plain numbers,
   ascending. Insertion-code ranges (`27A-29`), descending ranges and ranges
   spanning an insertion-code residue that is not listed are refused with the
   remedy -- previously kept verbatim (matching nothing) or expanded to nothing.
10. **Duplicate rows** (ADR-103). Identical requests merge only for the same
    sample id; distinct sample ids stay distinct samples; rows that would merge
    with different labels are refused; differing kept metadata is recorded per
    row; the summary says when one sample serves several rows.
11. **Cross sides from a wide sheet** (ADR-104). A side the row could narrow
    in several ways must be chosen with `--cross-side METHOD.SIDE=COLUMN|whole`
    (`cross_sides=` in Python), or the run is refused; unused populated
    columns are named per request, and the preview and summary show both
    sides' resolved selections. Previously such columns were silently unused.
12. **Processing rules are part of a calculation's identity** (ADR-105). Every
    completion record now holds the processing rules it was computed under
    (`request.processing`: the alternate-location and modified-residue rules),
    and every result records the modified-residue rule
    (`source.modified_residue_policy`). `--if-exists skip` keeps an existing
    output only when they match; an output written by an earlier 1.0.0rc4 build,
    which may have excluded `HETATM` TPO, SEP and PTR for the identical request,
    is refused and recomputed with `--if-exists overwrite`. Analyses compare the
    rule with the other processing settings; an unrecorded rule reads as
    *not recorded* and needs `--allow-mixed-processing`. Key spaces are
    unchanged; analysis processing fingerprints change.
13. **Partial results are explicit everywhere** (ADR-106). What a selection
    asked for and did not get is recorded structurally (`not_found`,
    `unsupported`, `declined`, missing atoms; `omissions` in each completion
    record). New exit code **4**: everything finished and nothing failed, but
    some result is partial -- returned by `run`, `batch` and `analyze`
    (previously 0). A resumed batch that keeps a partial output says so. Every
    analysis refuses partial samples, naming each and what it lacks, unless
    `--allow-partial` (`allow_partial=True`) -- an acknowledgement separate
    from `--allow-mixed-processing` that never relaxes the key-space check; the
    report then records each sample's omissions.

### Known limitations

See [docs/KNOWN_LIMITATIONS.md](https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md). In particular:
candidate-strategy extensions (such as CA-kNN) are not part of 1.0; mmCIF input
is not supported, and mmCIF-only entries cannot be fetched; Windows is not
supported; outputs written by 1.0.0rc3 or earlier are not read, and outputs
written by earlier 1.0.0rc4 builds cannot be resumed (they are recomputed); and modified
residues other than TPO, SEP and PTR are not read, while the Calpha methods
represent none (whether they should, and with which label, is undecided).

## [1.0.0rc3] — release candidate (not published)

The third release candidate. What it provided is described under 1.0.0rc4.

### Changed since 1.0.0rc2

1. **One sampling rule for every analysis, compared from structured records.**
   Feature matrices (and so similarity, clustering and the clustered heat map)
   now refuse a mixed sampling rule exactly as the common-key operations do,
   unless `allow_mixed_sampling=True`; the key space is never relaxed. The rule
   is the method, the sampling settings without the chain names, and the
   selection *policy* each result now records (`provenance.selection_request`):
   residue-type sets, entity vocabularies and the kind of narrowing are part of
   it; structures, chains, entity copies, residue numbers and atom tokens are
   targets and are not. rc2 decided cross sides by deleting `name=value` tokens
   from a free-text description — which also deleted a ligand's atom subset, a
   protein side's residue types and a pooled vocabulary — and never saw the
   arguments of one-sided methods, so `aa_ind` over LYS and over TYR compared as
   one rule. Every report records each sample's fingerprints, settings,
   selection rule and targets and its own warnings. `feature_matrix`'s
   `require_same_sampling` and `SamplingProfile.rule_fingerprint`/`rule_dict`
   are removed; the sampling fingerprint itself is unchanged. Results saved by
   rc2 have no selection record and are compared only with the opt-in (except
   `protein`-like methods whose selector takes no argument); the batch request
   schema is now 3, so `--if-exists skip` does not reuse an rc2 output as if it
   were current (it reports the row; `--if-exists overwrite` recomputes it).
2. **Download cancellation is prompt in every phase.** A connection being made
   is polled for cancellation, every connected socket is shut down on
   cancellation, and a response body is read in pieces with a cancellation
   check between them, so Ctrl-C/SIGTERM no longer waits for a stalled or
   slowly trickling response (rc2 waited for the whole body). Nothing is
   installed after a cancellation, and the report stays exact. Only a host-name
   lookup and a TLS handshake step in progress are waited for. `--timeout` is
   documented as what it is — a limit on each wait, not on a whole download;
   no total time limit exists.
3. **The figure needs matplotlib only.** The renderer used `seaborn.clustermap`,
   which needs scipy at run time, which neither seaborn nor the `plot` extra
   installed — so a clean `biotsr[plot]` could not draw — and it imported pandas
   without declaring it, switched the caller's matplotlib backend and closed
   every open figure. It now draws the heat map, both dendrograms, class and
   cluster stripes with legends and a labelled colour bar with matplotlib alone,
   on a standalone `Figure`, from BioTSR's own linkage. `plot = ["matplotlib>=3.5"]`.
4. **`maxclust` follows its definition.** Asking for *t* clusters gives the
   lowest cut that leaves at most *t*: exactly *t* with distinct merge heights,
   fewer (with a warning) only when merges tie. rc2 reproduced scipy 1.11.1's
   search, which could return fewer clusters than an exact cut allows even
   without ties, and two clusters for two samples and *t* = 1; scipy 1.18.1
   behaves as BioTSR now does. The agreement with each scipy release is stated
   precisely (ADR-060). The two illustrative results in the analysis record
   (1BNA ARI 1.0, 1UBQ 0.521) are unchanged.
5. **Explicit scale and exports.** Similarity and distance matrices are fractions
   unless explicitly `scaled(100)`; a scaled matrix records its scale, its
   diagonal and conversions follow it, a value outside the scale is refused, and
   clustering no longer multiplies by 100 behind the caller's back.
   `to_csv()` (square, sample ids on both axes, the measure and kind in the first
   header cell) and `to_long_csv()` export them; `ClusterAssignment.to_csv()`
   exports cluster and dendrogram position per sample. The default CSV headings
   `Drugs` and `Protein Name` are now `sample_id`.
6. **The number of clusters is never inferred.** `clustermap_spec` requires
   `n_clusters`; class labels are optional (unlabelled sets draw a cluster
   stripe only) and a mapping is refused in favour of a sequence in sample order.
7. **Invalid counts are refused by the pairwise measures** (negative,
   non-finite, non-numeric), an explicit zero is an absent key everywhere, the
   empty-sample rule is one rule for every measure, and a matrix over a set with
   an empty sample carries a warning. Cosine is capped at 1.
8. **altLoc provenance is accurate.** `amino_acid_atomic` results record
   `source.altloc_policy = "not_b"`, the rule its standard admission applies
   (ADR-004); every other method records `blank_or_a` as before. The admission
   itself is unchanged.
9. **Clearer refusals of an empty calculation.** A target excluded by the
   completeness rule is refused naming that rule and its reason, and a cross
   whose two sides have no atom in contact is refused naming the contact
   cutoff — both instead of "selection produced no atoms" with a list of
   chains.
10. **Documentation.** An executed analysis tutorial (08); a reference-evidence
    inventory (`docs/REFERENCE_EVIDENCE.md`); measured analysis run times;
    `--jobs` stated as parallelising structure computation, not analysis.

### Known limitations

See [docs/KNOWN_LIMITATIONS.md](https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md). In particular:
candidate-strategy extensions (such as CA-kNN) are not part of 1.0; mmCIF input
is not supported, and mmCIF-only entries cannot be fetched; and output file
names of cross methods are not valid on Windows, which is not supported.

## [1.0.0rc2] — release candidate (not published)

The second release candidate; not published.

### Changed since 1.0.0rc1

**Removed.** 1.0.0rc1 also shipped ways to reproduce historical programs'
outputs. They are not part of 1.0:

- the compatibility profiles of every method, and the cross methods'
  compatibility modes;
- the `--profile` and `--format` options of `biotsr run`, `batch` and
  `analyze`, and the `profile` manifest column (now refused, with a message
  saying to drop it or pass `--ignore-column profile`);
- the compatibility writers in `biotsr.io`, and `biotsr.profile_option` /
  `profile_options`;
- `mode=` on `biotsr.analyze_cross` (now refused with `InvalidRequest`);
- the historical atom lexicons (v1 and v2); `biotsr.ATOM_LEXICON_STANDARD_112`
  is the one packaged lexicon;
- from `biotsr.analysis`: the `LEGACY_STRING_SORTED` column ordering,
  `FeatureMatrix.to_legacy_csv`, `SimilarityMatrix.to_legacy_matrix_csv` and
  `SimilarityMatrix.quantized`.

A request that 1.0.0rc1 accepted under its standard profile and 1.0.0rc2 still
accepts gives the same keys, frequencies and fingerprints. With the
compatibility profiles gone, the historical metal-ion misreading (manganese or
zinc read as nitrogen through a narrow atom-name slice) is not reachable
anywhere: every admission that reads `HETATM` records uses the PDB
specification's atom-name columns, so metal ions are labelled as their
elements (Mn 100, Zn 105).

**Changed behaviour.**

1. **One profile per method.** `profile=` in Python accepts only the method's
   own standard profile or `None`; a profile name, another method's profile or
   a derived profile raises `InvalidRequest`. `run.json` no longer records an
   output format, and its request record is schema 2.
2. **Arguments are checked per method.** An argument the method does not read
   is refused by name, in Python and on the command line, instead of being
   ignored; `biotsr methods NAME` lists the method's target arguments.
   `--chain-2` is required for every cross method. `aa_ind` no longer accepts
   `residue_name` as an alias of `residue_types`.
3. **One molecule per single-entity method.** For `drug`, `carbohydrate`,
   `nucleotide`, their `*_atom_selection` variants and the entity side of the
   three entity cross methods, `--entity-id` may be omitted only when the chain
   holds exactly one copy of the entity; several copies are refused and listed,
   and an id that does not exist is refused. The copy used is recorded
   (provenance `target`, e.g. `BEN_1`).
4. **Empty and partial selections.** An explicit selection that matches nothing
   is refused with a specific message. Requested items that are missing are
   named in a `requested but not found: …` warning, printed by the CLI and
   recorded in `run.json` (`warnings`, per result).
5. **Fewer than three points** gives a `no triangles: …` warning; in a batch the
   row's status is `no_triangles`, which is not a failure.
6. **Common keys require one sampling rule.** `analyze common-keys` and the
   Python common-key functions refuse a set whose sampling rules differ (size
   filter, contact cutoff and rule, selector, completeness…), ignoring which
   structures and chains were sampled, unless `--allow-mixed-sampling` /
   `allow_mixed_sampling=True` is given; the report then carries a warning. The
   JSON report records each sample's key-space, sampling and sampling-rule
   fingerprints and its sampling settings; the CSV gains the fingerprint
   columns. The key-space check is never relaxed.
7. **`analyze common-keys --manifest` writes keys only** (no triplet files, and
   no `--no-triplets` option). With `--from-results`, flags that would compute
   (`--method`, `--chain`, modifiers, `--pdb-dir`, the column options, download
   options, `--jobs`, `--if-exists`, `--on-error`) are refused; hidden
   directories are never read; and a saved result whose files do not match the
   sizes recorded in its `run.json` is refused.
8. **Sample ids are validated.** A sample id must be one directory name
   directly inside `--output`: no path separators, no absolute path, no `.` or
   `..`, no leading dot, no control characters, none of the names the output
   root uses, at most 200 bytes. An invalid id fails its batch row, or exits 2
   for `biotsr run`; `--dry-run` reports it. A symbolic link at a sample's path
   is never written through, replaced or deleted.
9. **Interruption and worker failure.** SIGTERM is handled like Ctrl-C by every
   command (exit 130). A batch interrupted at any point after planning —
   downloading, checking existing outputs or computing — writes its summary
   (and, while downloading, `fetch_report.csv`), and `biotsr fetch` writes its
   report. A worker process that dies (for example, out of memory) no longer
   loses the summary: finished rows are kept, rows that were running are
   `failed` with a message naming the pool failure, rows never started are
   `not_started`, and the batch exits 3. Staging directories left by a process
   killed outright (SIGKILL, power loss) are removed by a later run on the same
   machine once that process is gone; ones from another machine are left and
   reported.
10. **CSV manifests.** A UTF-8 byte-order mark is accepted. Non-blank cells
    beyond the header, or under an empty header, are refused with the line and
    the remedy; duplicate header names are refused.
11. **Download settings** are validated identically by `fetch`, `batch`,
    `analyze` and Python: `--download-workers` 1–16, `--timeout` finite and
    greater than 0, `--retries` a whole number ≥ 0. Invalid values are refused,
    never replaced or clamped. `--timeout` applies per connection or per read,
    not per download. An obsolete entry (`OBSLTE` record) is downloaded with a
    note.
12. **Exit codes.** A batch exits 3 when any row failed, was interrupted by a
    pool failure or was never started. `--on-error stop` exits 3 after a
    computed row fails (2 if a row could not even be planned). `biotsr fetch`
    exits 0 only when every id was downloaded, refreshed or cached; a
    `refresh_failed` id (old file kept) counts as a failure (exit 3).
13. **`--jobs auto`** is capped by the scheduler allocation, CPU affinity,
    cgroup quota and the machine's CPU count, and never starts more workers
    than there are rows.

## [1.0.0rc1] — release candidate (not published)

The first release candidate. It also shipped compatibility profiles for
historical programs, which 1.0.0rc2 removes (see above).

[1.0.0rc2]: https://github.com/TSRMethod/BioTSR
[1.0.0rc1]: https://github.com/TSRMethod/BioTSR
