# BioTSR 1.0 — known limitations

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology. This page states what version 1.0 does **not**
do, or does only under stated conditions, so that nothing here has to be
discovered by accident. Each entry says where it is recorded in more detail.

## 1. The standard only

BioTSR 1.0 computes each method under its one standard profile, the accepted
TSR science as set by the TSR Research Group. It has no option to reproduce the
outputs of the historical TSR programs, their file layouts or their defects.
Where the standard differs from a historical program — `math.pi` where some
programs divided angles by the literal 3.14, a corrected atom-name reading, a
symmetric contact rule — the [methods records](methods/README.md) and the
[architecture decisions](architecture/README.md) say how and why. A result
computed by a historical program is therefore not expected to match a BioTSR
result key for key.

## 2. One PDB reader; PDB format only

BioTSR 1.0 ships exactly one structure reader, `pdb_legacy` (the name is kept
because results record it). mmCIF/PDBx input is not supported — the reader
parses fixed-column PDB records only, so convert mmCIF to PDB first. Since
1.0.0rc6 it reads the atom name from the full PDB field, columns 13–16, and
the element from columns 77–78 (ADR-108); the Cα test alone keeps columns
14–16, so an atom is a Cα only when it is named exactly `CA` there (calcium,
written from column 13, and `CA1` never are). Which reader and which input
rules produced a result are recorded and guarded, so results read under
different rules are not silently mixed (ADR-079, ADR-108).

How a file is read is fixed by each method's standard profile and stated in the
[user guide](USER_GUIDE.md), §12.5:

- **model 1 only** — reading stops at the first `ENDMDL`, or at a `MODEL`
  record numbered above 1; later models are never read;
- **`TER` ends protein and nucleotide reading** of the requested chain — for
  every protein side, including the protein side of
  `carbohydrate_protein_cross` and `nucleotide_protein_cross` since 1.0.0rc6
  (ADR-109); the `HETATM` families (`drug`, the carbohydrate methods) and the
  molecule side of the cross methods read on to the end of model 1, so a ligand
  or glycan written after the protein's `TER` stays selectable. A `TER` record
  is matched to a chain by its chain identifier: a `TER` with no chain
  identifier (a "bare" `TER`, written by some tools) ends no chain, and a `TER`
  naming another chain ends that chain only. BioTSR does not infer which chain a
  bare `TER` was meant to end; check such files before use;
- **alternate locations** — one rule for every method: atoms whose altLoc is
  blank or `A` are kept, and `B`, `C` and later conformers are excluded, with
  no fallback to another conformer when `A` is absent (ADR-095). Each excluded
  atom leaves a skip record (`altloc_rejected`) and the result a warning
  (`alternate locations: kept blank and A; excluded N atom record(s) ...`). A
  residue whose side chain exists only as `B` or `C` is then excluded by the
  method's completeness rule, with the reason recorded. Every result records
  the rule (`source.altloc_policy`: `blank_or_a`), and a different rule is a
  processing difference in analysis. Until 1.0.0rc3 `amino_acid_atomic`
  excluded only `B` and kept `C` and later conformers.
- **another conformer needs a prepared file.** There is no option to choose a
  conformer. To study a non-`A` conformation, prepare a PDB file whose altLoc
  indicators select it — for example keep only the wanted conformer and mark
  it `A` or blank.
- **modified amino acids stored as `HETATM`, atomic methods: three are read,
  the rest are named** (ADR-101). The amino-acid atomic methods (`amino_acid_atomic`,
  `amino_acid_atom_selection`, `aa_ind`, `aa_ind_mix`, `aa_ind_specific`) and
  the protein side of the atomic cross methods (`drug_protein_cross`,
  `carbohydrate_protein_cross`, `nucleotide_protein_cross`,
  `protein_protein_cross`) read exactly phosphothreonine `TPO`, phosphoserine
  `SEP` and phosphotyrosine `PTR`, the modified residues with defined
  heavy-atom counts (11, 10 and 16), labelled by element like every other atom, under the same
  alternate-location and completeness rules; none is converted to its
  unmodified residue. Every other modified residue stored as `HETATM` — the
  modified cysteine `CME` of 6G9M chain A (residue 161), `MSE`, … — is **not
  supported**: it is left out of the calculation, named in a warning
  (`not supported: …`) and in the side's selection record (`unsupported`),
  and a request of only such residues is refused with that reason. The
  `HETATM` entity methods (`drug`) do read such a residue, as a molecule.
- **the Cα methods represent eleven modified residues, and only those.**
  `protein`, `aa_selection` and `protein_protein_calpha_cross` label residues
  from the residue tables, which give TPO, SEP, PTR, NIY, M3L, TRF, CSO, LLP,
  CME and CSD codes of their own (24..33) and MSE methionine's code, at radix
  30 (ADR-107). A `HETATM` Cα of one of them is read when its element field is
  carbon or blank; no other atom of the residue is needed, as for an `ATOM`
  record. Any other modified residue — `HYP`, `KCX`, … — has no residue label
  and is left out and named (`not represented: …`), and so is a supported one
  without a readable Cα, with the reason. A free molecule of a supported
  residue's name written before its chain's `TER` record would also be read as
  a residue: the file does not distinguish it from a residue of the chain,
  and nothing checks the peptide bonds. An `ATOM` record is never
  element-checked, for any residue, so an `ATOM` CA that declares another
  element is read while the same `HETATM` CA is not.

In a batch, a calculation that leaves out something it asked for — not found,
not supported, not represented, or excluded by the completeness rule — has the
status `partial`, written like `ok`, counted separately and naming what is
missing. A partial result makes `biotsr run` and `biotsr batch` exit 4 when
nothing failed, and every analysis refuses it unless `--allow-partial`
(`allow_partial=True`) acknowledges it (ADR-106). A Cα calculation of a chain
that holds a modified residue outside the eleven is therefore partial.

**Cα results written before 1.0.0rc5 must be computed again.** They encoded the
20 regular residues at radix 20 and left out every modified residue, so every
Cα key now has another value, and a chain with a supported modified residue
now has more points and triangles. Their key space differs from a current
result's, so an analysis refuses to combine the two, and `--if-exists skip`
refuses to reuse them and leaves them untouched: recompute them with
`--if-exists overwrite` or into a new `--output`. Atomic results are
unaffected.

**Outputs written by earlier 1.0.0rc4 builds cannot be resumed.** Before
ADR-101 the amino-acid atomic methods, and the protein side of the atomic
cross methods, left out TPO, SEP and PTR written as `HETATM`, under the same
version number and the same key space. Their completion records do not state
the processing rules they were computed under (`request.processing`, added by
ADR-105), so `--if-exists skip` never keeps them as current: it reports
`… exists but cannot be reused: its completion record does not state the
processing rules it was computed under …`. Recompute them with
`--if-exists overwrite` (or into a new `--output`). They can still be analysed
as they are, deliberately: their `modified_residue_policy` is reported as *not
recorded*, so comparing them with current results needs
`--allow-mixed-processing`.

**Results written by 1.0.0rc5 or earlier record no input rules** (ADR-108).
The rules below changed in 1.0.0rc6 without the key space changing, so such a
result is not reused by `--if-exists skip` (it is reported as computed under
different processing rules and left untouched; recompute it with
`--if-exists overwrite` or into a new `--output`), and comparing it with an
rc6 result needs `--allow-mixed-processing`.

What BioTSR **refuses**, naming the records and what to do, rather than
computing a result that would look complete (ADR-108, ADR-109). Each check
applies to the records a calculation would use, never to unrelated records
elsewhere in the file:

- **a residue named `UNK`** in the requested selection — a whole chain that
  contains one (`protein`, a whole cross side, a per-residue sweep of a chain),
  a named UNK residue, residue type or atom, or a cross side or Cα-cross
  contact selection containing one. No calculation is done for UNK, and the
  rest is not computed without it; narrow the selection to leave it out. A
  chain that is not selected may contain UNK;
- **a selected atom with a blank element field** (columns 77–78) in an atomic
  calculation or in the contact atoms of any cross: it cannot be told from a
  hydrogen. BioTSR does not guess an element from the atom name. A `HETATM` Cα
  of a supported modified residue may still have a blank element (ADR-107);
- **two selected records of one atom** (repeated, or conflicting in any field),
  and **two different selected atoms at exactly the same coordinates**.
  Nothing is merged, chosen or repaired, and there is no distance tolerance:
  near-coincident atoms are accepted;
- **a residue number that cannot be read** (hybrid-36, written by some tools
  past 9,999 residues) on a record the selection would include. Hybrid-36 is
  not decoded in 1.0; renumber the file, or leave those residues out;
- **an atom the atom table cannot label**, naming the atom, its residue and its
  element. There is no option to drop it and continue; an atom the selection
  does not use (an unrequested atom of a narrowed ligand) is not labelled;
- **a protein–protein cross of a chain with itself** (`protein_protein_cross`,
  `protein_protein_calpha_cross`), and **any cross whose two sides share an
  atom** — for example `entity_name="TPO"` on a chain whose TPO is also a
  protein residue. A molecule side may share the protein's chain identifier;
- **SSE annotation of a chain whose `HELIX`/`SHEET` record cannot be read
  exactly**: an endpoint absent from the chain, written in two places, with an
  unreadable number, on another chain, or after the end
  ([SSE-TSR](methods/SSE_TSR.md), §2). The same calculation without `sse` is
  not affected.

Other fixed readings of unusual input, listed so that none is a surprise:

- **letter case** — element symbols are compared without regard to case (`h`
  is hydrogen); alternate-location indicators are compared exactly as written
  (an altLoc `a` is not `A`; blank and `A` are kept);
- **the Cα atom** is the atom named exactly `CA` in columns 14–15; `CA1` is not
  a Cα, and nothing infers one;
- **residue names differ between protein sides.** The protein side of
  `drug_protein_cross` and `protein_protein_cross` reads every `ATOM` residue of
  its chain (and TPO, SEP, PTR); the protein side of `carbohydrate_protein_cross`
  and `nucleotide_protein_cross` reads only the 20 regular amino acids and TPO,
  SEP, PTR, the vocabulary of their historical contact programs, and names any
  other residue of a whole chain as `not supported` (the result is partial).
  The amino-acid atomic methods read every `ATOM` residue a request names, but
  their completeness rule admits only the 23 residues with defined atom
  counts. Atomic support is not widened because a residue has a Cα label;
- **terminal residues** (ADR-111). Cα methods keep the first and last residues'
  Cα points, with no completeness or backbone requirement. The per-residue and
  residue-type atomic methods keep their exact-atom-count rule, so a
  C-terminal residue carrying `OXT` (1UBQ GLY 76) or a 5′ nucleotide without
  its phosphate is excluded and named. Cross methods apply no completeness
  rule: a terminal residue's atoms, `OXT` included, are used like any other;
- **`nucleotide_chain`** has no atom-count rule since 1.0.0rc7 (ADR-114): each
  standard nucleotide takes part with the heavy atoms the file has, at a chain
  end or inside the chain, and every nucleotide whose atoms differ from the
  standard set is named in a warning (`DC_1 lacks P, OP1, OP2`). The warning
  describes the file; nothing is reconstructed, and a chain with incomplete
  nucleotides is computed from the atoms present. The residue-specific
  `nucleotide` methods keep the exact count;
- **the narrowed Cα cross** (`calpha_residues`, `calpha_residues_2`) decides
  interface residues from the heavy atoms of the selected residues only, on
  both sides (ADR-109).

- **atomic labels come from the element field** (ADR-113). Every atomic method
  labels a heavy atom by its element field (columns 77–78, any letter case);
  the atom name never chooses a label. An element field that is wrong in the
  file gives a wrong label: BioTSR does not compare it with the atom name,
  with a residue's standard atom names or with bond lengths, and corrects
  nothing. A blank field, or an element the table does not list (`X`), refuses
  the calculation that would use the atom. Atomic results computed before
  1.0.0rc7 labelled atoms by name and are in another key space: they are not
  reused and never compared with current ones — recompute them.

## 3. Platforms: Linux and macOS

BioTSR is tested on Linux and macOS: keys, frequencies, bins, categories,
counts and fingerprints agree exactly on both
([platform validation](PLATFORM_VALIDATION.md)). The raw floating-point
geometry columns of a triplet file (angles, longest edge) are computed with the
host's C maths library and may differ in their last digits between platforms;
no class, key or count depends on that difference.

Windows is not a supported platform for 1.0: BioTSR is neither tested nor
validated there.

## 4. Not part of 1.0

- **Feature selection and key combining before a comparison.** The parallel
  TSR comparison program of the `hsp70_actin` work offered two options before
  computing a measure: keeping only keys whose frequencies vary little across
  samples (a median-absolute-deviation threshold together with a minimum
  number of absences), and summing adjacent keys. BioTSR computes every
  measure on the keys as they are, which is that program's default; results
  produced with either option cannot be reproduced with 1.0.
- **Candidate-strategy extensions** — for example CA-kNN and other
  nearest-neighbour candidate strategies. They are absent from the method
  registry rather than stubbed; every 1.0 method draws all C(n,3) triangles.
- **mmCIF-only entries.** `biotsr fetch` and `--download-missing` download
  PDB-format files only. Entries RCSB serves only as mmCIF (typically very large
  structures) are reported as `pdb_format_unavailable`; convert them to PDB
  format yourself. Retrieval is opt-in: nothing uses the network unless asked.
- **Automatic chain selection.** Every run names its chain — in the row, through
  `--column chain=HEADER`, or with `--chain` for rows that leave it blank. There
  is no option to run every chain of a structure; BioTSR never picks a chain.
- **Automatic choice of a cross side from a wide sheet.** When a row of a
  `--sheet tsr` sheet fills several columns that could narrow a side of a
  cross method (`aa_ind`, `aa_ind_specific`, `drug_atom_selected`, …), BioTSR
  neither picks one nor ignores them: the run is refused until each such side
  is chosen with `--cross-side METHOD.SIDE=COLUMN` or `=whole` (ADR-104).
- **Molecules on other chains.** `drug`, `carbohydrate`, `nucleotide`, their
  selections and the entity side of the cross methods search only the chain
  named. A name without `--entity-id` selects every copy of that name **on
  that chain**, in one combined calculation; with `--entity-id` it selects that
  one copy (ADR-096). Copies on other chains are never added.
- **Ordinal molecule identifiers.** `--entity-id` is a PDB residue identifier
  (`404`, `52A`, `-5`); "the second copy" cannot be asked for, and an
  identifier whose records appear in two separate places on the chain is
  refused as ambiguous.
- **Insertion codes in `aa_selection` ranges** (ADR-102). A range is two
  plain residue numbers, ascending, both ends included (`10-12`, `-5--3`), and
  expands to plain numbers only. A range with an insertion code at either end
  (`27A-29`), a descending range (`12-10`) and any other malformed item are
  refused, naming the item. A range that spans a residue with an insertion
  code in the structure — `9-11` on a chain that holds `10A` — is refused too,
  unless that residue is also listed: whether a numeric range should include
  it has no single answer, so BioTSR asks. List it to include it (`9-11;10A`),
  or list the residues one by one to leave it out (`9;10;11`).
- **SSE outside the protein Cα methods.** An atom, a ligand, a glycan or a
  nucleotide has no secondary structure of its own; `sse` is refused by name.
- **Mirror and SSE for cross methods.** No historical cross implementation has
  either, and neither is extended by analogy.
- **Similarity and clustering on the command line.** `biotsr analyze` offers
  common keys; feature matrices, similarity measures, clustering and the
  clustered heat map are available from Python (`biotsr.analysis`,
  [tutorial 08](tutorials/08_similarity_clustering.md)). The figure needs the
  `plot` extra (matplotlib).
- **Analysis size.** Similarity and clustering run in one process, in time
  quadratic in the number of samples, holding every sample's key counts in
  memory; `--jobs` does not parallelise them. Measured timings are in
  [methods/ANALYSIS.md](methods/ANALYSIS.md), §5.8; nothing above 400 samples
  (similarity) or 2,000 samples (clustering) was measured.
- **Results and outputs written by 1.0.0rc3 or earlier are not read.** Their
  per-sample directories with a `run.json` inside, their files named with bin
  counts (`1UBQ_A.keys_29_35`) and their schema-1 sidecars are refused by
  `biotsr.load_result`, `biotsr analyze --from-results` and `biotsr describe`,
  with a message saying they were written by a release candidate. No reader for that layout is
  provided: the selection record, the alternate-location rule of
  `amino_acid_atomic` and the meaning of an omitted molecule identifier all
  changed. Recompute them into a new `--output` folder.

## 5. Memory: not every run streams

`biotsr run`, `biotsr batch` and `batch.run(...)` write triplet files as the
triangles are computed, so a large chain's triplet list is never held in
memory. Two cases still hold every triplet record until the run finishes: the
per-residue sweep of `amino_acid_atomic` (no `--residue`), and the Python calls
`analyze`, `analyze_cross` and `analyze_entities` with the default
`keep_triplets=True`. Their memory grows with the number of triangles; use
`--no-triplets`, or `keep_triplets=False` with a `triplet_sink` (API guide),
for very large structures.

Two other costs do not depend on how a run is made. A run's key counts are held
until it finishes (they grow with the number of distinct keys, which is far
smaller than the number of triangles), and analysis — `biotsr analyze
common-keys` and `biotsr.analysis` — holds every sample's key counts at once,
so its memory grows with samples × distinct keys per sample. Streaming lowers
memory, not output size: a triplet file is as large either way, and it is by
far the largest output. The samples that `analyze common-keys` computes from a
manifest are written keys only, with no triplet file (CLI guide, *Memory and
disk*).

## 6. Interruption: what cannot be caught

Ctrl-C and SIGTERM are handled by every command: a batch interrupted at any
point after planning writes its summary, `biotsr fetch` writes its report, and
the exit status is 130. A worker process that dies is reported as a failure
with the summary written (exit 3). Three things are outside BioTSR's control:

- **SIGKILL and power loss** end the process before it can write anything. Its
  staging directory (`.biotsr-staging-<pid>-<host>-<token>`) and per-prefix
  locks stay behind, and a calculation interrupted while its files were being
  installed keeps an incomplete completion record (`"complete": false`) that
  every reader refuses — never a result that looks complete. A later run on the
  same machine removes the leftovers once that process is gone and recomputes
  an incomplete calculation; leftovers of a process on another machine (a
  shared file system) cannot be checked, so they are kept and reported.
- **No total time limit on a download.** `--timeout` limits each wait — for a
  connection, or for the next piece of a response — not a whole download: a
  server that keeps sending, however slowly, is never timed out, and nothing
  bounds how long one download (or its retries, after backoff waits of up to
  60 s) can take. Use Ctrl-C or SIGTERM, which stop it promptly.
- **Two steps of a download are not interruptible**: a host-name lookup (the
  system resolver's own timeouts apply) and a TLS handshake step in progress
  (bounded by `--timeout`). Ctrl-C waits for them; every other phase — a
  connection being made, a stalled or slowly arriving response, a backoff wait,
  queued ids — stops at once.

## 7. Cosmetic

`biotsr analyze common-keys --output DIR` writes into `DIR/analysis/common_keys/`,
so `--output analysis/` produces `analysis/analysis/common_keys/`. This is the
documented output contract and is left unchanged for 1.0.
