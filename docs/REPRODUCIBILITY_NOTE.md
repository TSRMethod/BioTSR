# Reproducibility note — which command produces which method

A short mapping from the scientific methods described in
[`methods/ANALYSIS.md`](methods/ANALYSIS.md),
[`methods/CROSS_TSR.md`](methods/CROSS_TSR.md) and
[`methods/SSE_TSR.md`](methods/SSE_TSR.md) to the commands that run them, and a
statement of what a saved result records.

**Those three pages are authoritative for the science.** Nothing here defines a
method, a convention or a parameter; this is a lookup table, and command syntax
is deliberately not part of the method descriptions.

---

## Environment

This note describes BioTSR 1.0.0.

```bash
pip install "biotsr @ git+https://github.com/TSRMethod/BioTSR@v1.0.0"   # or, once https://pypi.org/project/biotsr/ lists 1.0.0: pip install biotsr
biotsr --version        # biotsr 1.0.0
```

The calculation uses the standard library only (the installation's one
dependency, joblib, runs batch rows in parallel and changes no result). Every
number below is reproducible offline once BioTSR is installed, given the
structure files. Linux and macOS are supported; input is PDB format.

`python -m biotsr` is equivalent to `biotsr` and is what to use when the environment's
`bin` directory is not on `PATH`.

---

## One method, one command

| Method (as the Methods documents name it) | Command |
|---|---|
| Protein Cα TSR | `biotsr run -m protein --pdb P.pdb -c A -o OUT` |
| …with amino-acid grouping | add `--aa-grouping` |
| …with the mirror convention | add `--mirror` |
| …with size filtering at *d* Å | add `--size-filter d` |
| Cα TSR over selected residues | `biotsr run -m aa_selection ... --residues '10-12;48;70-71'` |
| Amino-acid atomic TSR, one residue | `biotsr run -m amino_acid_atomic ... --residue 15 --residue-name LEU` |
| Amino-acid atomic TSR, every residue of a chain | `biotsr run -m amino_acid_atomic ... ` (omit `--residue`) |
| Amino-acid atom selection | `biotsr run -m amino_acid_atom_selection ... --atoms 'LEU_8_CD1;LYS_48_NZ;HIS_68_NE2'` |
| Amino-acid type selection | `biotsr run -m aa_ind ... --residue-types TYR` (one type) or `-m aa_ind_mix ... --residue-types 'TYR;PHE'` |
| Specific amino-acid instances | `biotsr run -m aa_ind_specific ... --residues 'LYS_6;GLU18'` |
| Ligand TSR | `biotsr run -m drug ... --entity-name BEN --entity-id 1` |
| Carbohydrate TSR, one monomer | `biotsr run -m carbohydrate ... --entity-name NAG --entity-id 1` |
| Carbohydrate TSR, every copy of one sugar on the chain | `biotsr run -m carbohydrate ... --entity-name NAG` (omit `--entity-id`) |
| Carbohydrate TSR, whole chain | `biotsr run -m carbohydrate_chain ... ` |
| Nucleotide TSR, one residue | `biotsr run -m nucleotide ... --entity-name DC --entity-id 3` |
| DNA/RNA TSR, whole strand | `biotsr run -m nucleotide_chain ... ` |
| Ligand–protein cross TSR | `biotsr run -m drug_protein_cross ... --chain-2 B --entity-name BEN --entity-id 1 --contact-cutoff 3.5` |
| Carbohydrate–protein cross TSR | `biotsr run -m carbohydrate_protein_cross ... --chain-2 B --entity-name NAG --entity-id 1` |
| Nucleotide–protein cross TSR | `biotsr run -m nucleotide_protein_cross ... --chain-2 B --entity-name DG --entity-id 8` |
| Protein–protein cross TSR (all atom) | `biotsr run -m protein_protein_cross ... --chain-2 D --contact-cutoff 5.0` |
| Protein–protein cross TSR (Cα interface) | `biotsr run -m protein_protein_calpha_cross ... --chain-2 D --contact-cutoff 5.0` |
| SSE-TSR (18 categories) | add `--sse` to either Cα method |
| Common keys over a sample set | `biotsr analyze common-keys --manifest M.csv -o OUT` |

`biotsr methods NAME` prints any method's own description, the target arguments
it accepts and the standard profile it runs under. An argument a method does not
read is refused by name. For the single-entity methods and the entity side of
a cross, `--entity-id` names one copy; without it, every copy of
`--entity-name` on the named chain is analysed in one calculation (ADR-096).

`biotsr batch` also reads a wide dataset sheet (`--sheet tsr` or `--sheet
custom`), with the methods named on the command line (`--method`, repeatable)
rather than in the sheet: each method reads its own columns, and a row that
does not describe a method is failed or, with `--inapplicable skip`, recorded
as `inapplicable` (ADR-097; `CLI_GUIDE.md`, §7).

---

## Which convention ran

Every method runs under one profile, its standard profile; there is no profile
or output-format option (ADR-087). The method name therefore determines the
conventions, and the result records them in full: the profile's name in
the completion record `PREFIX.run.json` and the sidecar's provenance block, and every key-defining and
sampling choice in the two fingerprinted identities below.

Where a standard convention differs from a historical TSR program's — for
example the angle divisor (`math.pi` against the literal `3.14`), the top Cα
distance class, the grouped radix, the cysteine sulfur label or β-strand
identity — the difference is stated in the methods pages and the
[architecture decisions](architecture/DECISIONS.md). BioTSR does not reproduce
the historical programs' output files.

---

## What a reader of the output can recover

Every calculation writes its files into the one `--output` folder under one
prefix — `PREFIX.keys`, `PREFIX.meta.json`, `PREFIX.triplets` (unless
`--no-triplets`) and `PREFIX.run.json` — and its two records leave nothing about
the run to be reconstructed by inference:

* `PREFIX.meta.json` — the full `KeySpace` and `SamplingProfile` with their
  fingerprints, the profile that ran, the selection record (per side: chain,
  scope, arguments and every molecule or residue instance used), the source
  file's SHA-256, the parser backend, the alternate-location rule, the
  modified-residue rule (`source.modified_residue_policy`: which modified amino
  acids written as `HETATM` were read) and the input rules
  (`source.input_rules`: how atom records were read, ADR-108), the
  candidate ledger, the SSE assignment and coverage, every skip and every
  warning (result schema 2);
* `PREFIX.run.json` — the completion record: the command, the package version,
  the resolved profile, the SSE strand-identity convention, the arguments the
  user supplied, the source manifest rows for a batch, the normalised request
  together with the processing rules it was computed under
  (`request.processing`), and per result its files with their sizes, counts,
  fingerprints, warnings and `omissions` (what the selection asked for and did
  not get, `{}` when nothing). It is written last, and the command line never
  reads a calculation whose record is missing or incomplete. `--if-exists
  skip` reuses a calculation only when its recorded request, processing rules
  included, equals the new one; a record without processing rules, written
  before they were recorded, is recomputed with `--if-exists overwrite`, never
  reused.

The prefix is built from the request —
`structure_chain[_chain2]_method[_target][_modifiers]`, for example
`7LF9_D_carbohydrate_NAG_1` or `1UBQ_A_protein_aagroup_mirror_size25` — so the
same request always gets the same name; the bin counts are in the metadata, not
the name.

```bash
biotsr describe OUT/                  # one line per calculation
biotsr describe OUT/PREFIX            # how it was produced
biotsr describe OUT/PREFIX --json     # the completion record itself
```

Two fingerprints are recorded on every result:

* **`KeySpace`** — the angle convention, both bin schemes, the label scheme, the
  square-root and median conventions, the canonical ordering, the mirror
  convention, the key formula and the numerical precision. Two results with
  different `KeySpace` fingerprints have keys that do not mean the same thing,
  and the analysis layer refuses to combine them.
* **`SamplingProfile`** — the selector, the scope, the completeness policy, the
  candidate strategy, the size filter, the contact cutoff, filter rule and
  projection, and the degenerate policy, together with each side's chain and
  description.

Every analysis — feature matrices, common keys, similarity and clustering —
requires one `KeySpace`, decided from the key definitions and never from method
names. The processing settings among the sampling fields, with the
alternate-location, modified-residue and input rules, must also agree unless
mixed
processing is explicitly allowed (`--allow-mixed-processing`,
`allow_mixed_processing=True`), and the report then records every sample's
settings and a warning; a result that does not record the modified-residue
rule is reported as *not recorded*, never as the current rule. A partial
sample — one whose selection did not get everything it asked for — is
refused unless it is acknowledged separately (`--allow-partial`,
`allow_partial=True`), and the report then records what each lacks
(ADR-106). What was selected — structures, chains, residue types,
molecules, residues, atoms, one copy or all copies — never blocks a
comparison (ADR-094). The full `SamplingProfile`
fingerprint is reported but not used to decide comparability, because it
includes chain names and free-text side descriptions.

One convention sits in **neither**, by design: the β-strand identity that names
an SSE category (ADR-063). The completion record and the sidecar record which convention
ran, and `biotsr.sse.require_common_strand_identity` is what any category-level
comparison must call.

---

## A dataset

One manifest row per run, columns named after the flags:

```csv
sample_id,pdb_id,method,chain,label,sse
p53,1TUP,protein,A,tumour_suppressor,true
barnase,1BRS,protein,A,rnase,true
barstar,1BRS,protein,D,inhibitor,true
```

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/
```

`results/batch_summary.csv` has one row per **request**, listing the manifest
rows it came from, with its status: `ok`, `partial`, `no_triangles`,
`skipped`, `inapplicable`, `failed`, `interrupted` or `not_started`. Failures
additionally appear in `results/batch_failures.csv` with their reason, and the
process exits 3 when any row failed, was interrupted by a worker failure, or
never started. Otherwise it exits 4 when a result is `partial` — its selection
did not get everything it asked for, and its message says what is missing —
or is `no_triangles` with such a selection, or when `--if-exists skip` kept a
partial result, and 0 when every result is complete. No sample is ever silently omitted. A `sample_id` replaces the output prefix and
must be a valid file-name prefix inside `--output` (ADR-098).

Full details: `CLI_GUIDE.md`.

---

## Reproducing an interpretation without rerunning the pipeline

A saved native result carries the provenance needed to *interpret* it, so a
result can be re-analysed — months later, on another machine, by someone else —
without re-reading the structure file.

```bash
biotsr run -m protein --pdb 1TUP.pdb -c A -o results/ --no-triplets     # once
biotsr analyze common-keys --from-results results/ -o analysis/         # thereafter
```

or from Python:

```python
import biotsr
loaded = biotsr.load_result("results/1TUP_A_protein.keys")
```

What a loaded result restores, exactly:

* every key and its per-category frequency;
* **both fingerprints**, recomputed bit-for-bit from the recorded semantics, so
  the analysis layer's compatibility rules apply to a loaded result exactly as
  they apply to a freshly computed one. Loading launders nothing;
* the **profile** that produced it — its name, the profile it was derived from
  and which modifiers were active;
* the **SSE annotation provenance** — the assignment version, the β-strand and
  helix identity conventions, the `HELIX`/`SHEET` record counts and
  `has_records`, so a "this structure is all coil" reading can still be told from
  "this structure carried no annotation";
* the **source** — the structure's path, its SHA-256 as it was read, the parser
  backend, the model and the altLoc policy (`blank_or_a` for every method);
* the selection summary, the candidate ledger, the skip records and the warnings.

Two things a loaded result deliberately does not carry, and says so on the
object rather than pretending: the per-triangle records (`alt_loc` is per-atom
and is in no column of the triplet file) and the first-seen key order (which no
BioTSR output uses). Re-run the calculation for either; the keys are identical.

The format states its own version, `result_schema_version` (2), so a build
either reads a file or refuses it by name — it never reads one schema as
another (ADR-074, ADR-078). Results and output folders written by early
pre-release builds (per-sample directories, `*.keys_29_35` names, schema 1)
are refused with a message saying so; recompute them.

Anything **not** recorded is not recoverable, and that is the boundary to state
in a methods section: a saved result reproduces the *interpretation*, not the
computation. Reproducing the computation needs the structure file, whose SHA-256
the sidecar carries so the right one can be identified.

### Measuring it

`validation/scripts/performance_api.py` measures, on pinned structures:
compute, serialize and deserialize time; the peak Python allocation of each; and
the practical comparison — recompute-then-analyse against load-then-analyse. It
writes its measurements where `--json` points; no measurement is committed,
because timings belong to the machine they were taken on.

The ratio between the two paths is a property of the structure, not of the
package: the recompute side scales with C(n,3) in selected points and the load
side with the number of distinct keys. It is reported as a measurement on named
structures and is not a general speedup claim.
