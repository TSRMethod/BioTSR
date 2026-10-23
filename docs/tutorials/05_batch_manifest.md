# 05 — A dataset manifest

For more than a handful of structures, describe the dataset in a **manifest**:
a CSV where **one row is one run** and the columns are the `biotsr run` flags
without their dashes.

<!-- file: samples.csv -->
```csv
sample_id,pdb_id,method,chain,label
ubiquitin,1UBQ,protein,A,ubiquitin
barnase,1BRS,protein,A,barnase
barstar,1BRS,protein,D,barstar
```

`pdb_id` is resolved against `--pdb-dir` (`1UBQ.pdb`, `1ubq.ent`, …).

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --no-triplets
```

```text
manifest   samples.csv (3 source rows, 3 requests)
requests   3 ok, 0 skipped, 0 failed
summary    batch/batch_summary.csv
```

Every calculation writes into the one `--output` folder, under its own
**prefix** — here the row's `sample_id` — so `batch/` now holds
`ubiquitin.keys`, `ubiquitin.meta.json` and `ubiquitin.run.json` (the
completion record), the same for `barnase` and `barstar`, and
`batch_summary.csv`. `--no-triplets` skips the per-triangle files, which are
the large ones.

`batch_summary.csv` has one row per request, in order: `request`, `rows` (the
manifest line numbers it came from; line 1 is the header), `participant`,
`sample_id` (the output prefix), `status`, `method`, `profile`,
`structure_id`, `chain`, `chain_2`, `selection`, `n_results`, `n_triangles`,
`n_unique_keys`, `output` (folder and prefix), `message` and `chain_source`.

```python
import csv

with open("batch/batch_summary.csv") as fh:
    rows = list(csv.DictReader(fh))
assert [r["sample_id"] for r in rows] == ["ubiquitin", "barnase", "barstar"]
assert [r["rows"] for r in rows] == ["2", "3", "4"]
assert [r["n_triangles"] for r in rows] == ["70300", "204156", "105995"]
assert [r["output"] for r in rows] == ["batch/ubiquitin", "batch/barnase", "batch/barstar"]
```

Three rules keep a manifest honest:

1. a blank cell means *not set*, never a default;
2. an unknown column is an error that says how to give it a role;
3. by default a failing row is recorded in `batch_failures.csv` and the rest
   still run (`--on-error stop` halts instead).

A `sample_id` is one file-name stem inside `--output` (no `/`, no leading
`.`): it becomes the prefix of that row's files. Without a `sample_id` column
the prefix is built from the request (`1BRS_D_protein`, below). A list inside
one cell separates its items with `;`.

## A dataset sheet with columns of its own

An existing sheet rarely uses BioTSR's column names. Here `protein` (an
accepted historical name for `pdb_id`) and `chain` are understood as they are;
`group` is a class label and `group1` a numeric code for it. The method is
not in the sheet; it comes from `--method`:

<!-- file: sample_details.csv -->
```csv
protein,chain,group,group1
1UBQ,A,ubiquitin,0
1BRS,A,barnase,1
1BRS,D,barstar,2
```

Give each extra column a role — read `group` as the label, carry `group1` along
as metadata — and check the plan first. `--dry-run` computes and writes
nothing, and lists every request with its source row and output prefix:

```bash
biotsr batch --manifest sample_details.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --output sheet/ --no-triplets --dry-run
```

```text
chains     3 from the row, 0 from --chain
outcome    3 would run, 0 would be downloaded then run, 0 would be skipped (verified), 0 inapplicable, 0 chain not specified, 0 would fail
3 4 protein - D - 1BRS_D_protein [would run]
```

Then run it:

```bash
biotsr batch --manifest sample_details.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --output sheet/ --no-triplets
```

```text
requests   3 ok, 0 skipped, 0 failed
```

The kept column appears in the summary and in each completion record, and the
label never changes a calculation:

```python
import json

with open("sheet/batch_summary.csv") as fh:
    rows = list(csv.DictReader(fh))
assert [r["sample_id"] for r in rows] == ["1UBQ_A_protein", "1BRS_A_protein", "1BRS_D_protein"]
assert [r["group1"] for r in rows] == ["0", "1", "2"]
assert [r["chain_source"] for r in rows] == ["row", "row", "row"]
with open("sheet/1BRS_D_protein.run.json") as fh:
    run = json.load(fh)
assert run["complete"] is True and run["prefix"] == "1BRS_D_protein"
assert run["label"] == "barstar" and run["metadata"] == {"group1": "2"}
assert run["source"]["rows"] == [4]
assert run["request"]["modifiers"] == {}
```

Running the same command again with `--if-exists skip` keeps every calculation
whose completion record shows exactly the same request, and computes nothing:

```bash
biotsr batch --manifest sample_details.csv --method protein --column label=group --keep-column group1 --pdb-dir pdb/ --output sheet/ --no-triplets --if-exists skip
```

```text
requests   0 ok, 3 skipped, 0 failed
```

## A wide dataset sheet

A wide sheet has **one row per structure and the columns of many calculations
at once**: several chain columns (`aa_chain`, `drug_chain`, `carb_chain`,
`nt_chain`, `aa_chain_2`) and several alternative selections. This one is
`tests/data/manifests/wide_sheet.csv` in the repository, a synthetic sheet in
the TSR dataset-sheet layout. Its structures, 6G9M, 7LF9 and 2AC0 besides
3PTB and 1BRS, are under `tests/data/structures/` (see the
[tutorials README](README.md)):

<!-- file: wide_sheet.csv -->
```csv
protein,aa_chain,drug_chain,carb_chain,nt_chain,aa_chain_2,aa_group,aa_ind,aa_ind_mix,aa_ind_specific,aa_ind_2,aa_ind_mix_2,aa_ind_specific_2,drug_name,drug_id,carb_name,carb_id,nt_name,nt_id,aa_selection,aa_selection_2,aa_atom_selection,aa_atom_selection_2,drug_atom_selected,car_atom_selected,nt_atom_selected,comment_1
6G9M,A,A,,,,kinase,TYR,TYR;ARG,ARG_15;GLU12,,,,ESW,404,,,,,12-15;20,,TYR_30_OH;TYR_25_CZ,,C3;C5;O6,,,ligand row
3PTB,A,A,,,,protease,TYR,,,,,,BEN,1,,,,,,,,,,,,ligand row in a second structure
7LF9,B,,D,,,glycoprotein,TYR,,,,,,,,NAG,1,,,,,,,,C1;C2;N2,,one NAG of the glycan
7LF9,B,,D,,,glycoprotein,TYR,,,,,,,,NAG,2,,,,,,,,C1;C2;N2,,another NAG of the same glycan
7LF9,B,,D,,,glycoprotein,TYR,,,,,,,,NAG,5,,,,,,,,C1;C2;N2,,a NAG on the other branch
7LF9,B,,D,,,glycoprotein,TYR,,,,,,,,MAN,,,,,,,,,,,every MAN of the glycan: no identifier
1BRS,A,,,,D,barnase-barstar,LYS,LYS;ARG,LYS_19;TYR24,ASP,ASP;GLU,LYS1;ARG_11,,,,,,,1-30,1-20,,,,,,two protein chains
2AC0,A,,,E,,DNA-binding,TYR,,,,,,,,,,DT,7,,,,,,,P;C1';N1,nucleotide row
```

Which calculations to run is not in the sheet: you choose them with
`--method`, which is required with `--sheet` and repeatable, and each method
reads **only its own columns**. `--sheet tsr` knows the TSR layout:
`carbohydrate` reads `carb_chain`, `carb_name` and `carb_id`; `aa_ind` reads
`aa_chain` and `aa_ind`; the structure is `protein` and the label `aa_group`
(an annotation only). Preview first:

```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method carbohydrate --method aa_ind --inapplicable skip --pdb-dir pdb/ --output wide/ --no-triplets --dry-run
```

```text
manifest   wide_sheet.csv (8 source rows, 13 requests)
outcome    9 would run, 0 would be downloaded then run, 0 would be skipped (verified), 4 inapplicable, 0 chain not specified, 0 would fail
1 2 carbohydrate 1 - - - [inapplicable: the row has no carb_name / carb_id for carbohydrate]
5 4 carbohydrate 1 D entity_name=NAG; entity_id=1 7LF9_D_carbohydrate_NAG_1 [would run]
6 4;5;6;7 aa_ind 1 B residue_types=TYR 7LF9_B_aa_ind_TYR [would run]
9 7 carbohydrate 1 D entity_name=MAN 7LF9_D_carbohydrate_MAN_allcopies [would run]
```

Each line is one request: its number, the source row(s), the method, the
participant, the chains, the selection, and the output prefix with the
outcome. Four things to read in it:

* **Rows that do not describe a calculation.** The ligand rows (2 and 3) and
  the rows without a sugar (8 and 9) have no `carb_name` or `carb_id`, so they
  do not describe a `carbohydrate` calculation. By default such a row fails;
  `--inapplicable skip` records it as `inapplicable` instead. It is never
  silently left out.
* **Shared calculations.** Rows 4 to 7 all ask for TYR on 7LF9 chain B; that
  is one calculation, computed once, listing every source row (`4;5;6;7`).
* **One copy or all copies.** Rows 4, 5 and 6 name NAG 1, 2 and 5: three
  calculations, one copy each. Row 7 names MAN without an identifier: **every
  MAN on chain D** (MAN 4 and 7) in one calculation, `…_MAN_allcopies`. That is
  not the whole glycan, which is the `carbohydrate_chain` method.
* **Unused columns.** Columns that no selected method reads are listed, so a
  column you meant to use but no method reads is visible:

```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method carbohydrate --method aa_ind --inapplicable skip --pdb-dir pdb/ --output wide/ --no-triplets
```

```text
unused     aa_atom_selection, aa_atom_selection_2, aa_chain_2, aa_ind_2, aa_ind_mix, aa_ind_mix_2, aa_ind_specific, aa_ind_specific_2, aa_selection, aa_selection_2, car_atom_selected, drug_atom_selected, drug_chain, drug_id, drug_name, nt_atom_selected, nt_chain, nt_id, nt_name (read by none of the selected methods)
requests   9 ok, 0 skipped, 0 failed, 4 inapplicable
```

The `comment_1` column is ignored, as every `comment*` column is.

```python
with open("wide/batch_summary.csv") as fh:
    wide = list(csv.DictReader(fh))
assert len(wide) == 13
inapplicable = [(r["rows"], r["method"]) for r in wide if r["status"] == "inapplicable"]
assert inapplicable == [("2", "carbohydrate"), ("3", "carbohydrate"),
                        ("8", "carbohydrate"), ("9", "carbohydrate")]
by_prefix = {r["sample_id"]: r for r in wide if r["status"] == "ok"}
assert by_prefix["7LF9_B_aa_ind_TYR"]["rows"] == "4;5;6;7"
assert by_prefix["7LF9_D_carbohydrate_MAN_allcopies"]["n_triangles"] == "1540"
assert by_prefix["7LF9_D_carbohydrate_NAG_1"]["n_triangles"] == "364"

# The completion record names its source row, participant and the columns read.
with open("wide/7LF9_D_carbohydrate_MAN_allcopies.run.json") as fh:
    man = json.load(fh)
assert man["source"]["rows"] == [7] and man["source"]["participant"] == "1"
assert man["source"]["columns"] == {"chain": "carb_chain", "entity_name": "carb_name"}
assert man["results"][0]["target"] == "MAN_allcopies"
```

The folder itself lists its calculations, one line each — prefix, method,
triangles and whether the calculation is complete:

```bash
biotsr describe wide/
```

```text
7LF9_B_aa_ind_TYR aa_ind 95284 triangles complete (keys only)
7LF9_D_carbohydrate_MAN_allcopies carbohydrate 1540 triangles complete (keys only)
7LF9_D_carbohydrate_NAG_1 carbohydrate 364 triangles complete (keys only)
```

### The second participant

The amino-acid methods have a second participant in a TSR sheet, the `_2`
columns (`aa_chain_2`, `aa_ind_2`, …). `--participant 2` reads those instead
(`all` reads both). Only the 1BRS row fills `aa_ind_2`; it names ASP on chain
D. One `batch_summary.csv` describes one batch, so this run gets its own
folder:

```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method aa_ind --participant 2 --inapplicable skip --pdb-dir pdb/ --output wide_p2/ --no-triplets
```

```text
requests   1 ok, 0 skipped, 0 failed, 7 inapplicable
```

```python
with open("wide_p2/batch_summary.csv") as fh:
    (p2,) = [r for r in csv.DictReader(fh) if r["status"] == "ok"]
assert (p2["rows"], p2["participant"], p2["sample_id"]) == ("8", "2", "1BRS_D_aa_ind_ASP")
```

### Cross methods: choose each side

A cross method has two sides, and a TSR sheet can narrow each of them in
several ways: the protein side of `drug_protein_cross` by `aa_ind`,
`aa_ind_mix`, `aa_ind_specific` or `aa_atom_selection`, the drug by
`drug_atom_selected`. Those columns were written for other methods, and the
sheet does not say which one a cross calculation means. BioTSR does not pick
one, and it does not silently ignore them either: while a row populates a
side's columns, the run is refused until you choose.

```python
import biotsr.batch

try:
    biotsr.batch.run("wide_sheet.csv", "cross/", pdb_dir="pdb/", sheet="tsr",
                     methods=["drug_protein_cross"], inapplicable="skip", dry_run=True)
except biotsr.batch.BatchError as refused:
    message = str(refused)
assert "the drug side of drug_protein_cross: drug_atom_selected on line 2" in message
assert "the protein side of drug_protein_cross: aa_atom_selection on line 2" in message
```

`--cross-side METHOD.SIDE=COLUMN` narrows a side by that column, and
`--cross-side METHOD.SIDE=whole` uses the whole side. The preview prints both
sides of every request, and the populated columns the choice leaves unused:

```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method drug_protein_cross --cross-side drug_protein_cross.drug=whole --cross-side drug_protein_cross.protein=whole --inapplicable skip --pdb-dir pdb/ --output cross/ --no-triplets --dry-run
```

```text
drug, chain A: ESW copy 404, all atoms
protein, chain A: whole chain
not applied: drug_atom_selected (the drug side is used whole); aa_ind (the protein side is used whole); aa_ind_mix (the protein side is used whole); aa_ind_specific (the protein side is used whole); aa_atom_selection (the protein side is used whole)
```

The sides are `drug`, `carbohydrate` or `nucleotide` and `protein` for the
molecule-protein methods, and `1` (the `aa_*` columns) and `2` (the `aa_*_2`
columns) for `protein_protein_cross` and `protein_protein_calpha_cross`. A row
whose chosen column is blank does not describe the narrowed calculation: it is
inapplicable, never computed on the whole side. Computed:

<!-- not-executed -->
```bash
biotsr batch --manifest wide_sheet.csv --sheet tsr --method drug_protein_cross --cross-side drug_protein_cross.drug=whole --cross-side drug_protein_cross.protein=whole --inapplicable skip --pdb-dir pdb/ --output cross/ --no-triplets
```

This prints `requests 1 ok, 0 skipped, 0 failed, 1 partial, 6 inapplicable`
and exits `4`, because one result is partial (below, and *Partial results
and exit code 4*). The same from Python:

```python
report = biotsr.batch.run(
    "wide_sheet.csv", "cross/", pdb_dir="pdb/", sheet="tsr",
    methods=["drug_protein_cross"], inapplicable="skip", triplets=False,
    cross_sides=["drug_protein_cross.drug=whole", "drug_protein_cross.protein=whole"])
assert (report.n_ok, report.n_partial, report.n_inapplicable) == (1, 1, 6)
with open("cross/batch_summary.csv") as fh:
    cross = {r["rows"]: r for r in csv.DictReader(fh)}
assert cross["2"]["status"] == "partial"
assert cross["2"]["selection"] == ("drug, chain A: ESW copy 404, all atoms | "
                                   "protein, chain A: whole chain")
assert "not supported: CME_161" in cross["2"]["message"]
assert cross["3"]["status"] == "ok"
```

Row 2 is `partial`: chain A of 6G9M contains CME 161, a modified cysteine
stored as `HETATM`, which the atomic methods do not read (next section), so the
whole protein side is the chain without it, and the summary says so.

### Modified amino acids

RCSB files write modified amino acids as `HETATM` records. The amino-acid
atomic methods and the atomic cross sides read exactly the three
phosphorylated residues their definitions cover — phosphothreonine `TPO`,
phosphoserine `SEP` and phosphotyrosine `PTR`, which have exact atom counts
(11, 10, 16) and atom labels; any other modified residue, such as `CME`, is
named as **not supported**. The Calpha methods label residues instead, and
their residue labels cover eleven modified residues — `TPO`, `SEP`, `PTR`,
`NIY`, `M3L`, `TRF`, `CSO`, `LLP`, `CME`, `CSD`, and `MSE` with methionine's
code — so they read all of those; any other is named as **not represented**.
`SYNPHOS.pdb` is a synthetic structure written for this
(`tests/data/synthetic/`): chain A holds TPO 3, SEP 5, PTR 7 and CME 8.

<!-- not-executed -->
```bash
biotsr run --method aa_ind_mix --pdb SYNPHOS.pdb --chain A --residue-types "TYR;TPO;CME" --output modified/ --no-triplets
biotsr run --method protein --pdb SYNPHOS.pdb --chain A --output modified/ --no-triplets
```

The first warns `not supported: residue type CME (…)`, writes its result and
exits `4`, because the result is partial. The second reads all four modified
residues as Cα points, with no warning, and exits `0`. From Python:

```python
mixed = biotsr.analyze("SYNPHOS.pdb", method="aa_ind_mix", chain="A",
                       residue_types="TYR;TPO;CME", keep_triplets=False)
assert set(mixed.selection.residues_used) == {"3_TPO", "6_TYR"}
record = mixed.annotations["selection_request"]["sides"][0]
assert record["unsupported"] == ["CME"]
assert mixed.warnings[0].startswith("not supported: residue type CME (")
calpha = biotsr.analyze("SYNPHOS.pdb", method="protein", chain="A", keep_triplets=False)
assert calpha.warnings == ()
assert calpha.selection.residues_used[2:8] == ("3", "4", "5", "6", "7", "8")
```

A request of unsupported residues only is refused, naming them; in a batch, a
calculation that lacks something it named — not found, not supported, or
excluded by the completeness rule — is reported `partial`, never `ok`.

### Partial results and exit code 4

A result is **partial** when its selection did not get everything it asked
for: a requested residue, atom, molecule or residue type is not in the
structure, a modified residue is not supported or not represented, or a
residue or molecule copy you named is declined by the completeness rule. A
partial result is written and usable, but it is never presented as the
complete selection, and the command that produced it exits `4` instead of `0`
(a batch with a failed request still exits `3`).

Chain A of `SYNPHOS.pdb` has TYR, LYS and SER; chain B has TYR and LYS but no
SER:

<!-- file: partial.csv -->
```csv
pdb_id,chain,method,residue_types,sample_id
SYNPHOS,A,aa_ind_mix,TYR;LYS;SER,chainA
SYNPHOS,B,aa_ind_mix,TYR;LYS;SER,chainB
```

<!-- not-executed -->
```bash
biotsr batch --manifest partial.csv --pdb-dir pdb/ --output partial/ --no-triplets
```

`chainA` is `ok`; `chainB` is `partial`, with the message
`requested but not found: residue type SER`; the summary line reads
`requests 1 ok, 0 skipped, 0 failed, 1 partial`, and the command exits `4`. A
script tells the two outcomes apart by the exit status:

```python
import contextlib
import io

from biotsr.cli import EXIT_INCOMPLETE, EXIT_USAGE, main


def biotsr_cli(*argv):
    """Run one biotsr command in-process; return its exit status and output."""
    log = io.StringIO()
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        code = main(list(argv))
    return code, log.getvalue()


code, log = biotsr_cli("batch", "--manifest", "partial.csv", "--pdb-dir", "pdb/",
                       "--output", "partial/", "--no-triplets")
assert code == EXIT_INCOMPLETE == 4
assert "requests   1 ok, 0 skipped, 0 failed, 1 partial" in log
with open("partial/batch_summary.csv") as fh:
    status = {r["sample_id"]: (r["status"], r["message"]) for r in csv.DictReader(fh)}
assert status == {"chainA": ("ok", ""),
                  "chainB": ("partial", "requested but not found: residue type SER")}
with open("partial/chainB.run.json") as fh:
    assert json.load(fh)["results"][0]["omissions"]["primary"] == {"not_found": ["SER"]}
```

The completion record of each result lists its `omissions` (`{}` when it is
complete). Resuming keeps the partial result, and says so:

<!-- not-executed -->
```bash
biotsr batch --manifest partial.csv --pdb-dir pdb/ --output partial/ --no-triplets --if-exists skip
```

```python
code, log = biotsr_cli("batch", "--manifest", "partial.csv", "--pdb-dir", "pdb/",
                       "--output", "partial/", "--no-triplets", "--if-exists", "skip")
assert code == EXIT_INCOMPLETE
assert "requests   0 ok, 2 skipped (1 of them partial), 0 failed" in log
assert ("verified existing output of the same request; kept partial result: "
        "requested but not found: residue type SER") in log
```

An analysis refuses a partial sample until you acknowledge it with
`--allow-partial`, and then exits `4` ([tutorial 07](07_analysis.md)):

<!-- not-executed -->
```bash
biotsr analyze common-keys --from-results partial/ --output partial_ck/
biotsr analyze common-keys --from-results partial/ --output partial_ck/ --allow-partial
```

```python
code, log = biotsr_cli("analyze", "common-keys", "--from-results", "partial/",
                       "--output", "partial_ck/")
assert code == EXIT_USAGE == 2
assert "1 of the 2 samples are partial" in log and "chainB (not found: SER)" in log
code, log = biotsr_cli("analyze", "common-keys", "--from-results", "partial/",
                       "--output", "partial_ck/", "--allow-partial")
assert code == EXIT_INCOMPLETE
assert "common keys 258 over 2 samples" in log
```

### Two rows, one calculation: samples and computations

A sample is what one output folder entry holds: its prefix (or `sample_id`) and
its label. Two rows asking for **the same calculation of the same sample** are
computed once, and the summary lists both rows. Rows with **distinct
`sample_id` values** are distinct samples even when their calculation is
identical: each is computed and written under its own prefix, and an analysis
sees both. Rows that would be one calculation but carry **different labels**
are refused — one output cannot carry two labels.

<!-- file: duplicates.csv -->
```csv
pdb_id,chain,method,residue_types,sample_id,label
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,LYS,,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_a,ubiquitin
1UBQ,A,aa_ind,TYR,ubq_tyr_b,ubiquitin
```

```bash
biotsr batch --manifest duplicates.csv --pdb-dir pdb/ --output duplicates/ --no-triplets
```

```python
with open("duplicates/batch_summary.csv") as fh:
    dup = list(csv.DictReader(fh))
assert [(r["rows"], r["sample_id"]) for r in dup] == [
    ("2;3", "1UBQ_A_aa_ind_LYS"), ("4", "ubq_tyr_a"), ("5", "ubq_tyr_b")]
assert "one sample from 2 source requests" in dup[0]["message"]
```

### Residue identifiers and ranges

One grammar serves the command line, a manifest and Python. A residue is its
PDB number with an optional insertion code (`48`, `52A`, `-5`); `aa_selection`
lists also take ranges of plain numbers, ascending, both ends included
(`10-12`, `-5--3`). A range never includes an insertion-code residue: `9-11`
on SYNPHOS chain A, which has `10A` and `10B`, is refused rather than silently
leaving them out. List them to include them, or list the residues one by one
to leave them out:

```bash
biotsr run --method aa_selection --pdb SYNPHOS.pdb --chain A --residues "9-11;10A;10B" --output ranges/ --no-triplets
```

```python
from biotsr.domain.errors import InvalidRequest

listed = biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A",
                        residues="9-11;10A;10B", keep_triplets=False)
assert listed.selection.residues_used == ("9", "10", "10A", "10B", "11")
for refused, reason in (("9-11", "span residue(s) 10A, 10B"),
                        ("10A-11", "is not supported"), ("11-9", "is descending")):
    try:
        biotsr.analyze("SYNPHOS.pdb", method="aa_selection", chain="A",
                       residues=refused, keep_triplets=False)
    except InvalidRequest as exc:
        assert reason in str(exc)
    else:
        raise AssertionError(refused)
```

See the [CLI guide](../CLI_GUIDE.md), §7, for every column, the `--chain`
default, `--map` and `--sheet custom` for sheets in another layout,
downloading missing structures and running requests in parallel.
