# 07 — Common keys and similarity

The analysis layer works on results from any method, as long as they share a
**key space**; comparing integers from two different key spaces is refused.

## Compute once

<!-- file: samples.csv -->
```csv
sample_id,pdb_id,method,chain,label
ubiquitin,1UBQ,protein,A,ubiquitin
barnase,1BRS,protein,A,barnase
barstar,1BRS,protein,D,barstar
```

```bash
biotsr batch --manifest samples.csv --pdb-dir pdb/ --output batch/ --no-triplets
```

## Common keys, from the command line

```bash
biotsr analyze common-keys --from-results batch/ --output analysis_out/
```

```text
loaded     3 saved result(s) from 3 calculation(s) in batch/
common keys 9523 over 3 samples
```

`--from-results` reads every completed calculation in the folder (its
`PREFIX.run.json` completion record), and each sample is named by its prefix.
`analysis_out/analysis/common_keys/all_samples.csv` reports, per sample, its
method, total keys, the common keys, their summed frequency and the sample's
key-space, processing and sampling fingerprints. The command always writes
into an `analysis/common_keys/` directory under `--output`.

## What has to match, and what never blocks

Every analysis — common keys, feature matrices, similarity, clustering, on
computed or loaded results — applies one rule:

1. **The key definitions must match, always.** Labels, angle and distance
   classes, encoding, ordering and the mirror convention: one key space. No
   option relaxes this. Different methods that share a key space compare
   (`aa_ind`, `drug`, `carbohydrate` and the other atomic methods share one);
   the same method with `--mirror` or `--aa-grouping` does not.
2. **The processing settings must match, unless you acknowledge the
   difference.** A size filter, a contact cutoff or a completeness rule
   changes how the triangles were drawn: the keys mean the same, but the
   frequencies answer a different question. Such a set is refused unless you
   pass `--allow-mixed-processing` (`allow_mixed_processing=True` in Python);
   the report then names each differing setting and records every sample's
   settings.
3. **A partial sample must be acknowledged, separately.** A sample whose
   selection did not get everything it asked for — a requested residue type
   absent from the chain, say — is refused unless you pass `--allow-partial`
   (`allow_partial=True`). This is not the same acknowledgement as
   `--allow-mixed-processing`: neither implies the other, and neither relaxes
   rule 1.
4. **What was selected never blocks.** Structures, chains, amino-acid types,
   ligands, sugars, residue numbers, atom subsets, one copy or all copies of a
   molecule: these are what a dataset is made of. Two samples need not even
   share a key; an empty intersection is a valid result.

A plain and a size-filtered run of the same chain differ in processing:

<!-- file: filtered.csv -->
```csv
sample_id,pdb_id,method,chain,size_filter
ubq_all,1UBQ,protein,A,
ubq_25,1UBQ,protein,A,25.0
```

```bash
biotsr batch --manifest filtered.csv --pdb-dir pdb/ --output filtered/ --no-triplets
biotsr analyze common-keys --from-results filtered/ --output filtered_out/ --allow-mixed-processing
```

```text
common keys 37522 over 2 samples
warning: the samples were not all processed with the same settings (2 processing policies; differing in: size_filter).
```

Without `--allow-mixed-processing` the same command is refused (exit 2),
naming `size_filter (25.0 vs None)`. See the [user guide](../USER_GUIDE.md),
§17.

## Similarity, from Python

```python
import biotsr
from biotsr.analysis import AnalysisSample, distance_matrix, feature_matrix

paths = {
    "ubiquitin": "batch/ubiquitin.keys",
    "barnase": "batch/barnase.keys",
    "barstar": "batch/barstar.keys",
}
samples = [AnalysisSample.from_result(biotsr.load_result(p), name)
           for name, p in paths.items()]

fm = feature_matrix(samples)       # one column per key, in numeric order
dm = distance_matrix(fm)           # generalized Jaccard distance
assert dm.measure == "generalized_jaccard"
assert dm.sample_ids == ("ubiquitin", "barnase", "barstar")
# distance between samples 1 and 2 (barnase and barstar)
assert round(dm.value(1, 2), 4) == 0.8674
assert fm.comparison["processing_uniform"] and fm.warnings == ()
```

Binary Jaccard, exported matrices, the dendrogram and the clustered heat map,
all from saved results: [tutorial 08](08_similarity_clustering.md).

## Different selections compare freely

Lysines of ubiquitin against tyrosines of barnase: different structures,
chains and residue types, the same method and the same processing. They
compare without any option:

```python
from biotsr.analysis import IncompatibleProcessing

lys = biotsr.analyze("1UBQ.pdb", method="aa_ind", chain="A",
                     residue_types="LYS", keep_triplets=False)
tyr = biotsr.analyze("1BRS.pdb", method="aa_ind", chain="A",
                     residue_types="TYR", keep_triplets=False)
types = feature_matrix([AnalysisSample.from_result(lys, "ubiquitin_LYS"),
                        AnalysisSample.from_result(tyr, "barnase_TYR")])
assert types.warnings == ()
assert types.comparison["processing_uniform"] is True
```

## What is refused

A processing difference is refused until it is acknowledged, and then
recorded:

```python
plain = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", keep_triplets=False)
short = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", size_filter=25.0,
                       keep_triplets=False)
pair = [AnalysisSample.from_result(plain, "ubq_all"),
        AnalysisSample.from_result(short, "ubq_25")]
try:
    feature_matrix(pair)
except IncompatibleProcessing as error:
    assert "size_filter" in str(error)
else:
    raise AssertionError("a processing difference must be acknowledged")

mixed = feature_matrix(pair, allow_mixed_processing=True)
assert mixed.comparison["processing_differences"] == ["size_filter"]
assert "size_filter" in mixed.warnings[0]
```

A **partial** sample is refused too, for a different reason: it is not the
selection that was requested. Chain B of the synthetic `SYNPHOS.pdb` has no
SER, so a TYR, LYS and SER selection there is partial, while the same
selection on chain A is complete. `--allow-partial` analyses such a set as
computed; the command then exits `4`
([tutorial 05](05_batch_manifest.md) computes this pair with
`biotsr batch`):

<!-- not-executed -->
```bash
biotsr analyze common-keys --from-results partial/ --output partial_ck/ --allow-partial
```

In Python the acknowledgement is `allow_partial=True`, and the refusal is
`biotsr.analysis.PartialSamples`. `allow_mixed_processing=True` does not
stand in for it:

```python
from biotsr.analysis import PartialSamples

chains = [AnalysisSample.from_result(
              biotsr.analyze("SYNPHOS.pdb", method="aa_ind_mix", chain=chain,
                             residue_types="TYR;LYS;SER", keep_triplets=False),
              "chain" + chain)
          for chain in ("A", "B")]
assert chains[1].omissions["primary"] == {"not_found": ["SER"]}
try:
    feature_matrix(chains, allow_mixed_processing=True)
except PartialSamples as error:
    assert "chainB (not found: SER)" in str(error)
else:
    raise AssertionError("a partial sample must be acknowledged")

acknowledged = feature_matrix(chains, allow_partial=True)
assert acknowledged.comparison["partial_allowed"] is True
assert list(acknowledged.comparison["partial_samples"]) == ["chainB"]
assert [s["partial"] for s in acknowledged.samples] == [False, True]
```

Mirror-signed and unsigned keys live in different key spaces, so mixing them is
an error rather than a silent zero overlap, and no option allows it:

```python
mirrored = biotsr.analyze("1UBQ.pdb", method="protein", chain="A", mirror=True)
try:
    feature_matrix(samples + [AnalysisSample.from_result(mirrored, "mirrored")],
                   allow_mixed_processing=True)
except biotsr.IncompatibleKeySpace as error:
    assert "mirror_convention" in str(error)
else:
    raise AssertionError("mixing key spaces must be refused")
```
