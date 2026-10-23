# Licensing and provenance

BioTSR is licensed under the **Apache License, Version 2.0** ([`LICENSE`](../LICENSE)).
The licence was applied on 2026-09-25, after the TSR Research Group completed an
ownership and provenance review of everything the repository contains. This page
records what that review established about the origin of BioTSR's contents;
[`NOTICE`](../NOTICE) carries the required third-party notices verbatim.

Third-party material is **not relicensed**: each item below keeps its own terms.

## How provenance was established

* every tracked file was searched for licence files, copyright headers, SPDX
  identifiers and licence names;
* the historical TSR repositories were inspected for a licence, their remote and
  their authors;
* every line of `src/biotsr` was compared, whitespace-insensitively, with every
  source file of the historical repositories, to separate copied code from
  reimplemented behaviour;
* every bundled data resource was hashed and matched against the historical
  repositories;
* third-party runtime and development dependencies were listed.

No vendored third-party source was found.

## What BioTSR contains, and where it comes from

| component | origin | terms |
|---|---|---|
| `src/biotsr/**` code | written for BioTSR | Apache-2.0 |
| theta and distance bin tables (`binning/theta.py`, `binning/distance.py`) and the amino-acid label and group tables (`labels/schemes.py`) | values from `TSR-Package` | MIT (`NOTICE` item 1) |
| `labels/resources/element_labels_v1.csv` | the element codes of the `CarbohydratesTSR` atom lexicon (ADR-113) | MIT, as derived from it (`NOTICE` item 2) |
| secondary-structure category rules (`sse/`) | reimplemented from `TSR-Package`'s `SSE_TSR.py` | MIT (`NOTICE` item 1) |
| atom-count tables (`selection/completeness.py`) and residue sets (`cross/sides.py`) | transcribed from the historical amino-acid, nucleotide and carbohydrate repositories | attribution (`NOTICE` item 4) |
| `analysis/clustering.py` | follows the algorithms of SciPy's `nn_chain`, `maxclust` and `dendrogram` (agreement and its scope: `docs/methods/ANALYSIS.md`, §5.7); no SciPy code is imported or bundled | BSD-3-Clause notice (`NOTICE` item 3) |
| `analysis/similarity.py` | independent implementation of published formulas | Apache-2.0 |
| `tests/data/structures/*.pdb` | unmodified RCSB PDB entries ([manifest](../tests/data/structures/manifest.json)) | wwPDB data policy |
| atom records in `tests/unit/test_element_labels.py` | single records copied verbatim from RCSB entries 1OL6, 1EVZ, 1UFJ, 4X6J, 3FSF, 1RER, 1ZHM, 2AXT, 1BXQ | wwPDB data policy |
| `tests/data/synthetic/` | a synthetic structure and the script that writes it ([README](../tests/data/synthetic/README.md)); no real coordinates | Apache-2.0 |
| `tests/data/legacy_outputs/` | unedited outputs of earlier BioTSR builds on the synthetic structure and on an excerpt of 1BNA, kept to test which are reused and which never are ([README](../tests/data/legacy_outputs/before_adr101_hetatm_modified_residues/README.md), [README](../tests/data/legacy_outputs/before_adr107_calpha_radix20/README.md), [README](../tests/data/legacy_outputs/before_adr113_adr114/README.md)) | Apache-2.0 |
| `tests/data/reference/standard_baseline.json`, `calpha_base30.json`, `atomic_element_labels.json` | BioTSR's own standard results on those structures | Apache-2.0 |
| `tests/data/reference/historical/` | one column (the triangle keys) of `TSR-Package` output on 1UBQ | MIT (`NOTICE` item 1) |
| gate evidence under `validation/evidence/` (repository only) | outputs of BioTSR's own gate scripts | not part of the installed package |

`TSR-Package` shares part of its key-generation code with the earlier
`hsp70_actin` repository, and BioTSR's bin and label tables reach that lineage
through `TSR-Package`. BioTSR contains no `hsp70_actin` code; the lineage is
acknowledged in `NOTICE` item 5.

The atom-count tables and residue sets transcribed from the historical
repositories that carry no licence text of their own are attributed in
`NOTICE` item 4.

**Redistribution.** The project owner has confirmed proceeding with the
redistribution of these tables as part of BioTSR — the amino-acid and
nucleotide atom-count tables and the nucleotide residue vocabulary, the bin and
label tables of the `TSR-Package` lineage, and the standard atom lexicon, whose
element codes the packaged element table carries (ADR-100, ADR-113). That confirmation is the project owner's; it is not a permission
document from a third party, and none is claimed. The attributions and licence
notices above and in `NOTICE` are kept unchanged.

## What each distribution contains

| distribution | contents |
|---|---|
| wheel | `src/biotsr/**`, its element label CSV, `LICENSE`, `NOTICE` |
| source distribution | the wheel's contents plus `tests/` (with `tests/data/structures/`, `tests/data/synthetic/` and `tests/data/legacy_outputs/`), `README.md`, `CHANGELOG.md`, `CITATION.cff`, `pyproject.toml` |
| source repository | everything above, plus the documentation and the release and gate scripts and evidence under `validation/` |

BioTSR has one runtime dependency, joblib (BSD-3-Clause), which pip installs
alongside it. It, the optional `plot` extra (matplotlib) and the
development tools (pytest, build, twine) are installed under their own licences
and never bundled.
