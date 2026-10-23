# Validation material

The checks BioTSR runs on itself as a released artefact, the scripts that run
them, and the committed records they write. None of it ships in the wheel or the
sdist; the scientific tests that ship with the package live in `tests/` and run
on the structures under `tests/data/`.

Every method runs under its standard profile and writes the native result
format, so everything here checks standard behaviour only.

## `evidence/`

| path | what it records | written by |
|---|---|---|
| `cli/GATE_RESULTS.json` | the command-line contract: the CLI computes exactly what the Python API computes, keeps its promises about output prefixes, completion records, exit codes and the shared output folder, describes the methods from the registry, and refuses what it must | `scripts/gates_cli.py` |
| `api_hardening/GATE_RESULTS.json` | the saved-result guarantees: round trips for every representative family, refusal of corrupt or mismatched files, analyses on loaded results equal to analyses on computed ones, the public API surface, and filesystem robustness | `scripts/gates_api_hardening.py` |
| `cli/PERFORMANCE.json`, `api_hardening/PERFORMANCE.json` | run time and memory of one machine; written when measured, not committed | `scripts/performance_cli.py`, `scripts/performance_api.py` |
| `analysis/SCIPY_EQUIVALENCE_scipy-<version>.json` | the standard-library clustering compared with one scipy and scikit-learn release over randomised and tie-heavy matrices (1.11.1 and 1.18.1); the scope of the agreement is in `docs/methods/ANALYSIS.md`, §5.7 | `scripts/analysis_scipy_equivalence.py` |
| `release/RELEASE_ARTIFACTS.json` | the built wheel and sdist: contents, metadata, the sdist's own test run, and a clean install of each; written at release time, not committed | `scripts/release_artifacts.py` |
| `release/RELEASE_CANDIDATE_MANIFEST.json` | one release candidate: version, commit, label-table and key-space fingerprints, reference runs on 1UBQ, gate counts, artifact hashes; written at release time, not committed | `scripts/release_manifest.py` |

The two `GATE_RESULTS.json` files are deterministic: they record no timing and
no machine path, so re-running a gate script on an unchanged tree leaves them
byte-identical. The performance, release and SciPy records describe one
machine and one run, and are refreshed only when a release is prepared.

## `scripts/`

| script | purpose |
|---|---|
| `gates_cli.py` | command-line gates (equivalence with the API, contract, discovery, safety); writes `evidence/cli/GATE_RESULTS.json`, or `--json PATH` |
| `gates_api_hardening.py` | saved-result and API-surface gates (round trip, refusal, analysis, surface, robustness); writes `evidence/api_hardening/GATE_RESULTS.json`, or `--json PATH` |
| `capability_matrix.py` | generates `docs/CAPABILITY_MATRIX.md` from the method registry (family, modifiers, target arguments, standard profile); `--check` fails if the page is stale |
| `sanitize_tree.py` | removes workstation paths and account identifiers from tracked evidence; `--check` fails if anything would change |
| `release_artifacts.py` | builds the wheel and the sdist, inspects both (required and forbidden members, the packaged label table, metadata), runs the sdist's test suite, and installs each artifact into a clean environment |
| `clean_install_check.py` | the clean-environment checks `release_artifacts.py` runs against each artifact: console script, `python -m biotsr`, real workflows from outside the repository, no removed legacy module installed, `joblib` importable |
| `release_manifest.py` | writes the release-candidate manifest |
| `analysis_scipy_equivalence.py` | the clustering equivalence check; needs numpy, scipy and scikit-learn, which BioTSR does not depend on, and skips without them |
| `performance_cli.py`, `performance_api.py`, `performance_sse.py` | measure run time and memory; snapshots of one machine |

Every script runs from the repository root with the package importable, for
example:

```bash
python validation/scripts/gates_cli.py
python validation/scripts/gates_api_hardening.py
python validation/scripts/capability_matrix.py --check
python validation/scripts/sanitize_tree.py --check
python validation/scripts/release_artifacts.py --dist-dir dist --json dist/RELEASE_ARTIFACTS.json
```

`release_artifacts.py` needs `build` and `pytest` in the running interpreter;
it reads only the checked-out tree, plus the time of the `HEAD` commit for a
reproducible wheel. `docs/RELEASING.md` says which of these CI runs, and when.
