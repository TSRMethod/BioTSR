# Releasing BioTSR

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology. This page is the release engineering contract:
how versions are numbered, which Pythons are supported, what CI runs, what ships
in the package, and the steps to a release. The tick-list is
[`RELEASE_CHECKLIST.md`](../RELEASE_CHECKLIST.md).

## 1. Versioning

Versions follow [PEP 440](https://peps.python.org/pep-0440/). The current
version is `1.0.0rc7`, a release candidate; `1.0.0` has not been released.

* **One source of truth:** `__version__` in `src/biotsr/__init__.py`.
  `pyproject.toml` declares `dynamic = ["version"]` and reads it statically
  (`[tool.setuptools.dynamic] version = { attr = "biotsr.__version__" }`), so the
  wheel and sdist metadata, `biotsr --version` and the completion record's
  (`PREFIX.run.json`) `biotsr_version`
  agree by construction. A test fails if a second, static version appears.
* **The package version is not the result-schema version.** `RESULT_SCHEMA_VERSION`
  (currently `2`) versions the native file format and changes only when that
  format does (ADR-074, ADR-078). Schema 2 (1.0.0rc4) is the one-prefix layout
  with the per-instance selection record; files of schema 1, written by
  1.0.0rc3 and earlier release candidates, are refused rather than read.
  Neither version is part of any scientific fingerprint.
* **Release candidates** are `1.0.0rcN`. A change after a candidate has been
  reviewed produces the next candidate; only the version string changes between
  the candidate that is released and `1.0.0`. The tag is `v1.0.0`, on `main`,
  created by a person.
* **After 1.0.0:** a change to what an existing key *means* (any `KeySpace`
  field or default), to a standard profile's science, or a breaking change to
  the native format or the public API is a **major** release; new methods or
  options that leave existing results bit-identical are **minor**; fixes that
  change no validated output are **patch** releases.

## 2. Python and platform support

**Supported: CPython 3.9, 3.10, 3.11, 3.12 and 3.13, on Linux and macOS.**

A version is listed in `pyproject.toml`'s classifiers only when CI runs the full
default test suite on it. The floor (`requires-python = ">=3.9"`) is the oldest
version so tested, and a test asserts that the classifiers, the CI matrix and the
floor agree.

Integer keys, frequencies, bins, categories, selection and candidate counts and
both fingerprints are identical on Linux and macOS. The optional native triplet
file prints raw floating-point geometry computed through the host C maths
library, and the test suite compares it byte for byte only against references
produced on macOS arm64: [`PLATFORM_VALIDATION.md`](PLATFORM_VALIDATION.md).

Python 3.9 reached its upstream end of life in October 2025. It remains
supported for 1.0 because it is tested and costs nothing (the calculation uses
the standard library only, and the one dependency, joblib, supports 3.9), and
because HPC environments often lag; dropping it is a minor-release decision.

**Windows is not supported**: BioTSR is neither tested nor validated there
([known limitations](KNOWN_LIMITATIONS.md)).

**The `plot` extra** (`matplotlib>=3.5`) was checked locally with matplotlib
3.5.3 on Python 3.9 and 3.11.2 on Python 3.13; CI runs it with the newest
matplotlib on Python 3.13. Nothing else in BioTSR imports a plotting library.

## 3. CI

Two workflows under `.github/workflows/`.

### On pushes to `main` and release branches, and on pull requests — `ci.yml`

| job | what |
|---|---|
| `tests` | installs `.[dev]` and runs the default test suite (`pytest -q -rs`) on Python 3.9–3.13 on Ubuntu, and on Python 3.13 on macOS arm64, and once more with `.[dev,plot]` on Python 3.13 on Ubuntu so the figure tests and the tutorial's drawing steps run, with `BIOTSR_REQUIRE_EVIDENCE=1`, so that tracked evidence that is missing fails instead of skipping; then requires `git status` to be empty. Standard results are compared against `tests/data/reference/`: keys files by SHA-256 on every platform, triplet files by SHA-256 on macOS arm64 |
| `gates` | on Python 3.13, Ubuntu: `capability_matrix.py --check` (the capability matrix matches the method registry), `sanitize_tree.py --check` (no workstation path or account identifier in the tree), `gates_cli.py` and `gates_api_hardening.py`; then requires `git status` to be empty, so committed gate evidence regenerates byte-identically |
| `package` | on Python 3.12, Ubuntu: `release_artifacts.py` builds the wheel and sdist, inspects their contents and metadata, runs the unpacked sdist's own test suite, runs `twine check`, and installs each artifact into a clean environment outside the checkout for the clean-environment checks; the built artifacts are uploaded |

`ci.yml` can also be started by hand (`workflow_dispatch`). Nothing in it needs
history, tags, earlier artifacts, network access to the PDB archive, or files
outside the repository.

### On demand and for `v*` tags — `release-validation.yml`

Once, on Python 3.13 on macOS arm64, with the `plot` extra installed: the
**entire** test suite including the slow reference cases over the large
structures (`pytest -q -rs -m ""`, with
`BIOTSR_REQUIRE_EVIDENCE=1`), the capability-matrix and sanitizer checks, both
gate scripts and the clean-tree check, the artifact build and clean installs
(`release_artifacts.py`), and the release manifest (`release_manifest.py`). The
artifacts and the manifest are uploaded as a workflow artifact.

### Validation scripts

The tracked scripts under `validation/scripts/`:

| script | purpose | run by |
|---|---|---|
| `capability_matrix.py` | regenerates `docs/CAPABILITY_MATRIX.md` from the method registry; `--check` fails if it is stale | `ci.yml`, release validation |
| `sanitize_tree.py` | finds workstation paths and account identifiers in the tree; `--check` fails if any remain | `ci.yml`, release validation |
| `gates_cli.py` | command-line gates: CLI/API equivalence, the output contract, discovery, safety; writes `validation/evidence/cli/GATE_RESULTS.json` | `ci.yml`, release validation |
| `gates_api_hardening.py` | save/load round trip, reader refusals, analysis on loaded results, the public surface; writes `validation/evidence/api_hardening/GATE_RESULTS.json` | `ci.yml`, release validation |
| `release_artifacts.py` | builds, inspects and clean-installs the wheel and sdist | `ci.yml`, release validation |
| `clean_install_check.py` | the clean-environment checks `release_artifacts.py` runs against each installed artifact | via `release_artifacts.py` |
| `release_manifest.py` | writes the release manifest (§6) | release validation |
| `analysis_scipy_equivalence.py` | the scipy-equivalence evidence behind ADR-060 (`validation/evidence/analysis/SCIPY_EQUIVALENCE_scipy-<version>.json`, one per scipy release checked: 1.11.1 and 1.18.1) | by hand |
| `performance_api.py`, `performance_cli.py`, `performance_sse.py` | run-time and memory measurements; snapshots of one machine | by hand |

### Skips

A test that cannot run skips with a reason of one of two kinds, defined in
`tests/conftest.py`:

* `repository-only evidence not present` — tracked material outside the
  package (for example `docs/` or `validation/`) is missing. That is expected in
  an unpacked sdist, which carries `tests/` and `tests/data/` but not `docs/` or
  `validation/`. In the repository it is a defect, and
  `BIOTSR_REQUIRE_EVIDENCE=1` (set in CI) turns it into a failure.
* `optional dependency not installed` — a third-party package some tests use
  when it is installed. joblib is a runtime dependency, so the
  parallel-execution tests fail rather than skip without it.

`release_artifacts.py` requires every skip in the unpacked sdist's test run to
be of an expected kind, and requires a minimum number of tests to pass there, so
the sdist cannot pass by skipping.

## 4. Deterministic evidence

Re-running a gate on identical inputs must leave the tree clean, on any machine
and from any checkout path (ADR-081):

* committed gate reports carry results and **no timing**;
* committed sidecars record input paths **relative to the repository**; an
  input's identity is its `source.sha256`, and no fingerprint includes a path;
* the version string is provenance: `validation/evidence/api_hardening/GATE_RESULTS.json`
  records it, so bumping the version changes that report once, and the change is
  committed with the bump.

## 5. What ships where

The repository and the package are different things.

| | repository | sdist | wheel |
|---|---|---|---|
| `src/biotsr/` and the element label table | yes | yes | yes |
| `tests/`, with `tests/data/` (structures, reference results, manifests) | yes | yes (tests that need `docs/` or `validation/` skip, saying so) | no |
| `README.md`, `CHANGELOG.md`, `CITATION.cff`, `pyproject.toml` | yes | yes | README as the long description |
| `LICENSE` (Apache-2.0) and `NOTICE` (third-party notices) | yes | yes | yes |
| `docs/` | yes | no | no |
| `validation/` — validation scripts and committed evidence | yes | no | no |

`MANIFEST.in` states the sdist's contents; `release_artifacts.py` fails if an
artifact carries anything else, or any absolute workstation path.

## 6. Release manifest

`validation/scripts/release_manifest.py` writes one JSON document connecting a
commit to its artifacts: version, commit, a content hash of `src/biotsr`, Python
and platform, tool versions, the packaged label table's hash and content fingerprint,
the standard `KeySpace` fingerprint of every method, reference-run fingerprints,
gate counts, the pytest summary when a log is supplied, and the SHA-256 of every
built artifact. Release validation writes it to `dist/` and uploads it with the
wheel and sdist. A copy is committed under `validation/evidence/release/` only as
the release step in §8 describes, so that a committed record always names a
commit that exists in this repository.

## 7. Version-bearing files

The files that name the package version or its pre-release status, and so change
with the release bump, are listed as `VERSION_BEARING` in
`tests/unit/test_release_readiness.py`:

* `src/biotsr/__init__.py`, `CITATION.cff`, `CHANGELOG.md`, `README.md`;
* `docs/USER_GUIDE.md`, `docs/CLI_GUIDE.md`, `docs/API_GUIDE.md`,
  `docs/PACKAGE_OVERVIEW.md`, `docs/REPRODUCIBILITY_NOTE.md`.

While the version is a release candidate, each must name the version or its
status; once it is final, `test_every_version_bearing_file_follows_the_version`
fails on any pre-release wording left behind. The committed evidence that
records the version (`validation/evidence/api_hardening/GATE_RESULTS.json`, and
the release records of §8 step 3) is regenerated with the bump.

## 8. Producing 1.0.0

Nothing is published until every decision in `RELEASE_CHECKLIST.md` §A has been
made by the people with the authority to make it. Then:

1. Merge the reviewed release-candidate branch into `main`.
2. Set `__version__ = "1.0.0"`; set `version: "1.0.0"` and `date-released` in
   `CITATION.cff`; date the `CHANGELOG.md` entry as `## [1.0.0] — YYYY-MM-DD`
   and point its link at the tag; and remove the pre-release wording from every
   other version-bearing file (§7).
3. Run `gates_api_hardening.py` and `gates_cli.py`, regenerate the release
   evidence (`release_artifacts.py --json validation/evidence/release/RELEASE_ARTIFACTS.json`
   and `release_manifest.py --output validation/evidence/release/RELEASE_CANDIDATE_MANIFEST.json`),
   commit the regenerated files, and confirm `git status` is clean after a
   second gate run.
4. Push; wait for `ci.yml` to be green; run `release-validation.yml`.
5. Tag: `git tag -a v1.0.0 -m "BioTSR 1.0.0"` and push the tag (this runs
   release validation on the tag).
6. Create the GitHub Release from the tag, with the `CHANGELOG.md` entry as
   notes and the wheel, sdist and manifest from the tag's release-validation run
   attached.
7. Optional PyPI: `twine upload --repository testpypi dist/*`, install from
   TestPyPI in a clean environment and run the clean-environment checks, then
   `twine upload dist/*`.
8. Optional Zenodo: archive the release, then add the DOI to `CITATION.cff` and
   the README in a follow-up commit.
9. Update the group website once the release is public.
