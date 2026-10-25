# BioTSR release checklist

The steps from a release candidate to a published release. The procedure behind
each step is in [docs/RELEASING.md](docs/RELEASING.md). Copy this list into the
release's tracking issue and tick items there; this file itself stays unticked.

Two kinds of item:

* **automated** — a script, test or CI job decides it; tick it when the named
  check is green on the commit being released;
* **human** — a person with the authority decides it; no script can.

## A. Decisions — nothing is published until these are made

- [ ] The maintainers authorize publication: final version, tag, GitHub
      Release, and any PyPI or Zenodo upload *(human)*
- [ ] Licence and provenance review is current for everything the release
      ships; the project owner's confirmation that the historical tables
      listed in [docs/LICENSING.md](docs/LICENSING.md) and [`NOTICE`](NOTICE)
      may be redistributed on the stated terms is recorded (ADR-100), and
      material added since is reviewed the same way *(human)*
- [ ] No decision record the release depends on is still proposed: every one
      in [docs/architecture/DECISIONS.md](docs/architecture/DECISIONS.md) is
      accepted, superseded or historical *(human)*
- [ ] `CITATION.cff` and the package metadata confirmed: the software's
      authors and their order as the maintainers set them; ORCIDs,
      affiliations, a DOI or a preferred paper are added only when the
      maintainers supply them, never guessed *(human)*
- [ ] Commit author and committer metadata reviewed before any push that would
      publish it *(human)*
- [ ] If the repository is to be made public: approval of what it contains,
      including its history and its workflow-run records *(human)*

## B. Release candidate — built and evidenced

- [ ] Version set once, in `biotsr.__version__` *(automated: the packaging tests)*
- [ ] `ci.yml` green on the release commit: every `tests` matrix job (Linux
      and macOS arm64, Python 3.9–3.13 on each), `gates` and `package`
      *(automated)*
- [ ] `release-validation.yml` green on the release commit (macOS arm64: the
      full suite including slow reference cases, the gates, the artifacts and
      the manifest) *(automated)*
- [ ] README and user-guide examples executed by the test suite
      (`tests/integration/test_tutorials.py`) and documentation links resolve
      (`tests/unit/test_docs_links.py`) *(automated)*
- [ ] Capability matrix matches the registry (`capability_matrix.py --check`)
      *(automated)*
- [ ] Standard reference results unchanged
      (`tests/integration/test_standard_references.py`) *(automated)*
- [ ] Wheel and sdist built, inspected, and each installed cleanly; the unpacked
      sdist's own test suite passes (`release_artifacts.py`) *(automated)*
- [ ] Gate evidence regenerates with a clean tree (`gates_cli.py`,
      `gates_api_hardening.py`) *(automated)*
- [ ] HEAD free of workstation paths and account identifiers
      (`sanitize_tree.py --check`) *(automated)*
- [ ] Release manifest generated for the release commit by release
      validation (`release_manifest.py`) *(automated)*
- [ ] PyPI name availability re-verified *(automated check)*
- [ ] No workflow but `publish.yml` can publish, and it only by hand, from a
      release tag, for reviewed hashes, after an environment approval
      (`test_release_readiness.py`) *(automated)*

## C. Publication — only after every item in A is done

The procedure is [docs/RELEASING.md](docs/RELEASING.md), §8–§10.

- [ ] Final version in `src/biotsr/__init__.py`, `CITATION.cff` (with
      `date-released` once the date is fixed), `CHANGELOG.md` (dated) and every
      other version-bearing file *(human, then automated:
      `test_every_version_bearing_file_follows_the_version`)*
- [ ] Version-bearing gate evidence regenerated; clean tree *(automated)*
- [ ] `ci.yml` and `release-validation.yml` green on the release branch's
      head, and `main` fast-forwarded to that same commit *(human)*
- [ ] Release tag `vX.Y.Z` on that commit; release validation green on the tag;
      the SHA-256 of its wheel and sdist recorded *(human)*
- [ ] GitHub Release from `CHANGELOG.md`, with the tag's wheel, sdist,
      manifest and artifact report attached *(human)*
- [ ] One-time PyPI setup: `testpypi` and `pypi` environments with required
      reviewers, and pending publishers on TestPyPI and PyPI *(human)*
- [ ] Optional: PyPI upload with `publish.yml` — TestPyPI first, its hashes and
      a clean install checked, then PyPI *(human)*
- [ ] Optional: Zenodo archive and software DOI; DOI added to `CITATION.cff` and
      the README *(human)*
- [ ] Website update (`tsrmethod.github.io`), only once the release is public
      *(human)*
