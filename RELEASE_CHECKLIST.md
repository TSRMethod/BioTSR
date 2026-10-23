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
- [ ] `CITATION.cff` attribution confirmed: the agreed group attribution,
      "TSR Research Group"; individual authors, their order, affiliations or a
      preferred paper are added only if the maintainers choose to, never
      guessed *(human)*
- [ ] Commit author and committer metadata reviewed before any push that would
      publish it *(human)*
- [ ] If the repository is to be made public: approval of what it contains,
      including its history and its workflow-run records *(human)*

## B. Release candidate — built and evidenced

- [ ] Version set once, in `biotsr.__version__` *(automated: the packaging tests)*
- [ ] `ci.yml` green on the release commit: every `tests` matrix job (Linux
      Python 3.9–3.13, macOS arm64 Python 3.13), `gates` and `package`
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
- [ ] Release manifest regenerated for the release commit
      (`release_manifest.py`) *(automated)*
- [ ] PyPI name availability re-verified *(automated check)*

## C. Publication — only after every item in A is done

- [ ] Final version in `src/biotsr/__init__.py`, `CITATION.cff` (with
      `date-released`), `CHANGELOG.md` (dated) and every other version-bearing
      file *(human, then automated: `test_every_version_bearing_file_follows_the_version`)*
- [ ] Version-bearing gate evidence and the release manifest regenerated; clean
      tree *(automated)*
- [ ] Release tag on `main` *(human)*
- [ ] GitHub Release from `CHANGELOG.md`, with the wheel, sdist and manifest
      attached *(human)*
- [ ] Optional: PyPI upload — TestPyPI first, clean install from TestPyPI, then
      PyPI *(human)*
- [ ] Optional: Zenodo archive and software DOI; DOI added to `CITATION.cff` and
      the README *(human)*
- [ ] Website update (`tsrmethod.github.io`), only once the release is public
      *(human)*
