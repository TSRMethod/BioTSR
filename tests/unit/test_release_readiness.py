"""Release gates: the durable promises of BioTSR 1.0.

Each test is one promise the release makes, so that breaking it later fails the
suite rather than surprising a user.  The categories:

* **branding** -- one name everywhere, the canonical expansion of TSR;
* **packaging** -- one version source, honest licence metadata, declared
  resources, an sdist that ships what it should and nothing it should not;
* **documentation** -- the capability matrix is the registry's, the release
  documents exist and say what they must;
* **scientific invariance** -- fingerprints and label-table bytes are those validated
  before the release candidate;
* **metal ions** -- no admission that reads HETATM records respells an ion;
* **CI and determinism** -- the workflows exist, test what the classifiers claim,
  and the committed gate evidence carries no timing and no workstation path;
* **public tree** -- no identifying material left in HEAD.

Building and installing the actual wheel and sdist is too slow for a unit test
and needs the ``build`` tool; it is ``validation/scripts/release_artifacts.py``,
which CI runs on every push.
"""

from __future__ import annotations

import ast
import glob
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys

import pytest

import biotsr
from biotsr.methods.registry import get_method, list_methods

from ..conftest import (
    FIXTURE_PDB,
    OPTIONAL_DEPENDENCY,
    REPOSITORY_ONLY,
    SRC,
    WORKSPACE,
    _require,
)

PYPROJECT = os.path.join(WORKSPACE, "pyproject.toml")
SCRIPTS = os.path.join(WORKSPACE, "validation", "scripts")
EXPANSION = "Triangular Spatial Relationship"
CANONICAL = ("BioTSR is the unified software implementation of the Triangular Spatial "
             "Relationship (TSR) methodology.")
#: sha256 of https://www.apache.org/licenses/LICENSE-2.0.txt, unmodified.
APACHE_2_SHA256 = "cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30"


def _read(*parts):
    path = os.path.join(WORKSPACE, *parts)
    _require(path)
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _script(name):
    """Import a validation script as a module (they are not a package)."""
    path = os.path.join(SCRIPTS, name + ".py")
    _require(path)
    spec = importlib.util.spec_from_file_location("_release_script_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _one_line(text):
    return " ".join(text.split())


# ---------------------------------------------------------------------------
# Branding
# ---------------------------------------------------------------------------


def test_the_readme_opens_with_the_canonical_sentence():
    assert CANONICAL in _one_line(_read("README.md"))


@pytest.mark.parametrize("path", ["README.md", "pyproject.toml", "CITATION.cff",
                                  "CHANGELOG.md", "docs/CAPABILITY_MATRIX.md"])
def test_release_documents_expand_tsr_canonically(path):
    assert EXPANSION in _one_line(_read(*path.split("/")))


def test_no_public_surface_uses_a_retired_name():
    pattern = re.compile(r"tsr[-_ ]unified|unified[-_ ]tsr", re.I)
    offenders = []
    candidates = ["README.md", "pyproject.toml", "CITATION.cff", "CHANGELOG.md",
                  "RELEASE_CHECKLIST.md", "MANIFEST.in"]
    candidates += [os.path.relpath(p, WORKSPACE) for p in glob.glob(
        os.path.join(WORKSPACE, "docs", "*.md")) + glob.glob(
        os.path.join(WORKSPACE, "docs", "tutorials", "*.md"))]
    for root, _dirs, files in os.walk(os.path.join(SRC, "biotsr")):
        candidates += [os.path.relpath(os.path.join(root, f), WORKSPACE)
                       for f in files if f.endswith(".py")]
    for rel in candidates:
        path = os.path.join(WORKSPACE, rel)
        if os.path.exists(path) and pattern.search(open(path, encoding="utf-8").read()):
            offenders.append(rel)
    assert offenders == []


def test_distribution_import_and_command_are_all_biotsr():
    text = _read("pyproject.toml")
    assert re.search(r'^name = "biotsr"$', text, re.M)
    assert re.search(r'^biotsr = "biotsr\.cli:main"$', text, re.M)
    assert biotsr.__name__ == "biotsr"
    assert os.path.exists(os.path.join(SRC, "biotsr", "__main__.py"))


def test_no_current_document_points_at_the_retired_future_layout():
    """The package lives at the repository root; ``future/`` survives only in history."""
    offenders = []
    for rel in ["README.md", "CHANGELOG.md", "RELEASE_CHECKLIST.md", "pyproject.toml",
                "docs/API_GUIDE.md", "docs/CLI_GUIDE.md", "docs/USER_GUIDE.md",
                "docs/REPRODUCIBILITY_NOTE.md",
                "docs/PACKAGE_OVERVIEW.md", "docs/KNOWN_LIMITATIONS.md", "docs/RELEASING.md",
                "docs/CAPABILITY_MATRIX.md", ".github/workflows/ci.yml",
                ".github/workflows/release-validation.yml",
                ".github/workflows/publish.yml"]:
        path = os.path.join(WORKSPACE, rel)
        if os.path.exists(path) and re.search(r"future/(biotsr|tsr)|pip install future",
                                              open(path, encoding="utf-8").read()):
            offenders.append(rel)
    assert offenders == []
    assert not os.path.exists(os.path.join(WORKSPACE, "future", "biotsr", "pyproject.toml"))


# ---------------------------------------------------------------------------
# Packaging
# ---------------------------------------------------------------------------


def test_the_version_has_exactly_one_source():
    text = _read("pyproject.toml")
    project = text.split("[project]", 1)[1].split("\n[", 1)[0]
    assert not re.search(r"^version\s*=", project, re.M), "a second, static version"
    assert re.search(r'^dynamic = \["version"\]$', project, re.M)
    assert re.search(r'^version = \{ attr = "biotsr\.__version__" \}$', text, re.M)


def test_the_version_is_a_pep440_release_candidate_of_1_0_0():
    assert re.fullmatch(r"1\.0\.0(rc[1-9]\d*)?", biotsr.__version__), biotsr.__version__


def test_the_cli_reports_the_package_version():
    proc = subprocess.run([sys.executable, "-m", "biotsr", "--version"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True,
                          env=dict(os.environ, PYTHONPATH=SRC))
    assert proc.returncode == 0
    assert proc.stdout.strip() == "biotsr %s" % biotsr.__version__


def test_the_licence_metadata_matches_the_licence():
    """Apache-2.0: the canonical text in LICENSE, and the SPDX expression (PEP 639)
    in the metadata, with LICENSE and NOTICE declared as licence files. The
    deprecated licence classifier is not used alongside it."""
    text = _read("pyproject.toml")
    with open(_require(os.path.join(WORKSPACE, "LICENSE")), "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == APACHE_2_SHA256
    assert re.search(r'^license = "Apache-2\.0"$', text, re.M)
    assert re.search(r'^license-files = \["LICENSE", "NOTICE"\]$', text, re.M)
    assert "License :: " not in text
    assert re.search(r'requires = \["setuptools>=77', text)


def test_notice_carries_every_third_party_notice_the_audit_requires():
    """ADR-085: third-party material keeps its own terms, and its notices ship.

    TSR-Package and CarbohydratesTSR (MIT) and SciPy (BSD-3-Clause), verbatim,
    plus attribution for the historical repositories BioTSR's atom-count tables
    and residue sets come from (docs/LICENSING.md).
    """
    _require(os.path.join(WORKSPACE, "LICENSE"))
    text = _read("NOTICE")
    for needle in ("Copyright 2026 TSR Research Group",
                   "Copyright (c) 2024 Poorya Khajouie",
                   "Copyright (c) 2026 Krishna Rauniyar",
                   "Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.",
                   "element_labels_v1.csv", "selection/completeness.py",
                   "analysis/clustering.py", "not relicensed"):
        assert needle in text, needle
    assert text.count("Permission is hereby granted, free of charge") == 2
    assert text.count("Redistribution and use in source and binary forms") == 1
    assert 'license: Apache-2.0' in _read("CITATION.cff")


def test_joblib_is_the_only_runtime_dependency():
    """ADR-086: one runtime dependency, for parallel batches; the old extra stays
    as an alias so instructions that name it still install."""
    text = _read("pyproject.toml")
    assert '\ndependencies = ["joblib>=1.3"]\n' in text
    assert '\nparallel = ["joblib>=1.3"]\n' in text


def test_every_packaged_label_table_is_declared_and_present():
    text = _read("pyproject.toml")
    # Named exactly, so a stray table left beside it is never packaged.
    assert 'biotsr = ["labels/resources/element_labels_v1.csv"]' in text
    resources = os.path.join(SRC, "biotsr", "labels", "resources")
    assert sorted(f for f in os.listdir(resources) if f.endswith(".csv")) == [
        "element_labels_v1.csv"]


def test_the_sdist_manifest_keeps_the_repository_out():
    text = _read("MANIFEST.in")
    for line in ("prune validation", "prune legacy", "prune docs", "graft src/biotsr",
                 "graft tests"):
        assert line in text.splitlines(), line


def test_urls_are_the_canonical_ones():
    text = _read("pyproject.toml")
    assert 'Repository = "https://github.com/TSRMethod/BioTSR"' in text
    assert 'Homepage = "https://tsrmethod.github.io/"' in text


# ---------------------------------------------------------------------------
# Documentation
# ---------------------------------------------------------------------------


def test_the_capability_matrix_is_the_registry():
    module = _script("capability_matrix")
    assert module.render() == _read("docs", "CAPABILITY_MATRIX.md")


def test_the_capability_matrix_never_claims_what_a_method_refuses():
    for name in list_methods():
        method = get_method(name)
        for modifier in ("aa_grouping", "mirror", "size_filter", "sse"):
            supported = modifier in method.supported_modifiers
            try:
                method.validate_modifiers({modifier})
                accepted = True
            except ValueError:
                accepted = False
            assert supported == accepted, (name, modifier)


def test_the_citation_file_states_only_verified_metadata():
    text = _read("CITATION.cff")
    assert "cff-version: 1.2.0" in text
    assert 'title: "BioTSR"' in text
    assert 'version: "%s"' % biotsr.__version__ in text
    assert 'repository-code: "https://github.com/TSRMethod/BioTSR"' in text
    # Nothing unverified is invented: no DOI and no licence until they exist.
    assert not re.search(r"^\s*doi:", text, re.M)
    # An unreleased candidate has no release date.
    if "rc" in biotsr.__version__:
        assert not re.search(r"^date-released:", text, re.M)
    if not os.path.exists(os.path.join(WORKSPACE, "LICENSE")):
        assert not re.search(r"^license:", text, re.M)
    # The citation and the package metadata name the same authors, in the same
    # order, and the same licence.  A citation author is a person (given and
    # family names) or an entity (name).
    project = _read("pyproject.toml")
    authors = re.findall(r'\{ name = "([^"]+)" \}', project.split("authors = [", 1)[1]
                         .split("]", 1)[0])
    block = text.split("\nauthors:\n", 1)[1]
    block = block[:re.search(r"^\S", block, re.M).start()]
    cited = []
    for entry in re.split(r"^  - ", block, flags=re.M)[1:]:
        fields = dict(re.findall(r'^\s*([a-z-]+): "([^"]*)"$', entry, re.M))
        cited.append(fields.get("name")
                     or "%s %s" % (fields["given-names"], fields["family-names"]))
    assert authors and cited == authors
    assert re.search(r'^license = "Apache-2\.0"$', project, re.M)
    assert re.search(r"^license: Apache-2\.0$", text, re.M)


def test_the_changelog_has_an_entry_for_this_version():
    assert "## [%s]" % biotsr.__version__ in _read("CHANGELOG.md")


def test_the_release_checklist_claims_no_approval():
    """The public checklist is a template: every human approval stays unticked.

    Whether an approval has been given is recorded where the people who give it
    record it, never asserted by a file in the repository.
    """
    text = _read("RELEASE_CHECKLIST.md")
    assert "- [x]" not in text.lower()
    approvals = text.split("## A.", 1)[1].split("## B.", 1)[0]
    assert approvals.count("- [ ]") >= 4
    for item in ("authorize publication", "Licence and provenance", "CITATION.cff",
                 "proposed", "PyPI", "Zenodo"):
        assert item in text, item


# ---------------------------------------------------------------------------
# Scientific invariance -- the values validated before the release candidate
# ---------------------------------------------------------------------------

PINNED_LABEL_TABLE_SHA256 = {
    "element_labels_v1.csv": "d7d8e489ffbb5dab",
}
PINNED_LABEL_FINGERPRINTS = {
    "ELEMENT_LABELS_112": "17bee1b08fe364249b67d7821b0dc8d46fbdc3c1b68de4e77d5482746c024736",
}


@pytest.mark.parametrize("name", sorted(PINNED_LABEL_TABLE_SHA256))
def test_packaged_label_table_bytes_are_unchanged(name):
    path = os.path.join(SRC, "biotsr", "labels", "resources", name)
    with open(path, "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest().startswith(PINNED_LABEL_TABLE_SHA256[name])


@pytest.mark.parametrize("name", sorted(PINNED_LABEL_FINGERPRINTS))
def test_label_content_fingerprints_are_unchanged(name):
    from biotsr.labels import elements
    assert getattr(elements, name).content_fingerprint() == PINNED_LABEL_FINGERPRINTS[name]


def test_reference_fingerprints_are_unchanged(ubq_pdb):
    protein = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    atomic = biotsr.analyze(ubq_pdb, method="amino_acid_atomic", chain="A",
                            residue="15", residue_name="LEU", keep_triplets=False)
    assert (protein.n_triangles, protein.n_unique_keys) == (70300, 45349)
    # Radix 30 since 1.0.0rc5 (ADR-107); 15c1cef9317d0605 at radix 20 before.
    assert protein.key_space.fingerprint.startswith("da39d73510c54a4c")
    assert protein.sampling.fingerprint.startswith("866885436e9a02bb")
    # Atoms labelled by element since 1.0.0rc7 (ADR-113); 1187735d9085815e
    # with atom-name labels before.  The keys themselves did not change here.
    assert atomic.key_space.fingerprint.startswith("1ecbef1b83344df6")
    assert atomic.sampling.fingerprint.startswith("9f3b52e7dca90beb")


def test_no_deferred_method_is_registered():
    deferred = re.compile(r"knn", re.I)
    assert len(list_methods()) == 20
    assert [m for m in list_methods() if deferred.search(m)] == []


# ---------------------------------------------------------------------------
# Metal ions: no method respells an atom name (ADR-082 is historical, ADR-108
# removes the narrow slice from every atomic path; test_input_rules.py checks
# the labels)
# ---------------------------------------------------------------------------


def test_no_admission_can_choose_a_narrow_atom_name_slice():
    """Atom names are read once, by the parser, from the full field; no
    admission carries a slice choice any more, so none can respell an ion."""
    import dataclasses

    from biotsr.selection.atomic import AtomAdmission

    names = {f.name for f in dataclasses.fields(AtomAdmission)}
    assert "atom_name_source" not in names
    assert "use_legacy_tokens" not in names


def test_entity_names_as_a_string_is_split_not_iterated():
    """``--entity-names 'NAG;FUC'`` used to become ('N','A','G',';',...)."""
    path = os.path.join(FIXTURE_PDB, "1KX5.pdb")
    _require(path)
    as_string = biotsr.analyze(path, method="carbohydrate_chain", chain="I",
                               entity_names="MN", keep_triplets=False)
    as_list = biotsr.analyze(path, method="carbohydrate_chain", chain="I",
                             entity_names=["MN"], keep_triplets=False)
    assert as_string.n_triangles == as_list.n_triangles == 35
    assert as_string.frequencies == as_list.frequencies
    assert biotsr.list_entity_targets(path, chain="I", method="carbohydrate_chain",
                                      entity_names="MN") == \
        biotsr.list_entity_targets(path, chain="I", method="carbohydrate_chain",
                                   entity_names=["MN"])


def _fixture(pdb_id):
    from ..conftest import structure

    return structure(pdb_id)


def test_entity_names_nag_fuc_string_reproduces_the_5wt9_chain_result():
    """``"NAG;FUC"`` is two names: NAG 14 + NAG 14 + FUC 10 atoms -> C(38,3) = 8,436."""
    path = _fixture("5WT9")
    as_string = biotsr.analyze(path, method="carbohydrate_chain", chain="A",
                               entity_names="NAG;FUC", keep_triplets=False)
    as_tuple = biotsr.analyze(path, method="carbohydrate_chain", chain="A",
                              entity_names=("NAG", "FUC"), keep_triplets=False)
    assert as_string.counts.n_points == as_tuple.counts.n_points == 38
    assert as_string.n_triangles == as_tuple.n_triangles == 8436
    assert as_string.frequencies == as_tuple.frequencies
    assert biotsr.list_entity_targets(path, chain="A", method="carbohydrate_chain",
                                      entity_names="NAG;FUC") == \
        (("1", "NAG"), ("2", "NAG"), ("3", "FUC"))


@pytest.mark.parametrize("pdb_id,chain,names,n_points,n_triangles", [
    ("1KX5", "I", "MN", 7, 35),          # one name, not ("M", "N")
    ("5WT9", "A", "NAG;FUC", 38, 8436),  # two names
])
def test_the_entity_names_flag_and_manifest_column_split_the_same_way(
        tmp_path, pdb_id, chain, names, n_points, n_triangles):
    from biotsr.cli.manifest import read_manifest, resolve_structure_path
    from biotsr.cli.parser import build_parser
    from biotsr.cli.spec import build_spec, execute

    path = _fixture(pdb_id)
    flag = vars(build_parser().parse_args(
        ["run", "-m", "carbohydrate_chain", "--pdb", path, "-c", chain,
         "--entity-names", names, "-o", str(tmp_path / "out")]))
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("pdb_path,method,chain,entity_names\n%s,carbohydrate_chain,%s,%s\n"
                        % (path, chain, names))
    row = read_manifest(str(manifest))[0]
    assert row["entity_names"] == names
    row["pdb"] = resolve_structure_path(row, "", row.pop("_where"))[0]
    row.pop("_row")
    for values in (flag, row):
        result = execute(build_spec(values)).results[""]
        assert (result.counts.n_points, result.n_triangles) == (n_points, n_triangles)


def test_cross_methods_split_an_entity_names_string_the_same_way():
    """7LF9 C/A: only NAG reaches the protein, so a split string finds 1,026
    triangles; the old character-by-character reading would have found none."""
    path = _fixture("7LF9")
    kw = dict(method="carbohydrate_protein_cross", chain="C", chain_2="A",
              keep_triplets=False)
    as_string = biotsr.analyze_cross(path, entity_names="NAG;MAN", **kw)
    as_tuple = biotsr.analyze_cross(path, entity_names=("NAG", "MAN"), **kw)
    assert as_string.n_triangles == as_tuple.n_triangles == 1026
    assert as_string.frequencies == as_tuple.frequencies
    assert as_string.sampling.fingerprint == as_tuple.sampling.fingerprint


# ---------------------------------------------------------------------------
# CI and deterministic evidence
# ---------------------------------------------------------------------------


def _workflow(name):
    return _read(".github", "workflows", name)


def test_ci_tests_every_python_the_classifiers_claim():
    claimed = set(re.findall(r'"Programming Language :: Python :: (3\.\d+)"', _read("pyproject.toml")))
    matrix = re.search(r'python-version: \[([^\]]+)\]', _workflow("ci.yml")).group(1)
    tested = set(re.findall(r'"(3\.\d+)"', matrix))
    assert claimed == tested, (sorted(claimed), sorted(tested))
    floor = re.search(r'requires-python = ">=(3\.\d+)"', _read("pyproject.toml")).group(1)
    assert min(tested, key=lambda v: tuple(map(int, v.split(".")))) == floor


def test_ci_runs_the_suite_the_gates_and_the_artifact_checks():
    ci = _workflow("ci.yml")
    for needle in ("python -m pytest", "release_artifacts.py",
                   "capability_matrix.py --check", "sanitize_tree.py --check",
                   "gates_api_hardening.py", 'test -z "$(git status --porcelain)"'):
        assert needle in ci, needle
    release = _workflow("release-validation.yml")
    for needle in ('-m ""', "release_artifacts.py", "release_manifest.py",
                   "gates_cli.py", "gates_api_hardening.py"):
        assert needle in release, needle


@pytest.mark.parametrize("name", ["ci.yml", "release-validation.yml", "publish.yml"])
def test_workflows_are_valid_yaml_without_local_paths(name):
    text = _workflow(name)
    assert "/Users/" not in text and "C:\\" not in text
    try:
        import yaml
    except ImportError:
        pytest.skip("%s: PyYAML (the structure is checked textually above)"
                    % OPTIONAL_DEPENDENCY)
    doc = yaml.safe_load(text)
    assert doc["jobs"] and all("runs-on" in job for job in doc["jobs"].values())


#: Triggers that would start a workflow without a person choosing to.
_AUTOMATIC_TRIGGER = re.compile(
    r"^  (push|pull_request|pull_request_target|release|schedule|workflow_run"
    r"|repository_dispatch|create|registry_package|deployment):", re.M)


@pytest.mark.parametrize("name", ["ci.yml", "release-validation.yml"])
def test_the_test_workflows_cannot_publish(name):
    """Pushes, pull requests, tags and test runs never publish anything: the
    test workflows hold no identity token, secret or upload step."""
    text = _workflow(name)
    for needle in ("id-token", "secrets.", "pypi-publish", "twine upload",
                   "gh release", "zenodo", "contents: write", "packages: write"):
        assert needle not in text.lower(), (name, needle)
    assert re.search(r"^permissions:\n  contents: read$", text, re.M)


def test_publishing_is_by_hand_from_a_reviewed_tag_after_approval():
    """publish.yml: started by a person only, from a release tag, for files
    whose hashes were reviewed, after a reviewer approves the environment; it
    publishes through Trusted Publishing and holds no stored credential."""
    text = _workflow("publish.yml")
    on = text.split("\non:\n", 1)[1].split("\npermissions:", 1)[0]
    assert re.search(r"^  workflow_dispatch:$", on, re.M)
    assert not _AUTOMATIC_TRIGGER.search(on)
    assert re.search(r"^permissions: \{\}$", text, re.M)
    verify, publish = text.split("\n  verify:\n", 1)[1].split("\n  publish:\n", 1)
    # Only the publish job, after the verification and the approval, can ask
    # for an identity token.
    assert "id-token" not in verify and text.count("id-token: write") == 1
    assert "needs: verify" in publish and "id-token: write" in publish
    assert re.search(r"^    environment:\n      name: \$\{\{ inputs\.target \}\}$",
                     publish, re.M)
    assert "options: [testpypi, pypi]" in on and "default: testpypi" in on
    for needle in ('"$GITHUB_REF_TYPE" != tag', '"v$version"', "origin/main",
                   "required_reviewers", ".github/workflows/release-validation.yml",
                   "conclusion", "RELEASE_CANDIDATE_MANIFEST.json",
                   "$WHEEL_SHA256", "$SDIST_SHA256", "sha256sum --check --strict",
                   "twine check"):
        assert needle in verify, needle
    assert "sha256sum --check --strict" in publish
    # No build: the files published are those release validation tested.
    assert "python -m build" not in text and "release_artifacts.py" not in text
    # Trusted Publishing, by a pinned action; nothing stored.
    uses = re.findall(r"pypa/gh-action-pypi-publish@(\S+)", text)
    assert uses and all(re.fullmatch(r"[0-9a-f]{40}", ref) for ref in uses)
    code = "\n".join(line for line in text.splitlines()
                     if not line.lstrip().startswith("#"))
    for needle in ("secrets.", "password", "api-token", "user:"):
        assert needle not in code, needle


@pytest.mark.parametrize("rel", [
    "validation/evidence/cli/GATE_RESULTS.json",
    "validation/evidence/api_hardening/GATE_RESULTS.json",
])
def test_committed_gate_evidence_carries_no_timing(rel):
    assert "elapsed_seconds" not in _read(*rel.split("/"))


def test_timing_is_kept_out_of_version_control():
    assert "validation/**/*_RUNTIME.json" in _read(".gitignore").splitlines()


# ---------------------------------------------------------------------------
# Public tree
# ---------------------------------------------------------------------------


def test_head_carries_no_workstation_path_or_account_identifier():
    module = _script("sanitize_tree")
    try:
        tracked = module.tracked_files()
    except Exception:  # noqa: BLE001 -- not a git checkout (e.g. an unpacked sdist)
        pytest.skip("%s: not a git checkout" % REPOSITORY_ONLY)
    pending = []
    for rel in tracked:
        path = os.path.join(WORKSPACE, rel)
        if rel.endswith(".xlsx"):
            if module.sanitize_xlsx(path)[1]:
                pending.append(rel)
        elif rel.endswith(module.TEXT_SUFFIXES):
            with open(path, encoding="utf-8", newline="") as fh:
                if module.sanitize_text(fh.read())[1]:
                    pending.append(rel)
    assert pending == []


def test_no_stray_historical_output_at_the_root():
    name = "proteinNumKeysDist.csv"
    assert not os.path.exists(os.path.join(WORKSPACE, name)) or not _is_tracked(name), name


#: The research group's private testing inputs -- a dataset sheet and its column
#: dictionary -- by file name and by SHA-256 (as last committed), so neither a
#: renamed copy nor a re-added original can enter the tracked tree.  The
#: originals live outside the repository; the tests use a synthetic sheet.
#: Spreadsheets are never tracked; the private inputs are also recognised by the
#: SHA-256 of their contents, so a renamed copy is caught without naming them.
PRIVATE_INPUT_NAMES = re.compile(r"\.(xls[xmb]?|ods|numbers)$", re.I)
PRIVATE_INPUT_SHA256 = {
    "4f83327283d29f312da3b55a712671d8653723c572a219f97677088cfcc2e339",
    "77a649468e2735f5104e132ebc0f10c4ff9cf3e6f0a23edd8d316df4a8ca9352",
}


def _tracked_files():
    proc = subprocess.run(["git", "-C", WORKSPACE, "ls-files", "-z"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if proc.returncode != 0:          # not a git checkout (an unpacked sdist)
        return [os.path.relpath(os.path.join(d, f), WORKSPACE)
                for d, _dirs, files in os.walk(WORKSPACE) for f in files]
    return [p for p in proc.stdout.decode("utf-8").split("\0") if p]


def test_private_testing_inputs_are_not_in_the_tree():
    """Neither by name nor, renamed, by content."""
    named, copied = [], []
    for rel in _tracked_files():
        if PRIVATE_INPUT_NAMES.search(os.path.basename(rel)):
            named.append(rel)
        path = os.path.join(WORKSPACE, rel)
        if not os.path.isfile(path) or os.path.getsize(path) > 50_000:
            continue
        with open(path, "rb") as fh:
            if hashlib.sha256(fh.read()).hexdigest() in PRIVATE_INPUT_SHA256:
                copied.append(rel)
    assert named == [] and copied == [], (named, copied)


#: The repository's top level, entry by entry.  Anything else -- a results or
#: download directory, a backup, a local archive -- was added by accident.
TOP_LEVEL = {".gitattributes", ".github", ".gitignore", "CHANGELOG.md", "CITATION.cff",
             "LICENSE", "MANIFEST.in", "NOTICE", "README.md", "RELEASE_CHECKLIST.md",
             "docs", "pyproject.toml", "src", "tests", "validation"}

#: What running BioTSR writes.  Only the evidence under validation/ and the
#: fixtures under tests/data/ are tracked outputs.
RUN_OUTPUT = re.compile(
    r"(^|/)(run\.json|batch_summary\.csv|batch_failures\.csv|fetch_report\.csv"
    r"|\.biotsr-batch\.lock|\.biotsr-fetch/.*|\.[^/]+\.partial-[^/]+/.*)$"
    r"|\.(pdb|ent)\.gz$|\.triplets_|\.keys_")
ARCHIVE = re.compile(r"\.(bundle|tar|tgz|zip|7z|gz|bz2|xz)$")
PRIVATE_GUIDE = re.compile(r"(^|/)private/|\.private\.", re.I)
#: Well above the largest tracked reference (a 16.5 MiB key file), far below
#: GitHub's 100 MB limit.
MAX_TRACKED_BYTES = 25 * 1024 * 1024


def test_the_tracked_tree_holds_only_release_content():
    """A shared check, independent of any one clone's local excludes or hooks."""
    proc = subprocess.run(["git", "-C", WORKSPACE, "ls-files", "-z"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if proc.returncode != 0:
        pytest.skip("%s: not a git checkout" % REPOSITORY_ONLY)
    tracked = [p for p in proc.stdout.decode("utf-8").split("\0") if p]
    problems = sorted(
        {"unexpected top-level entry: " + p.split("/", 1)[0]
         for p in tracked if p.split("/", 1)[0] not in TOP_LEVEL}
        | {"run output outside the evidence: " + p for p in tracked
           if RUN_OUTPUT.search(p) and not p.startswith(("validation/", "tests/data/"))}
        | {"archive: " + p for p in tracked
           if ARCHIVE.search(p) and not p.startswith("tests/data/")}
        | {"private material: " + p for p in tracked if PRIVATE_GUIDE.search(p)}
        | {"larger than %d MiB: %s" % (MAX_TRACKED_BYTES >> 20, p) for p in tracked
           if os.path.getsize(os.path.join(WORKSPACE, p)) > MAX_TRACKED_BYTES}
    )
    assert problems == []


#: Files that define the patterns below rather than contain a path.
_PATH_PATTERN_DEFINITIONS = {"validation/scripts/sanitize_tree.py",
                             "validation/scripts/release_artifacts.py",
                             "tests/unit/test_release_readiness.py"}
_HOME_PATH = re.compile(r"/Users/(?!<)[^/\s\"'<>]+/|/home/(?!<)[A-Za-z][\w.-]*/"
                        r"|[A-Z]:\\\\?Users\\\\?")


def test_no_tracked_text_names_a_home_directory():
    """Anywhere in the tree, not only where sanitize_tree.py rewrites."""
    proc = subprocess.run(["git", "-C", WORKSPACE, "grep", "-I", "-l", "-E",
                           r"/Users/|/home/|Users\\"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if proc.returncode not in (0, 1):
        pytest.skip("%s: not a git checkout" % REPOSITORY_ONLY)
    offenders = []
    for rel in proc.stdout.decode("utf-8").split():
        if rel in _PATH_PATTERN_DEFINITIONS:
            continue
        with open(os.path.join(WORKSPACE, rel), encoding="utf-8", errors="replace") as fh:
            if _HOME_PATH.search(fh.read()):
                offenders.append(rel)
    assert offenders == []


def _is_tracked(rel):
    proc = subprocess.run(["git", "-C", WORKSPACE, "ls-files", "--error-unmatch", rel],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return proc.returncode == 0


def test_no_runtime_message_names_an_internal_development_stage():
    """Users have never heard of development stages; messages must not name one."""
    offenders = []
    for root, _dirs, files in os.walk(os.path.join(SRC, "biotsr")):
        for f in files:
            if not f.endswith(".py"):
                continue
            path = os.path.join(root, f)
            tree = ast.parse(open(path, encoding="utf-8").read())
            docstrings = set()
            for node in ast.walk(tree):
                if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                     ast.AsyncFunctionDef)) and node.body:
                    first = node.body[0]
                    if isinstance(first, ast.Expr) and isinstance(
                            getattr(first, "value", None), ast.Constant):
                        docstrings.add(id(first.value))
            for node in ast.walk(tree):
                if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                        and id(node) not in docstrings
                        and re.search(r"\bStage 2[A-K]\b", node.value)):
                    offenders.append("%s:%d" % (os.path.relpath(path, SRC), node.lineno))
    assert offenders == []


# ---------------------------------------------------------------------------
# Internal development material stays out of what is released
# ---------------------------------------------------------------------------

#: Material that belongs to the private development record, not to the release:
#: named or implied individuals, development-stage labels, internal risk ids,
#: internal document names, and candidate-strategy work that is not published.
#: Stage-numbered script, directory and key names (``stage2d_...``) count too:
#: validation material is named after what it validates.  Internal document
#: names are matched by the SHA-256 of their lower-cased spelling, so this test
#: does not itself publish them.
INTERNAL_MATERIAL = re.compile(
    r"supervisor|(?-i:\bStage ?[12][A-K]?\b)|(?<![a-z])stage_?[12][a-k]"
    r"|\brisk A\d|risk register|Future.standard"
    r"|docs/history|docs/research|validation/reports"
    r"|unpublished"
    r"|release-1\b|closeout|release review|\brc\d_[0-9a-f]{7}",
    re.I,
)

INTERNAL_DOCUMENT_SHA256 = {
    "01c6db03b3449cb46e638f3d9761404d0409b4584fe410f2e12b03e85c072d75",
    "0b978514d6028599be4e93869cd245d0e8f15beaf9e4a0c9834d025b28606563",
    "0f564c3901417f683befdec32dfe1cf24953f26155f9b6321fa200a68e5a4a5f",
    "33d23fadefd8d47558ee521575cde0719ccec9c942a4c60ea4f8b74e3983e546",
    "66460773a992db2423fdf243da7aa78740023a0ca55574407533aeed995804a9",
    "acc7e68bdce1d17ccde8e7ad4093ca19b77d04ccd390ade3d93b15b77971434a",
    "bfafa208b5c5437bc6f3ae58364e087458639c43c3e26dfecd3625cdd9b117ae",
}


def _names_an_internal_document(line):
    return any(hashlib.sha256(word.lower().encode()).hexdigest() in INTERNAL_DOCUMENT_SHA256
               for word in re.findall(r"[A-Za-z][A-Za-z_]+[A-Za-z]", line))

#: Candidate-strategy work that is not published, by the SHA-256 of its
#: lower-cased name, so this test does not itself publish the names.
UNPUBLISHED_NAME_SHA256 = {
    "6da08dcd036eb07a156659edf09657883b729602751654ed852f1b71f686bddb",
    "a695f06537aa3989261fcc449621cc7ad54b43f1d57fcd1004c71df2c8e2975e",
    "d31ca2c924f2f96c19e819059d520cfb46033ce2eb28a774e14aeb6315c2f5cc",
}


def _names_unpublished_work(line):
    return any(hashlib.sha256(word.lower().encode()).hexdigest() in UNPUBLISHED_NAME_SHA256
               for word in re.findall(r"[A-Za-z][A-Za-z0-9]+", line))

#: (file, matched text) pairs that are the guard rather than the material: the
#: release-artifact check refuses a re-added ``docs/history/`` by this pattern.
INTERNAL_MATERIAL_GUARDS = {
    ("validation/scripts/release_artifacts.py", "docs/history"),
}


def _released_text_files():
    """Every text file a user of the package, the sdist, the documentation or
    the validation record reads: code, comments, docstrings, prose, and the
    text fields of the committed evidence (not the key and triplet data)."""
    out = ["README.md", "CHANGELOG.md", "RELEASE_CHECKLIST.md", "NOTICE",
           "CITATION.cff", "MANIFEST.in", "pyproject.toml", ".gitignore",
           ".gitattributes", ".github/workflows/ci.yml",
           ".github/workflows/release-validation.yml", ".github/workflows/publish.yml"]
    suffixes = {"docs": (".md", ".py"), os.path.join("src", "biotsr"): (".md", ".py"),
                "tests": (".md", ".py", ".json", ".csv"),
                "validation": (".md", ".py", ".json", ".csv", ".txt")}
    for root, wanted in suffixes.items():
        for dirpath, _dirs, files in os.walk(os.path.join(WORKSPACE, root)):
            out += [os.path.relpath(os.path.join(dirpath, f), WORKSPACE)
                    for f in files if f.endswith(wanted)]
    this = os.path.relpath(os.path.abspath(__file__.replace(".pyc", ".py")), WORKSPACE)
    return sorted(p for p in out if p != this)


def test_released_files_carry_no_internal_development_material():
    """Documentation, package source (code, comments and docstrings), metadata
    and CI configuration name no individual, no development stage, no internal
    record and no unpublished method.  Scientific rationale stays; its
    provenance is cited as an ADR or a public page."""
    offenders = []
    for rel in _released_text_files():
        path = os.path.join(WORKSPACE, rel)
        if not os.path.exists(path):
            continue
        # Untracked files are not released: private notes,
        # caches, and locally regenerated evidence.
        if not _tracked(rel):
            continue
        for number, line in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
            if _names_unpublished_work(line) or _names_an_internal_document(line):
                offenders.append("%s:%d: %s" % (rel, number, line.strip()[:100]))
                continue
            for match in INTERNAL_MATERIAL.finditer(line):
                if (rel, match.group(0)) not in INTERNAL_MATERIAL_GUARDS:
                    offenders.append("%s:%d: %s" % (rel, number, line.strip()[:100]))
                    break
    assert offenders == [], "\n".join(offenders)


def _tracked(rel):
    import subprocess

    try:
        out = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=WORKSPACE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        return True
    return out.returncode == 0 or out.returncode == 128  # 128: not a git checkout


# ---------------------------------------------------------------------------
# The release bump reaches every file that names the version
# ---------------------------------------------------------------------------

#: Every file that names the package version or its pre-release status, and so
#: changes with the release bump (docs/RELEASING.md, section 7).  An unpacked
#: sdist has no docs/, so only the files present are checked.
VERSION_BEARING = (
    "src/biotsr/__init__.py", "CITATION.cff", "CHANGELOG.md", "README.md",
    "docs/USER_GUIDE.md", "docs/CLI_GUIDE.md", "docs/API_GUIDE.md",
    "docs/PACKAGE_OVERVIEW.md", "docs/REPRODUCIBILITY_NOTE.md",
)

_PRE_RELEASE_WORDING = re.compile(
    r"\b\d+\.\d+\.\d+(?:a|b|rc)\d+\b|once v?1\.0 is (?:released|published)"
    r"|not on PyPI yet|not yet published|not yet released|release candidate",
    re.I,
)


def test_every_version_bearing_file_follows_the_version():
    final = re.search(r"(?:a|b|rc)\d+$", biotsr.__version__) is None
    present = [rel for rel in VERSION_BEARING
               if os.path.exists(os.path.join(WORKSPACE, rel))]
    assert "README.md" in present and "CITATION.cff" in present
    offenders = []
    for rel in present:
        with open(os.path.join(WORKSPACE, rel), encoding="utf-8") as fh:
            text = fh.read()
        wording = _PRE_RELEASE_WORDING.search(text)
        if final and wording:
            offenders.append("%s: %r" % (rel, wording.group(0)))
        elif not final and not (wording or biotsr.__version__ in text):
            # The list is only useful while it is complete in both directions.
            offenders.append("%s: names neither the version nor its status" % rel)
    assert offenders == []
