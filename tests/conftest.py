"""Shared test fixtures and paths."""

from __future__ import annotations

import os
import sys

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(TESTS_DIR)
SRC = os.path.join(PACKAGE_ROOT, "src")
# The package project is the repository root, so the workspace -- where
# validation/ and docs/ live -- is the package root itself.
WORKSPACE = PACKAGE_ROOT

if SRC not in sys.path:
    sys.path.insert(0, SRC)

#: The small, unmodified RCSB structures the tests read. They ship in the source
#: distribution, so the computational tests run from an unpacked sdist too.
STRUCTURES = os.path.join(TESTS_DIR, "data", "structures")
#: Synthetic structures written by a script in that folder (its README).
SYNTHETIC = os.path.join(TESTS_DIR, "data", "synthetic")
FIXTURE_PDB = STRUCTURES
#: Small fixed reference data the tests compare against (standard results frozen
#: from a tested build, and compact historical facts); see its README.
REFERENCE_DATA = os.path.join(TESTS_DIR, "data", "reference")


#: Every skip for missing data starts with one of these, so a run can tell an
#: *expected* skip from an accidental one. The source distribution carries the
#: package, its tests and the data in ``tests/data``; the documentation and the
#: release records stay in the repository, so there those tests skip under
#: :data:`REPOSITORY_ONLY` and nothing else.
REPOSITORY_ONLY = "repository-only evidence not present"
#: A third-party package some tests use when it is installed.
OPTIONAL_DEPENDENCY = "optional dependency not installed"

#: ``BIOTSR_REQUIRE_EVIDENCE=1`` (set in repository CI) turns a missing piece of
#: tracked evidence into a failure. Without it, a move that updated the data but
#: not a path would turn the suite green with hundreds of silent skips.
REQUIRE_EVIDENCE = os.environ.get("BIOTSR_REQUIRE_EVIDENCE") == "1"


def _missing(path) -> str:
    return "%s: %s" % (REPOSITORY_ONLY, os.path.relpath(path, WORKSPACE))


def _require(path):
    """``path``, or skip -- or fail, under ``BIOTSR_REQUIRE_EVIDENCE=1``."""
    if not os.path.exists(path):
        if REQUIRE_EVIDENCE:
            pytest.fail(_missing(path))
        pytest.skip(_missing(path))
    return path


def evidence_available(path) -> bool:
    """For module-level parametrisation, where a test cannot yet skip.

    Missing evidence under ``BIOTSR_REQUIRE_EVIDENCE=1`` is a collection error,
    which fails the run as loudly as a failed test.
    """
    if os.path.exists(path):
        return True
    if REQUIRE_EVIDENCE:
        raise RuntimeError(_missing(path))
    return False


def evidence_cases(cases, path):
    """``cases``, or one explicitly skipped placeholder when ``path`` is absent.

    Keeps the skip reason under :data:`REPOSITORY_ONLY` instead of pytest's
    generic "got empty parameter set".
    """
    if evidence_available(path) and cases:
        return cases
    return [pytest.param(None, marks=pytest.mark.skip(reason=_missing(path)))]


def structure(pdb_id):
    """The path of one tracked test structure, e.g. ``structure("1UBQ")``."""
    return _require(os.path.join(STRUCTURES, "%s.pdb" % pdb_id))


@pytest.fixture(scope="session")
def ubq_pdb():
    return _require(os.path.join(FIXTURE_PDB, "1UBQ.pdb"))


# ---------------------------------------------------------------------------
# Atomic-family fixtures: one per family, all pinned with a recorded provenance in
# tests/data/structures/manifest.json.
# ---------------------------------------------------------------------------

ACCEPTANCE_PDB = STRUCTURES


@pytest.fixture(scope="session")
def mapk_path():
    """6G9M -- protein with one ESW ligand (404); six CYS."""
    return _require(os.path.join(ACCEPTANCE_PDB, "6G9M.pdb"))


@pytest.fixture(scope="session")
def trypsin_path():
    """3PTB -- the classic small-ligand fixture (BEN/1), historical drug baseline."""
    return _require(os.path.join(FIXTURE_PDB, "3PTB.pdb"))


@pytest.fixture(scope="session")
def glyco_path():
    """5WT9 -- historical carbohydrate baseline (NAG/1 monomer, NAG+NAG+FUC chain)."""
    return _require(os.path.join(FIXTURE_PDB, "5WT9.pdb"))


@pytest.fixture(scope="session")
def lf9_path():
    """7LF9 -- carbohydrate NAG 1/2/5 in chain D."""
    return _require(os.path.join(ACCEPTANCE_PDB, "7LF9.pdb"))


@pytest.fixture(scope="session")
def bdna_path():
    """1BNA -- the historical nucleotide baseline, B-DNA dodecamer."""
    return _require(os.path.join(FIXTURE_PDB, "1BNA.pdb"))


@pytest.fixture(scope="session")
def ac0_path():
    """2AC0 -- DNA with DT/DC in chain E and DA in chain F, p53 protein chains."""
    return _require(os.path.join(ACCEPTANCE_PDB, "2AC0.pdb"))


# ---------------------------------------------------------------------------
# SSE fixtures: structures chosen for what secondary structure they contain.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def barnase_path():
    """1BRS -- 21 HELIX and 23 SHEET records, the richest SSE fixture pinned."""
    return _require(os.path.join(FIXTURE_PDB, "1BRS.pdb"))


@pytest.fixture(scope="session")
def p53_path():
    """1TUP -- 10 HELIX, 33 SHEET; a beta-sandwich, so sheet-heavy categories."""
    return _require(os.path.join(FIXTURE_PDB, "1TUP.pdb"))


@pytest.fixture(scope="session")
def nucleosome_path():
    """1KX5 -- 38 HELIX records; the helix-heavy counterpart to 1TUP."""
    return _require(os.path.join(FIXTURE_PDB, "1KX5.pdb"))
