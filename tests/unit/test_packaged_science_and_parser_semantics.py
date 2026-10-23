"""Release consistency checks: branding, packaged science, parser semantics.

These are the assertions the release's consistency is judged by.  They are deliberately *specific*: each one names the fact it
protects, so that a future change that breaks it has to argue with the fact
rather than with a vague expectation.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
import sys

import pytest

import biotsr
from biotsr.labels.elements import (
    ELEMENT_LABELS_112,
    ELEMENT_LABELS_RESOURCE,
)
from biotsr.domain.interpretation import (
    LEGACY_COLUMNS_V1,
    SPEC_COLUMNS_V1,
    IncompatibleStructureInterpretation,
    interpretation_of_backend,
    require_common_structure_interpretation,
    structure_interpretation_of,
)

#: The project directory (the repository root): this file is ``tests/unit/<name>.py``.
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(HERE, "src")
PKG = os.path.join(SRC, "biotsr")


# ---------------------------------------------------------------------------
# Branding
# ---------------------------------------------------------------------------

def test_distribution_name_is_biotsr():
    text = open(os.path.join(HERE, "pyproject.toml"), encoding="utf-8").read()
    assert 'name = "biotsr"' in text


def test_import_namespace_is_biotsr():
    assert biotsr.__name__ == "biotsr"
    assert os.path.basename(os.path.dirname(biotsr.__file__)) == "biotsr"


def test_console_script_and_module_entry_are_biotsr():
    text = open(os.path.join(HERE, "pyproject.toml"), encoding="utf-8").read()
    assert 'biotsr = "biotsr.cli:main"' in text
    assert os.path.isfile(os.path.join(PKG, "__main__.py"))


def test_cli_identifies_itself_as_biotsr():
    from biotsr.cli.parser import build_parser

    parser = build_parser()
    assert parser.prog == "biotsr"
    assert "BioTSR" in (parser.description or "")


def test_version_flag_prints_biotsr():
    """`biotsr --version` must not announce itself as `biotsr`."""
    from biotsr import __version__
    from biotsr.cli.parser import build_parser

    parser = build_parser()
    action = next(
        a for a in parser._actions if "--version" in getattr(a, "option_strings", [])
    )
    assert action.version == "biotsr %s" % __version__


def test_describe_renders_the_command_as_biotsr():
    """`biotsr describe` prints the command that produced an output."""
    import io

    from biotsr.cli.commands import _describe_run
    from biotsr.cli.console import Console

    buf = io.StringIO()
    console = Console(stdout=buf, stderr=io.StringIO())
    _describe_run(
        {"command": "run", "biotsr_version": "0.0.0", "method": "protein",
         "profile": {}, "parameters": {}, "results": []},
        console,
    )
    text = buf.getvalue()
    assert "command          biotsr run" in text, text
    assert "\ntsr " not in text, text


def test_no_stale_tsr_unified_branding_in_the_package():
    """The retired product name must not survive anywhere a user can see it."""
    stale = ("TSR-Unified", "TSR Unified", "tsr-unified")
    offenders = []
    for dirpath, _dirs, files in os.walk(PKG):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(dirpath, name)
            text = open(path, encoding="utf-8").read()
            for needle in stale:
                if needle in text:
                    offenders.append("%s: %s" % (path, needle))
    assert not offenders, offenders


#: A development-stage label ("Stage" followed by a stage number and letter).
STAGE_LABEL = re.compile(r"\bStage ?[12][A-K]?\b")


def test_no_user_facing_development_stage_wording_survives():
    """A released package must not tell a user which development stage it was at."""
    offenders = []
    for dirpath, _dirs, files in os.walk(PKG):
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.join(dirpath, name)
            for lineno, line in enumerate(open(path, encoding="utf-8"), 1):
                stripped = line.strip()
                if not STAGE_LABEL.search(line):
                    continue
                # Docstrings and comments are development history and may say so;
                # a runtime *message* may not.
                if stripped.startswith("#") or stripped.startswith("*"):
                    continue
                if '"' in line or "'" in line:
                    offenders.append("%s:%d: %s" % (path, lineno, stripped))
    assert not offenders, offenders


# ---------------------------------------------------------------------------
# Scientific resources: the element label table ships, and is the approved one
# ---------------------------------------------------------------------------

def test_canonical_label_table_is_package_data_in_pyproject():
    text = open(os.path.join(HERE, "pyproject.toml"), encoding="utf-8").read()
    assert 'biotsr = ["labels/resources/element_labels_v1.csv"]' in text


def test_every_shipped_label_resource_is_present_on_disk():
    resources = os.path.join(PKG, "labels", "resources")
    present = sorted(f for f in os.listdir(resources) if f.endswith(".csv"))
    assert present == [
        "element_labels_v1.csv",
    ], present


def test_element_labels_load_through_importlib_resources_not_cwd(tmp_path, monkeypatch):
    """The canonical resource must load with the CWD somewhere else entirely."""
    monkeypatch.chdir(tmp_path)
    from biotsr.labels.elements import element_label_scheme

    scheme = element_label_scheme()
    assert scheme.resolve_atom("SG", "S") == 111


def test_element_label_content_is_stable():
    """The approved resource is pinned by hash.

    If this fails, the atomic labels changed.  That is allowed only
    deliberately, as a new scheme with a recorded decision (ADR-113), because
    every atomic key space changes with it.
    """
    path = os.path.join(PKG, "labels", "resources", ELEMENT_LABELS_RESOURCE)
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    assert digest == "d7d8e489ffbb5dab7a433a61d55dd1a1b1680ef56cbfd67764a7b8971c335d9f", digest


# --- names that once decided a label --------------------------------------

@pytest.mark.parametrize(
    "atom_name,element,expected,why",
    [
        ("CA", "C", 106, "protein alpha carbon"),
        ("CA", "CA", 99, "calcium, the element"),
        ("CL1", "CL", 93, "a ligand chlorine"),
        ("CL1", "C", 106, "a ligand carbon named CL1"),
        ("NE", "N", 108, "protein nitrogen NE"),
        ("NE", "NE", 51, "neon, the element"),
        ("SG", "S", 111, "protein sulfur gamma"),
        ("SG", "SG", 72, "seaborgium, the element"),
        ("ND", "N", 108, "protein ND"),
        ("ND", "ND", 50, "neodymium, the element"),
        ("OS", "O", 109, "protein OS"),
        ("OS", "OS", 54, "osmium, the element"),
        ("HO", "HO", 37, "holmium, the element"),
        ("PA", "P", 110, "ATP alpha phosphorus"),
        ("PB", "P", 110, "ATP beta phosphorus"),
    ],
)
def test_atom_name_element_ambiguities_resolve_by_the_element(atom_name, element, expected, why):
    assert ELEMENT_LABELS_112.resolve_atom(atom_name, element) == expected, why


# NOTE: "no CYS/SG residue-specific special case exists in the package"
# is asserted authoritatively by tests/unit/test_element_labels.py.
# It is deliberately NOT restated here: a second, cruder copy of that grep would
# either duplicate the guarantee or -- worse -- weaken it, because the residue
# name CYS legitimately appears in the amino-acid label table, in the per-residue
# expected-atom-count table and in the cross-entity side rules.


# ---------------------------------------------------------------------------
# Parser semantics
# ---------------------------------------------------------------------------

def test_only_one_parser_backend_ships():
    backends = os.path.join(PKG, "structure", "backends")
    modules = sorted(
        f for f in os.listdir(backends)
        if f.endswith(".py") and f != "__init__.py"
    )
    assert modules == ["pdb_legacy.py"], modules


def test_pdb_strict_is_not_implemented_and_is_refused():
    assert not os.path.exists(
        os.path.join(PKG, "structure", "backends", "pdb_strict.py")
    )
    with pytest.raises(biotsr.InvalidRequest):
        biotsr.load_structure(__file__, backend="pdb_strict")


# Every method runs under its own standard profile only. A profile *name* (or
# any non-profile value) used to fail with an AttributeError deep inside; it is
# refused up front, saying to omit ``profile=``. Cross methods have no ``mode=``.
_UBQ = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "structures", "1UBQ.pdb")


@pytest.mark.parametrize("call, kwargs, hint", [
    (biotsr.analyze, dict(method="protein", profile="legacy_tsr_package_v1"),
     "no profile to choose by name"),
    (biotsr.analyze, dict(method="protein", profile="standard"), "Omit profile="),
    (biotsr.analyze, dict(method="protein", profile="no_such_profile"), "Omit profile="),
    (biotsr.analyze, dict(method="protein", profile=3), "not int"),
    (biotsr.analyze_cross, dict(method="protein_protein_calpha_cross", chain_2="A",
                                profile="legacy"), "Omit profile="),
    (biotsr.analyze_cross, dict(method="protein_protein_calpha_cross", chain_2="A",
                                mode="legacy"), "Omit mode="),
])
def test_a_profile_name_is_refused_saying_what_to_pass(call, kwargs, hint):
    with pytest.raises(biotsr.InvalidRequest) as info:
        call(_UBQ, chain="A", **kwargs)
    assert hint in str(info.value)


def test_the_methods_own_default_profile_is_accepted():
    profile = biotsr.methods.get_method("protein").default_profile
    result = biotsr.analyze(_UBQ, method="protein", chain="A", profile=profile,
                            keep_triplets=False)
    assert result.n_triangles == 70300


def test_interpretation_is_coarser_than_the_backend_name():
    """Coordinate precision is not an atom-naming semantic."""
    assert interpretation_of_backend("pdb_legacy") == LEGACY_COLUMNS_V1
    assert interpretation_of_backend("pdb_legacy+float32_coords") == LEGACY_COLUMNS_V1


def test_two_backends_with_the_same_semantics_stay_comparable():
    """The property that makes this durable rather than a backend-name check."""
    assert interpretation_of_backend("pdb_legacy") == interpretation_of_backend(
        "pdb_legacy+float32_coords"
    )


def test_a_different_interpretation_is_a_different_semantic():
    assert interpretation_of_backend("pdb_strict") == SPEC_COLUMNS_V1
    assert SPEC_COLUMNS_V1 != LEGACY_COLUMNS_V1


def test_an_undeclared_backend_is_refused_rather_than_assumed():
    with pytest.raises(IncompatibleStructureInterpretation):
        interpretation_of_backend("something_new")


def test_structure_interpretation_is_not_in_either_identity():
    """Adding it to KeySpace would break every saved result; it is not there."""
    from biotsr.domain.keyspace import KeySpace, SamplingProfile
    from dataclasses import fields

    for cls in (KeySpace, SamplingProfile):
        names = {f.name for f in fields(cls)}
        assert "structure_interpretation" not in names
        assert "parser_backend" not in names
