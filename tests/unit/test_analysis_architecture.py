"""Structural guards on the analysis layer.

The layer's value is that it is *generic*: one implementation over
``key -> count`` plus two identities.  That property is easy to state and easy to
lose -- one ``from ..geometry`` import, one ``if result.method ==`` branch, and
the next family needs its own copy.  These tests read the source and refuse
those, so the guarantee survives people who have not read the architecture
document.
"""

from __future__ import annotations

import ast
import os

import pytest

SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "src",
    "biotsr",
)
ANALYSIS = os.path.join(SRC, "analysis")

#: The layer may reach into ``domain`` for the two identities and the result
#: type, and nowhere else in the package.
ALLOWED_PACKAGE_IMPORTS = {"domain"}

#: Engine-side concerns the analysis layer must never acquire.
FORBIDDEN_PACKAGE_IMPORTS = (
    "geometry",
    "encoding",
    "binning",
    "structure",
    "selection",
    "candidates",
    "methods",
    "modifiers",
    "profiles",
    "labels",
    "cross",
    "io",
    "api",
)


def _modules():
    for name in sorted(os.listdir(ANALYSIS)):
        if name.endswith(".py"):
            yield name, os.path.join(ANALYSIS, name)


def _tree(path):
    with open(path) as fh:
        return ast.parse(fh.read())


def test_the_analysis_package_has_modules():
    assert [n for n, _ in _modules()], "no analysis modules found"


@pytest.mark.parametrize("name", [n for n, _ in _modules()])
def test_no_analysis_module_reaches_into_the_engine(name):
    """No geometry, no encoding, no binning, no parser, no method registry."""
    path = os.path.join(ANALYSIS, name)
    with open(path) as fh:
        text = fh.read()
    offenders = [
        module for module in FORBIDDEN_PACKAGE_IMPORTS if ("from ..%s" % module) in text
    ]
    assert offenders == [], "%s imports %r from the engine" % (name, offenders)


@pytest.mark.parametrize("name", [n for n, _ in _modules()])
def test_analysis_imports_from_tsr_go_only_to_domain_or_siblings(name):
    """The only package-level dependency is ``domain``; everything else is a sibling."""
    tree = _tree(os.path.join(ANALYSIS, name))
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 2:
            top = (node.module or "").split(".")[0]
            if top not in ALLOWED_PACKAGE_IMPORTS:
                offenders.append(node.module)
        elif isinstance(node, ast.ImportFrom) and node.level and node.level > 2:
            offenders.append("level %d: %s" % (node.level, node.module))
    assert offenders == [], "%s reaches outside domain: %r" % (name, offenders)


@pytest.mark.parametrize("name", [n for n, _ in _modules()])
def test_no_analysis_module_branches_on_the_method_name(name):
    """The whole point: the layer cannot tell which engine made a result.

    Reading ``annotations["method"]`` as provenance is fine; comparing it is not.
    """
    tree = _tree(os.path.join(ANALYSIS, name))
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Compare):
            for side in [node.left] + list(node.comparators):
                if isinstance(side, ast.Attribute) and side.attr == "method":
                    offenders.append(ast.dump(node)[:120])
    assert offenders == [], "%s branches on result.method: %r" % (name, offenders)


@pytest.mark.parametrize("name", [n for n, _ in _modules()])
def test_no_analysis_module_parses_a_pdb_file(name):
    """PDB parsing lives in ``structure``; the analysis layer never sees a file."""
    with open(os.path.join(ANALYSIS, name)) as fh:
        text = fh.read().lower()
    for marker in ("atom  ", "hetatm", ".pdb", "line[30:38]", "mmcif"):
        assert marker not in text, "%s looks like it parses structures (%r)" % (name, marker)


@pytest.mark.parametrize("name", [n for n, _ in _modules()])
def test_no_analysis_module_recomputes_a_key_or_a_triangle(name):
    """Keys arrive already encoded; geometry is not re-derived here."""
    tree = _tree(os.path.join(ANALYSIS, name))
    forbidden = {
        "acos",
        "atan2",
        "degrees",
        "radians",
        "cos",
        "sin",
        "encode_key",
        "theta_class",
        "distance_class",
    }
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and node.attr in forbidden:
            offenders.append(node.attr)
        elif isinstance(node, ast.Name) and node.id in forbidden:
            offenders.append(node.id)
    assert offenders == [], "%s does geometry or encoding work: %r" % (name, offenders)


def test_the_comparison_rule_is_enforced_in_one_place():
    """Every entry point routes through ``check_comparable``.

    A second, independent copy of the rule is how the guarantee rots -- the
    feature matrix once allowed mixed sampling while common keys refused it --
    so the key-space check and the processing check live together in
    ``samples.py``, and every operation that compares samples calls that one
    function.  Similarity and clustering consume a feature matrix, so they
    inherit its check rather than repeating it.
    """
    with open(os.path.join(ANALYSIS, "samples.py")) as fh:
        text = fh.read()
    assert text.count("def require_common_key_space") == 1
    assert text.count("def check_comparable") == 1
    body = text[text.index("def check_comparable"):]
    body = body[:body.index("\ndef ")]
    assert "require_common_key_space(" in body
    assert "require_common_processing(" in body
    callers = []
    for name, path in _modules():
        if name in ("samples.py", "__init__.py", "render.py"):
            continue
        with open(path) as fh:
            if "check_comparable(" in fh.read():
                callers.append(name)
    assert sorted(callers) == ["common_keys.py", "features.py"]


def test_the_analysis_package_exports_a_stable_surface():
    import biotsr.analysis as analysis

    for name in analysis.__all__:
        assert hasattr(analysis, name), name
    for expected in (
        "AnalysisSample",
        "feature_matrix",
        "common_key_intersection",
        "similarity_matrix",
        "distance_matrix",
        "average_linkage",
        "cluster",
        "clustermap_spec",
    ):
        assert expected in analysis.__all__, expected


def test_importing_the_analysis_layer_does_not_import_the_renderer():
    """``biotsr.analysis`` must work with no matplotlib and no seaborn installed.

    Checked in a fresh interpreter: in this process a render test that ran
    earlier (with the plot extra installed) has imported the renderer on
    purpose, which says nothing about what importing the analysis layer does.
    """
    import subprocess
    import sys

    code = ("import sys, biotsr.analysis; "
            "assert 'biotsr.analysis.render' not in sys.modules; "
            "assert 'matplotlib' not in sys.modules; "
            "assert 'seaborn' not in sys.modules")
    src = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "src")
    env = dict(os.environ, PYTHONPATH=src + os.pathsep + os.environ.get("PYTHONPATH", ""))
    assert subprocess.call([sys.executable, "-c", code], env=env) == 0
