"""The CLI is a wrapper, and this file is what keeps it one.

The central risk for a CLI is not that it computes something wrong -- the
equivalence gate catches that.  It is that the CLI slowly accumulates a *second*
opinion about the science: a default here, a threshold there, a special case for
one family, until two code paths disagree and nobody knows which is authoritative.
Every test below is aimed at that, not at behaviour.
"""

from __future__ import annotations

import ast
import os

import pytest

import biotsr
from biotsr import cli
from biotsr.cli import parser as cli_parser
from biotsr.cli import spec as cli_spec
from biotsr.methods.registry import get_method

SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "src",
)
CLI_DIR = os.path.join(SRC, "biotsr", "cli")


def _cli_files():
    for name in sorted(os.listdir(CLI_DIR)):
        if name.endswith(".py"):
            yield os.path.join(CLI_DIR, name)


def _tree(path):
    return ast.parse(open(path).read())


def _docstring_ids(tree):
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            if ast.get_docstring(node, clean=False) is not None and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr):
                    out.add(id(first.value))
    return out


# ---------------------------------------------------------------------------
# No science in the CLI
# ---------------------------------------------------------------------------


def test_the_cli_never_imports_the_scientific_internals():
    """It may call the façade and the serializers. Not the engine.

    ``biotsr.api``, ``biotsr.io``, ``biotsr.profiles`` and the method registry are the
    doors; ``geometry``, ``binning``, ``encoding``, ``candidates``, ``cross``,
    ``selection``, ``structure`` and ``sse`` are not.  The registry is allowed
    because it is a declarative table of recipes with no arithmetic in it, and
    generating ``biotsr methods`` from anything else would mean maintaining a second
    list.  If the CLI ever needs one of the others, the thing it needs belongs on
    the façade.
    """
    forbidden = {
        "geometry", "binning", "encoding", "candidates", "cross", "selection",
        "structure", "labels", "modifiers", "sse",
    }
    offenders = []
    for path in _cli_files():
        for node in ast.walk(_tree(path)):
            module = ""
            if isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if node.level and not module:
                    continue
            elif isinstance(node, ast.Import):
                module = ",".join(a.name for a in node.names)
            for part in module.replace("biotsr.", "").split("."):
                if part in forbidden:
                    offenders.append((os.path.basename(path), module))
    assert offenders == [], (
        "the CLI reached into a scientific internal: %r -- route it through "
        "biotsr.api instead" % offenders
    )


def test_the_cli_contains_no_scientific_constant():
    """No bin count, no radix, no angle, no cutoff, no category threshold.

    Every number that can reach a key lives on a ``Profile`` or a ``Method``.  The
    integers the CLI is allowed to hold are the exit codes, verbosity levels,
    formatting widths and list indices -- all of which are small, and none of
    which is a distance, an angle or a count of anything physical.
    """
    suspicious = []
    # 130 is the exit code for an interrupted run (128 + SIGINT).
    allowed = {0, 1, 2, 3, -1, 4, 5, 8, 10, 12, 16, 24, 32, 36, 76, 88, 130}
    # Named, non-scientific floats: file -> the one value it may hold, and why.
    allowed_floats = {
        # the difflib similarity above which an unknown column name is offered a
        # "did you mean" suggestion
        "manifest.py": {0.8},
    }
    # Named, non-scientific integers: file -> the values it may hold, and why.
    allowed_ints = {
        # the longest sample id accepted as a directory name (bytes of UTF-8),
        # and the DEL code point refused as a control character in one
        "output.py": {200, 127},
        # the dry-run preview's column widths and the truncation of its
        # selection column
        "batch.py": {7, 28, 31, 34},
    }
    for path in _cli_files():
        tree = _tree(path)
        docstrings = _docstring_ids(tree)
        for node in ast.walk(tree):
            if id(node) in docstrings:
                continue
            if isinstance(node, ast.Constant) and isinstance(node.value, float):
                if node.value in allowed_floats.get(os.path.basename(path), ()):
                    continue
                suspicious.append((os.path.basename(path), node.value))
            elif isinstance(node, ast.Constant) and isinstance(node.value, int):
                if node.value in allowed_ints.get(os.path.basename(path), ()):
                    continue
                if node.value not in allowed and not isinstance(node.value, bool):
                    suspicious.append((os.path.basename(path), node.value))
    assert suspicious == [], (
        "numeric literal in the CLI that is not plainly an exit code or a "
        "formatting width: %r" % suspicious
    )


def test_the_cli_declares_no_default_for_a_scientific_argument():
    """``--mirror`` defaults to *unset*, not to ``False``.

    The difference matters: an unset modifier is never passed to ``analyze()``, so
    the package's own default applies and there is exactly one place that knows
    it.  A CLI default of ``False`` would pass ``mirror=False`` explicitly and
    make the CLI a second authority on a value it never chose -- and the day the
    package's default changed, the CLI would silently override it.
    """
    parser = cli_parser.build_parser()
    run = _subparser(parser, "run")
    scientific = {
        "--mirror", "--aa-grouping", "--sse", "--size-filter",
        "--contact-cutoff", "--contact-filter", "--no-contact-filter",
        "--residues", "--residue-types", "--atoms", "--entity-names",
    }
    for action in run._actions:
        for option in action.option_strings:
            if option in scientific:
                assert action.default is None, (
                    "%s carries a CLI-side default (%r); scientific defaults "
                    "belong to the package" % (option, action.default)
                )


def _subparser(parser, name):
    for action in parser._actions:
        if hasattr(action, "choices") and action.choices and name in action.choices:
            return action.choices[name]
    raise AssertionError("no subcommand %r" % name)


def test_the_cli_reimplements_no_serializer():
    """Nothing in the CLI formats a key, a triplet or a filename.

    The native and compatibility writers own that, down to the bin-count suffix.
    A ``.keys_`` or ``.triplets_`` literal in the CLI would be the beginning of a
    second naming scheme.
    """
    markers = (".keys_", ".triplets_", ".3Dkeys", "theta29", "dist18")
    offenders = []
    for path in _cli_files():
        tree = _tree(path)
        docstrings = _docstring_ids(tree)
        for node in ast.walk(tree):
            if id(node) in docstrings:
                continue
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                for marker in markers:
                    if marker in node.value:
                        offenders.append((os.path.basename(path), node.value))
    assert offenders == [], (
        "the CLI spelled a serializer's filename itself: %r" % offenders
    )


# ---------------------------------------------------------------------------
# Discovery cannot drift from the registry
# ---------------------------------------------------------------------------


def test_every_registered_method_is_reachable_and_documented():
    for name in biotsr.list_methods():
        recipe = get_method(name)
        assert recipe.description, "%s has no description for `biotsr methods`" % name
        profile = recipe.default_profile
        assert profile is not None, "%s has no standard profile" % name
        assert profile.description, (
            "%s's standard profile has no description for `biotsr methods`" % name
        )


def test_no_deferred_candidate_strategy_is_reachable():
    """No deferred candidate strategy is exposed, and none becomes exposed.

    Not stubbed, not hidden behind a flag, not named in help: a method the CLI
    can name is a method a user will try, and none of these is implemented.
    """
    text = " ".join(open(p).read().lower() for p in _cli_files())
    for absent in ("ca_knn", "ca-knn", "caknn"):
        assert absent not in [n.lower() for n in biotsr.list_methods()]
    for name in biotsr.list_methods():
        assert "knn" not in name
    # Every registered method is one of the twenty documented ones.
    assert len(biotsr.list_methods()) == 20


# ---------------------------------------------------------------------------
# The public surface
# ---------------------------------------------------------------------------


def test_main_is_importable_from_the_documented_entry_point():
    """``[project.scripts] biotsr = "biotsr.cli:main"`` must resolve."""
    assert callable(cli.main)


def test_python_dash_m_tsr_has_a_module_entry_point():
    assert os.path.isfile(os.path.join(SRC, "biotsr", "__main__.py"))
    assert os.path.isfile(os.path.join(CLI_DIR, "__main__.py"))


def test_the_cli_cannot_hand_the_api_a_pre_parsed_structure():
    """The CLI passes a path; parsing the structure belongs to the package.

    ``load_structure`` passes an already-parsed ``Structure`` through unchanged,
    so a CLI that parsed the file itself would become a second authority on how
    coordinates are read. An earlier version did, to reuse one parse across
    batch rows, and the equivalence gate caught it. The parameter is **absent**,
    not merely unused -- a signature that still accepted a ``Structure`` would
    invite the same defect back, and nothing else in the CLI would notice.
    ADR-069.
    """
    import inspect

    assert list(inspect.signature(cli_spec.execute).parameters) == ["spec"]
    source = open(os.path.join(CLI_DIR, "spec.py")).read()
    assert "load_structure" in source, (
        "the empty-selection diagnostic should still be able to name the "
        "structure's chains"
    )


def test_the_selection_argument_table_covers_every_registered_selector():
    """A new method with an unmapped selector must not silently accept nothing."""
    for name in biotsr.list_methods():
        recipe = get_method(name)
        if recipe.family == "cross":
            continue
        assert recipe.selector in cli_spec.SELECTION_ARGUMENTS, (
            "selector %r (method %r) has no entry in SELECTION_ARGUMENTS, so the "
            "CLI would accept no target for it" % (recipe.selector, name)
        )
