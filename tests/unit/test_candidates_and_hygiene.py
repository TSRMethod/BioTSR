"""Candidate enumeration, plus the package-level hygiene rules.

The hygiene tests are the ones that keep the architecture's promises honest:

* the core imports nothing outside the standard library;
* the legacy angle literal (``3.14``) is absent from the package;
* module dependency directions are not violated;
* no code builds a ``.3Dkeys_*`` name, and no ``legacy_*`` module ships.
"""

from __future__ import annotations

import ast
import math
import os

import pytest

from biotsr.candidates.exhaustive import count, exhaustive

SRC = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "src",
    "biotsr",
)


# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4, 5, 10, 76])
def test_exhaustive_yields_exactly_c_n_3(n):
    produced = list(exhaustive(n))
    assert len(produced) == count(n) == math.comb(n, 3)


def test_enumeration_is_lexicographic_and_ascending():
    got = [c.as_tuple() for c in exhaustive(5)]
    assert got == sorted(got)
    assert got[0] == (0, 1, 2)
    assert got[-1] == (2, 3, 4)
    for i, j, k in got:
        assert i < j < k


def test_enumeration_is_repeatable_and_free_of_duplicates():
    a = [c.as_tuple() for c in exhaustive(20)]
    b = [c.as_tuple() for c in exhaustive(20)]
    assert a == b
    assert len(set(a)) == len(a)


def test_no_set_is_built():
    """Order must be fixed by construction, never by hash traversal -- the direct
    fix for the confirmed nondeterminism in the legacy frequency-table builder."""
    source = open(os.path.join(SRC, "candidates", "exhaustive.py")).read()
    tree = ast.parse(source)
    for node in ast.walk(tree):
        assert not isinstance(node, ast.SetComp)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            assert node.func.id != "set"


def test_negative_input_is_rejected():
    with pytest.raises(ValueError):
        list(exhaustive(-1))


# ---------------------------------------------------------------------------
# Hygiene
# ---------------------------------------------------------------------------


def _python_files():
    for root, _dirs, files in os.walk(SRC):
        for name in sorted(files):
            if name.endswith(".py"):
                yield os.path.join(root, name)


#: The standard-library modules the package may import. An allow-list rather
#: than a deny-list, so a new third-party dependency cannot arrive unnoticed --
#: adding to it is a deliberate edit with a reason, which is the point.
#:
#: ``shutil`` and ``traceback`` were added for the CLI:
#: ``shutil.rmtree`` clears a staging directory when a write is abandoned, and
#: ``traceback.print_exc`` is what ``--debug`` prints for a genuinely internal
#: failure. Both are standard library on every supported interpreter and neither
#: is reachable from key generation.
STDLIB_ALLOWED = {
    "__future__",
    "argparse",
    "ast",
    "bisect",
    "collections",
    # `copy` is stdlib and is used in exactly one place: `TSRResult.metadata()`
    # deep-copies the annotation blocks it embeds, so a caller mutating the
    # returned payload cannot reach into the result it came from.
    "copy",
    "csv",
    "dataclasses",
    # the "did you mean" suggestion for a misspelt manifest column
    "difflib",
    # telling "a sample is already there" from other failures of a rename
    "errno",
    "functools",
    "gzip",
    "hashlib",
    "importlib",
    "io",
    "itertools",
    "json",
    "math",
    "os",
    "re",
    # a staging-directory token unique to one invocation
    "secrets",
    "shutil",
    # a batch treats SIGTERM (a scheduler's time limit) like Ctrl-C
    "signal",
    "struct",
    "sys",
    "tempfile",
    "traceback",
    "typing",
    "zipfile",
}

#: The one module that may open a network connection, and what it may use for
#: it. Retrieval is opt-in; nothing else in the package touches the network.
RETRIEVAL_MODULE = "retrieval.py"
RETRIEVAL_ALLOWED = {
    "concurrent", "datetime", "email", "errno", "functools", "http", "random",
    "selectors", "socket", "ssl", "threading", "time", "urllib", "zlib",
}

#: The one module that may import joblib, the parallel backend: a declared
#: runtime dependency (ADR-086) that ``--jobs 1`` never needs.
PARALLEL_MODULE = os.path.join("cli", "parallel.py")
PARALLEL_ALLOWED = {"joblib"}

#: The CLI output module reads the host name for the staging-directory token
#: (``<pid>-<hosttag8>-<hex8>``), falling back to ``socket.gethostname()`` where
#: ``os.uname`` does not exist. Only that call is allowed (checked below); it
#: opens no connection.
CLI_OUTPUT_MODULE = os.path.join("cli", "output.py")
CLI_OUTPUT_ALLOWED = {"socket"}


#: The single sanctioned third-party importer, and the only libraries it may use.
#: Drawing a figure needs matplotlib, the ``plot`` extra (ADR-059). Nothing
#: else in the package may import a third-party module, and nothing else may
#: import this module at all.
RENDER_MODULE = os.path.join("analysis", "render.py")
RENDER_ALLOWED = {"matplotlib"}


def test_core_imports_nothing_outside_the_standard_library():
    """Key generation *and the whole numerical analysis layer* must run where
    pandas/scipy/TensorFlow/MPI do not."""
    offenders = []
    extra = {
        RETRIEVAL_MODULE: RETRIEVAL_ALLOWED,
        PARALLEL_MODULE: PARALLEL_ALLOWED,
        CLI_OUTPUT_MODULE: CLI_OUTPUT_ALLOWED,
    }
    for path in _python_files():
        rel = os.path.relpath(path, SRC)
        if rel == RENDER_MODULE:
            continue
        allowed = STDLIB_ALLOWED | extra.get(rel, set())
        tree = ast.parse(open(path).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    top = alias.name.split(".")[0]
                    if top not in allowed and top != "biotsr":
                        offenders.append((path, alias.name))
            elif isinstance(node, ast.ImportFrom):
                if node.level:  # relative import inside the package
                    continue
                top = (node.module or "").split(".")[0]
                if top not in allowed and top != "biotsr":
                    offenders.append((path, node.module))
    assert offenders == [], "third-party imports found in the core: %r" % offenders


def test_the_cli_output_module_uses_socket_only_for_the_host_name():
    """The ``socket`` carve-out is for ``gethostname`` and nothing else."""
    tree = ast.parse(open(os.path.join(SRC, CLI_OUTPUT_MODULE)).read())
    used = {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "socket"
    }
    assert used == {"gethostname"}, used
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert (node.module or "").split(".")[0] != "socket"


def test_joblib_is_imported_only_inside_a_function_of_the_parallel_module():
    """``import biotsr`` and ``--jobs 1`` must work even where joblib is missing."""
    path = os.path.join(SRC, PARALLEL_MODULE)
    tree = ast.parse(open(path).read())
    for node in tree.body:  # module scope only
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [a.name for a in node.names] + [getattr(node, "module", "") or ""]
            assert not any(n.split(".")[0] == "joblib" for n in names), (
                "joblib imported at module scope in %s" % PARALLEL_MODULE)


def _docstrings(tree):
    """Every docstring node in a module, so prose can be excluded from code checks."""
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc is not None and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr):
                    out.add(id(first.value))
    return out


def test_the_legacy_angle_literal_is_absent_from_the_package():
    """`3.14` must not be reachable at all: the package ships the standard only.

    Checked against the AST, not the raw text: a module may still *describe* the
    historical divisor in prose, but no module may contain it as an executable
    float literal.
    """
    hits = []
    for path in _python_files():
        tree = ast.parse(open(path).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value.__class__ is float:
                if node.value == 3.14:
                    hits.append(os.path.relpath(path, SRC))
    assert sorted(set(hits)) == [], (
        "the legacy angle literal must not appear as a float literal anywhere in "
        "the package; found in %r" % sorted(set(hits))
    )


def test_no_module_builds_a_legacy_angle_convention():
    """A second, independent guard: the legacy convention is gone everywhere."""
    import biotsr.geometry.conventions as conventions

    hits = []
    for path in _python_files():
        text = open(path).read()
        if "LEGACY_ANGLE_CONVENTION" in text:
            hits.append(os.path.relpath(path, SRC))
    assert hits == []
    assert not hasattr(conventions, "LEGACY_ANGLE_CONVENTION")


def test_no_legacy_module_ships_with_the_package():
    """The quarantined legacy writers and profiles were removed, not hidden."""
    import importlib

    names = []
    for root, _dirs, files in os.walk(SRC):
        for name in files:
            if name.endswith(".py") and name.startswith("legacy"):
                names.append(os.path.relpath(os.path.join(root, name), SRC))
    assert names == [], names
    for module in (
        "biotsr.io.legacy_text",
        "biotsr.io.legacy_atomic",
        "biotsr.io.legacy_cross",
        "biotsr.io.legacy_sse",
        "biotsr.profiles.legacy",
        "biotsr.profiles.catalogue",
    ):
        with pytest.raises(ImportError):
            importlib.import_module(module)


def test_domain_is_a_dependency_sink():
    """`domain` may import the standard library and itself, nothing else in `biotsr`."""
    domain_dir = os.path.join(SRC, "domain")
    for path in sorted(os.listdir(domain_dir)):
        if not path.endswith(".py"):
            continue
        tree = ast.parse(open(os.path.join(domain_dir, path)).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level:
                assert node.level == 1, (
                    "%s reaches outside biotsr.domain (relative level %d, module %r)"
                    % (path, node.level, node.module)
                )


def test_geometry_never_imports_the_structure_layer():
    """The geometry engine never sees a PDB line."""
    forbidden = ("structure", "selection", "io", "methods", "binning", "encoding")
    geometry_dir = os.path.join(SRC, "geometry")
    for name in sorted(os.listdir(geometry_dir)):
        if not name.endswith(".py"):
            continue
        text = open(os.path.join(geometry_dir, name)).read()
        for module in forbidden:
            assert ("from ..%s" % module) not in text, "%s imports %s" % (name, module)


def test_encoding_never_imports_geometry_or_binning():
    text = open(os.path.join(SRC, "encoding", "key.py")).read()
    assert "from ..geometry" not in text
    assert "from ..binning" not in text
    assert "from ..structure" not in text
    assert "from ..selection" not in text


def test_no_code_builds_a_3dkeys_name():
    """``.3Dkeys_*`` must not be reachable at all.

    AST-based for the same reason as the angle-literal check: modules may mention
    the historical name in prose (docstrings), but no string constant in code may
    contain it. The legacy writers that built it are gone.
    """
    hits = []
    for path in _python_files():
        tree = ast.parse(open(path).read())
        docs = _docstrings(tree)
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and "3Dkeys" in node.value
                and id(node) not in docs
            ):
                hits.append(os.path.relpath(path, SRC))
    assert sorted(set(hits)) == []


def test_the_native_writer_names_files_by_prefix_not_by_a_legacy_extension():
    text = open(os.path.join(SRC, "io", "native.py")).read()
    assert 'KEYS_SUFFIX = ".keys"' in text
    assert "3Dkeys" not in text.split('"""', 2)[-1], (
        "the native writer must not emit a .3Dkeys_* name anywhere in its code"
    )


def test_packaged_resources_are_never_loaded_from_the_working_directory():
    """The label table must come from the package, not from wherever the process sits.

    The legacy amino-acid analyzer calls
    ``readDrugLexicalCsv("drug_atom_lexical_txt.csv")`` -- a bare relative path --
    so it only runs if the process has already ``cd``-ed into the package
    directory. Three repositories share that defect. This asserts the new loader
    resolves the resource against the package instead, by loading it from an
    unrelated working directory and checking it still works.
    """
    import os
    import subprocess
    import sys
    import tempfile

    src = os.path.dirname(SRC)
    script = (
        "import sys; sys.path.insert(0, %r);"
        "from biotsr.labels.elements import ELEMENT_LABELS_112 as L;"
        "print(L.n_labels, L.n_distinct_codes, L.resolve_atom('CA', 'C'))" % src
    )
    with tempfile.TemporaryDirectory() as elsewhere:
        proc = subprocess.run(
            [sys.executable, "-c", script], cwd=elsewhere,
            capture_output=True, text=True,
        )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.split() == ["112", "109", "106"]


def test_no_module_reads_a_table_by_bare_relative_path():
    """A second, static guard: no ``open("...csv")`` on a bare filename."""
    offenders = []
    for path in _python_files():
        tree = ast.parse(open(path).read())
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)):
                continue
            if node.func.id != "open" or not node.args:
                continue
            first = node.args[0]
            if isinstance(first, ast.Constant) and isinstance(first.value, str):
                if first.value.endswith(".csv") and os.sep not in first.value:
                    offenders.append((os.path.relpath(path, SRC), first.value))
    assert offenders == [], "bare relative resource path(s): %r" % offenders


def test_the_renderer_is_the_only_third_party_importer_and_is_narrowly_scoped():
    """``analysis/render.py`` may import matplotlib -- nothing else."""
    path = os.path.join(SRC, RENDER_MODULE)
    tree = ast.parse(open(path).read())
    offenders = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name.split(".")[0]
                if top not in STDLIB_ALLOWED and top not in RENDER_ALLOWED:
                    offenders.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                continue
            top = (node.module or "").split(".")[0]
            if top not in STDLIB_ALLOWED and top not in RENDER_ALLOWED:
                offenders.append(node.module)
    assert offenders == [], "render.py reaches beyond its carve-out: %r" % offenders


def test_nothing_imports_the_renderer():
    """``import biotsr`` and every analysis computation work with no plotting stack.

    Checked structurally rather than by uninstalling matplotlib: no other module
    in the package may name ``render`` in an import at all, at module scope or
    inside a function.
    """
    offenders = []
    for path in _python_files():
        rel = os.path.relpath(path, SRC)
        if rel == RENDER_MODULE:
            continue
        tree = ast.parse(open(path).read())
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and (node.module or "").endswith("render"):
                offenders.append((rel, node.module))
            elif isinstance(node, ast.ImportFrom) and node.level:
                for alias in node.names:
                    if alias.name == "render":
                        offenders.append((rel, alias.name))
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.split(".")[-1] == "render":
                        offenders.append((rel, alias.name))
    assert offenders == [], "the renderer must stay optional: %r" % offenders
