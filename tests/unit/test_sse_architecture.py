"""Structural guards on the SSE layer.

The whole claim of the SSE layer is "SSE is an annotation, not a second engine".
That is
cheap to say and easy to lose: one ``from ..geometry`` import, one copied
canonical-ordering branch, one private key formula, and there are two engines
again -- which is exactly how the legacy family ended up with
``SSE_TSR.py`` as a 452-line copy of ``TSR.py`` that duplicates its distance
function, its ordering, its binners, its key arithmetic and its PDB parser, and
then drifts (the copy's point-set rule silently drops insertion-code residues).

These tests read the source and refuse that outcome.
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
SSE = os.path.join(SRC, "sse")

#: The SSE layer may reach into ``domain`` (identities) and ``structure`` (the
#: HELIX/SHEET records it interprets). Nothing else.
ALLOWED_PACKAGE_IMPORTS = {"domain", "structure"}

#: Every engine-side concern the legacy copy duplicated.
FORBIDDEN_PACKAGE_IMPORTS = (
    "geometry",
    "encoding",
    "binning",
    "candidates",
    "methods",
    "modifiers",
    "profiles",
    "labels",
    "selection",
    "analysis",
    "cross",
    "io",
)


def _sse_files():
    return [
        os.path.join(SSE, name)
        for name in sorted(os.listdir(SSE))
        if name.endswith(".py")
    ]


def _package_imports(tree):
    """Which *sibling packages of ``biotsr.sse``* a module reaches into.

    ``from .states import ...`` (level 1) is internal to ``biotsr.sse`` and is not a
    reach; ``from ..geometry import ...`` (level 2) is.
    """
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.level >= 2 and node.module:
                found.add(node.module.split(".")[0])
            elif not node.level and node.module and node.module.startswith("biotsr."):
                found.add(node.module.split(".")[1])
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("biotsr."):
                    found.add(alias.name.split(".")[1])
    return found


def _code_tokens(path):
    """Every token of a file that is not a string literal or a comment."""
    import tokenize

    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(fh.readline):
            if tok.type in (tokenize.STRING, tokenize.COMMENT):
                continue
            yield tok


def test_the_sse_layer_imports_nothing_from_the_engine():
    offences = []
    for path in _sse_files():
        imports = _package_imports(ast.parse(open(path).read()))
        for name in sorted(imports & set(FORBIDDEN_PACKAGE_IMPORTS)):
            offences.append((os.path.basename(path), name))
    assert not offences, (
        "biotsr.sse must not reach into the engine; found %r. An annotation that "
        "imports geometry or encoding is on its way to becoming a second engine."
        % (offences,)
    )


def test_the_sse_layer_imports_only_what_it_is_allowed_to():
    for path in _sse_files():
        imports = _package_imports(ast.parse(open(path).read()))
        # `sse` itself is the package's own internal cross-imports.
        extra = imports - ALLOWED_PACKAGE_IMPORTS - {"sse"}
        assert not extra, (os.path.basename(path), sorted(extra))


def test_the_sse_layer_does_no_geometry():
    """No distances, no angles, no square roots, no pi -- in executable code.

    Prose in a docstring may name them; that is what the tokenizer filter is for.
    """
    banned = {"acos", "sqrt", "pi", "hypot", "degrees", "radians"}
    offences = []
    for path in _sse_files():
        for tok in _code_tokens(path):
            if tok.string in banned:
                offences.append((os.path.basename(path), tok.string, tok.start[0]))
            if tok.string == "3.14":
                offences.append((os.path.basename(path), tok.string, tok.start[0]))
    assert not offences, (
        "geometry appeared in biotsr.sse: %r. The annotation layer must read "
        "identities, never coordinates." % (offences,)
    )


def test_the_sse_layer_never_touches_a_coordinate():
    """It is handed ``AtomIdentity``; ``coord``/``Point`` must not appear."""
    offences = []
    for path in _sse_files():
        for tok in _code_tokens(path):
            if tok.string in ("coord", "Point", "coordinates", "squared_distance_to"):
                offences.append((os.path.basename(path), tok.string))
    assert not offences, offences


def test_the_sse_layer_contains_no_key_arithmetic():
    """The key formula's constants must not appear anywhere in executable code."""
    import tokenize

    banned_numbers = {29, 35, 20, 13, 112}
    offences = []
    for path in _sse_files():
        for tok in _code_tokens(path):
            if tok.type != tokenize.NUMBER:
                continue
            try:
                value = int(tok.string)
            except ValueError:
                continue
            if value in banned_numbers:
                offences.append((os.path.basename(path), tok.string, tok.start[0]))
    assert not offences, (
        "a radix or bin count appears in biotsr.sse: %r. SSE has no key space of its "
        "own and must not know the shape of one." % (offences,)
    )


def test_the_sse_layer_does_no_pdb_parsing():
    """It consumes ``SecondaryStructureRecord``; it never slices a line."""
    offences = []
    for path in _sse_files():
        for tok in _code_tokens(path):
            if tok.string in ("open", "readlines", "HETATM", "ATOM"):
                offences.append((os.path.basename(path), tok.string))
    assert not offences, offences


def test_the_engine_consults_the_annotator_after_the_key_is_final():
    """The ANNOTATION hook must sit below the KEY_SIGN hook, in source order.

    This is the structural reason key invariance holds. If the annotator were
    consulted earlier, nothing would stop a future edit from letting it influence
    the ordering, the bin or the sign.
    """
    text = open(os.path.join(SRC, "methods", "engine.py")).read()
    key_line = text.index("key = mirror.apply(")
    annotate_line = text.index("sse_annotator.annotate(")
    assert key_line < annotate_line


def test_the_engine_hands_the_annotator_identities_and_nothing_else():
    """Reading the call site: it receives ``ordered_identities``, full stop."""
    tree = ast.parse(open(os.path.join(SRC, "methods", "engine.py")).read())
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "annotate"
    ]
    assert len(calls) == 1
    call = calls[0]
    assert not call.keywords
    assert len(call.args) == 1
    assert isinstance(call.args[0], ast.Name)
    assert call.args[0].id == "ordered_identities"


def test_there_is_no_second_engine_loop():
    """No module outside ``methods/engine.py`` runs the ordering/geometry pipeline."""
    io_dir = os.path.join(SRC, "io")
    io_files = sorted(
        os.path.join(io_dir, name) for name in os.listdir(io_dir) if name.endswith(".py")
    )
    for path in _sse_files() + io_files:
        text = open(path).read()
        for forbidden in ("canonical_order", "triangle_geometry", "exhaustive("):
            assert forbidden not in text, (os.path.basename(path), forbidden)


def test_sse_adds_no_field_to_either_identity():
    """``KeySpace`` must not learn about secondary structure.

    ``SamplingProfile`` gained ``seq_label_policy`` with the SSE layer, which is NOT
    an SSE field: it is the point-set rule ``SSE_TSR.py`` applies (the rule that
    drops insertion-code residues), and it
    is a sampling fact whether or not anything is annotated.
    """
    from biotsr.domain.keyspace import KeySpace, SamplingProfile
    import dataclasses

    key_fields = {f.name for f in dataclasses.fields(KeySpace)}
    sampling_fields = {f.name for f in dataclasses.fields(SamplingProfile)}
    assert not any("sse" in f for f in key_fields)
    assert not any("sse" in f for f in sampling_fields)
    assert "seq_label_policy" in sampling_fields
