"""One resolved, validated description of one TSR run -- and how to execute it.

``RunSpec`` is the CLI's only model of a calculation.  A command line and a
manifest row both become one of these, so ``biotsr run`` and ``biotsr batch`` cannot
drift apart: there is exactly one place that decides what a set of user inputs
means, and exactly one place that turns it into an API call.

The rule this module is built around: **it validates the things the API cannot
see, and delegates everything else.**  A nonexistent chain, an unparsable residue
selection, an illegal modifier, an empty pool -- all of those already raise, with
a message written for a human, inside the package.  Catching them to re-word them
would produce a second, worse set of messages that drift out of date.  What the
API genuinely cannot catch is an argument it was never handed: ``--entity-name``
passed to ``protein`` is silently irrelevant to ``analyze()``, so it is rejected
here, naming what the method does take.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Dict, Mapping, Optional, Tuple

from .. import api as _api
from ..domain.errors import InvalidRequest
from ..methods.arguments import (
    CROSS_ARGUMENTS,
    SELECTION_ARGUMENTS,
    accepted_arguments,
    flag,
    refuse_unused_arguments,
)
from ..methods.registry import get_method
from .errors import CLIError

__all__ = [
    "RunSpec",
    "RunOutcome",
    "SELECTION_ARGUMENTS",
    "CROSS_ARGUMENTS",
    "MODIFIER_ARGUMENTS",
    "build_spec",
    "execute",
    "execute_streaming",
    "can_stream_triplets",
    "naming_request",
    "output_prefix",
    "sweeps_residues",
]

#: ``SELECTION_ARGUMENTS`` and ``CROSS_ARGUMENTS`` are re-exported from
#: ``biotsr.methods.arguments``, the one table the API and the CLI both read.

#: The three modifiers plus the annotation, as the registry names them.
MODIFIER_ARGUMENTS: Tuple[str, ...] = ("aa_grouping", "mirror", "size_filter", "sse")

def naming_request(spec: "RunSpec") -> Mapping[str, Any]:
    """The request :func:`biotsr.io.naming.request_prefix` names a calculation from."""
    return {
        "structure_id": spec.structure_id or _stem(spec.pdb_path),
        "method": spec.method,
        "chain": spec.chain,
        "chain_2": spec.cross.get("chain_2") if spec.is_cross else None,
        "sweep": sweeps_residues(spec),
        "selection": dict(spec.selection),
        "cross": _named_cross_settings(spec),
        "modifiers": dict(spec.modifiers),
    }


def _named_cross_settings(spec: "RunSpec") -> Mapping[str, Any]:
    """The cross settings a prefix names: a value equal to the method's own
    default (its reference contact cutoff, the symmetric contact filter) is the
    same calculation as leaving it out, so it is named the same way."""
    out = {k: v for k, v in spec.cross.items() if k != "chain_2"}
    if not spec.is_cross:
        return out
    if out.get("contact_cutoff") is not None and float(out["contact_cutoff"]) == \
            _api.default_contact_cutoff(spec.method):
        out.pop("contact_cutoff")
    if out.get("contact_filter") == "both_sides":
        out.pop("contact_filter")
    return out


def output_prefix(spec: "RunSpec") -> str:
    """The calculation's output prefix: its ``sample_id``, or the automatic name."""
    from ..io.naming import request_prefix

    return spec.sample_id or request_prefix(naming_request(spec))


def _stem(path: str) -> str:
    import os

    base = os.path.basename(path)
    return base[:-4] if base.lower().endswith(".pdb") else base


@dataclass(frozen=True)
class RunSpec:
    """Everything one calculation needs, already checked."""

    method: str
    pdb_path: str
    chain: str
    structure_id: Optional[str] = None
    sample_id: str = ""
    label: Optional[str] = None
    keep_triplets: bool = True
    #: ``{"mirror": True, "size_filter": 25.0, ...}`` -- only what the user set.
    modifiers: Mapping[str, Any] = field(default_factory=dict)
    #: Method-specific target arguments -- only what the user set.
    selection: Mapping[str, Any] = field(default_factory=dict)
    #: Cross-only arguments -- only what the user set.
    cross: Mapping[str, Any] = field(default_factory=dict)
    #: Where the row came from, for a batch failure report.  ``0`` for ``biotsr run``.
    source_row: int = 0

    @property
    def is_cross(self) -> bool:
        return get_method(self.method).family == "cross"

    def user_parameters(self) -> Mapping[str, Any]:
        """Exactly what the user specified -- never a default echoed back.

        Defaults belong to the package (ADR-069): recording ``mirror:
        false`` here would make the CLI look like the authority on a value it
        never chose.  What is absent was not asked for.
        """
        out: Dict[str, Any] = {
            "method": self.method,
            "pdb": self.pdb_path,
            "chain": self.chain,
        }
        if self.structure_id:
            out["structure_id"] = self.structure_id
        for source in (self.modifiers, self.selection, self.cross):
            for key in sorted(source):
                out[key] = source[key]
        return out


@dataclass(frozen=True)
class RunOutcome:
    """What one executed ``RunSpec`` produced."""

    spec: RunSpec
    #: ``{name: TSRResult}``.  One entry for a pooled method; one per admitted
    #: residue for the per-entity sweep.
    results: Mapping[str, Any]
    #: The profile that was actually used, read off the result rather than off
    #: the request -- ``Profile.with_modifiers`` decorates the name with the
    #: modifiers it resolved, e.g. ``standard_protein_calpha_v1+mirror``.
    profile: str
    #: The undecorated profile: the method's standard profile.
    base_profile: str = ""
    #: Warnings about the calculation as a whole rather than one result: the
    #: modified residues a per-residue sweep leaves out (ADR-101).
    warnings: Tuple[str, ...] = ()


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------


def _method_or_raise(name: str):
    if not name:
        raise CLIError(
            "no method given",
            "pass --method NAME; `biotsr methods` lists them",
        )
    try:
        return get_method(name)
    except KeyError:
        from ..api import list_methods

        raise CLIError(
            "unknown method %r" % name,
            "supported methods: %s" % ", ".join(list_methods()),
        )


def build_spec(
    values: Mapping[str, Any],
    *,
    source_row: int = 0,
    where: str = "",
) -> RunSpec:
    """Validate one set of user inputs and freeze it into a :class:`RunSpec`.

    ``values`` is a flat mapping with ``None`` for anything the user did not
    specify -- the shape ``argparse`` produces and the shape a manifest row is
    normalised into, so both go through this one function.
    """
    prefix = "%s: " % where if where else ""
    method_name = values.get("method") or ""
    recipe = _method_or_raise(method_name)
    is_cross = recipe.family == "cross"

    pdb_path = values.get("pdb")
    if not pdb_path:
        raise CLIError(
            prefix + "no structure given",
            "pass --pdb PATH (or, in a manifest, a pdb_path column)",
        )
    chain = values.get("chain")
    if not chain:
        raise CLIError(
            prefix + "no chain given",
            "pass --chain X; every TSR method is scoped to a chain and the "
            "package never guesses one",
        )

    modifiers: Dict[str, Any] = {}
    for name in MODIFIER_ARGUMENTS:
        value = values.get(name)
        if value is None or value is False:
            continue
        if name not in recipe.supported_modifiers:
            raise CLIError(
                prefix + "method %r does not support %r"
                % (method_name, name.replace("_", "-")),
                "it supports: %s"
                % (
                    ", ".join(
                        "--" + m.replace("_", "-")
                        for m in sorted(recipe.supported_modifiers)
                    )
                    or "no modifiers"
                ),
            )
        modifiers[name] = value

    if "size_filter" in modifiers and not _positive_finite(modifiers["size_filter"]):
        raise CLIError(
            prefix + "--size-filter must be a finite number greater than 0 (got %r)"
            % modifiers["size_filter"],
            "it is a maximum triangle edge length in angstroms; omit it to "
            "apply no size filter at all",
        )

    accepted = accepted_arguments(recipe)
    try:
        refuse_unused_arguments(recipe, values, cli=True)
    except InvalidRequest as exc:
        raise CLIError(prefix + str(exc), "`biotsr methods %s` describes it"
                       % method_name)

    selection: Dict[str, Any] = {}
    cross: Dict[str, Any] = {}
    target = cross if is_cross else selection
    for name in accepted:
        value = values.get(name)
        if value is None or value is False:
            continue
        target[name] = value

    if is_cross:
        if cross.get("contact_cutoff") is not None and not _positive_finite(
                cross["contact_cutoff"]):
            raise CLIError(
                prefix + "--contact-cutoff must be a finite number greater than 0 "
                "(got %r)" % cross["contact_cutoff"],
                "it is the distance in angstroms within which an atom is "
                "admitted to the contact pool; any positive value is accepted",
            )
        if cross.get("no_contact_filter") and cross.get("contact_cutoff") is not None:
            raise CLIError(
                prefix + "--no-contact-filter and --contact-cutoff contradict "
                "each other",
                "--no-contact-filter skips the contact scan entirely, so there "
                "is no cutoff to apply",
            )
        if cross.get("no_contact_filter") and cross.get("contact_filter"):
            raise CLIError(
                prefix + "--no-contact-filter and --contact-filter contradict "
                "each other",
                "--contact-filter chooses which side the contact scan narrows, "
                "and --no-contact-filter asks for no scan at all",
            )
        if not cross.get("chain_2"):
            raise CLIError(
                prefix + "cross method %r needs --chain-2" % method_name,
                "both chains of a cross calculation are named explicitly; the "
                "package never defaults the second one",
            )
    else:
        _, required = SELECTION_ARGUMENTS.get(recipe.selector, ((), ()))
        if selection.get("residue") and not selection.get("residue_name"):
            raise CLIError(
                prefix + "method %r was given --residue without --residue-name"
                % method_name,
                "name one residue with both --residue and --residue-name; "
                "--residue-name alone runs every residue of that name, and "
                "neither runs every residue of the chain",
            )
        missing = [
            n for n in required
            if n not in selection and not _sweeps_entities(recipe, values)
        ]
        if missing:
            raise CLIError(
                prefix + "method %r needs %s"
                % (method_name, " and ".join(flag(m) for m in missing)),
                _selection_hint(recipe),
            )

    # Residue identifiers are PDB numbers with an optional insertion code;
    # normalised here (``"007"`` -> ``"7"``) so the request, its output name
    # and its record agree, and refused here when they are not identifiers.
    from ..api import parse_residue_id

    for store in (selection, cross):
        for name in ("residue", "entity_id"):
            if store.get(name) not in (None, ""):
                try:
                    store[name] = parse_residue_id(store[name], flag(name))
                except InvalidRequest as exc:
                    raise CLIError(prefix + str(exc))
    for name in ("residue_name", "entity_name"):
        for store in (selection, cross):
            if store.get(name):
                store[name] = str(store[name]).strip().upper()

    sample_id = values.get("sample_id") or ""
    if sample_id:
        from .output import validate_prefix

        validate_prefix(str(sample_id), where=where)
    return RunSpec(
        method=method_name,
        pdb_path=str(pdb_path),
        chain=str(chain),
        structure_id=values.get("structure_id") or None,
        sample_id=str(sample_id),
        label=values.get("label") or None,
        keep_triplets=bool(values.get("keep_triplets", True)),
        modifiers=modifiers,
        selection=selection,
        cross=cross,
        source_row=source_row,
    )


def _positive_finite(value) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number > 0 and number != float("inf")


def _sweeps_entities(recipe, values) -> bool:
    """A per-residue method with no residue identifier runs one calculation per residue.

    ``amino_acid_atomic`` is the package's only ``scope="per_entity"`` method and
    ``analyze_entities()`` is the API call for exactly this: every residue of
    the chain, or with ``--residue-name`` every residue of that name.  Making
    the sweep the behaviour of an omitted ``--residue`` is not a silent choice
    of a scientific argument: it is the difference between "this residue" and
    "each of them", stated in the output prefix (``perresidue``) and the
    completion record.
    """
    return recipe.scope == "per_entity" and not values.get("residue")


def sweeps_residues(spec: "RunSpec") -> bool:
    """Whether ``spec`` is the per-residue sweep of ``amino_acid_atomic``."""
    return not spec.is_cross and _sweeps_entities(get_method(spec.method), spec.selection)


def _selection_hint(recipe) -> str:
    return {
        "protein_calpha_selection": "e.g. --residues '10-12;48;70-71'",
        "amino_acid_residue_atoms": (
            "e.g. --residue 15 --residue-name LEU; --residue-name LEU alone runs "
            "every leucine, and neither runs every residue of the chain"
        ),
        "amino_acid_atom_selection": "e.g. --atoms 'LYS_48_NZ;HIS_68_NE2'",
        "residue_type_atoms": "e.g. --residue-types 'TYR' or 'TYR;PHE'",
        "specific_residue_atoms": "e.g. --residues 'LYS_6;GLU18'",
        "entity_atoms": (
            "e.g. --entity-name ESW --entity-id 404 (omit --entity-id for every "
            "copy on the chain)"
        ),
        "entity_atom_selection": (
            "e.g. --entity-name ESW --entity-id 404 --atoms 'N4;C5;O6'"
        ),
    }.get(recipe.selector, "")


# ---------------------------------------------------------------------------
# Executing
# ---------------------------------------------------------------------------


def _require_points(spec: RunSpec, results) -> None:
    """Refuse to write a calculation that selected nothing.

    ``analyze()`` is right to return a zero-point result: an empty selection is a
    fact about the structure, the reason is recorded on the result, and a library
    must not decide for its caller that a fact is an error.  A command line is a
    different contract.  ``biotsr run --chain Z`` on a structure with no chain Z
    would otherwise write a directory of empty files and exit 0, and the next
    thing to read that directory has no way to tell it from a real result.

    So the CLI says so, names the chains that do exist, and exits 2.  Nothing is
    written, because the guard runs before the output directory is committed.
    """
    if any(r.counts.n_points for r in results.values()):
        return
    target = ", ".join(
        "%s=%r" % (k, v) for k, v in sorted(spec.selection.items())
    ) or ", ".join("%s=%r" % (k, v) for k, v in sorted(spec.cross.items()))
    skips = [s for r in results.values() for s in r.skipped]
    declined = [s.detail for s in skips if s.stage == "completeness" and s.detail]
    others = [s for s in skips if s.stage != "completeness"]
    if declined:
        # The target exists; the method's completeness rule declined it. Naming
        # the chains would send the user looking for a typo that is not there.
        # Atoms an admission rule dropped beforehand are counted, not hidden.
        raise CLIError(
            "nothing to compute: method %r, chain %r%s: the method's completeness "
            "rule excluded the selection (%s)%s"
            % (spec.method, spec.chain, (", " + target) if target else "",
               "; ".join(declined),
               "; %d other record(s) were not admitted (%s)"
               % (len(others), ", ".join(sorted({s.stage for s in others})))
               if others else ""),
            "this method needs every atom its rule lists; leave this target out "
            "of the dataset, or use an *_atom_selection method that names the "
            "atoms present",
        )
    empty = [w for r in results.values() for w in r.warnings
             if w.startswith("no triangles:") and "within" in w]
    if empty:
        # A cross whose sides both selected atoms, none of them in contact:
        # the selection was not empty, the contact pool was.
        raise CLIError(
            "nothing to compute: method %r, chain %r%s: %s"
            % (spec.method, spec.chain, (", " + target) if target else "",
               empty[0][len("no triangles: "):]),
            "check the two chains, or widen --contact-cutoff",
        )
    available = ""
    try:
        # Read the file a second time, only now, and only to name the chains that
        # do exist. Chain labels do not depend on coordinate precision, so this
        # cannot affect any number -- there are no numbers left to affect.
        parsed = _api.load_structure(spec.pdb_path)
        available = ", ".join(c.value for c in parsed.chains())
    except Exception:  # pragma: no cover -- diagnostics must never mask the error
        available = ""
    raise CLIError(
        "selection produced no atoms: method %r, chain %r%s"
        % (spec.method, spec.chain, (", " + target) if target else ""),
        (
            "chains in this structure: %s" % available
            if available
            else "check the chain and the target selection against the structure"
        ),
    )


def can_stream_triplets(spec: RunSpec) -> bool:
    """Whether :func:`execute_streaming` can hand this run's triplets to a sink.

    Every single-result method can; the per-residue sweep of
    ``amino_acid_atomic`` produces one small result per residue and keeps its
    records in memory as before.
    """
    return spec.is_cross or not _sweeps_entities(get_method(spec.method), spec.selection)


def execute_streaming(spec: RunSpec, triplet_sink) -> RunOutcome:
    """:func:`execute`, with every triplet record passed to ``triplet_sink`` as it
    is produced instead of being kept on the result.

    The records, their order and every other field of the result are exactly
    those :func:`execute` produces; only where the records go differs.  This is
    what keeps a large chain's triplet file from being held in memory before it
    is written.  Only for runs where :func:`can_stream_triplets` is true.
    """
    if not can_stream_triplets(spec):
        raise ValueError("this run cannot stream its triplets: %s" % spec.method)
    return _execute(spec, triplet_sink)


def execute(spec: RunSpec) -> RunOutcome:
    """Call the public API once (or once per entity) and return what it produced.

    Everything scientific happens inside ``analyze``/``analyze_cross``/
    ``analyze_entities``.  This function chooses which of the three to call and
    passes the arguments through; it computes nothing.  It hands the API a
    path, never a pre-parsed structure, so how the structure is read is decided
    in one place, by the package.
    """
    return _execute(spec, None)


def _outcome(spec: RunSpec, results) -> RunOutcome:
    first = next(iter(results.values()))
    return RunOutcome(
        spec=spec,
        results=results,
        profile=first.annotations.get("profile", ""),
        base_profile=first.annotations.get("base_profile", ""),
    )


def _execute(spec: RunSpec, triplet_sink) -> RunOutcome:
    streaming = {}
    if triplet_sink is not None:
        streaming = {"keep_triplets": False, "triplet_sink": triplet_sink}
    source = spec.pdb_path
    recipe = get_method(spec.method)

    if spec.is_cross:
        result = _api.analyze_cross(
            source,
            method=spec.method,
            chain=spec.chain,
            structure_id=spec.structure_id,
            **dict(dict(spec.cross, **spec.modifiers),
                   **(streaming or {"keep_triplets": spec.keep_triplets}))
        )
        _require_points(spec, {"": result})
        return _outcome(spec, {"": result})

    kwargs = dict(spec.selection)
    kwargs.update(spec.modifiers)

    if _sweeps_entities(recipe, spec.selection):
        results = _api.analyze_entities(
            source,
            method=spec.method,
            chain=spec.chain,
            structure_id=spec.structure_id,
            keep_triplets=spec.keep_triplets,
            **kwargs
        )
        if not results:
            raise CLIError(
                "no residue%s of chain %r was admitted by %r"
                % (" named %s" % kwargs["residue_name"] if kwargs.get("residue_name")
                   else "", spec.chain, spec.method),
                "the completeness policy declined every candidate; run with a "
                "named --residue to see the recorded reason",
            )
        _require_points(spec, results)
        outcome = _outcome(spec, results)
        if not kwargs.get("residue_name"):
            left_out = _api.unsupported_residues(source, chain=spec.chain,
                                                 method=spec.method)
            if left_out:
                outcome = replace(outcome, warnings=(
                    "not supported: %s (%s); not part of this per-residue sweep"
                    % (", ".join(left_out), _api.UNSUPPORTED_MODIFIED),))
        return outcome

    kwargs.update(streaming or {"keep_triplets": spec.keep_triplets})
    result = _api.analyze(
        source,
        method=spec.method,
        chain=spec.chain,
        structure_id=spec.structure_id,
        **kwargs
    )
    _require_points(spec, {"": result})
    return _outcome(spec, {"": result})
