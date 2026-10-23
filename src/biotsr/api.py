"""The public façade -- the only module most users touch.

Simple things are one call::

    from biotsr import analyze
    result = analyze("tests/data/structures/1UBQ.pdb", method="protein", chain="A")

Modifiers are optional keywords on the same call -- they compose freely, and
there is no per-combination class for any of the eight combinations::

    analyze(path, method="protein", chain="A",
            aa_grouping=True, mirror=True, size_filter=25.0)

Every method runs under its own standard profile -- the accepted scientific
conventions for that method -- and an argument the method does not read is
refused by name rather than ignored.
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple, Union

from .candidates.exhaustive import count as candidate_count
from .domain.identity import StructureId
from .domain.result import (
    SelectionSummary,
    SourceInfo,
    TSRResult,
    TripletRecord,
)
from .domain.errors import InvalidRequest
from .methods import builtin as _builtin_recipes  # noqa: F401  (registers the methods)
from .methods.arguments import refuse_unused_arguments
from .methods.engine import run_core
from .methods.registry import Method, get_method, list_methods
from .labels.schemes import LabelDomain, LabelScheme, UnlabelledPointPolicy
from .profiles.standard import Profile
from .selection.atomic import (
    SCOPE_PER_ENTITY,
    SELECTOR_ATOM_SELECTION,
    SELECTOR_ENTITY_ATOM_SELECTION,
    SELECTOR_ENTITY_ATOMS,
    SELECTOR_ENTITY_CHAIN_ATOMS,
    SELECTOR_RESIDUE_ATOMS,
    SELECTOR_RESIDUE_TYPE_ATOMS,
    SELECTOR_SPECIFIC_RESIDUE_ATOMS,
    entity_targets,
    hetero_amino_acid_residues,
    parse_atom_names,
    parse_atom_selection,
    parse_residue_id,
    parse_residue_refs,
    parse_residue_types,
    refuse_unreadable,
    residue_targets,
    select_atom_selection,
    select_entity_atoms,
    select_residue_atoms,
)
from .selection.base import (
    CALPHA_MODIFIED_RESIDUES,
    AltLocPolicy,
    InputRules,
    ModifiedResiduePolicy,
    altloc_policy_of,
)
from .selection.base import processing_identity as _processing_identity
from .selection import integrity
from .selection.completeness import MODIFIED_AMINO_ACIDS
from .sse.assignment import SSEAnnotator, build_assignment_map
from .selection.protein import (
    REASON_NOT_LABELLED,
    SELECTOR_CALPHA_SELECTION,
    SELECTOR_NAME,
    calpha_exclusions,
    check_range_insertion_codes,
    parse_residue_selection,
    select_calpha,
    select_calpha_residues,
)
from .cross.candidates import STRATEGY_MIXED_SIDE as STRATEGY_MIXED_SIDE_LEX
from .selection.truncation import LEGACY_SINGLE_PASS
from .structure.backends.pdb_legacy import read_pdb_legacy
from .structure.model import Structure

__all__ = [
    "analyze",
    "analyze_cross",
    "analyze_entities",
    "list_methods",
    "list_residue_targets",
    "list_entity_targets",
    "list_cross_methods",
    "load_structure",
    "parse_residue_id",
    "unsupported_residues",
    "processing_identity",
]


def load_structure(
    source: Union[str, Structure],
    *,
    structure_id: Optional[str] = None,
    backend: str = "pdb_legacy",
) -> Structure:
    """Read a structure from a path, or pass one through unchanged.

    Coordinates are kept at full double precision, which every method uses.
    """
    if isinstance(source, Structure):
        return source
    if backend != "pdb_legacy":
        raise InvalidRequest(
            "unknown parser backend %r; BioTSR includes one reader, "
            "'pdb_legacy', and it is the v1.0 standard. A reader using the PDB "
            "specification's atom-name columns ('pdb_strict') is not part of "
            "this release; see the known limitations, "
            "https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md."
            % backend
        )
    sid = StructureId.parse(structure_id) if structure_id else None
    return read_pdb_legacy(source, sid)


def analyze(
    source: Union[str, Structure],
    *,
    method: str = "protein",
    chain: str,
    profile: Optional[Profile] = None,
    structure_id: Optional[str] = None,
    keep_triplets: bool = True,
    triplet_sink: Optional[Callable[[TripletRecord], None]] = None,
    aa_grouping: bool = False,
    mirror: bool = False,
    size_filter: Optional[float] = None,
    sse: bool = False,
    residue: Optional[str] = None,
    residue_name: Optional[str] = None,
    atoms=None,
    entity_name: Optional[str] = None,
    entity_id: Optional[str] = None,
    entity_names=None,
    residue_types=None,
    residues=None,
    **modifiers,
) -> TSRResult:
    """Run one TSR calculation and return a structured result.

    ``modifiers`` is accepted only so unsupported ones fail loudly with a message
    naming what the method *does* support, rather than being silently ignored.

    Pass ``triplet_sink`` together with ``keep_triplets=False`` to stream triplet
    records straight into a serializer.  The result is otherwise identical --
    same keys, same counts, same provenance -- but ``result.triplets`` is ``None``
    because nothing was retained.

    The three modifiers:

    * ``aa_grouping`` -- swap the residue labels (30 codes) for the grouped
      labels (23 codes), both at radix 30. Protein Ca only; rejected in any
      other label domain.
    * ``mirror`` -- negate a key when its pre-fold angle exceeds 90 degrees.
      Default off, and never inferred: with ``mirror=False`` no key is ever
      negative.
    * ``size_filter`` -- keep only triangles whose ``max_distance`` is strictly
      below the cutoff, in angstroms.

    And one annotation, which is not a modifier in the same sense:

    * ``sse`` -- annotate each triangle with its secondary-structure category,
      one of the 18, read from the structure's own ``HELIX``/``SHEET`` records.
      **It does not change a single key.**  The keys, their signs, the candidate
      ledger and both identity fingerprints are bit-identical to the same call
      with ``sse=False``; what changes is the bucket each key's frequency is
      counted in, from the implicit ``"regular"``/``"mirror"`` to the SSE
      category.  Protein Ca only (``protein``, ``aa_selection``): the categories
      are defined over protein residues, so requesting it elsewhere raises
      (ADR-036).  It composes freely with all three modifiers.

    Atomic methods take a target instead of a whole chain:

    * ``amino_acid_atomic`` -- ``residue`` is the residue identifier ("15",
      "52A") and ``residue_name`` the three-letter name; both must match the
      structure.  Use :func:`analyze_entities` to run every residue of a chain,
      or every residue of one name.
    * ``amino_acid_atom_selection`` -- ``atoms`` is a ``RESIDUE_SEQNUM_ATOM``
      list or semicolon-separated string, e.g. ``"LYS_48_NZ;HIS_68_NE2"``.
    * ``aa_ind`` / ``aa_ind_mix`` -- ``residue_types`` names the amino-acid
      type(s), e.g. ``"TYR"`` or ``"TYR;PHE"``.
    * ``aa_ind_specific`` -- ``residues`` names residue instances, e.g.
      ``"LYS_6;GLU18"``.
    * ``aa_selection`` -- ``residues`` is a residue-number selection with ranges,
      e.g. ``"10-12;48;70-71"``.  A **Ca** method: radix 30, 35 distance
      classes.
    * ``drug`` / ``carbohydrate`` / ``nucleotide`` -- ``entity_name`` names the
      molecule and ``entity_id`` its residue identifier (number plus optional
      insertion code), e.g. ``entity_name="ESW", entity_id="404"``: that one
      copy.  Without ``entity_id``, every copy of that name on ``chain`` is used
      in one combined calculation (scope ``all_copies``; triangles may span the
      copies).  The ``*_atom_selection`` variants additionally take
      ``atoms="N4;C5;O6"``, applied to each selected copy.
    * ``carbohydrate_chain`` / ``nucleotide_chain`` -- whole-chain polymer scope;
      ``entity_names`` overrides the method's declared residue vocabulary.
    """
    recipe: Method = get_method(method)
    if recipe.family == "cross":
        raise InvalidRequest(
            "method %r is a cross method: it reads two sides, so it is run with "
            "analyze_cross(source, method=%r, chain=..., chain_2=...), not "
            "analyze()" % (recipe.name, recipe.name))
    _require_standard_profile(recipe, profile)
    refuse_unused_arguments(recipe, dict(
        modifiers, residue=residue, residue_name=residue_name, atoms=atoms,
        entity_name=entity_name, entity_id=entity_id, entity_names=entity_names,
        residue_types=residue_types, residues=residues,
    ))
    _refuse_cross_only(recipe, modifiers)
    requested = dict(modifiers)
    if aa_grouping:
        requested["aa_grouping"] = aa_grouping
    if mirror:
        requested["mirror"] = mirror
    if size_filter is not None:
        requested["size_filter"] = size_filter
    if sse:
        requested["sse"] = sse
    recipe.validate_modifiers(requested)

    base = recipe.default_profile
    active = base.with_modifiers(
        aa_grouping=aa_grouping,
        mirror=mirror,
        size_filter=size_filter,
        method=recipe.name,
    )

    structure = load_structure(source, structure_id=structure_id)
    _require_chain_present(structure, chain)
    selection, target_label, selection_warnings, request, explained = _select(
        recipe, active, structure, chain,
        residue=residue, residue_name=residue_name, atoms=atoms,
        entity_name=entity_name, entity_id=entity_id, entity_names=entity_names,
        residue_types=residue_types, residues=residues,
    )

    # The ANNOTATION hook. Built from the structure's own records and the same
    # chain the points came from; it is handed to the engine and nothing else in
    # the pipeline is told about it. The profile supplies one thing: which PDB
    # fields name an individual beta-strand (ADR-063).
    annotator = None
    sse_map = None
    sse_warnings = ()
    if sse:
        sse_map = build_assignment_map(
            structure, chain, strand_identity=active.sse_strand_identity
        )
        annotator = SSEAnnotator(assignments=sse_map)
        sse_warnings = sse_map.warnings

    core = run_core(
        selection.points,
        active,
        candidate_count=candidate_count(len(selection.points)),
        keep_triplets=keep_triplets,
        triplet_sink=triplet_sink,
        sse_annotator=annotator,
    )

    sampling = _sampling_profile(recipe, active)
    outcome_warnings = _selection_warnings(selection, core.counts, explained)

    return TSRResult(
        method=recipe.name,
        structure_id=structure.structure_id,
        chain=chain,
        keys=core.as_category_counts(),
        key_order=core.key_order,
        triplets=core.triplets,
        key_space=active.key_space(),
        sampling=sampling,
        selection=SelectionSummary(
            selector=selection.selector,
            chains_requested=selection.chains_requested,
            chains_used=selection.chains_used,
            n_points=len(selection.points),
            residues_used=selection.residues_used,
            atoms_missing=selection.atoms_missing,
        ),
        source=SourceInfo(
            structure_id=structure.structure_id,
            path=structure.source_path,
            sha256=structure.source_sha256,
            parser_backend=structure.parser_backend,
            model=1,
            # The policy that actually selected the atoms: the C-alpha selectors
            # apply blank-or-A themselves; an atomic method applies its
            # admission's rule, which for amino_acid_atomic is "not B".
            altloc_policy=(
                AltLocPolicy.BLANK_OR_A
                if recipe.selector in (SELECTOR_NAME, SELECTOR_CALPHA_SELECTION)
                else altloc_policy_of(active.atom_admission)
            ),
            modified_residue_policy=ModifiedResiduePolicy.for_label_domain(
                active.label_scheme.domain),
            input_rules=InputRules.CURRENT,
        ),
        counts=core.counts,
        skipped=tuple(selection.skipped) + tuple(core.skipped),
        warnings=tuple(structure.warnings) + tuple(selection.warnings)
        + tuple(selection_warnings) + tuple(sse_warnings) + outcome_warnings,
        annotations={
            "profile": active.name,
            "base_profile": base.name,
            "legacy_compatibility_profile": active.is_legacy_compatibility,
            "modifiers": active.active_modifiers(),
            "truncated_at_line": selection.truncated_at_line,
            "target": target_label,
            "selection_request": request,
            "sse": dict(sse_map.to_dict()) if sse_map is not None else None,
            "sse_coverage": (
                dict(
                    zip(
                        ("n_helix", "n_sheet", "n_coil"),
                        sse_map.coverage(selection.points.identities),
                    )
                )
                if sse_map is not None
                else None
            ),
        },
    )


def _require_chain_present(structure: Structure, chain: str, what: str = "chain") -> None:
    """Refuse a chain the structure does not contain.

    The gap this closes: ``analyze(path, method="protein",
    chain="ZZ")`` on a structure with only chain A returned a **complete,
    well-formed result with zero points, zero triangles, no warning and no skip
    record**.  It is indistinguishable from "this chain genuinely has no TSR
    keys", it carries both fingerprints, it serializes, and an analysis over a
    manifest containing one typo would silently include it.

    Every selector had the same hole, because every one of them filters
    ``structure.atoms`` by chain and an empty filter is not an error to a filter.
    The check therefore lives here, once, above all of them, rather than being
    added to each.

    **It is about the chain, not about the selection**, and that distinction is
    load-bearing. ``biotsr/cli/spec.py::_require_points`` refuses an empty
    calculation at the command line and its docstring gives the reason
    ``analyze`` must *not*: "an empty selection is a fact about the structure ...
    a library must not decide for its caller that a fact is an error." That is
    exactly right for a chain that **exists** and yields nothing -- ``protein``
    over 1BNA's DNA chain A selects no Ca and is correct to. It is not a fact
    about anything when the named chain is absent: there is no structure being
    described, only a request that misnames one. So this refuses case one and
    leaves case two to the CLI guard, which still handles it unchanged.

    The wording deliberately keeps the two phrases the command-line contract
    promises for a bad chain -- "no atoms" and "chains in this structure" --
    because moving *where* a request is refused should not change *what the user
    is told*, and the command-line validation asserts exactly that.

    Silence on a typo is precisely the property the package's skip ledger exists
    to eliminate: "nothing is dropped silently" has to include dropping
    everything.
    """
    present = structure.chains()
    if any(c.value == chain for c in present):
        return
    names = ", ".join(c.value for c in present) or "(none)"
    raise InvalidRequest(
        "%s %r is not in %s, so the selection produced no atoms; chains in this "
        "structure: %s. A run against a chain that is not there would otherwise "
        "produce an empty result carrying a valid fingerprint, which is "
        "indistinguishable from a real finding."
        % (what, chain, structure.structure_id.value, names)
    )


def _require_standard_profile(recipe: Method, profile) -> None:
    """Accept ``profile=None`` or the method's own standard profile, nothing else.

    Every public method runs under the accepted standard conventions for that
    method.  A profile built for another method, a derived profile, or a
    profile given by name would produce keys whose meaning nobody has defined
    for this method, so each is refused here, by name, before any structure is
    read -- through the Python API exactly as on the command line, which has no
    way to name a profile at all.
    """
    if profile is None:
        return
    standard = recipe.default_profile
    if isinstance(profile, str):
        raise InvalidRequest(
            "profile=%r: BioTSR runs method %r under its standard profile only "
            "(%s), so there is no profile to choose by name. Omit profile=."
            % (profile, recipe.name, standard.name)
        )
    if not isinstance(profile, Profile):
        raise InvalidRequest(
            "profile must be omitted (or None), not %s" % type(profile).__name__
        )
    if profile is standard or profile == standard:
        return
    raise InvalidRequest(
        "profile %r is not the standard profile of method %r (%s). Every method "
        "runs under its own standard profile; a profile from another method, or "
        "one built by changing a standard profile, defines keys this method does "
        "not. Omit profile=, and use the modifiers the method supports "
        "(%s) to vary a run."
        % (profile.name, recipe.name, standard.name,
           ", ".join(sorted(recipe.supported_modifiers)) or "it has none")
    )


def _refuse_cross_only(recipe: Method, modifiers) -> None:
    """``mode=`` is refused by name: there is one mode, the standard."""
    if "mode" in modifiers:
        raise InvalidRequest(
            "mode=%r: BioTSR 1.0 has one scientific mode for every method, the "
            "standard, so there is no mode argument. Omit mode=."
            % (modifiers["mode"],)
        )


def _selection_warnings(selection, counts, explained=()) -> tuple:
    """What a result must say about its own selection, so nobody infers it.

    * every explicitly requested residue, atom or entity that was not found;
    * the instances or residues a completeness rule excluded, when others
      were kept (a result that keeps none says so in its "no triangles" line);
    * how many atom records the alternate-location rule excluded;
    * a selection too small to form a triangle -- a result with no keys is
      not an ordinary calculation, and is labelled as such.
    """
    out = []
    missing = tuple(m for m in (getattr(selection, "atoms_missing", ()) or ())
                    if m not in explained)
    if missing:
        out.append("requested but not found: %s" % ", ".join(missing))
    rejected = [s.detail for s in selection.skipped if s.stage == "completeness"]
    n = len(selection.points)
    if rejected and n >= 3:
        out.append(
            "excluded by the method's completeness rule (%d): %s"
            % (len(rejected), "; ".join(rejected[:10])
               + ("; and %d more" % (len(rejected) - 10) if len(rejected) > 10 else ""))
        )
    import re

    altlocs = sorted({m.group(1) for m in (
        re.search(r"altloc=([^/\s)]+)", s.detail) for s in selection.skipped
        if s.reason == "altloc_rejected") if m})
    n_altloc = sum(1 for s in selection.skipped if s.reason == "altloc_rejected")
    if n_altloc:
        out.append(
            "alternate locations: kept blank and A; excluded %d atom record(s)%s"
            % (n_altloc, " with indicator %s" % ", ".join(altlocs) if altlocs else "")
        )
    if n < 3:
        out.append(
            "no triangles: the selection has %d point%s and a triangle needs three"
            "%s" % (n, "" if n == 1 else "s",
                    "; excluded by the method's completeness rule: %s"
                    % "; ".join(rejected) if rejected else "")
        )
    return tuple(out)


def _entity_instances(structure, chain, admission, name):
    """The residue identifiers of every copy of entity ``name`` in ``chain``, in order."""
    from .selection.atomic import _admitted_chain_atoms
    from .domain.identity import ChainId

    wanted = str(name).strip().upper()
    seen = []
    for atom, _cutoff in _admitted_chain_atoms(structure, ChainId(chain), admission):
        if atom.identity.residue_name.upper() != wanted:
            continue
        ok, _reason = admission.accepts(atom)
        if not ok:
            continue
        label = atom.identity.residue.legacy_seq_label
        if label not in seen:
            seen.append(label)
    return seen


def _split_identifiers(structure, chain, admission, name):
    """Identifiers of ``name`` whose records appear in two separate places.

    A residue identifier names one molecule only if its records are one
    contiguous block; when the same name and number reappear after another
    residue, the file holds two molecules BioTSR cannot tell apart, and a
    selection naming that identifier is refused rather than pooled.
    """
    from .selection.atomic import _admitted_chain_atoms
    from .domain.identity import ChainId

    wanted = str(name).strip().upper()
    previous = None
    seen = set()
    split = []
    for atom, _cutoff in _admitted_chain_atoms(structure, ChainId(chain), admission):
        ok, _reason = admission.accepts(atom)
        if not ok:
            continue
        key = (atom.identity.residue_name.upper(), atom.identity.residue.legacy_seq_label)
        if key == previous:
            continue
        if key in seen and key[0] == wanted and key[1] not in split:
            split.append(key[1])
        seen.add(key)
        previous = key
    return split


def _entity_names_in(structure, chain, admission):
    from .selection.atomic import entity_targets

    names = []
    for _seq, name in entity_targets(structure, chain, admission):
        if name not in names:
            names.append(name)
    return names


def _resolve_entity_instances(recipe_name, structure, chain, admission, entity_name,
                              entity_id, what="entity"):
    """Which copies of ``entity_name`` on ``chain`` a calculation uses, and its scope.

    * with ``entity_id`` -- that one copy (``one_instance``).  The identifier is
      the molecule's PDB residue number with an optional insertion code, never
      a position; it must exist on the named chain;
    * without -- every copy of that name on the named chain, in one combined
      calculation (``all_copies``), even when there is only one.  Triangles may
      then span the copies, by design.

    The chain is always the caller's: no other chain is searched.  Returns
    ``(identifiers, scope)``.
    """
    from .domain.selection_record import SCOPE_ALL_COPIES, SCOPE_ONE_INSTANCE

    copies = _entity_instances(structure, chain, admission, entity_name)
    split = _split_identifiers(structure, chain, admission, entity_name)
    name = str(entity_name).strip().upper()
    if entity_id not in (None, ""):
        label = parse_residue_id(entity_id, "%s id" % what)
        if label not in copies:
            raise InvalidRequest(
                "%s %s %s is not in chain %s of %s: %s"
                % (what, name, label, chain, structure.structure_id.value,
                   ("copies of %s there have ids %s" % (name, ", ".join(copies)))
                   if copies else _absent(structure, chain, admission, entity_name))
            )
        if label in split:
            raise InvalidRequest(
                "%s %s %s appears in two separate places in chain %s of %s, so the "
                "identifier does not name one molecule; BioTSR cannot tell the "
                "two apart and does not pool them. Renumber one of them in the "
                "structure file." % (what, name, label, chain,
                                     structure.structure_id.value)
            )
        return (label,), SCOPE_ONE_INSTANCE
    if not copies:
        raise InvalidRequest(
            "%s %s is not in chain %s of %s: %s"
            % (what, name, chain, structure.structure_id.value,
               _absent(structure, chain, admission, entity_name))
        )
    if split:
        raise InvalidRequest(
            "%s %s id%s %s appear%s in two separate places in chain %s of %s, so "
            "the copies of %s there cannot be told apart; BioTSR does not pool "
            "them. Renumber them in the structure file."
            % (what, name, "" if len(split) == 1 else "s", ", ".join(split),
               "s" if len(split) == 1 else "", chain, structure.structure_id.value,
               name)
        )
    return tuple(copies), SCOPE_ALL_COPIES


def _instance_tokens(residues_used):
    """``("1_NAG", "2_NAG")`` (``<seq>_<NAME>``) -> ``["NAG_1", "NAG_2"]``."""
    out = []
    for entry in residues_used:
        seq, _, name = str(entry).partition("_")
        out.append("%s_%s" % (name, seq))
    return out


def _absent(structure, chain, admission, entity_name) -> str:
    names = _entity_names_in(structure, chain, admission)
    return (
        "no %s there; the entities this method can see in that chain are: %s"
        % (entity_name, ", ".join(names) if names else "(none)")
    )


#: Why a modified residue stored as ``HETATM`` is not part of a calculation.
UNSUPPORTED_MODIFIED = (
    "a modified residue stored as HETATM; of the modified amino acids BioTSR "
    "reads only %s, which have defined atom counts" % ", ".join(MODIFIED_AMINO_ACIDS))
NOT_REPRESENTED_CALPHA = REASON_NOT_LABELLED


class _ContactOnlyLabels(LabelScheme):
    """Every atom resolves to one placeholder code: for atoms read only for
    their coordinates and residue identity (the Calpha cross interface)."""

    def resolve_atom(self, atom_name, element=""):
        return 1


#: The interface heavy atoms of ``protein_protein_calpha_cross`` (ADR-052).
CONTACT_ONLY_LABELS = _ContactOnlyLabels(
    name="contact_only", table={"*": 1}, n_labels=1,
    domain=LabelDomain.ATOM_ELEMENT,
    description="Placeholder labels for atoms used only for contact detection.")


def _calpha_unread(structure, chain, selection, requested=None):
    """``(tokens, notes)`` for the residues a Calpha selection leaves out.

    ``tokens`` are ``NAME_ID`` in file order (requested order with
    ``requested``); ``notes`` are warnings, one per reason, each naming its
    residues -- never "not found" for a residue that is there (ADR-107).
    """
    excluded = calpha_exclusions(structure, chain, selection.residues_used,
                                 requested=requested)
    tokens = ["%s_%s" % (name, seq) for seq, name, _ in excluded]
    by_reason = {}
    for seq, name, reason in excluded:
        by_reason.setdefault(reason, []).append("%s_%s" % (name, seq))
    notes = ["not represented: %s (%s)" % (", ".join(t), r) for r, t in by_reason.items()]
    return tokens, notes, [seq for seq, _, _ in excluded]


def processing_identity(method: str):
    """The processing rules a calculation of ``method`` is drawn under.

    ``{"altloc_policy": ..., "modified_residue_policy": ..., "input_rules":
    ...}``, plus ``"completeness_policy"`` for a method that has one -- what
    every completion record stores and ``--if-exists skip`` compares (ADR-105,
    ADR-114).  The modified-residue rule depends on whether the method labels
    residues (the Calpha methods) or atoms (ADR-107).
    """
    recipe = get_method(method)
    return _processing_identity(recipe.default_profile.label_scheme.domain,
                                recipe.completeness_policy)


def _unread_residues(structure, chain, admission):
    """``{seq_label: name}`` of the chain's ``HETATM`` amino-acid residues an
    atomic selection does not read (ADR-101).  The Calpha methods report
    theirs through :func:`_calpha_unread` (ADR-107)."""
    found = hetero_amino_acid_residues(structure, chain, admission.truncation)
    return {seq: name for seq, name in found
            if name not in admission.hetero_residue_names}


def _unsupported_warning(tokens, calpha=False) -> str:
    return "%s: %s (%s)" % (
        "not represented" if calpha else "not supported", ", ".join(tokens),
        NOT_REPRESENTED_CALPHA if calpha else UNSUPPORTED_MODIFIED)


def _refuse_named_unk(structure, chain, tokens) -> None:
    """Refuse a request that names UNK residues (``tokens``) explicitly."""
    integrity.refuse_unk(
        [t for t in tokens if t == integrity.UNK or t.startswith(integrity.UNK + "_")],
        "chain %s of %s" % (chain, structure.structure_id.value))


def _declined_tokens(selection, refs):
    """``NAME_ID`` of the named residues the completeness rule declined."""
    out = []
    for ref in refs:
        stem = "%s_%s:" % (ref.residue_seq, ref.residue_name.upper())
        if any(s.stage == "completeness" and s.detail.startswith(stem)
               for s in selection.skipped):
            out.append(ref.token())
    return out


def _select(recipe, profile, structure, chain, *, residue, residue_name, atoms,
            entity_name=None, entity_id=None, entity_names=None,
            residue_types=None, residues=None):
    """Route to the method's own selector. One place, no engine-side branching.

    Returns ``(selection, target_label, warnings, record, explained)``;
    ``record`` is the structured selection record
    (:mod:`biotsr.domain.selection_record`), built from the values actually
    used and the instances actually selected, and ``explained`` the requested
    items this function has already reported (not found, not supported,
    declined), which the generic "requested but not found" line leaves out.
    """
    from .domain import selection_record as sr

    warnings = ()

    def record(scope, **arguments):
        return sr.selection_record([sr.side_record(
            side="primary", selector=recipe.selector, chain=chain, scope=scope,
            **arguments)])

    if recipe.selector == SELECTOR_RESIDUE_ATOMS:
        if residue is None or residue_name is None:
            raise ValueError(
                "method %r needs a target residue: pass residue=<residue id> and "
                "residue_name=<three-letter name>, or use analyze_entities() to "
                "run every residue of the chain (optionally only residues of one "
                "name)" % recipe.name
            )
        label = parse_residue_id(residue, "residue")
        name = str(residue_name).strip().upper()
        _refuse_named_unk(structure, chain, ["%s_%s" % (name, label)])
        refuse_unreadable(structure, chain, profile.atom_admission, names=(name,))
        selection = select_residue_atoms(
            structure, chain, label, name,
            profile.label_scheme,
            admission=profile.atom_admission,
            completeness=profile.completeness,
        )
        if not selection.points:
            # Empty because the residue is not there (refused), or because the
            # completeness rule declined it (a result whose warning says so).
            present = dict(residue_targets(structure, chain, profile.atom_admission))
            if present.get(label) != name:
                unread = _unread_residues(structure, chain, profile.atom_admission)
                if unread.get(label) == name:
                    raise InvalidRequest(
                        "residue %s %s of chain %s of %s is not supported: it is %s"
                        % (label, name, chain, structure.structure_id.value,
                           UNSUPPORTED_MODIFIED))
                raise InvalidRequest(
                    "residue %s %s is not in chain %s of %s%s"
                    % (label, name, chain, structure.structure_id.value,
                       "; residue %s there is %s" % (label, present[label])
                       if label in present else "")
                )
        return (selection, "%s_%s" % (label, name), warnings,
                record(sr.SCOPE_ONE_RESIDUE, residue=label, residue_name=name,
                       instances=_instance_tokens(selection.residues_used),
                       excluded=[] if selection.points else ["%s_%s" % (name, label)],
                       declined=[] if selection.points else ["%s_%s" % (name, label)]),
                ())

    if recipe.selector == SELECTOR_ATOM_SELECTION:
        if atoms is None:
            raise ValueError(
                "method %r needs an atom selection: pass atoms='LYS_48_NZ;...' "
                "or a list of such tokens" % recipe.name
            )
        requests, warnings = parse_atom_selection(atoms)
        if not requests:
            raise ValueError(
                "atom selection %r produced no usable RESIDUE_ID_ATOM entries"
                % (atoms,)
            )
        _refuse_named_unk(structure, chain, [
            "%s_%s" % (r.residue_name, r.residue_seq) for r in requests])
        selection = select_atom_selection(
            structure, chain, requests, profile.label_scheme,
            admission=profile.atom_admission,
        )
        unread = _unread_residues(structure, chain, profile.atom_admission)
        unsupported = [r.token() for r in requests
                       if unread.get(r.residue_seq) == r.residue_name.upper()]
        not_found = [t for t in selection.atoms_missing if t not in unsupported]
        if not selection.points:
            raise InvalidRequest(
                "none of the requested atoms is in chain %s of %s: %s"
                % (chain, structure.structure_id.value, "; ".join(filter(None, (
                    ", ".join(not_found),
                    _unsupported_warning(unsupported) if unsupported else ""))))
            )
        if unsupported:
            warnings = tuple(warnings) + (_unsupported_warning(unsupported),)
        tokens = [r.token() for r in requests]
        return (selection, ";".join(tokens), warnings,
                record(sr.SCOPE_NAMED_ATOMS, atom_selection=tokens,
                       not_found=not_found, unsupported=unsupported),
                tuple(unsupported))

    if recipe.selector == SELECTOR_RESIDUE_TYPE_ATOMS:
        types, warnings = parse_residue_types(residue_types)
        _check_arity(recipe, types)
        _refuse_named_unk(structure, chain, types)
        refuse_unreadable(structure, chain, profile.atom_admission, names=types)
        absent = [t for t in types
                  if not _entity_instances(structure, chain, profile.atom_admission, t)]
        unread = _unread_residues(structure, chain, profile.atom_admission)
        unsupported = [t for t in absent if t in unread.values()]
        not_found = [t for t in absent if t not in unsupported]
        found_text = ("requested but not found: residue type%s %s"
                      % ("" if len(not_found) == 1 else "s", ", ".join(not_found)))
        unsupported_text = _unsupported_warning(
            ["residue type %s" % ", ".join(unsupported)]) if unsupported else ""
        if len(absent) == len(types):
            if not unsupported:
                raise InvalidRequest(
                    "no residue of type %s is in chain %s of %s"
                    % (", ".join(types), chain, structure.structure_id.value)
                )
            raise InvalidRequest(
                "nothing to compute in chain %s of %s: %s"
                % (chain, structure.structure_id.value, "; ".join(
                    filter(None, (found_text if not_found else "", unsupported_text)))))
        selection = select_entity_atoms(
            structure, chain, profile.label_scheme,
            admission=profile.atom_admission,
            entity_names=types,
            completeness=profile.completeness,
            selector=SELECTOR_RESIDUE_TYPE_ATOMS,
        )
        if not_found:
            warnings = tuple(warnings) + (found_text,)
        if unsupported:
            warnings = tuple(warnings) + (unsupported_text,)
        return (selection, ";".join(types), warnings,
                record(sr.SCOPE_RESIDUE_TYPES, residue_types=types,
                       not_found=not_found, unsupported=unsupported),
                ())

    if recipe.selector == SELECTOR_SPECIFIC_RESIDUE_ATOMS:
        refs, warnings = parse_residue_refs(residues)
        if not refs:
            raise ValueError(
                "method %r needs residues: pass residues='LYS_6;GLU18'"
                % recipe.name
            )
        _refuse_named_unk(structure, chain, [r.token() for r in refs])
        refuse_unreadable(structure, chain, profile.atom_admission,
                          names=[r.residue_name for r in refs])
        selection = select_entity_atoms(
            structure, chain, profile.label_scheme,
            admission=profile.atom_admission,
            entity_refs=refs,
            completeness=profile.completeness,
            selector=SELECTOR_SPECIFIC_RESIDUE_ATOMS,
        )
        tokens = [r.token() for r in refs]
        used = _instance_tokens(selection.residues_used)
        unread = _unread_residues(structure, chain, profile.atom_admission)
        unsupported = [r.token() for r in refs
                       if unread.get(r.residue_seq) == r.residue_name.upper()]
        declined = _declined_tokens(selection, refs)
        not_found = [t for t in tokens if t not in used and t not in unsupported
                     and t not in declined]
        notes = []
        if not_found:
            notes.append("requested but not found: %s" % ", ".join(not_found))
        if unsupported:
            notes.append(_unsupported_warning(unsupported))
        if declined:
            notes.append("requested but excluded by the method's completeness "
                         "rule: %s" % ", ".join(declined))
        if not selection.points and not used:
            if not unsupported and not declined:
                raise InvalidRequest(
                    "none of the requested residues is in chain %s of %s: %s"
                    % (chain, structure.structure_id.value, ", ".join(tokens))
                )
            raise InvalidRequest(
                "none of the requested residues can be used in chain %s of %s: %s"
                % (chain, structure.structure_id.value, "; ".join(notes)))
        warnings = tuple(warnings) + tuple(notes)
        return (selection, ";".join(tokens), warnings,
                record(sr.SCOPE_NAMED_RESIDUES, residues=tokens, instances=used,
                       excluded=[t for t in tokens if t not in used],
                       not_found=not_found, unsupported=unsupported,
                       declined=declined),
                tuple(t for t in tokens if t not in used))

    if recipe.selector in (SELECTOR_ENTITY_ATOMS, SELECTOR_ENTITY_ATOM_SELECTION):
        if not entity_name:
            raise ValueError(
                "method %r needs the molecule it applies to: pass "
                "entity_name=<name>, and entity_id=<residue id> for one copy "
                "(omit it for every copy on the chain), e.g. entity_name='ESW', "
                "entity_id='404'" % recipe.name
            )
        want_atoms = recipe.selector == SELECTOR_ENTITY_ATOM_SELECTION
        atom_names = parse_atom_names(atoms) if atoms is not None else ()
        if want_atoms and not atom_names:
            raise ValueError(
                "method %r needs an atom list: pass atoms='N4;C5;O6' -- use %r for "
                "the whole molecule instead"
                % (recipe.name, recipe.name.replace("_atom_selection", ""))
            )
        name = str(entity_name).strip().upper()
        _refuse_named_unk(structure, chain, [name])
        refuse_unreadable(structure, chain, profile.atom_admission, names=(name,))
        instances, scope = _resolve_entity_instances(
            recipe.name, structure, chain, profile.atom_admission,
            entity_name, entity_id,
        )
        selection = select_entity_atoms(
            structure, chain, profile.label_scheme,
            admission=profile.atom_admission,
            entity_names=(entity_name,),
            entity_seqs=instances,
            atom_names=atom_names or None,
            completeness=profile.completeness,
            selector=recipe.selector,
        )
        if (want_atoms and not selection.points
                and not any(s.stage == "completeness" for s in selection.skipped)):
            raise InvalidRequest(
                "none of the requested atoms (%s) is in %s %s of chain %s"
                % (", ".join(atom_names), name,
                   instances[0] if scope == sr.SCOPE_ONE_INSTANCE
                   else "(any copy: %s)" % ", ".join(instances), chain)
            )
        used = _instance_tokens(selection.residues_used)
        declined = ["%s_%s" % (name, i) for i in instances
                    if "%s_%s" % (name, i) not in used
                    and any(s.stage == "completeness" and s.detail.startswith(
                        "%s_%s:" % (i, name)) for s in selection.skipped)]
        target = ("%s_%s" % (name, instances[0]) if scope == sr.SCOPE_ONE_INSTANCE
                  else "%s_allcopies" % name)
        return (selection, target, warnings,
                record(scope, entity_name=name,
                       entity_id=instances[0] if scope == sr.SCOPE_ONE_INSTANCE else None,
                       atoms=atom_names, instances=used, excluded=declined,
                       # Declined counts as an omission only for the one copy
                       # the caller named; a copy declined among all copies is
                       # recorded in ``excluded``, as a residue declined among
                       # all residues of a type is in the skip records.
                       declined=declined if scope == sr.SCOPE_ONE_INSTANCE else (),
                       atom_subset_missing=dict(selection.atoms_missing_by_instance)),
                ())

    if recipe.selector == SELECTOR_ENTITY_CHAIN_ATOMS:
        # ``_split``, never ``tuple()``: the command line and a manifest hand this
        # over as the string "NAG;FUC", and tuple() of a string is its characters.
        names = _split(entity_names) if entity_names else recipe.default_entity_names
        _refuse_named_unk(structure, chain, [str(n).strip().upper() for n in names])
        refuse_unreadable(structure, chain, profile.atom_admission, names=names or None)
        if entity_names:
            absent = [n for n in names
                      if not _entity_instances(structure, chain,
                                               profile.atom_admission, n)]
            if len(absent) == len(names):
                raise InvalidRequest(
                    "none of the requested entities (%s) is in chain %s of %s: %s"
                    % (", ".join(names), chain, structure.structure_id.value,
                       _absent(structure, chain, profile.atom_admission,
                               "/".join(names)))
                )
            if absent:
                warnings = tuple(warnings) + (
                    "requested but not found: %s" % ", ".join(absent),)
        selection = select_entity_atoms(
            structure, chain, profile.label_scheme,
            admission=profile.atom_admission,
            entity_names=names or None,
            completeness=profile.completeness,
            selector=SELECTOR_ENTITY_CHAIN_ATOMS,
        )
        # The vocabulary is recorded only when the caller gave one; otherwise
        # it is the method's declared default, named by the method itself.
        return (selection, chain, warnings,
                record(sr.SCOPE_CHAIN_VOCABULARY,
                       entity_names=names if entity_names else None,
                       instances=_instance_tokens(selection.residues_used),
                       not_found=absent if entity_names else ()),
                ())

    if recipe.selector == SELECTOR_CALPHA_SELECTION:
        labels, warnings = parse_residue_selection(residues)
        if not labels:
            raise ValueError(
                "method %r needs a residue selection: pass "
                "residues='10-12;48;70-71'" % recipe.name
            )
        check_range_insertion_codes(structure, chain, residues, labels)
        selection = select_calpha_residues(
            structure, chain, labels, profile.label_scheme,
            altloc_policy=AltLocPolicy.BLANK_OR_A,
            truncation=LEGACY_SINGLE_PASS,
        )
        missing = list(selection.atoms_missing)
        unsupported, unread_notes, present = _calpha_unread(
            structure, chain, selection, requested=missing)
        not_found = [m for m in missing if m not in present]
        notes = []
        if not_found:
            notes.append("requested but not found: %s" % ", ".join(not_found))
        notes.extend(unread_notes)
        if not selection.points:
            raise InvalidRequest(
                "none of the requested residues has an admitted Calpha in chain "
                "%s of %s: %s" % (chain, structure.structure_id.value,
                                  "; ".join(notes) or ", ".join(labels))
            )
        warnings = tuple(warnings) + tuple(notes)
        # The selection as given (ranges kept), and the residues it expanded to.
        given = [t.strip() for t in str(residues).replace(",", ";").split(";")
                 if t.strip()] if isinstance(residues, str) else list(residues)
        return (selection, ";".join(labels), warnings,
                record(sr.SCOPE_RESIDUE_LIST, residues=given,
                       instances=list(selection.residues_used), excluded=missing,
                       not_found=not_found, unsupported=unsupported),
                tuple(missing))

    selection = select_calpha(
        structure, chain, profile.label_scheme,
        altloc_policy=AltLocPolicy.BLANK_OR_A,
        truncation=LEGACY_SINGLE_PASS,
        seq_label_policy=profile.seq_label_policy,
    )
    # The chain's modified residues are part of it but not of its Calpha
    # points: said, so a chain with a gap is never read as the whole chain.
    unsupported, unread_notes, _ = _calpha_unread(structure, chain, selection)
    if unsupported:
        warnings = tuple(warnings) + tuple(unread_notes)
    return (selection, chain, warnings,
            record(sr.SCOPE_WHOLE_CHAIN, unsupported=unsupported), ())


def _check_arity(recipe, types):
    """Enforce ``aa_ind`` = one type and ``aa_ind_mix`` = two or more.

    The two recipes are otherwise identical, so without this a one-type
    ``aa_ind_mix`` would quietly be an ``aa_ind`` under a different name in
    provenance -- and provenance is the whole point of having two names.
    """
    if not types:
        raise ValueError(
            "method %r needs residue types: pass residue_types='TYR' or "
            "residue_types='TYR;PHE'" % recipe.name
        )
    if len(types) < recipe.min_entity_types:
        raise ValueError(
            "method %r needs at least %d residue types; got %d (%s). Use %r for a "
            "single type."
            % (recipe.name, recipe.min_entity_types, len(types), ", ".join(types),
               "aa_ind")
        )
    if recipe.max_entity_types is not None and len(types) > recipe.max_entity_types:
        raise ValueError(
            "method %r takes at most %d residue type(s); got %d (%s). Use %r for "
            "two or more."
            % (recipe.name, recipe.max_entity_types, len(types), ", ".join(types),
               "aa_ind_mix")
        )


def list_residue_targets(source, *, chain: str, method: str = "amino_acid_atomic",
                         structure_id: Optional[str] = None):
    """Every ``(residue_seq_label, residue_name)`` an atomic method would consider.

    Offered, not yet admitted -- the completeness policy has its say later, and
    records a reason when it declines.
    """
    recipe = get_method(method)
    structure = load_structure(source, structure_id=structure_id)
    _require_chain_present(structure, chain)
    return residue_targets(structure, chain, recipe.default_profile.atom_admission)


def analyze_entities(
    source,
    *,
    method: str = "amino_acid_atomic",
    chain: str,
    profile: Optional[Profile] = None,
    structure_id: Optional[str] = None,
    keep_triplets: bool = True,
    include_rejected: bool = False,
    residue_name: Optional[str] = None,
    **modifiers,
):
    """Run a per-residue method over every residue of a chain -- one result each.

    ``residue_name`` narrows the sweep to residues of that name (``"LEU"``:
    every leucine of the chain, each its own calculation).  Each result is
    exactly the one-residue calculation ``analyze(..., residue=<id>,
    residue_name=<name>)``; nothing is pooled.

    Returns ``{"<seq>_<NAME>": TSRResult}`` in structure order.  Residues the
    completeness policy declines are **omitted by default** -- which is the legacy
    outcome, where a rejected residue produces no files at all -- but pass
    ``include_rejected=True`` to get their (empty, zero-triangle) results and read
    the recorded reason instead of inferring it from an absence.
    """
    recipe = get_method(method)
    if recipe.scope != SCOPE_PER_ENTITY:
        raise ValueError(
            "analyze_entities() is for per-entity methods; %r has scope %r"
            % (recipe.name, recipe.scope)
        )
    structure = load_structure(source, structure_id=structure_id)
    # Without this, a mistyped chain returns `{}` -- "this chain has no residues
    # this method can work on", which is a real answer for a real chain and a
    # silent lie for one that is not there.
    _require_chain_present(structure, chain)
    _require_standard_profile(recipe, profile)
    if "residue" in modifiers:
        raise InvalidRequest(
            "analyze_entities() runs every residue of the chain (or every residue "
            "of one name, with residue_name=), so residue= would contradict it; "
            "use analyze(method=%r, residue=..., residue_name=...) for one "
            "residue" % recipe.name
        )
    # A residue whose number cannot be read would be missing from the sweep.
    refuse_unreadable(structure, chain, recipe.default_profile.atom_admission,
                      names=(residue_name,) if residue_name else None)
    targets = residue_targets(structure, chain, recipe.default_profile.atom_admission)
    # A sweep over the whole chain includes its UNK residues, and one over a
    # residue name includes them when the name is UNK: either is refused, not
    # computed without them (ADR-109).
    _refuse_named_unk(structure, chain, [
        "%s_%s" % (name, seq) for seq, name in targets
        if not residue_name or name.upper() == str(residue_name).strip().upper()])
    if residue_name:
        wanted = str(residue_name).strip().upper()
        targets = tuple((seq, name) for seq, name in targets if name.upper() == wanted)
        if not targets:
            names = sorted({name for _seq, name in residue_targets(
                structure, chain, recipe.default_profile.atom_admission)})
            unread = _unread_residues(structure, chain,
                                      recipe.default_profile.atom_admission)
            if wanted in unread.values():
                raise InvalidRequest(
                    "residue %s of chain %s of %s is not supported: %s"
                    % (", ".join("%s_%s" % (wanted, seq) for seq, name
                                 in unread.items() if name == wanted),
                       chain, structure.structure_id.value, UNSUPPORTED_MODIFIED))
            raise InvalidRequest(
                "no residue named %s is in chain %s of %s; residue names there: %s"
                % (wanted, chain, structure.structure_id.value,
                   ", ".join(names) or "(none)")
            )

    results = {}
    for seq, name in targets:
        result = analyze(
            structure,
            method=method,
            chain=chain,
            profile=profile,
            structure_id=structure_id,
            keep_triplets=keep_triplets,
            residue=seq,
            residue_name=name,
            **modifiers,
        )
        if result.counts.n_points == 0 and not include_rejected:
            continue
        results["%s_%s" % (seq, name)] = result
    return results


def _sampling_profile(recipe: Method, profile: Profile):
    from .domain.keyspace import SamplingProfile

    return SamplingProfile(
        selector=recipe.selector,
        scope=recipe.scope,
        candidate_strategy=recipe.candidate_strategy,
        completeness_policy=recipe.completeness_policy,
        composition_rule=recipe.composition_rule,
        size_filter=profile.size_filter.cutoff if profile.size_filter else None,
        degenerate_policy=profile.degenerate_policy_label(),
        seq_label_policy=profile.seq_label_policy,
    )


def list_entity_targets(source, *, chain: str, method: str,
                        entity_names=None, structure_id: Optional[str] = None):
    """Every ``(seq_label, entity_name)`` a family's admission rules can see.

    The generalization of :func:`list_residue_targets` to the drug, carbohydrate
    and nucleotide families: it answers "what is actually in this chain that this
    method could work on?" without committing to a calculation.  Entities the
    completeness policy would later decline are still listed -- admission is a
    separate, recorded decision.
    """
    recipe = get_method(method)
    structure = load_structure(source, structure_id=structure_id)
    _require_chain_present(structure, chain)
    names = (
        _split(entity_names) if entity_names is not None
        else (recipe.default_entity_names or None)
    )
    return entity_targets(
        structure, chain, recipe.default_profile.atom_admission,
        entity_names=names,
    )


# ===========================================================================
# Cross TSR
# ===========================================================================

def default_contact_cutoff(method: str) -> float:
    """The reference contact cutoff of a cross method, in angstroms."""
    from .methods.builtin.cross import cross_recipe

    return float(cross_recipe(method).default_cutoff)


def unsupported_residues(source, *, chain: str, method: str) -> Tuple[str, ...]:
    """``NAME_ID`` of the chain's modified residues ``method`` cannot represent.

    The ``HETATM`` amino-acid residues its admission does not read -- those
    other than the supported modified residues of the residue labels for a
    Calpha method (ADR-107), those other than TPO, SEP and PTR for an atomic
    one (ADR-101).  A per-residue sweep (:func:`analyze_entities`) covers the
    residues the method reads; this names the ones it leaves out, so a caller
    can say so.
    """
    recipe = get_method(method)
    structure = load_structure(source)
    if recipe.selector in (SELECTOR_NAME, SELECTOR_CALPHA_SELECTION):
        selection = select_calpha(structure, chain, recipe.default_profile.label_scheme,
                                  unlabelled_policy=UnlabelledPointPolicy.DROP_AND_RECORD,
                                  refuse_unk=False)
        return tuple(_calpha_unread(structure, chain, selection)[0])
    return tuple("%s_%s" % (name, seq) for seq, name in _unread_residues(
        structure, chain, recipe.default_profile.atom_admission).items())


def cross_entity_role(method: str) -> Optional[str]:
    """The molecule side's role of a cross method ("drug", "carbohydrate",
    "nucleotide"), or ``None`` for the protein-protein methods."""
    from .methods.builtin.cross import cross_recipe

    recipe = cross_recipe(method)
    if recipe.entity_side is None:
        return None
    return recipe.side_a_role if recipe.entity_side == "a" else recipe.side_b_role


def list_cross_methods():
    """The registered cross methods, in name order."""
    from .methods.builtin.cross import CROSS_METHODS

    return tuple(sorted(m.name for m in CROSS_METHODS))


def _split(value):
    """Accept ``"A;B"``, ``["A", "B"]`` or ``None`` uniformly."""
    if value is None:
        return ()
    if isinstance(value, str):
        return tuple(v.strip() for v in value.split(";") if v.strip())
    return tuple(str(v).strip() for v in value if str(v).strip())


def _protein_side_selector(default, residue_types, residues, atom_selection,
                           calpha_residues):
    """Which selector a protein side uses, given the narrowing the caller asked for.

    The mapping follows the historical dataset-sheet format's column definitions
    (``biotsr.io.input_columns``), column by column::

        aa_ind / aa_ind_mix   -> residue types, all heavy atoms   (atomic, radix 112)
        aa_ind_specific       -> named residue instances, all heavy atoms   (atomic)
        aa_atom_selection     -> named atoms of named residues    (atomic)
        aa_selection          -> "using only Ca atoms"            (residues, radix 30)
        (nothing)             -> the family's default whole-chain selector
    """
    from .cross.sides import (
        SIDE_SELECTOR_ATOM_SELECTION,
        SIDE_SELECTOR_CALPHA_SELECTION,
        SIDE_SELECTOR_RESIDUE_TYPES,
        SIDE_SELECTOR_SPECIFIC_RESIDUES,
    )

    named = [
        (residue_types, SIDE_SELECTOR_RESIDUE_TYPES, "residue_types"),
        (residues, SIDE_SELECTOR_SPECIFIC_RESIDUES, "residues"),
        (atom_selection, SIDE_SELECTOR_ATOM_SELECTION, "atom_selection"),
        (calpha_residues, SIDE_SELECTOR_CALPHA_SELECTION, "calpha_residues"),
    ]
    given = [(sel, arg) for value, sel, arg in named if value]
    if not given:
        return default
    if len(given) > 1:
        raise ValueError(
            "a protein cross side takes at most one narrowing argument; got %s"
            % ", ".join(repr(arg) for _, arg in given)
        )
    return given[0][0]


def _cross_side_accounting(structure, spec, picked, admission):
    """``(not_found, unsupported, notes)`` of one cross side (ADR-101).

    What the side's request named that the structure does not contain, and
    the modified residues the side cannot represent -- for a whole-chain
    side, every one of the chain's; for a narrowed side, those it named.
    """
    from .cross import sides as cs

    calpha = spec.is_calpha
    not_found, unsupported = [], []
    if calpha:
        requested = (None if spec.selector == cs.SIDE_SELECTOR_CALPHA_CHAIN
                     else list(picked.atoms_missing))
        unsupported, unread_notes, present = _calpha_unread(
            structure, spec.chain, picked, requested=requested)
        if requested is not None:
            not_found = [m for m in requested if m not in present]
        notes = []
        if not_found:
            notes.append("requested but not found: %s (%s side, chain %s)"
                         % (", ".join(not_found), spec.role, spec.chain))
        notes.extend("%s (%s side, chain %s)" % (n, spec.role, spec.chain)
                     for n in unread_notes)
        return not_found, unsupported, notes
    unread = _unread_residues(structure, spec.chain, admission)
    outside = []
    if spec.selector == cs.SIDE_SELECTOR_PROTEIN_CHAIN:
        unsupported = ["%s_%s" % (name, seq) for seq, name in unread.items()]
        outside = _outside_vocabulary(structure, spec, admission)
        unsupported += outside
    elif spec.selector == cs.SIDE_SELECTOR_RESIDUE_TYPES:
        absent = [t for t in spec.residue_types
                  if not _entity_instances(structure, spec.chain, admission, t)]
        unsupported = [t for t in absent if t in unread.values()]
        not_found = [t for t in absent if t not in unsupported]
    elif spec.selector == cs.SIDE_SELECTOR_SPECIFIC_RESIDUES:
        refs, _ = parse_residue_refs(";".join(spec.residues))
        unsupported = [r.token() for r in refs
                       if unread.get(r.residue_seq) == r.residue_name.upper()]
        not_found = [t for t in picked.atoms_missing if t not in unsupported]
    elif spec.selector == cs.SIDE_SELECTOR_ATOM_SELECTION:
        requests, _ = parse_atom_selection(spec.atom_selection)
        unsupported = [r.token() for r in requests
                       if unread.get(r.residue_seq) == r.residue_name.upper()]
        not_found = [t for t in picked.atoms_missing if t not in unsupported]
    else:
        not_found = list(picked.atoms_missing)
    notes = []
    if not_found:
        notes.append("requested but not found: %s (%s side, chain %s)"
                     % (", ".join(not_found), spec.role, spec.chain))
    modified = [t for t in unsupported if t not in outside]
    if modified:
        notes.append("%s (%s side, chain %s)" % (
            _unsupported_warning(modified, calpha=calpha), spec.role, spec.chain))
    if outside:
        notes.append(
            "not supported: %s (residues outside this side's residue vocabulary, "
            "the 20 regular amino acids and %s; %s side, chain %s)"
            % (", ".join(outside), ", ".join(MODIFIED_AMINO_ACIDS), spec.role,
               spec.chain))
    return not_found, unsupported, notes


def _outside_vocabulary(structure, spec, admission):
    """``NAME_ID`` of the residues a vocabulary-limited whole protein side
    leaves out (the carbohydrate- and nucleotide-protein cross read the 23
    residue names of their historical contact programs; ADR-109)."""
    from .selection.truncation import truncation_line
    from .domain.identity import ChainId

    if not spec.entity_names:
        return []
    vocabulary = {n.upper() for n in spec.entity_names}
    chain_id = ChainId(spec.chain)
    cutoff = truncation_line(structure, chain_id, admission.truncation)
    out = []
    for atom in structure.atoms:
        if atom.chain != chain_id or (cutoff is not None and atom.file_order >= cutoff):
            continue
        if atom.hetero and atom.residue_name not in admission.hetero_residue_names:
            continue
        if atom.residue_name.upper() in vocabulary:
            continue
        token = "%s_%s" % (atom.residue_name, atom.identity.residue.legacy_seq_label)
        if token not in out:
            out.append(token)
    return out


def _cross_selection_record(spec_a, spec_b, selection_a, selection_b,
                            accounting=None):
    """Both sides' selection records (``biotsr.domain.selection_record``)."""
    from .cross import sides as cs
    from .domain import selection_record as sr

    def scope_of(spec):
        if spec.selector == cs.SIDE_SELECTOR_ENTITY:
            return sr.SCOPE_ONE_INSTANCE if spec.entity_id else sr.SCOPE_ALL_COPIES
        return {
            cs.SIDE_SELECTOR_ENTITY_CHAIN: sr.SCOPE_CHAIN_VOCABULARY,
            cs.SIDE_SELECTOR_RESIDUE_TYPES: sr.SCOPE_RESIDUE_TYPES,
            cs.SIDE_SELECTOR_SPECIFIC_RESIDUES: sr.SCOPE_NAMED_RESIDUES,
            cs.SIDE_SELECTOR_ATOM_SELECTION: sr.SCOPE_NAMED_ATOMS,
            cs.SIDE_SELECTOR_CALPHA_SELECTION: sr.SCOPE_RESIDUE_LIST,
        }.get(spec.selector, sr.SCOPE_WHOLE_CHAIN)

    sides = []
    accounting = accounting or {}
    for side, spec, picked in (("a", spec_a, selection_a), ("b", spec_b, selection_b)):
        not_found, unsupported, _notes = accounting.get(side, ((), (), ()))
        scope = scope_of(spec)
        names_listed = scope in (sr.SCOPE_ONE_INSTANCE, sr.SCOPE_ALL_COPIES,
                                 sr.SCOPE_NAMED_RESIDUES, sr.SCOPE_CHAIN_VOCABULARY)
        arguments = dict(
            entity_name=spec.entity_name, entity_id=spec.entity_id,
            residue_types=spec.residue_types,
            atoms=spec.atoms,
        )
        if scope == sr.SCOPE_CHAIN_VOCABULARY:
            arguments["entity_names"] = spec.entity_names
        if scope in (sr.SCOPE_NAMED_RESIDUES, sr.SCOPE_RESIDUE_LIST):
            arguments["residues"] = spec.residues
        if spec.atom_selection:
            arguments["atom_selection"] = [
                t.strip() for t in spec.atom_selection.split(";") if t.strip()]
        sides.append(sr.side_record(
            side=side, role=spec.role, selector=spec.selector, chain=spec.chain,
            scope=scope,
            instances=_instance_tokens(picked.residues_used) if names_listed else (),
            atom_subset_missing=dict(picked.atoms_missing_by_instance),
            not_found=not_found, unsupported=unsupported,
            **arguments
        ))
    return sr.selection_record(sides)


def _calpha_side_labels(spec):
    """The residue identifiers a narrowed Calpha cross side selects, or ``None``
    for a whole chain."""
    from .cross.sides import SIDE_SELECTOR_CALPHA_SELECTION

    if spec.selector != SIDE_SELECTOR_CALPHA_SELECTION:
        return None
    labels, _ = parse_residue_selection(";".join(spec.residues))
    return labels


def _refuse_cross_side_unk(structure, spec, picked, admission) -> None:
    """Refuse a cross side whose requested selection includes UNK (ADR-109).

    A whole protein side includes every residue of its chain that the side
    reads (up to the chain's ``TER``), whether or not the side's residue
    vocabulary would keep it; a narrowed side includes what it names.  The
    Calpha sides are checked by their selectors.
    """
    from .cross import sides as cs
    from .selection.truncation import truncation_line
    from .domain.identity import ChainId

    context = "the %s side, chain %s of %s" % (spec.role, spec.chain,
                                                structure.structure_id.value)
    if spec.selector == cs.SIDE_SELECTOR_PROTEIN_CHAIN:
        chain_id = ChainId(spec.chain)
        cutoff = truncation_line(structure, chain_id, admission.truncation)
        integrity.refuse_unk(integrity.unk_residues(
            a for a in structure.atoms
            if a.chain == chain_id and (cutoff is None or a.file_order < cutoff)
            and (not a.hetero or a.residue_name in admission.hetero_residue_names)),
            context)
        return
    named = list(spec.residue_types) + list(spec.entity_names)
    if spec.entity_name:
        named.append(spec.entity_name)
    if spec.atom_selection:
        named += [r.residue_name for r in parse_atom_selection(spec.atom_selection)[0]]
    if spec.selector == cs.SIDE_SELECTOR_SPECIFIC_RESIDUES:
        named += [r.residue_name for r in parse_residue_refs(";".join(spec.residues))[0]]
    integrity.refuse_unk(
        ["UNK"] if any(str(n).strip().upper() == integrity.UNK for n in named) else [],
        context)


def _refuse_shared_atoms(structure, method, spec_a, spec_b, selection_a,
                         selection_b) -> None:
    """Refuse two cross sides that select the same physical atom (ADR-109).

    Nothing is subtracted from either side: the caller's selections are left as
    they are, and the overlap is named so they can be corrected.
    """
    in_a = {p.identity.as_tuple(): p for p in selection_a.points}
    shared = [p.identity.describe() for p in selection_b.points
              if p.identity.as_tuple() in in_a]
    if not shared:
        return
    raise InvalidRequest(
        "the two sides of %r select %d of the same atom(s) in %s: %s. A cross "
        "calculation is between two disjoint sets of atoms; BioTSR does not "
        "remove the overlap for you. Change the %s side (%s) or the %s side "
        "(%s) so that they do not share atoms."
        % (method, len(shared), structure.structure_id.value,
           "; ".join(shared[:6]) + ("; and %d more" % (len(shared) - 6)
                                    if len(shared) > 6 else ""),
           spec_a.role, spec_a.describe(), spec_b.role, spec_b.describe()))


def analyze_cross(
    source: Union[str, Structure],
    *,
    method: str,
    chain: str,
    chain_2: Optional[str] = None,
    profile: Optional[Profile] = None,
    structure_id: Optional[str] = None,
    contact_cutoff: Optional[float] = None,
    contact_filter: Optional[str] = None,
    no_contact_filter: bool = False,
    entity_name: Optional[str] = None,
    entity_id: Optional[str] = None,
    entity_names=None,
    atoms=None,
    residue_types=None,
    residues=None,
    atom_selection=None,
    calpha_residues=None,
    residue_types_2=None,
    residues_2=None,
    atom_selection_2=None,
    calpha_residues_2=None,
    keep_triplets: bool = True,
    triplet_sink: Optional[Callable[[TripletRecord], None]] = None,
    size_filter: Optional[float] = None,
    **modifiers,
) -> TSRResult:
    """Run one cross-TSR calculation between two sides of a structure.

    Arguments are routed by **role**, matching the historical dataset columns,
    not by which side a family happens to list first:

    * for ``drug_protein_cross`` / ``carbohydrate_protein_cross`` /
      ``nucleotide_protein_cross``, ``chain`` is the entity chain
      (``drug_chain`` / ``carb_chain`` / ``nt_chain``) and ``chain_2`` the protein
      chain (``aa_chain``).  ``entity_name``/``entity_id`` identify one copy;
      ``entity_name`` alone uses every copy of that name on ``chain`` (scope
      ``all_copies``); ``atoms`` narrows each copy (``drug_atom_selected`` and
      friends), and ``entity_names`` instead pools every listed entity of the
      chain.  The
      protein side is narrowed with ``residue_types`` / ``residues`` /
      ``atom_selection``;
    * for ``protein_protein_cross``, ``chain`` is ``aa_chain`` and ``chain_2`` is
      ``aa_chain_2``, the unsuffixed narrowing arguments apply to the first side
      and the ``_2`` ones to the second;
    * ``protein_protein_calpha_cross`` is narrowed only by ``calpha_residues`` /
      ``calpha_residues_2`` (residue-number selections).

    **Both chains are always named** (``chain_2`` is required) and honoured:
    the calculation never scans a chain the caller did not name.  An argument
    the method does not read is refused by name.

    The science is the standard for every family: ``math.pi``, 29 theta
    classes, 18 atomic distance classes, radix 112 (radix 30 and 35 Calpha
    distance classes for ``protein_protein_calpha_cross``), full physical atom
    identity in the contact pool, symmetric atom-level contact filtering,
    exact mixed-side candidates, deterministic output.

    ``contact_cutoff`` overrides the family's reference cutoff (5 A for
    protein-protein, 3.5 A otherwise).  It is a **sampling** parameter: it
    decides which atoms enter the pool, and changing it leaves the ``KeySpace``
    fingerprint untouched.  It is emphatically *not* an edge-length constraint
    on the resulting triangles -- an atom admitted because it lies 3 A from one
    atom of the other side may well form a triangle with a 20 A edge.  To
    constrain triangle size, pass ``size_filter``, which is applied to
    ``max_distance`` after geometry and is OFF by default.

    ``contact_filter`` narrows only one side by the cutoff (``"side_a_only"``,
    ``"side_b_only"``) instead of the standard ``"both_sides"``; atomic cross
    methods only.  Like the cutoff it is a sampling decision and does not touch
    ``KeySpace``.  ``no_contact_filter=True`` skips the contact scan and uses
    both selections whole; meaningful when both sides are already explicitly
    enumerated.
    """
    from .cross import pools as cross_pools
    from .cross.runner import run_cross
    from .cross.sides import (
        SIDE_SELECTOR_ENTITY_CHAIN,
        SIDE_SELECTOR_ENTITY,
        SIDE_SELECTOR_PROTEIN_CHAIN,
        CrossSideSpec,
        select_side,
    )
    from .selection.atomic import select_entity_atoms
    from .methods.builtin.cross import cross_recipe

    recipe = get_method(method)
    if recipe.family != "cross":
        raise InvalidRequest(
            "analyze_cross() runs cross methods; %r is a one-sided method, so call "
            "analyze(method=%r, ...)" % (method, method)
        )
    cross = cross_recipe(method)
    _require_standard_profile(recipe, profile)
    _refuse_cross_only(recipe, modifiers)
    if method == "protein_protein_cross" and (calpha_residues or calpha_residues_2):
        raise InvalidRequest(
            "method %r does not use %s: it selects C-alpha residue points, which "
            "label residues "
            "in the residue label space; protein_protein_cross is fully "
            "atomic (atoms labelled by element, radix 112), and the two label spaces are never "
            "mixed (ADR-045). Use method='protein_protein_calpha_cross' for a "
            "C-alpha interface, or residues= / atom_selection= here."
            % (method, ", ".join(n for n, v in (("calpha_residues", calpha_residues),
                                               ("calpha_residues_2", calpha_residues_2))
                                 if v))
        )
    if method == "protein_protein_calpha_cross" and any(
            (residue_types, residues, atom_selection, residue_types_2, residues_2,
             atom_selection_2)):
        named = [n for n, v in (
            ("residue_types", residue_types), ("residues", residues),
            ("atom_selection", atom_selection), ("residue_types_2", residue_types_2),
            ("residues_2", residues_2), ("atom_selection_2", atom_selection_2)) if v]
        raise InvalidRequest(
            "method %r does not use %s: it is a C-alpha method, where a point is a "
            "residue, so it takes calpha_residues= / calpha_residues_2= (a "
            "residue-number selection such as '35-37;54'). For atomic "
            "protein-protein cross use method='protein_protein_cross'."
            % (method, ", ".join(named))
        )
    refuse_unused_arguments(recipe, dict(
        modifiers, chain_2=chain_2, contact_cutoff=contact_cutoff,
        contact_filter=contact_filter, no_contact_filter=no_contact_filter,
        entity_name=entity_name, entity_id=entity_id, entity_names=entity_names,
        atoms=atoms, residue_types=residue_types, residues=residues,
        atom_selection=atom_selection, calpha_residues=calpha_residues,
        residue_types_2=residue_types_2, residues_2=residues_2,
        atom_selection_2=atom_selection_2, calpha_residues_2=calpha_residues_2,
    ))
    if entity_names and (entity_name or entity_id not in (None, "") or atoms):
        raise InvalidRequest(
            "entity_names= pools every listed entity of the chain, and "
            "entity_name= / entity_id= / atoms= name one copy; they contradict "
            "each other. Pass one or the other."
        )
    selector_a, selector_b = cross.side_a_selector, cross.side_b_selector
    admission_a, admission_b = cross.side_a_admission, cross.side_b_admission

    requested = dict(modifiers)
    if size_filter is not None:
        requested["size_filter"] = size_filter
    # No `sse` here, deliberately: the 18 categories are defined over protein
    # residues. An atom inside a residue has no secondary structure of its own,
    # and a ligand, a glycan and a nucleotide have none at all, so no cross
    # method declares `sse` and passing it arrives in **modifiers, where it is
    # refused by name (ADR-036, ADR-047).
    recipe.validate_modifiers(requested)

    base = recipe.default_profile
    active = base.with_modifiers(size_filter=size_filter, method=recipe.name)

    structure = load_structure(source, structure_id=structure_id)
    # Both sides, for the same reason `analyze` checks one: a cross against a
    # chain that is not there produces an empty pool, and an empty pool produces
    # a result with zero triangles and a perfectly good fingerprint.
    _require_chain_present(structure, chain)
    if chain_2 is None:
        # Both participating chains are user-specified, for every cross method
        # (ADR-049). Defaulting the second chain to the first would be a guess --
        # and for a nucleotide or a glycan usually a wrong one, since the entity
        # rarely shares the protein's chain -- so it is asked for.
        if cross.entity_side is None:
            raise ValueError(
                "cross method %r is between two chains: pass chain_2=<the second "
                "chain>. There is no one-chain form of it -- defaulting chain_2 "
                "to chain=%r would cross a chain with itself (ADR-049)."
                % (method, chain)
            )
        raise ValueError(
            "cross method %r needs BOTH chains: pass chain_2=<the second chain>. "
            "The standard honours the chain named on each side and never scans "
            "an unrelated chain, so there is no chain it could default to "
            "(ADR-049). For %r, chain=%r is the %s side and chain_2 is the %s "
            "side." % (method, method, chain,
                       cross.side_b_role if cross.entity_side == "b"
                       else cross.side_a_role,
                       cross.side_a_role if cross.entity_side == "b"
                       else cross.side_b_role)
        )
    _require_chain_present(structure, chain_2, "chain_2")
    if cross.entity_side is None and chain == chain_2:
        # Both sides are protein: two different chains (ADR-109). A molecule
        # side may share its chain identifier with the protein, because there
        # the two sides are told apart by what they select.
        raise InvalidRequest(
            "method %r crosses two protein chains, and chain=%r and chain_2=%r "
            "name the same chain. Name two different chains; BioTSR 1.0 does not "
            "cross a protein chain with itself." % (method, chain, chain_2))

    # -- build the two side specs, routed by role ---------------------------
    entity_admission = admission_a if cross.entity_side == "a" else admission_b

    def entity_spec(entity_chain, role):
        if entity_names:
            names = _split(entity_names)
            if not any(_entity_instances(structure, entity_chain, entity_admission, n)
                       for n in names):
                raise InvalidRequest(
                    "none of the requested entities (%s) is in chain %s of %s: %s"
                    % (", ".join(names), entity_chain, structure.structure_id.value,
                       _absent(structure, entity_chain, entity_admission,
                               "/".join(names)))
                )
            return CrossSideSpec(
                selector=SIDE_SELECTOR_ENTITY_CHAIN,
                chain=entity_chain,
                role=role,
                entity_names=names,
            )
        if not entity_name:
            raise ValueError(
                "method %r needs the molecule it applies to: pass "
                "entity_name=<name>, and entity_id=<residue id> for one copy "
                "(omit it for every copy on the chain), e.g. entity_name='ESW', "
                "entity_id='404'; or entity_names=... for a whole-chain "
                "vocabulary" % method
            )
        from .domain.selection_record import SCOPE_ONE_INSTANCE

        refuse_unreadable(structure, entity_chain, entity_admission,
                          names=(entity_name,))
        instances, scope = _resolve_entity_instances(
            method, structure, entity_chain, entity_admission, entity_name,
            entity_id, what=role,
        )
        return CrossSideSpec(
            selector=SIDE_SELECTOR_ENTITY,
            chain=entity_chain,
            role=role,
            entity_name=str(entity_name).strip().upper(),
            # No identifier: every copy on the chain (``all_copies``).
            entity_id=instances[0] if scope == SCOPE_ONE_INSTANCE else None,
            atoms=_split(atoms),
        )

    def protein_spec(default_selector, protein_chain, role, entity_names_default,
                     r_types, r_specific, a_selection, ca_residues):
        selector = _protein_side_selector(
            default_selector, r_types, r_specific, a_selection, ca_residues
        )
        return CrossSideSpec(
            selector=selector,
            chain=protein_chain,
            role=role,
            entity_names=tuple(entity_names_default or ())
            if selector == default_selector else (),
            residue_types=_split(r_types),
            residues=_split(r_specific) or _split(ca_residues),
            atom_selection=a_selection or "",
        )

    entity_spec_used = None
    if cross.entity_side == "a":
        spec_a = entity_spec(chain, cross.side_a_role)
        entity_spec_used = spec_a
        spec_b = protein_spec(
            selector_b, chain_2, cross.side_b_role,
            cross.side_b_entity_names,
            residue_types, residues, atom_selection, calpha_residues,
        )
    elif cross.entity_side == "b":
        # The caller names the entity's chain first; the protein side is chain_2.
        spec_a = protein_spec(
            selector_a, chain_2, cross.side_a_role,
            cross.side_a_entity_names,
            residue_types, residues, atom_selection, calpha_residues,
        )
        spec_b = entity_spec(chain, cross.side_b_role)
        entity_spec_used = spec_b
    else:
        spec_a = protein_spec(
            selector_a, chain, cross.side_a_role,
            cross.side_a_entity_names,
            residue_types, residues, atom_selection, calpha_residues,
        )
        spec_b = protein_spec(
            selector_b, chain_2, cross.side_b_role,
            cross.side_b_entity_names,
            residue_types_2, residues_2, atom_selection_2, calpha_residues_2,
        )

    if method == "protein_protein_calpha_cross" and no_contact_filter:
        raise ValueError(
            "%r is defined by a contact cutoff: an interface residue is one "
            "with at least one non-hydrogen atom pair within it. Without a "
            "cutoff there is no interface and no method, so "
            "no_contact_filter=True is refused rather than reinterpreted."
            % method
        )

    # -- selection ----------------------------------------------------------
    from .cross.sides import SIDE_SELECTOR_CALPHA_SELECTION

    for spec in (spec_a, spec_b):
        if spec.selector == SIDE_SELECTOR_CALPHA_SELECTION:
            labels, _ = parse_residue_selection(";".join(spec.residues))
            check_range_insertion_codes(structure, spec.chain, list(spec.residues),
                                        labels)
    selection_a = select_side(
        structure, spec_a, active.label_scheme, admission=admission_a
    )
    selection_b = select_side(
        structure, spec_b, active.label_scheme, admission=admission_b
    )
    accounting = {
        "a": _cross_side_accounting(structure, spec_a, selection_a, admission_a),
        "b": _cross_side_accounting(structure, spec_b, selection_b, admission_b),
    }
    for side, spec, picked in (("a", spec_a, selection_a), ("b", spec_b, selection_b)):
        _refuse_cross_side_unk(structure, spec, picked,
                               admission_a if side == "a" else admission_b)
    _refuse_shared_atoms(structure, method, spec_a, spec_b, selection_a, selection_b)
    for side, spec, picked in (("a", spec_a, selection_a), ("b", spec_b, selection_b)):
        if not picked.points:
            notes = accounting[side][2]
            raise InvalidRequest(
                "the %s side of %r selected nothing in chain %s of %s%s"
                % (spec.role, method, spec.chain, structure.structure_id.value,
                   ": %s" % "; ".join(notes) if notes else
                   "; check the chain and the narrowing arguments for that side")
            )

    interface_a = interface_b = None
    if cross.interface_admission is not None:
        # Step 1 of the Ca rule: the heavy atoms -- side chains included -- that
        # decide which residues are interface residues. A separate selection from
        # the Ca points, because the two steps ask different questions: "does this
        # residue touch the other chain at all?" and "what point represents it?".
        #
        # These atoms are never labelled for a key: only their coordinates and
        # residue identities decide contact, and the points that reach the key
        # formula are the Ca points, labelled in the residue label space. So
        # they are read with a scheme that gives every atom the same placeholder
        # code -- an atom with no element label still makes contact, and
        # cannot stop a calculation whose keys never use its label.
        #
        # A narrowed side is narrowed here too: only the heavy atoms of its
        # selected residues can make a residue an interface residue, on either
        # side, so a partner residue touching only unselected residues is not
        # in the interface (ADR-109).
        interface_a, interface_b = (
            list(select_entity_atoms(
                structure, spec.chain, CONTACT_ONLY_LABELS,
                admission=cross.interface_admission,
                entity_seqs=_calpha_side_labels(spec),
                selector=SIDE_SELECTOR_PROTEIN_CHAIN,
                require_entity_name=False,
            ).points)
            for spec in (spec_a, spec_b))
        integrity.refuse_coincident_points(
            interface_a + interface_b,
            "%s of %s" % (method, structure.structure_id.value))

    # -- contact-pool policy ------------------------------------------------
    # A reference value, not a mathematical constant: 5 A for protein-protein,
    # 3.5 A for drug / carbohydrate / nucleotide-protein. Always overridable,
    # and the value actually used is recorded in SamplingProfile.contact_cutoff.
    cutoff = cross.default_cutoff if contact_cutoff is None else contact_cutoff
    if no_contact_filter:
        if contact_cutoff is not None:
            raise ValueError(
                "no_contact_filter=True and contact_cutoff=%r contradict each "
                "other: either scan for contacts at a cutoff, or use both "
                "selections whole" % (contact_cutoff,)
            )
        if contact_filter is not None:
            raise ValueError(
                "no_contact_filter=True and contact_filter=%r contradict each "
                "other: contact_filter chooses which side a contact scan narrows, "
                "and no_contact_filter asks for no scan" % (contact_filter,)
            )
        policy = cross_pools.ContactPoolPolicy(
            name="no_contact_filter",
            cutoff=None,
            contact_filter=cross_pools.ContactFilter.NONE,
        )
    else:
        policy = cross.standard_pool(cutoff)

    if contact_filter is not None:
        known = (
            cross_pools.ContactFilter.BOTH_SIDES,
            cross_pools.ContactFilter.SIDE_A_ONLY,
            cross_pools.ContactFilter.SIDE_B_ONLY,
        )
        if contact_filter not in known:
            raise ValueError(
                "unknown contact_filter %r; it is one of %s (use "
                "no_contact_filter=True for no scan at all)"
                % (contact_filter, ", ".join(repr(k) for k in known))
            )
        from dataclasses import replace as _replace

        policy = _replace(
            policy,
            contact_filter=contact_filter,
            name="%s+filter(%s)" % (policy.name, contact_filter),
        )

    integrity.refuse_coincident_points(
        list(selection_a.points) + list(selection_b.points),
        "%s of %s" % (method, structure.structure_id.value))
    run = run_cross(
        structure,
        method=method,
        spec_a=spec_a,
        spec_b=spec_b,
        points_a=list(selection_a.points),
        points_b=list(selection_b.points),
        selection_a=selection_a,
        selection_b=selection_b,
        profile=active,
        pool_policy=policy,
        interface_atoms_a=interface_a,
        interface_atoms_b=interface_b,
        candidate_strategy=STRATEGY_MIXED_SIDE_LEX,
        keep_triplets=keep_triplets,
        triplet_sink=triplet_sink,
        annotations={
            "base_profile": base.name,
            "selection_request": _cross_selection_record(
                spec_a, spec_b, selection_a, selection_b, accounting),
            # The molecule analysed: the copy named (NAME_ID), every copy on
            # the chain (NAME_allcopies), or a pooled side's vocabulary.
            "target": (
                "%s_%s" % (entity_spec_used.entity_name,
                           entity_spec_used.entity_id or "allcopies")
                if entity_spec_used is not None and entity_spec_used.entity_name
                else ";".join(entity_spec_used.entity_names)
                if entity_spec_used is not None else None
            ),
        },
    )
    result = run.result
    extra = [note for side in ("a", "b") for note in accounting[side][2]]
    if result.counts.n_triangles == 0:
        extra.append(
            "no triangles: no mixed-side triangle could be formed%s"
            % (" -- no atom pair across the two sides lies within %g A"
               % policy.cutoff if policy.cutoff is not None
               and not result.counts.n_points else "")
        )
    if extra:
        from dataclasses import replace as _replace

        result = _replace(result, warnings=tuple(result.warnings) + tuple(extra))
    return result
