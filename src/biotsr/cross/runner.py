"""The cross runner: two sides -> contact pool -> shared core -> Result.

This is the only module that knows a cross has two sides.  Everything it calls
downstream of pool construction is the same code path an ordinary one-sided
method uses -- ``run_core`` receives a ``PointSet`` and a ``Profile``, exactly as
before, plus the mixed-side candidate iterator and its exact count.

The flow, with no lossy intermediate:

    structure parsing
        -> side-A selection        (an existing selector, unchanged)
        -> side-B selection        (an existing selector, unchanged)
        -> contact-pool construction under full AtomIdentity
        -> mixed-side candidate enumeration
        -> shared ordering / geometry / binning / encoding
        -> Result

A contact CSV is never required to get here.  The legacy pipelines route the pool
through ``drug_protein_3A.py``'s CSV, which is exactly where physical identity is
lost -- the CSV keeps chain, residue and atom name, but ``cross_key.py`` then
re-keys on the atom name alone.  Here the pool is built in memory from
``AtomIdentity`` values (ADR-042).
"""

from __future__ import annotations

from typing import Callable, Optional, Sequence

from ..domain.keyspace import SamplingProfile
from ..domain.point import LabeledPoint
from ..domain.result import (
    SelectionSummary,
    SourceInfo,
    TSRResult,
    TripletRecord,
)
from ..methods.engine import run_core
from ..profiles.standard import Profile
from ..selection.base import AltLocPolicy, InputRules, ModifiedResiduePolicy
from ..structure.model import Structure
from .candidates import (
    COMPOSITION_RULE,
    STRATEGY_MIXED_SIDE,
    count,
    mixed_side,
)
from .pools import (
    ContactPoolPolicy,
    CrossPools,
    PoolProjection,
    build_contact_pools,
    build_interface_calpha_pools,
)
from .sides import CrossSideSpec

__all__ = ["CrossRun", "run_cross"]


_STRATEGIES = {
    STRATEGY_MIXED_SIDE: mixed_side,
}


class CrossRun:
    """A completed cross calculation: the result plus the pool ledger."""

    __slots__ = ("result", "pools", "selection_a", "selection_b")

    def __init__(self, result, pools, selection_a, selection_b):
        self.result = result
        self.pools = pools
        self.selection_a = selection_a
        self.selection_b = selection_b


def run_cross(
    structure: Structure,
    *,
    method: str,
    spec_a: CrossSideSpec,
    spec_b: CrossSideSpec,
    points_a: Sequence[LabeledPoint],
    points_b: Sequence[LabeledPoint],
    selection_a,
    selection_b,
    profile: Profile,
    pool_policy: ContactPoolPolicy,
    interface_atoms_a: Optional[Sequence[LabeledPoint]] = None,
    interface_atoms_b: Optional[Sequence[LabeledPoint]] = None,
    candidate_strategy: str = STRATEGY_MIXED_SIDE,
    keep_triplets: bool = True,
    triplet_sink: Optional[Callable[[TripletRecord], None]] = None,
    annotations: Optional[dict] = None,
) -> CrossRun:
    """Build the contact pool and run the shared core over mixed-side candidates."""
    if candidate_strategy not in _STRATEGIES:
        raise ValueError(
            "unknown cross candidate strategy %r; known strategies: %s"
            % (candidate_strategy, ", ".join(sorted(_STRATEGIES)))
        )

    if pool_policy.projection == PoolProjection.HEAVY_ATOM_INTERFACE_TO_CALPHA:
        # Two selections per side: the heavy atoms decide which residues are
        # interface residues, the Ca points are what the interface residues
        # contribute. Neither substitutes for the other.
        if interface_atoms_a is None or interface_atoms_b is None:
            raise ValueError(
                "pool policy %r classifies interface residues by heavy-atom "
                "contact, so it needs the heavy atoms of both sides; pass "
                "interface_atoms_a= and interface_atoms_b=" % (pool_policy.name,)
            )
        pools: CrossPools = build_interface_calpha_pools(
            interface_atoms_a, interface_atoms_b, points_a, points_b, pool_policy
        )
    else:
        if interface_atoms_a is not None or interface_atoms_b is not None:
            raise ValueError(
                "pool policy %r projects at the atom level -- the contacting atom "
                "IS the vertex -- so a separate interface selection has no meaning"
                % (pool_policy.name,)
            )
        pools = build_contact_pools(points_a, points_b, pool_policy)
    sides = pools.sides

    n_candidates = count(pools.side_a_pool, pools.side_b_pool)
    candidates = _STRATEGIES[candidate_strategy](sides)

    core = run_core(
        pools.points,
        profile,
        candidates=candidates,
        candidate_count=n_candidates,
        keep_triplets=keep_triplets,
        triplet_sink=triplet_sink,
    )

    sampling = SamplingProfile(
        selector="%s|%s" % (spec_a.selector, spec_b.selector),
        scope="cross_two_sided",
        candidate_strategy=candidate_strategy,
        completeness_policy=None,
        composition_rule=COMPOSITION_RULE,
        size_filter=profile.size_filter.cutoff if profile.size_filter else None,
        degenerate_policy=profile.degenerate_policy_label(),
        cross_method=method,
        side_a_chain=spec_a.chain_label,
        side_b_chain=spec_b.chain_label,
        side_a_selector=spec_a.describe(),
        side_b_selector=spec_b.describe(),
        contact_cutoff=pool_policy.cutoff,
        contact_filter=pool_policy.contact_filter,
        contact_filter_rule=pool_policy.contact_filter_rule,
        contact_projection=pool_policy.projection,
        pool_identity=pool_policy.identity,
        pool_order=pool_policy.order,
    )

    chains_requested = tuple(
        dict.fromkeys((spec_a.chain, spec_b.chain))
    )
    notes = {
        "profile": profile.name,
        "legacy_compatibility_profile": profile.is_legacy_compatibility,
        "modifiers": profile.active_modifiers(),
        "cross_pools": pools.to_dict(),
        "side_a": spec_a.describe(),
        "side_b": spec_b.describe(),
    }
    if annotations:
        notes.update(annotations)

    result = TSRResult(
        method=method,
        structure_id=structure.structure_id,
        chain="%s|%s" % (spec_a.chain, spec_b.chain),
        keys=core.as_category_counts(),
        key_order=core.key_order,
        triplets=core.triplets,
        key_space=profile.key_space(),
        sampling=sampling,
        selection=SelectionSummary(
            selector="%s|%s" % (spec_a.selector, spec_b.selector),
            chains_requested=chains_requested,
            chains_used=tuple(
                dict.fromkeys(
                    tuple(selection_a.chains_used) + tuple(selection_b.chains_used)
                )
            ),
            n_points=len(pools.points),
            residues_used=tuple(selection_a.residues_used)
            + tuple(selection_b.residues_used),
            atoms_missing=tuple(selection_a.atoms_missing)
            + tuple(selection_b.atoms_missing),
        ),
        source=SourceInfo(
            structure_id=structure.structure_id,
            path=structure.source_path,
            sha256=structure.source_sha256,
            parser_backend=structure.parser_backend,
            model=1,
            altloc_policy=AltLocPolicy.BLANK_OR_A,
            modified_residue_policy=ModifiedResiduePolicy.for_label_domain(
                profile.label_scheme.domain),
            input_rules=InputRules.CURRENT,
        ),
        counts=core.counts,
        skipped=tuple(selection_a.skipped)
        + tuple(selection_b.skipped)
        + tuple(pools.skipped)
        + tuple(core.skipped),
        warnings=tuple(structure.warnings)
        + tuple(selection_a.warnings)
        + tuple(selection_b.warnings),
        annotations=notes,
    )
    return CrossRun(result, pools, selection_a, selection_b)
