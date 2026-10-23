"""The ``protein`` method: ordinary whole-chain protein C-alpha TSR.

A declarative recipe -- no arithmetic, no modifier logic, no engine.
``supported_modifiers`` declares the three Ca modifiers and ``sse``, which is
where the 18 protein secondary-structure categories are defined and the only
family they are defined for (ADR-036).  Requesting it on any other method
still raises, naming what *is* supported rather than silently no-opping.

Note what is *not* here: there is no per-modifier code and no per-combination
class.  ``supported_modifiers`` is a declaration; the resolution happens in
``Profile.with_modifiers`` and the execution in the one shared engine loop.
"""

from __future__ import annotations

from ...candidates.exhaustive import STRATEGY_NAME
from ...profiles.standard import STANDARD_PROTEIN_CALPHA_V1
from ...selection import protein as protein_selection
from ..registry import Method, register_method

__all__ = ["PROTEIN_METHOD"]

PROTEIN_METHOD = register_method(
    Method(
        name="protein",
        family="calpha",
        selector=protein_selection.SELECTOR_NAME,
        scope=protein_selection.SCOPE,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_PROTEIN_CALPHA_V1,
        supported_modifiers=frozenset(
            {"aa_grouping", "mirror", "size_filter", "sse"}
        ),
        completeness_policy=None,  # TSR.py has none; none is invented here.
        composition_rule=None,
        description=(
            "Ordinary regular protein Ca TSR: every Ca atom of one chain, all "
            "C(n,3) triangles, residue labels for the 20 regular amino acids and "
            "the supported modified residues (radix 30). Supports the three "
            "Ca modifiers -- aa_grouping, mirror, size_filter -- composably, and "
            "sse annotation on top of any of them (which never changes a key)."
        ),
    )
)


# ---------------------------------------------------------------------------
# aa_selection -- residue-position selection, still a Ca method
# ---------------------------------------------------------------------------

from ...profiles.standard import STANDARD_CALPHA_SELECTION_V1  # noqa: E402

__all__ += ["CALPHA_SELECTION_METHOD"]

CALPHA_SELECTION_METHOD = register_method(
    Method(
        name="aa_selection",
        family="calpha",
        selector=protein_selection.SELECTOR_CALPHA_SELECTION,
        scope=protein_selection.SCOPE,
        candidate_strategy=STRATEGY_NAME,
        default_profile=STANDARD_CALPHA_SELECTION_V1,
        supported_modifiers=frozenset(
            {"aa_grouping", "mirror", "size_filter", "sse"}
        ),
        completeness_policy=None,
        composition_rule=None,
        description=(
            "aa_selection: the Ca atoms of an explicitly listed set of residues "
            "in aa_chain ('10-12;48;70-71'). A Ca method, so it uses the residue "
            "labels at radix 30 and the 35-class Ca distance scheme -- "
            "residue-selection strategy never changes the label space (ADR-032). "
            "Supports the three Ca modifiers, aa_grouping included, and sse: the "
            "selected residues keep their own secondary-structure states, and the "
            "representation stays Ca/radix-30 rather than becoming atomic."
        ),
    )
)
