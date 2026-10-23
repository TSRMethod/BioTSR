"""Method registry: named, declarative recipes.

A method is a *recipe*, not an engine -- it pins a family's defaults and contains
no arithmetic.  The registry is intentionally a dict
plus validation, not a plugin framework with entry points and discovery: the
evidence supports a small, closed, curated set of scientific methods.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, FrozenSet, Optional, Tuple

from ..profiles.standard import Profile

__all__ = ["Method", "register_method", "get_method", "list_methods"]


@dataclass(frozen=True)
class Method:
    """A named scientific recipe."""

    name: str
    family: str
    selector: str
    scope: str
    candidate_strategy: str
    default_profile: Profile
    supported_modifiers: FrozenSet[str] = frozenset()
    completeness_policy: Optional[str] = None
    composition_rule: Optional[str] = None
    #: For a polymer/chain-scoped method: the entity-name vocabulary that defines
    #: what belongs to this family. Empty means "no name filter". Carried on the
    #: recipe rather than inside the selector, so the rule is visible at the
    #: method's declaration site (it is the carbohydrate scoping decision the
    #: legacy audit flagged, made explicit rather than hidden in a JSON file).
    default_entity_names: Tuple[str, ...] = ()
    #: Minimum / maximum number of entity types this method accepts, for the
    #: recipes that differ only in arity (``aa_ind`` vs ``aa_ind_mix``).
    min_entity_types: int = 0
    max_entity_types: Optional[int] = None
    description: str = ""

    def validate_modifiers(self, requested) -> None:
        """Reject any modifier this method does not declare.

        Declaration is necessary but not sufficient: a declared modifier can
        still be rejected at resolution time if the label domain in play does not
        support it (e.g. ``aa_grouping`` outside the protein amino-acid domain).
        Both refusals name what is supported; neither is ever a silent no-op.
        """
        unknown = sorted(set(requested) - set(self.supported_modifiers))
        if unknown:
            raise ValueError(
                "method %r does not support %s; it supports %s"
                % (
                    self.name,
                    ", ".join(repr(u) for u in unknown),
                    ", ".join(sorted(self.supported_modifiers)) or "no modifiers yet",
                )
            )


_REGISTRY: Dict[str, Method] = {}


def register_method(method: Method) -> Method:
    if method.name in _REGISTRY:
        raise ValueError("method %r is already registered" % method.name)
    _REGISTRY[method.name] = method
    return method


def get_method(name: str) -> Method:
    try:
        return _REGISTRY[name]
    except KeyError:
        raise KeyError(
            "unknown method %r; registered methods: %s"
            % (name, ", ".join(sorted(_REGISTRY)) or "(none)")
        )


def list_methods():
    return tuple(sorted(_REGISTRY))
