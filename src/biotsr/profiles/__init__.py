"""Named bundles of scientific choices.

Every method runs under exactly one profile, its standard one
(``biotsr.get_method(name).default_profile``); profiles are how the package
states its conventions, not a user-facing choice.
"""

from .standard import STANDARD_PROTEIN_CALPHA_V1, Profile

__all__ = [
    "STANDARD_PROTEIN_CALPHA_V1",
    "Profile",
]
