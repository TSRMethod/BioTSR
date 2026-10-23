"""The standard profile for protein C-alpha TSR, and modifier resolution.

Everything here is the scientific default:

* ``math.pi``, never the legacy divisor;
* the proper 35-class Ca distance scheme (``d >= 132.45 -> 35``);
* the 29-class theta scheme, fed the folded theta;
* the residue labels -- the 20 regular amino acids and the supported modified
  residues, 30 codes at radix **30** -- unless AA grouping is requested, which
  swaps in the 23 grouped codes, also at radix **30** (ADR-107);
* mirror signing off, so no key is ever negated;
* no size filter;
* degenerate triples skipped and recorded, never forced to ``theta = 0``.

There is no configuration in which a profile produces the historical
programs' arithmetic.

**Modifiers are resolved into a profile, not into the engine.**
``Profile.with_modifiers(...)`` returns a new frozen profile whose label scheme,
mirror convention and size filter reflect the request, and whose ``name`` records
what was applied.  That is what makes composition free: the engine has one code
path and never asks which modifiers are active, and there is no combinatorial
family of per-combination classes.

Binning is fixed per family and is not a modifier: **Ca calculations always use
29 theta classes and 35 Ca distance classes**; every non-Ca atomic family uses 29
theta and 18 distance classes (ADR-021).  Neither may be substituted for the
other -- and because the scheme's content fingerprint enters ``KeySpace``, a
substitution cannot silently produce comparable keys.

``Profile`` is shared by every family.  The atomic families' profiles live in
``profiles/atomic.py`` and differ from the Ca ones only in configuration: a
different distance scheme, a different label scheme, an atom-admission policy and
a completeness policy.  No family has its own ``Profile`` subclass.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

from ..binning.base import LadderBinScheme
from ..binning.distance import DISTANCE_SCHEME_CALPHA_35
from ..binning.theta import THETA_SCHEME_29
from ..domain.keyspace import KeySpace
from ..encoding.key import KEY_FORMULA_VERSION, KeyFormulaSpec
from ..geometry.conventions import (
    DEFAULT_ACOS_TOLERANCE,
    MEDIAN_FROM_MIDPOINT,
    SQRT_POW_HALF,
    AngleConvention,
    DegeneratePolicy,
    MedianConvention,
    PI_CONVENTION,
    SqrtConvention,
)
from ..geometry.ordering import STANDARD_ORDERING, OrderingConvention
from ..labels.schemes import AMINO_ACID_30, AMINO_ACID_GROUP_23, LabelScheme
from ..modifiers.aa_grouping import apply_grouping
from ..selection.atomic import AtomAdmission
from ..selection.base import SeqLabelPolicy
from ..selection.completeness import CompletenessPolicy, NO_COMPLETENESS_POLICY
from ..modifiers.mirror import (
    MIRROR_DISABLED,
    MIRROR_NEGATE_WHEN_PREFOLD_GT_90,
    MirrorConvention,
)
from ..modifiers.size_filter import SizeFilter
from ..sse.assignment import STRAND_IDENTITY_SHEET_ID_AND_SERIAL
from ..structure.backends import pdb_legacy

__all__ = [
    "Profile",
    "STANDARD_PROTEIN_CALPHA_V1",
    "STANDARD_CALPHA_SELECTION_V1",
]


@dataclass(frozen=True)
class Profile:
    """A named bundle of every scientific choice a run makes.

    Selecting a profile is a single explicit act, and the choice propagates into
    the ``KeySpace`` fingerprint, so results from different profiles can never be
    compared by accident.
    """

    name: str
    angle_convention: AngleConvention
    theta_scheme: LadderBinScheme
    distance_scheme: LadderBinScheme
    label_scheme: LabelScheme
    parser_backend: str
    #: How precisely the reader keeps the file's coordinates. Not a tuning knob:
    #: ``Bio.PDB`` stores coordinates as float32, so the two TSR-Package cross
    #: programs -- the only implementations in the family that read through it --
    #: work from single-precision coordinates. It is a property of how the
    #: structure was READ, so it is recorded in provenance
    #: (``SourceInfo.parser_backend``) rather than in ``KeySpace``: it does not
    #: change what a key means, it changes the input the key is computed from, the
    #: same way a different structure file would.
    #:
    #: **Superseded (ADR-053).**  It is now
    #: key-affecting and enters ``KeySpace`` as ``numerical_precision``, because
    #: the difference was measured rather than assumed: on 1BRS A/D, 2 of
    #: 1,230,714 triangles change theta class between the two precisions. Two
    #: runs that differ only here can disagree on a key, which is exactly what
    #: ``KeySpace`` exists to make impossible to combine by accident. The
    #: standard is ``float64``; only the ``TSR-Package`` cross compatibility
    #: profiles ask for ``float32``.
    coordinate_precision: str = pdb_legacy.COORD_PRECISION_DOUBLE
    #: How a squared distance becomes a distance. Key-affecting (see
    #: geometry/conventions.py), so it enters KeySpace with the angle convention.
    sqrt_convention: SqrtConvention = SQRT_POW_HALF
    #: How the base-to-apex median is computed. Key-affecting; only the two
    #: TSR-Package cross programs deviate (ADR-039).
    median_convention: MedianConvention = MEDIAN_FROM_MIDPOINT
    #: Which lineage's spelling of the four-case canonical ordering applies.
    #: Key-affecting; only the two TSR-Package cross programs deviate, by
    #: comparing labels as strings and by a chained-inequality branch
    #: (ADR-038).
    ordering_convention: OrderingConvention = STANDARD_ORDERING
    grouped_label_scheme: LabelScheme = AMINO_ACID_GROUP_23
    mirror_convention: MirrorConvention = MIRROR_DISABLED
    size_filter: Optional[SizeFilter] = None
    aa_grouping: bool = False
    # Selection-side policies. They decide which atoms and which entities take
    # part, never what a key means, so they enter `SamplingProfile` and not
    # `KeySpace`. Protein Ca leaves both at their defaults: it admits atoms by
    # its own selector's rules and has no completeness concept at all.
    atom_admission: Optional[AtomAdmission] = None
    completeness: CompletenessPolicy = NO_COMPLETENESS_POLICY
    #: Which residues the Ca selector admits by the shape of their sequence
    #: label. Selection-side, so it enters `SamplingProfile`. Every residue is
    #: admitted (`accept_all`).
    seq_label_policy: str = SeqLabelPolicy.ACCEPT_ALL
    #: Which PDB fields name an individual β-strand when a run is annotated with
    #: `sse=True`: `(sheetID, serNum)` (ADR-063).
    #:
    #: It is in NEITHER identity, and that is the point: it moves no key, so it
    #: is not a `KeySpace` field, and it draws no triangle, so it is not a
    #: `SamplingProfile` field either -- an annotated run and its plain twin
    #: must keep sharing both fingerprints. It travels in the result's SSE
    #: **annotation provenance** (`annotations["sse"]["strand_identity"]`).
    sse_strand_identity: str = STRAND_IDENTITY_SHEET_ID_AND_SERIAL
    degenerate_policy: str = DegeneratePolicy.SKIP_AND_RECORD
    clamp_acos: bool = True
    acos_tolerance: float = DEFAULT_ACOS_TOLERANCE
    #: False for every profile BioTSR provides. Recorded in result provenance
    #: (``legacy_compatibility_profile``) so the result schema is unchanged.
    is_legacy_compatibility: bool = False
    description: str = ""

    # -- modifier resolution -------------------------------------------------

    def with_modifiers(
        self,
        *,
        aa_grouping: bool = False,
        mirror: bool = False,
        size_filter=None,
        method: str = "",
    ) -> "Profile":
        """Return a new profile with the requested modifiers applied.

        Each modifier attaches at its own hook (``modifiers/base.py``), so this
        is the only place their configuration meets.  Nothing here branches on
        *combinations* -- the eight combinations are eight independent yes/no
        answers, which is precisely why no combinatorial classes are needed.
        """
        changes = {}
        applied = []

        if aa_grouping:
            changes["label_scheme"] = apply_grouping(
                self.label_scheme, self.grouped_label_scheme, method=method
            )
            changes["aa_grouping"] = True
            applied.append("aa_grouping")

        if mirror:
            changes["mirror_convention"] = MIRROR_NEGATE_WHEN_PREFOLD_GT_90
            applied.append("mirror")

        resolved_filter = SizeFilter.from_value(size_filter)
        if resolved_filter is not None:
            changes["size_filter"] = resolved_filter
            applied.append("size_filter(%g)" % resolved_filter.cutoff)

        if not changes:
            return self
        changes["name"] = "%s+%s" % (self.name, "+".join(applied))
        return replace(self, **changes)

    # -- derived identities --------------------------------------------------

    def key_formula(self) -> KeyFormulaSpec:
        return KeyFormulaSpec(
            n_labels=self.label_scheme.n_labels,
            n_theta=self.theta_scheme.n_classes,
            n_distance=self.distance_scheme.n_classes,
            label_range=self.label_scheme.code_range,
        )

    def key_space(self) -> KeySpace:
        """Everything that changes what an integer key means.

        AA grouping enters here through a different label scheme -- another
        name and content fingerprint, at the same radix 30 -- and mirror enters
        as the sign convention.  The size filter deliberately does not: it changes which
        triangles were drawn, not what a key means.
        """
        return KeySpace(
            angle_convention=self.angle_convention.content_fingerprint(),
            sqrt_convention=self.sqrt_convention.content_fingerprint(),
            median_convention=self.median_convention.content_fingerprint(),
            numerical_precision=self.coordinate_precision,
            theta_scheme=self.theta_scheme.as_ref(),
            distance_scheme=self.distance_scheme.as_ref(),
            label_scheme=self.label_scheme.as_ref(),
            canonical_ordering_version=self.ordering_convention.version,
            key_formula_version=KEY_FORMULA_VERSION,
            mirror_convention=self.mirror_convention.name,
        )

    def degenerate_policy_label(self) -> str:
        """The sampling-profile label, including whether acos clamping was on.

        Clamping changes which candidates survive (never what a key means), so it
        belongs with the degenerate policy in ``SamplingProfile``, not in
        ``KeySpace``.
        """
        return "%s%s" % (
            self.degenerate_policy,
            "+acos_clamp" if self.clamp_acos else "",
        )

    def active_modifiers(self):
        """The modifier names in effect, for provenance. Deterministically ordered."""
        active = []
        if self.aa_grouping:
            active.append("aa_grouping")
        if self.mirror_convention.enabled:
            active.append("mirror")
        if self.size_filter is not None:
            active.append("size_filter")
        return tuple(active)


STANDARD_PROTEIN_CALPHA_V1 = Profile(
    name="standard_protein_calpha_v1",
    angle_convention=PI_CONVENTION,
    theta_scheme=THETA_SCHEME_29,
    distance_scheme=DISTANCE_SCHEME_CALPHA_35,
    label_scheme=AMINO_ACID_30,
    grouped_label_scheme=AMINO_ACID_GROUP_23,
    parser_backend=pdb_legacy.BACKEND_NAME,
    mirror_convention=MIRROR_DISABLED,
    size_filter=None,
    degenerate_policy=DegeneratePolicy.SKIP_AND_RECORD,
    clamp_acos=True,
    is_legacy_compatibility=False,
    description=(
        "Standard ordinary protein Ca TSR: math.pi, 35-class Ca distance "
        "scheme, 29-class theta, residue labels for the 20 regular amino acids "
        "and the supported modified residues (30 codes, radix 30), no mirror "
        "signing, no size filter. AA grouping swaps in the 23 grouped codes, "
        "also at radix 30. "
        "Structures are read by pdb_legacy, the one PDB reader BioTSR 1.0 ships; "
        "the reader affects parsing only, never key semantics, and is recorded "
        "in provenance."
    ),
)


# ---------------------------------------------------------------------------
# aa_selection -- residue-position selection in the Ca family
# ---------------------------------------------------------------------------

#: The standard for ``aa_selection``: a protein Cα point **is a residue**, so it is
#: labelled in the residue label space at radix 30, exactly like ordinary
#: whole-chain Cα TSR.  Candidate/residue-selection strategy never changes the
#: label space (ADR-032).
STANDARD_CALPHA_SELECTION_V1 = replace(
    STANDARD_PROTEIN_CALPHA_V1,
    name="standard_calpha_selection_v1",
    description=(
        "Standard aa_selection: the Ca atoms of an explicitly listed set of "
        "residues, labelled by residue identity in the residue label space at "
        "radix 30, 35-class Ca distance scheme, 29-class theta. Identical key "
        "space to ordinary protein Ca TSR -- only which residues take part differs, "
        "which is a SamplingProfile concern and not a KeySpace one."
    ),
)
