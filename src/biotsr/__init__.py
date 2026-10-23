"""BioTSR -- the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology, developed by the TSR Research Group.

The shared core; ordinary protein C-alpha TSR with its three composable
modifiers (AA grouping, mirror, size filtering); the shared atomic foundation, and
the drug, carbohydrate, nucleotide and cross families built on it; the analysis
layer -- feature matrices, common keys, similarity and hierarchical clustering --
in ``biotsr.analysis``; and secondary-structure annotation in ``biotsr.sse``.

SSE is an *annotation*, not a second engine: ``analyze(..., sse=True)`` produces
the same integer keys as the same call without it, and changes only which of the
18 categories each key's frequency is counted in.  It is offered by the two
C-alpha methods and refused, by name, everywhere the categories are undefined.

A result can be saved and **loaded back**: ``biotsr.io.save`` writes the native
format and ``biotsr.load_result`` reads it, restoring both scientific identities
with their exact fingerprints, so ``biotsr.analysis`` runs over previously computed
results without re-reading a single PDB file.  What a loaded result deliberately
does not carry -- triplet records and first-seen key order -- is stated on the
object and refused by the operations that would need it (ADR-075).

The ``biotsr`` command line is part of the package: ``biotsr ...``, or ``python -m biotsr``
where the environment's ``bin`` directory is not on ``PATH``.  ``biotsr.batch``
(``from biotsr import batch``; it is not imported for you) runs a whole dataset manifest from Python with the same engine, and
``biotsr.fetch_structures`` downloads PDB-format entries when asked -- the only
part of the package that uses the network.

Candidate-strategy extensions (for example CA-kNN) are not part of this release,
and are deliberately absent rather than stubbed.

The core imports nothing outside the Python standard library, and neither does
the analysis layer.  ``biotsr.analysis.render`` is the single exception and is never
imported for you: it draws the historical clustermap with matplotlib and seaborn,
and nothing else in the package depends on it.
"""

from . import analysis
from .analysis.samples import IncompatibleProcessing
from . import sse
from .api import (
    analyze,
    analyze_cross,
    analyze_entities,
    list_cross_methods,
    list_methods,
    list_entity_targets,
    list_residue_targets,
    load_structure,
)
from .domain import (
    AtomIdentity,
    CategoryCounts,
    ChainId,
    IncompatibleKeySpace,
    IncompatibleSampling,
    KeySpace,
    LabeledPoint,
    Point,
    PointSet,
    ResidueIdentity,
    SamplingProfile,
    StructureId,
    TSRResult,
    TriangleCandidate,
    TripletRecord,
)
from .domain.errors import (
    IncompatibleRepresentation,
    InvalidRequest,
    MalformedResult,
    TSRError,
    UnsupportedFormat,
)
from .domain.result import RESULT_SCHEMA_VERSION
from .geometry.triangle import GeometryResult
from .io.native_read import load_result
from .labels.elements import ELEMENT_LABELS_112
from .labels.schemes import (
    AMINO_ACID_30,
    AMINO_ACID_GROUP_23,
    LabelDomain,
    LabelScheme,
)
from .modifiers import (
    MIRROR_DISABLED,
    MIRROR_NEGATE_WHEN_PREFOLD_GT_90,
    IncompatibleModifier,
    MirrorConvention,
    SizeFilter,
)
from .methods.registry import get_method
from .profiles.atomic import (
    STANDARD_AMINOACID_ATOM_SELECTION_V1,
    STANDARD_AMINOACID_ATOMIC_V1,
)
from .profiles.cross import STANDARD_ATOMIC_CROSS_V1
from .profiles.standard import STANDARD_PROTEIN_CALPHA_V1, Profile
from .retrieval import FetchRecord, FetchReport, fetch_structures
from .sse.stratification import IncompatibleSSEStratification

__version__ = "1.0.0rc7"

__all__ = [
    "analysis",
    "sse",
    "analyze",
    "analyze_cross",
    "analyze_entities",
    "list_entity_targets",
    "list_residue_targets",
    "AMINO_ACID_30",
    "ELEMENT_LABELS_112",
    "STANDARD_AMINOACID_ATOMIC_V1",
    "STANDARD_AMINOACID_ATOM_SELECTION_V1",
    "STANDARD_ATOMIC_CROSS_V1",
    "AMINO_ACID_GROUP_23",
    "IncompatibleModifier",
    "LabelDomain",
    "LabelScheme",
    "MIRROR_DISABLED",
    "MIRROR_NEGATE_WHEN_PREFOLD_GT_90",
    "MirrorConvention",
    "SizeFilter",
    "load_structure",
    "load_result",
    "fetch_structures",
    "FetchReport",
    "FetchRecord",
    "RESULT_SCHEMA_VERSION",
    "get_method",
    "TSRError",
    "InvalidRequest",
    "IncompatibleRepresentation",
    "IncompatibleSSEStratification",
    "MalformedResult",
    "UnsupportedFormat",
    "list_cross_methods",
    "list_methods",
    "AtomIdentity",
    "CategoryCounts",
    "ChainId",
    "GeometryResult",
    "IncompatibleKeySpace",
    "IncompatibleSampling",
    "IncompatibleProcessing",
    "KeySpace",
    "LabeledPoint",
    "Point",
    "PointSet",
    "Profile",
    "ResidueIdentity",
    "SamplingProfile",
    "STANDARD_PROTEIN_CALPHA_V1",
    "StructureId",
    "TSRResult",
    "TriangleCandidate",
    "TripletRecord",
    "__version__",
]
