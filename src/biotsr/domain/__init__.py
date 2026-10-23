"""Pure value objects.  No I/O, no algorithms, no imports from the rest of ``biotsr``."""

from .candidate import TriangleCandidate
from .identity import AtomIdentity, ChainId, ResidueIdentity, StructureId
from .interpretation import (
    LEGACY_COLUMNS_V1,
    SPEC_COLUMNS_V1,
    STRUCTURE_INTERPRETATION,
    IncompatibleStructureInterpretation,
    interpretation_of_backend,
    require_common_structure_interpretation,
    structure_interpretation_of,
)
from .keyspace import (
    IncompatibleKeySpace,
    IncompatibleSampling,
    KeySpace,
    SamplingProfile,
    SchemeRef,
)
from .point import LabeledPoint, Point, PointSet
from .result import (
    DEFAULT_CATEGORY,
    CategoryCounts,
    ResultCounts,
    SelectionSummary,
    SkipRecord,
    SourceInfo,
    TSRResult,
    TripletRecord,
)

__all__ = [
    "AtomIdentity",
    "IncompatibleStructureInterpretation",
    "LEGACY_COLUMNS_V1",
    "SPEC_COLUMNS_V1",
    "STRUCTURE_INTERPRETATION",
    "interpretation_of_backend",
    "require_common_structure_interpretation",
    "structure_interpretation_of",
    "CategoryCounts",
    "ChainId",
    "DEFAULT_CATEGORY",
    "IncompatibleKeySpace",
    "IncompatibleSampling",
    "KeySpace",
    "LabeledPoint",
    "Point",
    "PointSet",
    "ResidueIdentity",
    "ResultCounts",
    "SamplingProfile",
    "SchemeRef",
    "SelectionSummary",
    "SkipRecord",
    "SourceInfo",
    "StructureId",
    "TSRResult",
    "TriangleCandidate",
    "TripletRecord",
]
