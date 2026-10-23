"""The analysis layer: feature matrices, common keys, similarity, clustering.

One implementation serves every TSR method family.  It consumes ``TSRResult``
objects -- or bare ``key -> count`` mappings -- and never asks which engine made
them, so regular Ca, selected Ca, grouped Ca, mirror, size-filtered, amino-acid
atomic, drug, carbohydrate, nucleotide and every cross family go through exactly
the same code.

What it *does* ask is whether the integers mean the same thing.  A shared
``KeySpace`` is a precondition of every operation that touches key coordinates
and is refused rather than coerced, always.  A difference in processing (size
filter, contact cutoff or rule, completeness rule, candidate strategy,
alternate-location rule) is refused unless the caller passes
``allow_mixed_processing=True``, and then recorded.  What was selected --
structures, chains, residue types, molecules, residue numbers, atoms, and the
method name -- never blocks a comparison (``samples.check_comparable``).

``render`` is not imported here.  It is the only module in the package that
touches matplotlib (the ``plot`` extra), and importing ``biotsr.analysis`` must
keep working where it is not installed.
"""

from __future__ import annotations

from .clustering import (
    AVERAGE,
    ClusterAssignment,
    Hierarchy,
    Linkage,
    adjusted_rand_index,
    average_linkage,
    cluster,
    cut_hierarchy,
    flat_clusters,
    hierarchy,
    leaf_order,
)
from .common_keys import (
    CommonKeyIntersection,
    KeyOccurrenceTable,
    OneVsAllOccurrence,
    common_key_intersection,
    common_key_intersection_by_group,
    key_occurrence_percent,
    key_occurrence_percent_one_vs_all,
)
from .features import (
    FEATURE_ORDERINGS,
    FIRST_ENCOUNTER,
    NUMERIC_ASCENDING,
    FeatureMatrix,
    feature_matrix,
)
from .samples import (
    PROCESSING_FIELDS,
    TARGET_FIELDS,
    AnalysisSample,
    IncompatibleProcessing,
    PartialSamples,
    check_comparable,
    describe_omissions,
    label_from_name,
    processing_fingerprint,
    processing_policy,
    require_common_key_space,
    require_common_processing,
    require_common_structure_interpretation,
    require_common_sampling,
    sampling_variation,
)
from .similarity import (
    DISTANCE,
    MEASURES,
    SIMILARITY,
    MeasureSpec,
    SimilarityMatrix,
    binary_jaccard,
    cosine,
    distance_matrix,
    generalized_jaccard,
    sarika,
    similarity_matrix,
    wu,
)
from .visualization import ClustermapSpec, clustermap_spec

__all__ = [
    # samples
    "AnalysisSample",
    "IncompatibleProcessing",
    "PartialSamples",
    "describe_omissions",
    "PROCESSING_FIELDS",
    "TARGET_FIELDS",
    "check_comparable",
    "processing_policy",
    "processing_fingerprint",
    "label_from_name",
    "require_common_key_space",
    "require_common_processing",
    "require_common_structure_interpretation",
    "require_common_sampling",
    "sampling_variation",
    # features
    "FeatureMatrix",
    "feature_matrix",
    "FEATURE_ORDERINGS",
    "NUMERIC_ASCENDING",
    "FIRST_ENCOUNTER",
    # common keys
    "CommonKeyIntersection",
    "KeyOccurrenceTable",
    "OneVsAllOccurrence",
    "common_key_intersection",
    "common_key_intersection_by_group",
    "key_occurrence_percent",
    "key_occurrence_percent_one_vs_all",
    # similarity
    "SIMILARITY",
    "DISTANCE",
    "MEASURES",
    "MeasureSpec",
    "SimilarityMatrix",
    "similarity_matrix",
    "distance_matrix",
    "generalized_jaccard",
    "binary_jaccard",
    "wu",
    "sarika",
    "cosine",
    # clustering
    "AVERAGE",
    "Linkage",
    "Hierarchy",
    "ClusterAssignment",
    "average_linkage",
    "hierarchy",
    "cut_hierarchy",
    "flat_clusters",
    "leaf_order",
    "cluster",
    "adjusted_rand_index",
    # visualization
    "ClustermapSpec",
    "clustermap_spec",
]
