"""The clustermap, described as data before it is drawn.

Visualization belongs to the analysis layer, but it must not *define* the
analysis.  So the figure is specified here as a plain, testable object -- which
matrix, which tree on rows and on columns, which leaf order, which palette slot
each row's annotation stripes take, the figure size -- and ``render.py`` turns
that specification into an image.  Everything a reader needs to check is
therefore checkable without matplotlib installed.

The default figure is the whole hierarchical result and needs **no number of
clusters**: the distance matrix as a heat map, ordered by the average-linkage
tree, with the *same* tree drawn on rows and columns (ADR-099).  Two stripes are
optional and kept distinct:

* **observed class** -- labels the caller supplies from their own table (a
  manifest's ``label`` column, a CSV), matched by sample id.  Nothing here
  derives a class from a file or sample name;
* **cluster** -- a flat partition, drawn only when the caller asks for one with
  an explicit cutting rule: ``n_clusters`` (``maxclust``) or ``height`` (cut
  every merge above it).  The rule is never inferred, and must not be chosen
  from the observed classes: a cut picked to match the classes agrees with them
  by construction.  The adjusted Rand index is reported only for a defined
  partition and fully labelled samples.

Stripes are coloured from matplotlib's ``tab10`` (classes) and ``Set2``
(clusters) over the *sorted* unique values (an evenly sampled ``turbo`` when
there are more values than a palette has colours).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Optional, Sequence, Tuple

from .clustering import (
    ClusterAssignment,
    Hierarchy,
    adjusted_rand_index,
    cut_hierarchy,
    hierarchy,
)
from .similarity import DISTANCE, SimilarityMatrix

__all__ = [
    "ClustermapSpec",
    "clustermap_spec",
    "DEFAULT_FIGSIZE",
    "ROW_NAMES_HEADING",
]

#: Inches.  Fits about 60 readable row labels; pass ``figsize`` for more.
DEFAULT_FIGSIZE = (10, 8)
#: The heading of :meth:`ClustermapSpec.to_row_names_csv`.
ROW_NAMES_HEADING = "sample_id"


@dataclass(frozen=True)
class ClustermapSpec:
    """Everything needed to draw the clustermap, and nothing computed by drawing it.

    ``matrix`` is drawn as given: its ``scale`` (fractions or percent) is the
    colour-bar scale, and ``matrix.value_label`` names it.  ``assignment`` is
    ``None`` unless a flat partition was asked for.
    """

    matrix: SimilarityMatrix
    hierarchy: Hierarchy
    class_labels: Tuple[Optional[str], ...]
    assignment: Optional[ClusterAssignment] = None
    class_palette: str = "tab10"
    cluster_palette: str = "Set2"
    row_cluster: bool = True
    col_cluster: bool = True
    figsize: Tuple[float, float] = DEFAULT_FIGSIZE
    cmap: str = "viridis"
    title: Optional[str] = None
    adjusted_rand: Optional[float] = None
    annotations: Mapping[str, Any] = field(default_factory=dict)

    @property
    def sample_ids(self) -> Tuple[str, ...]:
        return self.hierarchy.sample_ids

    @property
    def has_classes(self) -> bool:
        return any(c is not None for c in self.class_labels)

    @property
    def has_clusters(self) -> bool:
        return self.assignment is not None

    @property
    def unique_classes(self) -> Tuple[str, ...]:
        """Sorted, so a class keeps its colour whatever the sample order."""
        return tuple(sorted({c for c in self.class_labels if c is not None}))

    @property
    def unique_clusters(self) -> Tuple[int, ...]:
        return tuple(sorted(set(self._partition().labels)))

    def _partition(self) -> ClusterAssignment:
        if self.assignment is None:
            raise ValueError(
                "this clustermap has no flat partition: pass n_clusters= or "
                "height= to clustermap_spec to cut the tree by an explicit rule"
            )
        return self.assignment

    def class_palette_index(self) -> Mapping[str, int]:
        """Which slot of the class palette each class takes."""
        return {name: i for i, name in enumerate(self.unique_classes)}

    def cluster_palette_index(self) -> Mapping[int, int]:
        return {label: i for i, label in enumerate(self.unique_clusters)}

    def row_class_slots(self) -> Tuple[Optional[int], ...]:
        """Per row, in input order, the class stripe's palette slot (``None``: unlabelled)."""
        index = self.class_palette_index()
        return tuple(None if c is None else index[c] for c in self.class_labels)

    def row_cluster_slots(self) -> Tuple[int, ...]:
        index = self.cluster_palette_index()
        return tuple(index[label] for label in self._partition().labels)

    def ordered_sample_ids(self) -> Tuple[str, ...]:
        """Rows top to bottom in the drawn figure (columns left to right)."""
        return self.hierarchy.ordered_sample_ids()

    def cut_description(self) -> Optional[str]:
        """``"maxclust 4"`` or ``"height 0.5"``; ``None`` without a partition."""
        if self.assignment is None:
            return None
        if self.assignment.height is not None:
            return "height %g" % self.assignment.height
        return "maxclust %d" % self.assignment.n_clusters_requested

    def to_row_names_csv(self, heading: str = ROW_NAMES_HEADING) -> str:
        """The plotted order: one sample id per line, top row first, under ``heading``."""
        return "%s\n%s" % (
            heading,
            "".join("%s\n" % name for name in self.ordered_sample_ids()),
        )

    def to_assignments_csv(self) -> str:
        """``sample_id,cluster,leaf_position[,class_label]`` in input order."""
        return self._partition().to_csv(
            self.class_labels if self.has_classes else None)

    def to_grouped_csv(self) -> str:
        """Sample ids grouped under each cluster label."""
        out = []
        for label, names in self._partition().clusters().items():
            out.append("Cluster %d:\n" % label)
            out.extend("%s\n" % n for n in names)
            out.append("\n")
        return "".join(out)

    def provenance(self) -> Mapping[str, Any]:
        payload = dict(self.hierarchy.provenance())
        if self.assignment is not None:
            payload.update(self.assignment.provenance())
            payload["flat_partition"] = self.cut_description()
            payload.pop("flat_partition_note", None)
        payload.update(
            {
                "plotted_matrix": self.matrix.value_label,
                "plotted_scale": self.matrix.scale,
                "matrix": {k: v for k, v in self.matrix.provenance().items()
                           if k != "samples"},
                "row_cluster": self.row_cluster,
                "col_cluster": self.col_cluster,
                "row_linkage_is_col_linkage": True,
                "class_labels": list(self.class_labels) if self.has_classes else None,
                "class_labels_source": ("supplied by the caller"
                                        if self.has_classes else None),
                "class_palette": self.class_palette if self.has_classes else None,
                "class_palette_slots": dict(self.class_palette_index()),
                "cluster_palette": self.cluster_palette if self.has_clusters else None,
                "cluster_palette_slots": ({str(k): v for k, v in
                                           self.cluster_palette_index().items()}
                                          if self.has_clusters else None),
                "cmap": self.cmap,
                "figsize": list(self.figsize),
                "adjusted_rand_index": self.adjusted_rand,
                "visualization_affects_results": False,
            }
        )
        return payload


def clustermap_spec(
    matrix: SimilarityMatrix,
    class_labels: Optional[Sequence[Optional[str]]] = None,
    n_clusters: Optional[int] = None,
    *,
    height: Optional[float] = None,
    figsize: Tuple[float, float] = DEFAULT_FIGSIZE,
    cmap: str = "viridis",
    title: Optional[str] = None,
) -> ClustermapSpec:
    """Build the tree over ``matrix`` and describe the figure.

    With neither ``n_clusters`` nor ``height`` -- the default -- the figure is
    the tree, the dendrograms and the ordered heat map, and no partition is
    made.  ``n_clusters`` (``maxclust``) or ``height`` (cut merges above it)
    adds a flat partition and its stripe; at most one may be given.
    ``class_labels`` is optional: one label (or ``None``) per sample, in the
    matrix's sample order; with none the class stripe is left out.  The
    adjusted Rand index is computed only for a partition and when every sample
    has a class label.
    """
    if matrix.kind != DISTANCE:
        raise ValueError(
            "the clustermap plots a distance matrix; %r is a %s. Convert it "
            "explicitly with .to_distance()." % (matrix.measure, matrix.kind)
        )
    if class_labels is None:
        labels: Tuple[Optional[str], ...] = (None,) * matrix.n_samples
    else:
        if isinstance(class_labels, Mapping):
            raise TypeError(
                "class_labels must be a sequence in the matrix's sample order; "
                "for a mapping sample_id -> label pass "
                "[labels.get(s) for s in matrix.sample_ids]"
            )
        labels = tuple(None if c is None else str(c) for c in class_labels)
    if len(labels) != matrix.n_samples:
        raise ValueError(
            "%d class labels for %d samples" % (len(labels), matrix.n_samples)
        )
    if n_clusters is not None and height is not None:
        raise ValueError(
            "give one cutting rule, not two: n_clusters= (maxclust) or height= "
            "(cut every merge above it)"
        )
    tree = hierarchy(matrix)
    assignment = None
    ari = None
    if n_clusters is not None or height is not None:
        assignment = cut_hierarchy(tree, n_clusters=n_clusters, height=height)
        if all(c is not None for c in labels):
            ari = adjusted_rand_index(labels, assignment.labels)
    return ClustermapSpec(
        matrix=matrix,
        hierarchy=tree,
        class_labels=labels,
        assignment=assignment,
        figsize=tuple(figsize),
        cmap=cmap,
        title=title,
        adjusted_rand=ari,
    )
