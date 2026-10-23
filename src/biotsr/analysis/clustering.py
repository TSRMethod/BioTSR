"""Hierarchical clustering: the computation, with no plotting attached.

Every evidenced clustering path in the family is the same three calls, and the
three scripts that make them (``clusters/clustermap/clustermap.py``, its
``_chain_hexamer`` sibling and ``Kinases-and-Phosphatases-Clustering/clustermap_n.py``)
differ only in how they choose the cluster count and what they print:

    linkage = scipy.cluster.hierarchy.linkage(squareform(D), method='average')
    labels  = scipy.cluster.hierarchy.fcluster(linkage, t, criterion='maxclust')
    order   = seaborn's dendrogram_row.reordered_ind

so: **average linkage (UPGMA) on a precomputed distance matrix**, flat clusters
by ``maxclust``, and scipy's default dendrogram leaf order.  Here the tree and
its leaf order (:func:`hierarchy`) are the default result and need no cluster
count; a flat partition is optional and needs an explicit rule
(:func:`cut_hierarchy`: ``maxclust`` or a cut height).  No Ward, no
Euclidean, no scaling beyond the ``100 *`` percent conversion, no random
initialisation and therefore no seed -- the whole path is deterministic
(ADR-060; docs/methods/ANALYSIS.md).

This module reimplements those three in the standard library so key generation's
"no third-party imports" guarantee extends to the analysis layer.  It follows
the algorithms of scipy's ``nn_chain``, ``cluster_maxclust_monocrit`` and
``dendrogram``: the mergesort-stable ordering of ``Z`` by height, the
union-find relabelling, and the Lance-Williams update written in scipy's operand
order.  **Scope of the agreement:** ``validation/scripts/analysis_scipy_equivalence.py``
compares it with *one installed scipy version* over random matrices with heavy
ties, and the result is recorded with that version in
``validation/evidence/analysis/SCIPY_EQUIVALENCE_scipy-<version>.json``.  scipy has changed its
own handling of tied heights between releases, so agreement with every scipy
version is not claimed; with tied distances the partition, the leaf order and
the linkage rows can differ between scipy versions, and between this module and
a scipy version it was not checked against.

**Ties.**  Equal distances are resolved deterministically (the lowest index
wins in the nearest-neighbour search, and equal heights keep their merge order),
so the same matrix always gives the same linkage, partition and leaf order.
When merges tie at the cut, ``maxclust`` returns fewer clusters than requested;
the assignment records both numbers and a warning.

**Scale.**  The linkage is computed on the matrix exactly as given and records
its ``scale``: fractions (``1``) by default, percent (``100``) only if the
caller passed ``matrix.scaled(100)``.  Multiplying every distance by a constant
changes no partition or order in exact arithmetic; in floating point, two
distances one rounding unit apart can become equal and resolve as a tie, so
nothing is rescaled here behind the caller's back.

``maxclust`` follows its definition -- the lowest threshold that leaves at most
*t* clusters -- as later scipy releases do.  scipy 1.11.1 did not always find
that threshold (it could return fewer clusters than requested with distinct
merge heights, and two clusters for ``n == 2, t == 1``); see
:func:`flat_clusters` and ADR-060.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence, Tuple

from .similarity import DISTANCE, SimilarityMatrix

__all__ = [
    "AVERAGE",
    "MAXCLUST",
    "DISTANCE_CUT",
    "Linkage",
    "Hierarchy",
    "ClusterAssignment",
    "average_linkage",
    "hierarchy",
    "cut_hierarchy",
    "flat_clusters",
    "flat_clusters_at_height",
    "leaf_order",
    "cluster",
    "adjusted_rand_index",
]

#: The only linkage method any authoritative source uses.
AVERAGE = "average"
#: Flat-cut rules.  ``maxclust``: the lowest cut leaving at most *t* clusters.
#: ``distance``: cut every merge above a stated height.
MAXCLUST = "maxclust"
DISTANCE_CUT = "distance"


def _condensed_index(n: int, i: int, j: int) -> int:
    if i > j:
        i, j = j, i
    return n * i - (i * (i + 1)) // 2 + (j - i - 1)


@dataclass(frozen=True)
class Linkage:
    """A linkage matrix, in scipy's ``(n-1) x 4`` layout.

    Each row is ``(left, right, height, size)`` with cluster ids ``0..n-1`` for
    the original samples and ``n + k`` for the cluster formed by row ``k``.
    """

    rows: Tuple[Tuple[float, float, float, float], ...]
    sample_ids: Tuple[str, ...]
    method: str
    metric: str
    scale: float

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def heights(self) -> Tuple[float, ...]:
        return tuple(r[2] for r in self.rows)

    @property
    def is_monotonic(self) -> bool:
        h = self.heights
        return all(h[i] <= h[i + 1] for i in range(len(h) - 1))

    def provenance(self) -> Mapping[str, Any]:
        return {
            "algorithm": "agglomerative hierarchical clustering (nn-chain)",
            "linkage_method": self.method,
            "metric": self.metric,
            "scale": self.scale,
            "n_samples": self.n_samples,
            "monotonic": self.is_monotonic,
            "deterministic": True,
            "random_seed": None,
        }


def average_linkage(matrix: SimilarityMatrix) -> Linkage:
    """UPGMA over a precomputed distance matrix, on its values as given.

    ``matrix`` must be a distance; a similarity is refused rather than silently
    inverted, because ``1 - J`` is a modelling decision the caller makes with
    ``to_distance()`` (ADR-057).  Heights are in the matrix's own units
    (``matrix.scale``: 1 for fractions, 100 for percent).  Time is quadratic
    and the working copy of the condensed matrix is held in memory.
    """
    if matrix.kind != DISTANCE:
        raise ValueError(
            "average_linkage needs a distance matrix; %r is a %s. Convert it "
            "explicitly with .to_distance()." % (matrix.measure, matrix.kind)
        )
    n = matrix.n_samples
    if n < 2:
        raise ValueError("hierarchical clustering needs at least 2 samples, got %d" % n)
    rows = _nn_chain_average(list(matrix.condensed()), n)
    return Linkage(
        rows=tuple(tuple(r) for r in rows),
        sample_ids=matrix.sample_ids,
        method=AVERAGE,
        metric="precomputed:%s:%s" % (matrix.measure, matrix.kind),
        scale=matrix.scale,
    )


def _nn_chain_average(condensed: Sequence[float], n: int):
    """scipy's ``nn_chain`` specialised to average linkage.

    The Lance-Williams update is written ``(nx*d_xi + ny*d_yi) / (nx + ny)``
    because that is the operand order of scipy's ``_average``; any other
    grouping is a different sequence of IEEE roundings.
    """
    D = list(condensed)
    Z = [[0.0, 0.0, 0.0, 0.0] for _ in range(n - 1)]
    size = [1] * n
    chain = [0] * n
    chain_length = 0
    for k in range(n - 1):
        if chain_length == 0:
            chain_length = 1
            for i in range(n):
                if size[i] > 0:
                    chain[0] = i
                    break
        while True:
            x = chain[chain_length - 1]
            if chain_length > 1:
                y = chain[chain_length - 2]
                current_min = D[_condensed_index(n, x, y)]
            else:
                y = -1
                current_min = math.inf
            for i in range(n):
                if size[i] == 0 or x == i:
                    continue
                d = D[_condensed_index(n, x, i)]
                if d < current_min:
                    current_min = d
                    y = i
            if chain_length > 1 and y == chain[chain_length - 2]:
                break
            chain[chain_length] = y
            chain_length += 1
        chain_length -= 2
        if x > y:
            x, y = y, x
        nx, ny = size[x], size[y]
        Z[k][0] = float(x)
        Z[k][1] = float(y)
        Z[k][2] = current_min
        Z[k][3] = float(nx + ny)
        size[x] = 0
        size[y] = nx + ny
        for i in range(n):
            if size[i] == 0 or i == y:
                continue
            d_xi = D[_condensed_index(n, i, x)]
            d_yi = D[_condensed_index(n, i, y)]
            D[_condensed_index(n, i, y)] = (nx * d_xi + ny * d_yi) / (nx + ny)
    order = sorted(range(n - 1), key=lambda i: Z[i][2])  # stable, as mergesort
    Z = [Z[i] for i in order]
    _relabel(Z, n)
    return Z


def _relabel(Z, n: int) -> None:
    """scipy's ``label``: union-find renumbering after the height sort."""
    parent = list(range(2 * n - 1))
    sizes = [1.0] * (2 * n - 1)
    next_label = n

    def find(x: int) -> int:
        p = x
        while parent[x] != x:
            x = parent[x]
        while parent[p] != x:
            parent[p], p = x, parent[p]
        return x

    for i in range(n - 1):
        x_root = find(int(Z[i][0]))
        y_root = find(int(Z[i][1]))
        if x_root < y_root:
            Z[i][0], Z[i][1] = float(x_root), float(y_root)
        else:
            Z[i][0], Z[i][1] = float(y_root), float(x_root)
        parent[x_root] = next_label
        parent[y_root] = next_label
        total = sizes[x_root] + sizes[y_root]
        sizes[next_label] = total
        next_label += 1
        Z[i][3] = total


def _max_height_per_node(linkage: Linkage) -> Tuple[float, ...]:
    """scipy's ``get_max_dist_for_each_cluster``: the tallest merge in a subtree."""
    n = linkage.n_samples
    md = [0.0] * (n - 1)
    for i, row in enumerate(linkage.rows):
        highest = row[2]
        for child in (int(row[0]), int(row[1])):
            if child >= n:
                value = md[child - n]
                if value > highest:
                    highest = value
        md[i] = highest
    return tuple(md)


def _monocrit_labels(linkage: Linkage, mc: Sequence[float], cutoff: float):
    """scipy's ``cluster_monocrit``: flat labels, numbered in its traversal order."""
    n = linkage.n_samples
    Z = linkage.rows
    labels = [0] * n
    stack = [0] * n
    visited = bytearray(2 * n - 1)
    k = 0
    stack[0] = 2 * n - 2
    n_cluster = 0
    leader = -1
    while k >= 0:
        root = stack[k] - n
        i_lc = int(Z[root][0])
        i_rc = int(Z[root][1])
        if leader == -1 and mc[root] <= cutoff:
            leader = root
            n_cluster += 1
        if i_lc >= n and not visited[i_lc]:
            visited[i_lc] = 1
            k += 1
            stack[k] = i_lc
            continue
        if i_rc >= n and not visited[i_rc]:
            visited[i_rc] = 1
            k += 1
            stack[k] = i_rc
            continue
        if i_lc < n:
            if leader == -1:
                n_cluster += 1
            labels[i_lc] = n_cluster
        if i_rc < n:
            if leader == -1:
                n_cluster += 1
            labels[i_rc] = n_cluster
        if leader == root:
            leader = -1
        k -= 1
    return labels


def _count_clusters(linkage: Linkage, mc: Sequence[float], threshold: float, max_nc: int) -> int:
    """scipy's counting traversal inside ``cluster_maxclust_monocrit``."""
    n = linkage.n_samples
    Z = linkage.rows
    visited = bytearray(2 * n - 1)
    stack = [0] * n
    nc = 0
    k = 0
    stack[0] = 2 * n - 2
    while k >= 0:
        root = stack[k] - n
        i_lc = int(Z[root][0])
        i_rc = int(Z[root][1])
        if mc[root] <= threshold:
            nc += 1
            if nc > max_nc:
                break
            k -= 1
            visited[i_lc] = 1
            visited[i_rc] = 1
            continue
        if not visited[i_lc]:
            visited[i_lc] = 1
            if i_lc >= n:
                k += 1
                stack[k] = i_lc
                continue
            nc += 1
            if nc > max_nc:
                break
        if not visited[i_rc]:
            visited[i_rc] = 1
            if i_rc >= n:
                k += 1
                stack[k] = i_rc
                continue
            nc += 1
            if nc > max_nc:
                break
        k -= 1
    return nc


def flat_clusters(linkage: Linkage, n_clusters: int) -> Tuple[int, ...]:
    """Flat clusters by the ``maxclust`` criterion, 1-based labels.

    The definition: cut the tree at the **lowest** threshold that leaves at most
    ``n_clusters`` clusters, i.e. the smallest ``r`` among "below every merge"
    and the merge heights for which no flat cluster holds a merge above ``r``.
    When merge heights are distinct this gives exactly
    ``min(n_clusters, n)`` clusters; when the merges around the cut tie, every
    merge at that height goes together and fewer clusters form (the
    :func:`cluster` result then carries a warning).

    Labels are numbered in scipy's ``cluster_monocrit`` traversal order.  scipy
    1.11.1 searched only part of the candidate thresholds and could return
    fewer clusters than an exact cut allows even with distinct heights (and
    two clusters for ``n == 2, t == 1``); later scipy releases (1.18.1 checked)
    and this function both follow the definition above (ADR-060, amended).
    """
    if isinstance(n_clusters, bool) or not isinstance(n_clusters, int):
        raise TypeError("n_clusters must be an integer, got %r" % (n_clusters,))
    if n_clusters < 1:
        raise ValueError("n_clusters must be at least 1, got %d" % n_clusters)
    n = linkage.n_samples
    mc = _max_height_per_node(linkage)
    if n_clusters >= n:
        # Below every merge: each sample is its own cluster.
        return tuple(_monocrit_labels(linkage, mc, -math.inf))
    # The cluster count is non-increasing in the threshold, and the root's
    # height always leaves one cluster: find the lowest admissible height.
    candidates = sorted(set(mc))
    lower, upper = -1, len(candidates) - 1
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if _count_clusters(linkage, mc, candidates[middle], n_clusters) > n_clusters:
            lower = middle
        else:
            upper = middle
    return tuple(_monocrit_labels(linkage, mc, candidates[upper]))


def flat_clusters_at_height(linkage: Linkage, height: float) -> Tuple[int, ...]:
    """Flat clusters by the ``distance`` criterion, 1-based labels.

    Every merge above ``height`` is cut, so no flat cluster contains two
    samples joined higher than ``height`` (scipy's ``fcluster(...,
    criterion="distance")``).  ``height`` is in the matrix's own units.  Labels
    are numbered in the same traversal order as :func:`flat_clusters`.
    """
    if isinstance(height, bool) or not isinstance(height, (int, float)):
        raise TypeError("height must be a number, got %r" % (height,))
    if not height >= 0 or height == math.inf:
        raise ValueError("height must be a finite number of at least 0, got %r"
                         % (height,))
    return tuple(_monocrit_labels(linkage, _max_height_per_node(linkage), float(height)))


def leaf_order(linkage: Linkage) -> Tuple[int, ...]:
    """scipy's default dendrogram leaf order, which seaborn reports as ``reordered_ind``.

    Depth-first, left child before right child, with neither ``count_sort`` nor
    ``distance_sort`` -- seaborn calls ``dendrogram(..., no_plot=True)`` and takes
    ``['leaves']``, so the defaults are what the historical figures used.
    """
    n = linkage.n_samples
    leaves = []
    stack = [2 * n - 2]
    while stack:
        node = stack.pop()
        if node < n:
            leaves.append(node)
            continue
        row = linkage.rows[node - n]
        stack.append(int(row[1]))
        stack.append(int(row[0]))
    return tuple(leaves)


@dataclass(frozen=True)
class Hierarchy:
    """The average-linkage tree over a distance matrix, and its leaf order.

    This is the whole of the default workflow -- distance matrix, tree,
    dendrogram and ordered heat map -- and it needs no number of clusters.
    An uncut tree describes nested groupings at every height; it does not
    determine one final number of groups, so none is reported here.  A flat
    partition is a separate, optional step with an explicit rule
    (:func:`cut_hierarchy`).
    """

    linkage: Linkage
    leaves: Tuple[int, ...]
    sample_ids: Tuple[str, ...]
    measure: str
    key_space_fingerprint: str
    comparison: Mapping[str, Any]
    feature_ordering: str
    warnings: Tuple[str, ...] = ()
    samples: Tuple[Mapping[str, Any], ...] = ()

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def ordered_sample_ids(self) -> Tuple[str, ...]:
        return tuple(self.sample_ids[i] for i in self.leaves)

    def leaf_position(self) -> Tuple[int, ...]:
        """Per sample, in input order: its row in the dendrogram order (0 = top)."""
        position = [0] * len(self.leaves)
        for rank, index in enumerate(self.leaves):
            position[index] = rank
        return tuple(position)

    def to_linkage_csv(self) -> str:
        """The merges, one per line: ``step,left,right,height,size``.

        ``left``/``right`` are sample ids for a leaf and ``cluster<k>`` for the
        cluster formed at step *k* (1-based), so the tree can be read without
        knowing scipy's numbering.  ``height`` is in the matrix's units.
        """
        import csv
        import io

        n = self.n_samples
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["step", "left", "right", "height", "size"])

        def name(node):
            node = int(node)
            return self.sample_ids[node] if node < n else "cluster%d" % (node - n + 1)

        for step, (left, right, height, size) in enumerate(self.linkage.rows, 1):
            writer.writerow([step, name(left), name(right), repr(height), int(size)])
        return out.getvalue()

    def to_order_csv(self) -> str:
        """``sample_id,leaf_position``, top row first (the plotted order)."""
        return "sample_id,leaf_position\n" + "".join(
            "%s,%d\n" % (name, rank)
            for rank, name in enumerate(self.ordered_sample_ids()))

    def provenance(self) -> Mapping[str, Any]:
        payload = dict(self.linkage.provenance())
        payload.update(
            {
                "measure": self.measure,
                "feature_ordering": self.feature_ordering,
                "key_space_fingerprint": self.key_space_fingerprint,
                "comparison": dict(self.comparison),
                "preprocessing": "none: the matrix is clustered as given "
                                 "(scale %g)" % self.linkage.scale,
                "sample_order": list(self.sample_ids),
                "leaf_order": [self.sample_ids[i] for i in self.leaves],
                "flat_partition": None,
                "flat_partition_note": (
                    "none: an uncut hierarchy does not determine a number of "
                    "groups; a flat partition needs an explicit cutting rule "
                    "(cut_hierarchy)"
                ),
                "warnings": list(self.warnings),
                "samples": [dict(p) for p in self.samples],
            }
        )
        return payload


def hierarchy(matrix: SimilarityMatrix) -> Hierarchy:
    """Average linkage and the dendrogram leaf order -- no cluster count needed."""
    link = average_linkage(matrix)
    return Hierarchy(
        linkage=link,
        leaves=leaf_order(link),
        sample_ids=matrix.sample_ids,
        measure=matrix.measure,
        key_space_fingerprint=matrix.key_space_fingerprint,
        comparison=matrix.comparison,
        feature_ordering=matrix.feature_ordering,
        warnings=tuple(matrix.warnings),
        samples=matrix.samples,
    )


@dataclass(frozen=True)
class ClusterAssignment:
    """A flat partition cut from a hierarchy by an explicit rule.

    ``criterion`` is :data:`MAXCLUST` (``n_clusters_requested`` set) or
    :data:`DISTANCE_CUT` (``height`` set).  Exists independently of any figure.
    """

    linkage: Linkage
    labels: Tuple[int, ...]
    n_clusters_requested: Optional[int]
    leaves: Tuple[int, ...]
    sample_ids: Tuple[str, ...]
    measure: str
    key_space_fingerprint: str
    comparison: Mapping[str, Any]
    feature_ordering: str
    warnings: Tuple[str, ...] = ()
    samples: Tuple[Mapping[str, Any], ...] = ()
    criterion: str = MAXCLUST
    height: Optional[float] = None

    @property
    def n_clusters(self) -> int:
        return len(set(self.labels))

    def ordered_sample_ids(self) -> Tuple[str, ...]:
        return tuple(self.sample_ids[i] for i in self.leaves)

    def leaf_position(self) -> Tuple[int, ...]:
        """Per sample, in input order: its row in the dendrogram order (0 = top)."""
        position = [0] * len(self.leaves)
        for rank, index in enumerate(self.leaves):
            position[index] = rank
        return tuple(position)

    def to_csv(self, class_labels: Optional[Sequence[Optional[str]]] = None) -> str:
        """``sample_id,cluster,leaf_position[,class_label]``, one row per sample, input order.

        ``cluster`` is the 1-based flat-cluster label; ``leaf_position`` the row
        the sample occupies in the dendrogram (and in the clustermap), from 0 at
        the top.  ``class_label`` is written only when labels are given, and an
        unlabelled sample is written as an empty cell.
        """
        import csv
        import io

        if class_labels is not None and len(class_labels) != len(self.sample_ids):
            raise ValueError("%d class labels for %d samples"
                             % (len(class_labels), len(self.sample_ids)))
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        header = ["sample_id", "cluster", "leaf_position"]
        if class_labels is not None:
            header.append("class_label")
        writer.writerow(header)
        position = self.leaf_position()
        for i, name in enumerate(self.sample_ids):
            row = [name, self.labels[i], position[i]]
            if class_labels is not None:
                row.append("" if class_labels[i] is None else class_labels[i])
            writer.writerow(row)
        return out.getvalue()

    def clusters(self) -> Mapping[int, Tuple[str, ...]]:
        out: dict = {}
        for name, label in zip(self.sample_ids, self.labels):
            out.setdefault(label, []).append(name)
        return {k: tuple(v) for k, v in sorted(out.items())}

    def provenance(self) -> Mapping[str, Any]:
        payload = dict(self.linkage.provenance())
        payload.update(
            {
                "criterion": self.criterion,
                "n_clusters_requested": self.n_clusters_requested,
                "height": self.height,
                "n_clusters_formed": self.n_clusters,
                "measure": self.measure,
                "feature_ordering": self.feature_ordering,
                "key_space_fingerprint": self.key_space_fingerprint,
                "comparison": dict(self.comparison),
                "cut_rule_choice": (
                    "chosen by the caller; BioTSR never infers a number of "
                    "clusters or a cut height, and never takes it from known "
                    "class labels"
                ),
                "preprocessing": "none: the matrix is clustered as given "
                                 "(scale %g)" % self.linkage.scale,
                "sample_order": list(self.sample_ids),
                "leaf_order": [self.sample_ids[i] for i in self.leaves],
                "cluster_labels": list(self.labels),
                "warnings": list(self.warnings),
                "samples": [dict(p) for p in self.samples],
            }
        )
        return payload


def cut_hierarchy(tree: Hierarchy, *, n_clusters: Optional[int] = None,
                  height: Optional[float] = None) -> ClusterAssignment:
    """A flat partition of ``tree``, by exactly one explicit rule.

    * ``n_clusters=t`` -- ``maxclust``: the lowest cut leaving at most *t*
      clusters (:func:`flat_clusters`).  When merges tie at the cut, fewer form
      and the result carries a warning;
    * ``height=h`` -- ``distance``: cut every merge above *h*, in the matrix's
      units (:func:`flat_clusters_at_height`).

    The rule is the caller's; it is never inferred, and choosing it from known
    class labels would make any agreement with those labels a consequence of
    the choice.
    """
    if (n_clusters is None) == (height is None):
        raise ValueError(
            "a flat partition needs exactly one explicit cutting rule: "
            "n_clusters=<t> (maxclust) or height=<h> (cut merges above h). "
            "The tree and the ordered heat map need neither."
        )
    warnings = tuple(tree.warnings)
    if n_clusters is not None:
        if isinstance(n_clusters, bool) or not isinstance(n_clusters, int):
            raise TypeError("n_clusters must be an integer, got %r" % (n_clusters,))
        labels = flat_clusters(tree.linkage, n_clusters)
        formed = len(set(labels))
        if formed != min(n_clusters, tree.n_samples):
            warnings += (
                "%d cluster(s) were requested and %d formed: merges tie at the "
                "cut, so no threshold leaves exactly %d"
                % (n_clusters, formed, n_clusters),
            )
        criterion = MAXCLUST
    else:
        labels = flat_clusters_at_height(tree.linkage, height)
        criterion = DISTANCE_CUT
    return ClusterAssignment(
        linkage=tree.linkage,
        labels=labels,
        n_clusters_requested=n_clusters,
        leaves=tree.leaves,
        sample_ids=tree.sample_ids,
        measure=tree.measure,
        key_space_fingerprint=tree.key_space_fingerprint,
        comparison=tree.comparison,
        feature_ordering=tree.feature_ordering,
        warnings=warnings,
        samples=tree.samples,
        criterion=criterion,
        height=None if height is None else float(height),
    )


def cluster(matrix: SimilarityMatrix, n_clusters: int) -> ClusterAssignment:
    """Average linkage cut by ``maxclust`` into at most ``n_clusters``, in one call.

    Shorthand for ``cut_hierarchy(hierarchy(matrix), n_clusters=n_clusters)``.
    ``n_clusters`` is the caller's choice and is recorded as such.  When
    distances tie, fewer clusters than requested can be formed; the result then
    carries a warning and ``n_clusters`` (formed) differs from
    ``n_clusters_requested``.
    """
    return cut_hierarchy(hierarchy(matrix), n_clusters=n_clusters)


def adjusted_rand_index(labels_true: Sequence[Any], labels_pred: Sequence[Any]) -> float:
    """Agreement between a clustering and known classes, as the legacy scripts report it.

    ``clustermap.py`` prints ``sklearn.metrics.adjusted_rand_score(residue_types,
    cluster_labels)``, so the ARI is part of the evidenced numerical output and
    not a presentation detail.  This follows sklearn's own route -- the pair
    confusion matrix, with its two exact special cases -- rather than the
    textbook contingency formula, because the two disagree on degenerate inputs
    and sklearn's answer is the one in the frozen references.
    """
    if len(labels_true) != len(labels_pred):
        raise ValueError(
            "label vectors have different lengths: %d and %d"
            % (len(labels_true), len(labels_pred))
        )
    n = len(labels_true)
    if n == 0:
        return 1.0
    n_classes = len(set(labels_true))
    n_clusters = len(set(labels_pred))
    if n_classes == n_clusters == 1 or n_classes == n_clusters == n:
        return 1.0
    contingency: dict = {}
    a_counts: dict = {}
    b_counts: dict = {}
    for t, p in zip(labels_true, labels_pred):
        contingency[(t, p)] = contingency.get((t, p), 0) + 1
        a_counts[t] = a_counts.get(t, 0) + 1
        b_counts[p] = b_counts.get(p, 0) + 1
    sum_sq = sum(v * v for v in contingency.values())
    sum_a_sq = sum(v * v for v in a_counts.values())
    sum_b_sq = sum(v * v for v in b_counts.values())
    # sklearn's pair_confusion_matrix, scaled by 2 throughout as it is there.
    tp = sum_sq - n
    fp = sum_b_sq - sum_sq
    fn = sum_a_sq - sum_sq
    tn = n * n - sum_a_sq - sum_b_sq + sum_sq
    if fn == 0 and fp == 0:
        return 1.0
    num = 2.0 * (tp * tn - fn * fp)
    den = (tp + fn) * (fn + tn) + (tp + fp) * (fp + tn)
    return num / den
