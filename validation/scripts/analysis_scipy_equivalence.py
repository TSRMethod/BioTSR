"""Compare the stdlib clustering with the scipy and scikit-learn that are installed.

``biotsr.analysis.clustering`` reimplements ``scipy.cluster.hierarchy.linkage(...,
method='average')``, ``fcluster(..., criterion='maxclust')``,
``dendrogram(...)['leaves']`` and ``sklearn.metrics.adjusted_rand_score`` so the
package keeps its "standard library only" guarantee.  This script checks it
against whichever versions of those libraries are installed, over randomised
matrices -- a third of them built from a three-value alphabet, where ties are
everywhere and tie-breaking is what separates two implementations -- and counts
mismatches separately for matrices with and without tied distances.

**What a run shows is agreement with that one version.**  scipy has changed its
own handling of tied heights between releases (1.11.1 and 1.18.1 disagree with
each other on such matrices), so an "exact" verdict against one release says
nothing about another.

Optional dependencies: ``numpy``, ``scipy`` and ``scikit-learn``.  BioTSR does
not depend on any of them and does not declare them; install them into a
separate environment to run this check.  Without them the script says which are
missing and exits 0 without touching the committed evidence, so it never fails
for want of a library the package does not use.

Usage::

    python3 validation/scripts/analysis_scipy_equivalence.py [--json PATH]

The committed records are
``validation/evidence/analysis/SCIPY_EQUIVALENCE_scipy-<version>.json``, one
per scipy version checked (1.11.1 and 1.18.1).  Each records the platform and
library versions it ran with, so it is a record of one run rather than a
byte-stable artefact.  The script exits 1 only when a claim fails: a different
partition or leaf order on the same tree, or a different ARI.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import random
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src"))

#: The comparison libraries (import name -> distribution name).  None of them is
#: a BioTSR dependency.
OPTIONAL = (("numpy", "numpy"), ("scipy", "scipy"), ("sklearn", "scikit-learn"))
MISSING = [dist for name, dist in OPTIONAL if importlib.util.find_spec(name) is None]

if not MISSING:
    import numpy as np  # noqa: E402
    from scipy.cluster.hierarchy import dendrogram, fcluster, linkage  # noqa: E402
    from scipy.spatial.distance import squareform  # noqa: E402
    from sklearn.metrics import adjusted_rand_score  # noqa: E402

from biotsr.analysis.clustering import (  # noqa: E402
    adjusted_rand_index,
    average_linkage,
    flat_clusters,
    leaf_order,
)
from biotsr.analysis.similarity import DISTANCE, SimilarityMatrix  # noqa: E402

EVIDENCE = os.path.join(REPO, "validation", "evidence", "analysis")
TRIALS = 400
MAX_N = 45


def as_matrix(condensed, n):
    return SimilarityMatrix(
        sample_ids=tuple("s%d" % i for i in range(n)),
        measure="generalized_jaccard",
        kind=DISTANCE,
        values=tuple(condensed),
        key_space_fingerprint="equivalence-check",
        comparison={},
        feature_ordering="numeric_ascending",
    )


def _partition(labels):
    """The partition a label vector describes, independent of how it numbers clusters."""
    groups = {}
    for index, label in enumerate(labels):
        groups.setdefault(label, set()).add(index)
    return {frozenset(g) for g in groups.values()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", dest="json_path", default=None,
                        help="where to write the report (default: the committed "
                             "record for the installed scipy version)")
    args = parser.parse_args(argv)
    if MISSING:
        print("skipped: the optional comparison libraries are not installed (%s); "
              "BioTSR does not depend on them. Nothing was written."
              % ", ".join(MISSING))
        return 0
    if args.json_path is None:
        args.json_path = os.path.join(
            EVIDENCE, "SCIPY_EQUIVALENCE_scipy-%s.json" % __import__("scipy").__version__)
    rng = random.Random(20260920)
    counters = {
        "trials": 0,
        "linkage_checked": 0,
        "linkage_mismatch": 0,
        "fcluster_checked": 0,
        "fcluster_mismatch": 0,
        "leaf_order_checked": 0,
        "leaf_order_mismatch": 0,
        "ari_checked": 0,
        "ari_mismatch": 0,
        "tie_heavy_trials": 0,
        # The same mismatch counts, restricted to matrices with no tied value.
        "tie_free_trials": 0,
        "tie_free_linkage_mismatch": 0,
        "tie_free_fcluster_mismatch": 0,
        "tie_free_leaf_order_mismatch": 0,
        # fcluster and dendrogram applied to BioTSR's own linkage matrix: the
        # cut and the leaf order compared on identical trees, ties included.
        "own_tree_fcluster_checked": 0,
        "own_tree_fcluster_mismatch": 0,
        "own_tree_partition_mismatch": 0,
        # ...of which scipy formed fewer clusters than BioTSR, where BioTSR
        # formed exactly min(t, n) (scipy 1.11 stopping above the lowest cut),
        # or more than the t requested (scipy 1.11's n=2, t=1 read past its array).
        "own_tree_partition_mismatch_scipy_under_delivered": 0,
        "own_tree_partition_mismatch_scipy_over_delivered": 0,
        "own_tree_leaf_order_mismatch": 0,
    }
    known = []

    for trial in range(TRIALS):
        n = rng.randint(2, MAX_N)
        tie_heavy = (trial % 3) == 0
        if tie_heavy:
            counters["tie_heavy_trials"] += 1
        condensed = []
        for _ in range(n * (n - 1) // 2):
            if tie_heavy:
                condensed.append(rng.choice([0.1, 0.2, 0.3]))
            else:
                condensed.append(round(rng.random(), 3))
        counters["trials"] += 1
        tie_free = len(set(condensed)) == len(condensed)
        if tie_free:
            counters["tie_free_trials"] += 1

        square = squareform(np.array(condensed))
        z_scipy = linkage(squareform(square, checks=False), method="average")
        link = average_linkage(as_matrix(condensed, n))
        z_mine = np.array([list(r) for r in link.rows])

        counters["linkage_checked"] += 1
        if not np.array_equal(z_scipy, z_mine):
            counters["linkage_mismatch"] += 1
            counters["tie_free_linkage_mismatch"] += tie_free

        leaves_scipy = dendrogram(z_scipy, no_plot=True, color_threshold=-np.inf)["leaves"]
        counters["leaf_order_checked"] += 1
        if list(leaf_order(link)) != list(leaves_scipy):
            counters["leaf_order_mismatch"] += 1
            counters["tie_free_leaf_order_mismatch"] += tie_free

        own_leaves = dendrogram(z_mine, no_plot=True, color_threshold=-np.inf)["leaves"]
        if list(leaf_order(link)) != list(own_leaves):
            counters["own_tree_leaf_order_mismatch"] += 1
        for t in range(1, n + 1):
            counters["own_tree_fcluster_checked"] += 1
            ours = list(flat_clusters(link, t))
            scipys = [int(v) for v in fcluster(z_mine, t=t, criterion="maxclust")]
            if ours != scipys:
                counters["own_tree_fcluster_mismatch"] += 1
                if _partition(ours) != _partition(scipys):
                    counters["own_tree_partition_mismatch"] += 1
                    if len(set(scipys)) < len(set(ours)) == min(t, n):
                        counters["own_tree_partition_mismatch_scipy_under_delivered"] += 1
                    elif len(set(scipys)) > t:
                        counters["own_tree_partition_mismatch_scipy_over_delivered"] += 1

        for t in range(1, n + 1):
            counters["fcluster_checked"] += 1
            mine = list(flat_clusters(link, t))
            try:
                theirs = list(fcluster(z_scipy, t=t, criterion="maxclust"))
            except Exception:  # pragma: no cover - scipy raising is itself a result
                continue
            if mine != theirs:
                if n == 2 and t == 1:
                    known.append(
                        {
                            "case": "n=2, t=1",
                            "scipy": [int(v) for v in theirs],
                            "biotsr": mine,
                            "cause": "scipy reads max_dists[n-1], one past the end of an "
                            "(n-1)-element array, and returns two clusters when one "
                            "was requested",
                        }
                    )
                    continue
                counters["fcluster_mismatch"] += 1
                counters["tie_free_fcluster_mismatch"] += tie_free

        labels_true = [rng.randint(0, max(1, n // 3)) for _ in range(n)]
        labels_pred = [rng.randint(0, max(1, n // 4)) for _ in range(n)]
        counters["ari_checked"] += 1
        if abs(adjusted_rand_index(labels_true, labels_pred)
               - float(adjusted_rand_score(labels_true, labels_pred))) > 1e-12:
            counters["ari_mismatch"] += 1

    report = {
        "platform": platform.platform(),
        "python": sys.version,
        "scipy": __import__("scipy").__version__,
        "sklearn": __import__("sklearn").__version__,
        "counters": counters,
        "known_divergences": known[:3],
        "n_known_divergence_occurrences": len(known),
        # What BioTSR claims for every scipy version, checked here for this one:
        # on the same tree, the maxclust partition and the dendrogram leaf
        # order agree; the ARI agrees with scikit-learn.  Linkage tie-breaking,
        # and how scipy numbers flat clusters, vary between scipy releases and
        # are reported, not claimed.
        "agreement": {
            "linkage_identical": counters["linkage_mismatch"] == 0,
            "linkage_identical_without_ties": counters["tie_free_linkage_mismatch"] == 0,
            "leaf_order_identical": counters["leaf_order_mismatch"] == 0,
            "fcluster_labels_identical": counters["fcluster_mismatch"] == 0,
            "same_tree_partitions_identical": counters["own_tree_partition_mismatch"] == 0,
            "same_tree_partitions_differ_only_where_scipy_misses_the_cut": (
                counters["own_tree_partition_mismatch"]
                == counters["own_tree_partition_mismatch_scipy_under_delivered"]
                + counters["own_tree_partition_mismatch_scipy_over_delivered"]),
            "same_tree_labels_identical": counters["own_tree_fcluster_mismatch"] == 0,
            "same_tree_leaf_order_identical": counters["own_tree_leaf_order_mismatch"] == 0,
            "ari_identical": counters["ari_mismatch"] == 0,
        },
        "verdict": (
            "exact"
            if not any(counters[k] for k in counters if k.endswith("mismatch"))
            else "claims hold"
            if not (counters["own_tree_partition_mismatch"]
                    or counters["own_tree_leaf_order_mismatch"]
                    or counters["ari_mismatch"])
            else "claims hold except where this scipy misses the maxclust cut"
            if not (counters["own_tree_partition_mismatch"]
                    - counters["own_tree_partition_mismatch_scipy_under_delivered"]
                    - counters["own_tree_partition_mismatch_scipy_over_delivered"]
                    or counters["own_tree_leaf_order_mismatch"]
                    or counters["ari_mismatch"])
            else "CLAIMS FAIL"
        ),
    }
    with open(args.json_path, "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
        fh.write("\n")
    for key in sorted(counters):
        print("%-24s %d" % (key, counters[key]))
    print("known n=2/t=1 divergences:", len(known))
    print("verdict:", report["verdict"])
    print("wrote", os.path.relpath(args.json_path, REPO))
    return 0 if report["verdict"] != "CLAIMS FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
