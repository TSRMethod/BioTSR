#!/usr/bin/env python3
"""API hardening: what does saving, loading and re-analysing actually cost?

The one comparison that matters practically:

    recompute the structure + TSR + analysis
        versus
    load the saved TSR result + analysis

Nothing here optimises anything, and nothing here is a promise about other
machines or other structures.  These are measurements on pinned fixtures on one
laptop, reported so that a later claim has something to stand on.  In particular:
**a ratio measured on 1UBQ is not a speedup the package offers in general** --
the recompute side scales with C(n,3) in the number of selected points and the
load side scales with the number of distinct keys, so the ratio is a property of
the structure, not of the code.  Three structures of very different size are
measured for exactly that reason.

Memory is measured with ``tracemalloc``, which counts Python allocations rather
than RSS.  That is the right instrument for the question asked here -- "does the
result hold a pathological duplicate of something?" -- and the wrong one for
"how much RAM does a cluster job need", so the numbers are labelled as what they
are.

Usage::

    python3 validation/scripts/performance_api.py [--repeats N] [--json PATH]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
import time
import tracemalloc

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
FIXTURES = os.path.join(WORKSPACE, "tests", "data", "structures")
EVIDENCE = os.path.join(WORKSPACE, "validation", "evidence", "api_hardening")

sys.path.insert(0, SRC)

import biotsr  # noqa: E402
from biotsr.analysis import (  # noqa: E402
    AnalysisSample,
    common_key_intersection,
    distance_matrix,
    feature_matrix,
)
from biotsr.io import native, native_read  # noqa: E402

#: Small, medium and large by selected-point count, which is what the candidate
#: enumeration scales in.
#:
#: Each case is paired with the **next** one, cyclically, for the two-sample
#: comparison below.  Pairing every case with a fixed partner would have paired
#: the first one with itself, and a feature matrix of a sample against its own
#: copy is 0% sparse by construction -- a number that says nothing about the
#: representation and reads as though it did.
CASES = [
    ("1UBQ", "A"),
    ("1BRS", "A"),
    ("3PTB", "A"),
]


def timed(fn, repeats):
    """Best of ``repeats``.  Best, not mean: the distribution is one-sided.

    A slow run measures the machine (a scheduler decision, a page fault, another
    process); it does not measure the code.  The fastest run is the closest
    available estimate of what the code costs, and reporting a mean would mix in
    noise that is not a property of anything being compared.
    """
    best = None
    value = None
    for _ in range(repeats):
        start = time.perf_counter()
        value = fn()
        elapsed = time.perf_counter() - start
        best = elapsed if best is None else min(best, elapsed)
    return best, value


def peak_mb(fn):
    """Peak Python allocation during ``fn``, in MB, plus its return value."""
    tracemalloc.start()
    try:
        value = fn()
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak / (1024.0 * 1024.0), value


def fixture(pdb_id):
    path = os.path.join(FIXTURES, "%s.pdb" % pdb_id)
    if not os.path.exists(path):
        raise SystemExit("pinned fixture missing: %s" % path)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--json", dest="json_path", default=None)
    args = parser.parse_args()

    tmp = tempfile.mkdtemp(prefix="tsr-api-perf-")
    rows = []

    try:
        for index, (pdb_id, chain) in enumerate(CASES):
            partner_id, partner_chain = CASES[(index + 1) % len(CASES)]
            partner_path = fixture(partner_id)
            path = fixture(pdb_id)
            out = os.path.join(tmp, pdb_id)

            # --- compute -----------------------------------------------------
            compute_s, result = timed(
                lambda: biotsr.analyze(path, method="protein", chain=chain,
                                    keep_triplets=False),
                args.repeats,
            )
            compute_mb, _ = peak_mb(
                lambda: biotsr.analyze(path, method="protein", chain=chain,
                                    keep_triplets=False)
            )

            # --- serialize ---------------------------------------------------
            save_s, written = timed(
                lambda: native.save(result, out, triplets=False), args.repeats
            )
            save_mb, _ = peak_mb(
                lambda: native.save(result, out + "_mem", triplets=False)
            )

            # --- deserialize -------------------------------------------------
            load_s, loaded = timed(
                lambda: native_read.load_result(written["keys"]), args.repeats
            )
            load_mb, _ = peak_mb(lambda: native_read.load_result(written["keys"]))

            # --- the comparison that matters ---------------------------------
            # One analysis over two samples, reached two ways.  The analysis
            # itself is identical, so the difference is entirely in how the
            # samples were obtained.
            def recompute_then_analyse():
                a = biotsr.analyze(path, method="protein", chain=chain,
                                keep_triplets=False)
                b = biotsr.analyze(partner_path, method="protein",
                                chain=partner_chain, keep_triplets=False)
                samples = [AnalysisSample.from_result(a, sample_id="a"),
                           AnalysisSample.from_result(b, sample_id="b")]
                return common_key_intersection(samples).n_common_keys

            other = biotsr.analyze(partner_path, method="protein",
                                chain=partner_chain, keep_triplets=False)
            other_written = native.save(
                other, os.path.join(tmp, "partner_%s" % partner_id), triplets=False
            )

            def load_then_analyse():
                a = native_read.load_result(written["keys"])
                b = native_read.load_result(other_written["keys"])
                samples = [AnalysisSample.from_result(a, sample_id="a"),
                           AnalysisSample.from_result(b, sample_id="b")]
                return common_key_intersection(samples).n_common_keys

            recompute_s, n_recompute = timed(recompute_then_analyse, args.repeats)
            reload_s, n_reload = timed(load_then_analyse, args.repeats)
            assert n_recompute == n_reload, (
                "the two paths disagreed: %d vs %d" % (n_recompute, n_reload)
            )

            recompute_mb, _ = peak_mb(recompute_then_analyse)
            reload_mb, _ = peak_mb(load_then_analyse)

            # --- feature matrix, dense vs sparse ------------------------------
            samples = [
                AnalysisSample.from_result(loaded, sample_id="a"),
                AnalysisSample.from_result(
                    native_read.load_result(other_written["keys"]), sample_id="b"
                ),
            ]
            fm_mb, fm = peak_mb(lambda: feature_matrix(samples))
            dense_cells = len(fm.sample_ids) * len(fm.columns)
            stored_cells = sum(len(row) for row in fm.rows)
            dm_mb, _ = peak_mb(lambda: distance_matrix(fm))

            rows.append(
                {
                    "structure": pdb_id,
                    "chain": chain,
                    "partner": "%s_%s" % (partner_id, partner_chain),
                    "n_points": result.counts.n_points,
                    "n_candidates": result.counts.n_candidates,
                    "n_triangles": result.n_triangles,
                    "n_unique_keys": result.n_unique_keys,
                    "keys_file_bytes": os.path.getsize(written["keys"]),
                    "sidecar_bytes": os.path.getsize(written["metadata"]),
                    "seconds": {
                        "compute": compute_s,
                        "serialize": save_s,
                        "deserialize": load_s,
                        "recompute_then_analyse": recompute_s,
                        "load_then_analyse": reload_s,
                        "ratio_recompute_over_load": (
                            recompute_s / reload_s if reload_s else None
                        ),
                    },
                    "peak_python_mb": {
                        "compute": compute_mb,
                        "serialize": save_mb,
                        "deserialize": load_mb,
                        "recompute_then_analyse": recompute_mb,
                        "load_then_analyse": reload_mb,
                        "feature_matrix": fm_mb,
                        "distance_matrix": dm_mb,
                    },
                    "feature_matrix": {
                        "dense_cells": dense_cells,
                        "stored_cells": stored_cells,
                        "sparsity": (
                            1.0 - stored_cells / float(dense_cells)
                            if dense_cells else 0.0
                        ),
                    },
                }
            )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # --- report -------------------------------------------------------------
    print("%-8s %8s %10s %10s  %9s %9s %9s  %7s"
          % ("PDB", "points", "triangles", "keys",
             "compute", "save", "load", "ratio"))
    for row in rows:
        s = row["seconds"]
        print("%-8s %8d %10d %10d  %8.3fs %8.4fs %8.4fs  %6.1fx"
              % (row["structure"], row["n_points"], row["n_triangles"],
                 row["n_unique_keys"], s["compute"], s["serialize"],
                 s["deserialize"], s["ratio_recompute_over_load"]))

    print()
    print("peak Python allocation (tracemalloc, not RSS), MB")
    print("%-8s %9s %9s %9s  %11s %11s"
          % ("PDB", "compute", "save", "load", "recompute+an", "load+an"))
    for row in rows:
        m = row["peak_python_mb"]
        print("%-8s %9.1f %9.1f %9.1f  %11.1f %11.1f"
              % (row["structure"], m["compute"], m["serialize"], m["deserialize"],
                 m["recompute_then_analyse"], m["load_then_analyse"]))

    print()
    print("feature matrix storage (this case paired with the next, cyclically)")
    for row in rows:
        f = row["feature_matrix"]
        print("  %-8s + %-8s %7d dense cells, %7d stored (%.1f%% sparse)"
              % (row["structure"], row["partner"], f["dense_cells"],
                 f["stored_cells"], 100.0 * f["sparsity"]))

    print()
    print("The ratio is a property of the structure, not of the package: the")
    print("recompute side scales with C(n,3) in selected points and the load side")
    print("with the number of distinct keys. It is not a general speedup claim.")

    payload = {
        "repeats": args.repeats,
        "instrument": "time.perf_counter (best of N); tracemalloc peak, not RSS",
        "caveat": (
            "Measured on pinned fixtures on one machine. The recompute/load "
            "ratio is a property of the structure's size, not a general claim."
        ),
        "rows": rows,
    }
    # Timings belong to the machine they were taken on: nothing is written into
    # the repository unless --json names a file.
    if not args.json_path:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    path = os.path.abspath(args.json_path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("-> %s" % os.path.relpath(path, WORKSPACE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
