#!/usr/bin/env python3
"""Measure what SSE annotation costs, against the same run without it.

The question is narrow and answerable: for one structure, one
chain, one profile, how much slower and how much larger is a run with
``sse=True`` than the identical run with ``sse=False``?  Everything else is held
fixed, and the *outputs are compared* as part of the measurement -- a speed
number for a run that produced different keys would be meaningless.

Each configuration is run ``REPEATS`` times and the **minimum** is reported, not
the mean: the minimum is the measurement least contaminated by whatever else the
machine was doing, which is the standard practice ``timeit`` documents.

Nothing here is allowed to change scientific behaviour, and nothing here is
optimized in response to what it finds; it reports.

Usage::

    python3 validation/scripts/performance_sse.py [--json PATH]

The timings are a record of one machine, not evidence a rebuild can reproduce,
so nothing is written into the tree: the report is printed, and saved only where
``--json`` asks.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import platform
import sys
import time
import tracemalloc

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "src"))

import biotsr  # noqa: E402
from biotsr.sse import check_conservation  # noqa: E402

FIXTURES = os.path.join(REPO, "tests", "data", "structures")

REPEATS = 3

#: (label, structure, chain, modifier kwargs)
CONFIGURATIONS = [
    ("1UBQ_A regular", "1UBQ", "A", {}),
    ("1UBQ_A grouping", "1UBQ", "A", {"aa_grouping": True}),
    ("1UBQ_A mirror", "1UBQ", "A", {"mirror": True}),
    ("1UBQ_A size_filter_25", "1UBQ", "A", {"size_filter": 25.0}),
    ("1UBQ_A all three", "1UBQ", "A",
     {"aa_grouping": True, "mirror": True, "size_filter": 25.0}),
    ("1BRS_A regular", "1BRS", "A", {}),
    ("1TUP_A regular", "1TUP", "A", {}),
]


def _time_run(path, chain, sse, kwargs, keep_triplets):
    best = None
    result = None
    for _ in range(REPEATS):
        gc.collect()
        start = time.perf_counter()
        result = biotsr.analyze(
            path, method="protein", chain=chain, sse=sse,
            keep_triplets=keep_triplets, **kwargs
        )
        elapsed = time.perf_counter() - start
        best = elapsed if best is None else min(best, elapsed)
    return best, result


def _peak_memory(path, chain, sse, kwargs, keep_triplets):
    gc.collect()
    tracemalloc.start()
    biotsr.analyze(
        path, method="protein", chain=chain, sse=sse,
        keep_triplets=keep_triplets, **kwargs
    )
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", dest="json_path", default="",
                        help="also save the report here (nothing is written by default)")
    args = parser.parse_args(argv)
    rows = []
    for label, structure, chain, kwargs in CONFIGURATIONS:
        path = os.path.join(FIXTURES, "%s.pdb" % structure)
        if not os.path.exists(path):
            print("skipping %s: fixture absent" % label)
            continue

        entry = {
            "label": label,
            "structure": structure,
            "chain": chain,
            "modifiers": {k: v for k, v in kwargs.items()},
        }

        for mode, keep in (("with_triplets", True), ("keys_only", False)):
            t_plain, plain = _time_run(path, chain, False, kwargs, keep)
            t_sse, annotated = _time_run(path, chain, True, kwargs, keep)

            # The measurement is only meaningful if the two runs agree.
            assert annotated.frequencies == plain.frequencies, label
            assert annotated.key_space.fingerprint == plain.key_space.fingerprint
            assert annotated.sampling.fingerprint == plain.sampling.fingerprint
            check_conservation(annotated)

            entry.setdefault("counts", {
                "n_points": plain.selection.n_points,
                "n_candidates": plain.counts.n_candidates,
                "n_triangles": plain.n_triangles,
                "n_unique_keys": plain.n_unique_keys,
                "n_filtered_by_size": plain.counts.n_filtered_by_size,
                "n_sse_categories_observed": sum(
                    1 for n in annotated.category_totals().values() if n
                ),
            })
            entry[mode] = {
                "baseline_seconds": round(t_plain, 6),
                "sse_seconds": round(t_sse, 6),
                "overhead_ratio": round(t_sse / t_plain, 4),
                "overhead_percent": round(100 * (t_sse / t_plain - 1), 2),
                "microseconds_per_candidate_baseline": round(
                    1e6 * t_plain / max(plain.counts.n_candidates, 1), 4
                ),
                "microseconds_per_candidate_sse": round(
                    1e6 * t_sse / max(plain.counts.n_candidates, 1), 4
                ),
            }

        peak_plain = _peak_memory(path, chain, False, kwargs, True)
        peak_sse = _peak_memory(path, chain, True, kwargs, True)
        entry["memory_with_triplets"] = {
            "baseline_peak_bytes": peak_plain,
            "sse_peak_bytes": peak_sse,
            "ratio": round(peak_sse / peak_plain, 4) if peak_plain else None,
            "extra_bytes_per_triangle": (
                round((peak_sse - peak_plain) / entry["counts"]["n_triangles"], 2)
                if entry["counts"]["n_triangles"]
                else None
            ),
        }
        peak_plain_keys = _peak_memory(path, chain, False, kwargs, False)
        peak_sse_keys = _peak_memory(path, chain, True, kwargs, False)
        entry["memory_keys_only"] = {
            "baseline_peak_bytes": peak_plain_keys,
            "sse_peak_bytes": peak_sse_keys,
            "ratio": round(peak_sse_keys / peak_plain_keys, 4) if peak_plain_keys else None,
        }

        rows.append(entry)
        print(
            "%-26s cand=%-9d tri=%-9d  runtime x%.3f (triplets) / x%.3f (keys only)"
            " mem x%.3f"
            % (
                label,
                entry["counts"]["n_candidates"],
                entry["counts"]["n_triangles"],
                entry["with_triplets"]["overhead_ratio"],
                entry["keys_only"]["overhead_ratio"],
                entry["memory_with_triplets"]["ratio"],
            )
        )

    payload = {
        "generated_by": "validation/scripts/performance_sse.py",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "repeats_per_measurement": REPEATS,
        "statistic": "minimum of %d runs" % REPEATS,
        "note": (
            "Every pair is asserted to produce identical keys, identical "
            "fingerprints and conserved SSE frequencies before its timings are "
            "recorded. No configuration was optimized in response to these numbers."
        ),
        "runs": rows,
    }
    if args.json_path:
        os.makedirs(os.path.dirname(os.path.abspath(args.json_path)), exist_ok=True)
        with open(args.json_path, "w") as fh:
            json.dump(payload, fh, indent=2)
            fh.write("\n")
        print("\nwrote %s" % args.json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
