#!/usr/bin/env python3
"""How much does the command line cost?

The question is not "is the CLI fast" -- it computes nothing, so of
course it is.  It is whether the orchestration layer is *negligible* next to the
science, because if it is not, a batch of 5,000 samples pays for it 5,000 times.

Three numbers are measured separately, because they behave differently:

``interpreter``
    ``python -c "import biotsr"`` -- what every invocation pays before any argument
    is read.  Paid once per process, so it matters to a batch only if the batch
    is one process per sample (it is not; ``biotsr batch`` is one process).
``orchestration``
    Everything the CLI does around one calculation *within* a process: parse the
    arguments, resolve the method and profile, validate, name the output
    prefix, stage the files, install them and write the ``PREFIX.run.json``
    completion record.  Measured by running the same
    calculation through ``execute()`` and through the full ``main()`` and taking
    the difference, with serialization held constant.
``science``
    The ``analyze()``/``analyze_cross()`` call itself.

Nothing here optimises anything.  The command line does not touch the engine.

Usage::

    python3 validation/scripts/performance_cli.py [--repeats N] [--json PATH]
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
FIXTURES = os.path.join(WORKSPACE, "tests", "data", "structures")
EVIDENCE = os.path.join(WORKSPACE, "validation", "evidence", "cli")

sys.path.insert(0, SRC)

import biotsr  # noqa: E402
from biotsr.cli import main as cli_main  # noqa: E402
from biotsr.io import native  # noqa: E402

#: (tag, CLI flags, direct call). Deliberately spans three orders of magnitude of
#: computation, because the whole question is a ratio.
CASES = [
    (
        "drug_3PTB_BEN",
        ["run", "-m", "drug", "--pdb", os.path.join(FIXTURES, "3PTB.pdb"),
         "-c", "A", "--entity-name", "BEN", "--entity-id", "1", "--no-triplets"],
        lambda: biotsr.analyze(
            os.path.join(FIXTURES, "3PTB.pdb"), method="drug", chain="A",
            entity_name="BEN", entity_id="1", keep_triplets=False,
        ),
    ),
    (
        "protein_1UBQ_A",
        ["run", "-m", "protein", "--pdb", os.path.join(FIXTURES, "1UBQ.pdb"),
         "-c", "A", "--no-triplets"],
        lambda: biotsr.analyze(
            os.path.join(FIXTURES, "1UBQ.pdb"), method="protein", chain="A",
            keep_triplets=False,
        ),
    ),
    (
        "protein_1TUP_A",
        ["run", "-m", "protein", "--pdb", os.path.join(FIXTURES, "1TUP.pdb"),
         "-c", "A", "--no-triplets"],
        lambda: biotsr.analyze(
            os.path.join(FIXTURES, "1TUP.pdb"), method="protein", chain="A",
            keep_triplets=False,
        ),
    ),
    (
        "protein_1TUP_A_sse",
        ["run", "-m", "protein", "--pdb", os.path.join(FIXTURES, "1TUP.pdb"),
         "-c", "A", "--sse", "--no-triplets"],
        lambda: biotsr.analyze(
            os.path.join(FIXTURES, "1TUP.pdb"), method="protein", chain="A",
            sse=True, keep_triplets=False,
        ),
    ),
]


def _timings(fn, repeats):
    """Every run's elapsed time, in order."""
    out = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        out.append(time.perf_counter() - start)
    return out


def _best(fn, repeats):
    """The fastest of ``repeats`` runs.

    The minimum rather than the mean: the thing being measured is how long the
    work takes, and every sample is the work plus a non-negative amount of
    interference from the rest of the machine.
    """
    return min(_timings(fn, repeats))


def _interpreter_startup(repeats):
    """``import biotsr`` in a fresh process, and a bare interpreter for comparison."""
    env = dict(os.environ, PYTHONPATH=SRC)

    def run(code):
        subprocess.check_call(
            [sys.executable, "-c", code], env=env, stdout=subprocess.DEVNULL
        )

    bare = _best(lambda: run("pass"), repeats)
    imported = _best(lambda: run("import biotsr"), repeats)
    cli = _best(lambda: run("from biotsr.cli import main"), repeats)
    return {
        "bare_interpreter_s": bare,
        "import_tsr_s": imported,
        "import_tsr_cli_s": cli,
        "package_import_cost_s": imported - bare,
        "cli_import_cost_s": cli - imported,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--json", dest="json_path", default="")
    args = parser.parse_args()

    import tempfile

    print("interpreter startup (%d repeats)" % args.repeats)
    startup = _interpreter_startup(args.repeats)
    for key, value in sorted(startup.items()):
        print("  %-26s %8.4f s" % (key, value))
    print()

    print("%-22s %10s %10s %10s %9s %9s" % (
        "case", "science", "cli total", "overhead", "overhead%", "noise"))
    rows = []
    for tag, argv, direct in CASES:
        science_runs = _timings(direct, args.repeats)
        science = min(science_runs)
        # The spread of the science measurement is the noise floor of this case.
        # An orchestration figure smaller than it -- including a negative one --
        # means the CLI's cost is below what this method can resolve, which is
        # the answer to the question, not a failure of it.
        noise = max(science_runs) - science

        def one_cli():
            import contextlib
            import io as _io

            out = tempfile.mkdtemp(prefix="tsr-perf-")
            try:
                # The record path on stdout is the command's result, not noise --
                # it is swallowed here only so the table stays readable.
                with contextlib.redirect_stdout(_io.StringIO()):
                    code = cli_main(
                        argv + ["-o", out, "-q", "--if-exists", "overwrite"]
                    )
                assert code == 0, code
            finally:
                import shutil

                shutil.rmtree(out, ignore_errors=True)

        total = _best(one_cli, args.repeats)

        # Serialization, measured on its own so it is not charged to the CLI:
        # `biotsr run` writes files that a bare `analyze()` call does not, and that
        # cost belongs to `biotsr.io`, which the API user pays too.
        result = direct()

        def serialize():
            out = tempfile.mkdtemp(prefix="tsr-perf-io-")
            try:
                native.save(result, out, triplets=False)
            finally:
                import shutil

                shutil.rmtree(out, ignore_errors=True)

        io_cost = _best(serialize, args.repeats)
        overhead = total - science - io_cost
        rows.append(
            {
                "case": tag,
                "science_s": science,
                "science_runs_s": science_runs,
                "science_noise_s": noise,
                "serialization_s": io_cost,
                "cli_total_s": total,
                "orchestration_overhead_s": overhead,
                "orchestration_overhead_percent": 100.0 * overhead / total,
                "below_noise_floor": overhead < noise,
            }
        )
        print(
            "%-22s %10.4f %10.4f %10.4f %8.2f%%  %9.4f %s"
            % (tag, science, total, overhead, 100.0 * overhead / total, noise,
               "(below noise)" if overhead < noise else "")
        )

    worst = max(r["orchestration_overhead_s"] for r in rows)
    worst_share = max(r["orchestration_overhead_percent"] for r in rows)
    resolvable = [r for r in rows if not r["below_noise_floor"]]
    print()
    print("worst orchestration overhead : %.4f s per run" % worst)
    print("worst overhead share         : %.2f%% of one run" % worst_share)
    print(
        "cases below the noise floor  : %d of %d -- the CLI's cost is smaller "
        "than the run-to-run spread of the science it wraps"
        % (len(rows) - len(resolvable), len(rows))
    )
    print(
        "process startup (import biotsr) : %.4f s, paid once per process, not per "
        "sample" % startup["package_import_cost_s"]
    )

    payload = {
        "python": sys.version.split()[0],
        "repeats": args.repeats,
        "startup": startup,
        "cases": rows,
        "worst_orchestration_overhead_s": worst,
        "worst_orchestration_overhead_percent": worst_share,
        "n_cases_below_noise_floor": sum(1 for r in rows if r["below_noise_floor"]),
    }
    # Timings belong to the machine they were taken on: nothing is written into
    # the repository unless --json names a file.
    if not args.json_path:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    path = os.path.abspath(args.json_path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("-> %s" % os.path.relpath(path, WORKSPACE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
