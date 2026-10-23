"""Worker counts, and the one place that may import joblib.

``--jobs 1`` -- the default -- never comes here for execution: it is the plain
in-process loop, and needs nothing beyond the standard library.  ``--jobs N``
with N > 1 runs rows in N worker processes through joblib's process backend.
joblib is a standard dependency of BioTSR, installed with it (ADR-086); it is
imported here, inside a function, so that ``import biotsr`` and ``--jobs 1``
still work in an environment where it is missing.  Asking for parallel execution
there is an error with the installation command in it, never a silent fall-back
to one process.

How many CPUs this process may use
----------------------------------

On a workstation that is the machine's CPU count.  On a cluster it usually is
not: a job granted 8 cores on a 128-core node that starts 128 workers slows
itself down and its neighbours with it.  :func:`available_cpus` therefore reads,
the scheduler's allocation (Slurm, PBS/Torque, LSF, Grid Engine), the process's
CPU affinity mask, a cgroup CPU quota and the machine's own CPU count, and uses
the smallest of those that are present.  ``--jobs auto``
means exactly that number.  An explicit ``--jobs N`` is honoured as given, with
a warning when it exceeds it.

Nothing here changes a result: every row is computed by the same function,
with the same arguments, whichever process runs it.
"""

from __future__ import annotations

import math
import os
import re
from typing import Any, Callable, Iterable, Iterator, Optional, Tuple

from .errors import CLIError

__all__ = [
    "JOBS_AUTO",
    "available_cpus",
    "resolve_jobs",
    "parallel_map",
    "require_joblib",
]

JOBS_AUTO = "auto"

#: Scheduler variables that state how many CPUs this job was granted, in the
#: order they are consulted.  The first one set wins among these.
_SCHEDULER_VARIABLES = (
    ("SLURM_CPUS_PER_TASK", "Slurm"),
    ("SLURM_CPUS_ON_NODE", "Slurm"),
    ("PBS_NUM_PPN", "PBS/Torque"),
    ("NCPUS", "PBS"),
    ("LSB_DJOB_NUMPROC", "LSF"),
    ("NSLOTS", "Grid Engine"),
)

_LEADING_INT = re.compile(r"^\s*(\d+)")


def _scheduler_cpus() -> Optional[Tuple[int, str]]:
    for name, scheduler in _SCHEDULER_VARIABLES:
        match = _LEADING_INT.match(os.environ.get(name, ""))
        if match and int(match.group(1)) > 0:
            return int(match.group(1)), "%s (%s)" % (scheduler, name)
    return None


def _affinity_cpus() -> Optional[int]:
    getter = getattr(os, "sched_getaffinity", None)
    if getter is None:
        return None
    try:
        return len(getter(0)) or None
    except OSError:
        return None


def _cgroup_cpus() -> Optional[int]:
    """A cgroup v2 ``cpu.max`` or v1 ``cfs_quota_us`` limit, rounded up."""
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:
            quota, period = fh.read().split()[:2]
        if quota != "max":
            return max(1, int(math.ceil(float(quota) / float(period))))
        return None
    except (OSError, ValueError):
        pass
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fh:
            quota = float(fh.read().strip())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fh:
            period = float(fh.read().strip())
        if quota > 0 and period > 0:
            return max(1, int(math.ceil(quota / period)))
    except (OSError, ValueError):
        pass
    return None


def available_cpus() -> Tuple[int, str]:
    """``(count, where the count came from)`` for this process."""
    candidates = []
    scheduler = _scheduler_cpus()
    if scheduler:
        candidates.append(scheduler)
    affinity = _affinity_cpus()
    if affinity:
        candidates.append((affinity, "CPU affinity"))
    cgroup = _cgroup_cpus()
    if cgroup:
        candidates.append((cgroup, "cgroup CPU quota"))
    # The machine's own count bounds every other source: a scheduler variable
    # inherited from elsewhere (SLURM_CPUS_ON_NODE copied into an interactive
    # shell, say) must not start more workers than the machine has CPUs.
    machine = os.cpu_count()
    if machine:
        candidates.append((machine, "os.cpu_count()"))
    if candidates:
        return min(candidates, key=lambda c: c[0])
    return 1, "os.cpu_count() unavailable"


def resolve_jobs(value, warn: Optional[Callable[[str], None]] = None) -> int:
    """``"auto"``, ``"4"`` or ``4`` -> a worker count, or a :class:`CLIError`."""
    if value is None:
        return 1
    text = str(value).strip().lower()
    if text == JOBS_AUTO:
        return available_cpus()[0]
    try:
        jobs = int(text)
    except ValueError:
        raise CLIError(
            "--jobs %r is not a number" % value,
            "give a worker count such as --jobs 4, or --jobs auto for the CPUs "
            "this job was allocated",
        )
    if jobs < 1:
        raise CLIError(
            "--jobs must be at least 1 (got %d)" % jobs,
            "--jobs 1 runs the rows one after another in this process",
        )
    count, source = available_cpus()
    if jobs > count and warn is not None:
        warn(
            "--jobs %d asks for more workers than the %d CPU(s) available to "
            "this process (%s); they will compete for CPU time"
            % (jobs, count, source)
        )
    return jobs


def require_joblib():
    """The joblib module, or a :class:`CLIError` saying how to install it."""
    try:
        import joblib  # noqa: F401 -- imported here only (ADR-086)
    except ImportError:
        raise CLIError(
            "parallel execution (--jobs greater than 1) needs joblib, which "
            "BioTSR depends on but this environment cannot import",
            "reinstall BioTSR (pip install biotsr) or install joblib "
            "(pip install \"joblib>=1.3\") -- or run with --jobs 1, which "
            "needs nothing else",
        )
    return joblib


def parallel_map(function: Callable[[Any], Any], tasks: Iterable[Any],
                 jobs: int) -> Iterator[Any]:
    """``function(task)`` for every task in ``jobs`` processes, **in task order**.

    Results are yielded one at a time as soon as the next one in order is ready,
    and at most ``2 * jobs`` tasks are dispatched ahead of the consumer, so
    neither the task list nor the results pile up in memory.  Closing the
    iterator early (an error, Ctrl-C) stops dispatching and shuts the workers
    down.
    """
    joblib = require_joblib()
    runner = joblib.Parallel(
        n_jobs=jobs,
        backend="loky",
        return_as="generator",
        pre_dispatch="2*n_jobs",
        batch_size=1,
    )
    return runner(joblib.delayed(function)(task) for task in tasks)
