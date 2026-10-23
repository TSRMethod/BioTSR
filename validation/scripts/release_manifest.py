#!/usr/bin/env python3
"""The machine-readable release-candidate manifest.

One JSON document connecting the source tree to the built package and to the
evidence that validates it:

* **identity** -- BioTSR version, git commit and branch, whether the working
  tree was clean, and a content hash of ``src/biotsr`` that does not depend on
  git at all;
* **environment** -- Python, platform, and the versions of the (development-only)
  tools that built the artifacts; BioTSR's one runtime dependency is joblib;
* **science** -- the sha256 and content fingerprint of the packaged label table,
  the standard ``KeySpace`` fingerprint of every registered method, and the
  ``KeySpace``/``SamplingProfile`` fingerprints and counts of reference runs on
  the pinned 1UBQ structure;
* **validation** -- gate counts read from the committed gate evidence, and the
  pytest summary line when a log is supplied;
* **artifacts** -- file name, size and sha256 of every wheel and sdist in
  ``--dist-dir``.

Usage::

    python3 validation/scripts/release_manifest.py --dist-dir dist \\
        [--pytest-log LOG] [--output PATH]

Without ``--output`` the manifest is printed.  It is a *snapshot* of one machine
and one commit, not deterministic evidence, so CI writes it to ``dist/`` rather
than into the tree; the committed copy is
``validation/evidence/release/RELEASE_CANDIDATE_MANIFEST.json``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
sys.path.insert(0, SRC)

import biotsr  # noqa: E402
from biotsr.labels.elements import ELEMENT_LABELS_112  # noqa: E402
from biotsr.methods.registry import get_method, list_methods  # noqa: E402

EVIDENCE = os.path.join(WORKSPACE, "validation", "evidence")
FIXTURE = os.path.join(WORKSPACE, "tests", "data", "structures", "1UBQ.pdb")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git(*argv):
    try:
        return subprocess.run(["git", "-C", WORKSPACE] + list(argv), stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def source_tree_sha256():
    """sha256 over every file of ``src/biotsr`` (path + bytes), in path order."""
    h = hashlib.sha256()
    root = os.path.join(SRC, "biotsr")
    files = []
    for d, dirs, names in os.walk(root):
        dirs[:] = sorted(x for x in dirs if x != "__pycache__")
        files += [os.path.join(d, n) for n in names if not n.endswith((".pyc", ".pyo"))]
    for path in sorted(files, key=lambda p: os.path.relpath(p, SRC).replace(os.sep, "/")):
        rel = os.path.relpath(path, SRC).replace(os.sep, "/")
        h.update(rel.encode("utf-8") + b"\0")
        with open(path, "rb") as fh:
            h.update(fh.read())
        h.update(b"\0")
    return h.hexdigest(), len(files)


def tool_versions(names=("pip", "setuptools", "build", "twine", "wheel", "pytest")):
    out = {}
    for name in names:
        try:
            from importlib.metadata import version
            out[name] = version(name)
        except Exception:  # noqa: BLE001 -- absent is a fact, not an error
            out[name] = None
    return out


def lexicons():
    base = os.path.join(SRC, "biotsr", "labels", "resources")
    schemes = {
        "element_labels_v1.csv": ELEMENT_LABELS_112,
    }
    return {
        name: {
            "sha256": _sha256(os.path.join(base, name)),
            "scheme": scheme.name,
            "content_fingerprint": scheme.content_fingerprint(),
        }
        for name, scheme in sorted(schemes.items())
    }


def method_key_spaces():
    return {
        name: {
            "family": get_method(name).family,
            "standard_profile": get_method(name).default_profile.name,
            "key_space_fingerprint": get_method(name).default_profile.key_space().fingerprint,
        }
        for name in list_methods()
    }


def reference_runs():
    runs = {}
    for label, kwargs in (
        ("protein_1UBQ_A", dict(method="protein", chain="A")),
        ("protein_1UBQ_A_mirror_grouped", dict(method="protein", chain="A", mirror=True,
                                               aa_grouping=True)),
        ("amino_acid_atomic_1UBQ_A_LEU15", dict(method="amino_acid_atomic", chain="A",
                                                residue="15", residue_name="LEU")),
    ):
        r = biotsr.analyze(FIXTURE, keep_triplets=False, **kwargs)
        frequencies = r.frequencies          # a property that builds a new dict per access
        freq = "".join("%d\t%d\n" % (k, frequencies[k]) for k in sorted(frequencies))
        runs[label] = {
            "n_triangles": r.counts.n_triangles,
            "n_unique_keys": r.counts.n_unique_keys,
            "key_space_fingerprint": r.key_space.fingerprint,
            "sampling_fingerprint": r.sampling.fingerprint,
            "frequency_table_sha256": hashlib.sha256(freq.encode("ascii")).hexdigest(),
        }
    return runs


def _load(rel):
    path = os.path.join(WORKSPACE, rel)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


#: The committed gate evidence: what ``gates_cli.py`` and
#: ``gates_api_hardening.py`` write.
GATE_EVIDENCE = (
    ("cli", "validation/evidence/cli/GATE_RESULTS.json"),
    ("api_hardening", "validation/evidence/api_hardening/GATE_RESULTS.json"),
)


def gate_counts():
    out = {}
    for gates, rel in GATE_EVIDENCE:
        payload = _load(rel)
        if payload is None:
            out[gates] = None
            continue
        out[gates] = {"total": payload["n_gates"], "failed": payload["n_failed"],
                      "evidence": rel}
    return out


def _sdist_content_sha256(path):
    """sha256 over (member name, member sha256) of an sdist, in name order.

    The wheel is byte-reproducible under ``SOURCE_DATE_EPOCH``; an sdist's
    tar/gzip headers are not, so its *content* is what a rebuild can match.
    """
    import tarfile
    with tarfile.open(path) as tf:
        members = {m.name: hashlib.sha256(tf.extractfile(m).read()).hexdigest()
                   for m in tf.getmembers() if m.isfile()}
    listing = "".join("%s %s\n" % (k, members[k]) for k in sorted(members))
    return hashlib.sha256(listing.encode("utf-8")).hexdigest(), len(members)


def artifacts(dist):
    if not dist or not os.path.isdir(dist):
        return []
    out = []
    for f in sorted(os.listdir(dist)):
        if not f.endswith((".whl", ".tar.gz")):
            continue
        path = os.path.join(dist, f)
        entry = {"file": f, "bytes": os.path.getsize(path), "sha256": _sha256(path)}
        if f.endswith(".tar.gz"):
            entry["content_sha256"], entry["members"] = _sdist_content_sha256(path)
            entry["reproducibility"] = ("archive bytes vary between builds (tar/gzip "
                                        "headers); content_sha256 is reproducible")
        else:
            entry["reproducibility"] = ("byte-reproducible with SOURCE_DATE_EPOCH = "
                                        "the commit time")
        out.append(entry)
    return out


def pytest_summary(path):
    if not path or not os.path.exists(path):
        return None
    with open(path, errors="replace") as fh:
        lines = [l.strip() for l in fh if re.search(r"\d+ (passed|failed)", l)]
    return lines[-1].strip("= ") if lines else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dist-dir", default="")
    parser.add_argument("--pytest-log", default="")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    tree_hash, n_files = source_tree_sha256()
    status = _git("status", "--porcelain")
    manifest = {
        "manifest": "BioTSR release-candidate validation manifest",
        "generated_by": "validation/scripts/release_manifest.py",
        "identity": {
            "biotsr_version": biotsr.__version__,
            "result_schema_version": biotsr.RESULT_SCHEMA_VERSION,
            "git_commit": _git("rev-parse", "HEAD"),
            "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
            "git_tree_clean": (status == "") if status is not None else None,
            "source_tree_sha256": tree_hash,
            "source_tree_files": n_files,
        },
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "runtime_dependencies": tool_versions(("joblib",)),
            "development_tools": tool_versions(),
        },
        "science": {
            "lexicons": lexicons(),
            "standard_key_space_by_method": method_key_spaces(),
            "reference_runs": reference_runs(),
        },
        "validation": {
            "gates": gate_counts(),
            "pytest": pytest_summary(args.pytest_log),
        },
        "artifacts": artifacts(args.dist_dir),
    }
    text = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w") as fh:
            fh.write(text)
        print("-> %s" % args.output)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
