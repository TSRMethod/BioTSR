#!/usr/bin/env python3
"""Clean-install acceptance: install the package the way a user would, then use it.

Nothing else in this repository has ever installed ``biotsr``.  The test suite puts
``src/`` on ``sys.path`` and the validation scripts set ``PYTHONPATH``,
so every green run so far has proved only that the *source tree* works.  That is
exactly the arrangement that hides a missing entry point, an undeclared
dependency, package data that is not packaged, and a module that quietly depends
on being run from the checkout.

So this script:

1. builds a fresh virtual environment with no access to the source tree;
2. installs the package from the repository root through the ordinary packaging path;
3. runs the ``biotsr`` console script and ``python -m biotsr`` **from a directory
   outside the repository**, with the repository's ``src`` nowhere on the path;
4. checks the outputs against the same numbers a direct API call produces.

Usage::

    python3 validation/scripts/clean_install_check.py [--keep] [--json PATH] [--artifact FILE]

``--artifact`` installs a built wheel or sdist instead of the source directory
(``release_artifacts.py`` runs every check below against each release
artifact).  The JSON records no absolute path: the temporary
environment is reported as ``<clean-env>``.

Exit status is 0 when every check passes and 1 otherwise, so it can be a gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
PACKAGE = WORKSPACE  # the package project is the repository root
FIXTURES = os.path.join(WORKSPACE, "tests", "data", "structures")
LEXICONS = os.path.join(WORKSPACE, "src", "biotsr", "labels", "resources")


def source_version():
    """The version the source tree declares -- the one an install must report."""
    with open(os.path.join(WORKSPACE, "src", "biotsr", "__init__.py"), encoding="utf-8") as fh:
        return re.search(r'^__version__ = "([^"]+)"', fh.read(), re.M).group(1)

CHECKS = []


def check(name, category="install"):
    def decorate(fn):
        CHECKS.append((name, category, fn))
        return fn

    return decorate


class Env(object):
    """A virtual environment with the repository nowhere on its path."""

    def __init__(self, root, target=PACKAGE):
        self.root = root
        #: What ``pip install`` is given: the source directory, or a built artifact.
        self.target = target
        self.venv = os.path.join(root, "venv")
        self.cwd = os.path.join(root, "elsewhere")
        self.bootstrap_ok = False
        os.makedirs(self.cwd)

    @property
    def python(self):
        return os.path.join(self.venv, "bin", "python")

    @property
    def biotsr(self):
        return os.path.join(self.venv, "bin", "biotsr")

    def environ(self):
        """A clean environment: no PYTHONPATH, no inherited site directories."""
        env = dict(os.environ)
        env.pop("PYTHONPATH", None)
        env["PYTHONNOUSERSITE"] = "1"
        return env

    def run(self, argv, cwd=None):
        return subprocess.run(
            argv,
            cwd=cwd or self.cwd,
            env=self.environ(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

    def create(self):
        subprocess.check_call(
            [sys.executable, "-m", "venv", self.venv],
            stdout=subprocess.DEVNULL,
        )

    def _pip(self, *argv):
        return subprocess.run(
            [self.python, "-m", "pip"] + list(argv),
            env=self.environ(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

    def install(self):
        """Install the repository root the way a user would: ``pip install <dir>``.

        The venv module seeds an old ``setuptools`` and no ``wheel`` on some
        interpreters, and the PEP 517 build then fails with ``invalid command
        'bdist_wheel'`` before it ever reaches this package.  That is a property
        of the bootstrap, not of ``biotsr``, so the build backend is brought up to
        date first and the result is recorded -- the point of this script is to
        test the package's packaging, not pip's.
        """
        bootstrap = self._pip("install", "--upgrade", "pip", "setuptools", "wheel")
        self.bootstrap_ok = bootstrap.returncode == 0
        proc = self._pip("install", self.target)
        if proc.returncode != 0 and self.bootstrap_ok:
            proc = self._pip("install", "--no-build-isolation", self.target)
        return proc


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


@check("the console script exists on PATH")
def _console_script(env, state):
    return os.path.isfile(env.biotsr) and os.access(env.biotsr, os.X_OK), env.biotsr


@check("biotsr --help runs from outside the repository")
def _help(env, state):
    proc = env.run([env.biotsr, "--help"])
    return (
        proc.returncode == 0 and "COMMAND" in proc.stdout,
        "exit %d" % proc.returncode,
    )


@check("biotsr --version reports the packaged version")
def _version(env, state):
    proc = env.run([env.biotsr, "--version"])
    expected = "biotsr %s" % source_version()
    return proc.returncode == 0 and proc.stdout.strip() == expected, proc.stdout.strip()


@check("python -m biotsr is equivalent to the console script")
def _dash_m(env, state):
    a = env.run([env.biotsr, "methods"])
    b = env.run([env.python, "-m", "biotsr", "methods"])
    return (
        a.returncode == 0 and a.stdout == b.stdout,
        "%d bytes each" % len(a.stdout),
    )


@check("the source tree is not on the installed interpreter's path", "isolation")
def _not_from_source(env, state):
    proc = env.run(
        [env.python, "-c", "import biotsr, sys; sys.stdout.write(biotsr.__file__)"]
    )
    path = proc.stdout.strip()
    inside = os.path.abspath(path).startswith(os.path.abspath(PACKAGE))
    return (not inside) and "site-packages" in path, path


@check("packaged data (the element label table) is present and loadable", "packaging")
def _package_data(env, state):
    # The label table ships inside the wheel and is read through
    # importlib.resources. Three legacy repositories read their table by
    # relative path and therefore only work after a `cd`; this is the check
    # that says the package does not.
    proc = env.run(
        [env.python, "-c",
         "from biotsr import ELEMENT_LABELS_112 as L; "
         "print(L.n_labels, L.n_distinct_codes, len(L.table))"]
    )
    ok = proc.returncode == 0 and proc.stdout.split() == ["112", "109", "110"]
    return ok, ("radix %s, %s labels, %s elements" % tuple(proc.stdout.split())) if ok \
        else proc.stderr.strip()[:160]


@check("the installed label table is byte-identical to the source tree's", "packaging")
def _package_data_hashes(env, state):
    names = sorted(f for f in os.listdir(LEXICONS) if f.endswith(".csv"))
    code = (
        "import hashlib, json, importlib.resources as r\n"
        "base = r.files('biotsr.labels') / 'resources'\n"
        "print(json.dumps({n: hashlib.sha256((base / n).read_bytes()).hexdigest()"
        " for n in %r}))" % (names,)
    )
    proc = env.run([env.python, "-c", code])
    if proc.returncode != 0:
        return False, proc.stderr.strip()[-160:]
    installed = json.loads(proc.stdout)
    source = {}
    for n in names:
        with open(os.path.join(LEXICONS, n), "rb") as fh:
            source[n] = hashlib.sha256(fh.read()).hexdigest()
    return installed == source, "%d/%d sha256 match" % (
        sum(installed.get(n) == source[n] for n in names), len(names))


#: What the standard-only package no longer ships.  An install that still
#: carries any of these was built from a stale tree or a stale build directory.
REMOVED_MODULES = (
    "biotsr.io.legacy_text",
    "biotsr.io.legacy_atomic",
    "biotsr.io.legacy_cross",
    "biotsr.io.legacy_sse",
    "biotsr.profiles.legacy",
    "biotsr.profiles.catalogue",
    "biotsr.labels.lexicon",
)
REMOVED_FILES = (
    "io/legacy_text.py",
    "io/legacy_atomic.py",
    "io/legacy_cross.py",
    "io/legacy_sse.py",
    "profiles/legacy.py",
    "profiles/catalogue.py",
    "labels/resources/atom_lexicon_v1.csv",
    "labels/resources/atom_lexicon_v2.csv",
    "labels/resources/atom_lexicon_standard_v1.csv",
    "labels/lexicon.py",
)


@check("the installed package exposes no removed legacy module", "packaging")
def _no_legacy_modules(env, state):
    code = (
        "import importlib, json\n"
        "out = {}\n"
        "for name in %r:\n"
        "    try:\n"
        "        importlib.import_module(name)\n"
        "        out[name] = 'imported'\n"
        "    except ModuleNotFoundError:\n"
        "        out[name] = 'absent'\n"
        "print(json.dumps(out))" % (REMOVED_MODULES,)
    )
    proc = env.run([env.python, "-c", code])
    if proc.returncode != 0:
        return False, proc.stderr.strip()[-160:]
    found = sorted(n for n, v in json.loads(proc.stdout).items() if v != "absent")
    return not found, ("importable: %s" % ", ".join(found)) if found else (
        "%d removed modules, none importable" % len(REMOVED_MODULES))


@check("the installed package ships no removed legacy file or table", "packaging")
def _no_legacy_files(env, state):
    proc = env.run([env.python, "-c",
                    "import os, biotsr; print(os.path.dirname(biotsr.__file__))"])
    if proc.returncode != 0:
        return False, proc.stderr.strip()[-160:]
    root = proc.stdout.strip()
    present = sorted(rel for rel in REMOVED_FILES
                     if os.path.exists(os.path.join(root, *rel.split("/"))))
    tables = sorted(
        f for f in os.listdir(os.path.join(root, "labels", "resources"))
        if f.endswith(".csv"))
    ok = not present and tables == ["element_labels_v1.csv"]
    return ok, ("present: %s" % ", ".join(present)) if present else (
        "label tables: %s" % ", ".join(tables))


@check("joblib imports in the clean install", "packaging")
def _joblib_imports(env, state):
    proc = env.run([env.python, "-c", "import joblib; print(joblib.__version__)"])
    return proc.returncode == 0, (proc.stdout.strip() if proc.returncode == 0
                                  else proc.stderr.strip()[-160:])


@check("an ordinary install brings joblib, BioTSR's one dependency", "packaging")
def _dependencies(env, state):
    """ADR-086: joblib is the one runtime dependency, installed without asking.

    What BioTSR itself requires is checked, not the whole environment: joblib
    brings its own requirements, and the bootstrap's build tools bring theirs.
    """
    proc = env.run([env.python, "-m", "pip", "show", "biotsr"])
    requires = next((line.split(":", 1)[1] for line in proc.stdout.splitlines()
                     if line.startswith("Requires:")), "")
    declared = {r.strip().lower() for r in requires.split(",") if r.strip()}
    imported = env.run([env.python, "-c", "import joblib; print(joblib.__version__)"])
    return (declared == {"joblib"} and imported.returncode == 0,
            "requires %s; joblib %s" % (", ".join(sorted(declared)) or "nothing",
                                        imported.stdout.strip() or "not importable"))


#: The 1.0 methods (README, "Method families"): exactly these, nothing deferred.
RELEASED_METHODS = frozenset((
    "aa_ind", "aa_ind_mix", "aa_ind_specific", "aa_selection",
    "amino_acid_atom_selection", "amino_acid_atomic", "carbohydrate",
    "carbohydrate_atom_selection", "carbohydrate_chain", "carbohydrate_protein_cross",
    "drug", "drug_atom_selection", "drug_protein_cross", "nucleotide",
    "nucleotide_atom_selection", "nucleotide_chain", "nucleotide_protein_cross",
    "protein", "protein_protein_calpha_cross", "protein_protein_cross",
))


@check("biotsr methods lists exactly the released methods", "discovery")
def _methods(env, state):
    proc = env.run([env.biotsr, "methods", "--json"])
    records = json.loads(proc.stdout)
    names = {r["name"] for r in records}
    return names == RELEASED_METHODS, "%d methods" % len(names)


@check("a single-structure run works outside the repository", "workflow")
def _single(env, state):
    out = os.path.join(env.cwd, "single")
    proc = env.run(
        [env.biotsr, "run", "-m", "protein", "--pdb", os.path.join(FIXTURES, "1UBQ.pdb"),
         "-c", "A", "-o", out, "--no-triplets"]
    )
    if proc.returncode != 0:
        return False, proc.stderr.strip()[:200]
    # One shared folder: the calculation's files and its completion record
    # sit directly in --output under one prefix.
    payload = json.load(open(os.path.join(out, "1UBQ_A_protein.run.json")))
    state["single"] = payload
    counts = payload["results"][0]["counts"]
    return (
        payload["complete"] is True
        and counts["n_triangles"] == 70300 and counts["n_unique_keys"] == 45349,
        "%d triangles, %d keys" % (counts["n_triangles"], counts["n_unique_keys"]),
    )


@check("the installed CLI agrees with the in-repo API", "equivalence")
def _equivalence(env, state):
    """The numbers from the installed console script, against the source tree.

    Two different interpreters, two different copies of the package, one answer.
    """
    sys.path.insert(0, os.path.join(PACKAGE, "src"))
    import biotsr as source_tsr

    reference = source_tsr.analyze(
        os.path.join(FIXTURES, "1UBQ.pdb"), method="protein", chain="A",
        keep_triplets=False,
    )
    payload = state.get("single")
    if payload is None:
        return False, "the single-structure run did not produce a completion record"
    entry = payload["results"][0]
    return (
        entry["key_space_fingerprint"] == reference.key_space.fingerprint
        and entry["sampling_fingerprint"] == reference.sampling.fingerprint
        and entry["counts"]["n_triangles"] == reference.counts.n_triangles
        and entry["counts"]["n_unique_keys"] == reference.counts.n_unique_keys,
        entry["key_space_fingerprint"][:16],
    )


@check("the key file the installed CLI writes is byte-identical", "equivalence")
def _byte_identical(env, state):
    import hashlib

    sys.path.insert(0, os.path.join(PACKAGE, "src"))
    import biotsr as source_tsr
    from biotsr.io import native

    reference_dir = os.path.join(env.root, "reference")
    os.makedirs(reference_dir, exist_ok=True)
    result = source_tsr.analyze(
        os.path.join(FIXTURES, "1UBQ.pdb"), method="protein", chain="A",
        keep_triplets=False,
    )
    native.save(result, reference_dir, triplets=False)
    name = native.keys_filename(native.default_prefix(result))

    produced = os.path.join(env.cwd, "single", name)
    if not os.path.isfile(produced):
        return False, "no key file at %s" % produced
    a = hashlib.sha256(open(produced, "rb").read()).hexdigest()
    b = hashlib.sha256(open(os.path.join(reference_dir, name), "rb").read()).hexdigest()
    return a == b, a[:16]


@check("a batch run works outside the repository", "workflow")
def _batch(env, state):
    manifest = os.path.join(env.cwd, "samples.csv")
    open(manifest, "w").write(
        "sample_id,pdb_path,method,chain,label\n"
        "ubq,%s,protein,A,ubiquitin\n"
        "brs_a,%s,protein,A,barnase\n"
        "brs_b,%s,protein,B,barnase\n"
        % (
            os.path.join(FIXTURES, "1UBQ.pdb"),
            os.path.join(FIXTURES, "1BRS.pdb"),
            os.path.join(FIXTURES, "1BRS.pdb"),
        )
    )
    out = os.path.join(env.cwd, "batch")
    proc = env.run([env.biotsr, "batch", "-M", manifest, "-o", out, "--no-triplets"])
    ok = proc.returncode == 0 and os.path.isfile(
        os.path.join(out, "batch_summary.csv")
    )
    return ok, "exit %d" % proc.returncode


@check("a parallel batch (--jobs 2) matches the serial batch", "workflow")
def _parallel_batch(env, state):
    """Parallel execution works from an ordinary install, with no extra."""
    manifest = os.path.join(env.cwd, "samples.csv")
    serial, parallel = os.path.join(env.cwd, "batch"), os.path.join(env.cwd, "batch_jobs2")
    proc = env.run([env.biotsr, "batch", "-M", manifest, "-o", parallel,
                    "--no-triplets", "--jobs", "2"])
    if proc.returncode != 0:
        return False, "exit %d: %s" % (proc.returncode, proc.stderr.strip()[:200])
    differ, n = [], 0
    for sample in ("ubq", "brs_a", "brs_b"):
        for suffix in (".keys", ".meta.json"):
            name = sample + suffix
            n += 1
            a = open(os.path.join(serial, name), "rb").read()
            b = open(os.path.join(parallel, name), "rb").read()
            if a != b:
                differ.append(name)
    return n > 0 and not differ, "%d files byte-identical%s" % (
        n - len(differ), "; differ: %s" % differ if differ else "")


@check("a batch failure is attributable and exits 3", "workflow")
def _batch_failure(env, state):
    manifest = os.path.join(env.cwd, "bad.csv")
    open(manifest, "w").write(
        "sample_id,pdb_path,method,chain\n"
        "good,%s,protein,A\n"
        "bad,%s,protein,Q\n"
        % (os.path.join(FIXTURES, "1UBQ.pdb"), os.path.join(FIXTURES, "1UBQ.pdb"))
    )
    out = os.path.join(env.cwd, "batch_bad")
    proc = env.run([env.biotsr, "batch", "-M", manifest, "-o", out, "--no-triplets"])
    failures = os.path.join(out, "batch_failures.csv")
    return (
        proc.returncode == 3 and os.path.isfile(failures)
        and "bad" in open(failures).read(),
        "exit %d" % proc.returncode,
    )


@check("an analysis run works outside the repository", "workflow")
def _analysis(env, state):
    manifest = os.path.join(env.cwd, "samples.csv")
    out = os.path.join(env.cwd, "analysis")
    proc = env.run(
        [env.biotsr, "analyze", "common-keys", "-M", manifest, "-o", out,
         "--group-by", "label"]
    )
    root = os.path.join(out, "analysis", "common_keys")
    return (
        proc.returncode == 0 and os.path.isfile(os.path.join(root, "barnase.json")),
        "exit %d" % proc.returncode,
    )


@check("an SSE run works outside the repository", "workflow")
def _sse(env, state):
    out = os.path.join(env.cwd, "sse")
    proc = env.run(
        [env.biotsr, "run", "-m", "protein", "--pdb", os.path.join(FIXTURES, "1BRS.pdb"),
         "-c", "A", "--sse", "-o", out, "--no-triplets"]
    )
    if proc.returncode != 0:
        return False, proc.stderr.strip()[:200]
    payload = json.load(open(os.path.join(out, "1BRS_A_protein_sse.run.json")))
    return (
        payload["profile"]["sse_strand_identity"] == "sheet_id_and_serial",
        payload["profile"]["sse_strand_identity"],
    )


@check("a cross run works outside the repository", "workflow")
def _cross(env, state):
    out = os.path.join(env.cwd, "cross")
    proc = env.run(
        [env.biotsr, "run", "-m", "protein_protein_calpha_cross",
         "--pdb", os.path.join(FIXTURES, "1BRS.pdb"), "-c", "A", "--chain-2", "D",
         "--contact-cutoff", "8.0", "-o", out, "--no-triplets"]
    )
    return proc.returncode == 0, "exit %d" % proc.returncode


@check("a user error prints no traceback and exits 2", "errors")
def _user_error(env, state):
    proc = env.run(
        [env.biotsr, "run", "-m", "protein", "--pdb", "/nope.pdb", "-c", "A",
         "-o", os.path.join(env.cwd, "never")]
    )
    return (
        proc.returncode == 2 and "Traceback" not in proc.stderr
        and proc.stderr.startswith("biotsr: error:"),
        "exit %d" % proc.returncode,
    )


@check("nothing is written to the working directory", "isolation")
def _no_stray_writes(env, state):
    """Every path the CLI writes is one the user named."""
    expected = {
        "samples.csv", "bad.csv", "single", "batch", "batch_jobs2", "batch_bad",
        "analysis", "sse", "cross",
    }
    actual = set(os.listdir(env.cwd))
    stray = actual - expected
    return not stray, ", ".join(sorted(stray)) or "none"


@check("the CLI runs with the repository deleted from the path", "isolation")
def _no_repo_dependency(env, state):
    """A path-independence check with teeth: import from a directory that has
    a decoy ``biotsr`` package in it, and confirm the installed one still wins
    nothing from the checkout."""
    decoy = os.path.join(env.cwd, "decoy")
    os.makedirs(os.path.join(decoy, "validation"), exist_ok=True)
    proc = env.run([env.biotsr, "methods"], cwd=decoy)
    return proc.returncode == 0 and "protein" in proc.stdout, "exit %d" % proc.returncode


# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keep", action="store_true",
                        help="keep the temporary environment for inspection")
    parser.add_argument("--json", dest="json_path", default="")
    parser.add_argument("--artifact", default="",
                        help="install this wheel or sdist instead of the source directory")
    args = parser.parse_args()

    target = os.path.abspath(args.artifact) if args.artifact else PACKAGE
    root = tempfile.mkdtemp(prefix="tsr-clean-env-")
    results = []
    failed = 0
    try:
        env = Env(root, target)
        print("environment  %s" % env.venv)
        env.create()
        install = env.install()
        if install.returncode != 0:
            print("INSTALL FAILED")
            print(install.stdout[-4000:])
            print(install.stderr[-4000:])
            return 1
        print("installed    %s" % target)
        print()

        state = {}
        for name, category, fn in CHECKS:
            try:
                ok, detail = fn(env, state)
            except Exception as exc:  # noqa: BLE001 -- a check crashing is a failure
                ok, detail = False, "%s: %s" % (type(exc).__name__, exc)
            failed += not ok
            # Report the throwaway environment by role, never by its absolute
            # (machine- and user-specific) temporary path.
            detail = str(detail)
            for prefix in (os.path.realpath(root), root):
                detail = detail.replace(prefix, "<clean-env>")
            detail = detail.replace(WORKSPACE, "<workspace>")
            results.append(
                {"check": name, "category": category, "passed": bool(ok),
                 "detail": detail}
            )
            print("  %s %-58s %s" % ("PASS" if ok else "FAIL", name, detail))
    finally:
        if args.keep:
            print("\nkept %s" % root)
        else:
            shutil.rmtree(root, ignore_errors=True)

    print()
    print("%d checks, %d failed" % (len(results), failed))

    installed = (os.path.basename(target) if args.artifact
                 else os.path.relpath(PACKAGE, WORKSPACE))
    payload = {
        "package": installed,
        "python": sys.version.split()[0],
        "install_command": "pip install %s" % installed,
        "checks": results,
        "n_checks": len(results),
        "n_failed": failed,
    }
    # A record of one machine's install, not deterministic evidence: written
    # only where asked (release_artifacts.py passes a scratch path).
    if args.json_path:
        path = args.json_path
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print("-> %s" % path)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
