#!/usr/bin/env python3
"""Release check: build the release artifacts, inspect them, and install each one.

Every earlier clean-install check installed the package *directory*.  A user
installs a *built artifact* -- a wheel from PyPI, or an sdist that pip builds --
and the two can differ from the directory and from each other: a resource that
is in the tree but not declared as package data, a file the sdist forgot, a
stray file the wheel swept up.  So this script tests the artifacts themselves:

1. builds a wheel and an sdist from the repository root with ``python -m build``
   (an isolated PEP 517 build), with ``SOURCE_DATE_EPOCH`` pinned to the commit
   time so the wheel's bytes do not depend on when it was built;
2. inspects both archives: required files present, forbidden files absent, the
   packaged label table byte-identical to the source tree's, the metadata's
   name/version/licence/Python floor as declared, and no absolute workstation
   path inside any packaged text file;
3. unpacks the sdist and runs its own test suite there, outside the
   repository: it must pass, the computational tests must actually run on the
   structures the sdist carries, and every skip must be one of the expected
   kinds (repository-only evidence, untracked material, an optional
   dependency) -- so the sdist cannot "pass" by skipping;
4. runs ``twine check`` when twine is importable (a development tool only --
   never a BioTSR dependency);
5. installs **each artifact separately** into a fresh virtual environment with
   the repository nowhere on its path, and runs the clean-environment
   checks against it (``clean_install_check.py --artifact``): console script,
   ``python -m biotsr``, ``biotsr methods``, real single/batch/SSE/cross/
   analysis workflows run from outside the repository, the absence of every
   removed legacy module, and the installed CLI's fingerprints and key file
   against the source tree's.

Only the checked-out tree is read; the one thing taken from git is the time of
the ``HEAD`` commit (``git log -1``), for ``SOURCE_DATE_EPOCH``.

Usage::

    python3 validation/scripts/release_artifacts.py [--dist-dir DIR] [--json PATH]

Requires ``build`` and ``pytest`` in the running interpreter.  Exit 0 only when
every check passes.
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
import tarfile
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
WORKSPACE = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(WORKSPACE, "src")
LEXICONS = os.path.join(SRC, "biotsr", "labels", "resources")
CLEAN_ENV = os.path.join(HERE, "clean_install_check.py")

sys.path.insert(0, SRC)
import biotsr  # noqa: E402

#: The one packaged label table: element -> label, radix 112 (ADR-113).
REQUIRED_LEXICONS = (
    "element_labels_v1.csv",
)

#: Path fragments that must never appear inside a release artifact.
FORBIDDEN = (
    re.compile(r"(^|/)validation/"),
    re.compile(r"(^|/)legacy/"),
    re.compile(r"(^|/)docs/history/"),
    re.compile(r"(^|/)\.venv/|(^|/)venv/|(^|/)env/"),
    re.compile(r"(^|/)\.git/|(^|/)\.github/"),
    re.compile(r"__pycache__/|\.py[co]$"),
    re.compile(r"(^|/)\.DS_Store$|(^|/)\._"),
    re.compile(r"(^|/)\.pytest_cache/|(^|/)build/|(^|/)dist/"),
    re.compile(r"\.(pdb|ent|cif|xls|xlsx|xlsm|xlsb|ods|numbers|gz|zip)$"),
    # Private dataset tables and BioTSR's own run outputs.
    re.compile(r"(^|/)(sample_details[^/]*\.csv|column_means[^/]*)$"),
    re.compile(r"(^|/)(run\.json|batch_summary\.csv|batch_failures\.csv|fetch_report\.csv)$"),
    re.compile(r"\.(env|pem|key|netrc)$"),
    # The legacy writers, profiles, profile catalogue and lexicons that the
    # standard-only package removed; a stale tree or build directory would
    # otherwise carry them back in.
    re.compile(r"(^|/)biotsr/io/legacy_[^/]*\.py$"),
    re.compile(r"(^|/)biotsr/profiles/(legacy|catalogue)\.py$"),
    re.compile(r"(^|/)biotsr/labels/resources/atom_lexicon_v[12]\.csv$"),
    # The atom-name lexicon the element table replaced in 1.0.0rc7 (ADR-113).
    re.compile(r"(^|/)biotsr/labels/resources/atom_lexicon_standard_v1\.csv$"),
    re.compile(r"(^|/)biotsr/labels/lexicon\.py$"),
    re.compile(r"drug_atom_lexical"),
)

#: The exceptions, and only in the sdist: the unmodified PDB entries the
#: computational tests run on, and the one synthetic structure, by name, that
#: its generator in the same folder writes byte for byte.
SDIST_ALLOWED = re.compile(
    r"^tests/data/structures/[^/]+\.pdb$"
    r"|^tests/data/synthetic/SYNPHOS\.pdb$"
    # The one compressed test-reference file: one integer column, documented in
    # tests/data/reference/README.md.
    r"|^tests/data/reference/historical/1UBQ_A_TSR-Package_TSR_keys\.txt\.gz$"
)

#: A local absolute path in shipped text is a leak and a portability bug.
ABSOLUTE_PATH = re.compile(r"/Users/[A-Za-z0-9_.-]+/|/home/[A-Za-z0-9_.-]+/|[A-Z]:\\\\Users\\\\")

CHECKS = []


def record(name, ok, detail=""):
    # Only the checkout itself, not a sibling whose name merely starts with it.
    detail = re.sub(re.escape(WORKSPACE) + r"(?=[/\s\"']|$)", "<workspace>", str(detail))
    CHECKS.append({"check": name, "passed": bool(ok), "detail": detail})
    print("  %s %-64s %s" % ("PASS" if ok else "FAIL", name, detail))
    return ok


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    with open(path, "rb") as fh:
        return sha256_bytes(fh.read())


def source_date_epoch():
    try:
        out = subprocess.run(["git", "-C", WORKSPACE, "log", "-1", "--format=%ct"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             universal_newlines=True, check=True).stdout.strip()
        return out or None
    except (OSError, subprocess.CalledProcessError):
        return None


def build(outdir):
    env = dict(os.environ)
    epoch = source_date_epoch()
    if epoch:
        env["SOURCE_DATE_EPOCH"] = epoch
    # No --sdist/--wheel flags: build makes the sdist, then builds the wheel
    # FROM the unpacked sdist. That proves the sdist is buildable, and no stale
    # file in a local build/ directory can leak into the wheel.
    proc = subprocess.run(
        [sys.executable, "-m", "build", "--outdir", outdir, WORKSPACE],
        env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
    )
    wheels = sorted(f for f in os.listdir(outdir) if f.endswith(".whl")) if os.path.isdir(outdir) else []
    sdists = sorted(f for f in os.listdir(outdir) if f.endswith(".tar.gz")) if os.path.isdir(outdir) else []
    ok = proc.returncode == 0 and len(wheels) == 1 and len(sdists) == 1
    record("sdist builds from the repository root, and the wheel from the sdist", ok,
           ", ".join(wheels + sdists) if ok else proc.stdout[-600:])
    if not ok:
        return None, None, epoch
    return os.path.join(outdir, wheels[0]), os.path.join(outdir, sdists[0]), epoch


def _metadata_fields(text):
    fields = {}
    for line in text.split("\n\n", 1)[0].splitlines():
        if ":" in line and not line.startswith(" "):
            k, v = line.split(":", 1)
            fields.setdefault(k.strip(), []).append(v.strip())
    return fields


def _check_members(kind, names, read):
    bad = sorted(
        n for n in names
        if any(p.search(n) for p in FORBIDDEN)
        and not (kind == "sdist" and SDIST_ALLOWED.search(n))
    )
    if kind == "wheel":
        bad += sorted(n for n in names if n.startswith("tests/"))
    record("%s contains no forbidden file" % kind, not bad,
           "none" if not bad else ", ".join(bad[:5]))
    leaked = []
    for n in names:
        if n.endswith((".py", ".md", ".txt", ".csv", ".json", ".toml", ".cfg", ".in", "METADATA", "RECORD")):
            data = read(n)
            if data is not None and ABSOLUTE_PATH.search(data.decode("utf-8", "replace")):
                leaked.append(n)
    record("%s contains no absolute workstation path" % kind, not leaked,
           "none" if not leaked else ", ".join(leaked[:5]))


#: The licence and notice files the repository ships (Apache-2.0, ADR-085).
LICENCE_FILES = ("LICENSE", "NOTICE")


def _check_licence_files(kind, names, read):
    """Every licence file the repository has is in the artifact, byte for byte."""
    missing, differ = [], []
    for notice in LICENCE_FILES:
        path = os.path.join(WORKSPACE, notice)
        if not os.path.exists(path):
            continue
        found = [n for n in names if n == notice or n.endswith("/" + notice)]
        if not found:
            missing.append(notice)
        elif any(sha256_bytes(read(n)) != sha256_file(path) for n in found):
            differ.append(notice)
    present = [n for n in LICENCE_FILES if os.path.exists(os.path.join(WORKSPACE, n))]
    record("%s carries %s byte-identical to the repository"
           % (kind, " and ".join(present) or "no licence file (none in the repository)"),
           not missing and not differ,
           "missing %s; differ %s" % (missing, differ) if (missing or differ) else "sha256 match")


def _check_lexicons(kind, prefix, names, read):
    missing, differ = [], []
    for lex in REQUIRED_LEXICONS:
        member = prefix + "biotsr/labels/resources/" + lex
        if member not in names:
            missing.append(lex)
            continue
        if sha256_bytes(read(member)) != sha256_file(os.path.join(LEXICONS, lex)):
            differ.append(lex)
    resources = prefix + "biotsr/labels/resources/"
    extra = sorted(n[len(resources):] for n in names
                   if n.startswith(resources) and n.endswith(".csv")
                   and n[len(resources):] not in REQUIRED_LEXICONS)
    record("%s carries exactly the %d packaged label table(s), byte-identical to the source tree"
           % (kind, len(REQUIRED_LEXICONS)),
           not missing and not differ and not extra,
           "missing %s; differ %s; unexpected %s" % (missing, differ, extra)
           if (missing or differ or extra)
           else "%d/%d sha256 match" % (len(REQUIRED_LEXICONS), len(REQUIRED_LEXICONS)))


def _check_metadata(kind, text):
    f = _metadata_fields(text)
    expected_license = _declared_license()
    problems = []
    if f.get("Name", [""])[0] != "biotsr":
        problems.append("Name=%r" % f.get("Name"))
    if f.get("Version", [""])[0] != biotsr.__version__:
        problems.append("Version=%r != %r" % (f.get("Version"), biotsr.__version__))
    if f.get("Requires-Python", [""])[0] != ">=3.9":
        problems.append("Requires-Python=%r" % f.get("Requires-Python"))
    runtime = [r.replace(" ", "") for r in f.get("Requires-Dist", []) if "extra ==" not in r]
    if runtime != ["joblib>=1.3"]:     # ADR-086: joblib, and nothing else
        problems.append("runtime Requires-Dist=%r" % f.get("Requires-Dist"))
    if "parallel" not in f.get("Provides-Extra", []):
        problems.append("the biotsr[parallel] alias is not declared")
    lic = (f.get("License-Expression") or f.get("License") or [""])[0]
    if lic != expected_license:
        problems.append("licence %r != declared %r" % (lic, expected_license))
    urls = " ".join(f.get("Project-URL", []))
    for url in ("https://github.com/TSRMethod/BioTSR", "https://tsrmethod.github.io/"):
        if url not in urls:
            problems.append("missing Project-URL %s" % url)
    record("%s metadata: name, version, licence, Python floor, URLs, joblib the only runtime dep" % kind,
           not problems, "; ".join(problems) or "biotsr %s, %s" % (biotsr.__version__, lic))
    return f


def _declared_license():
    with open(os.path.join(WORKSPACE, "pyproject.toml"), encoding="utf-8") as fh:
        text = fh.read()
    m = re.search(r'^license\s*=\s*"([^"]+)"', text, re.M) or \
        re.search(r'^license\s*=\s*\{\s*text\s*=\s*"([^"]+)"', text, re.M)
    return m.group(1) if m else ""


def inspect_wheel(path):
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        read = lambda n: zf.read(n) if n in names else None  # noqa: E731
        top = sorted({n.split("/", 1)[0] for n in names})
        dist_info = "biotsr-%s.dist-info" % biotsr.__version__
        record("wheel top level is exactly biotsr/ and its dist-info",
               top == ["biotsr", dist_info], ", ".join(top))
        _check_members("wheel", names, read)
        _check_lexicons("wheel", "", names, read)
        f = _check_metadata("wheel", read(dist_info + "/METADATA").decode("utf-8"))
        ep = (read(dist_info + "/entry_points.txt") or b"").decode("utf-8")
        record("wheel declares the biotsr console script",
               re.search(r"^biotsr\s*=\s*biotsr\.cli:main\s*$", ep, re.M) is not None,
               ep.strip().replace("\n", " | "))
        licenses = [n for n in names if "/licenses/" in n or n.endswith(("LICENSE", "NOTICE"))]
        has_license = os.path.exists(os.path.join(WORKSPACE, "LICENSE"))
        record("wheel ships the licence file iff the repository has one",
               bool(licenses) == has_license,
               ", ".join(licenses) if licenses else "no LICENSE in repository or wheel")
        _check_licence_files("wheel", names, read)
        py = [n for n in names if n.endswith(".py")]
        src_py = sum(1 for d, _s, fs in os.walk(os.path.join(SRC, "biotsr"))
                     for x in fs if x.endswith(".py") and "__pycache__" not in d)
        record("wheel carries every source module", len(py) == src_py,
               "%d of %d modules" % (len(py), src_py))
        return {"file": os.path.basename(path), "sha256": sha256_file(path),
                "bytes": os.path.getsize(path), "members": len(names),
                "metadata_version": f.get("Metadata-Version", [""])[0]}


def _tracked_files():
    """The repository's tracked files, or ``None`` outside a git checkout."""
    try:
        out = subprocess.run(["git", "-C", WORKSPACE, "ls-files", "-z"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    return {p for p in out.decode("utf-8").split("\0") if p}


def _check_tracked(kind, rel):
    """Every source and test file in the artifact is a tracked file.

    A forbidden-name list catches what it names; this catches anything else a
    working tree may hold -- notes, a dataset table, a stray output -- because
    only what the repository tracks may ship.
    """
    tracked = _tracked_files()
    if tracked is None:
        record("%s carries only tracked source and test files" % kind, True,
               "not checked: no git checkout")
        return
    extra = sorted(n for n in rel if n.startswith(("src/biotsr/", "tests/"))
                   and n not in tracked)
    record("%s carries only tracked source and test files" % kind, not extra,
           "none untracked" if not extra else ", ".join(extra[:5]))


def inspect_sdist(path):
    with tarfile.open(path, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile()]
        names = [m.name for m in members]
        root = "biotsr-%s/" % biotsr.__version__

        def read(n):
            try:
                return tf.extractfile(n).read()
            except KeyError:
                return None

        record("sdist unpacks into one versioned directory",
               all(n.startswith(root) for n in names), root)
        rel = [n[len(root):] for n in names]
        required = ["pyproject.toml", "README.md", "src/biotsr/__init__.py", "PKG-INFO"]
        for notice in LICENCE_FILES:
            if os.path.exists(os.path.join(WORKSPACE, notice)):
                required.append(notice)
        for extra in ("CITATION.cff", "CHANGELOG.md"):
            if os.path.exists(os.path.join(WORKSPACE, extra)):
                required.append(extra)
        missing = [r for r in required if r not in rel]
        record("sdist carries pyproject, README, source, metadata%s"
               % "".join(", " + n for n in LICENCE_FILES if n in required),
               not missing, "missing %s" % missing if missing else ", ".join(required))
        _check_licence_files("sdist", names, read)
        _check_members("sdist", rel, lambda n: read(root + n))
        _check_tracked("sdist", rel)
        _check_lexicons("sdist", root + "src/", names, read)
        f = _check_metadata("sdist", read(root + "PKG-INFO").decode("utf-8"))
        tests = [r for r in rel if r.startswith("tests/")]
        structures = [t for t in tests if t.startswith("tests/data/structures/")
                      and t.endswith(".pdb")]
        record("sdist carries the test suite and its structures",
               any(t.endswith("conftest.py") for t in tests) and len(structures) >= 11,
               "%d test files, %d structures" % (len(tests), len(structures)))
        return {"file": os.path.basename(path), "sha256": sha256_file(path),
                "bytes": os.path.getsize(path), "members": len(names),
                "metadata_version": f.get("Metadata-Version", [""])[0]}


#: Every skip in an unpacked sdist's test run must start with one of these
#: (tests/conftest.py defines the first three; the opt-in check over an external
#: dataset sheet, which runs only when its environment variable names one,
#: defines the fourth in tests/integration/test_dataset_sheet_rows.py).
EXPECTED_SKIPS = (
    "repository-only evidence not present",
    "untracked optional material not present",
    "optional dependency not installed",
    "opt-in check not requested",
)
#: The unpacked sdist must actually run at least this many tests: its structures
#: are in tests/data, so the computational tests have what they need.
MIN_SDIST_PASSED = 1500


def test_sdist(path, scratch):
    """Unpack the sdist and run its test suite from inside it."""
    try:
        import pytest  # noqa: F401
    except ImportError:
        record("the unpacked sdist's test suite passes", False,
               "pytest is not installed in this interpreter")
        return {}
    target = os.path.join(scratch, "sdist-tests")
    with tarfile.open(path, "r:gz") as tf:
        # The "data" filter (refuse absolute paths, links out of the tree,
        # special files) where this Python has it; it is the 3.14 default.
        if hasattr(tarfile, "data_filter"):
            tf.extractall(target, filter="data")
        else:
            tf.extractall(target)
    root = os.path.join(target, "biotsr-%s" % biotsr.__version__)
    env = {k: v for k, v in os.environ.items()
           if k not in ("PYTHONPATH", "BIOTSR_REQUIRE_EVIDENCE",
                        "BIOTSR_REQUIRE_FROZEN_TEXT", "BIOTSR_FROZEN_TEXT_DIR")}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run([sys.executable, "-m", "pytest", "-q", "-rfEs", "-p", "no:cacheprovider"],
                          cwd=root, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    tail = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
    counts = {k: int(v) for v, k in re.findall(r"(\d+) (passed|failed|skipped|errors?)", tail)}
    unexpected = sorted({
        m.group(1).strip() for m in re.finditer(r"^SKIPPED \[\d+\] [^:]+:\d+: (.*)$",
                                                proc.stdout, re.M)
        if not m.group(1).strip().startswith(EXPECTED_SKIPS)
    })
    # Name what failed first: the end of the output is mostly skip lines.
    failed = re.findall(r"^(?:FAILED|ERROR) (\S+)", proc.stdout, re.M)
    record("the unpacked sdist's test suite passes", proc.returncode == 0,
           tail if proc.returncode == 0 else
           ("failed: %s | %s" % (", ".join(failed), tail) if failed
            else proc.stdout[-800:]).replace(root, "<sdist>"))
    record("the unpacked sdist runs its computational tests (at least %d pass)"
           % MIN_SDIST_PASSED, counts.get("passed", 0) >= MIN_SDIST_PASSED,
           "%d passed" % counts.get("passed", 0))
    record("every skip in the unpacked sdist is an expected kind", not unexpected,
           "; ".join(unexpected)[:600] if unexpected else
           "%d skipped, all repository-only, untracked, optional-dependency or opt-in"
           % counts.get("skipped", 0))
    return {"summary": tail, "counts": counts}


def twine_check(paths):
    try:
        import twine  # noqa: F401
    except ImportError:
        record("twine check", True, "skipped: twine not installed in this interpreter")
        return "skipped"
    proc = subprocess.run([sys.executable, "-m", "twine", "check", "--strict"] + list(paths),
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    # Report counts, not twine's own lines: they carry the absolute path of the
    # build directory (and colour codes), and this detail is committed evidence.
    passed = proc.stdout.count("PASSED")
    record("twine check --strict passes on both artifacts", proc.returncode == 0,
           "%d of %d artifacts PASSED" % (passed, len(paths)) if proc.returncode == 0
           else re.sub(r"\x1b\[[0-9;]*m", "", proc.stdout).replace(
               os.path.dirname(paths[0]), "<dist>")[-300:])
    return "passed" if proc.returncode == 0 else "failed"


def clean_install(kind, path, scratch):
    out = os.path.join(scratch, "clean_env_%s.json" % kind)
    proc = subprocess.run([sys.executable, CLEAN_ENV, "--artifact", path, "--json", out],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    print("\n".join("    " + l for l in proc.stdout.rstrip().splitlines()))
    payload = {}
    if os.path.exists(out):
        with open(out) as fh:
            payload = json.load(fh)
    n, failed = payload.get("n_checks", 0), payload.get("n_failed", -1)
    record("the %s installs into a clean environment and passes every clean-env check" % kind,
           proc.returncode == 0 and n > 0 and failed == 0, "%d checks, %d failed" % (n, failed))
    return {"n_checks": n, "n_failed": failed,
            "checks": payload.get("checks", [])}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dist-dir", default="", help="keep the built artifacts here")
    parser.add_argument("--json", dest="json_path", default="")
    args = parser.parse_args()

    scratch = tempfile.mkdtemp(prefix="biotsr-release-")
    outdir = os.path.abspath(args.dist_dir) if args.dist_dir else os.path.join(scratch, "dist")
    if args.dist_dir and os.path.isdir(outdir):
        for f in os.listdir(outdir):
            if f.startswith("biotsr-") and f.endswith((".whl", ".tar.gz")):
                os.remove(os.path.join(outdir, f))
    os.makedirs(outdir, exist_ok=True)
    report = {"biotsr_version": biotsr.__version__}
    try:
        print("build")
        wheel, sdist, epoch = build(outdir)
        report["source_date_epoch"] = epoch
        if wheel:
            print("\ninspect wheel")
            report["wheel"] = inspect_wheel(wheel)
            print("\ninspect sdist")
            report["sdist"] = inspect_sdist(sdist)
            print("\ntest suite inside the unpacked sdist")
            report["sdist_tests"] = test_sdist(sdist, scratch)
            print("\ntwine")
            report["twine_check"] = twine_check([wheel, sdist])
            for kind, path in (("wheel", wheel), ("sdist", sdist)):
                print("\nclean install from the %s" % kind)
                report["clean_install_" + kind] = clean_install(kind, path, scratch)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    failed = [c for c in CHECKS if not c["passed"]]
    report["checks"] = CHECKS
    report["n_checks"] = len(CHECKS)
    report["n_failed"] = len(failed)
    print("\n%d release-artifact checks, %d failed" % (len(CHECKS), len(failed)))
    if args.json_path:
        os.makedirs(os.path.dirname(os.path.abspath(args.json_path)), exist_ok=True)
        with open(args.json_path, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print("-> %s" % args.json_path)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
