# Platform validation: Linux and macOS

BioTSR is supported on **Linux and macOS**; Windows is not supported (the
package is not tested there). Continuous integration runs the test
suite on Linux for every supported Python and on macOS arm64 (see
[RELEASING.md](RELEASING.md)). This page states what is identical on every
platform, what may differ, and how the test suite compares results.

## What is identical on every platform

- integer keys, per triangle and as a key-frequency mapping, including every
  per-category count;
- theta class, distance class and mirror sign;
- SSE categories and their totals;
- selected atoms and residues, candidate and triangle counts, and the
  canonical vertex order;
- the `KeySpace` and `SamplingProfile` fingerprints;
- the native keys file, which holds only integer keys, category names and
  integer frequencies.

None of these involves a tolerance, and the test suite requires the reference
keys files to match by SHA-256 on every platform it runs on.

## What may differ, and why

The optional native **triplet** file prints raw floating-point geometry on every
line — `theta`, `theta_raw` and `max_distance` — using Python's shortest
round-trip spelling of each value (`repr`). CPython computes that geometry with
the host C maths library (`acos`, and `pow` for `x ** 0.5` and `x ** 2`), and
those functions are not correctly rounded on every platform: macOS arm64 returns
a value one unit in the last place away from the correctly rounded result for a
measurable fraction of `acos` inputs, while glibc on Linux returned the
correctly rounded value on every input probed. The final printed digits of those
values can therefore differ between macOS and Linux, at the level of 10⁻¹⁴ to
10⁻¹² degrees or ångströms.

This is not nondeterminism. Each platform is fully deterministic — the same
platform always prints the same bytes — and the difference lies in how the host
library rounds an intermediate value. The classes are separated from their
boundaries by many orders of magnitude more than this: a geometry change of
scientific size, such as dividing by `3.14` instead of `math.pi`, moves theta by
up to 0.09°.

## How the test suite compares results

`tests/data/reference/standard_baseline.json` holds the standard results of 38
calculations over the structures in `tests/data/structures/`, covering all 20
methods and the supported modifiers. They were computed on macOS arm64 with
CPython 3.13. `tests/integration/test_standard_references.py` recomputes every
case and requires:

- on **every platform**: the SHA-256 of the native keys file, the point,
  triangle and distinct-key counts, both identity fingerprints, the residues used
  and the category totals;
- on **macOS arm64 only**: the SHA-256 of the native triplet file
  (`triplets_sha256_macos_arm64`), for cases of up to 120,000 triangles. On
  other platforms the integer keys carry the comparison.

Cases with more than 250,000 triangles run only with `pytest -m slow` (and in
release validation, which runs on macOS arm64).

## What this means for users

Results computed on Linux and on macOS are scientifically identical: the same
keys, frequencies, categories, counts and fingerprints. If byte-for-byte
identity of a native triplet file with one produced on another machine matters
to you, compare on the same platform, or compare the keys files, which are
identical everywhere.
