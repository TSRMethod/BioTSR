# Changelog

All notable user-visible changes to BioTSR. Versions follow
[PEP 440](https://peps.python.org/pep-0440/); the release process is in
[docs/RELEASING.md](https://github.com/TSRMethod/BioTSR/blob/main/docs/RELEASING.md).

## [1.0.0]

The first release of BioTSR, the unified software implementation of the
Triangular Spatial Relationship (TSR) methodology: one Python package, one
command line and one shared scientific core — the TSR Engine — in place of a
family of independent historical TSR programs. Every method runs under one
standard profile, the accepted TSR science as set by the TSR Research Group,
and writes one self-describing output format.

### Platforms

- **CPython 3.9, 3.10, 3.11, 3.12 and 3.13, on Linux and macOS.** The full
  default test suite runs on every one of them on Linux (x86-64) and on macOS
  (arm64); the slow reference cases over large structures run on macOS arm64.
  Keys, frequencies and fingerprints are identical on both systems.
- **Windows is not supported**: it is neither tested nor validated.
- Input is PDB format; mmCIF is not read.

### Installation

- **From GitHub:** `pip install "biotsr @ git+https://github.com/TSRMethod/BioTSR@v1.0.0"`,
  or the wheel and source distribution attached to the
  [1.0.0 release](https://github.com/TSRMethod/BioTSR/releases/tag/v1.0.0); from a
  clone, `pip install .`.
- **From PyPI:** `pip install biotsr`, available when
  <https://pypi.org/project/biotsr/> lists version 1.0.0.
- One runtime dependency, joblib, runs batch rows in parallel; the calculation
  itself uses the standard library only. The optional `plot` extra
  (matplotlib) draws the clustered heat map.

### What 1.0.0 provides

- **Twenty methods in three families:** protein Cα; atomic — amino-acid,
  ligand, carbohydrate, nucleotide/DNA/RNA, and residue or atom selections;
  and five cross (two-sided) methods — drug–protein, carbohydrate–protein,
  nucleotide–protein, protein–protein all-atom and protein–protein Cα. The
  [capability matrix](https://github.com/TSRMethod/BioTSR/blob/main/docs/CAPABILITY_MATRIX.md)
  says what each supports.
- **One TSR Engine** — geometry, canonical vertex ordering, angle and distance
  binning, labelling and integer key encoding shared by every method — with
  composable modifiers: amino-acid grouping and SSE-TSR (18 secondary-structure
  categories) on the Cα methods, mirror signing and size filtering.
- **Labels.** Cα points are labelled by residue — the 20 regular amino acids
  and eleven modified residues — at radix 30, collision-free; atoms are
  labelled by their PDB element field at radix 112, from one packaged,
  validated element table.
- **Strict, recorded input rules.** Atom names and elements are read from their
  full PDB fields, nothing is guessed, and a selection with an unknown residue
  (`UNK`), an atom without an element, an atom written twice or an unreadable
  residue number is refused with the records named.
- **Two identities on every result** — a `KeySpace` fingerprint (what the keys
  mean) and a `SamplingProfile` fingerprint (how the triangles were drawn).
  Results whose keys are defined differently are refused, never combined;
  results processed under other settings, and partial results, are compared
  only when acknowledged (`--allow-mixed-processing`, `--allow-partial`).
- **Native output** — `PREFIX.keys`, `PREFIX.triplets`, `PREFIX.meta.json` and
  a completion record `PREFIX.run.json` written last — read back by
  `biotsr.load_result` with both fingerprints restored.
- **Datasets.** `biotsr batch` over manifests and wide dataset sheets, with
  `--dry-run`, verified resuming (`--if-exists skip`), parallel rows
  (`--jobs`), a per-request `batch_summary.csv`, and opt-in, verified download
  of PDB entries (`biotsr fetch`, `--download-missing`).
- **Analysis** over results of any method family: feature matrices, common
  keys, five similarity measures, the average-linkage tree, dendrogram and
  ordered heat map, and flat clusters by an explicit rule.
- **Licensed under the Apache License 2.0.** Third-party material keeps its
  own terms; its notices are in
  [`NOTICE`](https://github.com/TSRMethod/BioTSR/blob/main/NOTICE).

### Known limitations

The [known limitations](https://github.com/TSRMethod/BioTSR/blob/main/docs/KNOWN_LIMITATIONS.md)
state each one precisely. In short: only the standard science is computed —
the historical programs' outputs are not reproduced; only model 1 of a PDB
file is read, with alternate locations blank or `A`; a `TER` record without a
chain identifier ends no chain; atomic methods read three modified amino acids
written as `HETATM` (TPO, SEP, PTR) and name the others; triplet files are
compared byte for byte only on macOS arm64, where their reference geometry was
produced; and similarity and clustering are available from Python, not from
the command line.

### Results from the pre-release builds

1.0.0 reads and compares the results of the last pre-release build as its own.
It refuses — it does not silently combine — older results whose keys mean
something else: Cα results from builds that labelled residues at radix 20, and
atomic results from builds that labelled atoms by name, must be computed
again, and outputs in the earlier per-sample directory layout are not read.
What changed from one pre-release build to the next is recorded in
[docs/PRE_RELEASE_CHANGES.md](https://github.com/TSRMethod/BioTSR/blob/main/docs/PRE_RELEASE_CHANGES.md).

[1.0.0]: https://github.com/TSRMethod/BioTSR/releases/tag/v1.0.0
