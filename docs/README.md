# BioTSR documentation

BioTSR is the unified software implementation of the Triangular Spatial
Relationship (TSR) methodology. This page is the map of its documentation. If
you are new, start with the tutorials.

## Getting started

| | |
|---|---|
| [Tutorials](tutorials/README.md) | eight short, runnable walk-throughs of the established workflows, from one structure to a clustered heat map |
| [User guide](USER_GUIDE.md) | every method with a runnable example, and the dataset workflow: downloading structures, manifests, parallel batches, resume, memory |

## Reference

| | |
|---|---|
| [Command-line guide](CLI_GUIDE.md) | every `biotsr` command and option, the manifest format, downloading structures, parallel batches, exit codes, HPC use |
| [Python API guide](API_GUIDE.md) | run, inspect, save, load and analyse from Python |
| [Package overview](PACKAGE_OVERVIEW.md) | a longer tour of what the package contains, with Python examples per family |
| [Capability matrix](CAPABILITY_MATRIX.md) | which method supports which modifier and target argument, and its standard profile, generated from the method registry |

## Reproducibility

| | |
|---|---|
| [Reproducibility note](REPRODUCIBILITY_NOTE.md) | which command computes which published method |
| [Platform validation](PLATFORM_VALIDATION.md) | what is validated on Linux and macOS, and how |
| [Reference evidence](REFERENCE_EVIDENCE.md) | every reference calculation the tests compare against, what is compared, and which evidence is regression and which independent |
| [Known limitations](KNOWN_LIMITATIONS.md) | what version 1.0 does not do, stated precisely |

## Architecture

| | |
|---|---|
| [Architecture overview](ARCHITECTURE.md) | how the package is built: the pipeline, the TSR Engine, profiles, identities, serialization |
| [Architecture decisions](architecture/README.md) | ADR-001 to ADR-100: why BioTSR behaves as it does |

## Methods

| | |
|---|---|
| [Methods](methods/README.md) | cross TSR, SSE-TSR and analysis: the definitions, conventions and parameters the standard profiles implement |

## Development and provenance

| | |
|---|---|
| [Releasing](RELEASING.md) | versioning, Python support policy, CI, what ships in the package |
| [Licensing and provenance](LICENSING.md) | where each part of BioTSR comes from, under which terms, and what each distribution contains |

The reference outputs and the validation scripts that check the standard
results are under [`validation/`](../validation/) in the repository, and the test structures are under
[`tests/data/structures/`](../tests/data/structures/). They are not part of the
installed package.
