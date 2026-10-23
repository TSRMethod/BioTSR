"""The argument parser.  One executable, one vocabulary.

Every long flag is the ``analyze()``/``analyze_cross()`` keyword it maps to, with
underscores written as dashes and nothing abbreviated: ``size_filter=25.0``
becomes ``--size-filter 25.0``, everywhere, for every method.  A user who knows
the Python API knows the CLI, and vice versa -- which is the property the nine
legacy repositories most conspicuously lack, where the same concept is
``--dlen``, ``size_filter`` and ``distance_gap`` in three different programs.

Nothing here carries a scientific default.  ``--mirror`` defaults to *unset*, not
to ``False``: what the user did not ask for is not passed, and the package's own
default applies.  That is why ``--size-filter`` has no default rather than
``10000`` (the legacy sentinel for "off"), and why ``--contact-cutoff`` has none
rather than each family's historical value.
"""

from __future__ import annotations

import argparse

from ..retrieval import DEFAULT_RETRIES, DEFAULT_TIMEOUT, DEFAULT_WORKERS, MAX_WORKERS
from .commands import ON_ERROR_CONTINUE, ON_ERROR_STOP
from .output import IF_EXISTS_CHOICES, IF_EXISTS_FAIL

__all__ = ["build_parser"]

_EPILOG_MAIN = """\
examples:
  biotsr methods
  biotsr methods protein
  biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/
  biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/
  biotsr batch --manifest sheet.csv --sheet tsr --method carbohydrate \
      --method aa_ind --pdb-dir pdb/ --output results/
  biotsr fetch 1UBQ 4HHB --pdb-dir pdb/
  biotsr describe results/
  biotsr describe results/1UBQ_A_protein

exit codes:
  0    success
  1    unexpected internal error
  2    user, input or configuration error
  3    a batch finished with requests that failed, were stopped or never started
       (see batch_summary.csv and batch_failures.csv), or some structures could
       not be fetched (see fetch_report.csv)
  4    finished, nothing failed, but some result is partial: its selection
       lacks something it named (written and usable; an analysis then needs
       --allow-partial)
  130  interrupted (Ctrl-C or SIGTERM); a batch writes its summary first, and
       fetch its report

documentation:
  https://github.com/TSRMethod/BioTSR/blob/main/docs/USER_GUIDE.md
"""

_EPILOG_RUN = """\
examples:
  # ordinary protein C-alpha TSR
  biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/

  # the three C-alpha modifiers, composed
  biotsr run --method protein --pdb 1UBQ.pdb --chain A --output results/ \\
      --aa-grouping --mirror --size-filter 25.0

  # secondary-structure annotation (keys are unchanged; categories are added)
  biotsr run --method protein --pdb 1TUP.pdb --chain A --output results/ --sse

  # one ligand, named by its residue identifier
  biotsr run --method drug --pdb 3PTB.pdb --chain A --output results/ \\
      --entity-name BEN --entity-id 1

  # every NAG on chain D, in one combined calculation (prefix ..._NAG_allcopies)
  biotsr run --method carbohydrate --pdb 7LF9.pdb --chain D --output results/ \\
      --entity-name NAG

  # every residue of a chain, atomically, one calculation per residue
  biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A --output results/

  # every leucine of a chain, one calculation per leucine
  biotsr run --method amino_acid_atomic --pdb 1UBQ.pdb --chain A \\
      --residue-name LEU --output results/

  # a protein-protein interface
  biotsr run --method protein_protein_cross --pdb 1BRS.pdb \\
      --chain A --chain-2 D --contact-cutoff 5.0 --output results/
"""

_EPILOG_BATCH = """\
Every calculation writes PREFIX.keys, PREFIX.meta.json, PREFIX.triplets (unless
--no-triplets) and its completion record PREFIX.run.json into the --output
folder; batch_summary.csv lists every request with its source row(s).

1. A manifest: one row is one run, with the same names as the flags above,
   minus the dashes. A blank cell means "not set", never a default.

  sample_id,pdb_id,method,chain,mirror,size_filter
  ubq_regular,1UBQ,protein,A,,
  ubq_mirror,1UBQ,protein,A,true,25.0

   The method may instead be chosen on the command line, for every row:
   --method NAME (repeatable: each row then runs once per method). A manifest
   whose rows name methods together with --method is refused. Columns with
   other names are read with --column FIELD=HEADER, carried along with
   --keep-column or dropped with --ignore-column; any other unrecognised
   column is an error. Every request needs its chain (a chain column,
   --column chain=HEADER, or --chain ID for rows that leave it blank), and a
   cross method its chain_2; a request still without one is not run, gets the
   status chain_not_specified in batch_summary.csv, and the batch exits 4.

2. A wide dataset sheet (--sheet): one row per structure, several participant
   chains and alternative selections per row, and no method column. The
   methods come from --method (required, repeatable); each reads only its own
   columns. --sheet tsr reads the TSR dataset-sheet columns (protein, aa_chain,
   drug_chain, carb_chain, nt_chain, aa_chain_2, aa_ind, carb_name, carb_id,
   car_atom_selected ...); --sheet custom reads only what --map gives.
   --map METHOD[.PARTICIPANT].FIELD=HEADER maps any column, e.g.
   --map carbohydrate.chain=glycan_chain. Participant 2 (the *_2 columns of a
   second protein chain) runs with --participant 2 or all. A cross method
   reads its chains and molecule from the row; a side that the row could
   narrow in several ways (aa_ind, aa_ind_specific ... for the protein side)
   is chosen with --cross-side METHOD.SIDE=COLUMN, or =whole, and the run is
   refused until it is. A row that does not describe a method's calculation
   fails, or with --inapplicable skip is recorded as inapplicable. Rows asking
   for the same calculation of the same sample share it; rows with different
   labels for one calculation are refused (give them distinct sample ids).

examples:
  biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/
  # chains.csv (pdb_id,chain) has no method column: every row, once per method
  biotsr batch --manifest chains.csv --pdb-dir pdb/ --output results/ \\
      --method protein --method amino_acid_atomic --no-triplets

  # a TSR dataset sheet, two methods, preview first
  biotsr batch --manifest sample_details.csv --sheet tsr \\
      --method carbohydrate --method aa_ind --inapplicable skip \\
      --pdb-dir pdb/ --download-missing --output results/ --dry-run

  # drug-protein cross from the same sheet: the protein side narrowed by
  # aa_ind, the drug whole; the preview shows both sides of every request
  biotsr batch --manifest sample_details.csv --sheet tsr \
      --method drug_protein_cross --cross-side drug_protein_cross.protein=aa_ind \
      --cross-side drug_protein_cross.drug=whole --inapplicable skip \
      --pdb-dir pdb/ --output results/ --dry-run

  # the same, computed with four worker processes, keys only
  biotsr batch --manifest sample_details.csv --sheet tsr \\
      --method carbohydrate --method aa_ind --inapplicable skip \\
      --pdb-dir pdb/ --download-missing --output results/ --no-triplets --jobs 4

  # resume an interrupted run: calculations verified to be the same request
  # are kept, everything else (including half-installed ones) is computed
  biotsr batch --manifest samples.csv --pdb-dir pdb/ --output results/ \\
      --if-exists skip
"""

_EPILOG_FETCH = """\
Downloads PDB-format entries from RCSB over HTTPS as .pdb.gz, verifies the
transfer, the gzip integrity check and the decompressed file, and installs each
file atomically as <ID>.pdb. Files already in --pdb-dir are reused when valid and
never overwritten unless --refresh is given; a failed refresh keeps the old file.
Every id gets a row in fetch_report.csv.

examples:
  biotsr fetch 1UBQ 4HHB --pdb-dir pdb/
  biotsr fetch --id-file ids.txt --pdb-dir pdb/
  biotsr fetch --manifest sample_details.csv --pdb-dir pdb/ --download-workers 8
  biotsr fetch --manifest other.csv --column pdb_id=PDB --pdb-dir pdb/
  biotsr fetch --manifest sample_details.csv --pdb-dir pdb/ --dry-run
"""

_EPILOG_ANALYZE = """\
Reports the keys every sample shares. With --manifest the samples are computed
first, exactly as `biotsr batch` would compute them, keys only (no triplet
files), and read back from the outputs written; with --from-results an earlier
run's outputs are read and nothing is computed.

Either way the samples must share one key definition (key space: label
scheme, bins, encoding, mirror convention); that is never relaxed, and it is
decided from the key space, not from method names. Their processing settings
(size filter, contact cutoff and rule, completeness rule, candidate strategy,
alternate-location rule) must also agree, unless --allow-mixed-processing is
given; the report then records every sample's settings with a warning. What was
selected never prevents a comparison: structures, chains, residue types (LYS
versus TYR), ligands, carbohydrates, residue numbers, atom subsets and copies
may all differ, and samples need not share any key. all_samples.json (or
<group>.json) records every sample's method, selection, settings and
fingerprints.

examples:
  biotsr analyze common-keys --manifest samples.csv --pdb-dir pdb/ --output results/
  biotsr analyze common-keys --manifest samples.csv --pdb-dir pdb/ --output results/ \\
      --group-by label
"""


class _Formatter(argparse.RawDescriptionHelpFormatter):
    """Keep the epilog's examples verbatim and the option help wrapped."""

    def __init__(self, prog):
        super().__init__(prog, max_help_position=32, width=88)


def _add_common(parser, *, triplets=True):
    group = parser.add_argument_group("output")
    group.add_argument(
        "--output", "-o", required=True, metavar="DIR",
        help="the shared output folder: every calculation writes PREFIX.keys, "
             "PREFIX.meta.json, PREFIX.triplets and PREFIX.run.json into it. "
             "Required: the CLI never writes to the working directory.",
    )
    if triplets:
        group.add_argument(
            "--no-triplets", action="store_true",
            help="write the key frequencies but not the per-triangle file. The "
                 "keys are identical either way.",
        )
    group.add_argument(
        "--if-exists", choices=IF_EXISTS_CHOICES, default=IF_EXISTS_FAIL,
        help="what to do when a calculation's files are already there: 'fail' "
             "(default) changes nothing, 'overwrite' replaces them, 'skip' keeps "
             "them only when PREFIX.run.json is complete, records exactly the "
             "same request (structure content, method, chains, selections, "
             "modifiers, contact settings) and every file is present at the size "
             "it was written; otherwise that calculation fails. A calculation an "
             "interrupted run left incomplete is computed again under 'skip' and "
             "'overwrite'. Files BioTSR did not write are never replaced.",
    )

    group = parser.add_argument_group("reporting")
    group.add_argument("--quiet", "-q", action="store_true",
                       help="print errors only.")
    group.add_argument("--verbose", "-v", action="store_true",
                       help="print one line per calculation and the package's "
                            "own warnings.")
    group.add_argument("--debug", action="store_true",
                       help="print a traceback for unexpected internal errors. "
                            "Ordinary input errors never print one.")


def _add_science(parser, *, batch: bool):
    group = parser.add_argument_group("method")
    if batch:
        group.add_argument(
            "--method", "-m", metavar="NAME", action="append",
            help="a TSR method to run on every row; repeat it to run several "
                 "(each row then becomes one request per method). Required with "
                 "--sheet; for a manifest, use it instead of a method column, not "
                 "with one. `biotsr methods` lists them.",
        )
    else:
        group.add_argument(
            "--method", "-m", metavar="NAME",
            help="the TSR method to run (required). `biotsr methods` lists them; "
                 "each runs under its own standard profile.",
        )

    group = parser.add_argument_group("modifiers")
    group.add_argument("--aa-grouping", dest="aa_grouping", action="store_true",
                       default=None,
                       help="use the grouped residue labels (13 groups of the 20 "
                            "regular amino acids; MSE with MET; each other "
                            "supported modified residue on its own) instead of "
                            "the individual ones. Protein C-alpha methods only.")
    group.add_argument("--mirror", action="store_true", default=None,
                       help="negate a key whose pre-fold angle exceeds 90 degrees. "
                            "Off unless asked for; with it off no key is ever "
                            "negative.")
    group.add_argument("--size-filter", dest="size_filter", type=float,
                       metavar="ANGSTROMS",
                       help="keep only triangles whose longest edge is strictly "
                            "below this. Distinct from --contact-cutoff, which "
                            "selects atoms rather than triangles.")
    group.add_argument("--sse", action="store_true", default=None,
                       help="annotate each triangle with its secondary-structure "
                            "category. Does not change a single key: the keys, "
                            "their signs and both fingerprints are identical to "
                            "the same run without it. Protein C-alpha only.")

    if not batch:
        group = parser.add_argument_group("structure")
        group.add_argument("--pdb", metavar="PATH", help="the structure file to read.")
        group.add_argument("--chain", "-c", metavar="ID",
                           help="the chain to run on. Required: no method in this "
                                "package infers a chain.")
        group.add_argument("--chain-2", dest="chain_2", metavar="ID",
                           help="the second chain, for cross methods; for the two protein-protein "
                                "cross methods a different chain from --chain.")
        group.add_argument("--structure-id", dest="structure_id", metavar="ID",
                           help="override the identifier taken from the filename.")
        group.add_argument("--sample-id", dest="sample_id", metavar="NAME",
                           help="the output prefix -- the shared stem of this "
                                "calculation's files in --output, not a path. "
                                "Defaults to a deterministic name built from the "
                                "structure, chain(s), method, target and "
                                "modifiers, e.g. 7LF9_D_carbohydrate_NAG_1.")
        group.add_argument("--label", metavar="TEXT",
                           help="a free-text annotation carried into the "
                                "completion record.")

        group = parser.add_argument_group("target selection")
        group.add_argument("--residue", metavar="ID",
                           help="the residue's identifier in the structure file: "
                                "its number with an optional insertion code, e.g. "
                                "15 or 52A (amino_acid_atomic, with "
                                "--residue-name). Without it, one calculation is "
                                "run per residue of the chain.")
        group.add_argument("--residue-name", dest="residue_name", metavar="NAME",
                           help="three-letter residue name, e.g. LEU "
                                "(amino_acid_atomic). Alone, it runs one "
                                "calculation per residue of that name.")
        group.add_argument("--atoms", metavar="LIST",
                           help="atom selection, e.g. 'LYS_48_NZ;HIS_68_NE2' or "
                                "'N4;C5;O6' depending on the method.")
        group.add_argument("--residues", metavar="LIST",
                           help="residue selection, e.g. '10-12;48;52A' "
                                "(aa_selection: residue numbers with an optional "
                                "insertion code, and ranges of plain numbers, "
                                "ascending; a range spanning an insertion-code "
                                "residue is refused unless that residue is listed "
                                "too) or 'LYS_6;GLU18;GLY_52A' (aa_ind_specific: "
                                "NAME_ID or NAMEID, the identifier being the "
                                "residue number in the file).")
        group.add_argument("--residue-types", dest="residue_types", metavar="LIST",
                           help="amino-acid type(s), e.g. 'TYR' or 'TYR;PHE' "
                                "(aa_ind, aa_ind_mix). The modified residues TPO, "
                                "SEP and PTR are read from HETATM records; any "
                                "other modified residue is named as not supported.")
        group.add_argument("--entity-name", dest="entity_name", metavar="NAME",
                           help="molecule name, e.g. ESW or NAG (drug, "
                                "carbohydrate, nucleotide, and their cross "
                                "methods).")
        group.add_argument("--entity-id", dest="entity_id", metavar="ID",
                           help="the molecule's residue identifier in the "
                                "structure file, e.g. 404 or 52A: that one copy. "
                                "Omit it to use every copy of --entity-name on "
                                "the chain in one combined calculation "
                                "(all_copies; triangles may span copies).")
        group.add_argument("--entity-names", dest="entity_names", metavar="LIST",
                           help="override a chain-scoped method's residue "
                                "vocabulary.")

        group = parser.add_argument_group("cross methods")
        group.add_argument("--contact-cutoff", dest="contact_cutoff", type=float,
                           metavar="ANGSTROMS",
                           help="how close an atom must be to the other side to "
                                "enter the pool. Any positive value is accepted. "
                                "It selects atoms, not triangles: a triangle "
                                "built from admitted atoms may have a much longer "
                                "edge, and every implementation keeps it.")
        group.add_argument("--contact-filter", dest="contact_filter",
                           choices=("both_sides", "side_a_only", "side_b_only"),
                           help="which side(s) the cutoff narrows (default "
                                "both_sides, the standard symmetric rule). Atomic "
                                "cross methods only.")
        group.add_argument("--no-contact-filter", dest="no_contact_filter",
                           action="store_true", default=None,
                           help="skip the contact scan and use both selections "
                                "whole.")
        group.add_argument("--atom-selection", dest="atom_selection", metavar="LIST",
                           help="narrow the protein side to named atoms.")
        group.add_argument("--calpha-residues", dest="calpha_residues", metavar="LIST",
                           help="narrow the protein side to the C-alphas of named "
                                "residues.")
        group.add_argument("--residue-types-2", dest="residue_types_2", metavar="LIST",
                           help="second side: amino-acid type(s).")
        group.add_argument("--residues-2", dest="residues_2", metavar="LIST",
                           help="second side: residue selection.")
        group.add_argument("--atom-selection-2", dest="atom_selection_2",
                           metavar="LIST", help="second side: atom selection.")
        group.add_argument("--calpha-residues-2", dest="calpha_residues_2",
                           metavar="LIST",
                           help="second side: C-alpha residue selection.")


def _add_manifest(parser, manifest_required=True):
    group = parser.add_argument_group("input")
    group.add_argument("--manifest", "-M", required=manifest_required, metavar="CSV",
                       help="the dataset manifest: one row per run.")
    group.add_argument("--pdb-dir", dest="pdb_dir", metavar="DIR",
                       help="where to find structures named by a pdb_id column. "
                            "Not needed when the manifest gives pdb_path.")
    group.add_argument("--contact-cutoff", dest="contact_cutoff", type=float,
                       metavar="ANGSTROMS",
                       help="default contact cutoff for cross rows; a manifest "
                            "column overrides it.")
    group.add_argument("--contact-filter", dest="contact_filter",
                       choices=("both_sides", "side_a_only", "side_b_only"),
                       help="default contact filter for cross rows.")
    group.add_argument("--chain", "-c", metavar="ID",
                       help="the chain for rows that leave it blank or have no "
                            "chain column. A row's own chain always wins; the "
                            "summary records which one applied (chain_source).")

    group = parser.add_argument_group("manifest columns")
    group.add_argument("--column", action="append", metavar="FIELD=HEADER",
                       help="read field FIELD from the column named HEADER, e.g. "
                            "--column label=group or --column pdb_id=PDB. In a "
                            "sheet, for the fields every request shares (pdb_id, "
                            "pdb_path, structure_id, sample_id, label). "
                            "Repeatable.")
    group.add_argument("--keep-column", dest="keep_column", action="append",
                       metavar="HEADER",
                       help="carry the column HEADER through as metadata, into "
                            "the completion record and the batch summary. Never "
                            "used in a calculation. Repeatable.")
    group.add_argument("--ignore-column", dest="ignore_column", action="append",
                       metavar="HEADER",
                       help="drop the column HEADER. Repeatable.")

    group = parser.add_argument_group("wide dataset sheets")
    group.add_argument("--sheet", choices=("tsr", "custom"), default=None,
                       help="read --manifest as a wide sheet: one row per "
                            "structure, methods from --method, each reading its "
                            "own columns. 'tsr': the TSR dataset-sheet columns; "
                            "'custom': only what --map / --map-file give.")
    group.add_argument("--map", action="append", metavar="METHOD[.P].FIELD=HEADER",
                       help="read FIELD of METHOD (participant P, default 1) from "
                            "column HEADER, e.g. --map carbohydrate.chain=glycan, "
                            "--map aa_ind.2.residue_types=types_b, --map "
                            "protein_protein_cross.residue_types_2=aa_ind_2. "
                            "Repeatable; adds to or overrides the layout.")
    group.add_argument("--map-file", dest="map_file", metavar="CSV",
                       help="mappings as a CSV with columns "
                            "method,participant,field,column.")
    group.add_argument("--participant", choices=("1", "2", "all"), default=None,
                       help="which participant's columns one-sided methods read: "
                            "1 (default; aa_chain, aa_ind ...), 2 (aa_chain_2, "
                            "aa_ind_2 ...) or all. A cross method's participant "
                            "1 covers both its sides.")
    group.add_argument("--cross-side", dest="cross_side", action="append",
                       metavar="METHOD.SIDE=COLUMN",
                       help="how a cross method's side is selected when the row "
                            "could narrow it in several ways: by COLUMN, or "
                            "'whole' for the whole side. Sides and columns "
                            "(--sheet tsr): drug_protein_cross.drug "
                            "(drug_atom_selected) and .protein (aa_ind, "
                            "aa_ind_mix, aa_ind_specific, aa_atom_selection); "
                            "carbohydrate_protein_cross.carbohydrate "
                            "(car_atom_selected) and .protein; "
                            "nucleotide_protein_cross.nucleotide "
                            "(nt_atom_selected) and .protein; "
                            "protein_protein_cross.1 (the aa_* columns) and .2 "
                            "(the aa_*_2 columns); protein_protein_calpha_cross.1 "
                            "(aa_selection) and .2 (aa_selection_2). Required "
                            "for every side the row populates. Repeatable.")
    group.add_argument("--inapplicable", choices=("fail", "skip"), default="fail",
                       help="what a row that does not describe a method's "
                            "calculation becomes (a drug row for carbohydrate): "
                            "'fail' (default) records a failure naming the "
                            "missing columns; 'skip' records it as "
                            "'inapplicable'. Either way it is in the summary.")

    group = parser.add_argument_group("structures and execution")
    group.add_argument("--download-missing", dest="download_missing",
                       action="store_true",
                       help="download structures named by pdb_id that are not in "
                            "--pdb-dir before computing (HTTPS, verified; see "
                            "`biotsr fetch`). Off by default: BioTSR never uses "
                            "the network unless asked.")
    _add_download_options(group)
    group.add_argument("--jobs", "-j", metavar="N", default=None,
                       help="compute N requests at a time in worker processes, or "
                            "'auto' for the CPUs this job was allocated (the "
                            "scheduler's allocation, CPU affinity or cgroup "
                            "limit). Default 1: one request at a time in this "
                            "process. Results are identical either way.")


def _add_download_options(group):
    group.add_argument("--download-workers", dest="download_workers", type=int,
                       default=DEFAULT_WORKERS, metavar="N",
                       help="concurrent downloads, independent of --jobs "
                            "(default %d, from 1 to %d)." % (DEFAULT_WORKERS, MAX_WORKERS))
    group.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT,
                       metavar="SECONDS",
                       help="how long to wait for a connection, or for the next "
                            "piece of a response, before that attempt counts as "
                            "failed (default %g). It limits each wait, not a "
                            "whole download: a server that keeps sending, "
                            "however slowly, is never timed out."
                            % DEFAULT_TIMEOUT)
    group.add_argument("--retries", type=int, default=DEFAULT_RETRIES, metavar="N",
                       help="further attempts after a transient network error, "
                            "with exponential backoff of up to 60 s between "
                            "attempts (default %d)." % DEFAULT_RETRIES)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="biotsr",
        description="BioTSR -- the unified software implementation of the Triangular Spatial Relationship (TSR) methodology.",
        epilog=_EPILOG_MAIN,
        formatter_class=_Formatter,
    )
    from .. import __version__

    parser.add_argument(
        "--version", action="version", version="biotsr %s" % __version__
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    methods = sub.add_parser(
        "methods",
        help="list the supported TSR methods",
        description="List the supported TSR methods, their target arguments, "
                    "modifiers and standard profile. Generated from the method registry, so it cannot "
                    "disagree with what the package will actually run.",
        formatter_class=_Formatter,
    )
    methods.add_argument("method", nargs="?", metavar="NAME",
                         help="describe one method in full.")
    methods.add_argument("--json", action="store_true",
                         help="machine-readable output on stdout.")
    methods.set_defaults(command="methods")

    run = sub.add_parser(
        "run",
        help="run one TSR calculation on one structure",
        description="Run one TSR calculation and write its keys, triangles and "
                    "provenance into an output directory.",
        epilog=_EPILOG_RUN,
        formatter_class=_Formatter,
    )
    _add_science(run, batch=False)
    _add_common(run)
    run.set_defaults(command="run")

    batch = sub.add_parser(
        "batch",
        help="run every row of a dataset manifest or sheet",
        description="Run a whole dataset from one manifest or wide sheet, writing "
                    "every calculation into one shared output folder plus a "
                    "summary and, if anything failed, a failure report.",
        epilog=_EPILOG_BATCH,
        formatter_class=_Formatter,
    )
    _add_manifest(batch)
    _add_science(batch, batch=True)
    _add_common(batch)
    batch.add_argument("--on-error", dest="on_error",
                       choices=(ON_ERROR_CONTINUE, ON_ERROR_STOP),
                       default=ON_ERROR_CONTINUE,
                       help="'continue' (default) runs every request and records "
                            "each failure against its source row; 'stop' halts at "
                            "the first failure. Either way no failure is omitted.")
    batch.add_argument("--dry-run", dest="dry_run", action="store_true",
                       help="plan the batch and print every request -- source "
                            "row, method, participant, chains, selection, output "
                            "prefix and outcome -- plus unused columns, "
                            "structures present or to be downloaded and existing "
                            "outputs, without downloading, computing or writing "
                            "anything.")
    batch.set_defaults(command="batch")

    fetch = sub.add_parser(
        "fetch",
        help="download PDB-format structures from RCSB",
        description="Download PDB-format structures into a directory that "
                    "`biotsr batch --pdb-dir` reads, verifying each one, and "
                    "report every id.",
        epilog=_EPILOG_FETCH,
        formatter_class=_Formatter,
    )
    fetch.add_argument("ids", nargs="*", metavar="PDB_ID",
                       help="PDB ids, e.g. 1UBQ 4HHB (case-insensitive).")
    fetch.add_argument("--id-file", dest="id_file", metavar="FILE",
                       help="a text file of ids, whitespace- or comma-separated; "
                            "'#' starts a comment.")
    fetch.add_argument("--manifest", "-M", metavar="CSV",
                       help="read ids from one column of a CSV: pdb_id or an alias "
                            "such as protein. Other columns are not read or "
                            "checked, so any dataset sheet works.")
    fetch.add_argument("--column", action="append", metavar="pdb_id=HEADER",
                       help="with --manifest, the column that holds the ids.")
    fetch.add_argument("--pdb-dir", dest="pdb_dir", required=True, metavar="DIR",
                       help="where to put the files; also the cache that is "
                            "checked first. Required.")
    group = fetch.add_argument_group("downloading")
    _add_download_options(group)
    group.add_argument("--refresh", action="store_true",
                       help="download again even when a valid file is cached; "
                            "the cached file is replaced only after the new one "
                            "passes every check.")
    group.add_argument("--dry-run", dest="dry_run", action="store_true",
                       help="check the cache and report what would be "
                            "downloaded; no request is made and nothing is "
                            "written.")
    group.add_argument("--report", metavar="CSV", default=None,
                       help="where to write the report (default "
                            "DIR/fetch_report.csv).")
    group = fetch.add_argument_group("reporting")
    group.add_argument("--quiet", "-q", action="store_true", help="print errors only.")
    group.add_argument("--verbose", "-v", action="store_true",
                       help="print one line per id.")
    group.add_argument("--debug", action="store_true",
                       help="print a traceback for unexpected internal errors.")
    fetch.set_defaults(command="fetch")

    analyze = sub.add_parser(
        "analyze",
        help="analysis over a dataset",
        description="Analysis operations over a manifest of samples.",
        formatter_class=_Formatter,
    )
    analysis_sub = analyze.add_subparsers(dest="analysis", metavar="OPERATION")
    common = analysis_sub.add_parser(
        "common-keys",
        help="the keys every sample shares",
        description="Report the keys present in every sample of the manifest, "
                    "with each sample's total key count and the frequency it "
                    "contributes to the shared set.",
        epilog=_EPILOG_ANALYZE,
        formatter_class=_Formatter,
    )
    # `--manifest` stops being required here because `--from-results` is the
    # other way to name a dataset. Exactly one of the two must be given, which
    # `cmd_analyze` enforces with a message naming both -- argparse's own
    # "required" wording cannot express "one of these two".
    _add_manifest(common, manifest_required=False)
    common.add_argument(
        "--from-results", dest="from_results", metavar="DIR", default=None,
        help="analyse previously saved results instead of recomputing them: DIR "
             "is searched recursively for completed calculations (PREFIX.run.json "
             "completion records), so point it at an --output folder a previous "
             "run wrote. Each sample is named by its prefix. Mutually exclusive "
             "with --manifest. --method NAME (repeatable) keeps only the "
             "calculations of those methods, so one shared folder can be "
             "analysed one key space at a time. The structures are never "
             "re-read, so no --pdb-dir "
             "is needed and no scientific flag applies -- every sample keeps the "
             "key space and settings it was computed under, and comparing them "
             "is checked exactly as for a freshly computed set.")
    common.add_argument(
        "--allow-mixed-processing", dest="allow_mixed_processing",
        action="store_true",
        help="compare samples whose processing settings differ (for example one "
             "with a size filter and one without, or two contact cutoffs). Off "
             "by default, when such a set is refused. Never relaxes the key "
             "definition check. The report records every sample's settings and "
             "a warning.")
    common.add_argument(
        "--allow-partial", dest="allow_partial", action="store_true",
        help="analyse samples that are partial: their selection did not get "
             "everything it asked for (a requested residue, atom or molecule "
             "absent, a modified residue the method cannot read, one the "
             "completeness rule rejected). Off by default, when such a set is "
             "refused naming each partial sample and what it lacks. Separate "
             "from --allow-mixed-processing; never relaxes the key definition "
             "check. The report records each sample's omissions, and the exit "
             "code is 4 when any analysed sample is partial.")
    common.add_argument("--allow-mixed-sampling", dest="renamed_allow_mixed_sampling",
                        action="store_true", help=argparse.SUPPRESS)
    _add_science(common, batch=True)
    _add_common(common, triplets=False)
    common.add_argument("--group-by", dest="group_by", choices=("label",),
                        default=None,
                        help="report one intersection per value of the label "
                             "(a manifest's label column, or aa_group in a TSR "
                             "sheet) instead of one over everything.")
    common.add_argument("--on-error", dest="on_error",
                        choices=(ON_ERROR_CONTINUE, ON_ERROR_STOP),
                        default=ON_ERROR_STOP,
                        help="'stop' (default here) halts at the first failure: a "
                             "key intersection silently computed over a subset of "
                             "the dataset is a wrong answer, not a partial one.")
    common.set_defaults(command="analyze", analysis="common-keys")
    analyze.set_defaults(command="analyze", analysis=None)

    describe = sub.add_parser(
        "describe",
        help="print the provenance of an existing output",
        description="Read an output folder (one line per calculation), a "
                    "calculation's prefix or PREFIX.run.json (how it was "
                    "produced), or a .meta.json sidecar (what it means).",
        formatter_class=_Formatter,
    )
    describe.add_argument("path", metavar="PATH",
                          help="an output folder, a prefix (results/1UBQ_A_protein), "
                               "a .run.json, a .meta.json or a .keys file.")
    describe.add_argument("--json", action="store_true",
                          help="print the file itself, formatted.")
    describe.set_defaults(command="describe")

    return parser
