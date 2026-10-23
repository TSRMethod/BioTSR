"""Every example in ``docs/API_GUIDE.md``, executed.

The guide promises that a user needs nothing but it for ordinary use.  That
promise decays silently unless something runs the examples, so this module runs
each numbered section and asserts what the guide claims about it -- including the
claims that are easy to let rot, like "the same integer keys" and "byte-identical
reports".

Section numbers here match section numbers there.
"""

from __future__ import annotations

import glob
import os
import re

import pytest

import biotsr

from ..conftest import _require

#: ``tests/integration/`` -> ``tests/`` -> the repository root, whose ``docs/``
#: holds the guide beside the other user-facing documents.
GUIDE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "docs",
    "API_GUIDE.md",
)


# ---------------------------------------------------------------------------
# 1. Run a method from Python
# ---------------------------------------------------------------------------


def test_the_one_line_call(ubq_pdb):
    result = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    assert result.n_triangles > 0
    assert result.method == "protein"


def test_the_three_modifiers_compose_on_one_call(ubq_pdb):
    result = biotsr.analyze(
        ubq_pdb, method="protein", chain="A",
        aa_grouping=True, mirror=True, size_filter=25.0, keep_triplets=False,
    )
    assert result.key_space.label_scheme.name == "amino_acid_group_23"
    assert result.key_space.label_scheme.radix == 30
    assert any(k < 0 for k in result.keys)
    assert result.sampling.size_filter == 25.0


def test_sse_produces_the_same_integer_keys(p53_path):
    """The guide's central SSE claim, asserted exactly as it is written."""
    annotated = biotsr.analyze(p53_path, method="protein", chain="A", sse=True,
                            keep_triplets=False)
    plain = biotsr.analyze(p53_path, method="protein", chain="A", keep_triplets=False)
    assert annotated.frequencies == plain.frequencies
    assert annotated.key_space.fingerprint == plain.key_space.fingerprint
    assert annotated.sampling.fingerprint == plain.sampling.fingerprint
    # ... and the categories are what changed.
    assert set(annotated.category_totals()) != set(plain.category_totals())


def test_the_targeted_methods(ubq_pdb, trypsin_path):
    leu15 = biotsr.analyze(ubq_pdb, method="amino_acid_atomic", chain="A",
                        residue="15", residue_name="LEU", keep_triplets=False)
    assert leu15.n_triangles > 0
    ligand = biotsr.analyze(trypsin_path, method="drug", chain="A",
                         entity_name="BEN", entity_id="1", keep_triplets=False)
    assert ligand.n_triangles > 0


def test_the_cross_entry_point(barnase_path):
    cross = biotsr.analyze_cross(
        barnase_path, method="protein_protein_cross",
        chain="A", chain_2="D", contact_cutoff=3.2, keep_triplets=False,
    )
    assert cross.chain == "A|D"
    assert cross.n_triangles > 0


def test_the_method_listings():
    assert len(biotsr.list_methods()) == 20
    assert len(biotsr.list_cross_methods()) == 5


def test_an_unsupported_argument_raises_naming_what_is_supported(trypsin_path):
    with pytest.raises(ValueError) as exc:
        biotsr.analyze(trypsin_path, method="drug", chain="A", entity_name="BEN",
                    entity_id="1", sse=True)
    message = str(exc.value)
    assert "does not support" in message and "it supports" in message


def test_the_streaming_example(ubq_pdb, tmp_path):
    from biotsr.io.native import NativeTripletStreamWriter

    path = str(tmp_path / "1UBQ_A.triplets_29_35")
    with NativeTripletStreamWriter(path) as sink:
        result = biotsr.analyze(ubq_pdb, method="protein", chain="A",
                             keep_triplets=False, triplet_sink=sink.write)
    assert result.triplets is None
    assert sink.n_written == result.n_triangles
    # "otherwise identical": the keys and the provenance are unaffected.
    kept = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)
    assert result.frequencies == kept.frequencies


# ---------------------------------------------------------------------------
# 2. Inspect a TSRResult
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def result(ubq_pdb):
    return biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False)


def test_every_attribute_the_guide_lists_exists(result):
    assert isinstance(result.frequencies, dict)
    assert isinstance(result.keys, dict)
    assert isinstance(result.n_triangles, int)
    assert isinstance(result.n_unique_keys, int)
    assert result.category_totals() == {"regular": result.n_triangles}
    assert result.triplets is None  # keep_triplets=False
    assert isinstance(list(result.warnings), list)
    assert isinstance(list(result.skipped), list)


def test_the_ledger_balances(result):
    c = result.counts
    assert c.n_candidates == c.n_degenerate + c.n_filtered_by_size + c.n_triangles
    assert c.ledger_balances


def test_a_size_filter_changes_sampling_and_not_the_key_space(ubq_pdb, result):
    """The guide's example of what the two identities separate."""
    filtered = biotsr.analyze(ubq_pdb, method="protein", chain="A", size_filter=14.0,
                           keep_triplets=False)
    assert filtered.key_space.fingerprint == result.key_space.fingerprint
    assert filtered.sampling.fingerprint != result.sampling.fingerprint


def test_grouping_changes_the_key_space_at_the_same_radix(ubq_pdb, result):
    grouped = biotsr.analyze(ubq_pdb, method="protein", chain="A", aa_grouping=True,
                          keep_triplets=False)
    assert grouped.key_space.fingerprint != result.key_space.fingerprint
    assert grouped.key_space.label_scheme.radix == 30
    assert result.key_space.label_scheme.radix == 30
    assert grouped.key_space.label_scheme.name != result.key_space.label_scheme.name
    assert (
        grouped.key_space.label_scheme.content_fingerprint
        != result.key_space.label_scheme.content_fingerprint
    )


# ---------------------------------------------------------------------------
# 3. Save a result
# ---------------------------------------------------------------------------


def test_save_writes_the_three_files_the_guide_names(ubq_pdb, tmp_path):
    from biotsr.io import save

    full = biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=True)
    written = save(full, str(tmp_path / "results"))
    assert set(written) == {"keys", "metadata", "triplets"}
    assert os.path.basename(written["keys"]) == "1UBQ_A_protein.keys"
    assert os.path.basename(written["metadata"]) == "1UBQ_A_protein.meta.json"
    assert os.path.basename(written["triplets"]) == "1UBQ_A_protein.triplets"
    named = save(full, str(tmp_path / "results"), prefix="ubiquitin", triplets=False)
    assert os.path.basename(named["keys"]) == "ubiquitin.keys"


def test_two_runs_of_the_same_request_are_byte_identical(result, tmp_path):
    from biotsr.io import save

    a = save(result, str(tmp_path / "a"), triplets=False)
    b = save(result, str(tmp_path / "b"), triplets=False)
    for role in ("keys", "metadata"):
        with open(a[role], "rb") as fh:
            first = fh.read()
        with open(b[role], "rb") as fh:
            second = fh.read()
        assert first == second


def test_the_keys_file_has_the_columns_and_the_order_the_guide_states(
    result, tmp_path
):
    from biotsr.io import save

    written = save(result, str(tmp_path / "cols"), triplets=False)
    with open(written["keys"]) as fh:
        lines = fh.read().splitlines()
    assert lines[0] == "key\tcategory\tfrequency"
    rows = [line.split("\t") for line in lines[1:]]
    seen = [(int(k), c) for k, c, _f in rows]
    assert seen == sorted(seen)


# ---------------------------------------------------------------------------
# 4. Load the saved result
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def saved(result, tmp_path_factory):
    from biotsr.io import save

    out = str(tmp_path_factory.mktemp("guide_saved") / "1UBQ_A")
    return save(result, out, triplets=False)


def test_either_file_names_the_pair(saved):
    from_keys = biotsr.load_result(saved["keys"])
    from_sidecar = biotsr.load_result(saved["metadata"])
    assert from_keys.frequencies == from_sidecar.frequencies


def test_the_four_equalities_the_guide_asserts(result, saved):
    loaded = biotsr.load_result(saved["keys"])
    assert loaded.frequencies == result.frequencies
    assert loaded.key_space.fingerprint == result.key_space.fingerprint
    assert loaded.sampling.fingerprint == result.sampling.fingerprint
    assert loaded.counts == result.counts


def test_the_two_documented_absences(saved):
    loaded = biotsr.load_result(saved["keys"])
    assert loaded.triplets is None
    assert loaded.key_order == ()
    marker = loaded.annotations["loaded_from"]
    assert set(marker) == {
        "keys", "metadata", "result_schema_version",
        "triplets_loaded", "key_order_loaded",
    }
    assert marker["result_schema_version"] == 2
    assert marker["triplets_loaded"] is False
    assert marker["key_order_loaded"] is False


def test_the_two_exception_types_the_guide_names(result, tmp_path):
    import json

    from biotsr.io import save

    written = save(result, str(tmp_path / "bad"), triplets=False)
    with open(written["metadata"]) as fh:
        payload = json.load(fh)
    payload["result_schema_version"] = 99
    with open(written["metadata"], "w") as fh:
        json.dump(payload, fh)
    with pytest.raises(biotsr.UnsupportedFormat):
        biotsr.load_result(written["keys"])

    other = biotsr.analyze(result.source.path, method="protein", chain="A",
                        size_filter=12.0, keep_triplets=False)
    w2 = save(other, str(tmp_path / "other"), triplets=False)
    w1 = save(result, str(tmp_path / "mine"), triplets=False)
    with pytest.raises(biotsr.MalformedResult):
        biotsr.load_result(w1["keys"], metadata_path=w2["metadata"])


# ---------------------------------------------------------------------------
# 5. Verify provenance
# ---------------------------------------------------------------------------


def test_every_provenance_attribute_the_guide_lists(saved):
    loaded = biotsr.load_result(saved["keys"])
    assert loaded.method == "protein"
    assert loaded.annotations["profile"] == "standard_protein_calpha_v1"
    assert loaded.annotations["base_profile"] == "standard_protein_calpha_v1"
    assert loaded.annotations["legacy_compatibility_profile"] is False
    assert loaded.annotations["modifiers"] == ()
    assert loaded.source.path
    assert len(loaded.source.sha256) == 64
    assert loaded.source.parser_backend
    assert loaded.selection.chains_used == ("A",)
    assert loaded.key_space.numerical_precision == "float64"


def test_the_sse_provenance_example(p53_path, tmp_path):
    from biotsr.io import save
    from biotsr.sse import (
        category_totals_in_order,
        check_conservation,
        strand_identity_of,
    )

    annotated = biotsr.analyze(p53_path, method="protein", chain="A", sse=True,
                            keep_triplets=False)
    written = save(annotated, str(tmp_path / "sse"), triplets=False)
    loaded = biotsr.load_result(written["keys"])
    assert strand_identity_of(loaded) == "sheet_id_and_serial"
    assert loaded.annotations["sse"]["has_records"] is True
    totals = dict(category_totals_in_order(loaded))
    assert len(totals) == 18
    check_conservation(loaded)


def test_the_strand_identity_guard_example(p53_path):
    from biotsr.sse import (
        STRAND_IDENTITY_SHEET_ID_AND_SERIAL,
        IncompatibleSSEStratification,
        require_common_strand_identity,
    )

    a = biotsr.analyze(p53_path, method="protein", chain="A", sse=True,
                    keep_triplets=False)
    b = biotsr.analyze(p53_path, method="protein", chain="A", keep_triplets=False)
    assert require_common_strand_identity([a, a]) == STRAND_IDENTITY_SHEET_ID_AND_SERIAL
    # a result with no SSE annotation has no categories to compare
    with pytest.raises(IncompatibleSSEStratification):
        require_common_strand_identity([a, b])


# ---------------------------------------------------------------------------
# 6. Build a feature representation
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def samples(ubq_pdb, trypsin_path):
    from biotsr.analysis import AnalysisSample

    return [
        AnalysisSample.from_result(
            biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
        )
        for path in (ubq_pdb, trypsin_path)
    ]


def test_the_feature_matrix_attributes_the_guide_lists(samples):
    from biotsr.analysis import feature_matrix

    fm = feature_matrix(samples)
    assert fm.sample_ids == ("1UBQ_A", "3PTB_A")
    assert fm.columns == tuple(sorted(fm.columns))
    assert fm.index_of("1UBQ_A") == 0
    assert len(fm.row(0)) == len(fm.columns)
    assert fm.value(0, fm.columns[0]) >= 0
    assert fm.comparison["processing_uniform"] is True
    assert [p["sample_id"] for p in fm.samples] == ["1UBQ_A", "3PTB_A"]
    assert fm.to_header_csv().startswith("sample_id,")


def test_an_incompatible_key_space_is_refused_with_the_field_named(
    ubq_pdb, trypsin_path
):
    from biotsr.analysis import AnalysisSample, feature_matrix

    grouped = AnalysisSample.from_result(
        biotsr.analyze(ubq_pdb, method="protein", chain="A", aa_grouping=True,
                    keep_triplets=False),
        sample_id="grouped",
    )
    plain = AnalysisSample.from_result(
        biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False),
        sample_id="plain",
    )
    with pytest.raises(biotsr.IncompatibleKeySpace) as exc:
        feature_matrix([plain, grouped])
    assert "label_scheme" in str(exc.value)


def test_a_sampling_difference_is_refused_unless_allowed(ubq_pdb, trypsin_path):
    from biotsr.analysis import AnalysisSample, feature_matrix

    a = AnalysisSample.from_result(
        biotsr.analyze(ubq_pdb, method="protein", chain="A", keep_triplets=False),
        sample_id="unfiltered",
    )
    b = AnalysisSample.from_result(
        biotsr.analyze(ubq_pdb, method="protein", chain="A", size_filter=14.0,
                    keep_triplets=False),
        sample_id="filtered",
    )
    with pytest.raises(biotsr.IncompatibleSampling, match="size_filter"):
        feature_matrix([a, b])
    fm = feature_matrix([a, b], allow_mixed_processing=True)
    assert fm.comparison["processing_uniform"] is False
    assert any("size_filter" in w for w in fm.warnings)


# ---------------------------------------------------------------------------
# 7. Run an analysis
# ---------------------------------------------------------------------------


def test_the_common_key_example(samples):
    from biotsr.analysis import common_key_intersection

    intersection = common_key_intersection(samples)
    assert intersection.keys == tuple(sorted(intersection.keys))
    assert intersection.n_common_keys == len(intersection.keys)
    assert len(intersection.rows) == 2


def test_the_similarity_and_distance_example(samples):
    from biotsr.analysis import (
        MEASURES,
        distance_matrix,
        feature_matrix,
        generalized_jaccard,
        similarity_matrix,
    )

    assert 0.0 <= generalized_jaccard(samples[0].counts, samples[1].counts) <= 1.0
    fm = feature_matrix(samples)
    sm = similarity_matrix(fm)
    dm = distance_matrix(fm)
    # Condensed storage: one value for two samples, and the diagonal is a
    # definition rather than a computed number.
    assert len(sm.values) == 1
    assert sm.value(0, 1) == pytest.approx(1.0 - dm.value(0, 1))
    assert sm.square()[0][0] == 1.0
    assert dm.square()[0][0] == 0.0
    assert sorted(MEASURES) == [
        "binary_jaccard", "cosine", "generalized_jaccard", "sarika", "wu",
    ]


def test_the_clustering_example(ubq_pdb, trypsin_path, barnase_path):
    from biotsr.analysis import (
        AnalysisSample,
        clustermap_spec,
        distance_matrix,
        feature_matrix,
    )

    cases = [(ubq_pdb, "A", "small"), (trypsin_path, "A", "large"),
             (barnase_path, "A", "large"), (barnase_path, "B", "large")]
    built = [
        AnalysisSample.from_result(
            biotsr.analyze(path, method="protein", chain=chain, keep_triplets=False),
            sample_id="%d" % i,
            label=label,
        )
        for i, (path, chain, label) in enumerate(cases)
    ]
    spec = clustermap_spec(distance_matrix(feature_matrix(built)),
                           [s.label for s in built], n_clusters=2)
    assert spec.assignment.n_clusters_requested == 2
    assert len(spec.assignment.labels) == 4
    assert sorted(spec.assignment.leaves) == [0, 1, 2, 3]
    assert -1.0 <= spec.adjusted_rand <= 1.0


def test_label_from_name_takes_the_token_index():
    from biotsr.analysis import label_from_name

    assert label_from_name("1ABC_A_15_LEU_SITE", 3) == "LEU"
    assert label_from_name("1ABC_A_15_LEU_SITE", -1) == "SITE"


# ---------------------------------------------------------------------------
# "Compute once, analyse many times"
# ---------------------------------------------------------------------------


def test_the_compute_once_analyse_many_workflow(ubq_pdb, trypsin_path, tmp_path):
    from biotsr.analysis import AnalysisSample, common_key_intersection
    from biotsr.io import save

    root = tmp_path / "results"
    for path in (ubq_pdb, trypsin_path):
        r = biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
        save(r, str(root), triplets=False)          # one folder, one prefix each

    loaded = [
        AnalysisSample.from_result(biotsr.load_result(p))
        for p in sorted(glob.glob(str(root / "*.keys")))
    ]
    assert len(loaded) == 2
    from_disk = common_key_intersection(loaded)

    live = [
        AnalysisSample.from_result(
            biotsr.analyze(path, method="protein", chain="A", keep_triplets=False)
        )
        for path in (ubq_pdb, trypsin_path)
    ]
    assert from_disk.keys == common_key_intersection(live).keys


def test_the_cli_equivalent_of_that_workflow(ubq_pdb, trypsin_path, tmp_path):
    from biotsr.cli import EXIT_OK, main

    results = str(tmp_path / "results")
    for path in (ubq_pdb, trypsin_path):
        assert main(["run", "-m", "protein", "--pdb", path, "-c", "A",
                     "-o", results, "--no-triplets", "-q"]) == EXIT_OK
    assert main(["analyze", "common-keys", "--from-results", results,
                 "-o", str(tmp_path / "analysis"), "-q"]) == EXIT_OK
    assert os.path.isfile(
        str(tmp_path / "analysis" / "analysis" / "common_keys" / "all_samples.csv")
    )


# ---------------------------------------------------------------------------
# Errors, and the public surface
# ---------------------------------------------------------------------------


def test_the_error_tree_the_guide_draws():
    assert issubclass(biotsr.InvalidRequest, biotsr.TSRError)
    assert issubclass(biotsr.IncompatibleRepresentation, biotsr.TSRError)
    assert issubclass(biotsr.MalformedResult, biotsr.TSRError)
    assert issubclass(biotsr.UnsupportedFormat, biotsr.TSRError)
    for kind in (biotsr.IncompatibleKeySpace, biotsr.IncompatibleSampling,
                 biotsr.IncompatibleModifier, biotsr.IncompatibleSSEStratification):
        assert issubclass(kind, biotsr.IncompatibleRepresentation)
    assert issubclass(biotsr.TSRError, ValueError)


def test_every_import_path_in_the_public_surface_table_works():
    import biotsr.analysis
    import biotsr.cli
    import biotsr.io
    import biotsr.profiles
    import biotsr.sse

    assert biotsr.io.save and biotsr.io.load_result
    assert biotsr.profiles.Profile and biotsr.profiles.STANDARD_PROTEIN_CALPHA_V1
    assert biotsr.cli.main
    assert biotsr.analysis.feature_matrix
    assert biotsr.sse.require_common_strand_identity


def test_importing_the_analysis_layer_does_not_import_the_renderer():
    """The guide's claim that matplotlib is not required."""
    import subprocess
    import sys

    code = (
        "import sys, biotsr.analysis; "
        "assert 'biotsr.analysis.render' not in sys.modules; "
        "assert 'matplotlib' not in sys.modules"
    )
    src = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    env = dict(os.environ, PYTHONPATH=os.path.join(src, "src"))
    assert subprocess.call([sys.executable, "-c", code], env=env) == 0


# ---------------------------------------------------------------------------
# the guide itself
# ---------------------------------------------------------------------------


def test_the_guide_exists_and_has_the_seven_sections():
    with open(_require(GUIDE)) as fh:
        text = fh.read()
    for n, title in enumerate(
        [
            "Run a method from Python",
            "Inspect a `TSRResult`",
            "Save a result",
            "Load the saved result",
            "Verify provenance",
            "Build a feature representation",
            "Run an analysis",
        ],
        start=1,
    ):
        assert "## %d. %s" % (n, title) in text, title


def test_the_guide_names_no_method_the_registry_does_not_have():
    """A guide that documents a method which was renamed is worse than none."""
    with open(_require(GUIDE)) as fh:
        text = fh.read()
    registered = set(biotsr.list_methods())
    for quoted in set(re.findall(r'method="([a-z_]+)"', text)):
        assert quoted in registered, quoted
