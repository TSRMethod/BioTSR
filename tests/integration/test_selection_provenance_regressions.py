"""Smaller corrections made in the release candidates, kept as regressions.

* **altLoc provenance.**  rc2 recorded ``blank_or_a`` for every method while
  ``amino_acid_atomic`` admitted every altLoc except ``B``.  Since 1.0.0rc4
  every method applies the one rule, blank or ``A`` (ADR-095), and records it.
* **A target the completeness rule declines** was refused on the command line
  as "selection produced no atoms" with a list of the structure's chains, which
  sends the user hunting for a typo that is not there.
* **A completion record of another version** is never taken for this one's.
* **A cross with no atom in contact** is refused naming the cutoff.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

import biotsr
from biotsr.selection.base import AltLocPolicy

from ..conftest import SRC, structure

#: One target per method, enough to produce a result.
_TARGETS = {
    "protein": ("1UBQ", "A", {}),
    "aa_selection": ("1UBQ", "A", {"residues": "1-10"}),
    "amino_acid_atomic": ("1UBQ", "A", {"residue": "15", "residue_name": "LEU"}),
    "amino_acid_atom_selection": ("1UBQ", "A", {"atoms": "LEU_8_CD1;LEU_43_CD2;LYS_48_NZ"}),
    "aa_ind": ("1UBQ", "A", {"residue_types": "LYS"}),
    "aa_ind_mix": ("1UBQ", "A", {"residue_types": "LYS;ARG"}),
    "aa_ind_specific": ("1UBQ", "A", {"residues": "LYS_6;LYS_11"}),
    "drug": ("3PTB", "A", {"entity_name": "BEN", "entity_id": "1"}),
    "drug_atom_selection": ("3PTB", "A", {"entity_name": "BEN", "entity_id": "1",
                                          "atoms": "C1;C2;C3"}),
    "carbohydrate": ("5WT9", "A", {"entity_name": "NAG", "entity_id": "1"}),
    "carbohydrate_chain": ("5WT9", "A", {}),
    "carbohydrate_atom_selection": ("7LF9", "D", {"entity_name": "NAG", "entity_id": "1",
                                                  "atoms": "C1;C2;N2"}),
    "nucleotide": ("1BNA", "A", {"entity_name": "DG", "entity_id": "2"}),
    "nucleotide_chain": ("1BNA", "A", {}),
    "nucleotide_atom_selection": ("2AC0", "E", {"entity_name": "DT", "entity_id": "7",
                                                "atoms": "P;C1';N1"}),
}


@pytest.mark.parametrize("method", sorted(_TARGETS))
def test_every_method_records_the_one_altloc_rule_it_applies(method):
    pdb, chain, target = _TARGETS[method]
    result = biotsr.analyze(structure(pdb), method=method, chain=chain,
                            keep_triplets=False, **target)
    recipe = biotsr.get_method(method)
    assert result.source.altloc_policy == AltLocPolicy.BLANK_OR_A
    admission = recipe.default_profile.atom_admission
    if admission is not None:
        assert tuple(admission.accepted_alt_locs) == ("", "A")
    assert not hasattr(AltLocPolicy, "NOT_B")


def test_every_one_sided_method_is_covered():
    one_sided = {m for m in biotsr.list_methods()
                 if biotsr.get_method(m).family != "cross"}
    assert one_sided == set(_TARGETS)


def _three_conformers(tmp_path):
    """1UBQ with LEU 15 CD1 given altLocs A, B and C."""
    out = []
    with open(structure("1UBQ")) as fh:
        for line in fh.read().splitlines():
            if (line.startswith("ATOM") and line[21] == "A"
                    and line[22:26].strip() == "15" and line[12:16].strip() == "CD1"):
                for alt, occ in (("A", "0.40"), ("B", "0.30"), ("C", "0.30")):
                    out.append(line[:16] + alt + line[17:54] + "%6s" % occ + line[60:])
            else:
                out.append(line)
    path = tmp_path / "1UBQ.pdb"
    path.write_text("\n".join(out) + "\n")
    return str(path)


def test_a_third_conformer_is_excluded_and_the_residue_kept(tmp_path):
    """Until 1.0.0rc3 amino_acid_atomic admitted altLoc C as a ninth LEU atom and
    the exact-atom-count rule then excluded the residue. Now B and C are both
    excluded, so LEU 15 is computed from blank + A, exactly as the pooled
    selectors and the unmodified structure compute it."""
    path = _three_conformers(tmp_path)
    result = biotsr.analyze(path, method="amino_acid_atomic", chain="A",
                            residue="15", residue_name="LEU", keep_triplets=False)
    original = biotsr.analyze(structure("1UBQ"), method="amino_acid_atomic", chain="A",
                              residue="15", residue_name="LEU", keep_triplets=False)
    assert result.counts.n_points == 8
    assert result.frequencies == original.frequencies
    assert result.source.altloc_policy == "blank_or_a"
    assert any("excluded 2 atom record(s) with indicator B, C" in w
               for w in result.warnings)
    pooled = biotsr.analyze(path, method="aa_ind_specific", chain="A",
                            residues="LEU_15", keep_triplets=False)
    assert pooled.counts.n_points == 8 and pooled.source.altloc_policy == "blank_or_a"


def test_a_declined_target_is_refused_naming_the_completeness_rule(tmp_path):
    env = dict(os.environ, PYTHONPATH=SRC)
    done = subprocess.run(
        [sys.executable, "-m", "biotsr", "run", "--method", "nucleotide", "--pdb",
         structure("1BNA"), "--chain", "A", "--entity-name", "DC", "--entity-id", "1",
         "--output", str(tmp_path / "out"), "--no-triplets"],
        env=env, capture_output=True, text=True, timeout=120)
    assert done.returncode == 2
    assert "completeness rule excluded the selection" in done.stderr
    assert "expected exactly 19 atoms, got 16" in done.stderr
    assert "chains in this structure" not in done.stderr
    assert not (tmp_path / "out").exists()


def test_a_completion_record_of_another_version_is_never_taken_for_this_ones(tmp_path):
    """A record whose schema is not this version's is neither reused by
    --if-exists skip nor replaced by --if-exists overwrite: BioTSR cannot
    know which files it owns, so it names the version and changes nothing."""
    import json

    env = dict(os.environ, PYTHONPATH=SRC)
    manifest = tmp_path / "m.csv"
    manifest.write_text("sample_id,pdb_path,method,chain,residue_types\n"
                        "lys,%s,aa_ind,A,LYS\n" % structure("1UBQ"))
    out = tmp_path / "out"

    def batch(*extra):
        return subprocess.run(
            [sys.executable, "-m", "biotsr", "batch", "--manifest", str(manifest),
             "--output", str(out), "--no-triplets"] + list(extra),
            env=env, capture_output=True, text=True, timeout=120)

    assert batch().returncode == 0
    assert batch("--if-exists", "skip").returncode == 0      # same version: reused
    record = out / "lys.run.json"
    payload = json.loads(record.read_text())
    payload["schema"] = 99                                    # another version's
    record.write_text(json.dumps(payload))
    for policy in ("skip", "overwrite"):
        again = batch("--if-exists", policy)
        assert again.returncode == 3
        assert "different BioTSR version (schema 99" in again.stderr
    assert json.loads(record.read_text())["schema"] == 99


def test_a_cross_with_no_contact_is_refused_naming_the_cutoff(tmp_path):
    """Both sides select atoms (BEN's ring carbons, the protein chain), but none
    lies within 3.5 A of the other side: the pool, not the selection, is empty."""
    env = dict(os.environ, PYTHONPATH=SRC)
    done = subprocess.run(
        [sys.executable, "-m", "biotsr", "run", "--method", "drug_protein_cross",
         "--pdb", structure("3PTB"), "--chain", "A", "--chain-2", "A",
         "--entity-name", "BEN", "--atoms", "C1;C2;C3", "--output",
         str(tmp_path / "out"), "--no-triplets"],
        env=env, capture_output=True, text=True, timeout=120)
    assert done.returncode == 2
    assert "no atom pair across the two sides lies within 3.5 A" in done.stderr
    assert "--contact-cutoff" in done.stderr
    assert "selection produced no atoms" not in done.stderr
    assert not (tmp_path / "out").exists()
