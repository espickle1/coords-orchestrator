"""The build: what is written, what is flagged, and what is refused."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from coords.build import build_all
from coords.config import ConfigError, load_config
from coords.testing import make_toy_project

OUT_FILES = ["reference_maps.csv", "sequences.csv", "sequence_maps.csv", "structures.csv",
             "structure_chains.csv", "structure_maps.csv", "sites_resolved.csv", "manifest.json"]


def _csv(toy, name, **kw):
    return pd.read_csv(toy.maps / name, keep_default_na=False, na_values=[""], **kw)


def test_every_artifact_is_written(toy):
    for f in OUT_FILES:
        assert (toy.maps / f).is_file(), f
    assert (toy.maps / "reference" / "A" / "kin.fasta").is_file()
    assert (toy.maps / "reference" / "B" / "bind.fasta").is_file()


def test_manifest_records_inputs_aligner_and_counts(toy):
    m = json.loads((toy.maps / "manifest.json").read_text(encoding="utf-8"))
    assert m["schema_version"] == 1 and m["anchor"] == "A" and m["groups"] == ["A", "B"]
    assert m["families"] == ["bind", "kin"]
    assert m["aligner"]["matrix"] == "BLOSUM62" and m["aligner"]["biopython"]
    assert m["counts"]["sequences"] == 9 and m["counts"]["protein_chains"] == 7
    assert m["counts"]["chains_unregistered"] == 1 and m["counts"]["chains_low_identity"] == 1
    assert any(k.endswith("toy_a.cif") for k in m["inputs_sha256"])
    assert all(len(v) == 64 for v in m["inputs_sha256"].values())
    assert m["reference_to_reference"]["B->A.kin"]["identity"] > 0.6


def test_chain_assignment_is_by_alignment_not_by_description(toy):
    c = _csv(toy, "structure_chains.csv").set_index(["structure", "chain"])
    assert c.loc[("toy_a", "A"), "family"] == "kin" and c.loc[("toy_a", "B"), "family"] == "kin"
    assert c.loc[("toy_a", "C"), "family"] == "bind"
    assert c.loc[("toy_b", "A"), "group"] == "B"                     # best reference is the diverged group's
    assert c.loc[("toy_a2", "A"), "scheme_variant"] == "12col"


def test_an_unrelated_chain_is_unregistered_and_a_nucleic_acid_chain_is_not_listed(toy):
    c = _csv(toy, "structure_chains.csv")
    d = c[(c["structure"] == "toy_a") & (c["chain"] == "D")].iloc[0]
    assert d["quality"] == "unregistered" and pd.isna(d["family"])
    assert d["score_z"] != d["score_z"] or d["score_z"] < 8             # nothing like a family member
    assert "E" not in set(c[c["structure"] == "toy_a"]["chain"])        # the DNA chain
    sm = _csv(toy, "structure_maps.csv", dtype={"chain": str})
    assert (sm[(sm["structure"] == "toy_a") & (sm["chain"] == "D")]["status"] == "unregistered_chain").all()


def test_low_identity_chain_is_registered_but_flagged(toy):
    c = _csv(toy, "structure_chains.csv").set_index(["structure", "chain"])
    far = c.loc[("toy_far", "A")]
    assert far["family"] == "kin" and far["quality"] == "low_identity"
    assert 0.3 < far["anchor_identity"] < 0.5 and far["score_z"] > 8


def test_composed_group_chain_records_both_hops_and_the_direct_check(toy):
    b = _csv(toy, "structure_chains.csv").set_index(["structure", "chain"]).loc[("toy_b", "A")]
    assert b["hop_a_identity"] == 1.0                                  # chain IS the group B reference
    assert 0.6 < b["anchor_identity"] < 0.9                            # composed, through B -> A
    assert b["direct_agreement"] > 0.95                                # direct alignment agrees with the composition


def test_sequence_quality_flags(toy):
    s = _csv(toy, "sequences.csv").set_index("seq_id")
    assert s.loc["a_kin_frag", "quality"] == "partial"                 # 81 of 150 positions
    assert s.loc["a_kin_x", "n_unknown"] == 10 and s.loc["a_kin_x", "quality"] == "ok"
    assert s.loc["a_kin_ref", "quality"] == "ok" and s.loc["a_kin_ref", "identity"] == 1.0
    assert s.loc["a_kin_ins", "n_inserted"] == 2 and s.loc["a_kin_del", "n_gap"] == 1


def test_extra_input_columns_pass_through(toy):
    s = _csv(toy, "sequences.csv").set_index("seq_id")
    assert (s["source"] == "toy").all()


def test_fasta_paths_are_stored_relative_to_the_maps_directory(toy):
    s = _csv(toy, "sequences.csv").set_index("seq_id")
    assert s.loc["a_kin_del", "fasta"] == "../seqs.fasta"
    assert s.loc["a_kin_ins", "record"] == "insertion_variant" and s.loc["a_kin_del", "record"] == "a_kin_del"
    assert pd.isna(s.loc["a_kin_ref", "fasta"]) and s.loc["a_kin_ref", "sequence"] == toy.seqs["a_kin_ref"]


def test_hash_is_the_sha256_of_the_normalised_sequence(toy):
    from coords.seqio import sha256_seq
    s = _csv(toy, "sequences.csv").set_index("seq_id")
    for sid, seq in toy.seqs.items():
        assert s.loc[sid, "seq_hash"] == sha256_seq(seq)
    m = _csv(toy, "sequence_maps.csv")
    assert set(m["hash12"]) == {h[:12] for h in s["seq_hash"]}


def test_structures_table_records_a_relocatable_path_and_hash(toy):
    t = _csv(toy, "structures.csv").set_index("structure")
    assert t.loc["toy_a", "path"] == "../structures/toy_a.cif" and t.loc["toy_a", "entry_id"] == "TOYA"
    assert len(t.loc["toy_a", "sha256"]) == 64


def test_the_build_is_deterministic(toy, tmp_path):
    cfg = load_config(toy.config)
    cfg.out_dir = toy.root / "maps_again"      # same depth, so relative paths are the same
    build_all(cfg, verbose=False)
    for f in OUT_FILES:
        if f == "manifest.json":
            continue
        assert (cfg.out_dir / f).read_bytes() == (toy.maps / f).read_bytes(), f


# ---- refusals -------------------------------------------------------------------------------
@pytest.fixture()
def fresh(tmp_path):
    return make_toy_project(tmp_path / "p")


def test_duplicate_structure_stem_is_refused(fresh):
    (fresh.root / "structures" / "sub").mkdir()
    shutil.copy(fresh.root / "structures" / "toy_a.cif", fresh.root / "structures" / "sub" / "toy_a.cif")
    with pytest.raises(ConfigError, match="share the key 'toy_a'"):
        build_all(load_config(fresh.config), verbose=False)


def test_a_sequence_whose_hash_no_longer_matches_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    t.loc[t["seq_id"] == "a_kin_del", "seq_hash"] = "0" * 64
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="differs from the table's seq_hash"):
        build_all(load_config(fresh.config), verbose=False)


def test_a_sequence_with_no_reference_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    t.loc[0, "family"] = "nonexistent"
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="no reference for group 'A' family 'nonexistent'"):
        build_all(load_config(fresh.config), verbose=False)


def test_missing_fasta_record_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    t.loc[t["seq_id"] == "a_kin_ins", "record"] = "no_such_record"
    t["family"] = t["family"]
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="no record"):
        build_all(load_config(fresh.config), verbose=False)


def test_missing_required_column_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str).drop(columns=["family"])
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="missing column"):
        build_all(load_config(fresh.config), verbose=False)


def test_duplicate_seq_id_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    pd.concat([t, t.iloc[[0]]]).to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="duplicate seq_id"):
        build_all(load_config(fresh.config), verbose=False)


def test_a_build_without_sequences_structures_or_sites_still_works(tmp_path):
    t = make_toy_project(tmp_path / "p")
    (t.root / "coords.toml").write_text("""[project]
anchor = "A"
[references.A]
fasta = "refs/A.fasta"
""", encoding="utf-8")
    m = build_all(load_config(t.config), verbose=False)
    assert m["counts"]["sequences"] == 0 and m["counts"]["structures"] == 0 and m["counts"]["sites"] == 0
    from coords import Frame
    fr = Frame.load(t.root / "maps")
    assert fr.families == ["bind", "kin"] and len(fr.holes) == 0


@pytest.mark.parametrize("section, old, new, fragment", [
    ("structures", 'dir = "structures"', 'dir = "missing_dir"', r"\[structures\] dir not found"),
    ("sequences", 'table = "sequences.csv"', 'table = "missing.csv"', r"\[sequences\] table not found"),
    ("sites", 'table = "sites.csv"', 'table = "missing.csv"', r"\[sites\] table not found"),
])
def test_missing_inputs_are_named(fresh, section, old, new, fragment):
    fresh.config.write_text(fresh.config.read_text(encoding="utf-8").replace(old, new), encoding="utf-8")
    with pytest.raises(ConfigError, match=fragment):
        build_all(load_config(fresh.config), verbose=False)


def test_a_sequence_fasta_that_does_not_exist_is_named(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    t.loc[t["seq_id"] == "a_kin_del", "fasta"] = "gone.fasta"
    t.loc[t["seq_id"] == "a_kin_del", "seq_hash"] = ""
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="fasta not found"):
        build_all(load_config(fresh.config), verbose=False)


def test_a_sequence_with_neither_text_nor_fasta_is_refused(fresh):
    t = pd.read_csv(fresh.root / "sequences.csv", keep_default_na=False, dtype=str)
    t.loc[t["seq_id"] == "a_kin_ref", "sequence"] = ""
    t.to_csv(fresh.root / "sequences.csv", index=False)
    with pytest.raises(ConfigError, match="neither `sequence` nor `fasta`"):
        build_all(load_config(fresh.config), verbose=False)


def test_thresholds_and_aligner_come_from_the_config(fresh):
    fresh.config.write_text(fresh.config.read_text(encoding="utf-8") + """
[aligner]
matrix = "BLOSUM45"
open_gap = -12.0

[thresholds]
low_chain_identity = 0.95
""", encoding="utf-8")
    m = build_all(load_config(fresh.config), verbose=False)
    assert m["aligner"]["matrix"] == "BLOSUM45" and m["aligner"]["open_gap"] == -12.0
    assert m["thresholds"]["low_chain_identity"] == 0.95
    assert m["counts"]["chains_low_identity"] >= 2            # toy_b (composed 0.75) is now flagged too
