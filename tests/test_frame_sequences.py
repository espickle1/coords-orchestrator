"""f_i: reference position -> row of a per-residue table. Each case is one where the right answer differs
from the bare index, so code that passes the integer straight through fails it. Expected values come from
how the toy sequences were constructed (coords.testing), not from the library under test."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from coords import Frame, HoleError, Status
from coords.testing import make_toy_project

KIN = 150


def _read(toy, sid, row):
    return toy.seqs[sid][row]


def test_identity_sequence_reads_row_p_minus_1(fr, toy):
    assert fr.rows("kin", "a_kin_ref", [1, 75, 150]).tolist() == [0, 74, 149]


def test_deletion_upstream_shifts_every_later_row(fr, toy):
    ps = [10, 29, 31, 80, 150]                                     # position 30 is the deleted residue
    rows = fr.rows("kin", "a_kin_del", ps).tolist()
    assert rows == [9, 28, 29, 78, 148]
    for p, r in zip(ps, rows):                                      # and the residue read is the right one
        assert _read(toy, "a_kin_del", r) == toy.anchor_refs["kin"][p - 1]


def test_the_deleted_position_is_a_hole_not_the_next_row(fr):
    out = fr.rows("kin", "a_kin_del", [29, 30, 31])
    assert pd.isna(out[1]) and out[0] == 28 and out[2] == 29
    assert fr.holes["status"].tolist() == [Status.GAP_IN_SEQUENCE]
    assert fr.holes.iloc[0]["ref_pos"] == 30 and fr.holes.iloc[0]["seq_id"] == "a_kin_del"


def test_insertion_upstream_shifts_rows_up(fr, toy):
    ps = [5, 19, 20, 100, 150]                                     # WW inserted before position 20
    rows = fr.rows("kin", "a_kin_ins", ps).tolist()
    assert rows == [4, 18, 21, 101, 151]
    for p, r in zip(ps, rows):
        assert _read(toy, "a_kin_ins", r) == toy.anchor_refs["kin"][p - 1]


def test_the_inserted_residues_have_no_reference_position(fr):
    inv = fr.inverse_map("kin", "a_kin_ins")
    assert len(inv) == KIN + 2
    assert inv[19] == 0 and inv[20] == 0                            # the two W
    assert inv[18] == 19 and inv[21] == 20


def test_an_unsequenced_run_is_a_hole_but_its_flanks_keep_their_rows(fr):
    out = fr.rows("kin", "a_kin_x", [58, 59, 60, 65, 69, 70, 71])
    assert out[0] == 57 and out[1] == 58 and out[5] == 69 and out[6] == 70      # flanks unmoved
    assert all(pd.isna(out[i]) for i in (2, 3, 4))
    st = fr.holes["status"].tolist()
    assert st == [Status.UNKNOWN_RESIDUE] * 3


def test_a_fragment_has_holes_outside_its_span(fr, toy):
    assert fr.rows("kin", "a_kin_frag", [40, 41, 120]).tolist() == [0, 1, 80]    # span is positions 40..120
    out = fr.rows("kin", "a_kin_frag", [39, 121, 150])
    assert out.isna().all()
    assert set(fr.holes["status"]) == {Status.GAP_IN_SEQUENCE}


def test_seq_pos_is_one_based_and_matches_rows(fr):
    assert fr.seq_pos("kin", "a_kin_ins", 20) == 22
    assert fr.rows("kin", "a_kin_ins", [20])[0] == 21


def test_row_table_carries_the_residues_and_the_status(fr, toy):
    t = fr.row_table("kin", "a_kin_del", [29, 30, 31])
    assert t["status"].tolist() == [Status.OK, Status.GAP_IN_SEQUENCE, Status.OK]
    assert t["seq_aa"].tolist()[0] == t["ref_aa"].tolist()[0] and t["seq_aa"].tolist()[1] == ""
    assert t["row"].isna().tolist() == [False, True, False]
    assert t["seq_aa"].tolist()[2] == toy.anchor_refs["kin"][30]


def test_row_table_reads_residues_from_every_input_route(fr, toy):
    """inline sequence, a FASTA record named by seq_id, and a FASTA record named by a `record` column"""
    for sid in ("a_kin_ref", "a_kin_del", "a_kin_ins"):
        t = fr.row_table("kin", sid, [10])
        assert t.loc[0, "seq_aa"] == toy.anchor_refs["kin"][9], sid


# ---- a second group, composed through the reference-to-reference map ----------------------------
def _expected_b(p):
    """anchor position -> group B position, from how kin B was built (3 deleted after 40, GS after 101)"""
    return p if p <= 40 else None if p <= 43 else p - 3 if p <= 101 else p - 1


def test_group_b_rows_are_composed_through_the_anchor(fr, toy):
    for p in (1, 40, 44, 60, 101, 102, 150):
        assert fr.seq_pos("kin", "b_kin_ref", p) == _expected_b(p), p
    assert fr.holes.empty


def test_group_b_deletion_relative_to_the_anchor_is_a_hole(fr):
    out = fr.rows("kin", "b_kin_ref", [40, 41, 42, 43, 44])
    assert pd.isna(out[1]) and pd.isna(out[2]) and pd.isna(out[3])
    assert out[0] == 39 and out[4] == 40
    assert set(fr.holes["status"]) == {Status.GAP_IN_SEQUENCE}


def test_group_b_residues_read_are_the_aligned_ones(fr, toy):
    kin_a, kin_b = toy.anchor_refs["kin"], toy.group_b_refs["kin"]
    same = 0
    for p in range(1, KIN + 1):
        q = _expected_b(p)
        if q is not None and kin_a[p - 1] == kin_b[q - 1]:
            same += 1
            assert _read(toy, "b_kin_ref", q - 1) == kin_a[p - 1]
    assert same > 80


def test_group_b_insertion_has_no_anchor_position_and_to_anchor_says_so(fr):
    assert fr.to_anchor("B", "kin", 57) == 60                       # B position 57 is anchor position 60
    assert fr.to_anchor("B", "kin", 99) is None                     # the inserted G
    assert fr.holes.iloc[-1]["status"] == Status.NO_REF_POSITION
    arr = fr.own_to_anchor_array("B", "kin")
    assert arr[56] == 60 and arr[98] == 0 and len(arr) == len(fr.refs[("B", "kin")])
    assert fr.to_anchor("A", "kin", 17) == 17                       # the anchor group maps to itself


def test_inverse_map_round_trips_for_the_composed_group(fr, toy):
    inv = fr.inverse_map("kin", "b_kin_ref")
    assert len(inv) == len(toy.seqs["b_kin_ref"])
    for p in range(1, KIN + 1):
        q = _expected_b(p)
        if q is not None:
            assert inv[q - 1] == p
    assert inv[98] == 0 and inv[99] == 0                            # GS inserted after anchor position 101


def test_substituted_group_b_sequence_keeps_its_rows(fr):
    assert fr.seq_pos("kin", "b_kin_mut", 60) == _expected_b(60)


# ---- errors are errors, not holes ---------------------------------------------------------------
def test_position_outside_the_reference_is_an_error(fr):
    with pytest.raises(ValueError, match="outside the kin reference"):
        fr.rows("kin", "a_kin_ref", [KIN + 1])
    with pytest.raises(ValueError):
        fr.rows("kin", "a_kin_ref", [0])


def test_wrong_family_for_the_sequence_is_an_error(fr):
    with pytest.raises(ValueError, match="is a kin sequence, not bind"):
        fr.rows("bind", "a_kin_ref", [3])


def test_unknown_sequence_is_an_error(fr):
    with pytest.raises(KeyError, match="unknown seq_id"):
        fr.rows("kin", "nope", [3])


def test_on_hole_raise_raises_instead_of_recording(toy):
    fr = Frame.load(toy.maps, on_hole="raise")
    assert fr.rows("kin", "a_kin_del", [29])[0] == 28
    with pytest.raises(HoleError, match="gap_in_sequence"):
        fr.rows("kin", "a_kin_del", [30])
    with pytest.raises(HoleError, match="unknown_residue"):
        fr.rows("kin", "a_kin_x", [65])
    with pytest.raises(ValueError):
        Frame.load(toy.maps, on_hole="fill")


def test_holes_can_be_cleared_and_are_a_plain_table(fr):
    fr.rows("kin", "a_kin_del", [30])
    assert len(fr.holes) == 1
    fr.clear_holes()
    assert fr.holes.empty and list(fr.holes.columns) == ["status"]


# ---- hash checks and sequence text --------------------------------------------------------------
def test_check_hash_accepts_the_mapped_string_and_refuses_another(fr, toy):
    fr.check_hash("a_kin_ins", sequence=toy.seqs["a_kin_ins"])
    from coords.seqio import sha256_seq
    fr.check_hash("a_kin_ins", seq_hash=sha256_seq(toy.seqs["a_kin_ins"]))
    with pytest.raises(ValueError, match="different string"):
        fr.check_hash("a_kin_ins", sequence=toy.seqs["a_kin_ref"])
    with pytest.raises(ValueError, match="pass sequence= or seq_hash="):
        fr.check_hash("a_kin_ins")


def test_sequence_text_that_no_longer_matches_is_not_returned(tmp_path):
    t = make_toy_project(tmp_path / "p")
    from coords.build import build_all
    from coords.config import load_config
    build_all(load_config(t.config), verbose=False)
    fr = Frame.load(t.maps)
    assert fr.sequence("a_kin_del") == t.seqs["a_kin_del"]
    (t.root / "seqs.fasta").write_text(">a_kin_del\nAAAA\n>insertion_variant\nAAAA\n", encoding="utf-8")
    fr2 = Frame.load(t.maps)
    assert fr2.sequence("a_kin_del") is None                        # edited after the map was built
    assert fr2.sequence("a_kin_ref") == t.seqs["a_kin_ref"]         # inline text is stored with the map


def test_sequences_table_carries_passthrough_columns(fr):
    assert fr.sequences.loc["a_kin_ref", "source"] == "toy"
    assert fr.sequences.loc["a_kin_ref", "length"] == KIN
    assert fr.groups == ["A", "B"] and fr.families == ["bind", "kin"]


def test_frame_refuses_maps_of_another_schema(toy, tmp_path):
    import json, shutil
    d = tmp_path / "old"
    shutil.copytree(toy.maps, d)
    m = json.loads((d / "manifest.json").read_text(encoding="utf-8"))
    m["schema_version"] = 99
    (d / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(ValueError, match="schema_version"):
        Frame.load(d)
    with pytest.raises(FileNotFoundError, match="Run `python -m coords build"):
        Frame.load(tmp_path / "empty")
    with pytest.raises(ValueError, match="anchor"):
        Frame.load(toy.maps, anchor="B")
