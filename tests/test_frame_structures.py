"""g: reference position <-> structure residue, and the geometric query. Each case is one where the author
number is not the reference position. The expected author numbers come from coords.testing.toy_a_auth."""
from __future__ import annotations

import pandas as pd
import pytest

from coords import Frame, HoleError, ResKey, Status
from coords.testing import MISSING_LOOP, N_TERM_UNOBSERVED, toy_a_auth


# ---- author numbering is not the reference position ----------------------------------------------
@pytest.mark.parametrize("p", [10, 39, 40, 42, 100, 150])
def test_struct_returns_the_author_residue_one_per_copy(fr, p):
    auth, ins = toy_a_auth(p)
    keys = fr.struct("kin", "toy_a", p)
    assert keys == [ResKey(1, "A", auth, ins), ResKey(1, "B", auth, ins)]
    assert all(fr.ref_pos("toy_a", k) == p for k in keys)           # and back
    assert auth != p or p > 100                                      # the bare integer would be wrong


def test_the_offset_changes_inside_the_chain_so_no_scalar_offset_exists(fr):
    offsets = {toy_a_auth(p)[0] - p for p in range(10, 151)}
    assert offsets == {9, 8}
    assert {fr.struct("kin", "toy_a", p)[0].auth_seq_id - p for p in range(10, 151) if p not in MISSING_LOOP} == {9, 8}


def test_insertion_code_residues_are_distinct_from_their_parent(fr):
    assert fr.struct("kin", "toy_a", 40, chain="A") == [ResKey(1, "A", 49, "")]
    assert fr.struct("kin", "toy_a", 41, chain="A") == [ResKey(1, "A", 49, "A")]     # same number, other residue
    assert fr.struct("kin", "toy_a", 42, chain="A") == [ResKey(1, "A", 50, "")]
    assert fr.ref_pos("toy_a", ResKey(1, "A", 49, "")) == 40
    assert fr.ref_pos("toy_a", ResKey(1, "A", 49, "A")) == 41
    assert fr.ref_pos("toy_a", ResKey(1, "A", 49, "B")) is None
    assert fr.holes.iloc[-1]["status"] == Status.UNREGISTERED_CHAIN


def test_the_same_region_under_two_numberings_maps_to_one_position(fr):
    a = fr.struct("kin", "toy_a", 20, chain="A")[0]
    b = fr.struct("kin", "toy_a2", 20)[0]
    assert (a.auth_seq_id, b.auth_seq_id) == (29, 520)
    assert fr.ref_pos("toy_a", a) == fr.ref_pos("toy_a2", b) == 20


def test_a_family_numbered_from_one_is_the_identity_regression(fr):
    sm = fr.smaps
    c = sm[(sm["structure"] == "toy_a") & (sm["chain"] == "C") & sm["observed"]]
    assert len(c) == 90 and (c["auth_seq_id"] == c["ref_pos"]).all()
    assert (c["status"] == Status.OK).all()


# ---- holes ----------------------------------------------------------------------------------------
@pytest.mark.parametrize("p", list(MISSING_LOOP))
def test_a_missing_interior_loop_is_a_hole_on_that_chain_only(fr, p):
    assert fr.struct("kin", "toy_a", p, chain="A") == []
    assert fr.holes.iloc[-1]["status"] == Status.NOT_OBSERVED
    fr.clear_holes()
    keys = fr.struct("kin", "toy_a", p)                              # chain B has the residue
    assert [k.chain for k in keys] == ["B"]
    assert fr.holes["status"].tolist() == [Status.NOT_OBSERVED]      # chain A's copy is recorded, not returned


@pytest.mark.parametrize("p", [1, 5, 9])
def test_a_chain_not_starting_at_one_has_no_residues_before_its_first_modelled_one(fr, p):
    assert p in N_TERM_UNOBSERVED
    assert fr.struct("kin", "toy_a", p, chain="A") == []
    assert fr.struct("kin", "toy_a", 10, chain="A") == [ResKey(1, "A", 19, "")]


def test_a_position_the_chain_does_not_cover_is_not_on_structure(fr):
    assert fr.struct("kin", "toy_b", 42) == []                      # deleted in the B-group structure
    assert fr.holes.iloc[-1]["status"] == Status.NOT_ON_STRUCTURE
    assert fr.struct("bind", "toy_b", 3) == []                       # a family the structure does not contain


def test_unobserved_with_an_inferred_number_in_the_12col_variant_is_still_not_observed(fr):
    assert fr.struct("kin", "toy_a2", 1) == []                       # SEQRES 1 has an inferred number but no atoms
    assert fr.holes.iloc[-1]["status"] == Status.NOT_OBSERVED


def test_on_hole_raise_applies_to_structures(toy):
    fr = Frame.load(toy.maps, on_hole="raise")
    with pytest.raises(HoleError, match="not_observed"):
        fr.struct("kin", "toy_a", 62, chain="A")
    assert fr.struct("kin", "toy_a", 62, chain="B")


# ---- group B structure: composed through the reference-to-reference map -----------------------------
def test_a_group_b_structure_maps_through_the_anchor(fr):
    assert fr.struct("kin", "toy_b", 60) == [ResKey(1, "A", 157, "")]      # anchor 60 -> B 57 -> author 100 + 57
    assert fr.struct("kin", "toy_b", 40) == [ResKey(1, "A", 140, "")]
    assert fr.struct("kin", "toy_b", 102) == [ResKey(1, "A", 100 + 101, "")]  # past the GS insertion


def test_a_residue_inserted_relative_to_the_anchor_has_no_reference_position(fr):
    assert fr.ref_pos("toy_b", ResKey(1, "A", 199, "")) is None             # the inserted G of B
    assert fr.holes.iloc[-1]["status"] == Status.NO_REF_POSITION


def test_structure_maps_keep_the_group_reference_position_too(fr):
    sm = fr.smaps
    r = sm[(sm["structure"] == "toy_b") & (sm["auth_seq_id"] == 157)].iloc[0]
    assert (r["ref_pos"], r["group_ref_pos"]) == (60, 57)


# ---- chain bookkeeping ------------------------------------------------------------------------------
def test_chain_family_and_unregistered_chains(fr):
    assert fr.chain_family("toy_a", "A") == "kin" and fr.chain_family("toy_a", "C") == "bind"
    assert fr.chain_family("toy_a", "D") == ""                       # unrelated chain
    assert fr.chain_family("toy_a", "Z") == ""


def test_nucleic_acid_ligand_and_water_never_reach_the_atom_table(fr):
    a = fr.atoms("toy_a")
    assert set(a["chain"]) == {"A", "B", "C", "D"}
    assert not a["comp_id"].isin(["NAG", "HOH"]).any()


def test_unknown_structure_is_an_error(fr):
    with pytest.raises(FileNotFoundError, match="unknown structure"):
        fr.structure_path("nope")


def test_a_structure_file_that_moved_is_reported(toy, tmp_path):
    import shutil
    d = tmp_path / "m"
    shutil.copytree(toy.maps, d)                                      # relative paths now point nowhere
    with pytest.raises(FileNotFoundError, match="has moved"):
        Frame.load(d).atoms("toy_a")


# ---- neighbours ---------------------------------------------------------------------------------------
def _members(df):
    return set(zip(df["neighbor_reskey"].map(lambda k: k.chain), df["neighbor_ref_pos"]))


def test_neighbors_are_reference_positions_with_distances(fr):
    nb = fr.neighbors("kin", "toy_a", 10, 7.0)
    assert _members(nb) == {("A", 10), ("A", 11), ("B", 9), ("B", 10), ("B", 11)}
    row = nb[(nb["neighbor_ref_pos"] == 11) & (nb["neighbor_reskey"].map(lambda k: k.chain) == "A")].iloc[0]
    assert row["distance"] == pytest.approx(3.8) and not row["is_query"]
    assert nb[nb["is_query"]]["neighbor_ref_pos"].tolist() == [10]
    assert (nb["status"] == Status.OK).all()


def test_an_unobserved_neighbour_never_appears(fr):
    nb = fr.neighbors("kin", "toy_a", 10, 7.0, include_holes=True)
    assert ("A", 9) not in _members(nb)                               # chain A is not modelled before position 10
    assert ("B", 9) in _members(nb)                                   # chain B has it


def test_other_family_and_unregistered_neighbours_are_recorded_and_dropped(fr):
    full = fr.neighbors("kin", "toy_a", 10, 7.0, include_holes=True)
    by = full.groupby("status").size().to_dict()
    assert by[Status.OTHER_FAMILY] == 1 and by[Status.UNREGISTERED_CHAIN] == 3
    assert set(full[full["status"] == Status.OTHER_FAMILY]["neighbor_family"]) == {"bind"}
    assert set(fr.holes["status"]) >= {Status.OTHER_FAMILY, Status.UNREGISTERED_CHAIN}
    default = fr.neighbors("kin", "toy_a", 10, 7.0)
    assert set(default["neighbor_family"]) == {"kin"}


def test_cross_family_keeps_a_heteromer_partner_in_its_own_frame(fr):
    nb = fr.neighbors("kin", "toy_a", 10, 7.0, cross_family=True)
    c = nb[nb["neighbor_family"] == "bind"]
    assert len(c) == 1 and c.iloc[0]["neighbor_ref_pos"] == 10 and c.iloc[0]["status"] == Status.OK


def test_all_atom_centre_finds_contacts_a_ca_centre_misses(fr):
    ca = fr.neighbors("kin", "toy_a", 30, 4.5, center="CA")
    al = fr.neighbors("kin", "toy_a", 30, 4.5, center="all")
    assert _members(ca) == {("A", 29), ("A", 30), ("A", 31)}          # CA to CA across chains is 5.0 A
    assert ("B", 30) in _members(al) and _members(ca) < _members(al)
    assert al[al["neighbor_reskey"].map(lambda k: k.chain) == "B"]["distance"].min() < 4.5


def test_neighbors_of_the_other_copy(fr):
    only_b = fr.neighbors("kin", "toy_a", 30, 4.5, copy="B")
    assert set(only_b["query_reskey"]) == {ResKey(1, "B", 39, "")}
    both = fr.neighbors("kin", "toy_a", 30, 4.5, copy="all")
    assert set(both["query_reskey"]) == {ResKey(1, "A", 39, ""), ResKey(1, "B", 39, "")}


def test_alternate_locations_do_not_move_or_duplicate_a_residue(fr):
    nb = fr.neighbors("kin", "toy_a", 20, 4.0)                        # residue 20 has altlocs A and B (+50 A)
    assert _members(nb) >= {("A", 19), ("A", 20), ("A", 21)}
    assert nb[nb["neighbor_ref_pos"] == 21].iloc[0]["distance"] == pytest.approx(3.8)
    assert not nb.duplicated(["query_reskey", "neighbor_reskey"]).any()


def test_neighbors_cross_with_sequences_and_read_the_aligned_row(fr, toy):
    seqs = ["a_kin_ref", "a_kin_ins", "b_kin_ref"]
    nb = fr.neighbors("kin", "toy_a", 30, 4.5, center="all", seq_ids=seqs)
    assert set(nb["seq_id"]) == set(seqs)
    ref = nb[nb["seq_id"] == "a_kin_ref"]
    ins = nb[nb["seq_id"] == "a_kin_ins"]
    assert (ref["row"] == ref["neighbor_ref_pos"] - 1).all()
    assert (ins["row"] == ins["neighbor_ref_pos"] + 1).all()          # WW inserted upstream: +2 over p - 1
    assert (nb["row"] + 1 == nb["seq_pos"]).all()
    assert (ins["seq_aa"] == ins["ref_aa"]).all()                     # and it is the same residue
    assert fr.holes[fr.holes["status"] == Status.GAP_IN_SEQUENCE].empty


def test_a_neighbour_in_a_sequence_hole_is_recorded_and_not_read(fr):
    nb = fr.neighbors("kin", "toy_a", 30, 4.5, center="all", seq_ids=["a_kin_del"])
    assert 30 not in set(nb["neighbor_ref_pos"])                      # position 30 is deleted in a_kin_del
    held = fr.neighbors("kin", "toy_a", 30, 4.5, center="all", seq_ids=["a_kin_del"], include_holes=True)
    gap = held[held["status"] == Status.GAP_IN_SEQUENCE]
    assert set(gap["neighbor_ref_pos"]) == {30} and gap["row"].isna().all()


def test_cross_family_neighbours_read_their_own_familys_sequence(fr, toy):
    nb = fr.neighbors("kin", "toy_a", 10, 7.0, cross_family=True, seq_ids=["a_kin_ref", "a_bind"])
    bind = nb[nb["neighbor_family"] == "bind"]
    assert set(bind["seq_id"]) == {"a_bind"} and bind["row"].tolist() == [9]
    assert bind["seq_aa"].iloc[0] == toy.anchor_refs["bind"][9]
    kin = nb[nb["neighbor_family"] == "kin"]
    assert set(kin["seq_id"]) == {"a_kin_ref"}


def test_neighbors_of_an_unobserved_site_raise_under_raise(toy):
    fr = Frame.load(toy.maps, on_hole="raise")
    with pytest.raises(HoleError):
        fr.neighbors("kin", "toy_a", 62, 5.0, copy="A")


def test_bad_arguments(fr):
    with pytest.raises(ValueError, match="center"):
        fr.neighbors("kin", "toy_a", 10, 5.0, center="CB")
    with pytest.raises(KeyError, match="unknown seq_id"):
        fr.neighbors("kin", "toy_a", 10, 5.0, seq_ids=["nope"])
