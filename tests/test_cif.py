"""The mmCIF reader: both scheme variants, protein chains only, author numbering kept as deposited."""
from __future__ import annotations

import pandas as pd
import pytest

from coords.cif import read_structure, seqres_by_chain
from coords.testing import ChainSpec, random_protein, write_mmcif

SEQ = random_protein(30, 21)


def _chain(**kw):
    auth = [None, None, 3, 4, 5] + [6 + i for i in range(5, 30)]
    auth[10] = 10                                    # a number reused with an insertion code below
    ins = [""] * 30
    ins[11] = "A"
    auth[11] = 10                                    # 10 and 10A
    base = dict(chain="A", entity="1", seq=SEQ, auth=auth, ins=ins, description="kinase")
    base.update(kw)
    return ChainSpec(**base)


@pytest.fixture(params=["7col", "12col"])
def st(request, tmp_path):
    path = write_mmcif(tmp_path / "x.cif", "T001", [
        _chain(altloc_residues=[7], ligand=True),
        ChainSpec("D", "2", "ACGTACGT", list(range(1, 9)), (), "dna", poly_type="polydeoxyribonucleotide"),
    ], variant=request.param)
    s = read_structure(path)
    s.variant = request.param
    return s


def test_scheme_variant_is_detected(st):
    assert st.scheme_variant == st.variant


def test_entry_id_and_entities(st):
    assert st.entry_id == "T001"
    assert set(st.entities["entity_id"]) >= {"1", "2"}
    assert st.entities.set_index("entity_id").loc["1", "description"] == "kinase"


def test_unobserved_residues_are_flagged_in_both_variants(st):
    s = st.scheme[st.scheme["chain"] == "A"].sort_values("seq_id")
    assert s["observed"].tolist()[:5] == [False, False, True, True, True]
    assert int(s["observed"].sum()) == 28


def test_author_numbers_and_insertion_codes_are_kept_as_deposited(st):
    s = st.scheme[st.scheme["chain"] == "A"].sort_values("seq_id").reset_index(drop=True)
    assert (int(s.loc[10, "auth_seq_id"]), s.loc[10, "ins_code"]) == (10, "")
    assert (int(s.loc[11, "auth_seq_id"]), s.loc[11, "ins_code"]) == (10, "A")     # same number, different residue
    assert int(s.loc[2, "auth_seq_id"]) == 3                                      # not SEQRES index (3 is residue 3 by luck)
    assert int(s.loc[20, "auth_seq_id"]) == 6 + 20                                # offset differs from the SEQRES index


def test_unobserved_number_is_missing_in_7col_and_inferred_in_12col(st):
    s = st.scheme[st.scheme["chain"] == "A"].sort_values("seq_id").reset_index(drop=True)
    if st.variant == "7col":
        assert s.loc[0, "auth_seq_id"] is pd.NA or pd.isna(s.loc[0, "auth_seq_id"])
    else:
        assert int(s.loc[0, "auth_seq_id"]) == 1


def test_only_protein_atoms_survive(st):
    a = st.atoms
    assert set(a["chain"]) == {"A"}                                    # DNA chain D dropped
    assert not a["comp_id"].isin(["NAG", "HOH", "DA", "DC", "DG", "DT"]).any()
    assert not st.scheme[st.scheme["chain"] == "D"].empty              # but its SEQRES is read and typed
    assert st.scheme[st.scheme["chain"] == "D"]["poly_type"].iloc[0] == "polydeoxyribonucleotide"


def test_one_location_per_atom_blank_or_A_wins(st):
    ca = st.atoms[(st.atoms["atom"] == "CA") & (st.atoms["auth_seq_id"] == 6 + 7)]
    assert len(ca) == 1
    assert ca.iloc[0]["x"] == pytest.approx(3.8 * 7)                    # altloc A, not the +50 A copy B


def test_unobserved_residues_have_no_atoms(st):
    assert not st.atoms["auth_seq_id"].isin([1, 2]).any() or st.variant == "12col"
    n_ca = int((st.atoms["atom"] == "CA").sum())
    assert n_ca == 28


def test_seqres_by_chain_returns_the_deposited_sequence_including_unobserved(st):
    seq, g = seqres_by_chain(st)["A"]
    assert seq == SEQ and len(g) == 30


def test_semicolon_text_fields_and_quoted_tokens(tmp_path):
    text = """data_X
_entry.id X1
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'it's a "quoted" name'
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
_entity_poly.pdbx_seq_one_letter_code_can
1 'polypeptide(L)' A
;ACDEF
GHIK
;
loop_
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
1 A ALA 1 A 1 .
1 A CYS 2 A 2 .
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.0 0.0 0.0 1 ALA A CA 1
ATOM 2 C CA . CYS A 1 2 ? 3.8 0.0 0.0 2 CYS A CA 1
"""
    p = tmp_path / "q.cif"
    p.write_text(text, encoding="utf-8")
    s = read_structure(p)
    assert s.entities.iloc[0]["description"] == 'it\'s a "quoted" name'
    assert s.entities.iloc[0]["seq_can"] == "ACDEFGHIK"
    assert len(s.atoms) == 2 and s.structure == "q"


def test_file_without_a_scheme_is_an_error(tmp_path):
    p = tmp_path / "bad.cif"
    p.write_text("data_X\n_entry.id X\n", encoding="utf-8")
    with pytest.raises(ValueError, match="_atom_site"):
        read_structure(p)
