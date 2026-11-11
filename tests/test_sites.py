"""Sites: a table of positions whose numbering schemes differ, resolved once into the anchor frame."""
from __future__ import annotations

import pandas as pd
import pytest

from coords import Status
from coords.build import build_all
from coords.config import ConfigError, load_config
from coords.testing import make_toy_project


def _row(df, **kw):
    m = pd.Series(True, index=df.index)
    for k, v in kw.items():
        m &= df[k] == v
    assert m.sum() == 1, (kw, int(m.sum()))
    return df[m].iloc[0]


def test_scheme_offsets_are_added_to_the_position(fr):
    t = fr.sites("kin")
    assert _row(t, position=25, scheme="full")["ref_pos"] == 25
    assert _row(t, position=20, scheme="mature")["ref_pos"] == 25            # same residue, other numbering
    assert _row(t, position=10, scheme="mature", phenotype="z")["ref_pos"] == 15


def test_two_numberings_of_one_residue_resolve_to_the_same_reference_position(fr):
    t = fr.sites("kin")
    a, b = _row(t, position=25, scheme="full"), _row(t, position=20, scheme="mature")
    assert a["ref_pos"] == b["ref_pos"] and a["ref_aa"] == b["ref_aa"]


def test_rows_without_a_usable_scheme_are_excluded_not_guessed(fr):
    t = fr.sites("kin")
    unresolved = t[t["status"] == Status.UNRESOLVED_SCHEME]
    assert set(unresolved["phenotype"]) == {"blank", "t", "u"}               # blank scheme, "unresolved", unknown name
    assert unresolved["ref_pos"].isna().all() and (unresolved["ref_aa"] == "").all()


def test_a_position_off_the_reference_or_a_family_without_one_is_no_ref_position(fr):
    t = fr.sites()
    assert _row(t, phenotype="far")["status"] == Status.NO_REF_POSITION       # bind 1000 > 90
    assert _row(t, phenotype="nofam")["status"] == Status.NO_REF_POSITION     # no such family


def test_matches_ref_compares_the_listed_residue_with_the_anchor(fr, toy):
    t = fr.sites()
    assert bool(_row(t, position=25, scheme="full")["matches_ref"]) is True
    assert bool(_row(t, position=20, scheme="mature")["matches_ref"]) is True
    z = _row(t, phenotype="z")
    assert z["matches_ref"] == (toy.anchor_refs["kin"][14] == "W")
    assert pd.isna(_row(t, phenotype="blank")["matches_ref"])                 # no residue given


def test_extra_columns_are_kept(fr):
    assert set(fr.sites()["phenotype"]) >= {"x", "y", "z", "b"}


def test_every_unresolved_row_is_recorded_as_a_hole(fr):
    fr.sites()
    st = fr.holes["status"].value_counts().to_dict()
    assert st == {Status.UNRESOLVED_SCHEME: 3, Status.NO_REF_POSITION: 2}
    fr.clear_holes()
    fr.sites("bind")
    assert fr.holes["status"].tolist() == [Status.NO_REF_POSITION]            # only bind's own far row


def test_the_resolved_position_can_be_read_through_a_sequence(fr):
    s = fr.sites("kin")
    p = int(_row(s, position=20, scheme="mature")["ref_pos"])
    assert fr.rows("kin", "a_kin_ins", [p])[0] == p + 1                       # shifted by the upstream insertion


def test_site_columns_are_configurable(tmp_path):
    t = make_toy_project(tmp_path / "p")
    pd.DataFrame({"gene": ["kin"], "num": [3], "numbering": ["full"], "aa": ["K"]}).to_csv(t.root / "sites.csv", index=False)
    cfg = t.root / "coords.toml"
    cfg.write_text(cfg.read_text(encoding="utf-8").replace(
        '[sites]\ntable = "sites.csv"', '[sites]\ntable = "sites.csv"\nfamily_column = "gene"\n'
        'position_column = "num"\nscheme_column = "numbering"\nresidue_column = "aa"'), encoding="utf-8")
    build_all(load_config(cfg), verbose=False)
    out = pd.read_csv(t.maps / "sites_resolved.csv", keep_default_na=False, na_values=[""])
    assert out.loc[0, "family"] == "kin" and out.loc[0, "ref_pos"] == 3 and out.loc[0, "residue"] == "K"


def test_a_non_integer_position_is_refused(tmp_path):
    t = make_toy_project(tmp_path / "p")
    s = pd.read_csv(t.root / "sites.csv", dtype=str, keep_default_na=False)
    s.loc[0, "position"] = "12A"
    s.to_csv(t.root / "sites.csv", index=False)
    with pytest.raises(ConfigError, match="not an integer"):
        build_all(load_config(t.config), verbose=False)


def test_a_missing_site_column_is_refused(tmp_path):
    t = make_toy_project(tmp_path / "p")
    pd.read_csv(t.root / "sites.csv", dtype=str, keep_default_na=False).drop(columns=["scheme"]).to_csv(
        t.root / "sites.csv", index=False)
    with pytest.raises(ConfigError, match="missing column"):
        build_all(load_config(t.config), verbose=False)
