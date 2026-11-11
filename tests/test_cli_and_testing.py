"""The command line, the script helpers, and the public synthetic-data helpers."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from coords import Frame, Status, cli
from coords.testing import ChainSpec, homolog, make_toy_project, random_protein, synthetic_frame, write_mmcif


# ---- command line ----------------------------------------------------------------------------------
def test_build_and_info_through_main(tmp_path, capsys):
    t = make_toy_project(tmp_path / "p")
    cli.main(["build", str(t.config), "-q"])
    assert (t.maps / "manifest.json").is_file()
    with pytest.raises(SystemExit) as e:
        cli.main(["info", str(t.maps)])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert "anchor A" in out and "toy_a:D" in out and "unregistered" in out and "a_kin_frag" in out


def test_a_bad_config_exits_with_code_2_and_a_message(tmp_path, capsys):
    bad = tmp_path / "coords.toml"
    bad.write_text('[project]\nanchor = "Z"\n[references.A]\nfasta = "x.fa"\n', encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        cli.main(["build", str(bad)])
    assert e.value.code == 2
    assert "not one of the groups" in capsys.readouterr().err


def test_info_on_a_directory_that_was_never_built(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["info", str(tmp_path)])
    assert e.value.code == 2 and "coords build" in capsys.readouterr().err


def test_python_dash_m_runs(tmp_path):
    t = make_toy_project(tmp_path / "p")
    src = str(Path(__file__).resolve().parents[1] / "src")
    r = subprocess.run([sys.executable, "-m", "coords", "build", str(t.config), "-q"], capture_output=True,
                       text=True, env={**__import__("os").environ, "PYTHONPATH": src})
    assert r.returncode == 0, r.stderr
    assert (t.maps / "sequences.csv").is_file()


# ---- helpers a script uses -------------------------------------------------------------------------
def test_frame_args_and_hole_log(toy, tmp_path):
    ap = argparse.ArgumentParser()
    cli.add_frame_args(ap)
    args = ap.parse_args(["--coords-dir", str(toy.maps), "--on-hole", "skip"])
    fr = cli.frame_from_args(args)
    fr.rows("kin", "a_kin_del", [30, 30])
    fr.struct("kin", "toy_a", 62, chain="A")
    out = cli.write_holes(fr, tmp_path / "x" / "holes.csv")
    assert set(out["status"]) == {Status.GAP_IN_SEQUENCE, Status.NOT_OBSERVED}
    assert out[out["status"] == Status.GAP_IN_SEQUENCE]["n"].iloc[0] == 2          # deduplicated with a count
    assert (tmp_path / "x" / "holes.csv").is_file()
    raising = cli.frame_from_args(ap.parse_args(["--coords-dir", str(toy.maps), "--on-hole", "raise"]))
    assert raising.on_hole == "raise"
    with pytest.raises(SystemExit):
        ap.parse_args(["--on-hole", "fill"])


def test_an_empty_hole_log_is_a_valid_file(toy, tmp_path):
    out = cli.write_holes(Frame.load(toy.maps), tmp_path / "h.csv")
    assert out.empty and list(pd.read_csv(tmp_path / "h.csv").columns) == ["status", "n"]


# ---- public synthetic data -------------------------------------------------------------------------
def test_synthetic_frame_needs_no_files():
    ref = random_protein(80, 1)
    other = homolog(ref, 2, 0.2, deletions=[(30, 2)])
    fr = synthetic_frame({"fam": ref}, {"g2": {"fam": other}},
                         {"s1": ("anchor", "fam", ref), "s2": ("g2", "fam", other),
                          "s3": ("anchor", "fam", ref[:20] + "X" * 10 + ref[30:])})
    assert fr.rows("fam", "s1", [5, 50]).tolist() == [4, 49]
    assert fr.seq_pos("fam", "s2", 40) == 38                          # two residues deleted upstream, via the group map
    assert fr.rows("fam", "s2", [31]).isna().all()                    # a deleted position is a hole
    assert fr.rows("fam", "s3", [25]).isna().all() and fr.holes.iloc[-1]["status"] == Status.UNKNOWN_RESIDUE
    assert fr.sequence("s1") == ref and fr.groups == ["anchor", "g2"]


def test_random_protein_and_homolog_are_deterministic():
    assert random_protein(50, 3) == random_protein(50, 3) != random_protein(50, 4)
    s = random_protein(100, 5)
    assert homolog(s, 6, 0.3, deletions=[(40, 2)]) == homolog(s, 6, 0.3, deletions=[(40, 2)])
    h = homolog(s, 6, 0.5, deletions=[(40, 2)], insertions=[(70, "GG")], flank=6)
    assert len(h) == 100 - 2 + 2
    assert h[34:40] == s[34:40]                                       # indel flanks are never substituted


def test_write_mmcif_rejects_an_unknown_variant(tmp_path):
    with pytest.raises(ValueError, match="variant"):
        write_mmcif(tmp_path / "x.cif", "X", [ChainSpec("A", "1", "ACDE", [1, 2, 3, 4])], variant="9col")


def test_the_toy_project_is_self_consistent(toy):
    assert set(toy.seqs) >= {"a_kin_ref", "a_kin_del", "a_kin_ins", "a_kin_x", "a_kin_frag", "b_kin_ref"}
    assert len(toy.anchor_refs["kin"]) == 150 and len(toy.group_b_refs["kin"]) == 150 - 3 + 2
    assert toy.config.is_file() and (toy.root / "structures" / "toy_a.cif").is_file()
