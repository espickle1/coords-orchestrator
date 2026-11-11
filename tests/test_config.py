"""The config fails early and says what to change; nothing is guessed."""
from __future__ import annotations

from pathlib import Path

import pytest

from coords.config import ConfigError, load_config, load_references
from coords.seqio import normalize, read_fasta, sha256_seq, write_fasta


def _write(tmp_path: Path, toml: str, a="ACDEFGHIKL", b=None) -> Path:
    (tmp_path / "refs").mkdir(exist_ok=True)
    write_fasta(tmp_path / "refs" / "A.fasta", {"kin": a})
    write_fasta(tmp_path / "refs" / "B.fasta", {"kin": b or a})
    cfg = tmp_path / "coords.toml"
    cfg.write_text(toml, encoding="utf-8")
    return cfg


BASE = """[project]
anchor = "A"
[references.A]
fasta = "refs/A.fasta"
[references.B]
fasta = "refs/B.fasta"
"""


def test_valid_config_resolves_paths_against_its_own_directory(tmp_path):
    cfg = load_config(_write(tmp_path, BASE + '[structures]\ndir = "s"\n'))
    assert cfg.anchor == "A" and cfg.groups == ["A", "B"]
    assert cfg.references["B"].fasta == (tmp_path / "refs" / "B.fasta").resolve()
    assert cfg.structures.glob == "**/*.cif"
    assert cfg.out_dir == (tmp_path / "maps").resolve()
    assert cfg.sequences is None and cfg.sites is None


def test_anchor_is_always_first_in_group_order(tmp_path):
    toml = BASE.replace('anchor = "A"', 'anchor = "B"')
    assert load_config(_write(tmp_path, toml)).groups == ["B", "A"]


def test_sites_schemes_accept_integers_and_unresolved(tmp_path):
    cfg = load_config(_write(tmp_path, BASE + '[sites]\ntable = "t.csv"\n[sites.schemes]\nfull = 0\nmature = 16\ntag = "unresolved"\n'))
    assert dict(cfg.sites.schemes) == {"full": 0, "mature": 16, "tag": None}


@pytest.mark.parametrize("toml, fragment", [
    ('[references.A]\nfasta = "refs/A.fasta"\n', "anchor is required"),
    ('[project]\nanchor = "Z"\n[references.A]\nfasta = "refs/A.fasta"\n', "not one of the groups"),
    ('[project]\nanchor = "A"\n', "references"),
    ('[project]\nanchor = "A"\n[references.A]\nlabel = "x"\n', "needs fasta"),
    ('[project]\nanchor = "A"\n[references."a b"]\nfasta = "refs/A.fasta"\n', "must match"),
    ('[project]\nanchor = "A"\nbogus = 1\n[references.A]\nfasta = "refs/A.fasta"\n', "unknown key"),
    (BASE + '[sites]\ntable = "t.csv"\n[sites.schemes]\nbad = "later"\n', "integer offset"),
    (BASE + '[thresholds]\nnot_a_threshold = 1\n', "unknown key"),
    (BASE + '[structures]\n', "needs dir"),
    ("this is = not toml [", "not valid TOML"),
])
def test_bad_config_names_the_problem(tmp_path, toml, fragment):
    with pytest.raises(ConfigError, match=fragment):
        load_config(_write(tmp_path, toml))


def test_missing_config_file(tmp_path):
    with pytest.raises(ConfigError, match="not found"):
        load_config(tmp_path / "nope.toml")


def test_a_family_only_in_a_non_anchor_group_is_an_error(tmp_path):
    cfg = _write(tmp_path, BASE)
    write_fasta(tmp_path / "refs" / "B.fasta", {"kin": "ACDEFGHIKL", "extra": "MKTAYIAKQR"})
    with pytest.raises(ConfigError, match="anchor group 'A' does not"):
        load_references(load_config(cfg))


def test_a_family_name_that_cannot_be_a_filename_is_an_error(tmp_path):
    cfg = _write(tmp_path, BASE)
    (tmp_path / "refs" / "A.fasta").write_text(">a/b\nACDE\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="record id"):
        load_references(load_config(cfg))


def test_missing_reference_file(tmp_path):
    cfg = _write(tmp_path, BASE)
    (tmp_path / "refs" / "B.fasta").unlink()
    with pytest.raises(ConfigError, match="fasta not found"):
        load_references(load_config(cfg))


def test_seqio_normalisation_and_hash():
    assert normalize(" ac de*\n") == "ACDE"
    assert sha256_seq("ACDE") == sha256_seq(normalize("acde*"))
    assert len(sha256_seq("A")) == 64


def test_read_fasta_keys_by_first_token_and_keeps_first_duplicate(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">one some description\nAC\nDE\n>two\nKK\n>one\nWW\n", encoding="utf-8")
    assert read_fasta(p) == {"one": "ACDE", "two": "KK"}
