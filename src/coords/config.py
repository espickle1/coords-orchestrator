"""Project configuration (``coords.toml``) and the reference sequences it names.

A project says: which sequences define the frame (references, one FASTA per *group*, one record per
*family*), which group is the anchor, where the sequences to be read through the maps live, where the
structures live, and optionally a table of sites whose numbering schemes are to be resolved.

    [project]
    name   = "my-project"
    anchor = "A"              # the group whose references define the frame
    out_dir = "maps"          # relative to this file

    [references.A]            # one FASTA per group; the first header token of each record is its family
    fasta = "refs/A.fasta"
    label = "strain A"        # free text, recorded in the manifest

    [references.B]
    fasta = "refs/B.fasta"

    [sequences]               # optional: sequences to be read through the maps
    table = "sequences.csv"   # columns: seq_id, group, family, and sequence | fasta [+ record] [+ seq_hash]

    [structures]              # optional: mmCIF files (protein chains)
    dir  = "structures"
    glob = "**/*.cif"

    [sites]                   # optional: sites whose numbering scheme must be resolved
    table = "sites.csv"       # columns: family, position, scheme [+ residue] [+ anything else, kept]
    [sites.schemes]           # scheme name -> offset added to ``position``; "unresolved" excludes the row
    full_length = 0
    mature      = 16

    [aligner]                 # optional, defaults shown in coords.align.AlignerParams
    [thresholds]              # optional, defaults shown in Thresholds

A *group* is a set of homologous references (a strain, a subtype, an orthologue); a *family* is one
protein within it. Positions are always reference positions on the anchor group's reference of a family.
Paths are relative to the directory holding the config file.
"""
from __future__ import annotations

import re
import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple, Union

from .align import AlignerParams
from .seqio import read_fasta

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_+\-]*$")     # become directory names and CSV keys
UNRESOLVED = "unresolved"


class ConfigError(ValueError):
    """The configuration or the files it names are not usable. The message says what to change."""


@dataclass(frozen=True)
class Thresholds:
    """Quality gates. They flag, they do not delete: a flagged chain or sequence is still mapped."""

    min_chain_identity: float = 0.30    # best reference identity below this: the chain is unregistered
    low_chain_identity: float = 0.50    # composed identity to the anchor below this: chain flagged low_identity
    min_chain_coverage: float = 0.40    # fraction of a chain's SEQRES that must align, else unregistered
    min_chain_zscore: float = 8.0       # alignment score vs shuffled copies of the chain; unrelated chains sit near 0
    zscore_shuffles: int = 20           # shuffled copies used for that z-score
    min_seq_identity: float = 0.60      # sequence identity to its reference below this: unreliable
    min_seq_aligned: int = 30           # known (non-X) aligned residues below this: unreliable
    min_seq_coverage: float = 0.80      # known aligned residues / reference length below this: partial


@dataclass(frozen=True)
class ReferenceSpec:
    group: str
    fasta: Path
    label: str = ""


@dataclass(frozen=True)
class SequencesSpec:
    table: Path


@dataclass(frozen=True)
class StructuresSpec:
    dir: Path
    glob: str = "**/*.cif"


@dataclass(frozen=True)
class SitesSpec:
    table: Path
    family_column: str = "family"
    position_column: str = "position"
    scheme_column: str = "scheme"
    residue_column: str = "residue"
    schemes: Mapping[str, Optional[int]] = field(default_factory=dict)    # name -> offset, None = unresolved


@dataclass
class Config:
    path: Path
    base: Path
    name: str
    anchor: str
    out_dir: Path
    references: Dict[str, ReferenceSpec]
    sequences: Optional[SequencesSpec] = None
    structures: Optional[StructuresSpec] = None
    sites: Optional[SitesSpec] = None
    aligner: AlignerParams = AlignerParams()
    thresholds: Thresholds = Thresholds()

    @property
    def groups(self) -> List[str]:
        """Groups with the anchor first: on a tie a chain is assigned to the anchor."""
        return [self.anchor] + [g for g in self.references if g != self.anchor]


def _path(base: Path, value: str) -> Path:
    p = Path(value)
    return p if p.is_absolute() else (base / p).resolve()


def _section(raw: Mapping, key: str) -> Mapping:
    v = raw.get(key, {})
    if not isinstance(v, dict):
        raise ConfigError(f"[{key}] must be a table")
    return v


def _check_keys(section: str, given: Mapping, allowed: set) -> None:
    extra = sorted(set(given) - allowed)
    if extra:
        raise ConfigError(f"[{section}] has unknown key(s) {extra}; allowed: {sorted(allowed)}")


def load_config(path: Union[str, Path]) -> Config:
    path = Path(path).resolve()
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    try:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path}: not valid TOML: {e}") from e
    base = path.parent
    _check_keys("<top level>", raw, {"project", "references", "sequences", "structures", "sites",
                                     "aligner", "thresholds"})

    proj = _section(raw, "project")
    _check_keys("project", proj, {"name", "anchor", "out_dir"})
    anchor = proj.get("anchor")
    if not anchor:
        raise ConfigError("[project] anchor is required: the group whose references define the frame")

    refs_raw = _section(raw, "references")
    if not refs_raw:
        raise ConfigError("[references.<group>] is required: at least the anchor group")
    references: Dict[str, ReferenceSpec] = {}
    for group, spec in refs_raw.items():
        if not NAME_RE.match(group):
            raise ConfigError(f"group name {group!r} must match {NAME_RE.pattern} (it becomes a directory name)")
        if not isinstance(spec, dict) or "fasta" not in spec:
            raise ConfigError(f"[references.{group}] needs fasta = \"path\"")
        _check_keys(f"references.{group}", spec, {"fasta", "label"})
        references[group] = ReferenceSpec(group, _path(base, spec["fasta"]), str(spec.get("label", "")))
    if anchor not in references:
        raise ConfigError(f"anchor {anchor!r} is not one of the groups {sorted(references)}")

    sequences = structures = sites = None
    if "sequences" in raw:
        s = _section(raw, "sequences")
        _check_keys("sequences", s, {"table"})
        if "table" not in s:
            raise ConfigError("[sequences] needs table = \"path/to/sequences.csv\"")
        sequences = SequencesSpec(_path(base, s["table"]))
    if "structures" in raw:
        s = _section(raw, "structures")
        _check_keys("structures", s, {"dir", "glob"})
        if "dir" not in s:
            raise ConfigError("[structures] needs dir = \"path/to/structures\"")
        structures = StructuresSpec(_path(base, s["dir"]), str(s.get("glob", "**/*.cif")))
    if "sites" in raw:
        s = dict(_section(raw, "sites"))
        schemes_raw = s.pop("schemes", {})
        _check_keys("sites", s, {"table", "family_column", "position_column", "scheme_column", "residue_column"})
        if "table" not in s:
            raise ConfigError("[sites] needs table = \"path/to/sites.csv\"")
        schemes: Dict[str, Optional[int]] = {}
        for name, off in schemes_raw.items():
            if off == UNRESOLVED:
                schemes[name] = None
            elif isinstance(off, int) and not isinstance(off, bool):
                schemes[name] = off
            else:
                raise ConfigError(f"[sites.schemes] {name} = {off!r}: use an integer offset or \"{UNRESOLVED}\"")
        sites = SitesSpec(_path(base, s["table"]), **{k: v for k, v in s.items() if k != "table"},
                          schemes=schemes)

    def build(cls, section):
        s = _section(raw, section)
        _check_keys(section, s, {f.name for f in fields(cls)})
        try:
            return cls(**s)
        except TypeError as e:
            raise ConfigError(f"[{section}]: {e}") from e

    return Config(path=path, base=base, name=str(proj.get("name", path.parent.name)), anchor=anchor,
                  out_dir=_path(base, proj.get("out_dir", "maps")), references=references,
                  sequences=sequences, structures=structures, sites=sites,
                  aligner=build(AlignerParams, "aligner"), thresholds=build(Thresholds, "thresholds"))


def load_references(cfg: Config) -> Dict[Tuple[str, str], str]:
    """``{(group, family): sequence}`` for every reference, validated against the anchor.

    The anchor group defines the families. Every other group may cover a subset of them; a family a
    group has that the anchor lacks has no frame to map into and is an error.
    """
    refs: Dict[Tuple[str, str], str] = {}
    for group in cfg.groups:
        spec = cfg.references[group]
        if not spec.fasta.is_file():
            raise ConfigError(f"[references.{group}] fasta not found: {spec.fasta}")
        recs = read_fasta(spec.fasta)
        if not recs:
            raise ConfigError(f"{spec.fasta}: no sequences")
        for family, seq in recs.items():
            if not NAME_RE.match(family):
                raise ConfigError(f"{spec.fasta}: record id {family!r} must match {NAME_RE.pattern} "
                                  "(the first header token is the family name)")
            if not seq:
                raise ConfigError(f"{spec.fasta}: record {family!r} is empty")
            refs[(group, family)] = seq
    anchor_fams = {f for (g, f) in refs if g == cfg.anchor}
    for (g, f) in refs:
        if f not in anchor_fams:
            raise ConfigError(f"group {g!r} has family {f!r} but the anchor group {cfg.anchor!r} does not: "
                              "there is no frame to map it into. Add it to the anchor's FASTA.")
    return refs
