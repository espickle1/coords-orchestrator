# `coords.toml` reference

Paths are relative to the directory holding the config file. Unknown keys are an error, so a typo is
reported instead of ignored. Only `[project]` and `[references]` are required.

```toml
[project]
name    = "my-project"     # free text, recorded in the manifest (default: the directory name)
anchor  = "A"              # required: the group whose references define the frame
out_dir = "maps"           # where the maps are written (default "maps")

# A group is a set of homologous references (a strain, a subtype, an orthologue). One FASTA per group.
# The first header token of each record is the family name; it becomes a directory name, so it must match
# [A-Za-z0-9][A-Za-z0-9_+-]*. The anchor group defines the families; every other group may cover a subset.
[references.A]
fasta = "refs/A.fasta"
label = "anchor strain"    # free text, recorded in the manifest

[references.B]
fasta = "refs/B.fasta"

# Sequences you will read per-residue tables for (optional).
[sequences]
table = "sequences.csv"

# mmCIF structures (optional). The file stem is the structure key, so stems must be unique.
[structures]
dir  = "structures"
glob = "**/*.cif"          # default

# Sites given in different numbering schemes (optional).
[sites]
table           = "sites.csv"
family_column   = "family"      # these four are the defaults
position_column = "position"
scheme_column   = "scheme"
residue_column  = "residue"     # optional in the table

[sites.schemes]            # scheme name -> offset added to `position`
full_length = 0
mature      = 16
tag         = "unresolved" # or omit the scheme entirely: both exclude the row

# Aligner (optional). Defaults shown.
[aligner]
matrix     = "BLOSUM62"    # any name Bio.Align.substitution_matrices.load accepts
open_gap   = -10.0
extend_gap = -0.5
end_gap    = 0.0           # 0 = free end gaps, so partial sequences are tolerated
x_score    = 1.0           # score of X against anything: a positional placeholder (see design.md)

# Quality gates (optional). Defaults shown. They flag; they never delete.
[thresholds]
min_chain_identity = 0.30  # best-reference identity below this: the chain is unregistered
low_chain_identity = 0.50  # composed identity to the anchor below this: chain flagged low_identity
min_chain_coverage = 0.40  # fraction of a chain's SEQRES that must align, else unregistered
min_chain_zscore   = 8.0   # alignment score vs shuffled copies of the chain; unrelated chains sit near 0
zscore_shuffles    = 20    # shuffled copies used for that z-score
min_seq_identity   = 0.60  # a sequence below this identity to its reference is `unreliable`
min_seq_aligned    = 30    # fewer known (non-X) aligned residues than this is `unreliable`
min_seq_coverage   = 0.80  # known aligned residues / reference length below this is `partial`
```

## Input tables

**`sequences.csv`** (`[sequences] table`). One row per sequence.

| Column | |
|---|---|
| `seq_id` | Required, unique. |
| `group`, `family` | Required. There must be a reference for that pair. |
| `sequence` | The sequence text, **or** |
| `fasta` | A FASTA path (relative to the config file). `record` names the record; the default is `seq_id`, then `family`. |
| `seq_hash` | Optional. The sha256 of the normalised sequence; the build refuses a mismatch. |
| anything else | Carried through to `maps/sequences.csv`. |

The hash is the sha256 of the sequence upper-cased, with whitespace and a trailing `*` removed. Store it next to
a per-residue table you computed, then `Frame.check_hash` proves the table still describes the string the map
was built from.

**`sites.csv`** (`[sites] table`). One row per site: `family`, `position` (an integer in whatever numbering
the row's `scheme` uses), `scheme`, optionally `residue`; any other columns are kept. A row is resolved to
`ref_pos = position + offset`. A row whose scheme is empty, `unresolved`, or not listed in `[sites.schemes]` gets
`unresolved_scheme`; a resolved position outside the anchor reference, or a family with no reference, gets
`no_ref_position`. When `residue` is given, `matches_ref` says whether it equals the anchor's residue there.

## What a build writes

See [design.md](design.md#what-a-build-writes). `python -m coords info <out_dir>` summarises it and lists every
flagged chain and sequence.
